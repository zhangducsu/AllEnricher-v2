"""Core enrichment-analysis engine for AllEnricher.

This module implements over-representation analysis, GSEA, ssGSEA, multiple-
testing correction, input loading, result filtering, and database-level
workflow orchestration."""

from __future__ import annotations

import os
import csv
import logging
import collections
from abc import ABC, abstractmethod
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Set, Tuple, Any
from pathlib import Path
import math
import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests
from concurrent.futures import ThreadPoolExecutor, as_completed
from tqdm import tqdm

from allenricher.core.config import Config, EnrichmentMethod, CorrectionMethod

# Library modules obtain named loggers and leave handler configuration to the CLI.
logger = logging.getLogger(__name__)


def generate_term_url(term_id: str, database: str) -> str:
    """Build the canonical external URL for a database term.
    
    Args:
        term_id: Stable database identifier such as ``GO:0008150``.
        database: Database name.
    
    Returns:
        The external term URL, or an empty string for unsupported databases."""
    # Normalize database names for URL routing.
    database = database.upper()

    if database == "GO":
        return f"https://amigo.geneontology.org/amigo/term/{term_id}"
    elif database == "KEGG":
        return f"https://www.kegg.jp/entry/{term_id}"
    elif database == "REACTOME":
        return f"https://reactome.org/PathwayBrowser/#{term_id}"
    elif database == "WIKIPATHWAYS":
        return f"https://www.wikipathways.org/index.php/Pathway:{term_id}"
    elif database == "DO":
        return f"https://disease-ontology.org/?id={term_id}"
    elif database == "DISGENET":
        return f"https://www.disgenet.org/browser/0/10/{term_id}/"
    else:
        return ""


def add_result_term_metadata(
    frame: pd.DataFrame,
    term_data: Dict[str, Dict[str, Any]],
) -> pd.DataFrame:
    """Add stable term IDs, display names, and available hierarchy metadata."""
    result = frame.copy()
    term_lookup = {str(term_id): info for term_id, info in term_data.items()}

    id_column = next(
        (column for column in ("Term_ID", "pathway", "term_id", "ID", "id")
         if column in result.columns),
        None,
    )
    if id_column is None:
        term_ids = pd.Series(result.index.astype(str), index=result.index)
        result.index = term_ids.to_numpy()
        result.index.name = "Term_ID"
    else:
        term_ids = result[id_column].astype(str)

    names = []
    hierarchies = []
    for term_id in term_ids:
        info = term_lookup.get(str(term_id), {})
        raw_name = str(info.get("name") or term_id).strip()
        hierarchy = str(info.get("hierarchy") or "").strip()
        if "|" in raw_name:
            hierarchy = raw_name
            raw_name = raw_name.rsplit("|", 1)[-1].strip()
        names.append(raw_name or str(term_id))
        hierarchies.append(hierarchy)

    if id_column is not None:
        if "Term_ID" in result.columns:
            result["Term_ID"] = term_ids.to_numpy()
        else:
            result.insert(0, "Term_ID", term_ids.to_numpy())

    if "Term_Name" in result.columns:
        result["Term_Name"] = names
    else:
        position = result.columns.get_loc("Term_ID") + 1 if "Term_ID" in result.columns else 0
        result.insert(position, "Term_Name", names)

    if any(hierarchies):
        if "Hierarchy" in result.columns:
            result["Hierarchy"] = hierarchies
        else:
            result.insert(result.columns.get_loc("Term_Name") + 1, "Hierarchy", hierarchies)

    return result


@dataclass
class EnrichmentResult:
    """Statistics and metadata for one enriched term."""
    term_id: str
    term_name: str
    database: str
    pvalue: float
    adjusted_pvalue: float
    gene_count: int
    background_count: int = 0
    expected_count: float = 0.0
    rich_factor: float = 0.0
    gene_list: List[str] = field(default_factory=list)
    gene_ratio: str = ""
    background_ratio: str = ""
    term_url: str = ""

    # GSEA-specific fields
    nes: Optional[float] = None
    es: Optional[float] = None
    fdr: Optional[float] = None
    leading_edge: Optional[List[str]] = None
    set_size: Optional[int] = None
    rank_at_max: Optional[int] = None
    fwerp: Optional[float] = None          # Old version of single-circuit objects compatible fields
    tag_pct: str = ""
    gene_pct: str = ""
    
    def to_dict(self) -> Dict[str, Any]:
        """Serialize one result using the public ORA or GSEA column contract."""
        is_gsea = self.nes is not None
        
        result = collections.OrderedDict() if is_gsea else {}
        
        if is_gsea:
            result["Term_ID"] = self.term_id
            result["Term_Name"] = self.term_name
            result["Database"] = self.database
            result["setSize"] = self.set_size if self.set_size is not None else self.gene_count
            result["ES"] = round(self.es, 4) if self.es is not None else 0.0
            result["NES"] = round(self.nes, 4) if self.nes is not None else 0.0
            result["p_value"] = self.pvalue
            result["FDR"] = self.fdr if self.fdr is not None else self.adjusted_pvalue
            result["rank"] = self.rank_at_max if self.rank_at_max is not None else 0
            result["Tag %"] = self.tag_pct if self.tag_pct else ""
            result["Gene %"] = self.gene_pct if self.gene_pct else ""
            result["Lead_genes"] = ";".join(self.leading_edge) if self.leading_edge else ""
            result["matched_genes"] = ";".join(self.gene_list)
            result["Term_URL"] = self.term_url
        else:
            result["Term_ID"] = self.term_id
            result["Term_Name"] = self.term_name
            result["Database"] = self.database
            result["P_Value"] = self.pvalue
            result["Adjusted_P_Value"] = self.adjusted_pvalue
            result["Gene_Count"] = self.gene_count
            result["Background_Count"] = self.background_count
            result["Expected_Count"] = round(self.expected_count, 4)
            result["Rich_Factor"] = round(self.rich_factor, 4)
            result["Gene_Ratio"] = self.gene_ratio
            result["Background_Ratio"] = self.background_ratio
            result["Term_URL"] = self.term_url
            result["Genes"] = ";".join(self.gene_list)
            
        return result


class EnrichmentMethodBase(ABC):
    """Abstract interface implemented by enrichment methods."""
    
    @abstractmethod
    def calculate_pvalue(
        self,
        gene_count: int,
        background_count: int,
        gene_total: int,
        background_total: int
    ) -> float:
        """Calculate the one-sided enrichment P value for a contingency table."""
        pass
    
    @abstractmethod
    def calculate_enrichment(
        self,
        gene_set: Set[str],
        background_set: Set[str],
        term_genes: Set[str],
        term_name: str,
        term_id: str,
        database: str,
        background_total: Optional[int] = None
    ) -> Optional[EnrichmentResult]:
        """Analyze one term and return ``None`` when it has no query overlap."""
        pass


class FisherExactTest(EnrichmentMethodBase):
    """One-sided Fisher exact test for over-representation analysis."""
    
    def calculate_pvalue(
        self,
        gene_count: int,
        background_count: int,
        gene_total: int,
        background_total: int
    ) -> float:
        """Calculate the one-sided Fisher exact-test P value."""
        # Build the 2 x 2 contingency table.
        #                         In term    Outside term
        # Query genes                a            b
        # Background only            c            d
        
        a = gene_count
        b = gene_total - gene_count
        c = background_count - gene_count
        d = background_total - background_count - gene_total + gene_count
        
        # Clamp invalid negative cells defensively when external metadata are
        # inconsistent with the actual annotation universe.
        if a < 0: a = 0
        if b < 0: b = 0
        if c < 0: c = 0
        if d < 0: d = 0
        
        # A fully empty table contains no evidence of enrichment.
        if a + b + c + d == 0:
            return 1.0
        
        # Test whether the term is over-represented in the query genes.
        _, pvalue = stats.fisher_exact([[a, b], [c, d]], alternative='greater')
        return pvalue
    
    def calculate_enrichment(
        self,
        gene_set: Set[str],
        background_set: Set[str],
        term_genes: Set[str],
        term_name: str,
        term_id: str,
        database: str,
        background_total: Optional[int] = None
    ) -> Optional[EnrichmentResult]:
        """Calculate Fisher ORA statistics for one database term."""
        # Query genes annotated to this term.
        genes_in_term = gene_set & term_genes
        gene_count = len(genes_in_term)
        
        # Terms without query overlap do not produce result rows.
        if gene_count < 1:
            return None
        
        # Background genes annotated to this term.
        background_in_term = background_set & term_genes
        background_count = len(background_in_term)
        
        # Resolve query and background universe sizes.
        gene_total = len(gene_set)
        if background_total is None:
            background_total = len(background_set)
        
        # Expected overlap under random sampling from the background.
        expected_count = (background_count / background_total) * gene_total if background_total > 0 else 0
        
        # Rich factor is observed overlap divided by expected overlap.
        rich_factor = gene_count / expected_count if expected_count > 0 else 0
        
        # Calculate the one-sided enrichment P value.
        pvalue = self.calculate_pvalue(
            gene_count, background_count, gene_total, background_total
        )
        
        # Attach a canonical external term link when supported.
        term_url = generate_term_url(term_id, database)
        
        # adjusted_pvalue is initialized here and replaced during correction.
        return EnrichmentResult(
            term_id=term_id,
            term_name=term_name,
            database=database,
            pvalue=pvalue,
            adjusted_pvalue=pvalue,
            gene_count=gene_count,
            background_count=background_count,
            expected_count=expected_count,
            rich_factor=rich_factor,
            gene_list=list(genes_in_term),
            gene_ratio=f"{gene_count}/{gene_total}",
            background_ratio=f"{background_count}/{background_total}",
            term_url=term_url,
        )


class HypergeometricTest(EnrichmentMethodBase):
    """One-sided hypergeometric over-representation analysis.
    
    The survival function computes the probability of observing at least the
    measured overlap under random sampling without replacement."""
    
    def calculate_pvalue(
        self,
        gene_count: int,
        background_count: int,
        gene_total: int,
        background_total: int
    ) -> float:
        """Calculate the upper-tail hypergeometric P value."""
        # P(X >= k) = 1 - P(X < k) = 1 - P(X <= k-1)
        # The survival function evaluates the upper tail directly.
        
        M = background_total
        n = background_count
        N = gene_total
        k = gene_count
        
        # sf(k - 1, M, n, N) equals P(X >= k).
        pvalue = stats.hypergeom.sf(k - 1, M, n, N)
        return pvalue
    
    def calculate_enrichment(
        self,
        gene_set: Set[str],
        background_set: Set[str],
        term_genes: Set[str],
        term_name: str,
        term_id: str,
        database: str,
        background_total: Optional[int] = None
    ) -> Optional[EnrichmentResult]:
        """Calculate hypergeometric ORA statistics for one database term."""
        # Reuse the shared overlap and effect-size calculations.
        fisher = FisherExactTest()
        result = fisher.calculate_enrichment(
            gene_set, background_set, term_genes, term_name, term_id, database,
            background_total=background_total
        )
        
        if result is None:
            return None
        
        # Replace the temporary Fisher P value with the hypergeometric tail.
        effective_background = background_total if background_total is not None else len(background_set)
        pvalue = self.calculate_pvalue(
            result.gene_count,
            result.background_count,
            len(gene_set),
            effective_background
        )
        result.pvalue = pvalue
        result.adjusted_pvalue = pvalue
        
        return result


class GSEA(EnrichmentMethodBase):
    """Gene Set Enrichment Analysis for a preranked gene list.
    
    Running enrichment scores are calculated deterministically in Python. The
    production batch workflow delegates statistical inference to Bioconductor
    ``fgsea`` through :class:`EnrichmentAnalyzer`."""
    
    def __init__(self, permutations: int = 1000, min_size: int = 15, max_size: int = 500, seed: int = 42):
        """Configure gene-set size limits, permutation count, and random seed."""
        self.permutations = permutations
        self.min_size = min_size
        self.max_size = max_size
        self.seed = seed
    
    def calculate_pvalue(
        self,
        gene_count: int,
        background_count: int,
        gene_total: int,
        background_total: int
    ) -> float:
        """Compatibility stub; GSEA P values require a ranked gene list."""
        # Nominal significance requires the ranked list and is calculated later.
        return 1.0
    
    def calculate_enrichment_score(
        self,
        ranked_genes: List[str],
        gene_set: Set[str],
        gene_weights: Optional[Dict[str, float]] = None,
        return_rank: bool = False,
    ):
        """Calculate the running enrichment score and leading-edge statistics."""
        n = len(ranked_genes)
        nh = len(gene_set & set(ranked_genes))

        if nh == 0:
            return (0.0, [], 0) if return_rank else (0.0, [])

        hits = gene_set & set(ranked_genes)
        running_sum = 0.0
        max_es = min_es = 0.0
        max_rank = min_rank = 0

        nr = sum(abs(gene_weights.get(g, 1.0)) for g in hits) if gene_weights else nh
        hit_inc = 1.0 / nr if nr > 0 else 0
        miss_inc = 1.0 / (n - nh) if (n - nh) > 0 else 0

        for i, gene in enumerate(ranked_genes):
            if gene in gene_set:
                weight = abs(gene_weights.get(gene, 1.0)) if gene_weights else 1.0
                running_sum += hit_inc * weight
                if running_sum > max_es:
                    max_es = running_sum
                    max_rank = i + 1
            else:
                running_sum -= miss_inc
                if running_sum < min_es:
                    min_es = running_sum
                    min_rank = i + 1

        if abs(max_es) >= abs(min_es):
            es, rank_at_max = max_es, max_rank
            leading_edge = [g for g in ranked_genes[:rank_at_max] if g in hits]
        else:
            es, rank_at_max = min_es, min_rank
            leading_edge = [g for g in ranked_genes[rank_at_max:] if g in hits]
        # Protect the theoretical [-1, 1] range from floating-point drift.
        es = max(-1.0, min(1.0, float(es)))
        return (es, leading_edge, rank_at_max) if return_rank else (es, leading_edge)

    def _run_permutation_test(
        self,
        ranked_genes: List[str],
        gene_set: Set[str],
        observed_es: float,
        gene_weights: Optional[Dict[str, float]] = None,
        n_permutations: int = 1000,
        seed: int = 42
    ) -> float:
        """Estimate nominal significance from deterministic gene-label permutations."""
        rng = np.random.default_rng(seed)
        ranked_array = np.array(ranked_genes)
        count_ge = 0

        for _ in range(n_permutations):
            # Permute gene labels while preserving the ranked positions.
            permuted_genes = rng.permutation(ranked_array).tolist()
            # Only the permuted ES is required for the null distribution.
            permuted_es, _, _ = self.calculate_enrichment_score(
                permuted_genes, gene_set, gene_weights, return_rank=True
            )
            if ((observed_es >= 0 and permuted_es >= observed_es) or
                    (observed_es < 0 and permuted_es <= observed_es)):
                count_ge += 1

        # Add one to numerator and denominator to avoid a zero P value.
        pvalue = (count_ge + 1) / (n_permutations + 1)
        return pvalue

    def calculate_normalized_es(
        self,
        ranked_genes: List[str],
        gene_set: Set[str],
        gene_weights: Optional[Dict[str, float]] = None,
        return_rank: bool = False,
    ):
        """Normalize an enrichment score against permutations with the same sign."""
        # Calculate the observed ES and leading edge.
        es, leading_edge, rank_at_max = self.calculate_enrichment_score(
            ranked_genes, gene_set, gene_weights, return_rank=True
        )

        # A gene set without ranked-list overlap has a neutral score.
        if es == 0.0:
            result = (0.0, 0.0, 1.0, [])
            return result + (0,) if return_rank else result

        # Build sign-specific permutation null distributions.
        rng = np.random.default_rng(self.seed)
        ranked_array = np.array(ranked_genes)
        null_es_positive = []
        null_es_negative = []

        for _ in range(self.permutations):
            permuted_genes = rng.permutation(ranked_array).tolist()
            permuted_es, _, _ = self.calculate_enrichment_score(
                permuted_genes, gene_set, gene_weights, return_rank=True
            )

            # NES normalization compares scores with null values of the same sign.
            if permuted_es >= 0:
                null_es_positive.append(permuted_es)
            else:
                null_es_negative.append(permuted_es)

        # Mean absolute null ES for each direction.
        mean_pos = np.mean(null_es_positive) if null_es_positive else 1.0
        mean_neg = np.mean([abs(e) for e in null_es_negative]) if null_es_negative else 1.0

        # Normalize the observed ES by the matching null mean.
        if es > 0:
            nes = es / mean_pos if mean_pos > 0 else 0.0
        else:
            nes = es / mean_neg if mean_neg > 0 else 0.0

        # Calculate the sign-specific nominal P value.
        if es > 0:
            count_ge = sum(value >= es for value in null_es_positive)
            pvalue = (count_ge + 1) / (len(null_es_positive) + 1)
        else:
            count_ge = sum(value <= es for value in null_es_negative)
            pvalue = (count_ge + 1) / (len(null_es_negative) + 1)

        result = (es, nes, pvalue, leading_edge)
        return result + (rank_at_max,) if return_rank else result

    def calculate_enrichment(
        self,
        gene_set: Set[str],
        background_set: Set[str],
        term_genes: Set[str],
        term_name: str,
        term_id: str,
        database: str,
        ranked_genes: Optional[List[str]] = None,
        gene_weights: Optional[Dict[str, float]] = None
    ) -> Optional[EnrichmentResult]:
        """Analyze one gene set against a preranked gene list."""
        # Apply size limits after intersecting with the analyzed genes.
        overlap = gene_set & term_genes
        set_size = len(overlap)
        if set_size < self.min_size or set_size > self.max_size:
            return None

        # Legacy callers may omit ranks; preserve their background-order fallback.
        if ranked_genes is None:
            ranked_genes = list(background_set)

        # Calculate ES, NES, nominal P value, and leading-edge metadata.
        es, nes, pvalue, leading_edge, rank_at_max = self.calculate_normalized_es(
            ranked_genes, term_genes, gene_weights, return_rank=True
        )

        # The analyzer replaces this initial value during multiple-test correction.
        fdr = pvalue

        term_url = generate_term_url(term_id, database)

        return EnrichmentResult(
            term_id=term_id,
            term_name=term_name,
            database=database,
            pvalue=pvalue,
            adjusted_pvalue=pvalue,
            gene_count=set_size,
            gene_list=list(overlap),
            term_url=term_url,
            nes=nes,
            es=es,
            fdr=fdr,
            leading_edge=leading_edge,
            set_size=set_size,
            rank_at_max=rank_at_max,
        )

    def analyze_matrix(
        self,
        expression_matrix: pd.DataFrame,
        gene_sets: Dict[str, Set[str]],
        gene_weights_matrix: Optional[pd.DataFrame] = None
    ) -> pd.DataFrame:
        """Run GSEA for a gene-set collection and return a result table."""
        samples = expression_matrix.columns.tolist()
        results = {}

        for pathway_name, pathway_genes in gene_sets.items():
            nes_values = []
            for sample in samples:
                # Rank genes by expression within the current sample.
                sample_expr = expression_matrix[sample]
                ranked_genes = sample_expr.sort_values(ascending=False).index.tolist()

                # Use sample-specific weights when supplied.
                weights = None
                if gene_weights_matrix is not None and sample in gene_weights_matrix.columns:
                    weights = gene_weights_matrix[sample].to_dict()

                # Calculate the normalized enrichment score.
                _, nes, _, _ = self.calculate_normalized_es(
                    ranked_genes, pathway_genes, weights
                )
                nes_values.append(nes)

            results[pathway_name] = nes_values

        return pd.DataFrame(results, index=samples).T


class SSGSEA(EnrichmentMethodBase):
    """Single-sample Gene Set Enrichment Analysis.
    
    Scores summarize the relative enrichment of each gene set within one sample.
    The production matrix workflow uses Bioconductor GSVA when available."""

    def __init__(self, min_size: int = 1, max_size: Optional[int] = None):
        """Configure gene-set size limits and the ssGSEA weighting exponent."""
        self.method_name = "ssgsea"
        self.min_size = min_size
        self.max_size = max_size

    def calculate_pvalue(
        self,
        gene_count: int,
        background_count: int,
        gene_total: int,
        background_total: int
    ) -> float:
        """Compatibility stub; ssGSEA produces activity scores rather than P values."""
        # ssGSEA produces activity scores rather than inferential P values.
        return float('nan')

    def calculate_enrichment_score(
        self,
        ranked_genes: List[str],
        gene_set: Set[str],
        gene_weights: Optional[Dict[str, float]] = None
    ) -> Tuple[float, float, float, List[str]]:
        """Calculate an ssGSEA activity score for one ranked sample."""
        n = len(ranked_genes)
        nh = len(gene_set & set(ranked_genes))

        # A gene set without overlap has a neutral activity score.
        if nh == 0:
            return 0.0, 0.0, 0.0, []

        running_sum = 0.0
        max_es = 0.0
        min_es = 0.0
        leading_edge = []

        # Normalize hit and miss increments across the ranked list.
        nr = sum(abs(gene_weights.get(g, 1.0)) for g in gene_set if g in gene_weights) if gene_weights else nh
        hit_inc = 1.0 / nr if nr > 0 else 0
        miss_inc = 1.0 / (n - nh) if (n - nh) > 0 else 0

        # Traverse the ranked list and update the running score.
        for i, gene in enumerate(ranked_genes):
            if gene in gene_set:
                weight = abs(gene_weights.get(gene, 1.0)) if gene_weights else 1.0
                running_sum += hit_inc * weight
                if running_sum > max_es:
                    max_es = running_sum
                    leading_edge = ranked_genes[:i+1]
            else:
                running_sum -= miss_inc
                if running_sum < min_es:
                    min_es = running_sum

        # Retain only gene-set members in the leading-edge segment.
        leading_edge_genes = [g for g in leading_edge if g in gene_set]
        return max_es, min_es, max_es, leading_edge_genes

    def calculate_enrichment(
        self,
        gene_set: Set[str],
        background_set: Set[str],
        term_genes: Set[str],
        term_name: str,
        term_id: str,
        database: str,
        ranked_genes: Optional[List[str]] = None,
        gene_weights: Optional[Dict[str, float]] = None
    ) -> Optional[EnrichmentResult]:
        """Analyze one gene set for one ranked sample."""
        # Apply size limits after intersecting with the analyzed genes.
        overlap = gene_set & term_genes
        if len(overlap) < self.min_size or (
            self.max_size is not None and len(overlap) > self.max_size
        ):
            return None

        # Preserve the legacy background-order fallback when ranks are absent.
        if ranked_genes is None:
            ranked_genes = list(background_set)

        # Calculate the activity score and leading-edge genes.
        es, es_min, es_max, leading_edge = self.calculate_enrichment_score(
            ranked_genes, term_genes, gene_weights
        )

        # Normalize by the observed running-score range.
        denominator = abs(es_min) + abs(es_max)
        nes = es / denominator if denominator > 0 else 0.0

        # ssGSEA does not produce P values or FDR estimates.
        pvalue = float('nan')
        fdr = float('nan')

        term_url = generate_term_url(term_id, database)

        return EnrichmentResult(
            term_id=term_id,
            term_name=term_name,
            database=database,
            pvalue=pvalue,
            adjusted_pvalue=float('nan'),
            gene_count=len(overlap),
            background_count=len(term_genes),
            expected_count=0,
            rich_factor=0,
            gene_list=list(overlap),
            gene_ratio=f"{len(overlap)}/{len(gene_set)}",
            background_ratio=f"{len(term_genes)}/{len(background_set)}",
            term_url=term_url,
            nes=nes,
            es=es,
            fdr=fdr,
            leading_edge=leading_edge,
        )

    def analyze_matrix(
        self,
        expression_matrix: pd.DataFrame,
        gene_sets: Dict[str, Set[str]]
    ) -> pd.DataFrame:
        """Calculate a gene-set-by-sample activity matrix."""
        from allenricher.core.bioconductor import run_gsva

        return run_gsva(
            expression_matrix,
            gene_sets,
            method="ssgsea",
            tau=0.25,
            min_size=self.min_size,
            max_size=self.max_size,
        )


class EnrichmentAnalyzer:
    """Coordinate input loading, database analysis, correction, and output.
    
    The analyzer provides a shared workflow for ORA, GSEA, ssGSEA, and GSVA while
    preserving method-specific result schemas."""
    
    def __init__(self, config: Config):
        """Create an analyzer from a validated configuration."""
        self.config = config  # Global Configure Objects
        self._method: Optional[EnrichmentMethodBase] = None  # Examples of analytical methods that delayed initialization
        self.results: Dict[str, Any] = {}  # OLA is the list of result objects, GSEA is the fgsea official table

    @staticmethod
    def _is_tf_database(database: str) -> bool:
        return database.upper() in {'TRRUST', 'CHEA3', 'HTFTARGET', 'ANIMALTFDB'}

    def _filter_tf_gene_sets(
        self,
        database: str,
        method: str,
        gene_sets: Dict[str, Set[str]],
        universe: Set[str],
    ) -> Dict[str, Set[str]]:
        from allenricher.analysis.tf_enrichment import (
            filter_tf_gene_sets,
            resolve_tf_size_limits,
        )

        min_size, max_size = resolve_tf_size_limits(
            method,
            getattr(self.config, 'tf_min_size', None),
            getattr(self.config, 'tf_max_size', None),
        )
        filtered, stats = filter_tf_gene_sets(
            gene_sets, universe, min_size, max_size
        )
        logger.info(
            "TF %s/%s gene-set size filter: min=%d, max=%s, before=%d, "
            "after=%d, no_overlap=%d, below_min=%d, above_max=%d",
            database,
            method,
            min_size,
            'Inf' if max_size is None else max_size,
            stats['before'], stats['after'], stats['no_overlap'],
            stats['below_min'], stats['above_max'],
        )
        return filtered
    
    @property
    def method(self) -> EnrichmentMethodBase:
        """Return the configured enrichment-method instance."""
        if self._method is None:
            self._method = self._get_method()
        return self._method
    
    @method.setter
    def method(self, value: EnrichmentMethodBase):
        """Return the configured enrichment-method instance."""
        self._method = value
    
    def _get_method(self) -> EnrichmentMethodBase:
        """Instantiate the method selected in the configuration."""
        gsea_min_size = 15 if self.config.gsea_min_size is None else self.config.gsea_min_size
        gsea_max_size = 500 if self.config.gsea_max_size is None else self.config.gsea_max_size
        activity_min_size = 1 if self.config.gsea_min_size is None else self.config.gsea_min_size
        activity_max_size = self.config.gsea_max_size

        # Map of the name of the method to the example
        methods = {
            EnrichmentMethod.HYPERGEOMETRIC.value: HypergeometricTest(),  # Super-Geographic Test (ORA default method)
            EnrichmentMethod.GSEA.value: GSEA(  # GSEA Gene Set Enrichment Analysis
                permutations=self.config.gsea_permutations,
                min_size=gsea_min_size,
                max_size=gsea_max_size
            ),
            EnrichmentMethod.SSGSEA.value: SSGSEA(  # Single-sample gene-set enrichment analysis.
                min_size=activity_min_size,
                max_size=activity_max_size
            ),
        }

        # Delaying the import of GSVA to avoid recycling dependency (gsva.py imported the base class of enrichment.py)
        from allenricher.core.gsva import GSVA
        methods[EnrichmentMethod.GSVA.value] = GSVA(  # GSVA Gene Set Variation Analysis
            method=self.config.gsva_method,
            kcdf=self.config.gsva_kcdf,
            tau=self.config.gsva_tau,
            min_size=activity_min_size,
            max_size=activity_max_size
        )
        if self.config.method not in methods:
            raise ValueError(f"Unknown method: {self.config.method}")
        
        return methods[self.config.method]
    
    def load_gene_list(self, file_path: str) -> Set[str]:
        """Load gene identifiers from a text, CSV, TSV, or Excel file.
        
        Duplicate identifiers are removed while preserving first-occurrence order."""
        file_path = Path(file_path)
        
        if not file_path.exists():
            raise FileNotFoundError(f"Gene-list file does not exist: {file_path}")
        
        suffix = file_path.suffix.lower()
        
        genes = set()
        
        if suffix in ('.xlsx', '.xls'):
            genes = self._load_from_excel(file_path)
        elif suffix == '.csv':
            genes = self._load_from_csv(file_path)
        else:
            # Plain text and TSV inputs are detected from their contents.
            genes = self._load_from_text_auto(file_path)
        
        if not genes:
            raise ValueError(
                f"No valid gene identifiers were found in {file_path}. "
                "Use one identifier per row or a supported CSV/TSV layout."
            )
        
        logger.info("Loaded %d unique genes from %s", len(genes), file_path)
        return genes
    
    @staticmethod
    def load_ranked_gene_list(file_path: str) -> List[Tuple[str, float]]:
        """Load and validate a two-column ranked gene table for GSEA.
        
        The first column contains gene identifiers and the second contains numeric
        ranking statistics. Duplicate genes are rejected because ranks must be unique."""
        file_path_obj = Path(file_path)
        if not file_path_obj.exists():
            raise FileNotFoundError(f"Ranked gene table does not exist: {file_path}")
        
        df = pd.read_csv(file_path_obj, sep='\t')
        
        required_columns = {'gene', 'weight'}
        missing_columns = required_columns - set(df.columns)
        if missing_columns:
            raise ValueError(
                "The ranked gene table must contain 'gene' and 'weight' columns. "
                f"Available columns: {list(df.columns)}"
            )

        genes = df['gene'].astype('string').str.strip()
        weights = pd.to_numeric(df['weight'], errors='coerce')
        invalid_rows = genes.isna() | genes.eq('') | weights.isna() | ~np.isfinite(weights)
        if invalid_rows.any():
            rows = (invalid_rows[invalid_rows].index + 2).tolist()
            raise ValueError(
                "The ranked gene table contains empty gene IDs or non-finite weights "
                f"at rows: {rows}"
            )
        duplicated = genes[genes.duplicated(keep=False)].unique().tolist()
        if duplicated:
            raise ValueError(f"The ranked gene table contains duplicate gene IDs: {duplicated[:10]}")
        if weights.nunique() < 2:
            raise ValueError("The ranking statistic must contain at least two distinct values")

        ranked = pd.DataFrame({'gene': genes.astype(str), 'weight': weights.astype(float)})
        ranked = ranked.sort_values('weight', ascending=False, kind='mergesort')
        result = list(ranked.itertuples(index=False, name=None))
        
        logger.info("Loaded %d ranked genes from %s", len(result), file_path)
        return result
    
    def _load_from_excel(self, file_path: Path) -> Set[str]:
        """Load gene identifiers from the first column of an Excel worksheet."""
        genes = set()
        try:
            df = pd.read_excel(file_path, engine='openpyxl', header=None)
            for item in df.iloc[:, 0].dropna():
                gene = str(item).strip()
                if gene and not gene.startswith('#'):
                    genes.add(gene)
        except ImportError:
            raise ImportError(
                "Excel input requires openpyxl. Install it with: pip install openpyxl"
            )
        except Exception as e:
            raise ValueError(f"Could not read Excel file {file_path}: {e}")
        
        return genes
    
    def _load_from_csv(self, file_path: Path) -> Set[str]:
        """Load gene identifiers from the first column of a delimited text file."""
        genes = set()
        try:
            with open(file_path, 'r', encoding='utf-8') as f:
                reader = csv.reader(f)
                for row in reader:
                    if not row:
                        continue
                    gene = row[0].strip()
                    if gene and not gene.startswith('#'):
                        genes.add(gene)
        except Exception as e:
            raise ValueError(f"Could not read CSV file {file_path}: {e}")
        
        return genes
    
    def _load_from_text_auto(self, file_path: Path) -> Set[str]:
        """Detect a plain-text delimiter and load gene identifiers."""
        genes = set()
        
        with open(file_path, 'r', encoding='utf-8') as f:
            lines = f.readlines()
        
        # Ignore blank rows and comment lines.
        valid_lines = [line.strip() for line in lines if line.strip() and not line.strip().startswith('#')]
        
        if not valid_lines:
            return genes
        
        # A single row may contain several delimiter-separated identifiers.
        if len(valid_lines) == 1:
            line = valid_lines[0]
            # Prefer explicit delimiters before whitespace.
            for sep in [',', ';', '\t', ' ']:
                if sep in line:
                    parts = [p.strip() for p in line.split(sep) if p.strip()]
                    if len(parts) > 1:
                        genes.update(parts)
                        return genes
            genes.add(line)
            return genes
        
        # Detect a comma-delimited file and use its first column.
        first_line = valid_lines[0]
        if ',' in first_line:
            try:
                import io
                content = ''.join(lines)
                reader = csv.reader(io.StringIO(content))
                for row in reader:
                    if not row:
                        continue
                    gene = row[0].strip()
                    if gene and not gene.startswith('#'):
                        genes.add(gene)
                if len(genes) > 0:
                    return genes
            except Exception:
                pass
        
        # Otherwise parse one gene per row or the first column of a TSV.
        tab_count = sum(1 for line in valid_lines if '\t' in line)
        is_tsv = tab_count > 0 and tab_count >= len(valid_lines) * 0.5
        
        if is_tsv:
            # Skip a conventional header in the first TSV row.
            for i, line in enumerate(valid_lines):
                parts = line.split('\t')
                gene = parts[0].strip()
                if not gene:
                    continue
                if i == 0 and gene.lower() in ('gene', 'symbol', 'id', 'gene_id', 'gene_symbol', 'name', 'genes'):
                    continue
                genes.add(gene)
        else:
            for line in valid_lines:
                gene = line.strip()
                if gene and not gene.startswith('#'):
                    genes.add(gene)
        
        return genes
    
    def adjust_pvalues(
        self,
        results: List[EnrichmentResult],
        method: str = "BH"
    ) -> Any:
        """Apply the configured multiple-testing correction to valid P values."""
        if not results:
            return results
        
        # Extract raw P values from all terms.
        pvalues = [r.pvalue for r in results]
        
        # Activity methods may return NaN because they do not perform hypothesis tests.
        nan_count = sum(1 for p in pvalues if math.isnan(p))
        if nan_count > 0:
            logger.warning(
                "Detected %d/%d NaN p-values (possibly from ssGSEA); "
                "skipping multiple-testing correction",
                nan_count,
                len(pvalues),
            )
            return results
        
        # Map public correction names to statsmodels identifiers.
        correction_methods = {
            CorrectionMethod.BH.value: "fdr_bh",
            CorrectionMethod.BY.value: "fdr_by",
            CorrectionMethod.BONFERRONI.value: "bonferroni",
            CorrectionMethod.HOLM.value: "holm",
            CorrectionMethod.NONE.value: None,
        }
        
        # Preserve raw P values when correction is disabled.
        if method == CorrectionMethod.NONE.value:
            return results
        
        # Get the corresponding correction method identifier, use BH by default
        corr_method = correction_methods.get(method, "fdr_bh")
        
        # Call statsmodels' multi-test correction
        # Return value: reject (rejected or not, adjusted_pvalues, _, _
        _, adjusted, _, _ = multipletests(pvalues, method=corr_method)
        
        # Write back to each result object with the corrected p value
        for result, adj_p in zip(results, adjusted):
            result.adjusted_pvalue = adj_p
            # The FDR field of GSA should also be synchronized to the corrected value
            if result.fdr is not None:
                result.fdr = adj_p
        
        return results
    
    def filter_results(
        self,
        results: List[EnrichmentResult]
    ) -> List[EnrichmentResult]:
        """Filter result rows by the configured P- and Q-value thresholds."""
        filtered = []

        for result in results:
            if result.gene_count == 0:
                continue

            if not self.config.output_all:
                if (result.pvalue > self.config.pvalue_cutoff or
                        result.adjusted_pvalue > self.config.qvalue_cutoff):
                    continue

            set_size = result.set_size if result.set_size is not None else result.gene_count
            if set_size < self.config.min_genes:
                continue

            if self.config.max_genes != float('inf') and set_size > self.config.max_genes:
                continue

            filtered.append(result)

        return filtered
    
    def analyze_database(
        self,
        gene_set: Set[str],
        background_set: Set[str],
        term_data: Dict[str, Dict[str, Any]],
        database: str,
        ranked_gene_list: Optional[List[Tuple[str, float]]] = None
    ) -> List[EnrichmentResult]:
        """Analyze every term in one database with the selected method."""
        if self.config.method == EnrichmentMethod.GSEA.value:
            return self._run_gsea_with_fgsea(term_data, database, ranked_gene_list)
        if self.config.method in (EnrichmentMethod.SSGSEA.value, EnrichmentMethod.GSVA.value):
            raise ValueError(
                "ssGSEA and GSVA require an expression matrix; use "
                "analyze_activity_database() for activity analysis"
            )
        
        # ORA workflow
        results = []
        
        if not gene_set:
            logger.error("The query gene set is empty; check the input file and identifier type")
            return results
        
        if not background_set:
            logger.warning("The background gene set is empty; ORA cannot produce a valid result")
        
        if not term_data:
            logger.warning("Database %s contains no gene sets and will be skipped", database)
            return results
        
        # The M/N/n/k of ORA must come from the same statistical aggregate.
        universe = set(background_set)
        query_genes = set(gene_set) & universe
        background_total = len(universe)
        
        if background_total == 0:
            logger.warning("The resolved background is empty; skipping ORA for %s", database)
            return results
        
        gene_total = len(query_genes)
        if gene_total == 0:
            logger.warning("Query genes do not overlap the resolved background; skipping ORA for %s", database)
            return results

        is_tf_database = self._is_tf_database(database)
        if is_tf_database:
            filtered_sets = self._filter_tf_gene_sets(
                database,
                'ora',
                {
                    term_id: set(term_info.get('genes', []))
                    for term_id, term_info in term_data.items()
                },
                universe,
            )
            term_data = {
                term_id: {**term_data[term_id], 'genes': genes}
                for term_id, genes in filtered_sets.items()
            }

        # Apply size filters to the annotation universe before testing.
        term_stats = {}
        
        for term_id, term_info in term_data.items():
            term_genes = set(term_info.get("genes", [])) & universe
            num_in_C = len(term_genes)
            if not is_tf_database:
                if num_in_C < self.config.min_genes:
                    continue
                if self.config.max_genes != float('inf') and num_in_C > self.config.max_genes:
                    continue

            genes_in_term = query_genes & term_genes
            num_in_O = len(genes_in_term)
            
            term_stats[term_id] = {
                "num_in_O": num_in_O,
                "num_in_C": num_in_C,
                "genes_in_term": genes_in_term,
                "term_name": term_info.get("name", term_id),
            }
        
        # Test every eligible term so multiple-testing correction uses the full family.
        for term_id, ts in tqdm(
            term_stats.items(),
            desc=f"Analyzing {database}",
            leave=False
        ):
            try:
                num_in_O = ts["num_in_O"]
                num_in_C = ts["num_in_C"]
                
                pvalue = self.method.calculate_pvalue(
                    gene_count=num_in_O,
                    background_count=num_in_C,
                    gene_total=gene_total,
                    background_total=background_total
                )
                
                expected = num_in_C / background_total * gene_total if background_total > 0 else 0
                rich_factor = num_in_O / expected if expected > 0 else 0
                
                term_url = generate_term_url(term_id, database)
                
                result = EnrichmentResult(
                    term_id=term_id,
                    term_name=ts["term_name"],
                    database=database,
                    pvalue=pvalue,
                    adjusted_pvalue=pvalue,
                    gene_count=num_in_O,
                    background_count=num_in_C,
                    expected_count=round(expected, 6),
                    rich_factor=round(rich_factor, 6),
                    gene_list=sorted(ts["genes_in_term"]),
                    gene_ratio=f"{num_in_O}/{gene_total}",
                    background_ratio=f"{num_in_C}/{background_total}",
                    term_url=term_url,
                    set_size=num_in_C,
                )
                results.append(result)
            except Exception as e:
                logger.warning("Skipping term %s in %s after calculation error: %s", term_id, database, e)
                continue
        
        return results
    
    def _run_gsea_with_fgsea(
        self,
        term_data: Dict[str, Dict[str, Any]],
        database: str,
        ranked_gene_list: Optional[List[Tuple[str, float]]],
    ) -> pd.DataFrame:
        """Run Bioconductor fgsea and preserve its official result columns."""
        from allenricher.core.bioconductor import FGSEA_COLUMNS, run_fgsea

        gene_sets_dict: Dict[str, Set[str]] = {}
        for term_id, term_info in term_data.items():
            term_genes = set(term_info.get("genes", []))
            if term_genes:
                gene_sets_dict[term_id] = term_genes

        if not gene_sets_dict:
            logger.warning("Database %s contains no valid gene sets and will be skipped", database)
            return pd.DataFrame(columns=FGSEA_COLUMNS)
        if not ranked_gene_list:
            raise ValueError("GSEA requires a ranked gene table with 'gene' and numeric 'weight' columns")

        backend_min_size = getattr(self.method, "min_size", 15)
        backend_max_size = getattr(self.method, "max_size", 500)
        if self._is_tf_database(database):
            gene_sets_dict = self._filter_tf_gene_sets(
                database,
                'gsea',
                gene_sets_dict,
                {str(gene).strip() for gene, _ in ranked_gene_list},
            )
            backend_min_size = 1
            backend_max_size = max(len(ranked_gene_list), 1)

        logger.info(
            "Running Bioconductor fgseaMultilevel with %d gene sets and %d ranked genes",
            len(gene_sets_dict),
            len(ranked_gene_list),
        )
        return run_fgsea(
            ranked_gene_list,
            gene_sets_dict,
            min_size=backend_min_size,
            max_size=backend_max_size,
            seed=getattr(self.method, "seed", 42),
        )

    def analyze_activity_database(
        self,
        expression_matrix: pd.DataFrame,
        gene_sets: Dict[str, Set[str]],
        database: str,
    ) -> pd.DataFrame:
        """Run ssGSEA or GSVA and preserve the official activity matrix."""
        if not self._is_tf_database(database):
            return self.method.analyze_matrix(expression_matrix, gene_sets)

        from allenricher.core.bioconductor import run_gsva

        filtered = self._filter_tf_gene_sets(
            database,
            self.config.method,
            gene_sets,
            set(map(str, expression_matrix.index)),
        )
        method = 'ssgsea' if self.config.method == 'ssgsea' else self.config.gsva_method
        return run_gsva(
            expression_matrix,
            filtered,
            method=method,
            kcdf=self.config.gsva_kcdf,
            tau=0.25 if method == 'ssgsea' else self.config.gsva_tau,
            min_size=1,
            max_size=max(len(expression_matrix.index), 1),
        )
    
    def _finalize_database_results(self, results: Any) -> Any:
        """Attach term metadata, filter results, and store one database table."""
        if isinstance(results, pd.DataFrame):
            if results.empty or self.config.output_all:
                return results
            return results.loc[
                (results["pval"] <= self.config.pvalue_cutoff)
                & (results["padj"] <= self.config.qvalue_cutoff)
            ].reset_index(drop=True)

        results = self.adjust_pvalues(results, self.config.correction)
        results = self.filter_results(results)
        results.sort(key=lambda item: item.adjusted_pvalue)
        return results

    def run_analysis(
        self,
        gene_set: Set[str],
        background_set: Set[str],
        database_data: Dict[str, Dict[str, Dict[str, Any]]],
        parallel: bool = True,
        ranked_gene_list: Optional[List[Tuple[str, float]]] = None
    ) -> Dict[str, pd.DataFrame]:
        """Run the configured method across all loaded databases."""
        self.results = {}

        if self.config.method in (EnrichmentMethod.SSGSEA.value, EnrichmentMethod.GSVA.value):
            raise ValueError(
                "ssGSEA and GSVA require an expression matrix; use "
                "analyze_activity_database() for activity analysis"
            )
        
        if parallel and self.config.n_jobs > 1 and len(database_data) > 1:
            logger.info(
                "Analyzing %d databases with %d worker threads",
                len(database_data),
                self.config.n_jobs,
            )
            with ThreadPoolExecutor(max_workers=self.config.n_jobs) as executor:
                futures = {
                    executor.submit(
                        self.analyze_database,
                        gene_set,
                        background_set,
                        term_data,
                        database,
                        ranked_gene_list
                    ): database
                    for database, term_data in database_data.items()
                }
                
                # Collect each database result as its worker completes.
                for future in as_completed(futures):
                    database = futures[future]
                    try:
                        results = self._finalize_database_results(future.result())
                        self.results[database] = results
                        logger.info(f"Completed {database}: {len(results)} enriched terms")
                    except Exception as e:
                        if self.config.method == EnrichmentMethod.GSEA.value:
                            raise
                        logger.error(f"Error analyzing {database}: {e}")
        else:
            # Analyze databases sequentially when parallelism is disabled.
            for database, term_data in database_data.items():
                try:
                    results = self.analyze_database(
                        gene_set, background_set, term_data, database,
                        ranked_gene_list=ranked_gene_list
                    )
                    results = self._finalize_database_results(results)
                    self.results[database] = results
                    logger.info(f"Completed {database}: {len(results)} enriched terms")
                except Exception as e:
                    if self.config.method == EnrichmentMethod.GSEA.value:
                        raise
                    logger.error(f"Error analyzing {database}: {e}")
        
        dataframes = self.to_dataframes()
        for database, frame in dataframes.items():
            decorated = add_result_term_metadata(frame, database_data.get(database, {}))
            dataframes[database] = decorated
            self.results[database] = decorated
        return dataframes
    
    def to_dataframes(self) -> Dict[str, pd.DataFrame]:
        """Return result tables keyed by database name."""
        dataframes = {}
        
        for database, results in self.results.items():
            if isinstance(results, pd.DataFrame):
                if not results.empty:
                    dataframes[database] = results.copy()
            elif results:
                # Convert ORA result objects to the public tabular schema.
                data = [r.to_dict() for r in results]
                df = pd.DataFrame(data)
                dataframes[database] = df
        
        return dataframes
    
    def save_results(self, output_dir: str, metadata: dict = None) -> None:
        """Write each database result table as TSV and optionally CSV."""
        output_path = Path(output_dir)
        output_path.mkdir(parents=True, exist_ok=True)
        
        for database, results in self.results.items():
            if isinstance(results, pd.DataFrame):
                if results.empty:
                    continue
                df = results.copy()
            elif results:
                df = pd.DataFrame([result.to_dict() for result in results])
            else:
                continue

            if not df.empty:
                output_file = output_path / f"{database}_enrichment.tsv"

                df.to_csv(output_file, sep='\t', index=False, lineterminator='\n')

                logger.info(f"Saved {database} results to {output_file}")

    def get_annotated_genes(self, gene_set_data: Dict[str, Dict[str, Any]]) -> Set[str]:
        """Return the union of genes annotated by loaded databases."""
        annotated = set()
        for term_data in gene_set_data.values():
            genes = term_data.get("genes", set())
            if isinstance(genes, set):
                annotated.update(genes)
            elif isinstance(genes, (list, tuple)):
                annotated.update(genes)
        return annotated

    def resolve_background(self, 
                           gene_set_data: Dict[str, Dict[str, Any]],
                           user_background: Optional[Set[str]] = None,
                           background_mode: str = "annotated") -> Set[str]:
        """Resolve the ORA background from custom, genome, or annotated genes."""
        if background_mode == "custom":
            if not user_background:
                raise ValueError("background_mode='custom' requires user_background to be provided")
            return user_background
        elif background_mode == "annotated":
            return self.get_annotated_genes(gene_set_data)
        elif background_mode == "genome":
            if not user_background:
                logger.warning(
                    "background_mode='genome' was requested without a genome background; "
                    "falling back to annotated genes"
                )
                return self.get_annotated_genes(gene_set_data)
            return user_background
        else:
            raise ValueError(f"Unknown background_mode: {background_mode}. Expected 'annotated', 'genome', or 'custom'")
