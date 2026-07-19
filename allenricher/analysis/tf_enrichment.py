"""Transcription-factor enrichment workflows.

TF-target gene sets from TRRUST, ChEA3, AnimalTFDB, and hTFtarget are analyzed
with ORA, GSEA, ssGSEA, or GSVA using method-specific size limits."""

import logging
from typing import Dict, List, Optional, Set, Tuple

import numpy as np
import pandas as pd
from scipy.stats import hypergeom
from statsmodels.stats.multitest import multipletests

logger = logging.getLogger(__name__)

TF_SIZE_DEFAULTS = {
    'ora': (3, None),
    'gsea': (15, 5000),
    'ssgsea': (1, None),
    'gsva': (1, None),
}


def resolve_tf_size_limits(
    method: str,
    min_size: Optional[int] = None,
    max_size: Optional[int] = None,
) -> Tuple[int, Optional[int]]:
    """Resolve method-specific TF target-set size limits.

    ``None`` represents an unbounded maximum."""
    method = 'ora' if method == 'hypergeometric' else method.lower()
    if method not in TF_SIZE_DEFAULTS:
        raise ValueError(f"Unsupported TF analysis: {method}")
    default_min, default_max = TF_SIZE_DEFAULTS[method]
    resolved_min = default_min if min_size is None else min_size
    resolved_max = default_max if max_size is None else max_size
    if resolved_min < 1:
        raise ValueError("TF min_size must be greater than zero")
    if resolved_max is not None and resolved_min > resolved_max:
        raise ValueError("TF min_size cannot exceed max_size")
    return resolved_min, resolved_max


def filter_tf_gene_sets(
    gene_sets: Dict[str, Set[str]],
    universe: Set[str],
    min_size: int,
    max_size: Optional[int],
) -> Tuple[Dict[str, Set[str]], Dict[str, int]]:
    """Filter TF target sets after intersection with the analysis universe."""
    filtered: Dict[str, Set[str]] = {}
    stats = {
        'before': len(gene_sets),
        'after': 0,
        'no_overlap': 0,
        'below_min': 0,
        'above_max': 0,
    }
    for term_id, targets in gene_sets.items():
        overlap = set(targets) & universe
        if not overlap:
            stats['no_overlap'] += 1
        elif len(overlap) < min_size:
            stats['below_min'] += 1
        elif max_size is not None and len(overlap) > max_size:
            stats['above_max'] += 1
        else:
            filtered[term_id] = overlap
    stats['after'] = len(filtered)
    return filtered, stats

class TFEnrichmentAnalyzer:
    """Analyze TF-target gene sets against query, ranked, or expression data.

    The analyzer accepts the standardized database object returned by
    ``DatabaseManager`` and preserves TF library, context, evidence, and regulatory-
    mode metadata in downstream tables."""

    def __init__(self, tf_database: Dict, background_size: Optional[int] = None):
        """Build TF-to-target mappings from a standardized TF database."""
        if 'gene2tf' not in tf_database:
            raise ValueError(
                "tf_database must contain a 'gene2tf' matrix. Load the database "
                "through DatabaseManager before constructing the analyzer."
            )

        gene2tf_df = tf_database['gene2tf'].copy()
        self.tf_info = tf_database.get('tf_info', pd.DataFrame())
        self.edges = tf_database.get('edges', pd.DataFrame())
        # Convert the gene-by-TF incidence matrix into TF-to-target sets.
        tf_columns = [col for col in gene2tf_df.columns if col != 'Gene']

        # Normalize serialized numeric cells before testing incidence values.
        for col in tf_columns:
            gene2tf_df[col] = pd.to_numeric(gene2tf_df[col], errors='coerce').fillna(0).astype(int)

        self.tf_to_targets: Dict[str, Set[str]] = {}

        for tf in tf_columns:
            mask = gene2tf_df[tf] == 1
            targets = set(
                gene2tf_df.loc[mask, 'Gene'].dropna().astype(str)
                if 'Gene' in gene2tf_df.columns
                else gene2tf_df.index[mask].astype(str)
            )
            targets = {target.strip() for target in targets if target.strip()}
            if targets:
                self.tf_to_targets[tf] = targets

        # The default universe contains genes participating in at least one
        # TF-target relation; zero-only matrix rows must not expand it.
        self.background_genes = set().union(*self.tf_to_targets.values()) if self.tf_to_targets else set()
        if not self.background_genes:
            raise ValueError("The gene2tf matrix does not contain any TF-target relationships")
        if background_size is not None and background_size != len(self.background_genes):
            raise ValueError(
                "background_size does not define the TF ORA universe. Provide the "
                "actual detected genes through ora(background_genes=...)."
            )
        self.background_size = len(self.background_genes)

        self.term_metadata: Dict[str, Dict[str, str]] = {}
        if not self.tf_info.empty:
            id_column = 'Term_ID' if 'Term_ID' in self.tf_info.columns else 'TF'
            if id_column in self.tf_info.columns:
                for _, row in self.tf_info.iterrows():
                    term_id = str(row[id_column])
                    self.term_metadata[term_id] = {
                        str(column): str(value) if pd.notna(value) else ''
                        for column, value in row.items()
                    }

        logger.info(
            "Initialized TF enrichment analyzer with %d TFs and %d background genes",
            len(self.tf_to_targets), self.background_size
        )

    def _record_size_filter(
        self,
        method: str,
        stats: Dict[str, int],
        min_size: int,
        max_size: Optional[int],
    ) -> None:
        logger.info(
            "TF %s gene-set size filter: min=%d, max=%s, before=%d, after=%d, "
            "no_overlap=%d, below_min=%d, above_max=%d",
            method,
            min_size,
            'Inf' if max_size is None else max_size,
            stats['before'],
            stats['after'],
            stats['no_overlap'],
            stats['below_min'],
            stats['above_max'],
        )

    @staticmethod
    def _filter_values(value) -> Set[str]:
        if value is None:
            return set()
        values = value if isinstance(value, (list, tuple, set)) else str(value).split(',')
        return {str(item).strip().lower() for item in values if str(item).strip()}

    def _metadata(self, term_id: str) -> Dict[str, str]:
        item = dict(self.term_metadata.get(term_id, {}))
        item.setdefault('Term_ID', term_id)
        item.setdefault('Term_Name', term_id)
        item.setdefault('TF', term_id)
        item.setdefault('Library', item.get('Source', 'unknown'))
        item.setdefault('Context', '')
        item.setdefault('Evidence_Type', '')
        item.setdefault('Inference_Type', 'direct')
        return item

    def metadata_frame(self) -> pd.DataFrame:
        return pd.DataFrame([self._metadata(term_id) for term_id in self.tf_to_targets])

    def _scoped_gene_sets(
        self,
        library=None,
        tissue=None,
        regulation: str = 'all',
    ) -> Dict[str, Set[str]]:
        library_filter = self._filter_values(library)
        tissue_filter = self._filter_values(tissue)
        regulation = (regulation or 'all').lower()
        if regulation not in {'all', 'activation', 'repression', 'unknown'}:
            raise ValueError(f"Unsupported TF regulatory mode: {regulation}")

        regulation_targets: Dict[str, Set[str]] = {}
        if regulation != 'all':
            if self.edges.empty or not {'TF', 'Target', 'Mode'}.issubset(self.edges.columns):
                raise ValueError(
                    "The current TRRUST snapshot does not contain edge-level Mode metadata; "
                    "rebuild the local TRRUST database"
                )
            wanted = regulation.lower()
            for tf, target, mode in self.edges[['TF', 'Target', 'Mode']].itertuples(index=False, name=None):
                if str(mode).strip().lower() == wanted:
                    regulation_targets.setdefault(str(tf), set()).add(str(target))

        scoped: Dict[str, Set[str]] = {}
        for term_id, original_targets in self.tf_to_targets.items():
            item = self._metadata(term_id)
            tf = item['TF']
            if library_filter and item['Library'].lower() not in library_filter:
                continue
            if tissue_filter and item['Context'].lower() not in tissue_filter:
                continue
            targets = regulation_targets.get(tf, set()) if regulation != 'all' else original_targets
            targets = set(targets) & self.background_genes
            if targets:
                scoped[term_id] = targets
        return scoped

    def _select_gene_sets(
        self,
        tf_list: Optional[List[str]] = None,
        library=None,
        tissue=None,
        regulation: str = 'all',
        min_size: int = 1,
        max_size: Optional[int] = None,
    ) -> Dict[str, Set[str]]:
        tf_filter = self._filter_values(tf_list)
        selected: Dict[str, Set[str]] = {}
        for term_id, targets in self._scoped_gene_sets(library, tissue, regulation).items():
            item = self._metadata(term_id)
            if (
                tf_filter
                and term_id.lower() not in tf_filter
                and item['TF'].lower() not in tf_filter
            ):
                continue
            if len(targets) < min_size or (max_size is not None and len(targets) > max_size):
                continue
            selected[term_id] = targets
        return selected

    def ora(
        self,
        gene_set: List[str],
        tf_list: Optional[List[str]] = None,
        min_overlap: int = 3,
        library=None,
        tissue=None,
        regulation: str = 'all',
        min_size: int = 3,
        max_size: Optional[int] = None,
        background_genes: Optional[Set[str]] = None,
    ) -> pd.DataFrame:
        """Run one-sided hypergeometric ORA for TF target sets.

        Size filters and the background are applied independently within each source
        library. Benjamini-Hochberg correction therefore uses the complete tested term
        family for that library."""
        regulation = (regulation or 'all').lower()
        input_gene_set = {
            str(gene).strip() for gene in gene_set
            if gene is not None and str(gene).strip()
        }
        if not input_gene_set:
            raise ValueError("TF ORA requires a non-empty query gene list")

        scoped_gene_sets = self._scoped_gene_sets(
            library=library,
            tissue=tissue,
            regulation=regulation,
        )
        candidate_gene_sets = self._select_gene_sets(
            tf_list=tf_list,
            library=library,
            tissue=tissue,
            regulation=regulation,
            min_size=1,
            max_size=None,
        )

        custom_background = None
        if background_genes is not None:
            custom_background = {
                str(gene).strip() for gene in background_genes
                if gene is not None and str(gene).strip()
            }
            if not custom_background:
                raise ValueError("The custom TF ORA background cannot be empty")

        library_backgrounds: Dict[str, Set[str]] = {}
        for term_id, targets in scoped_gene_sets.items():
            source = self._metadata(term_id)['Library']
            if custom_background is not None:
                library_backgrounds[source] = custom_background
            else:
                library_backgrounds.setdefault(source, set()).update(targets)

        gene_sets: Dict[str, Set[str]] = {}
        size_stats = {'before': 0, 'after': 0, 'no_overlap': 0, 'below_min': 0, 'above_max': 0}
        for source, universe in library_backgrounds.items():
            source_sets = {
                term_id: targets for term_id, targets in candidate_gene_sets.items()
                if self._metadata(term_id)['Library'] == source
            }
            filtered, source_stats = filter_tf_gene_sets(
                source_sets, universe, min_size, max_size
            )
            gene_sets.update(filtered)
            for key in size_stats:
                size_stats[key] += source_stats[key]
        self._record_size_filter('ora', size_stats, min_size, max_size)

        for source, universe in sorted(library_backgrounds.items()):
            mapped = len(input_gene_set & universe)
            logger.info(
                "TF ORA background [%s]: %s; universe=%d, input=%d, mapped=%d (%.1f%%)",
                source,
                "custom detected-gene background" if custom_background is not None else "library annotation background",
                len(universe),
                len(input_gene_set),
                mapped,
                100.0 * mapped / len(input_gene_set),
            )

        results = []
        for term_id, targets in gene_sets.items():
            item = self._metadata(term_id)
            universe = library_backgrounds[item['Library']]
            gene_set_set = input_gene_set & universe
            n_gene_set = len(gene_set_set)
            overlap = gene_set_set & targets
            n_overlap = len(overlap)

            n_targets = len(targets)

            # One-sided hypergeometric probability P(X >= observed overlap).
            pvalue = hypergeom.sf(
                n_overlap - 1,
                len(universe),
                n_targets,
                n_gene_set
            )

            mode = (
                regulation.lower()
                if regulation != 'all'
                else item.get('Mode', item.get('mode', 'unknown')).lower()
            )
            a = n_overlap
            b = n_gene_set - n_overlap
            c = n_targets - n_overlap
            d = len(universe) - a - b - c
            odds_ratio = (a * d) / (b * c) if b > 0 and c > 0 else float('inf')

            results.append({
                'Term_ID': term_id,
                'Term_Name': item['Term_Name'],
                'TF': item['TF'],
                'Library': item['Library'],
                'Context': item['Context'],
                'Evidence_Type': item['Evidence_Type'],
                'Inference_Type': item['Inference_Type'],
                'Overlap': n_overlap,
                'TF_Targets': n_targets,
                'Target_Set_Size': n_targets,
                'Input_GeneSet_Size': len(input_gene_set),
                'GeneSet_Size': n_gene_set,
                'Background_Size': len(universe),
                'Mapping_Rate': n_gene_set / len(input_gene_set),
                'Overlap_Genes': ','.join(sorted(overlap)),
                'Pvalue': pvalue,
                'Odds_Ratio': odds_ratio,
                'Mode': mode,
            })

        if not results:
            logger.warning("TF ORA produced no tested TF terms (min_overlap=%d)", min_overlap)
            return pd.DataFrame(
                columns=['Term_ID', 'Term_Name', 'TF', 'Library', 'Context',
                         'Evidence_Type', 'Inference_Type', 'Overlap',
                         'TF_Targets', 'Target_Set_Size', 'Input_GeneSet_Size',
                         'GeneSet_Size', 'Background_Size', 'Mapping_Rate',
                         'Overlap_Genes', 'Pvalue', 'FDR', 'Odds_Ratio', 'Mode']
            )

        result_df = pd.DataFrame(results)

        # Match Enrichr/ChEA3 practice by correcting each evidence library independently.
        result_df['FDR'] = result_df.groupby('Library', dropna=False)['Pvalue'].transform(
            lambda values: multipletests(values.to_numpy(), method='fdr_bh')[1]
            if len(values) > 1 else values.to_numpy()
        )

        # min_overlap controls reported rows, not the multiple-testing universe.
        # Every source term that passed the source/context/size filters must
        # contribute to its library's BH denominator, including zero-overlap terms.
        result_df = result_df[result_df['Overlap'] >= min_overlap].copy()
        if result_df.empty:
            logger.warning(
                "TF ORA found no terms with at least %d overlapping genes", min_overlap
            )
            return pd.DataFrame(
                columns=['Term_ID', 'Term_Name', 'TF', 'Library', 'Context',
                         'Evidence_Type', 'Inference_Type', 'Overlap',
                         'TF_Targets', 'Target_Set_Size', 'Input_GeneSet_Size',
                         'GeneSet_Size', 'Background_Size', 'Mapping_Rate',
                         'Overlap_Genes', 'Pvalue', 'FDR', 'Odds_Ratio', 'Mode']
            )

        # Present the most significant TFs first.
        result_df = result_df.sort_values('Pvalue').reset_index(drop=True)

        logger.info("TF ORA completed with %d reported TF terms", len(result_df))
        return result_df

    def gsea(
        self,
        ranked_genes: List[Tuple[str, float]],
        tf_list: Optional[List[str]] = None,
        n_permutations: int = 1000,
        seed: int = 42,
        library=None,
        tissue=None,
        regulation: str = 'all',
        min_size: int = 15,
        max_size: int = 5000,
    ) -> pd.DataFrame:
        """Run preranked GSEA for selected TF target sets."""
        from allenricher.core.bioconductor import run_fgsea

        gene_sets = self._select_gene_sets(
            tf_list=tf_list,
            library=library,
            tissue=tissue,
            regulation=regulation,
            min_size=1,
            max_size=None,
        )
        gene_sets, stats = filter_tf_gene_sets(
            gene_sets, {str(gene).strip() for gene, _ in ranked_genes}, min_size, max_size
        )
        self._record_size_filter('gsea', stats, min_size, max_size)
        result = run_fgsea(ranked_genes, gene_sets, min_size=1,
                           max_size=max(len(ranked_genes), 1), seed=seed)
        logger.info("TF GSEA analysis completed: %d TFs", len(result))
        return result

    def ssgsea(
        self,
        expression_df: pd.DataFrame,
        tf_list: Optional[List[str]] = None,
        min_size: int = 1,
        max_size: Optional[int] = None,
        library=None,
        tissue=None,
        regulation: str = 'all',
    ) -> pd.DataFrame:
        """Calculate sample-level ssGSEA scores for TF target sets."""
        from allenricher.core.bioconductor import run_gsva

        gene_sets = self._select_gene_sets(
            tf_list=tf_list,
            library=library,
            tissue=tissue,
            regulation=regulation,
            min_size=1,
            max_size=None,
        )
        gene_sets, stats = filter_tf_gene_sets(
            gene_sets, set(map(str, expression_df.index)), min_size, max_size
        )
        self._record_size_filter('ssgsea', stats, min_size, max_size)
        result_df = run_gsva(
            expression_df,
            gene_sets,
            method="ssgsea",
            tau=0.25,
            min_size=1,
            max_size=max(len(expression_df.index), 1),
        )
        result_df.index.name = "TF"
        logger.info(
            "TF ssGSEA analysis completed: %d TFs x %d samples",
            len(result_df),
            len(result_df.columns),
        )
        return result_df

    def _calculate_ssgsea_score(
        self,
        ranked_genes: List[str],
        gene_set: Set[str],
        gene_weights: Optional[Dict[str, float]] = None
    ) -> Tuple[float, float, float]:
        """Calculate one ssGSEA score from a sample-specific ranked gene list."""
        n = len(ranked_genes)
        nh = len(gene_set & set(ranked_genes))

        if nh == 0:
            return 0.0, 0.0, 0.0

        # Calculate the sum of weights
        if gene_weights:
            nr = sum(abs(gene_weights.get(g, 1.0)) for g in gene_set if g in gene_weights)
        else:
            nr = nh

        hit_inc = 1.0 / nr if nr > 0 else 0
        miss_inc = 1.0 / (n - nh) if (n - nh) > 0 else 0

        running_sum = 0.0
        max_es = 0.0
        min_es = 0.0

        for gene in ranked_genes:
            if gene in gene_set:
                weight = abs(gene_weights.get(gene, 1.0)) if gene_weights else 1.0
                running_sum += hit_inc * weight
                if running_sum > max_es:
                    max_es = running_sum
            else:
                running_sum -= miss_inc
                if running_sum < min_es:
                    min_es = running_sum

        return max_es, min_es, max_es

    def gsva(
        self,
        expression_df: pd.DataFrame,
        tf_list: Optional[List[str]] = None,
        method: str = "gsva",
        kcdf: str = "Gaussian",
        tau: float = 1.0,
        min_size: int = 1,
        max_size: Optional[int] = None,
        library=None,
        tissue=None,
        regulation: str = 'all',
    ) -> pd.DataFrame:
        """Calculate GSVA activity scores for selected TF target sets."""
        if expression_df.empty:
            logger.warning("The expression matrix is empty; returning an empty result")
            return pd.DataFrame()

        selected = self._select_gene_sets(
            tf_list=tf_list,
            library=library,
            tissue=tissue,
            regulation=regulation,
            min_size=1,
            max_size=None,
        )

        gene_sets, stats = filter_tf_gene_sets(
            selected, set(map(str, expression_df.index)), min_size, max_size
        )
        self._record_size_filter('gsva', stats, min_size, max_size)

        if not gene_sets:
            logger.warning("No TF target sets passed the GSVA size filter")
            return pd.DataFrame()

        from allenricher.core.bioconductor import run_gsva

        result_df = run_gsva(
            expression_df,
            gene_sets,
            method=method,
            kcdf=kcdf,
            tau=tau,
            min_size=1,
            max_size=max(len(expression_df.index), 1),
        )

        # Name the result index explicitly.
        if not result_df.empty:
            result_df.index.name = "TF"

        logger.info(
            "TF GSVA analysis completed: %d TFs x %d samples",
            len(gene_sets),
            len(expression_df.columns),
        )
        return result_df

    def get_activators(self, result_df: pd.DataFrame) -> pd.DataFrame:
        """Return significantly enriched activating TFs."""
        if result_df.empty or 'Mode' not in result_df.columns:
            return pd.DataFrame()

        return result_df[result_df['Mode'] == 'activator'].reset_index(drop=True)

    def get_repressors(self, result_df: pd.DataFrame) -> pd.DataFrame:
        """Return significantly enriched repressing TFs."""
        if result_df.empty or 'Mode' not in result_df.columns:
            return pd.DataFrame()

        return result_df[result_df['Mode'] == 'repressor'].reset_index(drop=True)

    def _build_tf_mode_map(self) -> Dict[str, str]:
        """Return regulatory modes keyed by stable TF term ID."""
        mode_map: Dict[str, str] = {}

        if self.tf_info.empty:
            return mode_map

        # TTRUST format: Term_ID/TF, Mode, Target_Set_Size_Sizep
        mode_col = 'Mode' if 'Mode' in self.tf_info.columns else 'mode'
        if mode_col in self.tf_info.columns:
            tf_col = 'Term_ID' if 'Term_ID' in self.tf_info.columns else 'TF'
            for _, row in self.tf_info.iterrows():
                tf_name = str(row[tf_col])
                mode = str(row[mode_col]).lower()
                mode_map[tf_name] = mode

        # ChEA3 format: TF, lib_country, target_country (no mode column)
        # Default for all TFs to 'unknown'

        return mode_map

    @staticmethod
    def _calculate_es(
        ranked_genes: List[str],
        gene_set: Set[str],
        n_genes: int,
        nh: int
    ) -> float:
        """Calculate the running enrichment score for one TF target set."""
        if nh == 0 or n_genes == nh:
            return 0.0

        hit_inc = 1.0 / nh
        miss_inc = 1.0 / (n_genes - nh)

        running_sum = 0.0
        max_es = 0.0

        for gene in ranked_genes:
            if gene in gene_set:
                running_sum += hit_inc
                if running_sum > max_es:
                    max_es = running_sum
            else:
                running_sum -= miss_inc

        return max_es
