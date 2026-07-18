"""Gene Set Variation Analysis support for AllEnricher.

The public production workflow delegates GSVA, ssGSEA, PLAGE, and z-score
activity calculations to Bioconductor GSVA. The Python implementation retained
here supports compatibility and focused unit testing."""

import logging
from abc import ABC, abstractmethod
from dataclasses import dataclass
from typing import TYPE_CHECKING, Dict, List, Optional, Set, Tuple, Any

import numpy as np
import pandas as pd
from scipy.stats import norm

logger = logging.getLogger(__name__)

# Import types only during static checking to avoid a runtime import cycle.
if TYPE_CHECKING:
    from allenricher.core.enrichment import EnrichmentMethodBase, EnrichmentResult


class GSVA(ABC):
    """Calculate sample-level gene-set activity scores.
    
    Supported method variants are ``gsva``, ``plage``, and ``zscore``. Gene-set
    size limits are applied after identifiers are matched to the expression matrix."""

    def __init__(
        self,
        method: str = "gsva",
        kcdf: str = "Gaussian",
        tau: float = 1.0,
        min_size: int = 1,
        max_size: Optional[int] = None
    ):
        """Configure the GSVA method, kernel, weighting parameter, and size limits."""
        valid_methods = ("gsva", "plage", "zscore")
        if method not in valid_methods:
            raise ValueError(
                f"Unsupported GSVA method: '{method}'. Valid values: {valid_methods}"
            )
        valid_kcdf = ("Gaussian", "Poisson")
        if kcdf not in valid_kcdf:
            raise ValueError(
                f"Unsupported GSVA kernel: '{kcdf}'. Valid values: {valid_kcdf}"
            )

        self.method = method
        self.kcdf = kcdf
        self.tau = tau
        self.min_size = min_size
        self.max_size = max_size

    def calculate_pvalue(
        self,
        gene_count: int,
        background_count: int,
        gene_total: int,
        background_total: int
    ) -> float:
        """Return NaN because GSVA activity scores are not hypothesis-test P values."""
        return float('nan')

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
    ) -> Optional["EnrichmentResult"]:
        """Return a compatibility result for one gene set.
        
        Matrix-level activity analysis should use :meth:`analyze_matrix`."""
        # Import lazily to avoid a module-level cycle.
        from allenricher.core.enrichment import EnrichmentResult, generate_term_url

        # Apply size limits after matching identifiers.
        overlap = gene_set & term_genes
        if len(overlap) < self.min_size or (
            self.max_size is not None and len(overlap) > self.max_size
        ):
            return None

        term_url = generate_term_url(term_id, database)

        # This compatibility object is not the primary matrix-level GSVA output.
        return EnrichmentResult(
            term_id=term_id,
            term_name=term_name,
            database=database,
            pvalue=float('nan'),
            adjusted_pvalue=float('nan'),
            gene_count=len(overlap),
            background_count=len(term_genes),
            expected_count=0,
            rich_factor=0,
            gene_list=list(overlap),
            gene_ratio=f"{len(overlap)}/{len(gene_set)}",
            background_ratio=f"{len(term_genes)}/{len(background_set)}",
            term_url=term_url,
            nes=float('nan'),
            es=float('nan'),
            fdr=float('nan'),
            leading_edge=None
        )

    def _compute_ecdf_kde(
        self,
        expression_matrix: np.ndarray,
        kcdf: str = "Gaussian",
        tau: float = 1.0
    ) -> Tuple[np.ndarray, np.ndarray]:
        """Estimate lower- and upper-tail expression distributions for each gene.
        
        Args:
            expression_matrix: Numeric array with genes in rows and samples in columns.
            kcdf: ``Gaussian`` for continuous values or ``Poisson`` for counts.
            tau: Gaussian bandwidth multiplier.
        
        Returns:
            Lower- and upper-tail score matrices with the input shape."""
        n_genes, n_samples = expression_matrix.shape

        # Allocate lower- and upper-tail score matrices.
        ecdf_down = np.zeros_like(expression_matrix, dtype=float)
        ecdf_up = np.zeros_like(expression_matrix, dtype=float)

        for i in range(n_genes):
            gene_expr = expression_matrix[i, :]

            if kcdf == "Gaussian":
                # Gaussian smoothing for continuous expression values.
                std = np.std(gene_expr, ddof=1 if n_samples > 1 else 0)
                if not np.isfinite(std) or std == 0:
                    # Constant expression is neutral across samples.
                    ecdf_down[i, :] = 0.5
                    ecdf_up[i, :] = 0.5
                    continue

                # Standardize this gene across samples.
                mean = np.mean(gene_expr)
                z_scores = (gene_expr - mean) / std

                for j in range(n_samples):
                    # Estimate each sample relative to all observations for this gene.
                    z_j = z_scores[j]

                    # Gaussian weights use tau as the bandwidth multiplier.
                    kernel_weights = np.exp(-((z_scores - z_j) ** 2) / (2 * tau ** 2))
                    kernel_weights = kernel_weights / np.sum(kernel_weights)

                    # Lower and upper cumulative probabilities around the current value.
                    ecdf_down[i, j] = np.sum(kernel_weights * (z_scores <= z_j))
                    ecdf_up[i, j] = np.sum(kernel_weights * (z_scores > z_j))

            elif kcdf == "Poisson":
                # Empirical cumulative probabilities for count data.
                for j in range(n_samples):
                    ecdf_down[i, j] = np.sum(gene_expr <= gene_expr[j]) / n_samples
                    ecdf_up[i, j] = np.sum(gene_expr > gene_expr[j]) / n_samples

        return ecdf_down, ecdf_up

    def _compute_gsva_score(
        self,
        gene_set_genes: List[str],
        gene_to_idx: Dict[str, int],
        ecdf_down: np.ndarray,
        ecdf_up: np.ndarray,
        n_genes_total: int
    ) -> np.ndarray:
        """Calculate standard GSVA random-walk scores for one gene set."""
        n_samples = ecdf_down.shape[1]

        # Locate gene-set members represented in the expression matrix.
        gene_indices = [gene_to_idx[g] for g in gene_set_genes if g in gene_to_idx]
        gene_index_set = set(gene_indices)

        if len(gene_indices) == 0:
            return np.zeros(n_samples)

        scores = np.zeros(n_samples)

        # Count genes inside and outside the set.
        n_genes_in_set = len(gene_index_set)
        n_genes_not_in_set = n_genes_total - n_genes_in_set

        if n_genes_in_set == 0 or n_genes_not_in_set == 0:
            return np.zeros(n_samples)

        for j in range(n_samples):
            # Rank genes by their relative expression distribution in this sample.
            ecdf_diff = ecdf_up[:, j] - ecdf_down[:, j]

            sorted_indices = np.argsort(-ecdf_diff)

            running_sum = 0.0
            max_deviation = 0.0

            # Traverse genes from highest to lowest relative expression.
            for idx in sorted_indices:
                if idx in gene_index_set:
                    running_sum += 1.0 / n_genes_in_set
                else:
                    running_sum -= 1.0 / n_genes_not_in_set

                abs_deviation = abs(running_sum)
                if abs_deviation > max_deviation:
                    max_deviation = abs_deviation

            scores[j] = max_deviation

        return scores

    def _compute_plage_score(
        self,
        gene_set_genes: List[str],
        gene_to_idx: Dict[str, int],
        expression_matrix: np.ndarray
    ) -> np.ndarray:
        """Calculate PLAGE scores from the first right singular vector."""
        # Locate gene-set members represented in the expression matrix.
        gene_indices = [gene_to_idx[g] for g in gene_set_genes if g in gene_to_idx]

        if len(gene_indices) == 0:
            n_samples = expression_matrix.shape[1]
            return np.zeros(n_samples)

        # Extract matched genes.
        sub_matrix = expression_matrix[gene_indices, :]  # (n_genes_in_set, n_samples)

        # Standardize each gene across samples.
        means = np.mean(sub_matrix, axis=1, keepdims=True)
        stds = np.std(sub_matrix, axis=1, ddof=1 if sub_matrix.shape[1] > 1 else 0, keepdims=True)
        stds[~np.isfinite(stds) | (stds == 0)] = 1.0
        z_matrix = (sub_matrix - means) / stds

        # Use NumPy SVD to avoid an additional scikit-learn dependency.
        n_components = min(len(gene_indices), z_matrix.shape[1])
        if n_components == 0:
            return np.zeros(expression_matrix.shape[1])

        try:
            # Samples are rows and matched genes are columns for SVD.
            U, S, Vt = np.linalg.svd(z_matrix.T, full_matrices=False)
            scores = U[:, 0] * S[0]
        except np.linalg.LinAlgError:
            # Fall back to mean standardized expression if SVD does not converge.
            scores = np.mean(z_matrix, axis=0)

        return scores

    def _compute_zscore(
        self,
        gene_set_genes: List[str],
        gene_to_idx: Dict[str, int],
        expression_matrix: np.ndarray
    ) -> np.ndarray:
        """Calculate combined standardized-expression scores for one gene set."""
        n_samples = expression_matrix.shape[1]

        # Locate gene-set members represented in the expression matrix.
        gene_indices = [gene_to_idx[g] for g in gene_set_genes if g in gene_to_idx]

        if len(gene_indices) == 0:
            return np.zeros(n_samples)

        # Extract matched genes.
        sub_matrix = expression_matrix[gene_indices, :]  # (n_genes_in_set, n_samples)

        # Standardize each gene across samples.
        means = np.mean(sub_matrix, axis=1, keepdims=True)
        stds = np.std(sub_matrix, axis=1, ddof=1 if sub_matrix.shape[1] > 1 else 0, keepdims=True)
        stds[~np.isfinite(stds) | (stds == 0)] = 1.0
        z_matrix = (sub_matrix - means) / stds

        # Aggregate standardized expression across matched genes.
        scores = np.mean(z_matrix, axis=0)

        return scores

    def analyze_matrix(
        self,
        expression_matrix: pd.DataFrame,
        gene_sets: Dict[str, Set[str]]
    ) -> pd.DataFrame:
        """Calculate the gene-set-by-sample activity matrix.
        
        Args:
            expression_matrix: Expression matrix with gene IDs as the index.
            gene_sets: Mapping from gene-set IDs to gene identifiers.
        
        Returns:
            Activity scores with gene sets in rows and samples in columns."""
        from allenricher.core.bioconductor import run_gsva

        return run_gsva(
            expression_matrix,
            gene_sets,
            method=self.method,
            kcdf=self.kcdf,
            tau=self.tau,
            min_size=self.min_size,
            max_size=self.max_size,
        )
