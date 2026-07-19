"""Extended tests for GSEA and ssGSEA analysis behavior.

Coverage includes ES and NES calculations, permutation P values, leading-edge genes,
empty inputs, expression-matrix interfaces, score ranges, and single-sample analysis.
"""

import pytest
import pandas as pd
import numpy as np

from allenricher.core.enrichment import GSEA, SSGSEA


# ============================================================
# Synthetic test fixtures
# ============================================================

def _make_ranked_genes_with_enrichment_top(gene_set_size=15, total_genes=100, seed=42):
    """Place one gene set at the top of a synthetic ranked list."""
    rng = np.random.default_rng(seed)
    gene_set = {f"GENE_{i}" for i in range(gene_set_size)}
    other_genes = [f"GENE_{i}" for i in range(gene_set_size, total_genes)]
    rng.shuffle(other_genes)
    # Gene-set members occupy the leading positions.
    ranked_genes = list(gene_set) + other_genes
    return ranked_genes, gene_set


def _make_ranked_genes_with_enrichment_bottom(gene_set_size=15, total_genes=100, seed=42):
    """Place one gene set at the bottom of a synthetic ranked list."""
    rng = np.random.default_rng(seed)
    gene_set = {f"GENE_{i}" for i in range(total_genes - gene_set_size, total_genes)}
    other_genes = [f"GENE_{i}" for i in range(total_genes - gene_set_size)]
    rng.shuffle(other_genes)
    # Gene-set members occupy the trailing positions.
    ranked_genes = other_genes + list(gene_set)
    return ranked_genes, gene_set


def _make_expression_matrix(n_genes=50, n_samples=3, seed=42):
    """Create a synthetic gene-by-sample expression matrix."""
    rng = np.random.default_rng(seed)
    genes = [f"GENE_{i}" for i in range(n_genes)]
    samples = [f"SAMPLE_{i}" for i in range(n_samples)]
    data = rng.random((n_genes, n_samples))
    return pd.DataFrame(data, index=genes, columns=samples)


def _make_gene_sets(n_pathways=2, genes_per_pathway=10, total_genes=50, seed=42):
    """Create synthetic pathway gene sets for tests."""
    rng = np.random.default_rng(seed)
    gene_pool = [f"GENE_{i}" for i in range(total_genes)]
    gene_sets = {}
    for i in range(n_pathways):
        rng.shuffle(gene_pool)
        gene_sets[f"PATHWAY_{i}"] = set(gene_pool[:genes_per_pathway])
    return gene_sets


# ============================================================
# GSEA Test
# ============================================================

class TestGSEAEnrichmentScore:
    """Tests for GSEA enrichment-score calculations."""

    def test_gsea_enrichment_score_basic(self):
        """Return an ES near one when every hit precedes every miss."""
        gsea = GSEA(permutations=100, seed=42)

        # All five gene-set members precede the remaining genes.
        gene_set = {"A", "B", "C", "D", "E"}
        ranked_genes = ["A", "B", "C", "D", "E"] + [f"X_{i}" for i in range(15)]

        es, leading_edge = gsea.calculate_enrichment_score(ranked_genes, gene_set)

        # Concentrating all hits at the top should produce an ES near one.
        assert es > 0.9, f"ES should be close to 1.0, actually {es}"
        # Every hit belongs to the leading edge in this construction.
        assert set(leading_edge) == gene_set

    def test_gsea_nes_positive(self):
        """Return a positive NES for enrichment at the top of the ranking."""
        gsea = GSEA(permutations=100, seed=42)

        ranked_genes, gene_set = _make_ranked_genes_with_enrichment_top(
            gene_set_size=15, total_genes=100, seed=42
        )

        es, nes, pvalue, leading_edge = gsea.calculate_normalized_es(
            ranked_genes, gene_set
        )

        assert nes > 0, f"NES should be positive when gathering for the rich, actually{nes}"
        assert es > 0, f"ES should be positive when it is going to be rich, actually is{es}"

    def test_gsea_nes_negative(self):
        """When the negatives are concentrated NES should be negative

Attention.: Current GSEA It's... ES Calculate only trace maximum value (Heading),
So when the gene pool is at the bottom, ES Close. 0.
This test validates the NES behaviour when ES is 0 or close to 0.
        """
        gsea = GSEA(permutations=100, seed=42)

        ranked_genes, gene_set = _make_ranked_genes_with_enrichment_bottom(
            gene_set_size=15, total_genes=100, seed=42
        )

        es, nes, pvalue, leading_edge = gsea.calculate_normalized_es(
            ranked_genes, gene_set
        )

        # When the gene set is at the bottom, as only the maximum trace is currently being achieved,
        # ES should be close to 0 (miss first reduced, hit later accumulated but not above previous peak)
        assert es < 0, f"The Es with a rich bottom should be negative, actually {es}"
        assert nes < 0, f"The NES with a rich bottom should be negative, actually{nes}"
        assert leading_edge

    def test_gsea_permutation_pvalue(self):
        """Replace the check p values should be within [0, 1]"""
        gsea = GSEA(permutations=100, seed=42)

        ranked_genes, gene_set = _make_ranked_genes_with_enrichment_top(
            gene_set_size=15, total_genes=100, seed=42
        )

        es, nes, pvalue, leading_edge = gsea.calculate_normalized_es(
            ranked_genes, gene_set
        )

        assert 0.0 <= pvalue <= 1.0, f"p Value should be within [0, 1] and actual{pvalue}"

    def test_gsea_leading_edge(self):
        """The frontier gene should all belong to the gene pool."""
        gsea = GSEA(permutations=100, seed=42)

        ranked_genes, gene_set = _make_ranked_genes_with_enrichment_top(
            gene_set_size=15, total_genes=100, seed=42
        )

        es, nes, pvalue, leading_edge = gsea.calculate_normalized_es(
            ranked_genes, gene_set
        )

        # The frontier gene should all belong to the gene pool.
        for gene in leading_edge:
            assert gene in gene_set, f"Precipitous genes{gene}Not a gene set."

    def test_gsea_empty_gene_set(self):
        """The empty gene set should return ES=0"""
        gsea = GSEA(permutations=100, seed=42)

        ranked_genes = [f"GENE_{i}" for i in range(50)]
        gene_set = set()  # Empty gene set.

        es, leading_edge = gsea.calculate_enrichment_score(ranked_genes, gene_set)

        assert es == 0.0, f"The Es for the empty gene set should be 0, actually{es}"
        assert leading_edge == [], f"The frontier gene of the empty gene set should be an empty list"


# ============================================================
# GSEA Expression Matrix Analysis Test
# ============================================================

class TestGSEAAnalyzeMatrix:
    """GSEA expression matrix analysis interface test"""

    def test_gsea_analyze_matrix(self):
        """Express matrix analysis output shape correct"""
        gsea = GSEA(permutations=50, seed=42)

        expression_matrix = _make_expression_matrix(n_genes=50, n_samples=3, seed=42)
        gene_sets = _make_gene_sets(n_pathways=2, genes_per_pathway=10, total_genes=50, seed=42)

        result = gsea.analyze_matrix(expression_matrix, gene_sets)

        # Confirm output shape: Line = Line = Line = Sample
        assert result.shape == (2, 3), f"Output shape should be (2, 3), actual{result.shape}"
        # Validation of line names and listings
        assert list(result.index) == ["PATHWAY_0", "PATHWAY_1"]
        assert list(result.columns) == ["SAMPLE_0", "SAMPLE_1", "SAMPLE_2"]


# ============================================================
# SGSEA test
# ============================================================

class TestSSGSEAAnalyzeMatrix:
    """ssGSEA expression matrix analysis interface test"""

    def test_ssgsea_analyze_matrix(self):
        """ssGSEA expression matrix analysis output shape is correct"""
        ssgsea = SSGSEA(min_size=1, max_size=500)

        expression_matrix = _make_expression_matrix(n_genes=50, n_samples=3, seed=42)
        gene_sets = _make_gene_sets(n_pathways=2, genes_per_pathway=10, total_genes=50, seed=42)

        result = ssgsea.analyze_matrix(expression_matrix, gene_sets)

        # Confirm output shape: Line = Line = Line = Sample
        assert result.shape == (2, 3), f"Output shape should be (2, 3), actual{result.shape}"
        # Validation of line names and listings
        assert list(result.index) == ["PATHWAY_0", "PATHWAY_1"]
        assert list(result.columns) == ["SAMPLE_0", "SAMPLE_1", "SAMPLE_2"]

    def test_ssgsea_nes_range(self):
        """sGSEA NES should be within [1, 1]"""
        ssgsea = SSGSEA(min_size=1, max_size=500)

        expression_matrix = _make_expression_matrix(n_genes=50, n_samples=5, seed=42)
        gene_sets = _make_gene_sets(n_pathways=3, genes_per_pathway=10, total_genes=50, seed=42)

        result = ssgsea.analyze_matrix(expression_matrix, gene_sets)

        # Verify all NES values in [-1, 1]
        for col in result.columns:
            for val in result[col]:
                assert -1.0 <= val <= 1.0, f"sGSEA NES should be within [-1, 1] range, actually{val}"

    def test_ssgsea_single_sample(self):
        """Single sample analysis should be working."""
        ssgsea = SSGSEA(min_size=1, max_size=500)

        # There's only one sample of the matrix of expression.
        expression_matrix = _make_expression_matrix(n_genes=50, n_samples=1, seed=42)
        gene_sets = _make_gene_sets(n_pathways=2, genes_per_pathway=10, total_genes=50, seed=42)

        result = ssgsea.analyze_matrix(expression_matrix, gene_sets)

        # Confirm output shape: Line = Roads, Column = 1
        assert result.shape == (2, 1), f"The output shape should be (2, 1) and the actual {result.shape}"
        assert list(result.columns) == ["SAMPLE_0"]


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
