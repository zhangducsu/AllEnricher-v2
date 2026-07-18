"""Tests for Gene Set Variation Analysis (GSVA).

Coverage includes output shape, empty input, gene-set size filters, supported methods,
non-overlapping identifiers, single-sample input, and finite activity scores.
"""

import pytest
import pandas as pd
import numpy as np
from typing import Set

from allenricher.core.gsva import GSVA


class TestGSVA:
    """Tests for the GSVA activity-matrix interface."""

    @pytest.fixture
    def expression_matrix(self):
        """Create a reproducible 100-gene by 3-sample expression matrix."""
        np.random.seed(42)
        n_genes = 100
        n_samples = 3
        gene_names = [f"GENE_{i:04d}" for i in range(1, n_genes + 1)]
        sample_names = ["Sample_A", "Sample_B", "Sample_C"]
        data = np.random.randn(n_genes, n_samples)
        return pd.DataFrame(data, index=gene_names, columns=sample_names)

    @pytest.fixture
    def gene_sets(self):
        """Create two 15-gene pathway fixtures."""
        return {
            "Pathway_1": {f"GENE_{i:04d}" for i in range(1, 16)},      # GENE_0001 ~ GENE_0015
            "Pathway_2": {f"GENE_{i:04d}" for i in range(16, 31)},     # GENE_0016 ~ GENE_0030
        }

    def test_gsva_basic(self, expression_matrix, gene_sets):
        """Return a two-pathway by three-sample activity matrix."""
        gsva = GSVA(method="gsva", min_size=10, max_size=500)
        result = gsva.analyze_matrix(expression_matrix, gene_sets)

        # Two pathways are scored across three samples.
        assert result.shape == (2, 3), f"Expected shape (2, 3) actual{result.shape}"
        # Preserve pathway and sample labels.
        assert list(result.index) == ["Pathway_1", "Pathway_2"]
        assert list(result.columns) == ["Sample_A", "Sample_B", "Sample_C"]

    def test_gsva_empty_input(self):
        """Return an empty table for an empty expression matrix."""
        gsva = GSVA(method="gsva")
        empty_matrix = pd.DataFrame()
        gene_sets = {"Pathway_1": {"GENE_0001"}}

        result = gsva.analyze_matrix(empty_matrix, gene_sets)

        assert result.empty, "An empty expression matrix should return an empty DataFrame"

    def test_gsva_small_gene_set(self, expression_matrix):
        """Exclude gene sets smaller than ``min_size`` after intersection."""
        # Create a gene set with only 5 genes (less than the default men_size=10)
        small_gene_sets = {
            "Small_Pathway": {f"GENE_{i:04d}" for i in range(1, 6)},  # Only 5 genes
        }

        gsva = GSVA(method="gsva", min_size=10, max_size=500)
        result = gsva.analyze_matrix(expression_matrix, small_gene_sets)

        assert result.empty, "Gene sets smaller than min_size should be excluded"

    def test_gsva_large_gene_set(self, expression_matrix):
        """Exclude gene sets larger than ``max_size`` after intersection."""
        # Create a gene set with 600 genes (greater than the default max_size=500)
        # Only 100 genes in the matrix are expressed, so smaller max_sizes are needed
        large_gene_sets = {
            "Large_Pathway": {f"GENE_{i:04d}" for i in range(1, 101)},  # 100 genes
        }

        gsva = GSVA(method="gsva", min_size=10, max_size=50)  # max_size=50 < 100
        result = gsva.analyze_matrix(expression_matrix, large_gene_sets)

        assert result.empty, "Gene sets larger than max_size should be excluded"

    def test_gsva_plage_method(self, expression_matrix, gene_sets):
        """PLAGE method test: Validate PLAGE method to produce properly"""
        gsva = GSVA(method="plage", min_size=10, max_size=500)
        result = gsva.analyze_matrix(expression_matrix, gene_sets)

        # Authenticate output shape
        assert result.shape == (2, 3), f"PLAGE Output Shape Expectation (2, 3), is{result.shape}"
        # Validation results are limited
        assert np.all(np.isfinite(result.values)), "The PLAGE results should be all in limited values"

    def test_gsva_zscore_method(self, expression_matrix, gene_sets):
        """Z-score Method Test: Verify that Z-score is capable of being exported properly"""
        gsva = GSVA(method="zscore", min_size=10, max_size=500)
        result = gsva.analyze_matrix(expression_matrix, gene_sets)

        # Authenticate output shape
        assert result.shape == (2, 3), f"Z-score Output Shape Expectation (2, 3), actually {result.shape}"
        # Validation results are limited
        assert np.all(np.isfinite(result.values)), "Z-score results should be all limited"

    def test_gsva_no_overlap(self, expression_matrix):
        """No intersection test: DataFrame will be returned when the gene pool and the expression matrix do not intersect"""
        # Create a gene pool that does not overlap with the expression matrix
        no_overlap_sets = {
            "No_Overlap_Pathway": {"UNKNOWN_GENE_1", "UNKNOWN_GENE_2", "UNKNOWN_GENE_3"},
        }

        gsva = GSVA(method="gsva", min_size=1, max_size=500)  # Lower Min_size to exclude size filters
        result = gsva.analyze_matrix(expression_matrix, no_overlap_sets)

        assert result.empty, "The uninterrupted genome should be skipped and returned to empty DataFrame"

    def test_gsva_single_sample(self):
        """Single sample test: a matrix of expressions containing only one sample"""
        np.random.seed(42)
        n_genes = 100
        gene_names = [f"GENE_{i:04d}" for i in range(1, n_genes + 1)]
        data = np.random.randn(n_genes, 1)
        single_sample_matrix = pd.DataFrame(data, index=gene_names, columns=["Only_Sample"])

        gene_sets = {
            "Pathway_1": {f"GENE_{i:04d}" for i in range(1, 16)},
        }

        gsva = GSVA(method="gsva", min_size=10, max_size=500)
        result = gsva.analyze_matrix(single_sample_matrix, gene_sets)

        # Verify output shape: 1 route x 1 sample
        assert result.shape == (1, 1), f"Single sample output shape expectation (1, 1) actual{result.shape}"
        # Validation results are limited
        assert np.all(np.isfinite(result.values)), "The results of the single sample should be all limited"

    def test_gsva_result_range(self, expression_matrix, gene_sets):
        """Results for reasonableness test: Outcome value should be limited (non-NAN, non-Inf)"""
        gsva = GSVA(method="gsva", min_size=10, max_size=500)
        result = gsva.analyze_matrix(expression_matrix, gene_sets)

        # Verify that all values are limited
        assert np.all(np.isfinite(result.values)), "GSVA results should not contain nN or inf"

        # Validation results are within reasonable range (GSVA scores are usually close to [-1, 1], but the range depends on data)
        # It's not extreme.
        assert np.all(np.abs(result.values) < 1e6), "GSVA Results should not be extremes"


class TestGSVAInit:
    """GSVA Initialization Parameter Validation Test"""

    def test_invalid_method(self):
        """Reject an unknown GSVA method name."""
        with pytest.raises(ValueError, match="Unsupported GSVA method"):
            GSVA(method="invalid_method")

    def test_invalid_kcdf(self):
        """Test invalid nuclear function type to be thrown away ValueError"""
        with pytest.raises(ValueError, match="Unsupported GSVA kernel"):
            GSVA(method="gsva", kcdf="Invalid")

    def test_default_parameters(self):
        """Test Default Parameter Value"""
        gsva = GSVA()
        assert gsva.method == "gsva"
        assert gsva.kcdf == "Gaussian"
        assert gsva.tau == 1.0
        assert gsva.min_size == 1
        assert gsva.max_size is None


class TestGSVACalculatePvalue:
    """GSVA calculate_pvalue method test"""

    def test_returns_nan(self):
        """The calculate_pvalue of GSVA returns the n-N"""
        gsva = GSVA()
        result = gsva.calculate_pvalue(10, 50, 100, 1000)
        assert np.isnan(result), "The calculate_pvalue of GSVA returns the n-N"


class TestGSVACalculateEnrichment:
    """GSVA calculate_enrichment method test (compatible base-type interface)"""

    def test_returns_result_for_valid_input(self):
        """Valid input should return EnrichmentResult"""
        gsva = GSVA(min_size=2, max_size=100)
        gene_set = {"GENE_1", "GENE_2", "GENE_3"}
        background_set = {"GENE_1", "GENE_2", "GENE_3", "GENE_4", "GENE_5"}
        term_genes = {"GENE_1", "GENE_2", "GENE_3"}

        result = gsva.calculate_enrichment(
            gene_set=gene_set,
            background_set=background_set,
            term_genes=term_genes,
            term_name="Test_Pathway",
            term_id="TEST:001",
            database="TEST"
        )

        assert result is not None
        assert result.term_id == "TEST:001"
        assert result.gene_count == 3

    def test_returns_none_for_small_gene_set(self):
        """The gene set is too small to return to None"""
        gsva = GSVA(min_size=10, max_size=500)
        gene_set = {"GENE_1"}
        background_set = {"GENE_1", "GENE_2"}
        term_genes = {"GENE_1"}

        result = gsva.calculate_enrichment(
            gene_set=gene_set,
            background_set=background_set,
            term_genes=term_genes,
            term_name="Small_Pathway",
            term_id="TEST:002",
            database="TEST"
        )

        assert result is None


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
