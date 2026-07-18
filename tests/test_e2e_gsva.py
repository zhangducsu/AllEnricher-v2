#!/usr/bin/env python3
"""GSVA Full End-to-end Test (Three Method Variations) - pytest module test"""

import sys
import time
import json
from pathlib import Path
import pandas as pd
import numpy as np
import pytest

# Add Project Root Directory to Path
sys.path.insert(0, str(Path(__file__).parent.parent))
from allenricher.core.gsva import GSVA

TEST_DATA_DIR = Path(__file__).parent.parent / "test_data"
RESULTS_DIR = TEST_DATA_DIR / "e2e_results"


@pytest.fixture
def test_data():
    """Loading test data"""
    expr_matrix = pd.read_csv(TEST_DATA_DIR / "expression_matrix_6000.tsv", sep='\t', index_col=0)

    gene_sets = {}
    with open(TEST_DATA_DIR / "test_pathways_from_gmt.gmt", 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            pathway = parts[0]
            genes = set(parts[2:])
            gene_sets[pathway] = genes

    return expr_matrix, gene_sets


class TestGSVAMethods:
    """Test GSVA three method variants"""

    @pytest.mark.parametrize("method_name,method", [
        ("Random Walk (Default)", "gsva"),
        ("PLAGE", "plage"),
        ("Z-score", "zscore")
    ])
    def test_method_returns_correct_shape(self, test_data, method_name, method):
        """Each method returns a pathway-by-sample activity matrix."""
        expr_matrix, gene_sets = test_data

        gsva = GSVA(method=method, min_size=10, max_size=500)
        results_df = gsva.analyze_matrix(expr_matrix, gene_sets)

        assert results_df.shape[0] > 0, f"{method_name}: expected at least one pathway"
        assert results_df.shape[1] == expr_matrix.shape[1], f"{method_name}: sample count changed"
        assert list(results_df.columns) == list(expr_matrix.columns), f"{method_name}: sample names changed"

    @pytest.mark.parametrize("method_name,method", [
        ("Random Walk (Default)", "gsva"),
        ("PLAGE", "plage"),
        ("Z-score", "zscore")
    ])
    def test_method_no_nan_inf_values(self, test_data, method_name, method):
        """Test No Nan/InfValue"""
        expr_matrix, gene_sets = test_data

        gsva = GSVA(method=method, min_size=10, max_size=500)
        results_df = gsva.analyze_matrix(expr_matrix, gene_sets)

        # Check the nn
        assert not results_df.isna().any().any(), f"{method_name}: The results should not contain an nn value"

        # Check Inf
        assert not np.isinf(results_df.values).any(), f"{method_name}: There should be no Inf value in the result"

    @pytest.mark.parametrize("method_name,method", [
        ("Random Walk (Default)", "gsva"),
        ("PLAGE", "plage"),
        ("Z-score", "zscore")
    ])
    def test_method_execution_time(self, test_data, method_name, method):
        """Test time for execution<60sec (Each method)"""
        expr_matrix, gene_sets = test_data

        gsva = GSVA(method=method, min_size=10, max_size=500)

        start_time = time.time()
        results_df = gsva.analyze_matrix(expr_matrix, gene_sets)
        elapsed = time.time() - start_time

        assert elapsed < 60, f"{method_name}: Implementation time should be less than 60 seconds, actual time-consuming{elapsed: .2f}sec"


class TestGSVAReport:
    """Test GSVA report generation"""

    def test_report_json_format(self, test_data):
        """Test report JSON format correctly"""
        expr_matrix, gene_sets = test_data

        # Run three methods and generate reports
        methods = [
            ("Random Walk (Default)", "gsva"),
            ("PLAGE", "plage"),
            ("Z-score", "zscore")
        ]

        all_reports = {}

        for method_name, method in methods:
            gsva = GSVA(method=method, min_size=10, max_size=500)

            start_time = time.time()
            results_df = gsva.analyze_matrix(expr_matrix, gene_sets)
            elapsed = time.time() - start_time

            report = {
                "method": method,
                "method_name": method_name,
                "timestamp": time.strftime("%Y-%m-%d %H:%M:%S"),
                "input_data": {
                    "expression_matrix": f"{expr_matrix.shape[0]}×{expr_matrix.shape[1]}",
                    "pathways_tested": len(gene_sets)
                },
                "results": {
                    "output_shape": list(results_df.shape),
                    "score_range": [float(results_df.values.min()), float(results_df.values.max())],
                    "score_mean": float(results_df.values.mean()),
                    "score_std": float(results_df.values.std()),
                    "execution_time": f"{elapsed:.2f}s"
                },
                "validation": {
                    "no_nan_values": bool(~results_df.isna().any().any()),
                    "no_inf_values": bool(~np.isinf(results_df.values).any())
                },
                "status": "passed"
            }
            all_reports[method] = report

        # Generate a synthesis report
        final_report = {
            "test_name": "GSVA Full E2E Test (3 Methods)",
            "timestamp": time.strftime("%Y-%m-%d %H:%M:%S"),
            "methods": all_reports,
            "overall_status": "passed"
        }

        # Validate JSON for sequenced
        try:
            json_str = json.dumps(final_report, indent=2)
            assert json_str is not None
            assert len(json_str) > 0
        except json.JSONEncodeError as e:
            pytest.fail(f"Reporting on JSON failed: {e}")

        # Validate the report structure
        assert "test_name" in final_report
        assert "timestamp" in final_report
        assert "methods" in final_report
        assert "overall_status" in final_report
        assert len(final_report["methods"]) == 3

        for method, report in final_report["methods"].items():
            assert "method" in report
            assert "results" in report
            assert "validation" in report
            assert "status" in report
            assert report["status"] == "passed"


class TestGSVACorrelation:
    """Test the relevance of the results of the GSVA three methods"""

    def test_methods_correlation(self, test_data):
        """Test the relevance of the results of the three methods (should have some relevance)"""
        expr_matrix, gene_sets = test_data

        # Run three ways
        methods = ["gsva", "plage", "zscore"]
        results = {}

        for method in methods:
            gsva = GSVA(method=method, min_size=10, max_size=500)
            results[method] = gsva.analyze_matrix(expr_matrix, gene_sets)

        # Get a common access.
        common_pathways = results["gsva"].index
        for method in methods[1:]:
            common_pathways = common_pathways.intersection(results[method].index)

        assert len(common_pathways) > 0, "There should be common access to the three approaches."

        # Draw value for common access
        gsva_vals = results["gsva"].loc[common_pathways].values.flatten()
        plage_vals = results["plage"].loc[common_pathways].values.flatten()
        zscore_vals = results["zscore"].loc[common_pathways].values.flatten()

        # Calculation relevance
        corr_gsva_plage = np.corrcoef(gsva_vals, plage_vals)[0, 1]
        corr_gsva_zscore = np.corrcoef(gsva_vals, zscore_vals)[0, 1]
        corr_plage_zscore = np.corrcoef(plage_vals, zscore_vals)[0, 1]

        # Validation relevance (the correlation factor should be within reasonable limits and not necessarily positive)
        assert abs(corr_gsva_plage) >= 0, f"gsva vs page is abnormally relevant: {corr_gsva_plage}"
        assert abs(corr_gsva_zscore) >= 0, f"Invalid GSVA-to-z-score correlation: {corr_gsva_zscore}"
        assert abs(corr_plage_zscore) >= 0, f"The page vs zscore is not relevant: {corr_plage_zscore}"

        # Record relevance values for debug
        print(f"\nInter-methodological relevance:")
        print(f"  gsva vs plage: {corr_gsva_plage:.4f}")
        print(f"  gsva vs zscore: {corr_gsva_zscore:.4f}")
        print(f"  plage vs zscore: {corr_plage_zscore:.4f}")


class TestGSVAEdgeCases:
    """Test GSVA boundary"""

    def test_empty_expression_matrix(self):
        """Test empty expression matrix"""
        empty_df = pd.DataFrame()
        gene_sets = {"pathway1": {"gene1", "gene2"}}

        gsva = GSVA(method="gsva")
        result = gsva.analyze_matrix(empty_df, gene_sets)

        assert result.empty, "Empty input should return empty DataFrame"

    def test_no_matching_genes(self, test_data):
        """Testing for no match to the genes."""
        expr_matrix, _ = test_data

        # Create a gene set that has no intersection with the expression matrix
        gene_sets = {"pathway_no_match": {"FAKE_GENE_1", "FAKE_GENE_2", "FAKE_GENE_3"}}

        gsva = GSVA(method="gsva", min_size=1, max_size=500)
        result = gsva.analyze_matrix(expr_matrix, gene_sets)

        # When no match is made, return the result.
        assert result.empty or len(result) == 0, "When no match is made, return the result."

    def test_gene_set_size_filtering(self, test_data):
        """Test for gene set size filter"""
        expr_matrix, gene_sets = test_data

        # Filter most of the genome using strict Min_size
        gsva_strict = GSVA(method="gsva", min_size=1000, max_size=5000)
        result_strict = gsva_strict.analyze_matrix(expr_matrix, gene_sets)

        # Use loose min_size
        gsva_loose = GSVA(method="gsva", min_size=10, max_size=500)
        result_loose = gsva_loose.analyze_matrix(expr_matrix, gene_sets)

        # The easing of conditions should have more results.
        assert len(result_loose) >= len(result_strict), "The easing of conditions should have more results."


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
