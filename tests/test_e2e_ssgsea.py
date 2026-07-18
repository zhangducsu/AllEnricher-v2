#!/usr/bin/env python3
"""SGSEA E2E test - Unit test"""

import sys
import time
import json
from pathlib import Path
import pandas as pd
import numpy as np
import unittest

# Add Project Root Directory to Path
sys.path.insert(0, str(Path(__file__).parent.parent))
from allenricher.core.enrichment import SSGSEA

TEST_DATA_DIR = Path(__file__).parent.parent / "test_data"
RESULTS_DIR = TEST_DATA_DIR / "e2e_results"


class TestSsGSEAEndToEnd(unittest.TestCase):
    """SGSSEA End-to-end Test Class"""
    
    @classmethod
    def setUpClass(cls):
        """Test Class Initialisation, Load Test Data"""
        cls.expr_matrix = pd.read_csv(
            TEST_DATA_DIR / "expression_matrix_6000.tsv", 
            sep='\t', 
            index_col=0
        )
        
        # Read Test Channel
        cls.gene_sets = {}
        with open(TEST_DATA_DIR / "test_pathways_from_gmt.gmt", 'r') as f:
            for line in f:
                parts = line.strip().split('\t')
                pathway = parts[0]
                genes = set(parts[2:])
                cls.gene_sets[pathway] = genes
        
        # Create ssGSEA analyser
        cls.ssgsea = SSGSEA(min_size=10, max_size=500)
        expression_genes = set(cls.expr_matrix.index.astype(str))
        cls.expected_pathways = {
            pathway
            for pathway, genes in cls.gene_sets.items()
            if cls.ssgsea.min_size
            <= len(genes & expression_genes)
            <= cls.ssgsea.max_size
        }
        
        # Run analysis and time recording
        cls.start_time = time.time()
        cls.results_df = cls.ssgsea.analyze_matrix(cls.expr_matrix, cls.gene_sets)
        cls.elapsed_time = time.time() - cls.start_time
        
        # Ensure that results directory exists
        RESULTS_DIR.mkdir(exist_ok=True)
        
        # Save Results
        cls.results_df.to_csv(RESULTS_DIR / "ssgsea_results.csv")
    
    def test_01_output_shape(self):
        """Test result matrix shape correct"""
        expected_pathways = len(self.expected_pathways)
        expected_samples = self.expr_matrix.shape[1]
        
        self.assertEqual(
            self.results_df.shape[0], 
            expected_pathways,
            f"The number of routes does not match: expect {expected_pathways}, actual {self.results_df.shape[0]}"
        )
        self.assertEqual(
            self.results_df.shape[1], 
            expected_samples,
            f"Samples do not match: Expect {expected_samples}, actual {self.results_df.shape[1]}"
        )
        print(f"* Output matrix shape is correct: {self.results_df.shape}")
    
    def test_02_score_range(self):
        """Test scores are within [-1, 1] (sGSEA range)"""
        min_score = self.results_df.values.min()
        max_score = self.results_df.values.max()
        
        self.assertGreaterEqual(min_score, -1.0, f"Minimum score {min_score} less than -1")
        self.assertLessEqual(max_score, 1.0, f"Max. score.{max_score}greater than 1")
        print(f"* The score range is correct: [ =]{min_score: .3f}, {max_score: .3f}]")
    
    def test_03_no_nan_values(self):
        """Test No NaN Value"""
        has_nan = self.results_df.isna().any().any()
        self.assertFalse(has_nan, "The result contains the nn value")
        print("* No nn value")
    
    def test_04_no_inf_values(self):
        """No Inf for Test"""
        has_inf = np.isinf(self.results_df.values).any()
        self.assertFalse(has_inf, "The result contains Inf values")
        print("*No Inf value")
    
    def test_05_execution_time(self):
        """Test Implementation Time<30sec"""
        max_time = 30.0
        self.assertLess(
            self.elapsed_time, 
            max_time,
            f"Implementation time{self.elapsed_time: .2f}s Over{max_time}s"
        )
        print(f"* Implementation time: {self.elapsed_time: .2f}s (limitation: {max_time}s)")
    
    def test_06_report_json_format(self):
        """Test report JSON format correctly"""
        # Generate test report
        sample_means = self.results_df.mean(axis=0)
        pathway_means = self.results_df.mean(axis=1)
        
        report = {
            "test_name": "ssGSEA Full E2E Test",
            "timestamp": time.strftime("%Y-%m-%d %H:%M:%S"),
            "input_data": {
                "expression_matrix": f"{self.expr_matrix.shape[0]}x{self.expr_matrix.shape[1]}",
                "pathways_tested": len(self.gene_sets)
            },
            "results": {
                "output_shape": list(self.results_df.shape),
                "score_range": [
                    float(self.results_df.values.min()), 
                    float(self.results_df.values.max())
                ],
                "score_mean": float(self.results_df.values.mean()),
                "score_std": float(self.results_df.values.std()),
                "sample_means": sample_means.to_dict(),
                "pathway_means": pathway_means.to_dict(),
                "execution_time": f"{self.elapsed_time:.2f}s"
            },
            "validation": {
                "all_scores_in_range": bool(
                    (self.results_df.values >= -1).all() and 
                    (self.results_df.values <= 1).all()
                ),
                "no_nan_values": bool(~self.results_df.isna().any().any()),
                "no_inf_values": bool(~np.isinf(self.results_df.values).any())
            },
            "status": "passed"
        }
        
        # Save Report
        report_path = RESULTS_DIR / "ssgsea_report.json"
        with open(report_path, 'w') as f:
            json.dump(report, f, indent=2)
        
        # Validate JSON format
        self.assertTrue(report_path.exists(), "The report document was not generated")
        
        # Verify that JSON can be loaded correctly
        with open(report_path, 'r') as f:
            loaded_report = json.load(f)
        
        # Validate the report structure
        required_keys = ["test_name", "timestamp", "input_data", "results", "validation", "status"]
        for key in required_keys:
            self.assertIn(key, loaded_report, f"The report lacks the necessary field: {key}")
        
        print(f"* The report is in the correct form and is maintained in: {report_path}")
    
    def test_07_pathway_names_match(self):
        """Test the pass name."""
        result_pathways = set(self.results_df.index)
        expected_pathways = self.expected_pathways
        
        self.assertEqual(
            result_pathways, 
            expected_pathways,
            "The traffic name in the result does not match input"
        )
        print("* The name of the circuit matches")
    
    def test_08_sample_names_match(self):
        """Test sample name matches"""
        result_samples = set(self.results_df.columns)
        expected_samples = set(self.expr_matrix.columns)
        
        self.assertEqual(
            result_samples, 
            expected_samples,
            "The sample name in the result does not match input"
        )
        print("* Sample name matches")


if __name__ == "__main__":
    # Create Test Package
    loader = unittest.TestLoader()
    suite = loader.loadTestsFromTestCase(TestSsGSEAEndToEnd)
    
    # Run test using TextTestRunner
    runner = unittest.TextTestRunner(verbosity=2)
    result = runner.run(suite)
    
    # Output Summary
    print("\n" + "=" * 60)
    print("SDSSEA E2E module test test summary")
    print("=" * 60)
    print(f"Total number of tests: {result.testsRun}")
    print(f"Adopted: {result.testsRun - len(result.failures) - len(result.errors)}")
    print(f"Failed: {len(result.failures)}")
    print(f"Error: {len(result.errors)}")
    
    if result.wasSuccessful():
        print("\nAll ssGSEA E2E checks passed.")
    else:
        print("\nTest failed, please see the detailed output")
    
    print("=" * 60)
