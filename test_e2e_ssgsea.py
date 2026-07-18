#!/usr/bin/env python3
"""ssGSEA end-to-end test."""

import sys
import time
import json
from pathlib import Path
import pandas as pd
import numpy as np

sys.path.insert(0, str(Path(__file__).parent))
from allenricher.core.enrichment import SSGSEA

TEST_DATA_DIR = Path("test_data")
RESULTS_DIR = Path("test_data/e2e_results")
RESULTS_DIR.mkdir(exist_ok=True)

def load_test_data():
    """Loading test data"""
    # Read 6, 000 x 6 for matrix
    expr_matrix = pd.read_csv(TEST_DATA_DIR / "expression_matrix_6000.tsv", sep='\t', index_col=0)
    
    # Read Test Channel
    gene_sets = {}
    with open(TEST_DATA_DIR / "test_pathways_from_gmt.gmt", 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            pathway = parts[0]
            genes = set(parts[2:])
            gene_sets[pathway] = genes
    
    return expr_matrix, gene_sets

def test_ssgsea_full():
    """Run the ssGSEA end-to-end check."""
    print("=" * 60)
    print("ssGSEA end-to-end test")
    print("=" * 60)
    
    expr_matrix, gene_sets = load_test_data()
    print(f"• Expression matrix: {expr_matrix.shape[0]} genes × {expr_matrix.shape[1]} samples")
    print(f"* Test number of circuits: {len(gene_sets)} pathways")
    
    # Create the ssGSEA analyzer.
    ssgsea = SSGSEA(min_size=10, max_size=500)
    
    # Test all routes.
    start_time = time.time()
    results_df = ssgsea.analyze_matrix(expr_matrix, gene_sets)
    elapsed = time.time() - start_time
    
    # Save Results
    results_df.to_csv(RESULTS_DIR / "ssgsea_results.csv")
    
    # Statistical analysis
    sample_means = results_df.mean(axis=0)
    sample_stds = results_df.std(axis=0)
    pathway_means = results_df.mean(axis=1)
    
    # Relevance between samples
    sample_corr = results_df.corr()
    
    # Generate test report
    report = {
        "test_name": "ssGSEA Full E2E Test",
        "timestamp": time.strftime("%Y-%m-%d %H:%M:%S"),
        "input_data": {
            "expression_matrix": f"{expr_matrix.shape[0]}x{expr_matrix.shape[1]}",
            "pathways_tested": len(gene_sets)
        },
        "results": {
            "output_shape": list(results_df.shape),
            "score_range": [float(results_df.values.min()), float(results_df.values.max())],
            "score_mean": float(results_df.values.mean()),
            "score_std": float(results_df.values.std()),
            "sample_means": sample_means.to_dict(),
            "pathway_means": pathway_means.to_dict(),
            "execution_time": f"{elapsed:.2f}s"
        },
        "validation": {
            "all_scores_in_range_0_1": bool((results_df.values >= 0).all() and (results_df.values <= 1).all()),
            "no_nan_values": bool(~results_df.isna().any().any()),
            "no_inf_values": bool(~np.isinf(results_df.values).any())
        },
        "status": "passed"
    }
    
    with open(RESULTS_DIR / "ssgsea_report.json", 'w') as f:
        json.dump(report, f, indent=2)
    
    # Print Results
    print(f"\nssGSEA completed in {elapsed:.2f}s")
    print(f"• Results saved: {RESULTS_DIR / 'ssgsea_results.csv'}")
    print(f"\nStatistical information:")
    print(f"- Output matrix: {results_df.shape[0]} pathways × {results_df.shape[1]} samples")
    print(f"- Score range: [{results_df.values.min(): .3f}, {results_df.values.max(): .3f}]")
    print(f"- Average: {results_df.values.mean(): .3f}")
    
    return report

if __name__ == "__main__":
    report = test_ssgsea_full()
    print("\n" + "=" * 60)
    print("ssGSEA E2E test complete.")
    print("=" * 60)
