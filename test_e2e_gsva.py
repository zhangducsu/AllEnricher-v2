#!/usr/bin/env python3
"""GSVA Full-Stand-Endpoint Test (Three Method Variations)"""

import sys
import time
import json
from pathlib import Path
import pandas as pd
import numpy as np

sys.path.insert(0, str(Path(__file__).parent))
from allenricher.core.gsva import GSVA

TEST_DATA_DIR = Path("test_data")
RESULTS_DIR = Path("test_data/e2e_results")
RESULTS_DIR.mkdir(exist_ok=True)


def load_test_data():
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


def test_gsva_method(method_name: str, method: str) -> dict:
    """Test GSVA single method"""
    print(f"\nTest GSVA (Standard){method_name})...")

    expr_matrix, gene_sets = load_test_data()

    # Create GSVA Analyzer
    gsva = GSVA(method=method, min_size=10, max_size=500)

    # Implementation analysis
    start_time = time.time()
    results_df = gsva.analyze_matrix(expr_matrix, gene_sets)
    elapsed = time.time() - start_time

    # Save Results
    results_df.to_csv(RESULTS_DIR / f"gsva_{method}_results.csv")

    # Statistical analysis
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

    print(f"* Implementation time: {elapsed: .2f}s")
    print(f"* Score range: [{results_df.values.min(): .3f}, {results_df.values.max(): .3f}]")

    return report


def test_gsva_full():
    """GSVA Full Test - Three Methods"""
    print("=" * 60)
    print("GSVA Full-Stand-Endpoint Test (Three Method Variations)")
    print("=" * 60)

    # Test three ways.
    methods = [
        ("Random Walk (Default)", "gsva"),
        ("PLAGE", "plage"),
        ("Z-score", "zscore")
    ]

    all_reports = {}

    for method_name, method in methods:
        try:
            report = test_gsva_method(method_name, method)
            all_reports[method] = report
        except Exception as e:
            print(f"Error: {e}")
            all_reports[method] = {
                "method": method,
                "method_name": method_name,
                "status": "failed",
                "error": str(e)
            }

    # Inter-methodological comparison
    print("\nInter-method comparison...")
    comparison = compare_methods(all_reports)

    # Generate a synthesis report
    final_report = {
        "test_name": "GSVA Full E2E Test (3 Methods)",
        "timestamp": time.strftime("%Y-%m-%d %H:%M:%S"),
        "methods": all_reports,
        "comparison": comparison,
        "overall_status": "passed" if all(r.get("status") == "passed" for r in all_reports.values()) else "failed"
    }

    with open(RESULTS_DIR / "gsva_report.json", 'w') as f:
        json.dump(final_report, f, indent=2)

    print(f"\n* Report preservation: {RESULTS_DIR / 'gsva_report.json'}")

    return final_report


def compare_methods(reports: dict) -> dict:
    """Comparison of results of the three methods"""
    comparison = {
        "execution_times": {},
        "score_ranges": {},
        "correlations": {}
    }

    for method, report in reports.items():
        if report.get("status") == "passed":
            comparison["execution_times"][method] = report["results"]["execution_time"]
            comparison["score_ranges"][method] = report["results"]["score_range"]

    # Relevance between methods of calculation
    try:
        # Results of loading of three methods
        gsva_df = pd.read_csv(RESULTS_DIR / "gsva_gsva_results.csv", index_col=0)
        plage_df = pd.read_csv(RESULTS_DIR / "gsva_plage_results.csv", index_col=0)
        zscore_df = pd.read_csv(RESULTS_DIR / "gsva_zscore_results.csv", index_col=0)

        # Make sure the traffic sequence is consistent.
        common_pathways = gsva_df.index.intersection(plage_df.index).intersection(zscore_df.index)
        gsva_vals = gsva_df.loc[common_pathways].values.flatten()
        plage_vals = plage_df.loc[common_pathways].values.flatten()
        zscore_vals = zscore_df.loc[common_pathways].values.flatten()

        # Calculation relevance
        comparison["correlations"] = {
            "gsva_vs_plage": float(np.corrcoef(gsva_vals, plage_vals)[0, 1]),
            "gsva_vs_zscore": float(np.corrcoef(gsva_vals, zscore_vals)[0, 1]),
            "plage_vs_zscore": float(np.corrcoef(plage_vals, zscore_vals)[0, 1])
        }
    except Exception as e:
        comparison["correlations"] = {"error": str(e)}

    return comparison


if __name__ == "__main__":
    report = test_gsva_full()
    print("\n" + "=" * 60)
    print("GSVA E2E test complete!")
    print("=" * 60)
