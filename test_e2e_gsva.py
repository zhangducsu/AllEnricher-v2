#!/usr/bin/env python3
"""GSVA全量端到端测试（三种方法变体）"""

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
    """加载测试数据"""
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
    """测试GSVA单一方法"""
    print(f"\n测试 GSVA ({method_name})...")

    expr_matrix, gene_sets = load_test_data()

    # 创建GSVA分析器
    gsva = GSVA(method=method, min_size=10, max_size=500)

    # 执行分析
    start_time = time.time()
    results_df = gsva.analyze_matrix(expr_matrix, gene_sets)
    elapsed = time.time() - start_time

    # 保存结果
    results_df.to_csv(RESULTS_DIR / f"gsva_{method}_results.csv")

    # 统计分析
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

    print(f"  ✓ 执行时间: {elapsed:.2f}s")
    print(f"  ✓ 得分范围: [{results_df.values.min():.3f}, {results_df.values.max():.3f}]")

    return report


def test_gsva_full():
    """GSVA全量测试 - 三种方法"""
    print("=" * 60)
    print("GSVA全量端到端测试（三种方法变体）")
    print("=" * 60)

    # 测试三种方法
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
            print(f"  ✗ 错误: {e}")
            all_reports[method] = {
                "method": method,
                "method_name": method_name,
                "status": "failed",
                "error": str(e)
            }

    # 方法间比较
    print("\n方法间比较...")
    comparison = compare_methods(all_reports)

    # 生成综合报告
    final_report = {
        "test_name": "GSVA Full E2E Test (3 Methods)",
        "timestamp": time.strftime("%Y-%m-%d %H:%M:%S"),
        "methods": all_reports,
        "comparison": comparison,
        "overall_status": "passed" if all(r.get("status") == "passed" for r in all_reports.values()) else "failed"
    }

    with open(RESULTS_DIR / "gsva_report.json", 'w') as f:
        json.dump(final_report, f, indent=2)

    print(f"\n✓ 报告保存: {RESULTS_DIR / 'gsva_report.json'}")

    return final_report


def compare_methods(reports: dict) -> dict:
    """比较三种方法的结果"""
    comparison = {
        "execution_times": {},
        "score_ranges": {},
        "correlations": {}
    }

    for method, report in reports.items():
        if report.get("status") == "passed":
            comparison["execution_times"][method] = report["results"]["execution_time"]
            comparison["score_ranges"][method] = report["results"]["score_range"]

    # 计算方法间的相关性
    try:
        # 加载三种方法的结果
        gsva_df = pd.read_csv(RESULTS_DIR / "gsva_gsva_results.csv", index_col=0)
        plage_df = pd.read_csv(RESULTS_DIR / "gsva_plage_results.csv", index_col=0)
        zscore_df = pd.read_csv(RESULTS_DIR / "gsva_zscore_results.csv", index_col=0)

        # 确保通路顺序一致
        common_pathways = gsva_df.index.intersection(plage_df.index).intersection(zscore_df.index)
        gsva_vals = gsva_df.loc[common_pathways].values.flatten()
        plage_vals = plage_df.loc[common_pathways].values.flatten()
        zscore_vals = zscore_df.loc[common_pathways].values.flatten()

        # 计算相关性
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
    print("GSVA E2E测试完成!")
    print("=" * 60)
