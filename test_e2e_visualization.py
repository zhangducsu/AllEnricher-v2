#!/usr/bin/env python3
"""可视化集成端到端测试"""

import sys
import time
import json
from pathlib import Path
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))

from allenricher.visualization.gsea_plots import (
    plot_gsea_enrichment, plot_gsea_nes_barplot, plot_gsea_dotplot
)
from allenricher.visualization.gsva_plots import (
    plot_pathway_heatmap, plot_group_comparison, plot_sample_correlation
)
from allenricher.visualization.common_plots import (
    plot_enrichment_network, plot_volcano, plot_method_comparison
)

TEST_DATA_DIR = Path("test_data")
RESULTS_DIR = Path("test_data/e2e_results")
OUTPUT_DIR = Path("test_data/e2e_results/plots")
OUTPUT_DIR.mkdir(exist_ok=True)


def load_gene_sets_from_gmt(gmt_file):
    """从GMT文件加载基因集"""
    gene_sets = {}
    with open(gmt_file, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 3:
                gene_sets[parts[0]] = set(parts[2:])
    return gene_sets


def test_gsea_visualizations():
    """测试GSEA可视化"""
    print("\n" + "=" * 60)
    print("测试GSEA可视化")
    print("=" * 60)

    # 加载GSEA结果
    gsea_results = pd.read_csv(RESULTS_DIR / "gsea_results.csv")
    ranked_genes_df = pd.read_csv(TEST_DATA_DIR / "ranked_genes.tsv", sep='\t')

    tests = []

    # 1. 富集曲线图（前3个通路）
    try:
        start = time.time()
        gene_sets = load_gene_sets_from_gmt(TEST_DATA_DIR / "test_pathways_from_gmt.gmt")

        for i, row in gsea_results.head(3).iterrows():
            pathway_name = row['pathway']
            # 尝试匹配通路名称
            gene_set = set()
            for gs_name, gs_genes in gene_sets.items():
                if pathway_name.replace('_', ' ').lower() in gs_name.lower() or \
                   gs_name.lower() in pathway_name.replace('_', ' ').lower():
                    gene_set = gs_genes
                    break

            gene_weights = dict(zip(ranked_genes_df['gene'], ranked_genes_df['weight']))

            fig = plot_gsea_enrichment(
                ranked_genes=ranked_genes_df['gene'].tolist(),
                gene_weights=gene_weights,
                gene_set=gene_set,
                es=row['es'],
                nes=row['nes'],
                pvalue=row['pvalue'],
                title=f"{row['pathway']} (NES={row['nes']:.2f})",
                output_file=str(OUTPUT_DIR / f"gsea_enrichment_{i}.png")
            )
        tests.append({"name": "GSEA Enrichment Plot", "status": "passed", "time": time.time()-start})
        print("✓ GSEA富集曲线图")
    except Exception as e:
        tests.append({"name": "GSEA Enrichment Plot", "status": "failed", "error": str(e)})
        print(f"✗ GSEA富集曲线图: {e}")

    # 2. NES条形图
    try:
        start = time.time()
        fig = plot_gsea_nes_barplot(gsea_results, output_file=str(OUTPUT_DIR / "gsea_nes_barplot.png"))
        tests.append({"name": "GSEA NES Barplot", "status": "passed", "time": time.time()-start})
        print("✓ GSEA NES条形图")
    except Exception as e:
        tests.append({"name": "GSEA NES Barplot", "status": "failed", "error": str(e)})
        print(f"✗ GSEA NES条形图: {e}")

    # 3. 气泡图
    try:
        start = time.time()
        fig = plot_gsea_dotplot(gsea_results, output_file=str(OUTPUT_DIR / "gsea_dotplot.png"))
        tests.append({"name": "GSEA Dotplot", "status": "passed", "time": time.time()-start})
        print("✓ GSEA气泡图")
    except Exception as e:
        tests.append({"name": "GSEA Dotplot", "status": "failed", "error": str(e)})
        print(f"✗ GSEA气泡图: {e}")

    return tests


def test_gsva_visualizations():
    """测试GSVA/ssGSEA可视化"""
    print("\n" + "=" * 60)
    print("测试GSVA/ssGSEA可视化")
    print("=" * 60)

    tests = []

    # 加载结果
    ssgsea_results = pd.read_csv(RESULTS_DIR / "ssgsea_results.csv", index_col=0)
    gsva_results = pd.read_csv(TEST_DATA_DIR / "gsva_results.csv", index_col=0)

    # 1. 热图
    try:
        start = time.time()
        fig = plot_pathway_heatmap(
            ssgsea_results,
            output_file=str(OUTPUT_DIR / "ssgsea_heatmap.png"),
            figsize=(10, 8),
            dpi=150
        )
        tests.append({"name": "ssGSEA Heatmap", "status": "passed", "time": time.time()-start})
        print("✓ ssGSEA热图")
    except Exception as e:
        tests.append({"name": "ssGSEA Heatmap", "status": "failed", "error": str(e)})
        print(f"✗ ssGSEA热图: {e}")

    # 2. 样本相关性热图
    try:
        start = time.time()
        fig = plot_sample_correlation(
            ssgsea_results,
            output_file=str(OUTPUT_DIR / "ssgsea_correlation.png"),
            dpi=150
        )
        tests.append({"name": "Sample Correlation", "status": "passed", "time": time.time()-start})
        print("✓ 样本相关性热图")
    except Exception as e:
        tests.append({"name": "Sample Correlation", "status": "failed", "error": str(e)})
        print(f"✗ 样本相关性热图: {e}")

    # 3. GSVA热图
    try:
        start = time.time()
        fig = plot_pathway_heatmap(
            gsva_results,
            output_file=str(OUTPUT_DIR / "gsva_heatmap.png"),
            figsize=(10, 8),
            dpi=150
        )
        tests.append({"name": "GSVA Heatmap", "status": "passed", "time": time.time()-start})
        print("✓ GSVA热图")
    except Exception as e:
        tests.append({"name": "GSVA Heatmap", "status": "failed", "error": str(e)})
        print(f"✗ GSVA热图: {e}")

    # 4. 组间比较图
    try:
        start = time.time()
        # 创建模拟分组
        groups = {
            "Group_A": list(ssgsea_results.columns[:3]),
            "Group_B": list(ssgsea_results.columns[3:])
        }
        fig = plot_group_comparison(
            ssgsea_results,
            groups=groups,
            plot_type="box",
            output_file=str(OUTPUT_DIR / "group_comparison.png"),
            dpi=150
        )
        tests.append({"name": "Group Comparison", "status": "passed", "time": time.time()-start})
        print("✓ 组间比较图")
    except Exception as e:
        tests.append({"name": "Group Comparison", "status": "failed", "error": str(e)})
        print(f"✗ 组间比较图: {e}")

    return tests


def test_common_visualizations():
    """测试通用可视化"""
    print("\n" + "=" * 60)
    print("测试通用可视化")
    print("=" * 60)

    tests = []

    # 加载基因集
    gene_sets = load_gene_sets_from_gmt(TEST_DATA_DIR / "test_pathways_from_gmt.gmt")

    # 1. 网络图
    try:
        start = time.time()
        gsea_results = pd.read_csv(RESULTS_DIR / "gsea_results.csv")
        fig = plot_enrichment_network(
            gene_sets,
            results_df=gsea_results,
            output_file=str(OUTPUT_DIR / "enrichment_network.png"),
            dpi=150
        )
        tests.append({"name": "Enrichment Network", "status": "passed", "time": time.time()-start})
        print("✓ 通路网络图")
    except Exception as e:
        tests.append({"name": "Enrichment Network", "status": "failed", "error": str(e)})
        print(f"✗ 通路网络图: {e}")

    # 2. 火山图
    try:
        gsea_results = pd.read_csv(RESULTS_DIR / "gsea_results.csv")
        start = time.time()
        fig = plot_volcano(
            gsea_results,
            nes_col="nes",
            pvalue_col="pvalue",
            output_file=str(OUTPUT_DIR / "gsea_volcano.png"),
            dpi=150
        )
        tests.append({"name": "Volcano Plot", "status": "passed", "time": time.time()-start})
        print("✓ 火山图")
    except Exception as e:
        tests.append({"name": "Volcano Plot", "status": "failed", "error": str(e)})
        print(f"✗ 火山图: {e}")

    # 3. 方法比较图
    try:
        start = time.time()
        # 加载GSVA和ssGSEA结果进行比较
        gsva_results = pd.read_csv(TEST_DATA_DIR / "gsva_results.csv", index_col=0)
        ssgsea_results = pd.read_csv(RESULTS_DIR / "ssgsea_results.csv", index_col=0)

        # 计算每个通路的均值进行比较
        gsva_mean = gsva_results.mean(axis=1)
        ssgsea_mean = ssgsea_results.mean(axis=1)

        fig = plot_method_comparison(
            gsva_mean,
            ssgsea_mean,
            method_a_name="GSVA",
            method_b_name="ssGSEA",
            output_file=str(OUTPUT_DIR / "method_comparison.png"),
            dpi=150
        )
        tests.append({"name": "Method Comparison", "status": "passed", "time": time.time()-start})
        print("✓ 方法比较图")
    except Exception as e:
        tests.append({"name": "Method Comparison", "status": "failed", "error": str(e)})
        print(f"✗ 方法比较图: {e}")

    return tests


def main():
    """主函数"""
    print("=" * 60)
    print("可视化集成端到端测试")
    print("=" * 60)

    all_tests = []

    # 测试GSEA可视化
    all_tests.extend(test_gsea_visualizations())

    # 测试GSVA可视化
    all_tests.extend(test_gsva_visualizations())

    # 测试通用可视化
    all_tests.extend(test_common_visualizations())

    # 生成报告
    report = {
        "test_name": "Visualization Integration E2E Test",
        "timestamp": time.strftime("%Y-%m-%d %H:%M:%S"),
        "tests": all_tests,
        "summary": {
            "total": len(all_tests),
            "passed": sum(1 for t in all_tests if t["status"] == "passed"),
            "failed": sum(1 for t in all_tests if t["status"] == "failed"),
            "total_time": sum(t.get("time", 0) for t in all_tests)
        },
        "output_dir": str(OUTPUT_DIR),
        "status": "passed" if all(t["status"] == "passed" for t in all_tests) else "failed"
    }

    with open(RESULTS_DIR / "visualization_report.json", 'w') as f:
        json.dump(report, f, indent=2)

    print("\n" + "=" * 60)
    print("可视化测试完成!")
    print(f"总计: {report['summary']['total']} 个测试")
    print(f"通过: {report['summary']['passed']} 个")
    print(f"失败: {report['summary']['failed']} 个")
    print(f"总时间: {report['summary']['total_time']:.2f}s")
    print("=" * 60)

    return report


if __name__ == "__main__":
    main()
