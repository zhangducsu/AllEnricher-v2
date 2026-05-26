"""
通用可视化组件单元测试
======================

测试覆盖范围：
- plot_enrichment_network: 通路关系网络图
- plot_upset: UpSet 图（基因集交集）
- plot_volcano: 火山图
- plot_method_comparison: 方法间比较散点图
- 边界情况：空数据、单通路、无交集等
- 使用 test_data/gene_sets.gmt 中的基因集数据
"""

import os
from pathlib import Path

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pytest

# 使用非交互后端，避免弹窗
matplotlib.use("Agg")

from allenricher.visualization.common_plots import (
    plot_enrichment_network,
    plot_method_comparison,
    plot_upset,
    plot_volcano,
)


# ============================================================
# 测试数据加载工具
# ============================================================

TEST_DATA_DIR = Path(__file__).parent.parent / "test_data"


def _load_gene_sets() -> dict:
    """
    从 test_data/gene_sets.gmt 加载基因集

    Returns:
        {pathway_name: set(genes)}
    """
    gmt_path = TEST_DATA_DIR / "gene_sets.gmt"
    gene_sets = {}
    with open(gmt_path, "r") as f:
        for line in f:
            parts = line.strip().split("\t")
            name = parts[0]
            genes = set(parts[2:])
            gene_sets[name] = genes
    return gene_sets


def _make_volcano_df(n_pathways: int = 50, seed: int = 42) -> pd.DataFrame:
    """
    构造火山图测试数据

    Returns:
        包含 pathway, nes, pvalue 列的 DataFrame
    """
    rng = np.random.default_rng(seed)
    pathways = [f"Pathway_{i}" for i in range(n_pathways)]
    nes = rng.normal(0, 1.5, size=n_pathways)
    # 让部分通路显著
    nes[:5] = rng.uniform(2.0, 4.0, size=5)
    nes[5:10] = rng.uniform(-4.0, -2.0, size=5)
    pvalue = rng.uniform(1e-10, 0.1, size=n_pathways)
    # 显著通路 pvalue 更小
    pvalue[:10] = rng.uniform(1e-10, 0.01, size=10)

    return pd.DataFrame({
        "pathway": pathways,
        "nes": nes,
        "pvalue": pvalue,
    })


def _make_method_results(n_pathways: int = 30, seed: int = 42) -> tuple:
    """
    构造两种方法的模拟结果

    Returns:
        (results_a: pd.Series, results_b: pd.Series)
    """
    rng = np.random.default_rng(seed)
    pathways = [f"Pathway_{i}" for i in range(n_pathways)]
    base = rng.normal(0, 1.0, size=n_pathways)
    noise_a = rng.normal(0, 0.3, size=n_pathways)
    noise_b = rng.normal(0, 0.3, size=n_pathways)

    results_a = pd.Series(base + noise_a, index=pathways, name="Method_A")
    results_b = pd.Series(base + noise_b, index=pathways, name="Method_B")
    return results_a, results_b


# ============================================================
# plot_enrichment_network 测试
# ============================================================

class TestPlotEnrichmentNetwork:
    """测试 plot_enrichment_network"""

    def test_returns_figure(self):
        """测试返回 matplotlib Figure 对象"""
        gene_sets = _load_gene_sets()
        fig = plot_enrichment_network(gene_sets)
        assert isinstance(fig, matplotlib.figure.Figure)
        plt.close(fig)

    def test_with_results_df(self):
        """测试带 results_df 参数"""
        gene_sets = _load_gene_sets()
        results_df = pd.DataFrame({
            "pathway": list(gene_sets.keys()),
            "nes": np.random.randn(len(gene_sets)),
            "pvalue": np.random.uniform(0.001, 0.05, len(gene_sets)),
        })
        fig = plot_enrichment_network(gene_sets, results_df=results_df)
        assert isinstance(fig, matplotlib.figure.Figure)
        plt.close(fig)

    def test_top_n_filtering(self):
        """测试 top_n 参数"""
        gene_sets = _load_gene_sets()
        fig = plot_enrichment_network(gene_sets, top_n=5)
        assert isinstance(fig, matplotlib.figure.Figure)
        plt.close(fig)

    def test_min_overlap(self):
        """测试 min_overlap 参数"""
        gene_sets = _load_gene_sets()
        fig = plot_enrichment_network(gene_sets, min_overlap=10)
        assert isinstance(fig, matplotlib.figure.Figure)
        plt.close(fig)

    def test_circular_layout(self):
        """测试 circular 布局"""
        gene_sets = _load_gene_sets()
        fig = plot_enrichment_network(gene_sets, layout="circular")
        assert isinstance(fig, matplotlib.figure.Figure)
        plt.close(fig)

    def test_kamada_kawai_layout(self):
        """测试 kamada_kawai 布局"""
        gene_sets = _load_gene_sets()
        fig = plot_enrichment_network(gene_sets, layout="kamada_kawai")
        assert isinstance(fig, matplotlib.figure.Figure)
        plt.close(fig)

    def test_save_to_file(self, tmp_path):
        """测试保存到临时文件"""
        gene_sets = _load_gene_sets()
        output_file = str(tmp_path / "network.png")
        fig = plot_enrichment_network(gene_sets, output_file=output_file)
        assert os.path.exists(output_file)
        assert os.path.getsize(output_file) > 0
        plt.close(fig)

    def test_single_pathway(self):
        """测试单通路边界情况"""
        gene_sets = {"Pathway_A": {"GENE1", "GENE2", "GENE3"}}
        fig = plot_enrichment_network(gene_sets)
        assert isinstance(fig, matplotlib.figure.Figure)
        plt.close(fig)

    def test_empty_gene_sets(self):
        """测试空基因集边界情况"""
        gene_sets = {}
        fig = plot_enrichment_network(gene_sets)
        assert isinstance(fig, matplotlib.figure.Figure)
        plt.close(fig)

    def test_custom_title(self):
        """测试自定义标题"""
        gene_sets = _load_gene_sets()
        fig = plot_enrichment_network(gene_sets, title="Custom Network Title")
        assert isinstance(fig, matplotlib.figure.Figure)
        plt.close(fig)


# ============================================================
# plot_upset 测试
# ============================================================

class TestPlotUpset:
    """测试 plot_upset"""

    def test_returns_figure(self):
        """测试返回 matplotlib Figure 对象"""
        gene_sets = _load_gene_sets()
        fig = plot_upset(gene_sets)
        assert isinstance(fig, matplotlib.figure.Figure)
        plt.close(fig)

    def test_top_n_filtering(self):
        """测试 top_n 参数"""
        gene_sets = _load_gene_sets()
        fig = plot_upset(gene_sets, top_n=5)
        assert isinstance(fig, matplotlib.figure.Figure)
        plt.close(fig)

    def test_save_to_file(self, tmp_path):
        """测试保存到临时文件"""
        gene_sets = _load_gene_sets()
        output_file = str(tmp_path / "upset.png")
        fig = plot_upset(gene_sets, output_file=output_file)
        assert os.path.exists(output_file)
        assert os.path.getsize(output_file) > 0
        plt.close(fig)

    def test_empty_gene_sets(self):
        """测试空基因集边界情况"""
        gene_sets = {}
        fig = plot_upset(gene_sets)
        assert isinstance(fig, matplotlib.figure.Figure)
        plt.close(fig)

    def test_single_gene_set(self):
        """测试单基因集边界情况"""
        gene_sets = {"Pathway_A": {"GENE1", "GENE2", "GENE3"}}
        fig = plot_upset(gene_sets)
        assert isinstance(fig, matplotlib.figure.Figure)
        plt.close(fig)

    def test_no_overlapping_genes(self):
        """测试无重叠基因的情况"""
        gene_sets = {
            "Set_A": {"GENE1", "GENE2", "GENE3"},
            "Set_B": {"GENE4", "GENE5", "GENE6"},
            "Set_C": {"GENE7", "GENE8", "GENE9"},
        }
        fig = plot_upset(gene_sets)
        assert isinstance(fig, matplotlib.figure.Figure)
        plt.close(fig)

    def test_custom_title(self):
        """测试自定义标题"""
        gene_sets = _load_gene_sets()
        fig = plot_upset(gene_sets, title="Custom Upset Title")
        assert isinstance(fig, matplotlib.figure.Figure)
        plt.close(fig)

    def test_figure_has_two_axes(self):
        """测试图表包含两个子图（条形图 + 矩阵）"""
        gene_sets = _load_gene_sets()
        fig = plot_upset(gene_sets)
        assert len(fig.axes) == 2
        plt.close(fig)


# ============================================================
# plot_volcano 测试
# ============================================================

class TestPlotVolcano:
    """测试 plot_volcano"""

    def test_returns_figure(self):
        """测试返回 matplotlib Figure 对象"""
        df = _make_volcano_df()
        fig = plot_volcano(df)
        assert isinstance(fig, matplotlib.figure.Figure)
        plt.close(fig)

    def test_custom_columns(self):
        """测试自定义列名"""
        df = _make_volcano_df()
        df = df.rename(columns={"nes": "log2fc", "pvalue": "padj"})
        fig = plot_volcano(df, nes_col="log2fc", pvalue_col="padj")
        assert isinstance(fig, matplotlib.figure.Figure)
        plt.close(fig)

    def test_custom_thresholds(self):
        """测试自定义阈值"""
        df = _make_volcano_df()
        fig = plot_volcano(df, nes_threshold=2.0, pvalue_threshold=0.01)
        assert isinstance(fig, matplotlib.figure.Figure)
        plt.close(fig)

    def test_save_to_file(self, tmp_path):
        """测试保存到临时文件"""
        df = _make_volcano_df()
        output_file = str(tmp_path / "volcano.png")
        fig = plot_volcano(df, output_file=output_file)
        assert os.path.exists(output_file)
        assert os.path.getsize(output_file) > 0
        plt.close(fig)

    def test_empty_dataframe(self):
        """测试空 DataFrame 边界情况"""
        df = pd.DataFrame(columns=["pathway", "nes", "pvalue"])
        fig = plot_volcano(df)
        assert isinstance(fig, matplotlib.figure.Figure)
        plt.close(fig)

    def test_single_pathway(self):
        """测试单通路边界情况"""
        df = pd.DataFrame({
            "pathway": ["Only_Pathway"],
            "nes": [2.5],
            "pvalue": [0.001],
        })
        fig = plot_volcano(df)
        assert isinstance(fig, matplotlib.figure.Figure)
        plt.close(fig)

    def test_all_not_significant(self):
        """测试所有通路均不显著"""
        rng = np.random.default_rng(42)
        df = pd.DataFrame({
            "pathway": [f"PW_{i}" for i in range(20)],
            "nes": rng.normal(0, 0.5, size=20),
            "pvalue": rng.uniform(0.1, 1.0, size=20),
        })
        fig = plot_volcano(df)
        assert isinstance(fig, matplotlib.figure.Figure)
        plt.close(fig)

    def test_custom_title(self):
        """测试自定义标题"""
        df = _make_volcano_df()
        fig = plot_volcano(df, title="Custom Volcano Title")
        assert isinstance(fig, matplotlib.figure.Figure)
        plt.close(fig)

    def test_has_legend(self):
        """测试图表包含图例"""
        df = _make_volcano_df()
        fig = plot_volcano(df)
        ax = fig.axes[0]
        legend = ax.get_legend()
        assert legend is not None
        plt.close(fig)


# ============================================================
# plot_method_comparison 测试
# ============================================================

class TestPlotMethodComparison:
    """测试 plot_method_comparison"""

    def test_returns_figure(self):
        """测试返回 matplotlib Figure 对象"""
        results_a, results_b = _make_method_results()
        fig = plot_method_comparison(results_a, results_b)
        assert isinstance(fig, matplotlib.figure.Figure)
        plt.close(fig)

    def test_custom_names(self):
        """测试自定义方法名称"""
        results_a, results_b = _make_method_results()
        fig = plot_method_comparison(
            results_a, results_b,
            method_a_name="GSEA",
            method_b_name="ssGSEA",
        )
        assert isinstance(fig, matplotlib.figure.Figure)
        plt.close(fig)

    def test_save_to_file(self, tmp_path):
        """测试保存到临时文件"""
        results_a, results_b = _make_method_results()
        output_file = str(tmp_path / "comparison.png")
        fig = plot_method_comparison(
            results_a, results_b, output_file=output_file
        )
        assert os.path.exists(output_file)
        assert os.path.getsize(output_file) > 0
        plt.close(fig)

    def test_no_common_pathways(self):
        """测试无共有通路边界情况"""
        results_a = pd.Series({"A": 1.0, "B": 2.0}, name="Method_A")
        results_b = pd.Series({"C": 1.5, "D": 2.5}, name="Method_B")
        fig = plot_method_comparison(results_a, results_b)
        assert isinstance(fig, matplotlib.figure.Figure)
        plt.close(fig)

    def test_single_common_pathway(self):
        """测试仅一个共有通路"""
        results_a = pd.Series({"A": 1.0, "B": 2.0}, name="Method_A")
        results_b = pd.Series({"A": 1.5, "C": 2.5}, name="Method_B")
        fig = plot_method_comparison(results_a, results_b)
        assert isinstance(fig, matplotlib.figure.Figure)
        plt.close(fig)

    def test_custom_title(self):
        """测试自定义标题"""
        results_a, results_b = _make_method_results()
        fig = plot_method_comparison(
            results_a, results_b, title="Custom Comparison Title"
        )
        assert isinstance(fig, matplotlib.figure.Figure)
        plt.close(fig)

    def test_pearson_annotation(self):
        """测试图表包含 Pearson 相关系数标注"""
        results_a, results_b = _make_method_results()
        fig = plot_method_comparison(results_a, results_b)
        ax = fig.axes[0]
        # 检查文本标注存在（包含 "Pearson"）
        texts = [t.get_text() for t in ax.texts]
        assert any("Pearson" in t for t in texts)
        plt.close(fig)

    def test_with_real_gene_set_data(self):
        """测试使用真实基因集数据构造的方法结果"""
        gene_sets = _load_gene_sets()
        rng = np.random.default_rng(42)
        pathways = list(gene_sets.keys())
        results_a = pd.Series(rng.normal(0, 1, len(pathways)), index=pathways)
        results_b = pd.Series(rng.normal(0, 1, len(pathways)), index=pathways)

        fig = plot_method_comparison(results_a, results_b)
        assert isinstance(fig, matplotlib.figure.Figure)
        plt.close(fig)


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
