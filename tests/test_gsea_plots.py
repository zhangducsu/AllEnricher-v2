"""
GSEA 可视化模块单元测试
========================

测试覆盖范围：
- plot_gsea_enrichment 生成正确的 Figure 对象
- plot_gsea_nes_barplot 生成正确的 Figure 对象
- plot_gsea_dotplot 生成正确的 Figure 对象
- 使用 E2E 测试数据（test_data/ranked_genes.tsv 和 test_data/gene_sets.gmt）
- 验证图表保存到临时文件
"""

import os
import pytest
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from pathlib import Path

# 使用非交互后端，避免弹窗
matplotlib.use("Agg")

from allenricher.visualization.gsea_plots import (
    plot_gsea_enrichment,
    plot_gsea_nes_barplot,
    plot_gsea_dotplot,
)
from allenricher.visualization.plot_config import PlotConfig


# ============================================================
# 测试数据加载工具
# ============================================================

TEST_DATA_DIR = Path(__file__).parent.parent / "test_data"


def _load_ranked_genes() -> tuple:
    """
    从 test_data/ranked_genes.tsv 加载排序基因列表和权重

    Returns:
        (ranked_genes: list[str], gene_weights: dict[str, float])
    """
    tsv_path = TEST_DATA_DIR / "ranked_genes.tsv"
    df = pd.read_csv(tsv_path, sep="\t")
    ranked_genes = df["gene"].tolist()
    gene_weights = dict(zip(df["gene"], df["weight"]))
    return ranked_genes, gene_weights


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


def _make_gsea_results_df(n_pathways: int = 10, seed: int = 42) -> pd.DataFrame:
    """
    构造 GSEA 结果 DataFrame（用于条形图和气泡图测试）

    Returns:
        包含 pathway, nes, pvalue, gene_count 列的 DataFrame
    """
    rng = np.random.default_rng(seed)
    pathways = [f"Pathway_{i}" for i in range(n_pathways)]
    nes = rng.uniform(-3.0, 3.0, size=n_pathways)
    pvalue = rng.uniform(1e-10, 0.05, size=n_pathways)
    gene_count = rng.integers(10, 100, size=n_pathways)

    return pd.DataFrame({
        "pathway": pathways,
        "nes": nes,
        "pvalue": pvalue,
        "gene_count": gene_count,
    })


# ============================================================
# 测试类
# ============================================================

class TestPlotGseaEnrichment:
    """测试 plot_gsea_enrichment"""

    def test_returns_figure_object(self):
        """测试返回 matplotlib Figure 对象"""
        ranked_genes, gene_weights = _load_ranked_genes()
        gene_sets = _load_gene_sets()
        # 取第一个基因集
        pathway_name = list(gene_sets.keys())[0]
        gene_set = gene_sets[pathway_name]

        fig = plot_gsea_enrichment(
            ranked_genes=ranked_genes,
            gene_weights=gene_weights,
            gene_set=gene_set,
            es=0.65,
            nes=1.82,
            pvalue=0.003,
        )

        assert isinstance(fig, matplotlib.figure.Figure)
        plt.close(fig)

    def test_figure_has_three_axes(self):
        """测试图表包含三个子图（三面板）"""
        ranked_genes, gene_weights = _load_ranked_genes()
        gene_sets = _load_gene_sets()
        pathway_name = list(gene_sets.keys())[0]
        gene_set = gene_sets[pathway_name]

        fig = plot_gsea_enrichment(
            ranked_genes=ranked_genes,
            gene_weights=gene_weights,
            gene_set=gene_set,
            es=0.65,
            nes=1.82,
            pvalue=0.003,
        )

        assert len(fig.axes) == 3
        plt.close(fig)

    def test_custom_title(self):
        """测试自定义标题"""
        ranked_genes, gene_weights = _load_ranked_genes()
        gene_sets = _load_gene_sets()
        pathway_name = list(gene_sets.keys())[0]
        gene_set = gene_sets[pathway_name]

        fig = plot_gsea_enrichment(
            ranked_genes=ranked_genes,
            gene_weights=gene_weights,
            gene_set=gene_set,
            es=0.65,
            nes=1.82,
            pvalue=0.003,
            title="Custom Title Test",
        )

        # suptitle 应包含自定义标题
        assert fig._suptitle is not None
        assert "Custom Title Test" in fig._suptitle.get_text()
        plt.close(fig)

    def test_save_to_file(self, tmp_path):
        """测试图表保存到临时文件"""
        ranked_genes, gene_weights = _load_ranked_genes()
        gene_sets = _load_gene_sets()
        pathway_name = list(gene_sets.keys())[0]
        gene_set = gene_sets[pathway_name]

        output_file = str(tmp_path / "gsea_enrichment.png")
        fig = plot_gsea_enrichment(
            ranked_genes=ranked_genes,
            gene_weights=gene_weights,
            gene_set=gene_set,
            es=0.65,
            nes=1.82,
            pvalue=0.003,
            output_file=output_file,
        )

        assert os.path.exists(output_file)
        assert os.path.getsize(output_file) > 0
        plt.close(fig)

    def test_empty_gene_set(self):
        """测试空基因集不报错"""
        ranked_genes, gene_weights = _load_ranked_genes()

        fig = plot_gsea_enrichment(
            ranked_genes=ranked_genes[:50],
            gene_weights=gene_weights,
            gene_set=set(),
            es=0.0,
            nes=0.0,
            pvalue=1.0,
        )

        assert isinstance(fig, matplotlib.figure.Figure)
        plt.close(fig)

    def test_all_gene_sets(self):
        """测试所有 E2E 基因集均可正常绘图"""
        ranked_genes, gene_weights = _load_ranked_genes()
        gene_sets = _load_gene_sets()

        for pathway_name, gene_set in gene_sets.items():
            fig = plot_gsea_enrichment(
                ranked_genes=ranked_genes,
                gene_weights=gene_weights,
                gene_set=gene_set,
                es=0.5,
                nes=1.5,
                pvalue=0.01,
            )
            assert isinstance(fig, matplotlib.figure.Figure)
            plt.close(fig)


class TestPlotGseaNesBarplot:
    """测试 plot_gsea_nes_barplot"""

    def test_returns_figure_object(self):
        """测试返回 matplotlib Figure 对象"""
        df = _make_gsea_results_df()

        fig = plot_gsea_nes_barplot(df)

        assert isinstance(fig, matplotlib.figure.Figure)
        plt.close(fig)

    def test_top_n_filtering(self):
        """测试 top_n 参数过滤"""
        df = _make_gsea_results_df(n_pathways=30)

        fig = plot_gsea_nes_barplot(df, top_n=10)
        ax = fig.axes[0]

        # Y 轴标签数量应 <= top_n
        ytick_labels = [t.get_text() for t in ax.get_yticklabels()]
        assert len(ytick_labels) <= 10
        plt.close(fig)

    def test_save_to_file(self, tmp_path):
        """测试图表保存到临时文件"""
        df = _make_gsea_results_df()
        output_file = str(tmp_path / "nes_barplot.png")

        fig = plot_gsea_nes_barplot(df, output_file=output_file)

        assert os.path.exists(output_file)
        assert os.path.getsize(output_file) > 0
        plt.close(fig)

    def test_custom_title(self):
        """测试自定义标题"""
        df = _make_gsea_results_df()

        fig = plot_gsea_nes_barplot(df, title="Custom NES Title")
        ax = fig.axes[0]

        assert ax.get_title() == "Custom NES Title"
        plt.close(fig)

    def test_with_real_gsea_data(self):
        """测试使用真实 GSEA 结果格式数据"""
        df = pd.DataFrame({
            "pathway": ["Cell_Cycle", "Apoptosis", "Immune_Response"],
            "nes": [2.5, -1.8, 1.2],
            "pvalue": [0.001, 0.01, 0.03],
            "gene_count": [45, 30, 25],
        })

        fig = plot_gsea_nes_barplot(df)

        assert isinstance(fig, matplotlib.figure.Figure)
        plt.close(fig)


class TestPlotGseaDotplot:
    """测试 plot_gsea_dotplot"""

    def test_returns_figure_object(self):
        """测试返回 matplotlib Figure 对象"""
        df = _make_gsea_results_df()

        fig = plot_gsea_dotplot(df)

        assert isinstance(fig, matplotlib.figure.Figure)
        plt.close(fig)

    def test_top_n_filtering(self):
        """测试 top_n 参数过滤"""
        df = _make_gsea_results_df(n_pathways=30)

        fig = plot_gsea_dotplot(df, top_n=10)
        ax = fig.axes[0]

        ytick_labels = [t.get_text() for t in ax.get_yticklabels()]
        assert len(ytick_labels) <= 10
        plt.close(fig)

    def test_save_to_file(self, tmp_path):
        """测试图表保存到临时文件"""
        df = _make_gsea_results_df()
        output_file = str(tmp_path / "gsea_dotplot.png")

        fig = plot_gsea_dotplot(df, output_file=output_file)

        assert os.path.exists(output_file)
        assert os.path.getsize(output_file) > 0
        plt.close(fig)

    def test_custom_title(self):
        """测试自定义标题"""
        df = _make_gsea_results_df()

        fig = plot_gsea_dotplot(df, title="Custom Dot Plot Title")
        ax = fig.axes[0]

        assert ax.get_title() == "Custom Dot Plot Title"
        plt.close(fig)

    def test_has_colorbar(self):
        """测试图表包含颜色条"""
        df = _make_gsea_results_df()

        fig = plot_gsea_dotplot(df)

        # 查找 colorbar（通过 ScalarMappable 检查）
        # fig.axes 包含主轴和 colorbar 轴
        assert len(fig.axes) >= 2
        plt.close(fig)

    def test_with_real_gsea_data(self):
        """测试使用真实 GSEA 结果格式数据"""
        df = pd.DataFrame({
            "pathway": ["Cell_Cycle", "Apoptosis", "Immune_Response"],
            "nes": [2.5, -1.8, 1.2],
            "pvalue": [0.001, 0.01, 0.03],
            "gene_count": [45, 30, 25],
        })

        fig = plot_gsea_dotplot(df)

        assert isinstance(fig, matplotlib.figure.Figure)
        plt.close(fig)


class TestPlotConfig:
    """测试 PlotConfig 配置类"""

    def test_default_values(self):
        """测试默认配置值"""
        config = PlotConfig()
        assert config.figure_format == "png"
        assert config.figure_dpi == 300
        assert config.color_palette == "RdBu_r"
        assert config.font_family == "sans-serif"
        assert config.font_size == 10
        assert config.top_n_pathways == 20

    def test_custom_values(self):
        """测试自定义配置值"""
        config = PlotConfig(
            figure_format="pdf",
            figure_dpi=600,
            font_size=12,
        )
        assert config.figure_format == "pdf"
        assert config.figure_dpi == 600
        assert config.font_size == 12


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
