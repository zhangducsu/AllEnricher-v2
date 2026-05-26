"""
可视化集成端到端测试

测试所有可视化函数返回正确的Figure对象
"""

import pytest
import pandas as pd
import numpy as np
from pathlib import Path
import matplotlib.pyplot as plt
import matplotlib.figure

# 导入可视化模块
from allenricher.visualization.gsea_plots import (
    plot_gsea_enrichment, plot_gsea_nes_barplot, plot_gsea_dotplot
)
from allenricher.visualization.gsva_plots import (
    plot_pathway_heatmap, plot_group_comparison, plot_sample_correlation,
    plot_pathway_dotplot
)
from allenricher.visualization.common_plots import (
    plot_enrichment_network, plot_volcano, plot_method_comparison, plot_upset
)


# 测试数据路径
TEST_DATA_DIR = Path(__file__).parent.parent / "test_data"
RESULTS_DIR = TEST_DATA_DIR / "e2e_results"
OUTPUT_DIR = RESULTS_DIR / "plots"


@pytest.fixture
def gsea_results():
    """加载GSEA结果"""
    return pd.read_csv(RESULTS_DIR / "gsea_results.csv")


@pytest.fixture
def ranked_genes():
    """加载排序基因列表"""
    return pd.read_csv(TEST_DATA_DIR / "ranked_genes.tsv", sep='\t')


@pytest.fixture
def ssgsea_results():
    """加载ssGSEA结果"""
    return pd.read_csv(RESULTS_DIR / "ssgsea_results.csv", index_col=0)


@pytest.fixture
def gsva_results():
    """加载GSVA结果"""
    return pd.read_csv(TEST_DATA_DIR / "gsva_results.csv", index_col=0)


@pytest.fixture
def gene_sets():
    """加载基因集"""
    gene_sets = {}
    gmt_file = TEST_DATA_DIR / "test_pathways_from_gmt.gmt"
    with open(gmt_file, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 3:
                gene_sets[parts[0]] = set(parts[2:])
    return gene_sets


class TestGSEAPlots:
    """测试GSEA可视化函数"""

    def test_plot_gsea_enrichment_returns_figure(self, ranked_genes, gsea_results, gene_sets):
        """测试plot_gsea_enrichment返回Figure对象"""
        row = gsea_results.iloc[0]
        gene_set = list(gene_sets.values())[0] if gene_sets else set()
        gene_weights = dict(zip(ranked_genes['gene'], ranked_genes['weight']))

        fig = plot_gsea_enrichment(
            ranked_genes=ranked_genes['gene'].tolist()[:500],
            gene_weights=gene_weights,
            gene_set=gene_set,
            es=row['es'],
            nes=row['nes'],
            pvalue=row['pvalue'],
            title="Test Enrichment"
        )

        assert isinstance(fig, matplotlib.figure.Figure)
        plt.close(fig)

    def test_plot_gsea_nes_barplot_returns_figure(self, gsea_results):
        """测试plot_gsea_nes_barplot返回Figure对象"""
        fig = plot_gsea_nes_barplot(gsea_results, top_n=10)

        assert isinstance(fig, matplotlib.figure.Figure)
        plt.close(fig)

    def test_plot_gsea_dotplot_returns_figure(self, gsea_results):
        """测试plot_gsea_dotplot返回Figure对象"""
        fig = plot_gsea_dotplot(gsea_results, top_n=10)

        assert isinstance(fig, matplotlib.figure.Figure)
        plt.close(fig)

    def test_plot_gsea_nes_barplot_with_output(self, gsea_results, tmp_path):
        """测试plot_gsea_nes_barplot保存文件"""
        output_file = tmp_path / "test_barplot.png"
        fig = plot_gsea_nes_barplot(gsea_results, top_n=10, output_file=str(output_file))

        assert output_file.exists()
        plt.close(fig)

    def test_plot_gsea_dotplot_with_output(self, gsea_results, tmp_path):
        """测试plot_gsea_dotplot保存文件"""
        output_file = tmp_path / "test_dotplot.png"
        fig = plot_gsea_dotplot(gsea_results, top_n=10, output_file=str(output_file))

        assert output_file.exists()
        plt.close(fig)


class TestGSVAPlots:
    """测试GSVA/ssGSEA可视化函数"""

    def test_plot_pathway_heatmap_returns_figure(self, ssgsea_results):
        """测试plot_pathway_heatmap返回Figure对象"""
        fig = plot_pathway_heatmap(ssgsea_results, cluster_rows=False, cluster_cols=False)

        assert isinstance(fig, matplotlib.figure.Figure)
        plt.close(fig)

    def test_plot_sample_correlation_returns_figure(self, ssgsea_results):
        """测试plot_sample_correlation返回Figure对象"""
        fig = plot_sample_correlation(ssgsea_results, method="pearson")

        assert isinstance(fig, matplotlib.figure.Figure)
        plt.close(fig)

    def test_plot_group_comparison_returns_figure(self, ssgsea_results):
        """测试plot_group_comparison返回Figure对象"""
        groups = {
            "Group_A": list(ssgsea_results.columns[:3]),
            "Group_B": list(ssgsea_results.columns[3:])
        }
        fig = plot_group_comparison(
            ssgsea_results,
            groups=groups,
            pathways=list(ssgsea_results.index[:4]),
            plot_type="box"
        )

        assert isinstance(fig, matplotlib.figure.Figure)
        plt.close(fig)

    def test_plot_pathway_dotplot_returns_figure(self, ssgsea_results):
        """测试plot_pathway_dotplot返回Figure对象"""
        fig = plot_pathway_dotplot(ssgsea_results, top_n=10)

        assert isinstance(fig, matplotlib.figure.Figure)
        plt.close(fig)

    def test_plot_pathway_heatmap_with_output(self, ssgsea_results, tmp_path):
        """测试plot_pathway_heatmap保存文件"""
        output_file = tmp_path / "test_heatmap.png"
        fig = plot_pathway_heatmap(
            ssgsea_results,
            cluster_rows=False,
            cluster_cols=False,
            output_file=str(output_file),
            dpi=150
        )

        assert output_file.exists()
        plt.close(fig)

    def test_plot_sample_correlation_spearman(self, ssgsea_results):
        """测试plot_sample_correlation使用spearman方法"""
        fig = plot_sample_correlation(ssgsea_results, method="spearman")

        assert isinstance(fig, matplotlib.figure.Figure)
        plt.close(fig)

    def test_plot_group_comparison_violin(self, ssgsea_results):
        """测试plot_group_comparison使用violin类型"""
        groups = {
            "Group_A": list(ssgsea_results.columns[:3]),
            "Group_B": list(ssgsea_results.columns[3:])
        }
        fig = plot_group_comparison(
            ssgsea_results,
            groups=groups,
            pathways=list(ssgsea_results.index[:4]),
            plot_type="violin"
        )

        assert isinstance(fig, matplotlib.figure.Figure)
        plt.close(fig)

    def test_plot_group_comparison_bar(self, ssgsea_results):
        """测试plot_group_comparison使用bar类型"""
        groups = {
            "Group_A": list(ssgsea_results.columns[:3]),
            "Group_B": list(ssgsea_results.columns[3:])
        }
        fig = plot_group_comparison(
            ssgsea_results,
            groups=groups,
            pathways=list(ssgsea_results.index[:4]),
            plot_type="bar"
        )

        assert isinstance(fig, matplotlib.figure.Figure)
        plt.close(fig)


class TestCommonPlots:
    """测试通用可视化函数"""

    def test_plot_enrichment_network_returns_figure(self, gene_sets, gsea_results):
        """测试plot_enrichment_network返回Figure对象"""
        fig = plot_enrichment_network(
            gene_sets,
            results_df=gsea_results,
            top_n=10,
            min_overlap=1
        )

        assert isinstance(fig, matplotlib.figure.Figure)
        plt.close(fig)

    def test_plot_volcano_returns_figure(self, gsea_results):
        """测试plot_volcano返回Figure对象"""
        fig = plot_volcano(
            gsea_results,
            nes_col="nes",
            pvalue_col="pvalue"
        )

        assert isinstance(fig, matplotlib.figure.Figure)
        plt.close(fig)

    def test_plot_method_comparison_returns_figure(self, gsva_results, ssgsea_results):
        """测试plot_method_comparison返回Figure对象"""
        gsva_mean = gsva_results.mean(axis=1)
        ssgsea_mean = ssgsea_results.mean(axis=1)

        fig = plot_method_comparison(
            gsva_mean,
            ssgsea_mean,
            method_a_name="GSVA",
            method_b_name="ssGSEA"
        )

        assert isinstance(fig, matplotlib.figure.Figure)
        plt.close(fig)

    def test_plot_upset_returns_figure(self, gene_sets):
        """测试plot_upset返回Figure对象"""
        fig = plot_upset(gene_sets, top_n=10)

        assert isinstance(fig, matplotlib.figure.Figure)
        plt.close(fig)

    def test_plot_enrichment_network_without_results(self, gene_sets):
        """测试plot_enrichment_network不带results_df"""
        fig = plot_enrichment_network(gene_sets, top_n=10)

        assert isinstance(fig, matplotlib.figure.Figure)
        plt.close(fig)

    def test_plot_volcano_with_thresholds(self, gsea_results):
        """测试plot_volcano使用自定义阈值"""
        fig = plot_volcano(
            gsea_results,
            nes_col="nes",
            pvalue_col="pvalue",
            nes_threshold=0.5,
            pvalue_threshold=0.1
        )

        assert isinstance(fig, matplotlib.figure.Figure)
        plt.close(fig)

    def test_plot_volcano_with_output(self, gsea_results, tmp_path):
        """测试plot_volcano保存文件"""
        output_file = tmp_path / "test_volcano.png"
        fig = plot_volcano(
            gsea_results,
            nes_col="nes",
            pvalue_col="pvalue",
            output_file=str(output_file),
            dpi=150
        )

        assert output_file.exists()
        plt.close(fig)


class TestEdgeCases:
    """测试边界情况"""

    def test_empty_gsea_results(self):
        """测试空的GSEA结果"""
        empty_df = pd.DataFrame(columns=["pathway", "es", "nes", "pvalue", "gene_count"])
        fig = plot_gsea_nes_barplot(empty_df, top_n=10)

        assert isinstance(fig, matplotlib.figure.Figure)
        plt.close(fig)

    def test_single_pathway_gsea(self, gsea_results):
        """测试单个通路的GSEA结果"""
        single_df = gsea_results.head(1)
        fig = plot_gsea_nes_barplot(single_df, top_n=10)

        assert isinstance(fig, matplotlib.figure.Figure)
        plt.close(fig)

    def test_empty_gene_sets(self):
        """测试空的基因集"""
        empty_sets = {}
        fig = plot_enrichment_network(empty_sets, top_n=10)

        assert isinstance(fig, matplotlib.figure.Figure)
        plt.close(fig)

    def test_single_gene_set(self, gene_sets):
        """测试单个基因集"""
        if gene_sets:
            single_set = {list(gene_sets.keys())[0]: list(gene_sets.values())[0]}
            fig = plot_enrichment_network(single_set, top_n=10)

            assert isinstance(fig, matplotlib.figure.Figure)
            plt.close(fig)

    def test_plot_upset_empty(self):
        """测试空的UpSet图"""
        empty_sets = {}
        fig = plot_upset(empty_sets, top_n=10)

        assert isinstance(fig, matplotlib.figure.Figure)
        plt.close(fig)


class TestIntegration:
    """集成测试 - 测试完整工作流"""

    def test_full_visualization_pipeline(self, gsea_results, ssgsea_results, gsva_results, gene_sets, tmp_path):
        """测试完整的可视化管道"""
        output_dir = tmp_path / "viz_output"
        output_dir.mkdir(exist_ok=True)

        # 1. GSEA可视化 (注意：plot_gsea_nes_barplot和plot_gsea_dotplot不支持dpi参数)
        fig1 = plot_gsea_nes_barplot(gsea_results, output_file=str(output_dir / "gsea_bar.png"))
        fig2 = plot_gsea_dotplot(gsea_results, output_file=str(output_dir / "gsea_dot.png"))
        fig3 = plot_volcano(gsea_results, output_file=str(output_dir / "gsea_volcano.png"), dpi=150)

        # 2. GSVA可视化
        fig4 = plot_pathway_heatmap(
            ssgsea_results,
            cluster_rows=False,
            cluster_cols=False,
            output_file=str(output_dir / "ssgsea_heatmap.png"),
            dpi=150
        )
        fig5 = plot_sample_correlation(ssgsea_results, output_file=str(output_dir / "correlation.png"), dpi=150)

        # 3. 通用可视化
        fig6 = plot_enrichment_network(
            gene_sets,
            results_df=gsea_results,
            output_file=str(output_dir / "network.png"),
            dpi=150
        )

        # 验证所有文件都已创建
        expected_files = [
            "gsea_bar.png",
            "gsea_dot.png",
            "gsea_volcano.png",
            "ssgsea_heatmap.png",
            "correlation.png",
            "network.png"
        ]

        for filename in expected_files:
            assert (output_dir / filename).exists(), f"文件 {filename} 未创建"

        # 清理
        for fig in [fig1, fig2, fig3, fig4, fig5, fig6]:
            plt.close(fig)
