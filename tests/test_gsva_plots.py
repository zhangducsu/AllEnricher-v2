"""
ssGSEA/GSVA 可视化模块单元测试

测试覆盖范围：
- plot_pathway_heatmap: 生成正确的 Figure 对象，支持分组注释
- plot_group_comparison: 三种 plot_type（box/violin/bar）均生成正确 Figure
- plot_pathway_dotplot: 生成正确的 Figure 对象，支持分组
- plot_sample_correlation: 生成正确的 Figure 对象，支持 pearson/spearman
- 图表保存到临时文件验证
"""

import os
from pathlib import Path

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pytest

from allenricher.visualization.gsva_plots import (
    plot_group_comparison,
    plot_pathway_dotplot,
    plot_pathway_heatmap,
    plot_sample_correlation,
)

# 使用非交互式后端，避免测试时弹出窗口
matplotlib.use("Agg")


class TestPlotPathwayHeatmap:
    """plot_pathway_heatmap 测试"""

    @pytest.fixture
    def scores_df(self):
        """创建模拟通路活性得分矩阵（5 通路 x 6 样本）"""
        np.random.seed(42)
        pathways = [
            "HSA_Cell_Cycle",
            "HSA_DNA_Repair",
            "HSA_PI3K_AKT",
            "HSA_MAPK_Signaling",
            "HSA_Apoptosis",
        ]
        samples = [f"Normal_{i}" for i in range(1, 4)] + [f"Disease_{i}" for i in range(1, 4)]
        data = np.random.randn(len(pathways), len(samples))
        return pd.DataFrame(data, index=pathways, columns=samples)

    @pytest.fixture
    def annotation_col(self):
        """创建样本分组注释"""
        samples = [f"Normal_{i}" for i in range(1, 4)] + [f"Disease_{i}" for i in range(1, 4)]
        return pd.DataFrame({
            "Group": ["Normal"] * 3 + ["Disease"] * 3,
        }, index=samples)

    def test_returns_figure(self, scores_df):
        """测试返回 matplotlib Figure 对象"""
        fig = plot_pathway_heatmap(scores_df)
        assert fig is not None
        assert isinstance(fig, matplotlib.figure.Figure)
        plt.close(fig)

    def test_with_annotation(self, scores_df, annotation_col):
        """测试带分组注释的热图"""
        fig = plot_pathway_heatmap(scores_df, annotation_col=annotation_col)
        assert fig is not None
        assert isinstance(fig, matplotlib.figure.Figure)
        plt.close(fig)

    def test_no_clustering(self, scores_df):
        """测试不聚类"""
        fig = plot_pathway_heatmap(scores_df, cluster_rows=False, cluster_cols=False)
        assert fig is not None
        plt.close(fig)

    def test_save_to_file(self, scores_df, tmp_path):
        """测试保存到文件"""
        output_file = str(tmp_path / "heatmap.png")
        fig = plot_pathway_heatmap(scores_df, output_file=output_file)
        assert os.path.exists(output_file)
        assert os.path.getsize(output_file) > 0
        plt.close(fig)

    def test_save_pdf(self, scores_df, tmp_path):
        """测试保存为 PDF"""
        output_file = str(tmp_path / "heatmap.pdf")
        fig = plot_pathway_heatmap(scores_df, output_file=output_file)
        assert os.path.exists(output_file)
        assert os.path.getsize(output_file) > 0
        plt.close(fig)

    def test_custom_figsize(self, scores_df):
        """测试自定义 figsize"""
        fig = plot_pathway_heatmap(scores_df, figsize=(12, 10))
        assert fig is not None
        plt.close(fig)


class TestPlotGroupComparison:
    """plot_group_comparison 测试"""

    @pytest.fixture
    def scores_df(self):
        """创建模拟通路活性得分矩阵"""
        np.random.seed(42)
        pathways = [f"Pathway_{i}" for i in range(1, 16)]
        samples = [f"Normal_{i}" for i in range(1, 4)] + [f"Disease_{i}" for i in range(1, 4)]
        data = np.random.randn(len(pathways), len(samples))
        return pd.DataFrame(data, index=pathways, columns=samples)

    @pytest.fixture
    def groups(self):
        """创建分组字典"""
        return {
            "Normal": [f"Normal_{i}" for i in range(1, 4)],
            "Disease": [f"Disease_{i}" for i in range(1, 4)],
        }

    def test_box_plot(self, scores_df, groups):
        """测试箱线图类型"""
        fig = plot_group_comparison(scores_df, groups, plot_type="box")
        assert fig is not None
        assert isinstance(fig, matplotlib.figure.Figure)
        plt.close(fig)

    def test_violin_plot(self, scores_df, groups):
        """测试小提琴图类型"""
        fig = plot_group_comparison(scores_df, groups, plot_type="violin")
        assert fig is not None
        assert isinstance(fig, matplotlib.figure.Figure)
        plt.close(fig)

    def test_bar_plot(self, scores_df, groups):
        """测试柱状图类型"""
        fig = plot_group_comparison(scores_df, groups, plot_type="bar")
        assert fig is not None
        assert isinstance(fig, matplotlib.figure.Figure)
        plt.close(fig)

    def test_specific_pathways(self, scores_df, groups):
        """测试指定通路列表"""
        pathways = ["Pathway_1", "Pathway_2", "Pathway_3"]
        fig = plot_group_comparison(scores_df, groups, pathways=pathways)
        assert fig is not None
        plt.close(fig)

    def test_invalid_plot_type(self, scores_df, groups):
        """测试无效 plot_type 应抛出 ValueError"""
        with pytest.raises(ValueError, match="不支持的 plot_type"):
            plot_group_comparison(scores_df, groups, plot_type="scatter")

    def test_save_to_file(self, scores_df, groups, tmp_path):
        """测试保存到文件"""
        output_file = str(tmp_path / "group_comparison.png")
        fig = plot_group_comparison(scores_df, groups, output_file=output_file)
        assert os.path.exists(output_file)
        assert os.path.getsize(output_file) > 0
        plt.close(fig)


class TestPlotPathwayDotplot:
    """plot_pathway_dotplot 测试"""

    @pytest.fixture
    def scores_df(self):
        """创建模拟通路活性得分矩阵"""
        np.random.seed(42)
        pathways = [f"Pathway_{i}" for i in range(1, 26)]
        samples = [f"Normal_{i}" for i in range(1, 4)] + [f"Disease_{i}" for i in range(1, 4)]
        data = np.random.randn(len(pathways), len(samples))
        return pd.DataFrame(data, index=pathways, columns=samples)

    @pytest.fixture
    def groups(self):
        """创建分组字典"""
        return {
            "Normal": [f"Normal_{i}" for i in range(1, 4)],
            "Disease": [f"Disease_{i}" for i in range(1, 4)],
        }

    def test_returns_figure(self, scores_df):
        """测试返回 Figure 对象（无分组）"""
        fig = plot_pathway_dotplot(scores_df)
        assert fig is not None
        assert isinstance(fig, matplotlib.figure.Figure)
        plt.close(fig)

    def test_with_groups(self, scores_df, groups):
        """测试带分组的气泡图"""
        fig = plot_pathway_dotplot(scores_df, groups=groups)
        assert fig is not None
        plt.close(fig)

    def test_top_n(self, scores_df):
        """测试 top_n 参数"""
        fig = plot_pathway_dotplot(scores_df, top_n=5)
        assert fig is not None
        plt.close(fig)

    def test_save_to_file(self, scores_df, groups, tmp_path):
        """测试保存到文件"""
        output_file = str(tmp_path / "dotplot.png")
        fig = plot_pathway_dotplot(scores_df, groups=groups, output_file=output_file)
        assert os.path.exists(output_file)
        assert os.path.getsize(output_file) > 0
        plt.close(fig)


class TestPlotSampleCorrelation:
    """plot_sample_correlation 测试"""

    @pytest.fixture
    def scores_df(self):
        """创建模拟通路活性得分矩阵"""
        np.random.seed(42)
        pathways = [f"Pathway_{i}" for i in range(1, 11)]
        samples = [f"Normal_{i}" for i in range(1, 4)] + [f"Disease_{i}" for i in range(1, 4)]
        data = np.random.randn(len(pathways), len(samples))
        return pd.DataFrame(data, index=pathways, columns=samples)

    @pytest.fixture
    def annotation_col(self):
        """创建样本分组注释"""
        samples = [f"Normal_{i}" for i in range(1, 4)] + [f"Disease_{i}" for i in range(1, 4)]
        return pd.DataFrame({
            "Group": ["Normal"] * 3 + ["Disease"] * 3,
        }, index=samples)

    def test_pearson(self, scores_df):
        """测试 Pearson 相关性"""
        fig = plot_sample_correlation(scores_df, method="pearson")
        assert fig is not None
        assert isinstance(fig, matplotlib.figure.Figure)
        plt.close(fig)

    def test_spearman(self, scores_df):
        """测试 Spearman 相关性"""
        fig = plot_sample_correlation(scores_df, method="spearman")
        assert fig is not None
        plt.close(fig)

    def test_with_annotation(self, scores_df, annotation_col):
        """测试带分组注释"""
        fig = plot_sample_correlation(scores_df, annotation_col=annotation_col)
        assert fig is not None
        plt.close(fig)

    def test_invalid_method(self, scores_df):
        """测试无效 method 应抛出 ValueError"""
        with pytest.raises(ValueError, match="不支持的 method"):
            plot_sample_correlation(scores_df, method="kendall")

    def test_save_to_file(self, scores_df, tmp_path):
        """测试保存到文件"""
        output_file = str(tmp_path / "correlation.png")
        fig = plot_sample_correlation(scores_df, output_file=output_file)
        assert os.path.exists(output_file)
        assert os.path.getsize(output_file) > 0
        plt.close(fig)


class TestE2EWithRealData:
    """使用 E2E 测试数据的端到端测试"""

    @pytest.fixture
    def ssgsea_scores(self):
        """加载 ssGSEA 测试结果"""
        test_data_dir = Path(__file__).parent.parent / "test_data"
        csv_path = test_data_dir / "ssgsea_results.csv"
        if not csv_path.exists():
            pytest.skip(f"E2E 测试数据不存在: {csv_path}")
        df = pd.read_csv(csv_path, index_col=0)
        return df

    @pytest.fixture
    def gsva_scores(self):
        """加载 GSVA 测试结果"""
        test_data_dir = Path(__file__).parent.parent / "test_data"
        csv_path = test_data_dir / "gsva_results.csv"
        if not csv_path.exists():
            pytest.skip(f"E2E 测试数据不存在: {csv_path}")
        df = pd.read_csv(csv_path, index_col=0)
        return df

    @pytest.fixture
    def groups(self):
        """E2E 测试数据的分组"""
        return {
            "Normal": ["Normal_1", "Normal_2", "Normal_3"],
            "Disease": ["Disease_1", "Disease_2", "Disease_3"],
        }

    @pytest.fixture
    def annotation_col(self):
        """E2E 测试数据的分组注释"""
        samples = [f"Normal_{i}" for i in range(1, 4)] + [f"Disease_{i}" for i in range(1, 4)]
        return pd.DataFrame({
            "Group": ["Normal"] * 3 + ["Disease"] * 3,
        }, index=samples)

    def test_ssgsea_heatmap(self, ssgsea_scores, annotation_col, tmp_path):
        """ssGSEA 数据热图"""
        output_file = str(tmp_path / "ssgsea_heatmap.png")
        fig = plot_pathway_heatmap(
            ssgsea_scores, annotation_col=annotation_col, output_file=output_file
        )
        assert isinstance(fig, matplotlib.figure.Figure)
        assert os.path.exists(output_file)
        plt.close(fig)

    def test_gsva_heatmap(self, gsva_scores, annotation_col, tmp_path):
        """GSVA 数据热图"""
        output_file = str(tmp_path / "gsva_heatmap.png")
        fig = plot_pathway_heatmap(
            gsva_scores, annotation_col=annotation_col, output_file=output_file
        )
        assert isinstance(fig, matplotlib.figure.Figure)
        assert os.path.exists(output_file)
        plt.close(fig)

    def test_ssgsea_group_comparison(self, ssgsea_scores, groups, tmp_path):
        """ssGSEA 数据组间比较"""
        output_file = str(tmp_path / "ssgsea_group_comparison.png")
        fig = plot_group_comparison(
            ssgsea_scores, groups, plot_type="box", output_file=output_file
        )
        assert isinstance(fig, matplotlib.figure.Figure)
        assert os.path.exists(output_file)
        plt.close(fig)

    def test_ssgsea_dotplot(self, ssgsea_scores, groups, tmp_path):
        """ssGSEA 数据气泡图"""
        output_file = str(tmp_path / "ssgsea_dotplot.png")
        fig = plot_pathway_dotplot(
            ssgsea_scores, groups=groups, output_file=output_file
        )
        assert isinstance(fig, matplotlib.figure.Figure)
        assert os.path.exists(output_file)
        plt.close(fig)

    def test_ssgsea_correlation(self, ssgsea_scores, annotation_col, tmp_path):
        """ssGSEA 数据相关性热图"""
        output_file = str(tmp_path / "ssgsea_correlation.png")
        fig = plot_sample_correlation(
            ssgsea_scores, annotation_col=annotation_col, output_file=output_file
        )
        assert isinstance(fig, matplotlib.figure.Figure)
        assert os.path.exists(output_file)
        plt.close(fig)


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
