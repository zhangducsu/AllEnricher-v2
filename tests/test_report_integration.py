"""
报告生成器集成测试 - GSEA/GSVA 可视化嵌入
===========================================

测试 ReportGenerator 的 GSEA/GSVA 可视化嵌入功能：
- 新增参数不破坏现有功能（向后兼容）
- GSEA 结果嵌入报告
- GSVA 结果嵌入报告
- 无 GSEA/GSVA 结果时报告正常生成
- 报告 HTML 包含 base64 图片标签
- plot_types 过滤功能
"""

import sys
import os
from pathlib import Path
from unittest.mock import patch, MagicMock

import pytest
import pandas as pd

# 确保项目根目录在 sys.path 中
sys.path.insert(0, str(Path(__file__).parent.parent))

from allenricher.report.generator import ReportGenerator


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

@pytest.fixture
def output_dir(tmp_path):
    """创建临时输出目录"""
    return str(tmp_path / "report_output")


@pytest.fixture
def sample_results():
    """创建示例富集结果"""
    results = {
        "GO": pd.DataFrame({
            "Term_ID": ["GO:0005576", "GO:0005737"],
            "Term_Name": ["extracellular region", "cytoplasm"],
            "Gene_Count": [10, 8],
            "Rich_Factor": [0.05, 0.04],
            "P_Value": [1e-5, 1e-3],
            "Adjusted_P_Value": [1e-3, 0.05],
            "Genes": ["TP53,BRCA1,EGFR", "AKT1,MAPK1,STAT3"],
            "Term_URL": ["http://example.com/GO:0005576", ""],
        }),
        "KEGG": pd.DataFrame({
            "Term_ID": ["hsa04110", "hsa04010"],
            "Term_Name": ["Cell Cycle", "MAPK signaling pathway"],
            "Gene_Count": [8, 12],
            "Rich_Factor": [0.03, 0.05],
            "P_Value": [1e-4, 1e-3],
            "Adjusted_P_Value": [1e-2, 0.05],
            "Genes": ["CDK1,CDC20", "MAPK1,EGFR"],
            "Term_URL": ["", ""],
        }),
    }
    return results


@pytest.fixture
def gsea_results_df():
    """创建示例 GSEA 分析结果"""
    return pd.DataFrame({
        "pathway": ["HSA_Cell_Cycle", "HSA_Apoptosis", "HSA_PI3K_AKT", "HSA_MAPK"],
        "nes": [2.5, -1.8, 1.5, -2.1],
        "pvalue": [0.001, 0.005, 0.01, 0.002],
        "gene_count": [45, 30, 25, 35],
        "es": [0.65, -0.45, 0.38, -0.52],
    })


@pytest.fixture
def gsva_results_df():
    """创建示例 ssGSEA/GSVA 活性矩阵"""
    return pd.DataFrame({
        "Normal_1": [0.1, 0.2, 0.8, 0.3, 0.5],
        "Normal_2": [0.15, 0.25, 0.7, 0.35, 0.45],
        "Normal_3": [0.12, 0.18, 0.75, 0.28, 0.52],
        "Disease_1": [0.9, 0.8, 0.2, 0.7, 0.9],
        "Disease_2": [0.85, 0.75, 0.25, 0.65, 0.88],
        "Disease_3": [0.88, 0.82, 0.18, 0.72, 0.85],
    }, index=["HSA_Cell_Cycle", "HSA_Apoptosis", "HSA_PI3K_AKT", "HSA_MAPK", "HSA_DNA_Repair"])


@pytest.fixture
def gsva_groups():
    """创建示例分组"""
    return {
        "Normal": ["Normal_1", "Normal_2", "Normal_3"],
        "Disease": ["Disease_1", "Disease_2", "Disease_3"],
    }


@pytest.fixture
def ranked_genes():
    """创建示例排序基因列表"""
    return ["GENE001", "GENE002", "GENE003", "GENE004", "GENE005"]


@pytest.fixture
def gene_weights():
    """创建示例基因权重"""
    return {"GENE001": 5.0, "GENE002": 4.0, "GENE003": 3.0, "GENE004": 2.0, "GENE005": 1.0}


@pytest.fixture
def gene_sets():
    """创建示例基因集"""
    return {
        "HSA_Cell_Cycle": {"GENE001", "GENE002", "GENE003"},
        "HSA_Apoptosis": {"GENE002", "GENE004"},
    }


def _write_minimal_png(path):
    """写入一个最小 PNG 文件，用于报告嵌入测试。"""
    path.write_bytes(
        b"\x89PNG\r\n\x1a\n"
        b"\x00\x00\x00\rIHDR\x00\x00\x00\x01\x00\x00\x00\x01"
        b"\x08\x06\x00\x00\x00\x1f\x15\xc4\x89"
        b"\x00\x00\x00\nIDATx\x9cc`\x00\x00\x00\x02\x00\x01"
        b"\xe2!\xbc3\x00\x00\x00\x00IEND\xaeB`\x82"
    )


# ===========================================================================
# 1. 向后兼容测试
# ===========================================================================

class TestBackwardCompatibility:
    """测试新增参数不破坏现有功能"""

    def test_generate_without_new_params(self, output_dir, sample_results, tmp_path):
        """不传新增参数时，generate 方法应正常工作"""
        gen = ReportGenerator(output_dir)
        output_file = str(tmp_path / "report.html")
        gene_list = ["TP53", "BRCA1", "EGFR"]

        result = gen.generate(
            results=sample_results,
            output_file=output_file,
            gene_list=gene_list,
        )

        assert os.path.exists(result)
        with open(result, 'r', encoding='utf-8') as f:
            html = f.read()
        assert "AllEnricher" in html
        assert "GO" in html
        assert "KEGG" in html
        # 不应包含 GSEA/GSVA 区域
        assert "gsea-plots" not in html
        assert "gsva-plots" not in html

    def test_generate_no_results_backward_compat(self, output_dir, tmp_path):
        """无结果时向后兼容"""
        gen = ReportGenerator(output_dir)
        output_file = str(tmp_path / "report.html")

        result = gen.generate(
            results={},
            output_file=output_file,
        )

        assert os.path.exists(result)
        with open(result, 'r', encoding='utf-8') as f:
            html = f.read()
        assert "AllEnricher" in html

    def test_generate_all_filtered_out_backward_compat(self, output_dir, tmp_path):
        """所有结果被过滤掉时向后兼容"""
        gen = ReportGenerator(output_dir)
        output_file = str(tmp_path / "report.html")

        results = {
            "GO": pd.DataFrame({
                "Term_ID": ["GO:0005576"],
                "Term_Name": ["test"],
                "Gene_Count": [1],
                "Rich_Factor": [0.01],
                "P_Value": [0.5],
                "Adjusted_P_Value": [0.9],
                "Genes": ["GENE1"],
            }),
        }

        result = gen.generate(
            results=results,
            output_file=output_file,
            pvalue_cutoff=0.05,
            qvalue_cutoff=0.05,
        )

        assert os.path.exists(result)


# ===========================================================================
# 2. GSEA 结果嵌入测试
# ===========================================================================

class TestGSEAIntegration:
    """测试 GSEA 结果嵌入报告"""

    @patch("allenricher.visualization.gsea_plots.plot_gsea_nes_barplot")
    @patch("allenricher.visualization.gsea_plots.plot_gsea_dotplot")
    def test_gsea_results_embedded_in_report(
        self, mock_dotplot, mock_barplot, output_dir, sample_results,
        gsea_results_df, tmp_path
    ):
        """GSEA 结果应正确嵌入到 HTML 报告中"""
        # 模拟图表生成返回 Figure 对象
        mock_fig = MagicMock()
        mock_barplot.return_value = mock_fig
        mock_dotplot.return_value = mock_fig

        gen = ReportGenerator(output_dir)
        output_file = str(tmp_path / "report_gsea.html")

        result = gen.generate(
            results=sample_results,
            output_file=output_file,
            gene_list=["TP53"],
            gsea_results=gsea_results_df,
            plot_types=["nes_barplot", "dotplot"],  # 排除 enrichment（需要更多参数）
        )

        assert os.path.exists(result)
        with open(result, 'r', encoding='utf-8') as f:
            html = f.read()

        # 验证 GSEA 区域存在
        assert "gsea-plots" in html
        assert "GSEA Visualization" in html
        # 验证导航栏包含 GSEA 链接
        assert 'href="#gsea-plots"' in html

    @patch("allenricher.visualization.gsea_plots.plot_gsea_enrichment")
    @patch("allenricher.visualization.gsea_plots.plot_gsea_nes_barplot")
    @patch("allenricher.visualization.gsea_plots.plot_gsea_dotplot")
    def test_gsea_enrichment_plots(
        self, mock_dotplot, mock_barplot, mock_enrichment,
        output_dir, sample_results, gsea_results_df,
        ranked_genes, gene_weights, gene_sets, tmp_path
    ):
        """测试富集曲线图生成"""
        mock_fig = MagicMock()
        mock_enrichment.return_value = mock_fig
        mock_barplot.return_value = mock_fig
        mock_dotplot.return_value = mock_fig

        gen = ReportGenerator(output_dir)
        output_file = str(tmp_path / "report_gsea_full.html")

        result = gen.generate(
            results=sample_results,
            output_file=output_file,
            gsea_results=gsea_results_df,
            gsea_ranked_genes=ranked_genes,
            gsea_gene_weights=gene_weights,
            gsea_gene_sets=gene_sets,
            analysis_method="gsea",
        )

        assert os.path.exists(result)
        with open(result, 'r', encoding='utf-8') as f:
            html = f.read()

        # 验证 enrichment 图表被调用（仅对 gene_sets 中存在的通路）
        assert mock_enrichment.called

    @patch("allenricher.visualization.gsea_plots.plot_gsea_nes_barplot")
    @patch("allenricher.visualization.gsea_plots.plot_gsea_dotplot")
    def test_gsea_base64_images(
        self, mock_dotplot, mock_barplot, output_dir, sample_results,
        gsea_results_df, tmp_path
    ):
        """报告应包含 base64 编码的图片"""
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt

        def create_and_save_fig(*args, **kwargs):
            fig, ax = plt.subplots(figsize=(6, 4))
            ax.barh([0, 1], [1, 2])
            # 模拟 _save_figure 行为：如果传了 output_file 就保存
            output_file = kwargs.get("output_file")
            if output_file:
                fig.savefig(output_file, dpi=150, bbox_inches="tight")
            return fig

        mock_barplot.side_effect = create_and_save_fig
        mock_dotplot.side_effect = create_and_save_fig

        gen = ReportGenerator(output_dir)
        output_file = str(tmp_path / "report_gsea_b64.html")

        result = gen.generate(
            results=sample_results,
            output_file=output_file,
            gsea_results=gsea_results_df,
            plot_types=["nes_barplot", "dotplot"],
        )

        assert os.path.exists(result)
        with open(result, 'r', encoding='utf-8') as f:
            html = f.read()

        # 验证包含 base64 图片标签
        assert "data:image/png;base64," in html

        plt.close("all")


# ===========================================================================
# 3. GSVA 结果嵌入测试
# ===========================================================================

class TestGSVAIntegration:
    """测试 GSVA 结果嵌入报告"""

    @patch("allenricher.visualization.gsva_plots.plot_sample_correlation")
    @patch("allenricher.visualization.gsva_plots.plot_pathway_heatmap")
    def test_gsva_results_embedded_in_report(
        self, mock_heatmap, mock_corr, output_dir, sample_results,
        gsva_results_df, tmp_path
    ):
        """GSVA 结果应正确嵌入到 HTML 报告中"""
        mock_fig = MagicMock()
        mock_heatmap.return_value = mock_fig
        mock_corr.return_value = mock_fig

        gen = ReportGenerator(output_dir)
        output_file = str(tmp_path / "report_gsva.html")

        result = gen.generate(
            results=sample_results,
            output_file=output_file,
            gene_list=["TP53"],
            gsva_results=gsva_results_df,
            plot_types=["heatmap", "correlation"],  # 排除 group_comparison
        )

        assert os.path.exists(result)
        with open(result, 'r', encoding='utf-8') as f:
            html = f.read()

        # 验证 GSVA 区域存在
        assert "gsva-plots" in html
        assert "ssGSEA/GSVA Visualization" in html
        # 验证导航栏包含 GSVA 链接
        assert 'href="#gsva-plots"' in html

    @patch("allenricher.visualization.gsva_plots.plot_group_comparison")
    @patch("allenricher.visualization.gsva_plots.plot_sample_correlation")
    @patch("allenricher.visualization.gsva_plots.plot_pathway_heatmap")
    def test_gsva_with_groups(
        self, mock_heatmap, mock_corr, mock_group,
        output_dir, sample_results, gsva_results_df, gsva_groups, tmp_path
    ):
        """测试带分组的 GSVA 可视化"""
        mock_fig = MagicMock()
        mock_heatmap.return_value = mock_fig
        mock_corr.return_value = mock_fig
        mock_group.return_value = mock_fig

        gen = ReportGenerator(output_dir)
        output_file = str(tmp_path / "report_gsva_groups.html")

        result = gen.generate(
            results=sample_results,
            output_file=output_file,
            gsva_results=gsva_results_df,
            gsva_groups=gsva_groups,
            analysis_method="gsva",
        )

        assert os.path.exists(result)
        with open(result, 'r', encoding='utf-8') as f:
            html = f.read()

        # group_comparison 应被调用（有分组信息）
        assert mock_group.called

    @patch("allenricher.visualization.gsva_plots.plot_sample_correlation")
    @patch("allenricher.visualization.gsva_plots.plot_pathway_heatmap")
    def test_gsva_base64_images(
        self, mock_heatmap, mock_corr, output_dir, sample_results,
        gsva_results_df, tmp_path
    ):
        """GSVA 报告应包含 base64 编码的图片"""
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt

        def create_and_save_fig(*args, **kwargs):
            fig, ax = plt.subplots(figsize=(6, 4))
            ax.imshow([[1, 2], [3, 4]])
            # 模拟保存行为：如果传了 output_file 就保存
            output_file = kwargs.get("output_file")
            if output_file:
                fig.savefig(output_file, dpi=150, bbox_inches="tight", facecolor="white")
            return fig

        mock_heatmap.side_effect = create_and_save_fig
        mock_corr.side_effect = create_and_save_fig

        gen = ReportGenerator(output_dir)
        output_file = str(tmp_path / "report_gsva_b64.html")

        result = gen.generate(
            results=sample_results,
            output_file=output_file,
            gsva_results=gsva_results_df,
            plot_types=["heatmap", "correlation"],
        )

        assert os.path.exists(result)
        with open(result, 'r', encoding='utf-8') as f:
            html = f.read()

        assert "data:image/png;base64," in html

        plt.close("all")


# ===========================================================================
# 4. 无 GSEA/GSVA 结果测试
# ===========================================================================

class TestNoGseaGsvaResults:
    """测试无 GSEA/GSVA 结果时报告正常生成"""

    def test_no_gsea_gsva_params(self, output_dir, sample_results, tmp_path):
        """不传 GSEA/GSVA 参数时报告应正常生成"""
        gen = ReportGenerator(output_dir)
        output_file = str(tmp_path / "report_no_extra.html")

        result = gen.generate(
            results=sample_results,
            output_file=output_file,
            gene_list=["TP53"],
        )

        assert os.path.exists(result)
        with open(result, 'r', encoding='utf-8') as f:
            html = f.read()

        assert "gsea-plots" not in html
        assert "gsva-plots" not in html
        assert "GO" in html

    def test_none_gsea_results(self, output_dir, sample_results, tmp_path):
        """gsea_results=None 时报告应正常生成"""
        gen = ReportGenerator(output_dir)
        output_file = str(tmp_path / "report_none_gsea.html")

        result = gen.generate(
            results=sample_results,
            output_file=output_file,
            gsea_results=None,
            gsva_results=None,
        )

        assert os.path.exists(result)
        with open(result, 'r', encoding='utf-8') as f:
            html = f.read()

        assert "gsea-plots" not in html
        assert "gsva-plots" not in html

    def test_empty_gsea_results(self, output_dir, sample_results, tmp_path):
        """gsea_results 为空 DataFrame 时报告应正常生成"""
        gen = ReportGenerator(output_dir)
        output_file = str(tmp_path / "report_empty_gsea.html")

        empty_df = pd.DataFrame(columns=["pathway", "nes", "pvalue", "gene_count"])

        result = gen.generate(
            results=sample_results,
            output_file=output_file,
            gsea_results=empty_df,
        )

        assert os.path.exists(result)
        with open(result, 'r', encoding='utf-8') as f:
            html = f.read()

        assert "gsea-plots" not in html


# ===========================================================================
# 5. 已生成 GSEA R 图嵌入测试
# ===========================================================================

class TestGeneratedGseaPlotEmbedding:
    """测试已生成的 GSEA R 图会进入主报告 Visualization 区域"""

    def test_generated_r_gsea_plots_are_embedded(self, output_dir):
        gen = ReportGenerator(output_dir)
        gsea_plot_dir = Path(output_dir) / "gsea_plots"
        gsea_plot_dir.mkdir(parents=True, exist_ok=True)

        plot_names = [
            "KEGG_nes_barplot.png",
            "KEGG_dotplot.png",
            "KEGG_barplot.png",
            "KEGG_ridgeplot.png",
            "KEGG_emapplot.png",
            "KEGG_cnetplot.png",
            "KEGG_circos.png",
            "KEGG_enrichment2.png",
            "KEGG_heatmap.png",
            "hsa04110_enrichment.png",
        ]
        for plot_name in plot_names:
            _write_minimal_png(gsea_plot_dir / plot_name)

        results = {
            "KEGG": pd.DataFrame({
                "Term_ID": ["hsa04110"],
                "Term_Name": ["Cell cycle"],
                "NES": [2.1],
                "ES": [0.7],
                "P_Value": [0.001],
                "Adjusted_P_Value": [0.01],
                "Set_Size": [32],
                "Lead_genes": ["CDK1;CCNB1"],
                "matched_genes": ["CDK1;CCNB1;CDC20"],
            })
        }

        html = gen._generate_plot_section(results)

        assert html.count("data:image/png;base64,") == len(plot_names)
        for caption in [
            "GSEA signed NES ranking",
            "GSEA pathway summary",
            "GSEA top pathway bar plot",
            "Running enrichment score distribution",
            "Pathway overlap map",
            "Pathway-gene network",
            "Pathway-gene circos overview",
            "Multi-pathway enrichment trajectories",
            "Expression heatmap for selected genes",
            "Single-pathway enrichment trajectory",
        ]:
            assert caption in html


# ===========================================================================
# 6. plot_types 过滤测试
# ===========================================================================

class TestPlotTypesFilter:
    """测试 plot_types 过滤功能"""

    @patch("allenricher.visualization.gsea_plots.plot_gsea_nes_barplot")
    @patch("allenricher.visualization.gsea_plots.plot_gsea_dotplot")
    def test_gsea_plot_types_filter(
        self, mock_dotplot, mock_barplot, output_dir, sample_results,
        gsea_results_df, tmp_path
    ):
        """只生成指定的图表类型"""
        mock_fig = MagicMock()
        mock_barplot.return_value = mock_fig
        mock_dotplot.return_value = mock_fig

        gen = ReportGenerator(output_dir)
        output_file = str(tmp_path / "report_filter.html")

        # 只请求 nes_barplot，不请求 dotplot
        result = gen.generate(
            results=sample_results,
            output_file=output_file,
            gsea_results=gsea_results_df,
            plot_types=["nes_barplot"],
        )

        assert os.path.exists(result)
        mock_barplot.assert_called_once()
        mock_dotplot.assert_not_called()

    @patch("allenricher.visualization.gsva_plots.plot_sample_correlation")
    @patch("allenricher.visualization.gsva_plots.plot_pathway_heatmap")
    def test_gsva_plot_types_filter(
        self, mock_heatmap, mock_corr, output_dir, sample_results,
        gsva_results_df, tmp_path
    ):
        """只生成指定的 GSVA 图表类型"""
        mock_fig = MagicMock()
        mock_heatmap.return_value = mock_fig
        mock_corr.return_value = mock_fig

        gen = ReportGenerator(output_dir)
        output_file = str(tmp_path / "report_gsva_filter.html")

        # 只请求 correlation，不请求 heatmap
        result = gen.generate(
            results=sample_results,
            output_file=output_file,
            gsva_results=gsva_results_df,
            plot_types=["correlation"],
        )

        assert os.path.exists(result)
        mock_corr.assert_called_once()
        mock_heatmap.assert_not_called()

    @patch("allenricher.visualization.gsea_plots.plot_gsea_nes_barplot")
    @patch("allenricher.visualization.gsea_plots.plot_gsea_dotplot")
    def test_default_plot_types(
        self, mock_dotplot, mock_barplot, output_dir, sample_results,
        gsea_results_df, tmp_path
    ):
        """不指定 plot_types 时应使用默认值"""
        mock_fig = MagicMock()
        mock_barplot.return_value = mock_fig
        mock_dotplot.return_value = mock_fig

        gen = ReportGenerator(output_dir)
        output_file = str(tmp_path / "report_default_types.html")

        result = gen.generate(
            results=sample_results,
            output_file=output_file,
            gsea_results=gsea_results_df,
            # 不传 plot_types，应使用默认值 ["enrichment", "nes_barplot", "dotplot"]
        )

        assert os.path.exists(result)
        # enrichment 需要额外参数所以不会调用，但 nes_barplot 和 dotplot 应该被调用
        mock_barplot.assert_called_once()
        mock_dotplot.assert_called_once()


# ===========================================================================
# 7. GSEA + GSVA 同时存在测试
# ===========================================================================

class TestCombinedGseaGsva:
    """测试 GSEA 和 GSVA 结果同时存在"""

    @patch("allenricher.visualization.gsva_plots.plot_sample_correlation")
    @patch("allenricher.visualization.gsva_plots.plot_pathway_heatmap")
    @patch("allenricher.visualization.gsea_plots.plot_gsea_nes_barplot")
    @patch("allenricher.visualization.gsea_plots.plot_gsea_dotplot")
    def test_both_gsea_and_gsva(
        self, mock_dotplot, mock_barplot, mock_heatmap, mock_corr,
        output_dir, sample_results, gsea_results_df, gsva_results_df, tmp_path
    ):
        """GSEA 和 GSVA 结果同时存在时应都嵌入报告"""
        mock_fig = MagicMock()
        mock_barplot.return_value = mock_fig
        mock_dotplot.return_value = mock_fig
        mock_heatmap.return_value = mock_fig
        mock_corr.return_value = mock_fig

        gen = ReportGenerator(output_dir)
        output_file = str(tmp_path / "report_combined.html")

        result = gen.generate(
            results=sample_results,
            output_file=output_file,
            gsea_results=gsea_results_df,
            gsva_results=gsva_results_df,
            plot_types=["nes_barplot", "dotplot", "heatmap", "correlation"],
        )

        assert os.path.exists(result)
        with open(result, 'r', encoding='utf-8') as f:
            html = f.read()

        # 验证两个区域都存在
        assert "gsea-plots" in html
        assert "gsva-plots" in html
        # 验证导航栏包含两个链接
        assert 'href="#gsea-plots"' in html
        assert 'href="#gsva-plots"' in html
        # 验证 GSEA/GSVA 区域在 plots 和 tables 之间
        plots_pos = html.find('id="plots"')
        gsea_pos = html.find('id="gsea-plots"')
        gsva_pos = html.find('id="gsva-plots"')
        tables_pos = html.find('id="GO-table"')
        assert plots_pos < gsea_pos < tables_pos or plots_pos < gsva_pos < tables_pos


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
