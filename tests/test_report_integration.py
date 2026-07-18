"""
Report Generator Integration Test - GSEA/GSVA Visualize embedding
===========================================

Test ReportGenerator The... GSEA/GSVA Visualize embedded functionality: 
- Add parameters without compromising existing functionality (Recompatible)
- GSEA The result embedded in the report
- GSVA The result embedded in the report
- None GSEA/GSVA Report normal production at the end of the process
- Report HTML Organisation base64 Picture Tab
- plot_types Filter
"""

import sys
import os
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch, MagicMock

import pytest
import pandas as pd

# Ensure that root directory in sys.path
sys.path.insert(0, str(Path(__file__).parent.parent))

from allenricher.report.generator import ReportGenerator


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

@pytest.fixture
def output_dir(tmp_path):
    """Create temporary output directory"""
    return str(tmp_path / "report_output")


@pytest.fixture
def sample_results():
    """Create example rich result"""
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
    """Create Example GSEA Analysis"""
    return pd.DataFrame({
        "pathway": ["HSA_Cell_Cycle", "HSA_Apoptosis", "HSA_PI3K_AKT", "HSA_MAPK"],
        "nes": [2.5, -1.8, 1.5, -2.1],
        "pvalue": [0.001, 0.005, 0.01, 0.002],
        "gene_count": [45, 30, 25, 35],
        "es": [0.65, -0.45, 0.38, -0.52],
    })


@pytest.fixture
def gsva_results_df():
    """Create Example ssGSEA/GSVAActive matrix"""
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
    """Create Example Group"""
    return {
        "Normal": ["Normal_1", "Normal_2", "Normal_3"],
        "Disease": ["Disease_1", "Disease_2", "Disease_3"],
    }


@pytest.fixture
def ranked_genes():
    """Create a list of genes for example sorting"""
    return ["GENE001", "GENE002", "GENE003", "GENE004", "GENE005"]


@pytest.fixture
def gene_weights():
    """Create an example of gene weight"""
    return {"GENE001": 5.0, "GENE002": 4.0, "GENE003": 3.0, "GENE004": 2.0, "GENE005": 1.0}


@pytest.fixture
def gene_sets():
    """Create an example gene set"""
    return {
        "HSA_Cell_Cycle": {"GENE001", "GENE002", "GENE003"},
        "HSA_Apoptosis": {"GENE002", "GENE004"},
    }


# ===========================================================================
# 1. Retrocompatibility tests
# ===========================================================================

class TestBackwardCompatibility:
    """Test for additional parameters without destroying existing functionality"""

    def test_generate_without_new_params(self, output_dir, sample_results, tmp_path):
        """The generate method should work normally without sending additional parameters"""
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
        assert 'href="https://github.com/zhangducsu/AllEnricher-v2"' in html
        assert "BMC Bioinformatics. 2020;21:106." in html
        # Should not contain GSEA/GSVARegional
        assert "gsea-plots" not in html
        assert "gsva-plots" not in html
        assert ".header-content" in html
        assert ".header .meta" in html
        assert "display: block" in html

    def test_generate_no_results_backward_compat(self, output_dir, tmp_path):
        """Compatibility with no result"""
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
        assert 'href="https://github.com/zhangducsu/AllEnricher-v2"' in html
        assert "BMC Bioinformatics. 2020;21:106." in html
        assert "Materials and Methods Writing Reference" in html

    @pytest.mark.parametrize(
        "method, method_text",
        [
            ("hypergeometric", "over-representation analysis (ORA)"),
            ("gsea", "gene set enrichment analysis (GSEA)"),
            ("ssgsea", "single-sample gene set enrichment analysis (ssGSEA)"),
            ("gsva", "gene set variation analysis (GSVA)"),
        ],
    )
    def test_methods_reference_is_embedded_for_each_analysis_method(
        self, output_dir, tmp_path, method, method_text
    ):
        results = {
            "GO": pd.DataFrame(
                {
                    "Term_ID": ["GO:0001"],
                    "Term_Name": ["Example pathway"],
                    "P_Value": [0.01],
                    "Adjusted_P_Value": [0.02],
                    "Sample_1": [0.4],
                }
            )
        }
        metadata = {
            "allenricher_version": "2.1-test",
            "analysis_method": method,
            "species": "hsa",
            "databases": ["GO"],
            "parameters": {"ssgsea_tau": 0.25},
        }
        output_file = tmp_path / f"{method}_methods_report.html"

        ReportGenerator(output_dir).generate(
            results,
            str(output_file),
            analysis_method=method,
            metadata=metadata,
            pvalue_cutoff=1.0,
            qvalue_cutoff=1.0,
        )

        html = output_file.read_text(encoding="utf-8")
        assert "Materials and Methods Writing Reference" in html
        assert "AllEnricher version 2.1-test" in html
        assert method_text in html

    def test_tf_report_template_links_repository_and_shows_citation(self):
        template = (
            Path(__file__).parent.parent
            / "allenricher"
            / "report"
            / "templates"
            / "tf_report.html"
        ).read_text(encoding="utf-8")

        assert 'href="https://github.com/zhangducsu/AllEnricher-v2"' in template
        assert "BMC Bioinformatics. 2020;21:106." in template

    def test_generate_all_filtered_out_backward_compat(self, output_dir, tmp_path):
        """All results are filtered back in compatibility."""
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
# 2. GSEA results embedded testing
# ===========================================================================

class TestGSEAIntegration:
    """Test GSEA results embedded in report"""

    def test_gsea_results_embedded_in_report(
        self, output_dir, sample_results, gsea_results_df, tmp_path
    ):
        """GSEA results should be correctly embedded in HTML reports"""
        gen = ReportGenerator(output_dir)
        output_file = str(tmp_path / "report_gsea.html")

        result = gen.generate(
            results=sample_results,
            output_file=output_file,
            gene_list=["TP53"],
            gsea_results=gsea_results_df,
            plot_types=["enrichment"],
        )

        assert os.path.exists(result)
        with open(result, 'r', encoding='utf-8') as f:
            html = f.read()

        # Verify GSEA Area Exists
        assert "gsea-plots" in html
        assert "GSEA Visualization" in html
        # Verify that the navigation bar contains GSEA links
        assert 'href="#gsea-plots"' in html

    @patch("allenricher.visualization.gsea_plots.plot_gsea_enrichment")
    def test_gsea_enrichment_plots(
        self, mock_enrichment,
        output_dir, sample_results, gsea_results_df,
        ranked_genes, gene_weights, gene_sets, tmp_path
    ):
        """Test a curved image generation"""
        mock_fig = MagicMock()
        mock_enrichment.return_value = mock_fig

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

        # Verify that the enrichment chart is called (for the access that exists in the gene_sets only)
        assert mock_enrichment.called

# ===========================================================================
# 3. GSVA Results embedded for testing
# ===========================================================================

class TestGSVAIntegration:
    """Test GSVA Results embedded in report"""

    @patch("allenricher.visualization.gsva_plots.plot_sample_correlation")
    @patch("allenricher.visualization.gsva_plots.plot_pathway_heatmap")
    def test_gsva_results_embedded_in_report(
        self, mock_heatmap, mock_corr, output_dir, sample_results,
        gsva_results_df, tmp_path
    ):
        """GSVA results should be correctly embedded in HTML reports"""
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
            analysis_method="gsva",
            plot_types=["heatmap", "correlation"],  # Exclude group_comparison
        )

        assert os.path.exists(result)
        with open(result, 'r', encoding='utf-8') as f:
            html = f.read()

        # Verify GSVA Area Exists
        assert "gsva-plots" in html
        assert "GSVA Visualization" in html
        assert mock_heatmap.call_args.kwargs["title"] == "GSVA Pathway Activity"
        # Validation Navigation Bar contains GSVA links
        assert 'href="#gsva-plots"' in html

    @patch("allenricher.visualization.gsva_plots.plot_group_comparison")
    @patch("allenricher.visualization.gsva_plots.plot_sample_correlation")
    @patch("allenricher.visualization.gsva_plots.plot_pathway_heatmap")
    def test_gsva_with_groups(
        self, mock_heatmap, mock_corr, mock_group,
        output_dir, sample_results, gsva_results_df, gsva_groups, tmp_path
    ):
        """Test GSVA visualization with grouping"""
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

        # Group_comparison should be called (group information available)
        assert mock_group.called

    @patch("allenricher.visualization.gsva_plots.plot_sample_correlation")
    @patch("allenricher.visualization.gsva_plots.plot_pathway_heatmap")
    def test_gsva_base64_images(
        self, mock_heatmap, mock_corr, output_dir, sample_results,
        gsva_results_df, tmp_path
    ):
        """GSVA report should contain picture of base64 code"""
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt

        def create_and_save_fig(*args, **kwargs):
            fig, ax = plt.subplots(figsize=(6, 4))
            ax.imshow([[1, 2], [3, 4]])
            # Simulate saving behaviour: Save if output_file
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
# 4. N/A GSEA/GSVAResults test
# ===========================================================================

class TestNoGseaGsvaResults:
    """Report normal generation when no GSEA/GSVA results are tested"""

    def test_no_gsea_gsva_params(self, output_dir, sample_results, tmp_path):
        """No GSEA/GSVAReport when parameters should be generated properly"""
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
        """Report should be generated normally when gsea_reults=Noone"""
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
        """The report should be generated normally when gsea_reults are empty"""
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
# 5. Plot_types filtering test
# ===========================================================================

class TestPlotTypesFilter:
    """Tests for filters"""

    def test_retired_gsea_plot_type_is_ignored(
        self, output_dir, sample_results, gsea_results_df, tmp_path
    ):
        """The NES barplot that is retired should no longer be called by the report generator."""
        gen = ReportGenerator(output_dir)
        output_file = str(tmp_path / "report_filter.html")

        # The old configuration is still readable but does not generate decommissioned charts.
        result = gen.generate(
            results=sample_results,
            output_file=output_file,
            gsea_results=gsea_results_df,
            plot_types=["nes_barplot"],
        )

        assert os.path.exists(result)

    @patch("allenricher.visualization.gsva_plots.plot_sample_correlation")
    @patch("allenricher.visualization.gsva_plots.plot_pathway_heatmap")
    def test_gsva_plot_types_filter(
        self, mock_heatmap, mock_corr, output_dir, sample_results,
        gsva_results_df, tmp_path
    ):
        """Generate only the specified GSVA chart type"""
        mock_fig = MagicMock()
        mock_heatmap.return_value = mock_fig
        mock_corr.return_value = mock_fig

        gen = ReportGenerator(output_dir)
        output_file = str(tmp_path / "report_gsva_filter.html")

        # Request only claretation, not headmap
        result = gen.generate(
            results=sample_results,
            output_file=output_file,
            gsva_results=gsva_results_df,
            plot_types=["correlation"],
        )

        assert os.path.exists(result)
        mock_corr.assert_called_once()
        mock_heatmap.assert_not_called()

    def test_default_plot_types(
        self, output_dir, sample_results, gsea_results_df, tmp_path
    ):
        """Default values should be used when not specifying plot_types"""
        gen = ReportGenerator(output_dir)
        output_file = str(tmp_path / "report_default_types.html")

        result = gen.generate(
            results=sample_results,
            output_file=output_file,
            gsea_results=gsea_results_df,
            # Do not send plot_types, use default enrichment.
        )

        assert os.path.exists(result)


# ===========================================================================
# 6. GSEA + GSVA concurrent presence of tests
# ===========================================================================

class TestCombinedGseaGsva:
    """Test GSEA and GSVA results exist simultaneously"""

    @patch("allenricher.visualization.gsva_plots.plot_sample_correlation")
    @patch("allenricher.visualization.gsva_plots.plot_pathway_heatmap")
    def test_both_gsea_and_gsva(
        self, mock_heatmap, mock_corr,
        output_dir, sample_results, gsea_results_df, gsva_results_df, tmp_path
    ):
        """GSEA and GSVA results should be embedded in the report when they exist"""
        mock_fig = MagicMock()
        mock_heatmap.return_value = mock_fig
        mock_corr.return_value = mock_fig

        gen = ReportGenerator(output_dir)
        output_file = str(tmp_path / "report_combined.html")

        result = gen.generate(
            results=sample_results,
            output_file=output_file,
            gsea_results=gsea_results_df,
            gsva_results=gsva_results_df,
            plot_types=["enrichment", "heatmap", "correlation"],
        )

        assert os.path.exists(result)
        with open(result, 'r', encoding='utf-8') as f:
            html = f.read()

        # Verify that both areas exist.
        assert "gsea-plots" in html
        assert "gsva-plots" in html
        # Verify that the navigation bar contains two links
        assert 'href="#gsea-plots"' in html
        assert 'href="#gsva-plots"' in html
        # Validate GSEA/GSVA between plots and tables
        plots_pos = html.find('id="plots"')
        gsea_pos = html.find('id="gsea-plots"')
        gsva_pos = html.find('id="gsva-plots"')
        tables_pos = html.find('id="GO-table"')
        assert plots_pos < gsea_pos < tables_pos or plots_pos < gsva_pos < tables_pos


class TestReportCompletenessContract:
    """The report must faithfully present the official results sheet and all the maps that have been generated by CLI."""

    def test_gsea_report_embeds_all_saved_plots_and_readable_term_metadata(
        self, output_dir, tmp_path
    ):
        plot_dir = Path(output_dir) / "gsea_plots"
        plot_dir.mkdir(parents=True)
        (plot_dir / "PATH_001_enrichment.png").write_bytes(b"png")
        (plot_dir / "GO_lollipop.png").write_bytes(b"png")
        (plot_dir / "GO_lollipop.svg").write_text("<svg/>", encoding="utf-8")
        (plot_dir / "GO_ridgeplot.png").write_bytes(b"png")

        results = {
            "GO": pd.DataFrame({
                "Term_ID": ["PATH:001"],
                "Term_Name": ["Readable pathway name"],
                "Hierarchy": ["Class A|Subclass B|Readable pathway name"],
                "pathway": ["PATH:001"],
                "pval": [0.2],
                "padj": [0.3],
                "log2err": [0.01],
                "ES": [0.5],
                "NES": [1.6],
                "size": [20],
                "leadingEdge": ["G1,G2,G3"],
            })
        }
        config = SimpleNamespace(pvalue_cutoff=1.0, qvalue_cutoff=1.0)
        output_file = tmp_path / "gsea_report.html"

        ReportGenerator(output_dir, config).generate(
            results,
            str(output_file),
            analysis_method="gsea",
        )
        html = output_file.read_text(encoding="utf-8")

        assert "GSEA Enrichment Analysis Report" in html
        assert "Readable pathway name (PATH:001) - GSEA Enrichment Plot" in html
        assert "Class A|Subclass B|Readable pathway name" in html
        assert "3 plot types" in html
        assert html.count("<img ") == 3
        assert ">PNG</a>" in html and ">SVG</a>" in html
        assert 'href="#GO-table">GO Results</a>' in html
        assert "fonts.googleapis.com" not in html
        assert 'class="no-results-box"' not in html
        assert "fullValue ? fullValue.dataset.full" in html

    def test_activity_report_uses_explicit_term_fields_safe_ids_and_method_label(
        self, output_dir, tmp_path
    ):
        results = {
            "TF <Custom>": pd.DataFrame({
                "Term_ID": ["TF:001"],
                "Term_Name": ["Regulatory <script>alert(1)</script> pathway"],
                "Hierarchy": ["TF class|Regulatory pathway"],
                "Sample_1": [0.25],
                "Sample_2": [-0.15],
            })
        }
        output_file = tmp_path / "activity_report.html"

        ReportGenerator(output_dir).generate(
            results,
            str(output_file),
            analysis_method="ssgsea",
            ai_interpretation={"TF <Custom>": "**Result** <script>bad()</script>"},
        )
        html = output_file.read_text(encoding="utf-8")

        assert "ssGSEA Pathway Activity Report" in html
        assert "GSVA Pathway Activity Report" not in html
        assert 'id="TF-Custom-table"' in html
        assert '<td>TF:001</td><td>Regulatory &lt;script&gt;alert(1)&lt;/script&gt; pathway</td>' in html
        assert "TF class|Regulatory pathway" in html
        assert "AI Interpretation" in html
        assert "<script>bad()" not in html
        assert "&lt;script&gt;bad()&lt;/script&gt;" in html


def test_report_keeps_analysis_successful_when_ai_interpretation_fails(
    output_dir, sample_results, tmp_path
):
    config = SimpleNamespace(pvalue_cutoff=1.0, qvalue_cutoff=1.0)
    output_file = tmp_path / "ai_failure_report.html"
    ReportGenerator(output_dir, config).generate(
        sample_results,
        str(output_file),
        analysis_method="hypergeometric",
        ai_interpretation_error={
            "error_code": "AI_INTERPRETATION_FAILED",
            "backend": "deepseek",
            "mode": "summary",
            "message": "Error: Missing credentials",
        },
    )
    html = output_file.read_text(encoding="utf-8")
    assert "Analysis Summary" in html
    assert "AI interpretation unavailable." in html
    assert "AI_INTERPRETATION_FAILED" in html
    assert "Error: Missing credentials" in html


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
