"""Tests for CLI integration with maintained R plotting workflows."""

from unittest.mock import patch
from pathlib import Path

import pandas as pd
import pytest

from allenricher.core.config import Config
from allenricher.cli import (
    _METHOD_PLOT_TYPES,
    _generate_plots,
    _normalize_ranked_genes,
    _parse_groups,
    _select_gsea_enrichment_terms,
    _write_running_es_file,
    create_parser,
)


def test_gsea_plot_types_exclude_retired_plots():
    """Verify that retired GSEA plots are not exposed by the CLI."""
    assert {"heatmap", "circos", "cnetplot", "nes_barplot", "dotplot"}.isdisjoint(
        _METHOD_PLOT_TYPES["gsea"]
    )
    assert "enrichment2" in _METHOD_PLOT_TYPES["gsea"]
    assert "lollipop" in _METHOD_PLOT_TYPES["gsea"]


def test_public_plot_inventory_matches_supported_workflows():
    """Keep the public method-to-plot mapping explicit and stable."""
    assert set(_METHOD_PLOT_TYPES["gsea"]) == {
        "enrichment", "enrichment2", "barplot", "lollipop", "ridgeplot", "emapplot",
    }
    assert set(_METHOD_PLOT_TYPES["ssgsea"]) == {
        "heatmap", "group_comparison", "correlation",
    }
    assert _METHOD_PLOT_TYPES["gsva"] == _METHOD_PLOT_TYPES["ssgsea"]


def test_groups_reject_sample_assigned_to_multiple_groups():
    with pytest.raises(ValueError, match="assigned to both"):
        _parse_groups("Control:S1,S2;Disease:S2,S3")


def test_retired_common_plots_generate_nothing(tmp_path):
    """Verify that retired common plots do not generate output."""
    generated = _generate_plots(
        method="gsea",
        results={"GO": pd.DataFrame({"Term_Name": ["A"]})},
        ranked_genes=None,
        gene_weights=None,
        gene_sets=None,
        expr_matrix=None,
        groups=None,
        plot_types=["volcano", "method_comparison", "network", "upset"],
        output_dir=str(tmp_path),
    )

    assert generated == []
    assert not (tmp_path / "common_plots").exists()


def test_r_plots_are_the_default_and_python_plots_is_explicit():
    parser = create_parser()
    assert parser.parse_args(["analyze", "-m", "gsea"]).use_r_plots is True
    assert parser.parse_args(["analyze", "-m", "gsea", "--python-plots"]).use_r_plots is False


def test_use_r_plots_option():
    """Verify that the CLI accepts the R plotting option."""
    parser = create_parser()
    args = parser.parse_args([
        "analyze", "-i", "genes.txt", "-m", "gsea", "--use-r-plots",
        "-pt", "enrichment,enrichment2,heatmap",
    ])
    assert args.use_r_plots is True
    assert args.plot_types == "enrichment,enrichment2,heatmap"


def test_use_r_plots_option_accepts_lollipop():
    """Lolipop should be parsed as GSEA R plot type."""
    parser = create_parser()
    args = parser.parse_args([
        "analyze", "-i", "genes.txt", "-m", "gsea", "--use-r-plots",
        "-pt", "lollipop,emapplot",
    ])
    assert args.plot_types == "lollipop,emapplot"


def test_role_specific_palette_options_are_parsed():
    args = create_parser().parse_args([
        "analyze", "-i", "genes.txt",
        "--categorical-palette", "okabe_ito",
        "--sequential-palette", "cividis",
        "--diverging-palette", "colorbrewer_prgn",
    ])
    assert args.categorical_palette == "okabe_ito"
    assert args.sequential_palette == "cividis"
    assert args.diverging_palette == "colorbrewer_prgn"


def test_style_choices_keep_aliases_but_remove_colorblind():
    parser = create_parser()
    for style in ("nature", "science", "presentation", "cell", "omicshare"):
        args = parser.parse_args(["analyze", "-i", "genes.txt", "--style", style])
        assert args.style == style
    with pytest.raises(SystemExit):
        parser.parse_args(["analyze", "-i", "genes.txt", "--style", "colorblind"])


def test_emapplot_filter_defaults_are_exposed():
    """The R emapplot filter parameters should have a clear default value."""
    parser = create_parser()
    args = parser.parse_args([
        "analyze", "-i", "genes.txt", "-m", "gsea", "--use-r-plots",
        "-pt", "emapplot",
    ])
    assert args.emapplot_qvalue == 0.05
    assert args.emapplot_min_count == 3
    assert args.emapplot_top_n == 30


def test_emapplot_filter_options_are_parsed():
    """The R emapplot filter parameters should be overwritten by CLI."""
    parser = create_parser()
    args = parser.parse_args([
        "analyze", "-i", "genes.txt", "-m", "gsea", "--use-r-plots",
        "-pt", "emapplot",
        "--emapplot-qvalue", "0.01",
        "--emapplot-min-count", "5",
        "--emapplot-top-n", "12",
    ])
    assert args.emapplot_qvalue == 0.01
    assert args.emapplot_min_count == 5
    assert args.emapplot_top_n == 12


def test_gsea_enrichment_top_up_down_defaults_are_exposed():
    parser = create_parser()
    args = parser.parse_args([
        "analyze", "-i", "genes.txt", "-m", "gsea", "-pt", "enrichment",
    ])
    assert args.gsea_enrichment_top_up == 5
    assert args.gsea_enrichment_top_down == 5


def test_gsea_multi_top_up_down_defaults_are_exposed():
    parser = create_parser()
    args = parser.parse_args([
        "analyze", "-i", "genes.txt", "-m", "gsea", "-pt", "enrichment2",
    ])
    assert args.gsea_multi_top_up == 3
    assert args.gsea_multi_top_down == 3


def test_gsea_multi_top_up_down_options_are_parsed():
    parser = create_parser()
    args = parser.parse_args([
        "analyze", "-i", "genes.txt", "-m", "gsea", "-pt", "enrichment2",
        "--gsea-multi-top-up", "2", "--gsea-multi-top-down", "4",
    ])
    assert args.gsea_multi_top_up == 2
    assert args.gsea_multi_top_down == 4


def test_select_gsea_enrichment_terms_uses_top_up_and_down():
    df = pd.DataFrame({
        "Term_ID": ["UP1", "UP2", "UP3", "DOWN1", "DOWN2", "DOWN3"],
        "NES": [2.5, 1.5, 3.0, -2.2, -1.2, -3.1],
        "FDR": [0.01, 0.02, 0.03, 0.01, 0.02, 0.03],
    })
    selected = _select_gsea_enrichment_terms(df, top_up=2, top_down=1)
    assert list(selected["Term_ID"]) == ["UP3", "UP1", "DOWN3"]


@patch("allenricher.visualization.gsea_plots.plot_gsea_multi_enrichment")
def test_python_multi_pathway_outputs_up_and_down_separately(mock_plot, tmp_path):
    results = pd.DataFrame({
        "Term_ID": ["UP1", "UP2", "UP3", "DOWN1", "DOWN2", "DOWN3"],
        "Term_Name": ["Up 1", "Up 2", "Up 3", "Down 1", "Down 2", "Down 3"],
        "NES": [3.0, 2.5, 2.0, -3.0, -2.5, -2.0],
        "FDR": [0.01] * 6,
    })
    gene_sets = {term_id: {"G1"} for term_id in results["Term_ID"]}

    files = _generate_plots(
        method="gsea",
        results={"TEST": results},
        ranked_genes=["G1", "G2"],
        gene_weights={"G1": 1.0, "G2": -1.0},
        gene_sets=gene_sets,
        expr_matrix=None,
        groups=None,
        plot_types=["enrichment2"],
        output_dir=str(tmp_path),
        plot_style="science",
        plot_palette="tol_vibrant",
    )

    assert mock_plot.call_count == 2
    assert all(call.kwargs["style"] == "science" for call in mock_plot.call_args_list)
    assert all(call.kwargs["palette"] == "tol_vibrant" for call in mock_plot.call_args_list)
    assert any(path.endswith("TEST_enrichment2_up.png") for path in files)
    assert any(path.endswith("TEST_enrichment2_down.png") for path in files)


@patch("allenricher.visualization.r_plotter.plot_gsea_enrichment2_r", return_value=True)
@patch("allenricher.visualization.r_plotter.check_r_environment", return_value=True)
def test_r_multi_pathway_outputs_up_and_down_separately(mock_r_env, mock_plot, tmp_path):
    results = pd.DataFrame({
        "Term_ID": ["UP1", "UP2", "UP3", "DOWN1", "DOWN2", "DOWN3"],
        "Term_Name": ["Up 1", "Up 2", "Up 3", "Down 1", "Down 2", "Down 3"],
        "NES": [3.0, 2.5, 2.0, -3.0, -2.5, -2.0],
        "FDR": [0.01] * 6,
    })
    gene_sets = {term_id: {"G1"} for term_id in results["Term_ID"]}

    files = _generate_plots(
        method="gsea",
        results={"TEST": results},
        ranked_genes=["G1", "G2"],
        gene_weights={"G1": 1.0, "G2": -1.0},
        gene_sets=gene_sets,
        expr_matrix=None,
        groups=None,
        plot_types=["enrichment2"],
        output_dir=str(tmp_path),
        use_r_plots=True,
    )

    assert mock_plot.call_count == 2
    assert any(path.endswith("TEST_enrichment2_up.png") for path in files)
    assert any(path.endswith("TEST_enrichment2_down.png") for path in files)


def test_normalize_ranked_genes_accepts_tuple_list():
    """Sorts the gene tuple list into a Gene list and weight dictionary."""
    genes, weights = _normalize_ranked_genes([("A", 2.0), ("B", -1.5)], None)
    assert genes == ["A", "B"]
    assert weights == {"A": 2.0, "B": -1.5}


def test_write_running_es_file_uses_real_ranked_data(tmp_path):
    """The running ES intermediate table should be derived from a true sequence of genes and genomes."""
    results = pd.DataFrame({
        "Term_ID": ["TERM1"],
        "Term_Name": ["Pathway 1"],
        "NES": [1.5],
    })
    out_file = tmp_path / "running_es.tsv"

    result = _write_running_es_file(
        out_file,
        results,
        [("A", 2.0), ("B", 1.0), ("C", -0.5)],
        None,
        {"TERM1": {"A", "C"}},
    )

    assert result == str(out_file)
    df = pd.read_csv(out_file, sep="\t")
    assert list(df["Gene"]) == ["A", "B", "C"]
    assert list(df["Hit"]) == [True, False, True]
    assert "Running_ES" in df.columns


def test_python_backend_generates_gsea_lollipop(tmp_path):
    """If not using R, Python back should also generate lolipop."""
    results = {
        "KEGG": pd.DataFrame({
            "Term_Name": ["Pathway A", "Pathway B", "Pathway C", "Pathway D"],
            "NES": [2.1, -1.8, 1.3, -1.1],
            "FDR": [0.001, 0.01, 0.02, 0.04],
            "setSize": [30, 26, 18, 12],
        })
    }

    generated = _generate_plots(
        method="gsea",
        results=results,
        ranked_genes=None,
        gene_weights=None,
        gene_sets=None,
        expr_matrix=None,
        groups=None,
        plot_types=["lollipop"],
        output_dir=str(tmp_path),
        use_r_plots=False,
    )

    expected = tmp_path / "gsea_plots" / "KEGG_lollipop.png"
    assert str(expected) in generated
    assert expected.exists()


def test_python_backend_generates_gsea_barplot_and_svg(tmp_path):
    results = {
        "KEGG": pd.DataFrame({
            "Term_Name": ["Pathway A", "Pathway B"],
            "NES": [2.1, -1.8],
            "FDR": [0.001, 0.01],
            "setSize": [30, 26],
        })
    }

    generated = _generate_plots(
        method="gsea", results=results, ranked_genes=None, gene_weights=None,
        gene_sets=None, expr_matrix=None, groups=None, plot_types=["barplot"],
        output_dir=str(tmp_path), plot_format="svg",
    )

    expected = tmp_path / "gsea_plots" / "KEGG_barplot.svg"
    assert generated == [str(expected)]
    assert expected.read_text(encoding="utf-8").lstrip().startswith("<?xml")


@patch("allenricher.visualization.gsea_plots.plot_gsea_barplot")
def test_python_gsea_plots_receive_database_term_names(mock_plot, tmp_path):
    results = {
        "KEGG": pd.DataFrame({
            "pathway": ["hsa00010", "hsa00020"],
            "pval": [0.001, 0.002],
            "padj": [0.01, 0.02],
            "log2err": [0.1, 0.1],
            "ES": [0.6, -0.5],
            "NES": [2.1, -1.8],
            "size": [30, 26],
            "leadingEdge": ["G1,G2", "G3,G4"],
        })
    }
    original_columns = results["KEGG"].columns.tolist()

    _generate_plots(
        method="gsea", results=results, ranked_genes=None, gene_weights=None,
        gene_sets=None, expr_matrix=None, groups=None, plot_types=["barplot"],
        output_dir=str(tmp_path),
        term_name_maps={
            "KEGG": {
                "hsa00010": "Metabolism|Carbohydrate metabolism|Glycolysis / Gluconeogenesis",
                "hsa00020": "Metabolism|Carbohydrate metabolism|Citrate cycle (TCA cycle)",
            }
        },
    )

    plotted = mock_plot.call_args.kwargs["results_df"]
    assert plotted["Term_Name"].str.contains("Glycolysis|Citrate cycle").all()
    assert results["KEGG"].columns.tolist() == original_columns


@patch("allenricher.visualization.r_plotter.plot_gsea_barplot_r", return_value=True)
@patch("allenricher.visualization.r_plotter.check_r_environment", return_value=True)
def test_r_gsea_sidecar_contains_names_but_official_frame_stays_unchanged(
    mock_r_env, mock_plot, tmp_path
):
    frame = pd.DataFrame({
        "pathway": ["hsa00010", "hsa00020"],
        "pval": [0.001, 0.002],
        "padj": [0.01, 0.02],
        "log2err": [0.1, 0.1],
        "ES": [0.6, -0.5],
        "NES": [2.1, -1.8],
        "size": [30, 26],
        "leadingEdge": ["G1,G2", "G3,G4"],
    })
    official_columns = frame.columns.tolist()

    _generate_plots(
        method="gsea", results={"KEGG": frame}, ranked_genes=None,
        gene_weights=None, gene_sets=None, expr_matrix=None, groups=None,
        plot_types=["barplot"], output_dir=str(tmp_path), use_r_plots=True,
        term_name_maps={
            "KEGG": {
                "hsa00010": "Metabolism|Carbohydrate metabolism|Glycolysis / Gluconeogenesis",
                "hsa00020": "Metabolism|Carbohydrate metabolism|Citrate cycle (TCA cycle)",
            }
        },
    )

    sidecar = pd.read_csv(tmp_path / "gsea_plots" / "KEGG_enrichment.tsv", sep="\t")
    assert sidecar["Term_Name"].str.contains("Glycolysis|Citrate cycle").all()
    assert sidecar["Database"].eq("KEGG").all()
    assert frame.columns.tolist() == official_columns
    assert mock_plot.call_count == 1


def test_python_backend_generates_each_requested_format(tmp_path):
    results = {
        "KEGG": pd.DataFrame({
            "Term_Name": ["Pathway A", "Pathway B"],
            "NES": [2.1, -1.8],
            "FDR": [0.001, 0.01],
            "setSize": [30, 26],
        })
    }
    generated = _generate_plots(
        method="gsea", results=results, ranked_genes=None, gene_weights=None,
        gene_sets=None, expr_matrix=None, groups=None, plot_types=["lollipop"],
        output_dir=str(tmp_path), plot_formats=["png", "pdf", "svg"],
    )
    assert {Path(path).suffix for path in generated} == {".png", ".pdf", ".svg"}
    assert all(Path(path).stat().st_size > 0 for path in generated)


def test_python_backend_generates_gsea_ridgeplot(tmp_path):
    """Python back should use ranked lights and gene sets to generate ridgeplot."""
    ranked_genes = [f"G{i}" for i in range(1, 21)]
    gene_weights = {gene: 3.0 - index * 0.3 for index, gene in enumerate(ranked_genes)}
    gene_sets = {
        "T1": set(ranked_genes[:8]),
        "T2": set(ranked_genes[6:16]),
    }
    results = {
        "KEGG": pd.DataFrame({
            "Term_ID": ["T1", "T2"],
            "Term_Name": ["Pathway A", "Pathway B"],
            "NES": [2.1, -1.8],
            "p_value": [0.001, 0.01],
        })
    }

    generated = _generate_plots(
        method="gsea",
        results=results,
        ranked_genes=ranked_genes,
        gene_weights=gene_weights,
        gene_sets=gene_sets,
        expr_matrix=None,
        groups=None,
        plot_types=["ridgeplot"],
        output_dir=str(tmp_path),
        use_r_plots=False,
    )

    expected = tmp_path / "gsea_plots" / "KEGG_ridgeplot.png"
    assert str(expected) in generated
    assert expected.exists() and expected.stat().st_size > 0


def test_r_emapplot_receives_filter_options(tmp_path):
    """_gendere_prots should pass the Eapplot special filter parameter to Rwrapper."""
    results = {
        "KEGG": pd.DataFrame({
            "Term_ID": ["T1", "T2"],
            "Term_Name": ["Pathway A", "Pathway B"],
            "NES": [2.1, -1.8],
            "FDR": [0.001, 0.02],
            "Gene_Count": [8, 6],
            "Genes": ["A/B/C/D", "B/C/E/F"],
        })
    }

    with patch("allenricher.visualization.r_plotter.check_r_environment", return_value=True), \
         patch("allenricher.visualization.r_plotter.plot_gsea_emapplot_r", return_value=True) as plot:
        generated = _generate_plots(
            method="gsea",
            results=results,
            ranked_genes=None,
            gene_weights=None,
            gene_sets=None,
            expr_matrix=None,
            groups=None,
            plot_types=["emapplot"],
            output_dir=str(tmp_path),
            use_r_plots=True,
            emapplot_qvalue=0.01,
            emapplot_min_count=5,
            emapplot_top_n=12,
            plot_style="science",
            plot_palette="tol_sunset",
        )

    expected = tmp_path / "gsea_plots" / "KEGG_emapplot.png"
    assert str(expected) in generated
    assert plot.call_args.kwargs == {
        "top_n": 12, "qvalue": 0.01, "min_count": 5, "dpi": 300,
        "style": "science", "palette": "tol_sunset",
    }


def test_ssgsea_r_heatmap_uses_activity_wrapper(tmp_path):
    scores = pd.DataFrame(
        {"S1": [0.8, -0.4], "S2": [0.6, -0.2], "S3": [-0.5, 0.9], "S4": [-0.3, 0.7]},
        index=["Pathway A", "Pathway B"],
    )
    groups = {"Control": ["S1", "S2"], "Disease": ["S3", "S4"]}

    with patch("allenricher.visualization.r_plotter.check_r_environment", return_value=True), \
         patch("allenricher.visualization.r_plotter.plot_activity_heatmap_r", return_value=True) as plot:
        generated = _generate_plots(
            method="ssgsea",
            results={"GO": scores},
            ranked_genes=None,
            gene_weights=None,
            gene_sets=None,
            expr_matrix=None,
            groups=groups,
            plot_types=["heatmap"],
            output_dir=str(tmp_path),
            use_r_plots=True,
        )

    assert plot.call_count == 1
    assert plot.call_args.kwargs == {
        "analysis_method": "ssgsea", "top_n": 40, "dpi": 300,
        "style": "nature", "palette": None,
    }
    assert generated == [str(tmp_path / "gsea_plots" / "activity_heatmap.png")]


@pytest.mark.parametrize(
    ("method", "expected_title"),
    [("ssgsea", "ssGSEA Pathway Activity"), ("gsva", "GSVA Pathway Activity")],
)
def test_python_activity_heatmap_title_matches_method(tmp_path, method, expected_title):
    scores = pd.DataFrame({"S1": [0.8, -0.4], "S2": [-0.5, 0.9]}, index=["A", "B"])
    with patch("allenricher.visualization.gsva_plots.plot_pathway_heatmap") as plot:
        _generate_plots(
            method=method,
            results={"GO": scores},
            ranked_genes=None,
            gene_weights=None,
            gene_sets=None,
            expr_matrix=None,
            groups=None,
            plot_types=["heatmap"],
            output_dir=str(tmp_path),
        )

    assert plot.call_args.kwargs["title"] == expected_title


def test_activity_heatmap_limits_display_rows_but_not_source_matrix(tmp_path):
    scores = pd.DataFrame(
        {
            "S1": range(40),
            "S2": [value * 2 for value in range(40)],
            "S3": [value * -1 for value in range(40)],
        },
        index=[f"P{value}" for value in range(40)],
    )
    with patch("allenricher.visualization.gsva_plots.plot_pathway_heatmap") as plot:
        _generate_plots(
            method="gsva",
            results={"GO": scores},
            ranked_genes=None,
            gene_weights=None,
            gene_sets=None,
            expr_matrix=None,
            groups=None,
            plot_types=["heatmap"],
            output_dir=str(tmp_path),
            activity_heatmap_top_n=12,
        )

    plotted_scores = plot.call_args.args[0]
    assert plotted_scores.shape == (12, 3)
    assert scores.shape == (40, 3)


def test_activity_heatmap_has_dedicated_default_row_limit():
    config = Config()
    assert config.activity_heatmap_top_n == 40
    assert config.top_terms == 20


def test_gsva_r_correlation_uses_correlation_wrapper(tmp_path):
    scores = pd.DataFrame(
        {"S1": [0.8, -0.4], "S2": [0.6, -0.2], "S3": [-0.5, 0.9]},
        index=["Pathway A", "Pathway B"],
    )
    groups = {"Control": ["S1", "S2"], "Disease": ["S3"]}

    with patch("allenricher.visualization.r_plotter.check_r_environment", return_value=True), \
         patch("allenricher.visualization.r_plotter.plot_sample_correlation_r", return_value=True) as plot:
        generated = _generate_plots(
            method="gsva",
            results={"GO": scores},
            ranked_genes=None,
            gene_weights=None,
            gene_sets=None,
            expr_matrix=None,
            groups=groups,
            plot_types=["correlation"],
            output_dir=str(tmp_path),
            use_r_plots=True,
        )

    assert plot.call_count == 1
    assert generated == [str(tmp_path / "gsea_plots" / "sample_correlation.png")]


def test_ssgsea_r_group_comparison_writes_plot_and_statistics(tmp_path):
    scores = pd.DataFrame(
        {
            "S1": [0.8, -0.4], "S2": [0.6, -0.2],
            "S3": [-0.5, 0.9], "S4": [-0.3, 0.7],
        },
        index=["Pathway A", "Pathway B"],
    )
    groups = {"Control": ["S1", "S2"], "Disease": ["S3", "S4"]}

    with patch("allenricher.visualization.r_plotter.check_r_environment", return_value=True), \
         patch("allenricher.visualization.r_plotter.plot_group_comparison_r", return_value=True) as plot:
        generated = _generate_plots(
            method="ssgsea",
            results={"GO": scores},
            ranked_genes=None,
            gene_weights=None,
            gene_sets=None,
            expr_matrix=None,
            groups=groups,
            plot_types=["group_comparison"],
            output_dir=str(tmp_path),
            use_r_plots=True,
        )

    plot_dir = tmp_path / "gsea_plots"
    plot.assert_called_once_with(
        str(plot_dir / "activity_scores.tsv"),
        str(plot_dir / "sample_metadata.tsv"),
        str(plot_dir / "group_comparison.png"),
        str(plot_dir / "group_comparison.statistics.tsv"),
        dpi=300,
        style="nature",
        palette=None,
    )
    assert generated == [str(plot_dir / "group_comparison.png")]


def test_activity_plots_are_generated_for_each_database(tmp_path):
    scores = pd.DataFrame(
        {
            "S1": [0.8, -0.4], "S2": [0.6, -0.2],
            "S3": [-0.5, 0.9], "S4": [-0.3, 0.7],
        },
        index=["Pathway A", "Pathway B"],
    )
    groups = {"Control": ["S1", "S2"], "Disease": ["S3", "S4"]}

    with patch("allenricher.visualization.r_plotter.check_r_environment", return_value=True), \
         patch("allenricher.visualization.r_plotter.plot_group_comparison_r", return_value=True) as plot:
        generated = _generate_plots(
            method="ssgsea",
            results={"GO": scores, "PFAM": scores},
            ranked_genes=None,
            gene_weights=None,
            gene_sets=None,
            expr_matrix=None,
            groups=groups,
            plot_types=["group_comparison"],
            output_dir=str(tmp_path),
            use_r_plots=True,
        )

    plot_dir = tmp_path / "gsea_plots"
    assert plot.call_count == 2
    assert generated == [
        str(plot_dir / "GO_group_comparison.png"),
        str(plot_dir / "PFAM_group_comparison.png"),
    ]
    assert [call.args[3] for call in plot.call_args_list] == [
        str(plot_dir / "GO_group_comparison.statistics.tsv"),
        str(plot_dir / "PFAM_group_comparison.statistics.tsv"),
    ]


@patch("allenricher.visualization.gsea_plots.plot_gsea_enrichment")
def test_gsea_enrichment_uses_database_specific_gene_sets(mock_plot, tmp_path):
    frame = pd.DataFrame({
        "Term_ID": ["SHARED"],
        "Term_Name": ["Shared term"],
        "NES": [2.0],
        "FDR": [0.01],
        "ES": [0.6],
        "pval": [0.001],
        "size": [1],
    })
    gene_sets_by_database = {
        "GO": {"SHARED": {"G1"}},
        "CustomDB": {"SHARED": {"G2"}},
    }

    generated = _generate_plots(
        method="gsea",
        results={"GO": frame.copy(), "CustomDB": frame.copy()},
        ranked_genes=["G1", "G2", "G3"],
        gene_weights={"G1": 2.0, "G2": 1.5, "G3": -1.0},
        gene_sets={"SHARED": {"WRONG"}},
        gene_sets_by_database=gene_sets_by_database,
        expr_matrix=None,
        groups=None,
        plot_types=["enrichment"],
        output_dir=str(tmp_path),
    )

    used_sets = [call.kwargs["gene_set"] for call in mock_plot.call_args_list]
    assert {"G1"} in used_sets
    assert {"G2"} in used_sets
    assert {"WRONG"} not in used_sets
    assert any(path.endswith("GO_SHARED_enrichment.png") for path in generated)
    assert any(path.endswith("CustomDB_SHARED_enrichment.png") for path in generated)


@patch("allenricher.visualization.r_plotter.plot_gsea_enrichment_r", return_value=True)
@patch("allenricher.visualization.r_plotter.check_r_environment", return_value=True)
def test_r_gsea_running_es_is_database_scoped(mock_r_env, mock_plot, tmp_path):
    frame = pd.DataFrame({
        "Term_ID": ["SHARED"],
        "Term_Name": ["Shared term"],
        "NES": [2.0],
        "FDR": [0.01],
        "ES": [0.6],
        "pval": [0.001],
        "size": [1],
    })
    gene_sets_by_database = {
        "GO": {"SHARED": {"G1"}},
        "CustomDB": {"SHARED": {"G2"}},
    }

    _generate_plots(
        method="gsea",
        results={"GO": frame.copy(), "CustomDB": frame.copy()},
        ranked_genes=["G1", "G2", "G3"],
        gene_weights={"G1": 2.0, "G2": 1.5, "G3": -1.0},
        gene_sets={"SHARED": {"WRONG"}},
        gene_sets_by_database=gene_sets_by_database,
        expr_matrix=None,
        groups=None,
        plot_types=["enrichment"],
        output_dir=str(tmp_path),
        use_r_plots=True,
    )

    go_running = pd.read_csv(tmp_path / "gsea_plots" / "GO_running_es.tsv", sep="\t")
    custom_running = pd.read_csv(tmp_path / "gsea_plots" / "CustomDB_running_es.tsv", sep="\t")
    assert go_running.loc[go_running["Hit"], "Gene"].tolist() == ["G1"]
    assert custom_running.loc[custom_running["Hit"], "Gene"].tolist() == ["G2"]
    assert mock_plot.call_count == 2
