import shutil

import pandas as pd
import pandas.testing as pdt
import pytest

from allenricher.core.bioconductor import FGSEA_COLUMNS, run_fgsea, windows_to_wsl_path
from allenricher.core.config import Config
from allenricher.core.enrichment import (
    EnrichmentAnalyzer,
    SSGSEA,
    add_result_term_metadata,
)
from allenricher.core.gsva import GSVA
from allenricher.report.generator import ReportGenerator


def _fgsea_frame() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "pathway": ["UP", "DOWN"],
            "pval": [0.001, 0.02],
            "padj": [0.002, 0.03],
            "log2err": [0.2, 0.3],
            "ES": [0.8, -0.7],
            "NES": [2.1, -1.9],
            "size": [12, 10],
            "leadingEdge": ["G1;G2", "G99;G100"],
        }
    )


def test_windows_paths_translate_for_wsl_without_touching_values():
    assert windows_to_wsl_path(r"D:\Project Folder\result.tsv") == "/mnt/d/Project Folder/result.tsv"
    assert windows_to_wsl_path("Gaussian") == "Gaussian"


def test_gsea_export_adds_term_metadata_and_preserves_fgsea_values(monkeypatch, tmp_path):
    expected = _fgsea_frame()
    monkeypatch.setattr(
        "allenricher.core.bioconductor.run_fgsea",
        lambda *args, **kwargs: expected.copy(),
    )
    analyzer = EnrichmentAnalyzer(
        Config(method="gsea", gsea_min_size=1, gsea_max_size=100, n_jobs=1)
    )
    results = analyzer.run_analysis(
        {"G1"},
        {"G1", "G2", "G99", "G100"},
        {
            "TEST": {
                "UP": {
                    "name": "Metabolism|Carbohydrate|Up pathway",
                    "hierarchy": "Metabolism|Carbohydrate|Up pathway",
                    "genes": {"G1", "G2"},
                },
                "DOWN": {"name": "Down", "genes": {"G99", "G100"}},
            }
        },
        parallel=False,
        ranked_gene_list=[("G1", 2.0), ("G2", 1.0), ("G99", -1.0), ("G100", -2.0)],
    )

    exported = results["TEST"]
    assert exported["Term_ID"].tolist() == ["UP", "DOWN"]
    assert exported["Term_Name"].tolist() == ["Up pathway", "Down"]
    assert exported["Hierarchy"].tolist() == [
        "Metabolism|Carbohydrate|Up pathway", "",
    ]
    pdt.assert_frame_equal(exported[FGSEA_COLUMNS], expected)
    analyzer.save_results(str(tmp_path), metadata={"allenricher_version": "test"})
    assert (tmp_path / "TEST_enrichment.tsv").read_text(encoding="utf-8").startswith("Term_ID\tTerm_Name\tHierarchy\t")
    saved = pd.read_csv(tmp_path / "TEST_enrichment.tsv", sep="\t")
    pdt.assert_frame_equal(saved[FGSEA_COLUMNS], expected)


def test_ora_result_is_pure_reproducible_tsv(tmp_path):
    analyzer = EnrichmentAnalyzer(Config(method="hypergeometric", n_jobs=1))
    analyzer.run_analysis(
        {"G1", "G2"},
        {"G1", "G2", "G3", "G4"},
        {"TEST": {"T1": {
            "name": "Parent|Child|Term one",
            "hierarchy": "Parent|Child|Term one",
            "genes": {"G1", "G2", "G3"},
        }}},
        parallel=False,
    )
    analyzer.save_results(
        str(tmp_path),
        metadata={"allenricher_version": "test", "analysis_date": "changing-value"},
    )

    result_path = tmp_path / "TEST_enrichment.tsv"
    text = result_path.read_text(encoding="utf-8")
    assert text.startswith("Term_ID\t")
    assert "changing-value" not in text
    assert not text.startswith("#")
    saved = pd.read_csv(result_path, sep="\t")
    assert saved.loc[0, "Term_Name"] == "Term one"
    assert saved.loc[0, "Hierarchy"] == "Parent|Child|Term one"


def test_activity_export_keeps_term_ids_names_hierarchy_and_scores():
    official = pd.DataFrame({"S1": [0.4], "S2": [-0.2]}, index=["P1"])
    official.index.name = "Pathway"

    exported = add_result_term_metadata(
        official,
        {"P1": {"name": "Class A|Subclass B|Pathway one"}},
    )

    assert exported.index.name == "Term_ID"
    assert exported.index.tolist() == ["P1"]
    assert exported.loc["P1", "Term_Name"] == "Pathway one"
    assert exported.loc["P1", "Hierarchy"] == "Class A|Subclass B|Pathway one"
    assert exported.loc["P1", "S1"] == pytest.approx(0.4)


def test_gsea_significant_filter_uses_official_p_columns(monkeypatch):
    expected = _fgsea_frame()
    monkeypatch.setattr(
        "allenricher.core.bioconductor.run_fgsea",
        lambda *args, **kwargs: expected.copy(),
    )
    analyzer = EnrichmentAnalyzer(
        Config(
            method="gsea",
            gsea_min_size=1,
            gsea_max_size=100,
            output_all=False,
            pvalue_cutoff=0.01,
            qvalue_cutoff=0.01,
            n_jobs=1,
        )
    )
    results = analyzer.run_analysis(
        {"G1"},
        {"G1", "G2"},
        {"TEST": {"UP": {"name": "Up", "genes": {"G1", "G2"}}}},
        parallel=False,
        ranked_gene_list=[("G1", 1.0), ("G2", -1.0)],
    )
    assert results["TEST"]["pathway"].tolist() == ["UP"]


def test_ssgsea_and_gsva_keep_official_matrix(monkeypatch):
    expected = pd.DataFrame({"S1": [0.4], "S2": [-0.2]}, index=["Pathway_A"])
    expected.index.name = "Pathway"
    calls = []

    def fake_run_gsva(*args, **kwargs):
        calls.append(kwargs)
        return expected.copy()

    monkeypatch.setattr("allenricher.core.bioconductor.run_gsva", fake_run_gsva)
    expression = pd.DataFrame({"S1": [1.0], "S2": [2.0]}, index=["G1"])
    gene_sets = {"Pathway_A": {"G1"}}

    pdt.assert_frame_equal(SSGSEA(min_size=1).analyze_matrix(expression, gene_sets), expected)
    pdt.assert_frame_equal(GSVA(min_size=1).analyze_matrix(expression, gene_sets), expected)
    assert calls[0]["method"] == "ssgsea"
    assert calls[1]["method"] == "gsva"


@pytest.mark.parametrize("method", ["ssgsea", "gsva"])
def test_matrix_methods_reject_generic_gene_list_entrypoint(method):
    analyzer = EnrichmentAnalyzer(Config(method=method, n_jobs=1))

    with pytest.raises(ValueError, match="analyze_activity_database"):
        analyzer.run_analysis({"G1"}, {"G1"}, {}, parallel=False)


def test_tf_gsea_and_ssgsea_reuse_official_backends(monkeypatch):
    from allenricher.analysis.tf_enrichment import TFEnrichmentAnalyzer

    database = {
        "gene2tf": pd.DataFrame(
            {"Gene": ["G1", "G2", "G3"], "TF_A": [1, 1, 1]}
        )
    }
    analyzer = TFEnrichmentAnalyzer(database)
    fgsea_expected = _fgsea_frame().iloc[[0]].assign(pathway="TF_A")
    activity_expected = pd.DataFrame({"S1": [0.5]}, index=["TF_A"])
    activity_expected.index.name = "Pathway"

    monkeypatch.setattr(
        "allenricher.core.bioconductor.run_fgsea",
        lambda *args, **kwargs: fgsea_expected.copy(),
    )
    monkeypatch.setattr(
        "allenricher.core.bioconductor.run_gsva",
        lambda *args, **kwargs: activity_expected.copy(),
    )

    pdt.assert_frame_equal(
        analyzer.gsea([("G1", 2.0), ("G2", 1.0), ("G3", -1.0)]),
        fgsea_expected,
    )
    expression = pd.DataFrame({"S1": [3.0, 2.0, 1.0]}, index=["G1", "G2", "G3"])
    tf_activity = analyzer.ssgsea(expression, min_size=1)
    assert tf_activity.index.name == "TF"
    assert tf_activity.loc["TF_A", "S1"] == pytest.approx(0.5)


def test_tf_gsea_size_filter_is_applied_after_rank_intersection(monkeypatch):
    from allenricher.analysis.tf_enrichment import TFEnrichmentAnalyzer

    database = {
        "gene2tf": pd.DataFrame({
            "Gene": ["G1", "G2", "G3", "G4"],
            "TF_A": [1, 1, 1, 1],
        })
    }
    captured = {}

    def fake_fgsea(ranked_genes, gene_sets, **kwargs):
        captured["ranked_genes"] = ranked_genes
        captured["gene_sets"] = gene_sets
        captured["options"] = kwargs
        return pd.DataFrame(columns=FGSEA_COLUMNS)

    monkeypatch.setattr("allenricher.core.bioconductor.run_fgsea", fake_fgsea)
    TFEnrichmentAnalyzer(database).gsea(
        [("G1", 2.0), ("G2", -1.0)], min_size=1, max_size=2
    )

    assert captured["gene_sets"]["TF_A"] == {"G1", "G2"}
    assert captured["options"]["min_size"] == 1
    assert captured["options"]["max_size"] == 2


@pytest.mark.parametrize(
    ("method", "expected"),
    [
        ("ora", (3, None)),
        ("gsea", (15, 5000)),
        ("ssgsea", (1, None)),
        ("gsva", (1, None)),
    ],
)
def test_tf_method_specific_size_defaults(method, expected):
    from allenricher.analysis.tf_enrichment import resolve_tf_size_limits

    assert resolve_tf_size_limits(method) == expected


def test_generic_tf_gsea_uses_post_intersection_tf_defaults(monkeypatch):
    captured = {}

    def fake_fgsea(ranked_genes, gene_sets, **kwargs):
        captured["gene_sets"] = gene_sets
        captured["options"] = kwargs
        return pd.DataFrame(columns=FGSEA_COLUMNS)

    monkeypatch.setattr("allenricher.core.bioconductor.run_fgsea", fake_fgsea)
    ranked = [(f"G{i}", float(5100 - i)) for i in range(5100)]
    term_data = {
        "SMALL": {"name": "small", "genes": [f"G{i}" for i in range(14)]},
        "VALID": {"name": "valid", "genes": [f"G{i}" for i in range(15)]},
        "LARGE": {"name": "large", "genes": [f"G{i}" for i in range(5001)]},
    }
    analyzer = EnrichmentAnalyzer(Config(method="gsea", n_jobs=1))

    analyzer.analyze_database(set(), set(), term_data, "TRRUST", ranked)

    assert set(captured["gene_sets"]) == {"VALID"}
    assert captured["options"]["min_size"] == 1
    assert captured["options"]["max_size"] == 5100


def test_generic_tf_ssgsea_keeps_single_gene_sets(monkeypatch):
    captured = {}

    def fake_gsva(expression, gene_sets, **kwargs):
        captured["gene_sets"] = gene_sets
        captured["options"] = kwargs
        return pd.DataFrame({"S1": [1.0]}, index=["ONE"])

    monkeypatch.setattr("allenricher.core.bioconductor.run_gsva", fake_gsva)
    expression = pd.DataFrame({"S1": [1.0, 2.0]}, index=["G1", "G2"])
    analyzer = EnrichmentAnalyzer(Config(method="ssgsea", n_jobs=1))

    analyzer.analyze_activity_database(
        expression,
        {"ONE": {"G1"}, "NONE": {"OUTSIDE"}},
        "TRRUST",
    )

    assert captured["gene_sets"] == {"ONE": {"G1"}}
    assert captured["options"]["min_size"] == 1
    assert captured["options"]["max_size"] == 2


@pytest.mark.parametrize(
    ("method", "expected"),
    [
        ("gsea", (15, 500)),
        ("ssgsea", (1, None)),
        ("gsva", (1, None)),
    ],
)
def test_generic_method_specific_size_defaults(method, expected):
    configured = EnrichmentAnalyzer(Config(method=method)).method
    assert (configured.min_size, configured.max_size) == expected


def test_explicit_size_limits_override_method_defaults():
    configured = EnrichmentAnalyzer(
        Config(method="ssgsea", gsea_min_size=4, gsea_max_size=120)
    ).method
    assert (configured.min_size, configured.max_size) == (4, 120)


def test_tf_ora_uses_the_annotated_gene_universe():
    from scipy.stats import hypergeom
    from allenricher.analysis.tf_enrichment import TFEnrichmentAnalyzer

    database = {
        "gene2tf": pd.DataFrame(
            {
                "Gene": ["G1", "G2", "G3", "G4"],
                "TF_A": [1, 0, 0, 0],
                "TF_B": [0, 1, 1, 1],
            }
        )
    }
    result = TFEnrichmentAnalyzer(database).ora(
        ["G1", "OUTSIDE"], min_overlap=1, min_size=1
    )

    assert result.loc[0, "GeneSet_Size"] == 1
    assert result.loc[0, "Pvalue"] == pytest.approx(hypergeom.sf(0, 4, 1, 1))


@pytest.mark.parametrize(
    ("ranked", "message"),
    [
        ([(" ", 2.0), ("G2", -1.0)], "empty gene ID"),
        ([("G1", 2.0), (" G1 ", -1.0)], "duplicate gene IDs"),
        ([("G1", 1.0), ("G2", 1.0)], "two distinct values"),
    ],
)
def test_run_fgsea_direct_api_rejects_invalid_ranked_vectors(ranked, message):
    with pytest.raises(ValueError, match=message):
        run_fgsea(ranked, {"TF1": {"G1"}}, min_size=1)


def test_report_consumes_official_tables(tmp_path):
    generator = ReportGenerator(str(tmp_path))
    gsea_html = generator._generate_tables({"TEST": _fgsea_frame()}, analysis_method="gsea")
    assert all(f">{column}<" in gsea_html for column in FGSEA_COLUMNS)

    activity = pd.DataFrame({"S1": [0.4], "S2": [-0.2]}, index=["Pathway_A"])
    activity.index.name = "Pathway"
    activity_html = generator._generate_tables({"TEST": activity}, analysis_method="ssgsea")
    assert "Pathway_A" in activity_html
    assert ">S1<" in activity_html
    assert "Gene Count" not in activity_html


def test_python_and_running_es_consumers_accept_fgsea_table(tmp_path):
    import matplotlib.pyplot as plt

    from allenricher.cli import _write_running_es_file
    from allenricher.visualization.gsea_plots import plot_gsea_lollipop
    from allenricher.visualization.plotter import Plotter

    frame = _fgsea_frame()
    output = tmp_path / "lollipop.png"
    figure = plot_gsea_lollipop(frame, output_file=str(output), dpi=100)
    plt.close(figure)
    assert output.stat().st_size > 0

    barplot = tmp_path / "barplot.png"
    Plotter(str(tmp_path)).plot_barplot(frame, "TEST", barplot.name, top_n=2)
    assert barplot.stat().st_size > 0

    running_es = _write_running_es_file(
        tmp_path / "running_es.tsv",
        frame,
        [("G1", 2.0), ("G2", 1.0), ("G99", -1.0), ("G100", -2.0)],
        {"G1": 2.0, "G2": 1.0, "G99": -1.0, "G100": -2.0},
        {"UP": {"G1", "G2"}, "DOWN": {"G99", "G100"}},
    )
    assert running_es is not None
    assert set(pd.read_csv(running_es, sep="\t")["Term_ID"]) == {"UP", "DOWN"}


def test_plot_term_names_are_added_only_to_display_copy():
    from allenricher.cli import _with_plot_term_names

    official = _fgsea_frame().assign(pathway=["hsa00010", "hsa00020"])
    display = _with_plot_term_names(
        official,
        {
            "hsa00010": "Metabolism|Carbohydrate metabolism|Glycolysis / Gluconeogenesis",
            "hsa00020": "Metabolism|Carbohydrate metabolism|Citrate cycle (TCA cycle)",
        },
        "KEGG",
    )

    assert official.columns.tolist() == FGSEA_COLUMNS
    assert "Term_Name" not in official.columns
    assert display["Term_ID"].tolist() == ["hsa00010", "hsa00020"]
    assert display["Term_Name"].str.contains("Glycolysis|Citrate cycle").all()
    assert display["Database"].tolist() == ["KEGG", "KEGG"]


def test_python_gsea_barplot_uses_nes_and_pathway_names(tmp_path):
    import matplotlib.pyplot as plt

    from allenricher.cli import _with_plot_term_names
    from allenricher.visualization.gsea_plots import plot_gsea_barplot

    official = _fgsea_frame().assign(pathway=["hsa00010", "hsa00020"])
    display = _with_plot_term_names(
        official,
        {
            "hsa00010": "Metabolism|Carbohydrate metabolism|Glycolysis / Gluconeogenesis",
            "hsa00020": "Metabolism|Carbohydrate metabolism|Citrate cycle (TCA cycle)",
        },
        "KEGG",
    )
    output = tmp_path / "gsea_barplot.png"
    figure = plot_gsea_barplot(
        display, database="KEGG", output_file=str(output), dpi=100
    )
    labels = " ".join(text.get_text().replace("\n", " ") for text in figure.findobj(plt.Text))

    assert output.stat().st_size > 0
    assert "KEGG GSEA NES Ranking" in labels
    assert "Glycolysis / Gluconeogenesis" in labels
    assert "Citrate cycle (TCA cycle)" in labels
    assert "hsa00010" not in labels
    assert "Rich Factor" not in labels
    assert figure.axes[0].get_xlabel() == "Normalized enrichment score (NES)"
    plt.close(figure)


@pytest.mark.skipif(shutil.which("Rscript") is None, reason="Rscript is not installed")
def test_r_plot_consumer_accepts_fgsea_table(tmp_path):
    from allenricher.visualization.r_plotter import plot_gsea_barplot_r

    table = tmp_path / "fgsea.tsv"
    output = tmp_path / "barplot.png"
    _fgsea_frame().to_csv(table, sep="\t", index=False)
    assert plot_gsea_barplot_r(str(table), str(output), top_n=2, dpi=100)
    assert output.stat().st_size > 0


@pytest.mark.skipif(shutil.which("Rscript") is None, reason="Rscript is not installed")
def test_fgsea_runtime_returns_official_columns_and_signed_nes():
    ranked = [(f"G{i}", float(51 - i)) for i in range(1, 101)]
    gene_sets = {
        "UP": {f"G{i}" for i in range(1, 16)},
        "DOWN": {f"G{i}" for i in range(86, 101)},
    }
    result = run_fgsea(ranked, gene_sets, min_size=5, max_size=50)
    assert result.columns.tolist() == FGSEA_COLUMNS
    nes = result.set_index("pathway")["NES"]
    assert nes["UP"] > 0
    assert nes["DOWN"] < 0
