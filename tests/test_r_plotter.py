"""R Draw bridge layer test."""

from pathlib import Path
import shutil
from unittest.mock import MagicMock, patch

import pytest

from allenricher.visualization import r_plotter


def _script_path(name: str) -> Path:
    return r_plotter.R_SCRIPTS_DIR / name


def _successful_r_run(cmd, **kwargs):
    output = Path(cmd[cmd.index("--output") + 1])
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_bytes(b"plot")
    return MagicMock(returncode=0, stdout="ok", stderr="")


def test_run_r_script_resolves_only_file_arguments(tmp_path):
    """All true file parameters should be turned in absolute ways, and the gene_set_ids keep a normal string."""
    tsv = tmp_path / "result.tsv"
    expr = tmp_path / "expr.tsv"
    running_es = tmp_path / "running_es.tsv"
    scores = tmp_path / "scores.tsv"
    metadata = tmp_path / "metadata.tsv"
    statistics = tmp_path / "statistics.tsv"
    output = tmp_path / "plot.png"
    for path in (tsv, expr, running_es, scores, metadata, statistics):
        path.write_text("x\n", encoding="utf-8")

    with patch("allenricher.visualization.r_plotter.subprocess.run", side_effect=_successful_r_run) as run:
        ok = r_plotter.run_r_script(
            "gsea_barplot.R",
            {
                "tsv": str(tsv),
                "expr": str(expr),
                "running_es": str(running_es),
                "scores": str(scores),
                "metadata": str(metadata),
                "statistics": str(statistics),
                "gene_set_ids": "TERM:1,TERM:2",
                "top_n": "20",
                "qvalue": "0.05",
                "min_count": "3",
            },
            str(output),
        )

    assert ok is True
    cmd = run.call_args.args[0]
    assert cmd[:2] == ["Rscript", str(_script_path("gsea_barplot.R"))]
    assert cmd[cmd.index("--tsv") + 1] == str(tsv.resolve())
    assert cmd[cmd.index("--expr") + 1] == str(expr.resolve())
    assert cmd[cmd.index("--running_es") + 1] == str(running_es.resolve())
    assert cmd[cmd.index("--scores") + 1] == str(scores.resolve())
    assert cmd[cmd.index("--metadata") + 1] == str(metadata.resolve())
    assert cmd[cmd.index("--statistics") + 1] == str(statistics.resolve())
    assert cmd[cmd.index("--gene_set_ids") + 1] == "TERM:1,TERM:2"
    assert cmd[cmd.index("--top_n") + 1] == "20"
    assert cmd[cmd.index("--qvalue") + 1] == "0.05"
    assert cmd[cmd.index("--min_count") + 1] == "3"
    if r_plotter.os.name == "nt":
        env = run.call_args.kwargs["env"]
        assert all(env.get(key) != "C.UTF-8" for key in ("LANG", "LC_ALL", "LC_CTYPE"))
    assert run.call_args.kwargs["encoding"] == "utf-8"
    assert run.call_args.kwargs["errors"] == "replace"


def test_run_r_script_failure_logs_stdout_and_stderr(tmp_path, caplog):
    """returns False when R fails and records stdout/stderr- It's easy to locate."""
    output = tmp_path / "plot.png"
    completed = MagicMock(returncode=1, stdout="partial output", stderr="missing package")

    with patch("allenricher.visualization.r_plotter.subprocess.run", return_value=completed):
        ok = r_plotter.run_r_script("gsea_barplot.R", {}, str(output))

    assert ok is False
    assert "partial output" in caplog.text
    assert "missing package" in caplog.text


def test_run_r_script_rejects_missing_output_after_success(tmp_path, caplog):
    completed = MagicMock(returncode=0, stdout="", stderr="")
    with patch("allenricher.visualization.r_plotter.subprocess.run", return_value=completed):
        ok = r_plotter.run_r_script("gsea_barplot.R", {}, str(tmp_path / "missing.png"))
    assert ok is False
    assert "did not create a non-empty output" in caplog.text


def test_run_r_script_falls_back_to_wsl_when_windows_r_is_blocked(tmp_path):
    output = tmp_path / "plot.png"
    completed = MagicMock(returncode=0, stdout="", stderr="")

    def run(cmd, **kwargs):
        if cmd[0] == "Rscript":
            raise OSError(4551, "blocked")
        output.write_bytes(b"plot")
        return completed

    with patch("allenricher.visualization.r_plotter.subprocess.run", side_effect=run) as mocked:
        ok = r_plotter.run_r_script("gsea_barplot.R", {}, str(output))

    assert ok is True
    assert mocked.call_count == 2
    fallback = mocked.call_args_list[1].args[0]
    assert Path(fallback[0]).name.lower() in {"wsl", "wsl.exe"}
    assert fallback[1] == "Rscript"


def test_run_r_script_falls_back_when_windows_runtime_dll_is_blocked(tmp_path):
    if r_plotter.os.name != "nt":
        pytest.skip("Windows-specific fallback")
    output = tmp_path / "plot.png"

    def run(cmd, **kwargs):
        if cmd[0] == "Rscript":
            return MagicMock(
                returncode=1,
                stdout="",
                stderr=(
                    "package 'grDevices' in options(\"defaultPackages\") was not found\n"
                    "LoadLibrary failure: application control policy has blocked this file"
                ),
            )
        output.write_bytes(b"plot")
        return MagicMock(returncode=0, stdout="", stderr="")

    with patch("allenricher.visualization.r_plotter.subprocess.run", side_effect=run) as mocked:
        ok = r_plotter.run_r_script("gsea_barplot.R", {}, str(output))

    assert ok is True
    assert mocked.call_count == 2


def test_enrichment2_keeps_gene_set_ids_as_ids(tmp_path):
    """Gene_set_ids of enrichment2 should not be mistransmitted by Path.resolve()."""
    output = tmp_path / "plot.png"
    with patch("allenricher.visualization.r_plotter.subprocess.run", side_effect=_successful_r_run) as run:
        ok = r_plotter.plot_gsea_enrichment2_r(
            "result.tsv",
            ["GO:0001", "hsa00010"],
            str(output),
            running_es_path="running.tsv",
        )

    assert ok is True
    cmd = run.call_args.args[0]
    assert cmd[cmd.index("--gene_set_ids") + 1] == "GO:0001,hsa00010"
    assert cmd[cmd.index("--running_es") + 1].endswith("running.tsv")


def test_ridgeplot_accepts_running_es_file(tmp_path):
    """Riddlegeplot should be able to receive the real running-ES middle table."""
    output = tmp_path / "plot.png"
    with patch("allenricher.visualization.r_plotter.subprocess.run", side_effect=_successful_r_run) as run:
        ok = r_plotter.plot_gsea_ridgeplot_r(
            "result.tsv",
            str(output),
            top_n=15,
            running_es_path="running.tsv",
        )

    assert ok is True
    cmd = run.call_args.args[0]
    assert cmd[cmd.index("--running_es") + 1].endswith("running.tsv")
    assert cmd[cmd.index("--top_n") + 1] == "15"


def test_lollipop_targets_new_r_script(tmp_path):
    """Lollypop should call for a new gsea_lollipop.R script."""
    output = tmp_path / "plot.png"
    with patch("allenricher.visualization.r_plotter.subprocess.run", side_effect=_successful_r_run) as run:
        ok = r_plotter.plot_gsea_lollipop_r(
            "result.tsv", str(output), top_n=12,
            style="science", palette="tol_sunset",
        )

    assert ok is True
    cmd = run.call_args.args[0]
    assert cmd[:2] == ["Rscript", str(_script_path("gsea_lollipop.R"))]
    assert cmd[cmd.index("--top_n") + 1] == "12"
    assert cmd[cmd.index("--dpi") + 1] == "300"
    assert cmd[cmd.index("--style") + 1] == "science"
    assert cmd[cmd.index("--categorical_palette") + 1] == "tol_bright"
    assert cmd[cmd.index("--sequential_palette") + 1] == "colorbrewer_blues"
    assert cmd[cmd.index("--diverging_palette") + 1] == "tol_sunset"


def test_r_bridge_canonicalizes_legacy_styles(tmp_path):
    output = tmp_path / "plot.png"
    with patch("allenricher.visualization.r_plotter.subprocess.run", side_effect=_successful_r_run) as run:
        assert r_plotter.plot_gsea_lollipop_r("result.tsv", str(output), style="cell")
    cmd = run.call_args.args[0]
    assert cmd[cmd.index("--style") + 1] == "nature"


def test_r_bridge_passes_all_explicit_palette_roles(tmp_path):
    from allenricher.visualization.color_config import PaletteSelection

    output = tmp_path / "plot.png"
    selection = PaletteSelection(
        categorical="nature",
        sequential="viridis",
        diverging="colorbrewer_brbg",
    )
    with patch("allenricher.visualization.r_plotter.subprocess.run", side_effect=_successful_r_run) as run:
        ok = r_plotter.plot_gsea_lollipop_r(
            "result.tsv", str(output), palette=selection,
        )

    assert ok is True
    cmd = run.call_args.args[0]
    assert cmd[cmd.index("--categorical_palette") + 1] == "nature"
    assert cmd[cmd.index("--sequential_palette") + 1] == "viridis"
    assert cmd[cmd.index("--diverging_palette") + 1] == "colorbrewer_brbg"


def test_activity_heatmap_passes_analysis_method(tmp_path):
    output = tmp_path / "plot.png"
    with patch("allenricher.visualization.r_plotter.subprocess.run", side_effect=_successful_r_run) as run:
        ok = r_plotter.plot_activity_heatmap_r(
            "scores.tsv", "metadata.tsv", str(output), analysis_method="gsva"
        )

    assert ok is True
    cmd = run.call_args.args[0]
    assert cmd[cmd.index("--analysis_method") + 1] == "gsva"
    assert cmd[cmd.index("--top_n") + 1] == "40"


@pytest.mark.skipif(shutil.which("Rscript") is None, reason="Rscript is not installed")
def test_sample_correlation_r_accepts_duplicate_pathway_names(tmp_path):
    scores = tmp_path / "scores.tsv"
    metadata = tmp_path / "metadata.tsv"
    output = tmp_path / "correlation.png"
    scores.write_text(
        "Pathway\tS1\tS2\tS3\n"
        "Shared description\t1\t2\t4\n"
        "Shared description\t2\t1\t3\n"
        "Distinct description\t3\t4\t1\n",
        encoding="utf-8",
    )
    metadata.write_text(
        "Sample\tGroup\nS1\tControl\nS2\tControl\nS3\tTreatment\n",
        encoding="utf-8",
    )

    assert r_plotter.plot_sample_correlation_r(
        str(scores), str(metadata), str(output), dpi=72
    )
    assert output.stat().st_size > 0


def test_emapplot_passes_filter_args(tmp_path):
    """Emplot should pass to R scripts top_n, qvalue and min_count."""
    output = tmp_path / "plot.png"
    with patch("allenricher.visualization.r_plotter.subprocess.run", side_effect=_successful_r_run) as run:
        ok = r_plotter.plot_gsea_emapplot_r(
            "result.tsv",
            str(output),
            top_n=12,
            qvalue=0.01,
            min_count=5,
        )

    assert ok is True
    cmd = run.call_args.args[0]
    assert cmd[:2] == ["Rscript", str(_script_path("gsea_emapplot.R"))]
    assert cmd[cmd.index("--top_n") + 1] == "12"
    assert cmd[cmd.index("--qvalue") + 1] == "0.01"
    assert cmd[cmd.index("--min_count") + 1] == "5"
