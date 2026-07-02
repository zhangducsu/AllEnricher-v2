"""R 绘图桥接层测试。"""

from pathlib import Path
from unittest.mock import MagicMock, patch

from allenricher.visualization import r_plotter


def _script_path(name: str) -> Path:
    return r_plotter.R_SCRIPTS_DIR / name


def test_run_r_script_resolves_only_file_arguments(tmp_path):
    """tsv/expr/running_es 应转绝对路径，gene_set_ids 保持普通字符串。"""
    tsv = tmp_path / "result.tsv"
    expr = tmp_path / "expr.tsv"
    running_es = tmp_path / "running_es.tsv"
    output = tmp_path / "plot.png"
    for path in (tsv, expr, running_es):
        path.write_text("x\n", encoding="utf-8")

    completed = MagicMock(returncode=0, stdout="ok", stderr="")
    with patch("allenricher.visualization.r_plotter.subprocess.run", return_value=completed) as run:
        ok = r_plotter.run_r_script(
            "gsea_dotplot.R",
            {
                "tsv": str(tsv),
                "expr": str(expr),
                "running_es": str(running_es),
                "gene_set_ids": "TERM:1,TERM:2",
                "top_n": "20",
            },
            str(output),
        )

    assert ok is True
    cmd = run.call_args.args[0]
    assert cmd[:2] == ["Rscript", str(_script_path("gsea_dotplot.R"))]
    assert cmd[cmd.index("--tsv") + 1] == str(tsv.resolve())
    assert cmd[cmd.index("--expr") + 1] == str(expr.resolve())
    assert cmd[cmd.index("--running_es") + 1] == str(running_es.resolve())
    assert cmd[cmd.index("--gene_set_ids") + 1] == "TERM:1,TERM:2"
    assert cmd[cmd.index("--top_n") + 1] == "20"


def test_run_r_script_failure_logs_stdout_and_stderr(tmp_path, caplog):
    """R 失败时返回 False，并记录 stdout/stderr 便于定位问题。"""
    output = tmp_path / "plot.png"
    completed = MagicMock(returncode=1, stdout="partial output", stderr="missing package")

    with patch("allenricher.visualization.r_plotter.subprocess.run", return_value=completed):
        ok = r_plotter.run_r_script("gsea_dotplot.R", {}, str(output))

    assert ok is False
    assert "partial output" in caplog.text
    assert "missing package" in caplog.text


def test_enrichment2_keeps_gene_set_ids_as_ids(tmp_path):
    """enrichment2 的 gene_set_ids 不应被 Path.resolve() 误转成路径。"""
    output = tmp_path / "plot.png"
    completed = MagicMock(returncode=0, stdout="", stderr="")

    with patch("allenricher.visualization.r_plotter.subprocess.run", return_value=completed) as run:
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
    """ridgeplot 应能接收真实 running-ES 中间表。"""
    output = tmp_path / "plot.png"
    completed = MagicMock(returncode=0, stdout="", stderr="")

    with patch("allenricher.visualization.r_plotter.subprocess.run", return_value=completed) as run:
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


def test_heatmap_receives_gsea_tsv_for_pathway_gene_filtering(tmp_path):
    """heatmap 应接收 GSEA TSV，用于优先筛选 top 通路命中基因。"""
    expr = tmp_path / "expr.tsv"
    tsv = tmp_path / "gsea.tsv"
    output = tmp_path / "heatmap.png"
    expr.write_text("gene\tS1\nA\t1\n", encoding="utf-8")
    tsv.write_text("Term_ID\tTerm_Name\tGenes\tNES\nT1\tPathway\tA\t1.2\n", encoding="utf-8")
    completed = MagicMock(returncode=0, stdout="", stderr="")

    with patch("allenricher.visualization.r_plotter.subprocess.run", return_value=completed) as run:
        ok = r_plotter.plot_gsea_heatmap_r(str(expr), str(output), tsv_path=str(tsv), top_n=12)

    assert ok is True
    cmd = run.call_args.args[0]
    assert cmd[cmd.index("--expr") + 1] == str(expr.resolve())
    assert cmd[cmd.index("--tsv") + 1] == str(tsv.resolve())
    assert cmd[cmd.index("--top_n") + 1] == "12"
