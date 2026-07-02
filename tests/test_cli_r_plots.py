"""CLI 的 R GSEA 绘图集成测试。"""

from pathlib import Path

import pandas as pd

from allenricher.cli import (
    _METHOD_PLOT_TYPES,
    _generate_plots,
    _normalize_ranked_genes,
    _write_running_es_file,
    create_parser,
)


def test_gsea_r_plot_types_include_heatmap_and_enrichment2():
    """GSEA 图表白名单应包含已接入的 R 图表类型。"""
    assert "heatmap" in _METHOD_PLOT_TYPES["gsea"]
    assert "enrichment2" in _METHOD_PLOT_TYPES["gsea"]


def test_use_r_plots_option():
    """验证 --use-r-plots 参数可以被正确解析。"""
    parser = create_parser()
    args = parser.parse_args([
        "analyze", "-i", "genes.txt", "-m", "gsea", "--use-r-plots",
        "-pt", "enrichment,enrichment2,heatmap",
    ])
    assert args.use_r_plots is True
    assert args.plot_types == "enrichment,enrichment2,heatmap"


def test_normalize_ranked_genes_accepts_tuple_list():
    """排序基因 tuple 列表应转成 gene 列表和权重字典。"""
    genes, weights = _normalize_ranked_genes([("A", 2.0), ("B", -1.5)], None)
    assert genes == ["A", "B"]
    assert weights == {"A": 2.0, "B": -1.5}


def test_write_running_es_file_uses_real_ranked_data(tmp_path):
    """running ES 中间表应来自真实排序基因和基因集。"""
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


def test_gsea_r_heatmap_receives_result_tsv(tmp_path, monkeypatch):
    """CLI R heatmap 应拿到 GSEA 结果 TSV，用于筛选通路相关基因。"""
    calls = {}

    def fake_heatmap(expr_path, output_file, tsv_path="", top_n=12):
        calls["expr_path"] = expr_path
        calls["output_file"] = output_file
        calls["tsv_path"] = tsv_path
        calls["top_n"] = top_n
        return True

    monkeypatch.setattr("allenricher.visualization.r_plotter.check_r_environment", lambda: True)
    monkeypatch.setattr("allenricher.visualization.r_plotter.plot_gsea_heatmap_r", fake_heatmap)

    results = {
        "KEGG": pd.DataFrame({
            "Term_ID": ["TERM1"],
            "Term_Name": ["Pathway 1"],
            "NES": [1.8],
            "Genes": ["A;B"],
        })
    }
    expr = pd.DataFrame({"S1": [1.0, 2.0], "S2": [2.0, 1.0]}, index=["A", "B"])

    generated = _generate_plots(
        method="gsea",
        results=results,
        ranked_genes=None,
        gene_weights=None,
        gene_sets=None,
        expr_matrix=expr,
        groups=None,
        plot_types=["heatmap"],
        output_dir=str(tmp_path),
        use_r_plots=True,
    )

    assert len(generated) == 1
    assert calls["top_n"] == 12
    assert calls["tsv_path"].endswith("KEGG_enrichment.tsv")
    assert calls["expr_path"].endswith("KEGG_expression_matrix.tsv")
    assert Path(calls["tsv_path"]).exists()
