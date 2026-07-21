import pytest
from scipy.stats import hypergeom
from statsmodels.stats.multitest import multipletests

from allenricher.core.config import Config
from allenricher.core.enrichment import EnrichmentAnalyzer, GSEA


def _ora(**kwargs):
    config = Config(
        species="hsa",
        databases=["TEST"],
        method="hypergeometric",
        min_genes=1,
        output_all=True,
        **kwargs,
    )
    return EnrichmentAnalyzer(config)


def test_ora_uses_one_background_universe():
    results = _ora().run_analysis(
        {"A", "X", "Y"},
        {"A", "B", "C", "D"},
        {"TEST": {"T": {"name": "T", "genes": ["A", "X"]}}},
        parallel=False,
    )["TEST"]

    row = results.iloc[0]
    assert row["Gene_Count"] == 1
    assert row["Background_Count"] == 1
    assert row["Gene_Ratio"] == "1/1"
    assert row["P_Value"] == pytest.approx(hypergeom.sf(0, 4, 1, 1))


def test_ora_corrects_positive_overlap_hypotheses_like_v1():
    database = {
        "double": {"name": "double", "genes": ["A", "B"]},
        "single": {"name": "single", "genes": ["A", "C"]},
        "zero": {"name": "zero", "genes": ["D", "E"]},
    }
    results = _ora().run_analysis(
        {"A", "B"}, set("ABCDEFGHIJ"), {"TEST": database}, parallel=False
    )["TEST"]

    raw = [hypergeom.sf(1, 10, 2, 2), hypergeom.sf(0, 10, 2, 2)]
    expected = dict(zip(["double", "single"], multipletests(raw, method="fdr_bh")[1]))

    assert set(results["Term_ID"]) == {"double", "single"}
    for _, row in results.iterrows():
        assert row["Adjusted_P_Value"] == pytest.approx(expected[row["Term_ID"]])


def test_gsea_reports_signed_bottom_enrichment():
    ranked = [f"G{i}" for i in range(1, 11)]
    weights = {gene: 11 - index for index, gene in enumerate(ranked, 1)}

    es, leading_edge, rank = GSEA(min_size=1).calculate_enrichment_score(
        ranked, {"G8", "G9", "G10"}, weights, return_rank=True
    )

    assert es == pytest.approx(-1.0)
    assert leading_edge == ["G8", "G9", "G10"]
    assert rank == 7


def test_ranked_gene_loader_requires_real_finite_weights(tmp_path):
    ranked = tmp_path / "ranked.tsv"
    ranked.write_text("gene\tweight\nG1\t2.5\nG2\t-1.0\n", encoding="utf-8")

    assert EnrichmentAnalyzer.load_ranked_gene_list(str(ranked)) == [
        ("G1", 2.5),
        ("G2", -1.0),
    ]

    ranked.write_text("gene\nG1\nG2\n", encoding="utf-8")
    with pytest.raises(ValueError, match="gene.*weight"):
        EnrichmentAnalyzer.load_ranked_gene_list(str(ranked))

    ranked.write_text("gene\tweight\nG1\t1\nG2\tNaN\n", encoding="utf-8")
    with pytest.raises(ValueError, match="non-finite weights"):
        EnrichmentAnalyzer.load_ranked_gene_list(str(ranked))


def test_ranked_gene_loader_rejects_duplicate_genes_and_constant_scores(tmp_path):
    ranked = tmp_path / "ranked.tsv"
    ranked.write_text("gene\tweight\nG1\t2\nG1\t1\n", encoding="utf-8")
    with pytest.raises(ValueError, match="duplicate gene IDs"):
        EnrichmentAnalyzer.load_ranked_gene_list(str(ranked))

    ranked.write_text("gene\tweight\nG1\t1\nG2\t1\n", encoding="utf-8")
    with pytest.raises(ValueError, match="two distinct values"):
        EnrichmentAnalyzer.load_ranked_gene_list(str(ranked))
