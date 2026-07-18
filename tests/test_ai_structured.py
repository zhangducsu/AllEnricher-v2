"""Tests for method-aware AI evidence, validation, and report anchors."""

import json

import pandas as pd
import pytest

from allenricher.ai.interpreter import (
    build_structured_evidence,
    create_interpreter,
    validate_interpretation,
)
from allenricher.cli import create_parser
from allenricher.report.generator import ReportGenerator


def test_ora_evidence_keeps_statistics_and_genes():
    results = {
        "GO": pd.DataFrame({
            "Term_ID": ["GO:0001"],
            "Term_Name": ["Cell cycle"],
            "P_Value": [1e-5],
            "Adjusted_P_Value": [2e-4],
            "Gene_Count": [4],
            "Genes": ["A;B;C;D"],
            "EnrichFactor": [2.4],
        })
    }
    evidence = build_structured_evidence(results, method="hypergeometric")
    record = evidence["evidence"]["GO:R001"]
    assert record["values"]["adjusted_p_value"] == 2e-4
    assert record["values"]["enrich_factor"] == 2.4
    assert record["values"]["genes"] == ["A", "B", "C", "D"]
    assert record["raw"]["Term_Name"] == "Cell cycle"


def test_method_specific_evidence_limits_and_override():
    ora = pd.DataFrame({
        "Term_ID": [f"GO:{i:04d}" for i in range(20)],
        "Term_Name": [f"Term {i}" for i in range(20)],
        "P_Value": [0.001 + i / 10000 for i in range(20)],
        "Adjusted_P_Value": [0.002 + i / 10000 for i in range(20)],
        "Gene_Count": [5] * 20,
        "EnrichFactor": [1.0 + i / 100 for i in range(20)],
    })
    assert len(build_structured_evidence({"GO": ora}, "ora")["evidence"]) == 15
    assert len(build_structured_evidence({"GO": ora}, "ora", top_n=3)["evidence"]) == 3

    gsea = pd.DataFrame({
        "Term_ID": [f"R{i}" for i in range(24)],
        "Term_Name": [f"Pathway {i}" for i in range(24)],
        "pval": [0.01] * 24,
        "padj": [0.02] * 24,
        "ES": [1.0] * 12 + [-1.0] * 12,
        "NES": [2.0] * 12 + [-2.0] * 12,
        "size": [20] * 24,
        "leadingEdge": ["A;B"] * 24,
    })
    assert len(build_structured_evidence({"Reactome": gsea}, "gsea")["evidence"]) == 20
    assert len(build_structured_evidence({"Reactome": gsea}, "gsea", top_n=2)["evidence"]) == 4


def test_activity_evidence_uses_group_difference_not_sample_variance():
    activity = pd.DataFrame({
        "Term_ID": ["A", "B"],
        "Term_Name": ["Strong group shift", "Within-group variation"],
        "Control_1": [1.0, 0.0],
        "Control_2": [1.0, 10.0],
        "Disease_1": [-1.0, 5.0],
        "Disease_2": [-1.0, 5.0],
    })
    groups = {"Control": ["Control_1", "Control_2"], "Disease": ["Disease_1", "Disease_2"]}
    evidence = build_structured_evidence({"GO": activity}, "gsva", groups=groups, top_n=1)
    assert evidence["evidence"]["GSVA_GO:R001"]["term_id"] == "A"


def test_gsea_evidence_separates_nes_directions_and_leading_edge():
    results = {
        "Reactome": pd.DataFrame({
            "Term_ID": ["R1", "R2"],
            "Term_Name": ["Up pathway", "Down pathway"],
            "pval": [0.001, 0.002],
            "padj": [0.01, 0.02],
            "ES": [0.4, -0.3],
            "NES": [2.1, -1.8],
            "size": [30, 20],
            "leadingEdge": ["A;B", "C;D"],
        })
    }
    evidence = build_structured_evidence(results, method="gsea")
    records = list(evidence["evidence"].values())
    assert {record["values"]["direction"] for record in records} == {"positive", "negative"}
    assert records[0]["values"]["leading_edge_genes"]
    assert all(record["evidence_id"].startswith("GSEA_Reactome:R") for record in records)


def test_activity_evidence_computes_group_summary():
    results = {
        "GO": pd.DataFrame({
            "Term_ID": ["GO:1"],
            "Term_Name": ["Apoptosis"],
            "Control_1": [1.0],
            "Control_2": [1.2],
            "Disease_1": [-1.0],
            "Disease_2": [-1.1],
        })
    }
    groups = {"Control": ["Control_1", "Control_2"], "Disease": ["Disease_1", "Disease_2"]}
    evidence = build_structured_evidence(results, method="gsva", groups=groups)
    values = evidence["evidence"]["GSVA_GO:R001"]["values"]
    assert values["group_means"]["Control"] == pytest.approx(1.1)
    assert values["group_differences"]["Control vs Disease"] == pytest.approx(2.15)


@pytest.mark.parametrize("mode", ["summary", "reviewer", "caption"])
def test_mock_supports_all_structured_modes(mode):
    results = {"GO": pd.DataFrame({"Term_ID": ["GO:1"], "Term_Name": ["Signal"], "P_Value": [0.01], "Adjusted_P_Value": [0.02]})}
    output = create_interpreter("mock").interpret_structured_results(results, "hypergeometric", mode)
    assert output["profile"] == mode
    assert set(output["databases"]["GO"]) == {
        "core_themes", "key_evidence", "limitations", "computational_checks"
    }
    for item in output["databases"]["GO"]["key_evidence"]:
        assert set(item["evidence_ids"]) <= set(output["evidence"])


def test_invalid_json_and_unknown_evidence_are_rejected():
    results = {"GO": pd.DataFrame({"Term_ID": ["GO:1"], "Term_Name": ["Signal"], "P_Value": [0.01], "Adjusted_P_Value": [0.02]})}
    evidence = build_structured_evidence(results)
    with pytest.raises(ValueError, match="valid JSON"):
        validate_interpretation("not-json", evidence, "summary")
    invalid = {
        "schema_version": 1,
        "method": "ora",
        "profile": "summary",
        "databases": {
            "GO": {
                "core_themes": [{"text": "Unsupported claim", "evidence_ids": ["GO:R999"]}],
                "key_evidence": [], "limitations": [], "computational_checks": [],
            }
        },
    }
    with pytest.raises(ValueError, match="unknown evidence_id"):
        validate_interpretation(invalid, evidence, "summary")


def test_report_links_evidence_to_result_row_and_escapes_ai_text(tmp_path):
    results = {"GO": pd.DataFrame({"Term_ID": ["GO:1"], "Term_Name": ["Signal"], "P_Value": [0.01], "Adjusted_P_Value": [0.02]})}
    interpretation = create_interpreter("mock").interpret_structured_results(results, "hypergeometric", "summary")
    interpretation["databases"]["GO"]["core_themes"][0]["text"] = "<script>unsafe()</script>"
    output_file = tmp_path / "report.html"
    ReportGenerator(str(tmp_path)).generate(
        results, str(output_file), ai_interpretation=interpretation,
        analysis_method="hypergeometric", metadata={"database_version": "test"},
    )
    html = output_file.read_text(encoding="utf-8")
    assert 'id="evidence-GO-R001"' in html
    assert 'href="#evidence-GO-R001"' in html
    assert "&lt;script&gt;unsafe()&lt;/script&gt;" in html
    assert "<script>unsafe()</script>" not in html


def test_cli_exposes_ai_mode():
    args = create_parser().parse_args([
        "analyze", "-i", "genes.txt", "--ai", "mock", "--ai-mode", "reviewer", "--ai-top-n", "7"
    ])
    assert args.ai_mode == "reviewer"
    assert args.ai_top_n == 7
