"""Tests for method-aware AI evidence, validation, and report anchors."""

import json

import pandas as pd
import pytest

from allenricher.ai.interpreter import (
    AIInterpreter,
    build_interpretation_prompt,
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
    assert values["inference_scope"] == "descriptive_activity_pattern"


def test_evidence_relations_distinguish_shared_core_redundancy_and_conflict():
    ora = pd.DataFrame({
        "Term_ID": ["A", "B", "C"],
        "Term_Name": ["Cell cycle checkpoint", "Cell cycle transition", "Cell cycle checkpoint duplicate"],
        "Adjusted_P_Value": [0.001, 0.002, 0.003],
        "EnrichFactor": [3.0, 2.5, 2.0],
        "Genes": ["A;B;C;D", "A;B;C;E", "A;B;C;D"],
    })
    evidence = build_structured_evidence({"GO": ora}, "ora")
    relation_types = {relation["relation_type"] for relation in evidence["relations"]}
    assert "shared_core" in relation_types
    assert "redundant" in relation_types
    assert any(set(relation["shared_genes"]) >= {"A", "B", "C"} for relation in evidence["relations"])

    gsea = pd.DataFrame({
        "Term_ID": ["UP", "DOWN"],
        "Term_Name": ["Cell cycle activation", "Cell cycle suppression"],
        "padj": [0.01, 0.02],
        "NES": [2.0, -1.8],
        "leadingEdge": ["CDK1;CCNB1", "CDK1;CCNB1"],
    })
    gsea_evidence = build_structured_evidence({"Reactome": gsea}, "gsea")
    assert gsea_evidence["relations"][0]["relation_type"] == "conflicting"
    assert gsea_evidence["relations"][0]["direction_consistent"] is False

    broad_hierarchy = pd.DataFrame({
        "Term_ID": ["H1", "H2"],
        "Term_Name": ["Lipid transport", "DNA repair"],
        "Hierarchy": ["Biology|Cellular processes", "Biology|Cellular processes"],
        "Adjusted_P_Value": [0.01, 0.02],
    })
    assert build_structured_evidence({"GO": broad_hierarchy}, "ora")["relations"] == []


def test_prompt_requests_research_synthesis_and_activity_caution():
    evidence = build_structured_evidence({
        "GO": pd.DataFrame({
            "Term_ID": ["GO:1"], "Term_Name": ["Signal"],
            "P_Value": [0.01], "Adjusted_P_Value": [0.02],
        })
    })
    prompt = build_interpretation_prompt(evidence)
    assert "recurring biological themes" in prompt
    assert "research_summary is the first-read take-home block" in prompt
    assert "strongest supported program, exploratory signals" in prompt
    assert "disease- or infection-named pathway terms" in prompt
    assert "biological_meaning" in prompt
    assert "what the cited pathways, diseases, or TF families usually represent" in prompt
    assert "adjusted P/FDR/padj when supplied" in prompt
    assert "describe the signal as exploratory or low confidence" in prompt
    assert "metabolic and growth-control pathways contribute to the positive GSEA evidence" in prompt
    assert "Do not put generic method cautions in research_summary" in prompt
    assert "Do not write tautologies" in prompt
    assert "separate term meaning from analysis implication" in prompt
    assert "Do not repeat the same biological wording" in prompt
    assert "Do not count redundant terms as independent corroboration" in prompt
    assert "do not claim statistical significance" in prompt
    assert 'use only the phrases "positive target-set enrichment"' in prompt


def test_prompt_omits_raw_rows_sample_values_and_long_gene_lists():
    genes = [f"G{i}" for i in range(30)]
    results = {"GO": pd.DataFrame({
        "Term_ID": ["GO:1"],
        "Term_Name": ["Signal"],
        "Hierarchy": ["Process|Signal"],
        "P_Value": [0.01],
        "Adjusted_P_Value": [0.02],
        "Genes": [";".join(genes)],
        "Unused_Trace_Field": ["must not enter the prompt"],
    })}
    evidence = build_structured_evidence(results)
    prompt = build_interpretation_prompt(evidence)
    prompt_payload = json.loads(prompt.split("Structured evidence:\n", 1)[1])
    record = prompt_payload["evidence"]["GO:R001"]
    assert "raw" not in record
    assert "Unused_Trace_Field" not in json.dumps(prompt_payload)
    assert len(record["values"]["genes"]) == 20
    assert record["values"]["genes_total"] == 30


@pytest.mark.parametrize("mode", ["summary", "reviewer", "caption"])
def test_mock_supports_all_structured_modes(mode):
    results = {"GO": pd.DataFrame({"Term_ID": ["GO:1"], "Term_Name": ["Signal"], "P_Value": [0.01], "Adjusted_P_Value": [0.02]})}
    output = create_interpreter("mock").interpret_structured_results(results, "hypergeometric", mode)
    assert output["profile"] == mode
    assert output["research_summary"]
    assert output["overall_synthesis"]
    assert output["overall_synthesis"][0]["support_class"] == "single_signal"
    assert output["overall_synthesis"][0]["confidence"] == "exploratory"
    assert set(output["databases"]["GO"]) == {
        "core_themes", "biological_meaning", "key_evidence", "limitations", "computational_checks"
    }
    for item in output["databases"]["GO"]["key_evidence"]:
        assert set(item["evidence_ids"]) <= set(output["evidence"])


def test_validator_rejects_false_convergence_claim():
    results = {"GO": pd.DataFrame({
        "Term_ID": ["GO:1"], "Term_Name": ["Signal"],
        "P_Value": [0.01], "Adjusted_P_Value": [0.02],
    })}
    evidence = build_structured_evidence(results)
    payload = {
        "schema_version": 1,
        "method": "ora",
        "profile": "summary",
        "overall_synthesis": [{
            "text": "Multiple pathways converge.",
            "evidence_ids": ["GO:R001"],
            "support_class": "convergent",
        }],
        "databases": {
            "GO": {
                "core_themes": [{
                    "text": "A single signal is available.",
                    "evidence_ids": ["GO:R001"],
                    "support_class": "single_signal",
                }],
                "biological_meaning": [], "key_evidence": [], "limitations": [], "computational_checks": [],
            }
        },
    }
    with pytest.raises(ValueError, match="not supported by the cited evidence relations"):
        validate_interpretation(payload, evidence, "summary")


def test_validator_accepts_multi_record_semantic_convergence_with_moderate_confidence():
    results = {"GO": pd.DataFrame({
        "Term_ID": ["GO:1", "GO:2"],
        "Term_Name": ["Mitotic spindle", "Chromosome segregation"],
        "P_Value": [0.001, 0.002],
        "Adjusted_P_Value": [0.01, 0.02],
    })}
    evidence = build_structured_evidence(results)
    ids = list(evidence["evidence"])
    payload = {
        "schema_version": 1,
        "method": "ora",
        "profile": "summary",
        "overall_synthesis": [{
            "text": "Mitotic chromosome organization is repeatedly represented.",
            "evidence_ids": ids,
            "support_class": "convergent",
        }],
        "databases": {
            "GO": {
                "core_themes": [{
                    "text": "Mitotic chromosome organization is repeatedly represented.",
                    "evidence_ids": ids,
                    "support_class": "convergent",
                }],
                "biological_meaning": [], "key_evidence": [], "limitations": [], "computational_checks": [],
            }
        },
    }
    interpretation = validate_interpretation(payload, evidence, "summary")
    theme = interpretation["overall_synthesis"][0]
    assert theme["relationship_basis"] == "biological_semantics"
    assert theme["confidence"] == "moderate"


def test_validator_accepts_semantic_conflict_only_when_directions_are_mixed():
    results = {"Reactome": pd.DataFrame({
        "Term_ID": ["UP", "DOWN"],
        "Term_Name": ["Mitotic spindle", "Chromosome segregation"],
        "pval": [0.001, 0.002], "padj": [0.01, 0.02],
        "NES": [2.0, -1.8],
    })}
    evidence = build_structured_evidence(results, "gsea")
    assert evidence["relations"] == []
    ids = list(evidence["evidence"])
    item = {
        "text": "Related mitotic signals have opposite enrichment directions.",
        "evidence_ids": ids,
        "support_class": "conflicting",
    }
    payload = {
        "schema_version": 1, "method": "gsea", "profile": "summary",
        "overall_synthesis": [item],
        "databases": {"Reactome": {
            "core_themes": [item], "biological_meaning": [], "key_evidence": [],
            "limitations": [], "computational_checks": [],
        }},
    }
    theme = validate_interpretation(payload, evidence, "summary")["overall_synthesis"][0]
    assert theme["direction"] == "mixed"
    assert theme["relationship_basis"] == "biological_semantics"
    assert theme["confidence"] == "exploratory"


def test_validator_normalizes_gsea_direction_overclaim_outside_limitations():
    results = {"KEGG": pd.DataFrame({
        "Term_ID": ["hsa00010"],
        "Term_Name": ["Glycolysis"],
        "pval": [0.001], "padj": [0.01], "NES": [1.8],
    })}
    evidence = build_structured_evidence(results, "gsea")
    item = {
        "text": "Glycolysis upregulation is suggested.",
        "evidence_ids": ["GSEA_KEGG:R001"],
        "support_class": "single_signal",
    }
    payload = {
        "schema_version": 1, "method": "gsea", "profile": "summary",
        "overall_synthesis": [item],
        "databases": {"KEGG": {
            "core_themes": [], "biological_meaning": [], "key_evidence": [],
            "limitations": [{
                "text": "GSEA direction reflects ranked-list enrichment, not pathway activation.",
                "evidence_ids": ["GSEA_KEGG:R001"],
            }],
            "computational_checks": [],
        }},
    }
    output = validate_interpretation(payload, evidence, "summary")
    theme = output["overall_synthesis"][0]
    assert theme["text"] == "Glycolysis positive GSEA pattern is suggested."
    assert theme["auto_normalized_claims"] == ["gsea_direction_wording"]


def test_validator_normalizes_gsea_expression_overclaim_in_research_summary():
    results = {"KEGG": pd.DataFrame({
        "Term_ID": ["hsa04060"],
        "Term_Name": ["Cytokine signaling"],
        "pval": [0.001], "padj": [0.01], "NES": [-1.8],
    })}
    evidence = build_structured_evidence(results, "gsea")
    payload = {
        "schema_version": 1, "method": "gsea", "profile": "summary",
        "research_summary": [{
            "text": "Negative enrichment indicates lower expression of these genes.",
            "evidence_ids": ["GSEA_KEGG:R001"],
        }],
        "overall_synthesis": [],
        "databases": {"KEGG": {
            "core_themes": [], "biological_meaning": [], "key_evidence": [],
            "limitations": [], "computational_checks": [],
        }},
    }
    summary = validate_interpretation(payload, evidence, "summary")["research_summary"][0]
    assert summary["text"] == "Negative enrichment indicates the cited genes support a negative GSEA pattern."
    assert summary["auto_normalized_claims"] == ["gsea_direction_wording"]


def test_validator_preserves_gsea_activity_caution_wording():
    results = {"KEGG": pd.DataFrame({
        "Term_ID": ["hsa04060"],
        "Term_Name": ["Cytokine signaling"],
        "pval": [0.001], "padj": [0.01], "NES": [-1.8],
    })}
    evidence = build_structured_evidence(results, "gsea")
    payload = {
        "schema_version": 1, "method": "gsea", "profile": "summary",
        "research_summary": [],
        "overall_synthesis": [],
        "databases": {"KEGG": {
            "core_themes": [], "biological_meaning": [], "key_evidence": [],
            "limitations": [{
                "text": "Do not infer pathway activation or repression from these GSEA results.",
                "evidence_ids": ["GSEA_KEGG:R001"],
            }],
            "computational_checks": [],
        }},
    }
    limitation = validate_interpretation(payload, evidence, "summary")["databases"]["KEGG"]["limitations"][0]
    assert limitation["text"] == "Do not infer pathway activation or repression from these GSEA results."


def test_validator_normalizes_gsea_biological_shift_overclaim():
    results = {"KEGG": pd.DataFrame({
        "Term_ID": ["hsa04152"],
        "Term_Name": ["AMPK signaling"],
        "pval": [0.001], "padj": [0.01], "NES": [1.8],
    })}
    evidence = build_structured_evidence(results, "gsea")
    payload = {
        "schema_version": 1, "method": "gsea", "profile": "summary",
        "research_summary": [{
            "text": "Positive enrichment suggests a coordinated shift toward metabolic programs.",
            "evidence_ids": ["GSEA_KEGG:R001"],
        }],
        "overall_synthesis": [],
        "databases": {"KEGG": {
            "core_themes": [], "biological_meaning": [], "key_evidence": [],
            "limitations": [], "computational_checks": [],
        }},
    }
    summary = validate_interpretation(payload, evidence, "summary")["research_summary"][0]
    assert "GSEA pattern involving metabolic programs" in summary["text"]
    assert summary["auto_normalized_claims"] == ["gsea_direction_wording"]


def test_validator_downgrades_non_significant_adjusted_confidence():
    results = {"KEGG": pd.DataFrame({
        "Term_ID": ["A", "B"],
        "Term_Name": ["Insulin signaling", "AMPK signaling"],
        "pval": [0.001, 0.002], "padj": [0.21, 0.31],
        "NES": [1.8, 1.7],
    })}
    evidence = build_structured_evidence(results, "gsea")
    ids = list(evidence["evidence"])
    payload = {
        "schema_version": 1, "method": "gsea", "profile": "summary",
        "research_summary": [{
            "text": "Metabolic signaling has Confidence: moderate due to nominal P-values.",
            "evidence_ids": ids,
        }],
        "overall_synthesis": [{
            "text": "Metabolic signaling is represented.",
            "evidence_ids": ids,
            "support_class": "convergent",
        }],
        "databases": {"KEGG": {
            "core_themes": [{
                "text": "Metabolic signaling is represented.",
                "evidence_ids": ids,
                "support_class": "convergent",
            }],
            "biological_meaning": [], "key_evidence": [],
            "limitations": [], "computational_checks": [],
        }},
    }
    output = validate_interpretation(payload, evidence, "summary")
    summary = output["research_summary"][0]
    assert "Confidence: low because adjusted P/FDR values are not significant" in summary["text"]
    assert summary["auto_normalized_claims"] == ["non_significant_confidence"]
    assert output["overall_synthesis"][0]["confidence"] == "exploratory"
    assert output["databases"]["KEGG"]["core_themes"][0]["confidence"] == "exploratory"


def test_validator_removes_generic_method_caution_from_research_summary():
    results = {"KEGG": pd.DataFrame({
        "Term_ID": ["hsa04152"],
        "Term_Name": ["AMPK signaling"],
        "pval": [0.001], "padj": [0.01], "NES": [1.8],
    })}
    evidence = build_structured_evidence(results, "gsea")
    payload = {
        "schema_version": 1, "method": "gsea", "profile": "summary",
        "research_summary": [{
            "text": "Do not infer pathway activity or regulatory direction from these GSEA results.",
            "evidence_ids": ["GSEA_KEGG:R001"],
        }],
        "overall_synthesis": [],
        "databases": {"KEGG": {
            "core_themes": [], "biological_meaning": [], "key_evidence": [],
            "limitations": [{
                "text": "GSEA direction indicates where pathway genes fall in the ranked list.",
                "evidence_ids": ["GSEA_KEGG:R001"],
            }],
            "computational_checks": [],
        }},
    }
    output = validate_interpretation(payload, evidence, "summary")
    assert output["research_summary"] == []


def test_validator_rewrites_gsea_ranked_list_tautology():
    results = {"KEGG": pd.DataFrame({
        "Term_ID": ["hsa04152"],
        "Term_Name": ["AMPK signaling"],
        "pval": [0.001], "padj": [0.01], "NES": [1.8],
    })}
    evidence = build_structured_evidence(results, "gsea")
    payload = {
        "schema_version": 1, "method": "gsea", "profile": "summary",
        "research_summary": [{
            "text": "The strongest supported program is a coordinated enrichment of genes near the top of the ranked list for metabolic signaling.",
            "evidence_ids": ["GSEA_KEGG:R001"],
        }],
        "overall_synthesis": [],
        "databases": {"KEGG": {
            "core_themes": [], "key_evidence": [],
            "biological_meaning": [{
                "text": "AMPK pathways regulate energy balance. Their top-ranked enrichment suggests genes in these pathways are coordinately positioned near the top of the ranked list. Together, these terms indicate that this biological program contributes to the top-ranked signal.",
                "evidence_ids": ["GSEA_KEGG:R001"],
            }],
            "limitations": [], "computational_checks": [],
        }},
    }
    meaning = validate_interpretation(payload, evidence, "summary")["databases"]["KEGG"]["biological_meaning"][0]
    assert "coordinately positioned near the top" not in meaning["text"]
    assert "top-ranked signal" not in meaning["text"]
    assert "support the selected positive GSEA pattern" in meaning["text"]
    assert meaning["auto_normalized_claims"] == ["gsea_tautology"]
    summary = validate_interpretation(payload, evidence, "summary")["research_summary"][0]
    assert summary["text"].startswith("The strongest supported program involves metabolic signaling.")

    payload["research_summary"][0]["text"] = (
        "A separate exploratory signal shows genes from immune pathways concentrated near the bottom of the ranked list."
    )
    bottom_summary = validate_interpretation(payload, evidence, "summary")["research_summary"][0]
    assert bottom_summary["text"] == "A separate exploratory signal involves immune pathways."


def test_validator_links_named_tf_evidence_and_normalizes_tf_activity_overclaim():
    results = {"ChEA3": pd.DataFrame({
        "Term_ID": ["FOXO6_SET", "FOXO1_SET"],
        "Term_Name": ["FOXO6 [Coexpression]", "FOXO1 [Coexpression]"],
        "pval": [0.001, 0.002], "padj": [0.01, 0.02],
        "NES": [-1.8, 1.9],
    })}
    evidence = build_structured_evidence(results, "gsea")
    foxo6_id = next(
        evidence_id for evidence_id, record in evidence["evidence"].items()
        if record["term_name"].startswith("FOXO6")
    )
    foxo1_id = next(
        evidence_id for evidence_id, record in evidence["evidence"].items()
        if record["term_name"].startswith("FOXO1")
    )

    def payload(text, support_class="single_signal"):
        item = {"text": text, "evidence_ids": [foxo6_id], "support_class": support_class}
        return {
            "schema_version": 1, "method": "gsea", "profile": "summary",
            "overall_synthesis": [item],
            "databases": {"ChEA3": {
                "core_themes": [item], "biological_meaning": [], "key_evidence": [],
                "limitations": [], "computational_checks": [],
            }},
        }

    interpretation = validate_interpretation(
        payload("FOXO6 contrasts with the positive FOXO1 signal.", "conflicting"),
        evidence,
        "summary",
    )
    theme = interpretation["overall_synthesis"][0]
    assert theme["auto_linked_evidence_ids"] == [foxo1_id]
    assert set(theme["evidence_ids"]) == {foxo6_id, foxo1_id}
    normalized = validate_interpretation(
        payload("FOXO6 target enrichment indicates TF repression."), evidence, "summary"
    )
    assert "transcriptional regulation" in normalized["overall_synthesis"][0]["text"]
    assert normalized["overall_synthesis"][0]["auto_normalized_claims"] == ["gsea_direction_wording"]
    bad_meaning = payload("FOXO6 target set is negative.", "single_signal")
    bad_meaning["databases"]["ChEA3"]["biological_meaning"] = [{
        "text": "FOXO6 enrichment indicates TF repression.",
        "evidence_ids": [foxo6_id],
    }]
    normalized = validate_interpretation(bad_meaning, evidence, "summary")
    meaning = normalized["databases"]["ChEA3"]["biological_meaning"][0]
    assert "transcriptional regulation" in meaning["text"]
    assert meaning["auto_normalized_claims"] == ["gsea_direction_wording"]


def test_invalid_json_and_unknown_evidence_are_rejected():
    results = {"GO": pd.DataFrame({"Term_ID": ["GO:1"], "Term_Name": ["Signal"], "P_Value": [0.01], "Adjusted_P_Value": [0.02]})}
    evidence = build_structured_evidence(results)
    with pytest.raises(ValueError, match="valid JSON"):
        validate_interpretation("not-json", evidence, "summary")
    wrapped = (
        "Here is the corrected JSON:\n"
        + json.dumps({
            "schema_version": 1, "method": "ora", "profile": "summary",
            "overall_synthesis": [],
            "databases": {"GO": {
                "core_themes": [], "biological_meaning": [], "key_evidence": [],
                "limitations": [], "computational_checks": [],
            }},
        })
        + "\nDone."
    )
    assert validate_interpretation(wrapped, evidence, "summary")["databases"]["GO"]
    invalid = {
        "schema_version": 1,
        "method": "ora",
        "profile": "summary",
        "databases": {
            "GO": {
                "core_themes": [{"text": "Unsupported claim", "evidence_ids": ["GO:R999"]}],
                "biological_meaning": [], "key_evidence": [], "limitations": [], "computational_checks": [],
            }
        },
    }
    with pytest.raises(ValueError, match="unknown evidence_id"):
        validate_interpretation(invalid, evidence, "summary")


def test_interpreter_normalizes_tf_gsea_activity_wording_without_retry():
    results = {"ChEA3": pd.DataFrame({
        "Term_ID": ["FOXO6_SET"],
        "Term_Name": ["FOXO6 [Coexpression]"],
        "pval": [0.001], "padj": [0.01], "NES": [1.8],
    })}
    invalid = {
        "schema_version": 1, "method": "gsea", "profile": "summary",
        "overall_synthesis": [{
            "text": "FOXO6 activation is represented.",
            "evidence_ids": ["GSEA_ChEA3:R001"],
            "support_class": "single_signal",
        }],
        "databases": {"ChEA3": {
            "core_themes": [], "biological_meaning": [], "key_evidence": [],
            "limitations": [], "computational_checks": [],
        }},
    }
    class FakeBackend:
        calls = 0
        def _call_api(self, prompt):
            FakeBackend.calls += 1
            return json.dumps(invalid)

    previous = AIInterpreter.BACKENDS.get("fake")
    AIInterpreter.BACKENDS["fake"] = FakeBackend
    try:
        output = AIInterpreter("fake").interpret_structured_results(results, "gsea", "summary")
    finally:
        if previous is None:
            AIInterpreter.BACKENDS.pop("fake", None)
        else:
            AIInterpreter.BACKENDS["fake"] = previous
    assert FakeBackend.calls == 1
    assert output["overall_synthesis"][0]["text"] == "FOXO6 transcriptional regulation is represented."
    assert output["overall_synthesis"][0]["auto_normalized_claims"] == ["gsea_direction_wording"]


def test_interpreter_retries_unsupported_support_class_once():
    results = {"GO": pd.DataFrame({
        "Term_ID": ["A", "B"],
        "Term_Name": ["Mitotic spindle", "Chromosome segregation"],
        "P_Value": [0.001, 0.002],
        "Adjusted_P_Value": [0.01, 0.02],
    })}
    invalid = {
        "schema_version": 1, "method": "ora", "profile": "summary",
        "overall_synthesis": [{
            "text": "Mitotic chromosome terms are related.",
            "evidence_ids": ["GO:R001", "GO:R002"],
            "support_class": "shared_core",
        }],
        "databases": {"GO": {
            "core_themes": [], "biological_meaning": [], "key_evidence": [],
            "limitations": [], "computational_checks": [],
        }},
    }
    valid = {
        "schema_version": 1, "method": "ora", "profile": "summary",
        "overall_synthesis": [{
            "text": "Mitotic chromosome terms are related.",
            "evidence_ids": ["GO:R001", "GO:R002"],
            "support_class": "convergent",
        }],
        "databases": {"GO": {
            "core_themes": [], "biological_meaning": [], "key_evidence": [],
            "limitations": [], "computational_checks": [],
        }},
    }

    class FakeBackend:
        calls = 0
        def _call_api(self, prompt):
            FakeBackend.calls += 1
            return json.dumps(invalid if FakeBackend.calls == 1 else valid)

    previous = AIInterpreter.BACKENDS.get("fake")
    AIInterpreter.BACKENDS["fake"] = FakeBackend
    try:
        output = AIInterpreter("fake").interpret_structured_results(results, "ora", "summary")
    finally:
        if previous is None:
            AIInterpreter.BACKENDS.pop("fake", None)
        else:
            AIInterpreter.BACKENDS["fake"] = previous
    assert FakeBackend.calls == 2
    assert output["overall_synthesis"][0]["support_class"] == "convergent"


def test_interpreter_retries_missing_evidence_id_once():
    results = {"GO": pd.DataFrame({
        "Term_ID": ["GO:1"], "Term_Name": ["Signal"],
        "P_Value": [0.001], "Adjusted_P_Value": [0.01],
    })}
    invalid = {
        "schema_version": 1, "method": "ora", "profile": "summary",
        "overall_synthesis": [],
        "databases": {"GO": {
            "core_themes": [], "biological_meaning": [], "key_evidence": [],
            "limitations": [{"text": "A limitation is present.", "evidence_ids": []}],
            "computational_checks": [],
        }},
    }
    valid = {
        "schema_version": 1, "method": "ora", "profile": "summary",
        "overall_synthesis": [],
        "databases": {"GO": {
            "core_themes": [], "biological_meaning": [], "key_evidence": [],
            "limitations": [{"text": "A limitation is present.", "evidence_ids": ["GO:R001"]}],
            "computational_checks": [],
        }},
    }

    class FakeBackend:
        calls = 0
        def _call_api(self, prompt):
            FakeBackend.calls += 1
            return json.dumps(invalid if FakeBackend.calls == 1 else valid)

    previous = AIInterpreter.BACKENDS.get("fake")
    AIInterpreter.BACKENDS["fake"] = FakeBackend
    try:
        output = AIInterpreter("fake").interpret_structured_results(results, "ora", "summary")
    finally:
        if previous is None:
            AIInterpreter.BACKENDS.pop("fake", None)
        else:
            AIInterpreter.BACKENDS["fake"] = previous
    assert FakeBackend.calls == 2
    assert output["databases"]["GO"]["limitations"][0]["evidence_ids"] == ["GO:R001"]


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
    assert "<h2>AI Interpretation</h2>" in html
    assert 'class="ai-section-meta">Model: unknown | Profile: summary; requires expert review</p>' in html
    assert "AI Interpretation <span" not in html
    assert "Research summary" in html
    assert "Evidence pattern" in html
    assert "Biological meaning" in html
    assert "Support: single signal" in html
    assert "Confidence: exploratory" in html
    assert "&lt;script&gt;unsafe()&lt;/script&gt;" in html
    assert "<script>unsafe()</script>" not in html


def test_cli_exposes_ai_mode():
    args = create_parser().parse_args([
        "analyze", "-i", "genes.txt", "--ai", "mock", "--ai-mode", "reviewer", "--ai-top-n", "7"
    ])
    assert args.ai_mode == "reviewer"
    assert args.ai_top_n == 7
