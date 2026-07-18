from pathlib import Path

import pandas as pd

from allenricher.report.generator import ReportGenerator
from allenricher.report.methods_reference import build_methods_reference
from allenricher.cli import (
    _recorded_analysis_databases,
    _recorded_analysis_parameters,
    create_parser,
)
from allenricher.core.config import Config


def test_ora_reference_uses_only_recorded_values():
    content = build_methods_reference({
        "allenricher_version": "2.3.1",
        "analysis_method": "hypergeometric",
        "species": "hsa",
        "species_name": "Homo sapiens",
        "species_taxonomy_id": 9606,
        "databases": ["GO", "CUSTOM"],
        "database_versions": {"GO": "v20260701"},
        "parameters": {
            "background_mode": "custom",
            "correction": "BY",
            "pvalue_cutoff": 0.01,
            "qvalue_cutoff": 0.02,
            "gene_set_size_by_database": {
                "GO": {"min": 3, "max": "Inf"},
                "CUSTOM": {"min": 5, "max": 100},
            },
        },
    })

    paragraph = content["paragraphs"][0]
    assert "AllEnricher version 2.3.1" in paragraph
    assert "Homo sapiens (TaxID: 9606; hsa)" in paragraph
    assert "GO (version: v20260701)" in paragraph
    assert "CUSTOM (version: To be added)" in paragraph
    assert "custom" in paragraph and "BY" in paragraph
    assert "GO: 3-unbounded" in paragraph
    assert all(item["source"] not in {"KEGG", "Reactome"} for item in content["references"])
    assert content["references"][-1]["citation"] == "CUSTOM: To be added"


def test_english_gsea_reference_does_not_invent_unrecorded_parameters():
    content = build_methods_reference({
        "methods_language": "en",
        "allenricher_version": "2.3.1",
        "analysis_method": "gsea",
        "species": "mmu",
        "databases": ["TRRUST"],
        "database_versions": {"TRRUST": "v1"},
        "parameters": {
            "pvalue_cutoff": 0.05,
            "qvalue_cutoff": 0.05,
            "gene_set_size_by_database": {"TRRUST": {"min": 15, "max": 5000}},
        },
    })

    paragraph = content["paragraphs"][0]
    assert "gene set enrichment analysis (GSEA)" in paragraph
    assert "TRRUST: 15-5000" in paragraph
    assert "permutation" not in paragraph.lower()
    assert "fgsea" not in paragraph.lower()
    assert content["references"][1]["source"] == "TRRUST"


def test_disgenet_uses_the_frozen_v1_snapshot_version():
    content = build_methods_reference({
        "allenricher_version": "2.3.1",
        "analysis_method": "hypergeometric",
        "species": "hsa",
        "databases": ["DisGeNET"],
        "parameters": {
            "background_mode": "annotated",
            "correction": "BH",
            "pvalue_cutoff": 0.05,
            "qvalue_cutoff": 0.05,
            "gene_set_size_by_database": {
                "DisGeNET": {"min": 3, "max": "Inf"},
            },
        },
    })

    assert "DisGeNET (v20190612)" in content["paragraphs"][0]
    assert "DisGeNET (version)" not in content["paragraphs"][0]


def test_cli_records_method_specific_effective_size_limits():
    config = Config(method="gsea", databases=["GO", "TRRUST"])
    parameters = _recorded_analysis_parameters(config, "annotated")
    assert parameters["gene_set_size_by_database"] == {
        "GO": {"min": 15, "max": 500},
        "TRRUST": {"min": 15, "max": 5000},
    }
    args = create_parser().parse_args([
        "analyze", "-i", "genes.txt", "--methods-language", "en"
    ])
    assert args.methods_language == "en"


def test_cli_records_only_databases_that_contributed_results():
    assert _recorded_analysis_databases(
        ["GO"], "trrust", False, True
    ) == ["GO", "TRRUST"]
    assert _recorded_analysis_databases(
        ["GO"], "both", True, True
    ) == ["TRRUST", "ChEA3"]
    assert _recorded_analysis_databases(
        ["GO"], "trrust", False, False
    ) == ["GO"]


def test_html_report_contains_shared_methods_reference(tmp_path):
    output = tmp_path / "report.html"
    results = {
        "KEGG": pd.DataFrame({
            "Term_ID": ["hsa04110"],
            "Term_Name": ["Cell cycle"],
            "Gene_Count": [4],
            "P_Value": [0.001],
            "Adjusted_P_Value": [0.01],
        })
    }
    metadata = {
        "allenricher_version": "2.3.1",
        "analysis_method": "hypergeometric",
        "species": "hsa",
        "databases": ["KEGG"],
        "database_versions": {"KEGG": "v20260701"},
        "parameters": {
            "background_mode": "annotated",
            "correction": "BH",
            "pvalue_cutoff": 0.05,
            "qvalue_cutoff": 0.05,
            "gene_set_size_by_database": {"KEGG": {"min": 3, "max": "Inf"}},
        },
    }

    ReportGenerator(str(tmp_path)).generate(
        results,
        str(output),
        analysis_method="hypergeometric",
        metadata=metadata,
    )
    html = Path(output).read_text(encoding="utf-8")
    assert 'id="methods-reference"' in html
    assert "Materials and Methods Writing Reference" in html
    assert "doi:10.1093/nar/28.1.27" in html
    assert "WikiPathways" not in html
