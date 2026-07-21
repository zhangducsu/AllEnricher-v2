"""Real-world SCI E2E framework contract tests."""

from __future__ import annotations

import importlib.util
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
import yaml


ROOT = Path(__file__).parents[1]
SCRIPT = ROOT / "test_e2e_2026" / "18_real_world_sci" / "run_real_world_sci.py"
SPEC = importlib.util.spec_from_file_location("run_real_world_sci", SCRIPT)
assert SPEC and SPEC.loader
RUNNER = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(RUNNER)
PREPARE_SPEC = importlib.util.spec_from_file_location(
    "prepare_database_snapshot",
    ROOT / "test_e2e_2026" / "18_real_world_sci" / "prepare_database_snapshot.py",
)
assert PREPARE_SPEC and PREPARE_SPEC.loader
PREPARE = importlib.util.module_from_spec(PREPARE_SPEC)
PREPARE_SPEC.loader.exec_module(PREPARE)
ATLAS_SPEC = importlib.util.spec_from_file_location(
    "prepare_expression_atlas",
    ROOT / "test_e2e_2026" / "18_real_world_sci" / "prepare_expression_atlas.py",
)
assert ATLAS_SPEC and ATLAS_SPEC.loader
ATLAS = importlib.util.module_from_spec(ATLAS_SPEC)
ATLAS_SPEC.loader.exec_module(ATLAS)


def test_preregistered_matrix_has_exactly_108_cells():
    matrix = yaml.safe_load(RUNNER.DEFAULT_MATRIX.read_text(encoding="utf-8"))
    cells = sum(len(dataset["databases"]) * len(RUNNER.METHODS) for dataset in matrix["datasets"].values())

    assert cells == 108
    assert set(RUNNER.METHODS) == {"hypergeometric", "gsea", "ssgsea", "gsva"}


def test_offline_command_environment_blocks_network_proxies():
    env = RUNNER.command_environment(offline=True)

    assert env["ALLENRICHER_OFFLINE"] == "1"
    assert env["NO_PROXY"] == ""
    assert {env[name] for name in ("HTTP_PROXY", "HTTPS_PROXY", "ALL_PROXY")} == {
        "http://127.0.0.1:9"
    }


def test_paper_results_separate_main_matrix_from_python_appendix(tmp_path, monkeypatch):
    input_root = tmp_path / "inputs"
    dataset_root = input_root / "sample"
    dataset_root.mkdir(parents=True)
    (dataset_root / "provenance.json").write_text(
        """{
          "checks": {"query_genes": 10, "ranked_genes": 1000,
                     "reference_samples": 3, "test_samples": 3},
          "source_url": "https://example.test/study",
          "license": "CC BY 4.0"
        }""",
        encoding="utf-8",
    )
    monkeypatch.setattr(RUNNER, "INPUT_ROOT", input_root)
    matrix = {
        "datasets": {
            "sample": {
                "accession": "E-TEST-1",
                "latin_name": "Species example",
                "contrast": "test vs control",
            }
        }
    }
    main = [{"case_id": "main", "status": "PASS", "errors": []}]
    appendix = [{"case_id": "appendix", "status": "PASS", "errors": []}]

    RUNNER.write_paper_materials(tmp_path / "run", matrix, main, appendix)

    results = (tmp_path / "run" / "paper_materials" / "RESULTS.md").read_text(encoding="utf-8")
    assert "Main analysis matrix: 1/1 passed" in results
    assert "Python visualization appendix: 1/1 passed" in results


def test_gsva_command_uses_public_counts_groups_and_poisson_config(tmp_path):
    dataset = {"species": "hsa"}
    inputs = {
        "query": tmp_path / "query.txt",
        "ranked": tmp_path / "ranked.tsv",
        "expression": tmp_path / "expression.tsv",
        "groups": tmp_path / "groups.txt",
    }
    inputs["groups"].write_text("Control:S1;Treatment:S2\n", encoding="utf-8")
    config = tmp_path / "config.yaml"
    config.write_text("gsva_kcdf: Poisson\n", encoding="utf-8")

    command = RUNNER.build_command(
        "gsva", dataset, "GO", inputs, tmp_path / "output", config, use_r=True
    )

    assert command[command.index("-e") + 1] == str(inputs["expression"])
    assert command[command.index("--groups") + 1] == "Control:S1;Treatment:S2"
    assert command[command.index("--config") + 1] == str(config)
    assert "--use-r-plots" in command


def test_python_appendix_requests_only_retained_python_gsea_plots(tmp_path):
    inputs = {
        "query": tmp_path / "query.txt",
        "ranked": tmp_path / "ranked.tsv",
        "expression": tmp_path / "expression.tsv",
        "groups": tmp_path / "groups.txt",
    }
    command = RUNNER.build_command(
        "gsea", {"species": "hsa"}, "GO", inputs, tmp_path / "output",
        tmp_path / "config.yaml", use_r=False,
    )
    plots = command[command.index("-pt") + 1].split(",")

    assert set(plots) == {"lollipop", "ridgeplot", "enrichment", "enrichment2"}
    assert "barplot" not in plots
    assert "emapplot" not in plots
    assert "--use-r-plots" not in command


def test_ora_oracle_corrects_positive_overlap_hypotheses_like_v1(tmp_path):
    query = tmp_path / "query.txt"
    background = tmp_path / "background.txt"
    query.write_text("G1\nG2\n", encoding="utf-8")
    background.write_text("G1\nG2\nG3\nG4\nG5\nG6\n", encoding="utf-8")
    terms = {
        "HIT": {"name": "Hit term", "genes": {"G1", "G2", "G3"}},
        "MISS": {"name": "Miss term", "genes": {"G4", "G5", "G6"}},
    }
    actual = pd.DataFrame(
        {
            "Term_ID": ["HIT"],
            "Term_Name": ["Hit term"],
            "P_Value": [1 / 5],
            "Adjusted_P_Value": [1 / 5],
            "Gene_Count": [2],
            "Genes": ["G1;G2"],
        }
    )
    oracle = tmp_path / "oracle"
    oracle.mkdir()

    assert RUNNER.oracle_ora(actual, terms, query, background, oracle, "GO") == []
    saved = pd.read_csv(oracle / "independent_hypergeometric.tsv", sep="\t")
    assert saved["Term_ID"].tolist() == ["HIT"]
    assert np.isclose(saved.loc[0, "Adjusted_P_Value"], 1 / 5)


def test_single_and_multi_enrichment_plot_tokens_are_not_conflated():
    assert not RUNNER.has_plot_token(["go_enrichment2_up"], "enrichment")
    assert RUNNER.has_plot_token(["go_enrichment_up"], "enrichment")
    assert RUNNER.has_plot_token(["go_enrichment2_down"], "enrichment2")


def test_case_validator_rejects_report_that_omits_generated_plot(tmp_path, monkeypatch):
    output = tmp_path / "output"
    output.mkdir()
    pd.DataFrame({
        "Term_ID": ["GO:1"],
        "Term_Name": ["Readable pathway"],
        "P_Value": [0.01],
        "Adjusted_P_Value": [0.02],
        "Gene_Count": [2],
        "Genes": ["G1;G2"],
    }).to_csv(output / "GO_enrichment.tsv", sep="\t", index=False)
    (output / "analysis_metadata.json").write_text(
        '{"analysis_method":"hypergeometric","database_versions":{"GO":"v1"}}',
        encoding="utf-8",
    )
    (output / "GO_barplot.png").write_bytes(b"plot")
    (output / "report.html").write_text(
        '<html><body>ORA Enrichment Analysis Report Term ID Term Name '
        'Readable pathway</body></html>',
        encoding="utf-8",
    )
    query = tmp_path / "query.txt"
    background = tmp_path / "background.txt"
    query.write_text("G1\nG2\n", encoding="utf-8")
    background.write_text("G1\nG2\nG3\n", encoding="utf-8")
    monkeypatch.setattr(RUNNER, "oracle_ora", lambda *_args, **_kwargs: [])

    errors = RUNNER.validate_case_output(
        output,
        "GO",
        "hypergeometric",
        {"GO:1": {"name": "Readable pathway", "genes": {"G1", "G2"}}},
        {"query": query, "background": background},
        tmp_path / "oracle",
        requested_plot_tokens=("barplot",),
    )

    assert "report embeds 0/1 generated plot types" in errors


def test_contact_sheet_copy_names_stay_below_windows_max_path(tmp_path, monkeypatch):
    run_dir = tmp_path / "run"
    output = run_dir / "cases" / "drosophila_fas_p218__PUBLIC_GO_CUSTOM__hypergeometric" / "output"
    output.mkdir(parents=True)
    (output / "PUBLIC_GO_CUSTOM_barplot.png").write_bytes(b"plot")
    monkeypatch.setattr(RUNNER.subprocess, "run", lambda *_args, **_kwargs: type("Result", (), {"returncode": 0, "stderr": ""})())

    errors = RUNNER.create_contact_sheets(
        run_dir,
        {"datasets": {"drosophila_fas_p218": {}}},
    )

    copied = next((run_dir / "visual_review" / "drosophila_fas_p218" / "hypergeometric").glob("*.png"))
    assert copied.name == "PUBLIC_GO_CUSTOM_01_PUBLIC_GO_CUSTOM_barplot.png"
    old_name = "drosophila_fas_p218__PUBLIC_GO_CUSTOM__hypergeometric_01_PUBLIC_GO_CUSTOM_barplot.png"
    assert len(copied.name) < len(old_name)
    assert errors  # Other methods intentionally have no fixture images.


def test_public_go_custom_equivalence_ignores_only_ora_source_metadata(tmp_path):
    dataset = "example"
    matrix = {"datasets": {dataset: {"databases": ["GO", "PUBLIC_GO_CUSTOM"]}}}
    for method in RUNNER.METHODS:
        if method == "hypergeometric":
            standard = pd.DataFrame({
                "Term_ID": ["GO:1"], "Term_Name": ["DNA Repair"], "Database": ["GO"],
                "P_Value": [0.01], "Term_URL": ["https://example.org/GO:1"],
            })
            custom = standard.assign(Database="PUBLIC_GO_CUSTOM", Term_URL="")
        elif method == "gsea":
            standard = custom = pd.DataFrame({
                "Term_ID": ["GO:1"], "Term_Name": ["DNA Repair"],
                "pathway": ["GO:1"], "NES": [1.5],
            })
        else:
            standard = custom = pd.DataFrame({
                "Term_ID": ["GO:1"], "Term_Name": ["DNA Repair"], "S1": [0.5],
            })
        for database, frame in (("GO", standard), ("PUBLIC_GO_CUSTOM", custom)):
            output = tmp_path / "cases" / f"{dataset}__{database}__{method}" / "output"
            output.mkdir(parents=True)
            frame.to_csv(output / f"{database}_enrichment.tsv", sep="\t", index=False)

    assert RUNNER.compare_public_go_custom(tmp_path, matrix) == []


def test_offline_database_prepare_only_validates_frozen_snapshot(tmp_path, monkeypatch):
    database_root = tmp_path / "database"
    (database_root / "organism").mkdir(parents=True)
    (database_root / "SOURCE_MANIFEST.json").write_text("{}\n", encoding="utf-8")
    (database_root / "DATABASE_AUDIT.json").write_text("{}\n", encoding="utf-8")
    matrix_path = tmp_path / "matrix.yaml"
    matrix_path.write_text("datasets: {}\n", encoding="utf-8")
    monkeypatch.setattr(PREPARE, "validate_snapshot", lambda *_args: {"validated": True})
    monkeypatch.setattr(PREPARE, "source_manifest", lambda *_args: {"refreshed": True})
    monkeypatch.setattr(PREPARE, "fetch_tf_sources", lambda *_args: (_ for _ in ()).throw(AssertionError("network")))

    assert PREPARE.prepare(matrix_path, database_root, tmp_path / "basic", offline=True) == {
        "validated": True
    }
    assert yaml.safe_load((database_root / "DATABASE_AUDIT.json").read_text(encoding="utf-8")) == {
        "validated": True
    }


def test_runner_rejects_stale_database_audit(tmp_path):
    audit_path = tmp_path / "DATABASE_AUDIT.json"
    audit_path.write_text('{"old_dataset": {}}\n', encoding="utf-8")

    try:
        RUNNER.validate_database_audit({"datasets": {"current_dataset": {}}}, audit_path)
    except ValueError as exc:
        assert "missing: current_dataset" in str(exc)
        assert "stale: old_dataset" in str(exc)
    else:
        raise AssertionError("stale audit was accepted")


def test_database_prepare_resume_requires_a_loadable_valid_species_snapshot(tmp_path, monkeypatch):
    spec = {"species": "hsa", "databases": ["GO"]}
    monkeypatch.setattr(PREPARE, "validate_snapshot", lambda *_args: {"ready": True})
    assert PREPARE.species_snapshot_ready(tmp_path, spec)

    def fail_validation(*_args):
        raise FileNotFoundError("missing database")

    monkeypatch.setattr(PREPARE, "validate_snapshot", fail_validation)
    assert not PREPARE.species_snapshot_ready(tmp_path, spec)


def test_ssgsea_oracle_uses_official_alpha_quarter(tmp_path, monkeypatch):
    expression = tmp_path / "expression.tsv"
    expression.write_text("gene\tS1\nG1\t1\n", encoding="utf-8")
    expected = pd.DataFrame({"S1": [0.5]}, index=["P1"])
    expected.index.name = "Term_ID"
    captured = {}

    def fake_run_gsva(*_args, **kwargs):
        captured.update(kwargs)
        return expected.copy()

    monkeypatch.setattr(RUNNER, "run_gsva", fake_run_gsva)
    actual = expected.reset_index().assign(Term_Name="Pathway A")
    (tmp_path / "oracle").mkdir()

    assert not RUNNER.oracle_activity(
        actual,
        {"P1": {"name": "Pathway A", "genes": ["G1"]}},
        expression,
        "ssgsea",
        tmp_path / "oracle",
        "GO",
    )
    assert captured["tau"] == 0.25


def test_expression_atlas_retrieval_timestamp_is_stable_for_cached_input(tmp_path):
    original = "2026-07-15T00:00:00+00:00"
    (tmp_path / "provenance.json").write_text(
        '{"accession":"E-TEST-1","retrieved_at":"' + original + '"}\n', encoding="utf-8"
    )

    assert ATLAS.stable_retrieved_at(tmp_path, "E-TEST-1") == original
    assert ATLAS.stable_retrieved_at(tmp_path, "E-OTHER") != original


def test_expression_atlas_sample_selection_applies_all_contrast_conditions():
    design = pd.DataFrame(
        {
            "Run": ["HR0_1", "HR0_2", "HR0_3", "HR1_1", "HR1_2", "HR1_3", "LR0", "LR1", "OTHER"],
            "Analysed": ["Yes"] * 9,
            "Factor Value[time]": [
                "0 week", "0 week", "0 week", "1 week", "1 week", "1 week", "0 week", "1 week", "2 week"
            ],
            "Factor Value[phenotype]": [
                "high risk", "high risk", "high risk", "high risk", "high risk", "high risk",
                "low risk", "low risk", "high risk"
            ],
        }
    )

    reference, test = ATLAS._pick_samples(
        design,
        "Factor Value[time]",
        "0 week",
        "1 week",
        {"Factor Value[phenotype]": "high risk"},
    )

    assert reference == ["HR0_1", "HR0_2", "HR0_3"]
    assert test == ["HR1_1", "HR1_2", "HR1_3"]


def test_expression_atlas_aggregates_technical_runs_by_biological_unit():
    expression = pd.DataFrame(
        {
            "LR1": [1, 2], "LR2": [3, 4], "LR3": [5, 6], "LR4": [7, 8],
            "HR1": [2, 1], "HR2": [4, 3], "HR3": [6, 5], "HR4": [8, 7],
        },
        index=["G1", "G2"],
    )
    design = pd.DataFrame(
        {
            "Run": ["LR1", "LR2", "LR3", "LR4", "HR1", "HR2", "HR3", "HR4"],
            "Sample Characteristic[individual]": ["1", "1", "2", "3", "6", "6", "7", "8"],
        }
    )

    aggregated, reference, test = ATLAS._aggregate_technical_replicates(
        expression,
        design,
        ["LR1", "LR2", "LR3", "LR4"],
        ["HR1", "HR2", "HR3", "HR4"],
        "Sample Characteristic[individual]",
    )

    assert reference == ["Control_1", "Control_2", "Control_3"]
    assert test == ["Treatment_6", "Treatment_7", "Treatment_8"]
    assert aggregated.loc["G1", "Control_1"] == 4
    assert aggregated.loc["G2", "Treatment_6"] == 4


def test_expression_atlas_query_ids_use_one_consistent_analytics_mapping():
    analytics = pd.DataFrame({
        "Gene ID": ["ID1", "ID2"],
        "Gene Name": ["GENE1", pd.NA],
    })
    analytics["Analysis Gene"] = ATLAS._analysis_gene_ids(analytics)
    query = pd.DataFrame({
        "Gene ID": ["ID1", "ID2"],
        # Atlas filtered downloads may put Gene ID here when no symbol exists.
        "Gene Name": ["GENE1", "ID2"],
    })

    mapped = ATLAS._map_query_to_analytics_ids(query, analytics)

    assert mapped.iloc[0] == "GENE1"
    assert mapped.iloc[1] == "ID2"


def test_real_world_ora_command_uses_expression_background(tmp_path):
    inputs = {
        "query": tmp_path / "query.txt",
        "background": tmp_path / "background.txt",
        "ranked": tmp_path / "ranked.tsv",
        "expression": tmp_path / "expression.tsv",
        "groups": tmp_path / "groups.txt",
    }
    command = RUNNER.build_command(
        "hypergeometric", {"species": "hsa"}, "GO", inputs,
        tmp_path / "output", tmp_path / "config.yaml",
    )

    assert command[command.index("-b") + 1] == str(inputs["background"])


def test_frozen_real_world_inputs_share_one_valid_gene_universe():
    matrix = yaml.safe_load(RUNNER.DEFAULT_MATRIX.read_text(encoding="utf-8"))
    for case_id in matrix["datasets"]:
        converted = RUNNER.INPUT_ROOT / case_id / "converted"
        query = {
            line.strip()
            for line in (converted / "query_genes.txt").read_text(encoding="utf-8").splitlines()
            if line.strip()
        }
        background = {
            line.strip()
            for line in (converted / "background_genes.txt").read_text(encoding="utf-8").splitlines()
            if line.strip()
        }
        ranked = pd.read_csv(converted / "ranked_genes.tsv", sep="\t")
        expression = pd.read_csv(converted / "expression_counts.tsv", sep="\t", index_col=0)
        ranked_genes = set(ranked["gene"].astype(str))
        weights = pd.to_numeric(ranked["weight"], errors="raise")

        assert ranked["gene"].astype(str).is_unique, case_id
        assert expression.index.astype(str).is_unique, case_id
        assert np.isfinite(weights).all(), case_id
        assert weights.is_monotonic_decreasing, case_id
        assert (weights > 0).any() and (weights < 0).any(), case_id
        assert query <= ranked_genes, case_id
        assert query <= background, case_id
        assert ranked_genes == background == set(expression.index.astype(str)), case_id


def test_real_world_matrix_uses_independent_cattle_and_yeast_fixtures():
    matrix = yaml.safe_load(RUNNER.DEFAULT_MATRIX.read_text(encoding="utf-8"))
    cattle = matrix["datasets"]["cattle_metabolic_risk_week1"]
    yeast = matrix["datasets"]["yeast_tsa1_deletion_stress"]

    assert cattle["species"] == "bta"
    assert cattle["accession"] == "E-MTAB-5838"
    assert cattle["factor_column"] == "Factor Value[phenotype]"
    assert cattle["sample_filters"] == {"Factor Value[time]": "1 week"}
    assert cattle["biological_unit_column"] == "Sample Characteristic[individual]"
    assert yeast["accession"] == "E-MTAB-9355"
    assert yeast["sample_filters"]["Factor Value[compound]"] == (
        "azetidine-2-carboxylic acid 5 millimolar"
    )
    assert "cattle_metabolic_high_risk_week1" not in matrix["datasets"]
    assert "pig_cartilage_bear_day7_female" not in matrix["datasets"]
    assert "yeast_sec66_deletion" not in matrix["datasets"]


def test_expression_atlas_download_is_atomic_and_bypasses_environment_proxy(tmp_path, monkeypatch):
    class Response:
        headers = {"Content-Length": "6"}

        def __enter__(self):
            return self

        def __exit__(self, *_args):
            return False

        def raise_for_status(self):
            return None

        def iter_content(self, _chunk_size):
            yield b"abc"
            yield b"def"

    class Session:
        def __init__(self):
            self.trust_env = True

        def get(self, *_args, **_kwargs):
            assert self.trust_env is False
            return Response()

        def close(self):
            return None

    monkeypatch.setattr(ATLAS.requests, "Session", Session)
    target = tmp_path / "source" / "analytics.tsv"

    ATLAS.download("https://example.test/analytics", target, offline=False)

    assert target.read_bytes() == b"abcdef"
    assert not target.with_suffix(".tsv.part").exists()


def test_expression_atlas_download_removes_partial_file_after_retries(tmp_path, monkeypatch):
    calls = 0

    class Response:
        headers = {}

        def __enter__(self):
            return self

        def __exit__(self, *_args):
            return False

        def raise_for_status(self):
            return None

        def iter_content(self, _chunk_size):
            yield b"partial"
            raise OSError("connection interrupted")

    class Session:
        trust_env = True

        def get(self, *_args, **_kwargs):
            nonlocal calls
            calls += 1
            return Response()

        def close(self):
            return None

    monkeypatch.setattr(ATLAS.requests, "Session", Session)
    monkeypatch.setattr(ATLAS.time, "sleep", lambda _seconds: None)
    target = tmp_path / "source" / "analytics.tsv"

    with pytest.raises(RuntimeError, match="download failed after 4 attempts"):
        ATLAS.download("https://example.test/analytics", target, offline=False)

    assert calls == 4
    assert not target.exists()
    assert not target.with_suffix(".tsv.part").exists()


def test_animaltfdb_gene_info_subset_keeps_only_requested_taxa(tmp_path):
    import gzip

    source = tmp_path / "gene_info.gz"
    target = tmp_path / "selected.gz"
    with gzip.open(source, "wt", encoding="utf-8") as handle:
        handle.write("#tax_id\tGeneID\tSymbol\tLocusTag\tSynonyms\tdbXrefs\n")
        handle.write("9606\t1\tH1\t-\t-\tEnsembl:ENSG1\n")
        handle.write("7227\t2\tF1\t-\t-\tFlyBase:FBgn1\n")
        handle.write("10090\t3\tM1\t-\t-\tEnsembl:ENSMUSG1\n")

    PREPARE.prepare_animaltfdb_gene_info(source, target, {9606, 7227})
    with gzip.open(target, "rt", encoding="utf-8") as handle:
        text = handle.read()

    assert "9606\t" in text
    assert "7227\t" in text
    assert "10090\t" not in text

    PREPARE.prepare_animaltfdb_gene_info(source, target, {9606, 7227, 10090})
    with gzip.open(target, "rt", encoding="utf-8") as handle:
        rebuilt = handle.read()

    assert "10090\t" in rebuilt
