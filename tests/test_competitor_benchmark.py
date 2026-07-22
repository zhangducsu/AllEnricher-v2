"""Deterministic checks for the preregistered competitor benchmark."""

from __future__ import annotations

import importlib.util
import re
from pathlib import Path

import pandas as pd
import yaml


ROOT = Path(__file__).resolve().parents[1]
SUITE = ROOT / "test_e2e_2026" / "19_competitor_benchmark"
SPEC = importlib.util.spec_from_file_location("competitor_benchmark", SUITE / "run_benchmark.py")
benchmark = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(benchmark)
FIGURE_SPEC = importlib.util.spec_from_file_location(
    "competitor_publication", SUITE / "make_publication_figures.py"
)
publication = importlib.util.module_from_spec(FIGURE_SPEC)
assert FIGURE_SPEC.loader is not None
FIGURE_SPEC.loader.exec_module(publication)


def test_matrix_and_gmt_contract(tmp_path: Path) -> None:
    config = yaml.safe_load((SUITE / "benchmark_matrix.yaml").read_text(encoding="utf-8"))
    cases = benchmark.expected_result_sets(config)
    assert len(cases) == 52
    assert sum(case["tool"] == "g:Profiler" for case in cases) == 4

    source = tmp_path / "source.gmt"
    source.write_text("T1\tTerm one\tA\tB\tC\nT2\tTerm two\tB\tD\n", encoding="utf-8")
    target = tmp_path / "normalized.gmt"
    terms = benchmark.normalize_gmt(source, target, {"A", "B", "D"}, 2, None)
    assert list(terms) == ["T1", "T2"]
    assert target.read_text(encoding="utf-8") == "T1\tTerm one\tA\tB\nT2\tTerm two\tB\tD\n"


def test_case_level_metrics_do_not_treat_terms_as_replicates() -> None:
    rows = []
    for tool, q_values in (("AllEnricher", [0.01, 0.2]), ("clusterProfiler", [0.01, 0.2])):
        for term_id, q_value in zip(("T1", "T2"), q_values):
            rows.append({
                "tool": tool, "tool_version": "test", "dataset": "D", "species": "hsa",
                "database": "GO", "method": "ORA", "term_id": term_id,
                "term_name": term_id, "p_value": q_value / 2,
                "adjusted_p_value": q_value, "es": None, "nes": None,
                "leading_edge": "", "overlap_count": 1, "term_size": 3, "status": "PASS",
            })
    frame = pd.DataFrame(rows, columns=benchmark.NORMALIZED_COLUMNS)
    metric, _ = benchmark.compare_pair(
        frame[frame.tool == "AllEnricher"], frame[frame.tool == "clusterProfiler"],
        {"tool": "clusterProfiler", "dataset": "D", "database": "GO", "method": "ORA"},
    )
    assert metric["term_jaccard"] == 1
    assert metric["spearman"] > 0.999999
    assert metric["max_abs_p_diff"] == 0
    assert metric["max_abs_q_diff"] == 0

    expected = [{"tool": "clusterProfiler", "dataset": "D", "database": "GO", "method": "ORA"}]
    acceptance = {
        "ora_max_abs_p_diff": 1e-10,
        "ora_max_abs_q_diff": 1e-10,
        "clusterprofiler_gsea_max_abs_nes_diff": 1e-8,
        "clusterprofiler_gsea_max_abs_q_diff": 1e-8,
    }
    metrics, _ = benchmark.build_metrics(frame, expected, acceptance)
    assert metrics.iloc[0]["status"] == "PASS"

    frame.loc[(frame.tool == "clusterProfiler") & (frame.term_id == "T2"), "adjusted_p_value"] = 0.3
    metrics, _ = benchmark.build_metrics(frame, expected, acceptance)
    assert metrics.iloc[0]["status"] == "FAIL_NUMERIC"
    assert "BH-FDR" in metrics.iloc[0]["reason"]


def test_gsea_metrics_report_nominal_p_differences_and_valid_pair_counts() -> None:
    rows = []
    values = {
        "AllEnricher": ([0.01, 0.20, 0.50], [0.03, 0.30, 0.70], [1.5, -1.2, 0.8]),
        "clusterProfiler": ([0.02, 0.19, 0.70], [0.04, 0.31, 0.90], [1.4, -1.1, None]),
    }
    for tool, (p_values, q_values, nes_values) in values.items():
        for term_id, p_value, q_value, nes in zip(("T1", "T2", "T3"), p_values, q_values, nes_values):
            rows.append({
                "tool": tool, "tool_version": "test", "dataset": "D", "species": "hsa",
                "database": "GO", "method": "GSEA", "term_id": term_id,
                "term_name": term_id, "p_value": p_value,
                "adjusted_p_value": q_value, "es": nes, "nes": nes,
                "leading_edge": "A;B", "overlap_count": None, "term_size": 15, "status": "PASS",
            })
    frame = pd.DataFrame(rows, columns=benchmark.NORMALIZED_COLUMNS)
    metric, detail = benchmark.compare_pair(
        frame[frame.tool == "AllEnricher"], frame[frame.tool == "clusterProfiler"],
        {"tool": "clusterProfiler", "dataset": "D", "database": "GO", "method": "GSEA"},
    )

    assert abs(metric["max_abs_p_diff"] - 0.20) < 1e-12
    assert abs(detail["median_abs_p_diff"] - 0.01) < 1e-12
    assert detail["p_spearman"] > 0.999999
    assert detail["valid_p_pairs"] == 3
    assert detail["valid_q_pairs"] == 3
    assert detail["valid_nes_pairs"] == 2


def test_metric_points_dodge_species_and_databases() -> None:
    frame = pd.DataFrame([
        {"comparator": "clusterProfiler", "dataset": dataset, "database": database, "spearman": 1.0}
        for dataset in ("human_airway_dex", "cattle_metabolic_risk_week1")
        for database in ("GO", "KEGG")
    ])
    figure, axis = publication.plt.subplots()
    try:
        publication.plot_metric_points(axis, frame, ["spearman"], "Agreement")
        points = axis.collections
        x_positions = [round(float(point.get_offsets()[0, 0]), 6) for point in points]
        assert len(x_positions) == 4
        assert len(set(x_positions)) == 4
        assert tuple(points[0].get_facecolors()[0, :3]) != (1.0, 1.0, 1.0)
        assert tuple(points[1].get_facecolors()[0, :3]) != (1.0, 1.0, 1.0)
        assert tuple(points[2].get_facecolors()[0, :3]) == (1.0, 1.0, 1.0)
        assert tuple(points[3].get_facecolors()[0, :3]) == (1.0, 1.0, 1.0)
        assert points[0].get_edgecolors().size == 0
        assert points[1].get_edgecolors().size == 0
        assert tuple(points[2].get_edgecolors()[0, :3]) != (0.0, 0.0, 0.0)
        assert all(point.get_alpha() == 0.78 for point in points)
    finally:
        publication.plt.close(figure)


def test_publication_rejects_stale_ora_fdr_metrics() -> None:
    config = {
        "datasets": {"D": {}}, "databases": ["GO"],
        "acceptance": {"ora_max_abs_p_diff": 1e-10, "ora_max_abs_q_diff": 1e-10},
    }
    metrics = pd.DataFrame([{
        "comparator": "clusterProfiler", "dataset": "D", "database": "GO", "method": "ORA",
        "reference_terms": 2, "comparator_terms": 2, "term_jaccard": 1,
        "max_abs_p_diff": 0, "max_abs_q_diff": 0.5, "status": "PASS",
    }])
    try:
        publication.validate_publication_metrics(metrics, config)
    except ValueError as exc:
        assert "D/GO" in str(exc)
    else:
        raise AssertionError("stale ORA FDR metrics were accepted for publication")

def test_publication_exports_are_portable_and_leave_figure1_to_author() -> None:
    gmt = publication.source_gmt_path(
        r"D:\archive\database_snapshot\organism\v20260715\hsa\hsa.GO.gmt.gz"
    )
    assert gmt == "database_snapshot/organism/v20260715/hsa/hsa.GO.gmt.gz"

    command = publication.portable_command(
        r"D:\R\bin\Rscript.exe D:\repo\AllEnricher-v2\test_e2e_2026\99_runs\20260720_RUN\raw\case\result.tsv"
    )
    assert command == "Rscript {ARCHIVE_ROOT}/raw/case/result.tsv"
    assert not any(token in command for token in ("D:\\", "/mnt/", "file://"))

    manifest = {
        "created_at": "2026-07-20T15:02:31+00:00",
        "tool_versions": {
            "AllEnricher": "2.1.0", "clusterProfiler": "4.20.0",
            "WebGestaltR": "1.0.0", "g:Profiler": "test",
            "getENRICH": "5235d2d9eb3234b41b9c6f507e399f8d0e6b80bc",
            "GSVA": "2.6.2",
        },
    }
    summary = publication.tool_summary(manifest)
    assert summary.shape[0] == 6
    templates = "\n".join(summary["command_template"])
    assert "run_getenrich_case.R" not in templates
    assert "direct_gsva_oracle.R" not in templates
    assert (ROOT / "test_e2e_2026" / "18_real_world_sci" / "run_real_world_sci.py").is_file()
    assert (SUITE / "run_benchmark.py").is_file()
    assert (SUITE / "competitor_methods.R").is_file()
    assert list(summary.columns) == [
        "tool", "version", "role", "execution_environment", "access_date",
        "commit_sha", "command_template", "evidence_archive_path", "source_url",
    ]
    assert not hasattr(publication, "figure1")
    assert not hasattr(publication, "supplement_s3")
    assert "Figure_1_workflow_validation" not in publication.generate.__code__.co_consts
    assert "Figure_S3_activity_method_error" not in publication.generate.__code__.co_consts

def test_generated_publication_tables_preserve_records_without_path_leaks() -> None:
    paper = ROOT.parent / "Paper"
    supplementary = paper / "supplementary"
    table_s1 = pd.read_csv(supplementary / "Table_S1_datasets_inputs_databases.tsv", sep="\t")
    table_s2 = pd.read_csv(supplementary / "Table_S2_capability_evidence.tsv", sep="\t")
    table_s3 = pd.read_csv(supplementary / "Table_S3_versions_commands_access.tsv", sep="\t")
    table_s4 = pd.read_csv(
        supplementary / "Table_S4_case_metrics_failures.tsv", sep="\t", keep_default_na=False
    )
    commands = pd.read_csv(
        supplementary / "source_data" / "Data_S1_full_case_commands.tsv", sep="\t"
    )

    assert len(table_s1) == 8
    assert set(table_s1["annotation_distribution"]) == {
        "CC BY 4.0", "Not redistributed; hash only",
    }
    frozen_inputs = pd.read_csv(
        ROOT / "test_e2e_2026" / "99_runs" / "20260720_2300_COMPETITOR_OFFLINE_REPLAY_FINAL" / "input_statistics.tsv",
        sep="\t", dtype=str, keep_default_na=False,
    )
    published_inputs = pd.read_csv(
        supplementary / "Table_S1_datasets_inputs_databases.tsv",
        sep="\t", dtype=str, keep_default_na=False,
    )
    unchanged_columns = [column for column in frozen_inputs.columns if column != "source_gmt"]
    pd.testing.assert_frame_equal(
        frozen_inputs[unchanged_columns], published_inputs[unchanged_columns], check_dtype=False
    )
    assert len(commands) == 56
    assert commands["case_id"].is_unique
    assert list(table_s3.columns) == [
        "tool", "version", "role", "execution_environment", "access_date",
        "commit_sha", "command_template", "evidence_archive_path", "source_url",
    ]
    assert set(table_s2["value"]) <= {"Yes", "Partial", "No", "N/A"}
    assert len(table_s2) == len(publication.CAPABILITIES) * len(publication.FEATURES)
    assert list(dict.fromkeys(table_s2["feature"])) == publication.FEATURES
    assert set(table_s2["group"]) == set(publication.FEATURE_GROUPS)
    assert table_s2["definition"].str.len().gt(0).all()
    assert table_s2["evidence"].str.len().gt(0).all()
    assert "Machine-readable output" not in set(table_s2["feature"])
    workflow = table_s4[table_s4["record_type"] == "workflow_audit"]
    assert len(workflow) == 4
    assert set(workflow["status"]) == {"N/A", "INCOMPARABLE"}
    assert workflow["reason"].str.len().gt(0).all()

    publication_files = [
        supplementary / "Table_S1_datasets_inputs_databases.tsv",
        supplementary / "Table_S2_capability_evidence.tsv",
        supplementary / "Table_S3_versions_commands_access.tsv",
        supplementary / "Table_S4_case_metrics_failures.tsv",
        supplementary / "source_data" / "Data_S1_full_case_commands.tsv",
        paper / "AllEnricher_Bioinformatics_Application_Note_DRAFT.md",
        paper / "oup_submission" / "AllEnricher_Application_Note.tex",
    ]
    forbidden = re.compile(r"(?i)(?:[A-Z]:\\|/mnt/|file://|zhang_i5edc0)")
    leaked = [str(path) for path in publication_files if forbidden.search(path.read_text(encoding="utf-8"))]
    assert leaked == []