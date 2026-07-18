#!/usr/bin/env python3
"""Run the preregistered real-world 108-cell AllEnricher E2E matrix."""

from __future__ import annotations

import argparse
import hashlib
import json
import mimetypes
import os
import platform
import shlex
import shutil
import subprocess
import sys
import time
from datetime import datetime
from html import escape
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import yaml
from scipy.stats import hypergeom
from statsmodels.stats.multitest import multipletests

from allenricher.core.bioconductor import FGSEA_COLUMNS, run_fgsea, run_gsva
from allenricher.database.manager import DatabaseManager


PROJECT_ROOT = Path(__file__).resolve().parents[2]
E2E_ROOT = PROJECT_ROOT / "test_e2e_2026"
SCRIPT_ROOT = Path(__file__).resolve().parent
INPUT_ROOT = E2E_ROOT / "00_input_data" / "real_world_sci"
DATABASE_ROOT = INPUT_ROOT / "database_snapshot"
DEFAULT_MATRIX = SCRIPT_ROOT / "case_matrix.yaml"
METHODS = ("hypergeometric", "gsea", "ssgsea", "gsva")
TF_DATABASES = {"TRRUST", "ChEA3", "hTFtarget", "AnimalTFDB"}
ACTIVITY_PLOTS = ("heatmap", "group_comparison", "correlation")
GSEA_R_PLOTS = ("barplot", "lollipop", "ridgeplot", "emapplot", "enrichment", "enrichment2")
PYTHON_APPENDIX_PLOTS = {
    "gsea": ("lollipop", "ridgeplot", "enrichment", "enrichment2"),
    "ssgsea": ACTIVITY_PLOTS,
    "gsva": ACTIVITY_PLOTS,
}


def validate_database_audit(matrix: dict, audit_path: Path) -> None:
    if not audit_path.is_file():
        raise FileNotFoundError(
            "database snapshot is incomplete; run prepare_database_snapshot.py first"
        )
    audit = json.loads(audit_path.read_text(encoding="utf-8"))
    if not isinstance(audit, dict):
        raise ValueError("database audit must be a JSON object")
    expected = set(matrix["datasets"])
    actual = set(audit)
    if actual != expected:
        missing = ", ".join(sorted(expected - actual)) or "none"
        stale = ", ".join(sorted(actual - expected)) or "none"
        raise ValueError(
            "database audit does not match the current matrix "
            f"(missing: {missing}; stale: {stale}); run prepare_database_snapshot.py --offline"
        )


def write_json(path: Path, value: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def file_record(path: Path, base: Path) -> dict[str, Any]:
    return {
        "path": path.relative_to(base).as_posix(),
        "size": path.stat().st_size,
        "sha256": sha256(path),
        "suffix": path.suffix.lower(),
        "mime": mimetypes.guess_type(path.name)[0] or "application/octet-stream",
    }


def manifest(root: Path, base: Path | None = None) -> list[dict[str, Any]]:
    base = base or root
    if not root.exists():
        return []
    return [file_record(path, base) for path in sorted(root.rglob("*")) if path.is_file()]


def capture(command: list[str], timeout: int = 60) -> str:
    try:
        result = subprocess.run(
            command,
            cwd=PROJECT_ROOT,
            capture_output=True,
            text=True,
            encoding="utf-8",
            errors="replace",
            timeout=timeout,
        )
        return (result.stdout + result.stderr).strip()
    except Exception as exc:
        return f"unavailable: {type(exc).__name__}: {exc}"


def command_environment(offline: bool = False) -> dict[str, str]:
    env = os.environ.copy()
    env.update({"PYTHONUTF8": "1", "PYTHONIOENCODING": "utf-8", "MPLBACKEND": "Agg"})
    if offline:
        blocked_proxy = "http://127.0.0.1:9"
        env.update({
            "ALLENRICHER_OFFLINE": "1",
            "HTTP_PROXY": blocked_proxy,
            "HTTPS_PROXY": blocked_proxy,
            "ALL_PROXY": blocked_proxy,
            "NO_PROXY": "",
        })
    return env


def collect_environment(run_dir: Path, matrix_path: Path) -> None:
    rscript = shutil.which("Rscript") or r"D:\AppGallery\App\R-4.6.1\bin\x64\Rscript.exe"
    environment = {
        "timestamp": datetime.now().astimezone().isoformat(),
        "platform": platform.platform(),
        "python": sys.version,
        "pandas": pd.__version__,
        "numpy": np.__version__,
        "git_commit": capture(["git", "rev-parse", "HEAD"]),
        "git_status": capture(["git", "status", "--short"]),
        "R": capture([rscript, "--version"]),
        "pytest": capture([sys.executable, "-m", "pytest", "--version"]),
        "matrix": file_record(matrix_path, matrix_path.parent),
        "database_manifest": (
            file_record(DATABASE_ROOT / "SOURCE_MANIFEST.json", DATABASE_ROOT)
            if (DATABASE_ROOT / "SOURCE_MANIFEST.json").is_file() else None
        ),
    }
    write_json(run_dir / "00_metadata" / "environment.json", environment)
    help_text = capture([sys.executable, "-m", "allenricher", "analyze", "--help"])
    (run_dir / "00_metadata" / "parser_snapshot.txt").write_text(help_text + "\n", encoding="utf-8")


def case_inputs(case_id: str, database: str) -> dict[str, Path]:
    converted = INPUT_ROOT / case_id / "converted"
    return {
        "query": converted / "query_genes.txt",
        "background": converted / "background_genes.txt",
        "ranked": converted / "ranked_genes.tsv",
        "expression": converted / "expression_counts.tsv",
        "groups": converted / "groups.txt",
        "provenance": INPUT_ROOT / case_id / "provenance.json",
        "database_audit": DATABASE_ROOT / "DATABASE_AUDIT.json",
        "database_manifest": DATABASE_ROOT / "SOURCE_MANIFEST.json",
    }


def validate_fixture(case_id: str, paths: dict[str, Path]) -> list[str]:
    errors = [f"missing input: {name}" for name, path in paths.items() if not path.is_file()]
    if errors:
        return errors
    provenance = json.loads(paths["provenance"].read_text(encoding="utf-8"))
    checks = provenance.get("checks", {})
    if float(checks.get("query_mapping_rate", 0)) < 0.70:
        errors.append(f"{case_id}: query mapping rate below 70%")
    if float(checks.get("query_in_rank_rate", 0)) != 1.0:
        errors.append(f"{case_id}: query is not a subset of ranked genes")
    if float(checks.get("query_in_background_rate", 0)) != 1.0:
        errors.append(f"{case_id}: query is not a subset of expression background")
    if int(checks.get("ranked_genes", 0)) < 1000:
        errors.append(f"{case_id}: fewer than 1000 ranked genes")
    if not checks.get("rank_has_positive") or not checks.get("rank_has_negative"):
        errors.append(f"{case_id}: ranking lacks positive or negative values")
    return errors


def load_terms(species: str, database: str) -> dict[str, dict[str, Any]]:
    manager = DatabaseManager(str(DATABASE_ROOT), species)
    manager.load_database(database)
    loaded = manager.get_all_term_data()
    if database not in loaded:
        raise KeyError(f"database manager did not expose {database}: {list(loaded)}")
    return loaded[database]


def read_ranked(path: Path) -> list[tuple[str, float]]:
    frame = pd.read_csv(path, sep="\t")
    return list(zip(frame.iloc[:, 0].astype(str), pd.to_numeric(frame.iloc[:, 1], errors="raise")))


def read_expression(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", index_col=0)


def activity_gene_sets(terms: dict[str, dict[str, Any]]) -> dict[str, set[str]]:
    return {
        str(term_id): set(map(str, info.get("genes", [])))
        for term_id, info in terms.items()
    }


def compare_numeric(
    actual: pd.DataFrame,
    expected: pd.DataFrame,
    keys: list[str],
    numeric: list[str],
    absolute: float,
    relative: float,
) -> list[str]:
    errors: list[str] = []
    actual_keyed = actual.sort_values(keys).reset_index(drop=True)
    expected_keyed = expected.sort_values(keys).reset_index(drop=True)
    if actual_keyed[keys].astype(str).to_dict("records") != expected_keyed[keys].astype(str).to_dict("records"):
        return ["oracle keys differ from CLI output"]
    for column in numeric:
        left = pd.to_numeric(actual_keyed[column], errors="coerce").to_numpy(dtype=float)
        right = pd.to_numeric(expected_keyed[column], errors="coerce").to_numpy(dtype=float)
        if not np.allclose(left, right, atol=absolute, rtol=relative, equal_nan=True):
            delta = float(np.nanmax(np.abs(left - right))) if len(left) else 0.0
            errors.append(f"oracle mismatch in {column}; max abs delta={delta:.3g}")
    return errors


def oracle_ora(
    actual: pd.DataFrame,
    terms: dict[str, dict[str, Any]],
    query_path: Path,
    background_path: Path,
    oracle_dir: Path,
    database: str,
) -> list[str]:
    query = {line.strip() for line in query_path.read_text(encoding="utf-8").splitlines() if line.strip()}
    universe = {
        line.strip()
        for line in background_path.read_text(encoding="utf-8").splitlines()
        if line.strip()
    }
    query &= universe
    rows = []
    min_size, max_size = 3, None
    for term_id, info in terms.items():
        genes = set(map(str, info.get("genes", []))) & universe
        if len(genes) < min_size or (max_size is not None and len(genes) > max_size):
            continue
        overlap = query & genes
        rows.append(
            {
                "Term_ID": str(term_id),
                "P_Value": float(hypergeom.sf(len(overlap) - 1, len(universe), len(genes), len(query))),
                "Gene_Count": len(overlap),
            }
        )
    expected = pd.DataFrame(rows)
    if expected.empty:
        return ["ORA oracle produced no eligible hypotheses"]
    expected["Adjusted_P_Value"] = multipletests(expected["P_Value"], method="fdr_bh")[1]
    expected = expected[expected["Gene_Count"] > 0].copy()
    expected.to_csv(oracle_dir / "independent_hypergeometric.tsv", sep="\t", index=False)
    errors = compare_numeric(
        actual,
        expected,
        ["Term_ID"],
        ["P_Value", "Adjusted_P_Value", "Gene_Count"],
        absolute=1e-12,
        relative=1e-8,
    )
    if actual.get("Term_Name", pd.Series(dtype=str)).astype(str).eq(actual["Term_ID"].astype(str)).all():
        errors.append("all ORA display names equal term IDs")
    return errors


def oracle_gsea(
    actual: pd.DataFrame,
    terms: dict[str, dict[str, Any]],
    ranked_path: Path,
    oracle_dir: Path,
    database: str,
) -> list[str]:
    min_size = 15
    max_size = 5000 if database in TF_DATABASES else 500
    expected = run_fgsea(
        read_ranked(ranked_path),
        {str(term): set(map(str, info.get("genes", []))) for term, info in terms.items()},
        min_size=min_size,
        max_size=max_size,
        seed=42,
    )
    expected.to_csv(oracle_dir / "official_fgsea.tsv", sep="\t", index=False)
    required = {"Term_ID", "Term_Name", *FGSEA_COLUMNS}
    if missing := required - set(actual.columns):
        return [f"GSEA table lacks columns {sorted(missing)}"]
    errors = compare_numeric(
        actual[FGSEA_COLUMNS],
        expected,
        ["pathway"],
        ["pval", "padj", "log2err", "ES", "NES", "size"],
        absolute=1e-12,
        relative=1e-8,
    )
    if not actual["Term_ID"].astype(str).equals(actual["pathway"].astype(str)):
        errors.append("GSEA Term_ID differs from fgsea pathway")
    if actual["Term_Name"].astype(str).eq(actual["Term_ID"].astype(str)).all():
        errors.append("all GSEA display names equal term IDs")
    if any(info.get("hierarchy") or "|" in str(info.get("name") or "") for info in terms.values()):
        if "Hierarchy" not in actual.columns:
            errors.append("hierarchical GSEA table lacks Hierarchy column")
    if actual.empty:
        errors.append("fgsea result is empty")
    return errors


def oracle_activity(
    actual: pd.DataFrame,
    terms: dict[str, dict[str, Any]],
    expression_path: Path,
    method: str,
    oracle_dir: Path,
    database: str,
) -> list[str]:
    expression = read_expression(expression_path)
    expected = run_gsva(
        expression,
        activity_gene_sets(terms),
        method=method,
        kcdf="Poisson" if method == "gsva" else "Gaussian",
        tau=0.25 if method == "ssgsea" else 1.0,
        min_size=1,
        max_size=None,
    )
    expected.to_csv(oracle_dir / f"official_{method}.tsv", sep="\t", index=True, index_label="Term_ID")
    actual_indexed = actual.set_index("Term_ID")[expected.columns]
    errors = []
    if actual_indexed.index.tolist() != expected.index.tolist() or actual_indexed.columns.tolist() != expected.columns.tolist():
        errors.append("GSVA activity matrix labels/order differ from official oracle")
        return errors
    left = actual_indexed.to_numpy(dtype=float)
    right = expected.to_numpy(dtype=float)
    if not np.allclose(left, right, atol=1e-8, rtol=0, equal_nan=True):
        errors.append(f"GSVA activity matrix max abs delta={np.nanmax(np.abs(left - right)):.3g}")
    if actual.empty:
        errors.append(f"{method} result is empty")
    return errors


def log_errors(stdout: str, stderr: str) -> list[str]:
    combined = stdout + "\n" + stderr
    markers = (
        "Traceback (most recent call last)",
        " - ERROR - ",
        " - WARNING - ",
        "RuntimeWarning",
        "Warning message:",
        "findfont",
        "empty distance matrix",
        "\ufffd",
    )
    return [f"log contains {marker!r}" for marker in markers if marker in combined]


def expected_plot_tokens(method: str) -> tuple[str, ...]:
    if method == "hypergeometric":
        return ("barplot", "lollipop")
    if method == "gsea":
        return GSEA_R_PLOTS
    return ACTIVITY_PLOTS


def has_plot_token(names: list[str], token: str) -> bool:
    if token == "enrichment":
        return any("enrichment" in name and "enrichment2" not in name for name in names)
    return any(token in name for name in names)


def validate_case_output(
    output: Path,
    database: str,
    method: str,
    terms: dict[str, dict[str, Any]],
    inputs: dict[str, Path],
    oracle_dir: Path,
    requested_plot_tokens: tuple[str, ...] | None = None,
) -> list[str]:
    errors: list[str] = []
    result_path = output / f"{database}_enrichment.tsv"
    metadata_path = output / "analysis_metadata.json"
    if not result_path.is_file() or not result_path.stat().st_size:
        return [f"missing or empty result: {result_path.name}"]
    if not metadata_path.is_file() or not metadata_path.stat().st_size:
        errors.append("missing analysis_metadata.json")
    else:
        metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
        if metadata.get("analysis_method") != method:
            errors.append(f"metadata method={metadata.get('analysis_method')!r}, expected {method!r}")
        if not metadata.get("database_versions", {}).get(database):
            errors.append(f"metadata lacks version for {database}")
    if result_path.read_text(encoding="utf-8", errors="strict").startswith("#"):
        errors.append("formal result table contains metadata comments")
    actual = pd.read_csv(result_path, sep="\t")
    oracle_dir.mkdir(parents=True, exist_ok=True)
    common_required = {"Term_ID", "Term_Name"}
    if missing := common_required - set(actual.columns):
        errors.append(f"result table lacks columns {sorted(missing)}")
    elif database not in TF_DATABASES:
        invalid_names = (
            actual["Term_Name"].fillna("").astype(str).str.strip().eq("")
            | actual["Term_Name"].astype(str).str.casefold().eq(
                actual["Term_ID"].astype(str).str.casefold()
            )
        )
        if invalid_names.any():
            errors.append(f"{int(invalid_names.sum())} terms lack concrete names")
    has_hierarchy = any(
        info.get("hierarchy") or "|" in str(info.get("name") or "")
        for info in terms.values()
    )
    if has_hierarchy and "Hierarchy" not in actual.columns:
        errors.append("hierarchical result table lacks Hierarchy column")

    if method == "hypergeometric":
        required = {"Term_ID", "Term_Name", "P_Value", "Adjusted_P_Value", "Gene_Count", "Genes"}
        if missing := required - set(actual.columns):
            errors.append(f"ORA table lacks columns {sorted(missing)}")
        elif actual.empty:
            errors.append("ORA result is empty")
        else:
            errors.extend(oracle_ora(
                actual, terms, inputs["query"], inputs["background"], oracle_dir,
                database,
            ))
    elif method == "gsea":
        errors.extend(oracle_gsea(
            actual, terms, inputs["ranked"], oracle_dir, database
        ))
    else:
        required = {"Term_ID", "Term_Name"}
        if missing := required - set(actual.columns):
            errors.append(f"activity table lacks columns {sorted(missing)}")
        else:
            errors.extend(oracle_activity(
                actual, terms, inputs["expression"], method, oracle_dir, database
            ))
            if actual["Term_Name"].astype(str).eq(actual["Term_ID"].astype(str)).all():
                errors.append("all activity display names equal term IDs")
            if any(info.get("hierarchy") or "|" in str(info.get("name") or "") for info in terms.values()):
                if "Hierarchy" not in actual.columns:
                    errors.append("hierarchical activity table lacks Hierarchy column")

    plots = [path for path in output.rglob("*") if path.suffix.lower() in {".png", ".pdf", ".svg"}]
    names = [path.stem.lower() for path in plots]
    for token in requested_plot_tokens or expected_plot_tokens(method):
        if not has_plot_token(names, token):
            errors.append(f"missing requested plot type: {token}")
    errors.extend(f"empty plot: {path.name}" for path in plots if path.stat().st_size == 0)
    report_path = output / "report.html"
    if not report_path.is_file():
        errors.append("missing report.html")
    else:
        report_html = report_path.read_text(encoding="utf-8", errors="strict")
        if "\ufffd" in report_html:
            errors.append("report contains Unicode replacement characters")
        plot_stems = {
            path.relative_to(output).with_suffix("").as_posix()
            for path in plots
        }
        embedded_media = report_html.count("<img ") + report_html.count("<object ")
        if embedded_media < len(plot_stems):
            errors.append(
                f"report embeds {embedded_media}/{len(plot_stems)} generated plot types"
            )
        if not actual.empty and "Term_Name" in actual.columns:
            term_names = actual["Term_Name"].dropna().astype(str).str.strip()
            term_names = term_names[term_names.ne("")]
            if not term_names.empty and escape(term_names.iloc[0], quote=True) not in report_html:
                errors.append("report does not expose result term names")
        method_label = {
            "hypergeometric": "ORA Enrichment Analysis Report",
            "gsea": "GSEA Enrichment Analysis Report",
            "ssgsea": "ssGSEA Pathway Activity Report",
            "gsva": "GSVA Pathway Activity Report",
        }[method]
        if method_label not in report_html:
            errors.append(f"report title does not identify {method}")
    return errors


def shell_command(command: list[str]) -> str:
    return "#!/usr/bin/env bash\nset -euo pipefail\ncd " + shlex.quote(str(PROJECT_ROOT)) + "\n" + shlex.join(command) + "\n"


def build_command(
    method: str,
    dataset: dict[str, Any],
    database: str,
    inputs: dict[str, Path],
    output: Path,
    config_path: Path,
    use_r: bool = True,
) -> list[str]:
    command = [
        sys.executable,
        "-m",
        "allenricher",
        "analyze",
        "-i",
        str(inputs["query"]),
        "-s",
        dataset["species"],
        "-d",
        database,
        "-m",
        method,
        "--database-dir",
        str(DATABASE_ROOT),
        "-o",
        str(output),
        "-p",
        "1",
        "-q",
        "1",
        "--plot-format",
        "png",
        "--plot-dpi",
        "150",
        "--style",
        "nature",
        "-j",
        "1",
    ]
    if method == "gsea":
        plots = GSEA_R_PLOTS if use_r else PYTHON_APPENDIX_PLOTS["gsea"]
        command.extend(["-r", str(inputs["ranked"]), "-pt", ",".join(plots)])
        if use_r:
            command.append("--use-r-plots")
        command.extend(
            [
                "--emapplot-qvalue", "1.0", "--emapplot-min-count", "1", "--emapplot-top-n", "30",
                "--gsea-enrichment-top-up", "5", "--gsea-enrichment-top-down", "5",
                "--gsea-multi-top-up", "3", "--gsea-multi-top-down", "3",
            ]
        )
    elif method == "hypergeometric":
        command.extend(["-b", str(inputs["background"])])
    elif method in {"ssgsea", "gsva"}:
        groups = inputs["groups"].read_text(encoding="utf-8").strip()
        command.extend(
            [
                "-e", str(inputs["expression"]), "--groups", groups,
                "-pt", ",".join(ACTIVITY_PLOTS), "--config", str(config_path),
            ]
        )
        if use_r:
            command.append("--use-r-plots")
    return command


def run_case(
    run_dir: Path,
    case_name: str,
    command: list[str],
    inputs: dict[str, Path],
    database: str,
    method: str,
    terms: dict[str, dict[str, Any]],
    env: dict[str, str],
    timeout: int,
    requested_plot_tokens: tuple[str, ...] | None = None,
) -> dict[str, Any]:
    case_dir = run_dir / "cases" / case_name
    output = case_dir / "output"
    oracle_dir = case_dir / "oracle"
    output.mkdir(parents=True, exist_ok=True)
    (case_dir / "command.txt").write_text(subprocess.list2cmdline(command) + "\n", encoding="utf-8")
    (case_dir / "run.sh").write_text(shell_command(command), encoding="utf-8", newline="\n")
    write_json(
        case_dir / "input_manifest.json",
        [
            {"role": role, "absolute_path": str(path.resolve()), **file_record(path, path.parent)}
            for role, path in inputs.items()
        ],
    )
    started = time.perf_counter()
    try:
        result = subprocess.run(
            command,
            cwd=PROJECT_ROOT,
            env=env,
            capture_output=True,
            text=True,
            encoding="utf-8",
            errors="replace",
            timeout=timeout,
        )
        code, stdout, stderr = result.returncode, result.stdout, result.stderr
    except subprocess.TimeoutExpired as exc:
        code = 124
        stdout = exc.stdout if isinstance(exc.stdout, str) else ""
        stderr = (exc.stderr if isinstance(exc.stderr, str) else "") + f"\nTIMEOUT after {timeout}s\n"
    duration = round(time.perf_counter() - started, 3)
    (case_dir / "stdout.log").write_text(stdout, encoding="utf-8")
    (case_dir / "stderr.log").write_text(stderr, encoding="utf-8")
    (case_dir / "exit_code.txt").write_text(f"{code}\n", encoding="ascii")
    errors = [] if code == 0 else [f"exit code {code}"]
    errors.extend(log_errors(stdout, stderr))
    try:
        errors.extend(
            validate_case_output(
                output, database, method, terms, inputs, oracle_dir, requested_plot_tokens
            )
        )
    except Exception as exc:
        errors.append(f"validator error: {type(exc).__name__}: {exc}")
    output_records = manifest(output, case_dir)
    write_json(case_dir / "output_manifest.json", output_records)
    record = {
        "case_id": case_name,
        "dataset": case_name.split("__", 1)[0],
        "database": database,
        "method": method,
        "status": "PASS" if not errors else "FAIL",
        "exit_code": code,
        "duration_seconds": duration,
        "output_files": len(output_records),
        "errors": errors,
    }
    (case_dir / "review_notes.md").write_text(
        f"# {case_name}\n\n- Status: {record['status']}\n- Duration: {duration}s\n"
        f"- Database: {database}\n- Method: {method}\n- Errors: {errors or 'none'}\n",
        encoding="utf-8",
    )
    print(f"{record['status']:4} {case_name} ({duration:.1f}s)", flush=True)
    return record


def compare_offline(primary: Path, current: Path, records: list[dict[str, Any]]) -> list[dict[str, Any]]:
    comparison = []
    for record in records:
        case_id = record["case_id"]
        if "__" not in case_id:
            continue
        left_dir = primary / "cases" / case_id / "output"
        right_dir = current / "cases" / case_id / "output"
        left = {item["path"].split("output/", 1)[-1]: item["sha256"] for item in manifest(left_dir, primary / "cases" / case_id)}
        right = {item["path"].split("output/", 1)[-1]: item["sha256"] for item in manifest(right_dir, current / "cases" / case_id)}
        table_names = sorted(name for name in set(left) | set(right) if name.endswith("_enrichment.tsv"))
        mismatches = [name for name in table_names if left.get(name) != right.get(name)]
        comparison.append({"case_id": case_id, "status": "PASS" if not mismatches else "FAIL", "mismatches": mismatches})
    write_json(current / "OFFLINE_REPRODUCIBILITY.json", comparison)
    return comparison


def compare_public_go_custom(run_dir: Path, matrix: dict[str, Any]) -> list[str]:
    records = []
    errors = []
    for dataset_id, dataset in matrix["datasets"].items():
        if "PUBLIC_GO_CUSTOM" not in dataset["databases"]:
            continue
        for method in METHODS:
            standard = run_dir / "cases" / f"{dataset_id}__GO__{method}" / "output" / "GO_enrichment.tsv"
            custom = (
                run_dir / "cases" / f"{dataset_id}__PUBLIC_GO_CUSTOM__{method}" /
                "output" / "PUBLIC_GO_CUSTOM_enrichment.tsv"
            )
            mismatch = ""
            if not standard.is_file() or not custom.is_file():
                mismatch = "missing result table"
            else:
                left = pd.read_csv(standard, sep="\t").sort_values("Term_ID").reset_index(drop=True)
                right = pd.read_csv(custom, sep="\t").sort_values("Term_ID").reset_index(drop=True)
                if method == "hypergeometric":
                    left = left.drop(columns=["Database", "Term_URL"], errors="ignore")
                    right = right.drop(columns=["Database", "Term_URL"], errors="ignore")
                try:
                    pd.testing.assert_frame_equal(left, right, check_exact=False, rtol=1e-12, atol=1e-12)
                except AssertionError as exc:
                    mismatch = str(exc).splitlines()[0]
            status = "PASS" if not mismatch else "FAIL"
            records.append({"dataset": dataset_id, "method": method, "status": status, "error": mismatch})
            if mismatch:
                errors.append(f"PUBLIC_GO_CUSTOM mismatch {dataset_id}/{method}: {mismatch}")
    pd.DataFrame(records).to_csv(run_dir / "PUBLIC_GO_CUSTOM_EQUIVALENCE.tsv", sep="\t", index=False)
    return errors


def create_contact_sheets(run_dir: Path, matrix: dict[str, Any]) -> list[str]:
    errors: list[str] = []
    for dataset_id in matrix["datasets"]:
        for method in METHODS:
            review_dir = run_dir / "visual_review" / dataset_id / method
            review_dir.mkdir(parents=True, exist_ok=True)
            case_dirs = sorted((run_dir / "cases").glob(f"{dataset_id}__*__{method}"))
            case_dirs.extend(
                sorted((run_dir / "cases").glob(f"APPENDIX__{dataset_id}__*__{method}__python"))
            )
            for case_dir in case_dirs:
                images = sorted((case_dir / "output").rglob("*.png"))
                selected: list[Path] = []
                for token in expected_plot_tokens(method):
                    match = next((path for path in images if token in path.stem.lower() and path not in selected), None)
                    if match:
                        selected.append(match)
                parts = case_dir.name.split("__")
                case_label = f"PYTHON_{parts[2]}" if parts[0] == "APPENDIX" else parts[1]
                for index, source in enumerate(selected, 1):
                    shutil.copy2(source, review_dir / f"{case_label}_{index:02d}_{source.name}")
            images = list(review_dir.glob("*.png"))
            if not images:
                errors.append(f"no images for contact sheet {dataset_id}/{method}")
                continue
            output = run_dir / "visual_review" / f"CONTACT_SHEET_{dataset_id}_{method}.png"
            command = [
                sys.executable,
                str(E2E_ROOT / "create_visual_contact_sheet.py"),
                str(review_dir),
                str(output),
                "--columns", "5" if len(images) > 30 else "3",
                "--title", f"{dataset_id} - {method}",
            ]
            completed = subprocess.run(
                command, cwd=PROJECT_ROOT, capture_output=True, text=True,
                encoding="utf-8", errors="replace", timeout=600,
            )
            if completed.returncode:
                errors.append(f"contact sheet failed {dataset_id}/{method}: {completed.stderr.strip()}")
    return errors


def write_paper_materials(
    run_dir: Path,
    matrix: dict[str, Any],
    records: list[dict[str, Any]],
    appendix_records: list[dict[str, Any]] | None = None,
) -> None:
    paper = run_dir / "paper_materials"
    paper.mkdir(parents=True, exist_ok=True)
    appendix_records = appendix_records or []
    data_rows = []
    for dataset_id, spec in matrix["datasets"].items():
        provenance = json.loads((INPUT_ROOT / dataset_id / "provenance.json").read_text(encoding="utf-8"))
        data_rows.append(
            {
                "dataset": dataset_id,
                "accession": spec["accession"],
                "species": spec["latin_name"],
                "contrast": spec["contrast"],
                "query_genes": provenance["checks"]["query_genes"],
                "ranked_genes": provenance["checks"]["ranked_genes"],
                "samples": provenance["checks"]["reference_samples"] + provenance["checks"]["test_samples"],
                "source_url": provenance["source_url"],
                "license": provenance["license"],
            }
        )
    pd.DataFrame(data_rows).to_csv(paper / "DATA_SOURCES.tsv", sep="\t", index=False)
    pd.DataFrame([*records, *appendix_records]).drop(columns=["errors"], errors="ignore").to_csv(
        paper / "RUNTIME_PERFORMANCE.tsv", sep="\t", index=False
    )
    (paper / "METHODS.md").write_text(
        "# Methods\n\n"
        "Four preregistered Expression Atlas differential RNA-seq experiments were tested without "
        "local differential-expression modelling. Atlas-filtered gene lists were used for ORA, full "
        "Atlas log2 fold changes for fgsea ranking, and public count matrices for ssGSEA/GSVA. Each "
        "species-database pair was analyzed independently with hypergeometric ORA, Bioconductor "
        "fgseaMultilevel, GSVA::ssgseaParam, and GSVA::gsvaParam (Poisson kcdf). Numerical results "
        "were compared with independent SciPy/statsmodels or direct Bioconductor reruns.\n",
        encoding="utf-8",
    )
    passed = sum(record["status"] == "PASS" for record in records)
    appendix_passed = sum(record["status"] == "PASS" for record in appendix_records)
    appendix_result = (
        f" Python visualization appendix: {appendix_passed}/{len(appendix_records)} passed."
        if appendix_records else ""
    )
    (paper / "RESULTS.md").write_text(
        f"# Results\n\nMain analysis matrix: {passed}/{len(records)} passed."
        f"{appendix_result} "
        "See `MATRIX_COVERAGE.tsv`, per-case logs, numerical oracles, and visual contact sheets.\n",
        encoding="utf-8",
    )
    (paper / "LIMITATIONS.md").write_text(
        "# Limitations\n\nThese E2E analyses validate software behavior on four public bulk RNA-seq "
        "experiments but do not establish biological generalizability. DisGeNET uses the frozen "
        "AllEnricher-v1 free snapshot because later full releases are not freely downloadable. "
        "Expression Atlas upstream processing and gene identifier mappings remain external inputs.\n",
        encoding="utf-8",
    )
    (paper / "REFERENCES.bib").write_text(
        "@misc{expression_atlas, title={EMBL-EBI Expression Atlas}, "
        "url={https://www.ebi.ac.uk/gxa/}, note={Accessed 2026-07-15}}\n"
        "@article{fgsea, title={Fast gene set enrichment analysis}, author={Korotkevich, Gennady and others}, year={2021}}\n"
        "@article{gsva, title={GSVA: gene set variation analysis for microarray and RNA-seq data}, "
        "author={Hanzelmann, Sonja and Castelo, Robert and Guinney, Justin}, journal={BMC Bioinformatics}, year={2013}, volume={14}, pages={7}}\n",
        encoding="utf-8",
    )


def write_summary(run_dir: Path, records: list[dict[str, Any]], expected: int, extras: list[str]) -> None:
    passed = sum(record["status"] == "PASS" for record in records)
    failed = sum(record["status"] == "FAIL" for record in records)
    summary = {"expected_main_cases": expected, "executed": len(records), "pass": passed, "fail": failed, "issues": extras, "cases": records}
    write_json(run_dir / "E2E_SUMMARY.json", summary)
    lines = [
        "# Real-world SCI E2E summary", "",
        f"- Expected main matrix: {expected}", f"- Executed: {len(records)}",
        f"- PASS: {passed}", f"- FAIL: {failed}", f"- Framework issues: {len(extras)}", "",
        "| Case | Database | Method | Status | Seconds | Errors |", "|---|---|---|---|---:|---|",
    ]
    lines.extend(
        f"| {row['case_id']} | {row['database']} | {row['method']} | {row['status']} | "
        f"{row['duration_seconds']} | {'; '.join(row['errors']) or '-'} |" for row in records
    )
    if extras:
        lines.extend(["", "## Framework issues", "", *(f"- {item}" for item in extras)])
    (run_dir / "E2E_SUMMARY.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def run_python_appendix(
    run_dir: Path,
    matrix: dict[str, Any],
    config_path: Path,
    env: dict[str, str],
    timeout: int,
) -> list[dict[str, Any]]:
    records: list[dict[str, Any]] = []
    for dataset_id, dataset in matrix["datasets"].items():
        database = "GO"
        inputs = case_inputs(dataset_id, database)
        terms = load_terms(dataset["species"], database)
        for method in ("gsea", "ssgsea", "gsva"):
            case_name = f"APPENDIX__{dataset_id}__{database}__{method}__python"
            case_dir = run_dir / "cases" / case_name
            command = build_command(
                method, dataset, database, inputs, case_dir / "output", config_path, use_r=False
            )
            records.append(
                run_case(
                    run_dir,
                    case_name,
                    command,
                    inputs,
                    database,
                    method,
                    terms,
                    env,
                    timeout,
                    requested_plot_tokens=PYTHON_APPENDIX_PLOTS[method],
                )
            )
    pd.DataFrame(records).to_csv(run_dir / "PYTHON_APPENDIX_COVERAGE.tsv", sep="\t", index=False)
    return records


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--matrix", type=Path, default=DEFAULT_MATRIX)
    parser.add_argument("--mode", choices=("primary", "offline"), default="primary")
    parser.add_argument("--output", type=Path)
    parser.add_argument("--compare-to", type=Path)
    parser.add_argument("--case", help="Exact case id or substring")
    parser.add_argument("--timeout", type=int, default=3600)
    parser.add_argument("--fail-fast", action="store_true")
    parser.add_argument("--skip-visual-audit", action="store_true")
    args = parser.parse_args()

    matrix_path = args.matrix.resolve()
    matrix = yaml.safe_load(matrix_path.read_text(encoding="utf-8"))
    expected = sum(len(spec["databases"]) * len(METHODS) for spec in matrix["datasets"].values())
    if expected != 108:
        parser.error(f"preregistered matrix must contain 108 cells, found {expected}")
    try:
        validate_database_audit(matrix, DATABASE_ROOT / "DATABASE_AUDIT.json")
    except (OSError, ValueError) as exc:
        parser.error(str(exc))
    run_dir = (args.output or (
        E2E_ROOT / "99_runs" / f"{datetime.now():%Y%m%d_%H%M%S_%f}_REAL_WORLD_SCI_{args.mode.upper()}"
    )).resolve()
    run_dir.mkdir(parents=True, exist_ok=True)
    collect_environment(run_dir, matrix_path)
    shutil.copy2(matrix_path, run_dir / "00_metadata" / "case_matrix.yaml")
    config_path = run_dir / "00_metadata" / "runtime_config.yaml"
    config_path.write_text("gsva_kcdf: Poisson\n", encoding="utf-8")

    env = command_environment(offline=args.mode == "offline")
    records: list[dict[str, Any]] = []
    fixture_errors: list[str] = []
    for dataset_id, dataset in matrix["datasets"].items():
        inputs = case_inputs(dataset_id, "")
        fixture_errors.extend(validate_fixture(dataset_id, inputs))
        for database in dataset["databases"]:
            try:
                terms = load_terms(dataset["species"], database)
            except Exception as exc:
                fixture_errors.append(f"{dataset_id}/{database}: database load failed: {exc}")
                if args.fail_fast:
                    break
                continue
            for method in METHODS:
                case_name = f"{dataset_id}__{database}__{method}"
                if args.case and args.case not in case_name:
                    continue
                case_dir = run_dir / "cases" / case_name
                command = build_command(method, dataset, database, inputs, case_dir / "output", config_path)
                record = run_case(
                    run_dir, case_name, command, inputs, database, method, terms, env, args.timeout
                )
                records.append(record)
                if record["status"] == "FAIL" and args.fail_fast:
                    break
            if args.fail_fast and records and records[-1]["status"] == "FAIL":
                break
        if args.fail_fast and records and records[-1]["status"] == "FAIL":
            break

    coverage = pd.DataFrame(records)
    coverage.to_csv(run_dir / "MATRIX_COVERAGE.tsv", sep="\t", index=False)
    if not args.case and len(records) != expected:
        fixture_errors.append(f"matrix incomplete: executed {len(records)}/{expected}")

    if args.mode == "offline":
        if not args.compare_to:
            fixture_errors.append("offline mode requires --compare-to PRIMARY_RUN")
        else:
            comparison = compare_offline(args.compare_to.resolve(), run_dir, records)
            fixture_errors.extend(
                f"offline mismatch: {item['case_id']} {item['mismatches']}"
                for item in comparison if item["status"] == "FAIL"
            )

    if not args.case:
        fixture_errors.extend(compare_public_go_custom(run_dir, matrix))

    appendix_records: list[dict[str, Any]] = []
    if args.mode == "primary" and not args.case:
        appendix_records = run_python_appendix(run_dir, matrix, config_path, env, args.timeout)
        fixture_errors.extend(
            f"Python appendix failed: {record['case_id']}"
            for record in appendix_records if record["status"] == "FAIL"
        )

    if not args.case and not args.skip_visual_audit:
        audit = subprocess.run(
            [sys.executable, str(E2E_ROOT / "audit_e2e_outputs.py"), str(run_dir)],
            cwd=PROJECT_ROOT, capture_output=True, text=True,
            encoding="utf-8", errors="replace", timeout=1800,
        )
        (run_dir / "visual_audit.stdout.log").write_text(audit.stdout, encoding="utf-8")
        (run_dir / "visual_audit.stderr.log").write_text(audit.stderr, encoding="utf-8")
        if audit.returncode:
            fixture_errors.append(f"visual audit exit code {audit.returncode}")
        fixture_errors.extend(create_contact_sheets(run_dir, matrix))

    write_paper_materials(run_dir, matrix, records, appendix_records)
    write_summary(run_dir, records, expected, fixture_errors)
    write_json(run_dir / "run_manifest.json", manifest(run_dir, run_dir))
    print(f"RUN_DIR={run_dir}")
    return 1 if fixture_errors or any(record["status"] == "FAIL" for record in records) else 0


if __name__ == "__main__":
    raise SystemExit(main())
