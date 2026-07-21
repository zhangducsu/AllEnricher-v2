#!/usr/bin/env python3
"""Run the preregistered AllEnricher competitor benchmark."""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import math
import os
import platform
import shutil
import subprocess
import sys
import urllib.request
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable

import numpy as np
import pandas as pd
import yaml
from scipy.stats import spearmanr

from allenricher import __version__ as allenricher_version
from allenricher.core.bioconductor import windows_to_wsl_path
from allenricher.core.config import Config
from allenricher.core.enrichment import EnrichmentAnalyzer


PROJECT_ROOT = Path(__file__).resolve().parents[2]
SCRIPT_ROOT = Path(__file__).resolve().parent
DEFAULT_CONFIG = SCRIPT_ROOT / "benchmark_matrix.yaml"
R_ADAPTER = SCRIPT_ROOT / "competitor_methods.R"
NORMALIZED_COLUMNS = [
    "tool", "tool_version", "dataset", "species", "database", "method",
    "term_id", "term_name", "p_value", "adjusted_p_value", "es", "nes",
    "leading_edge", "overlap_count", "term_size", "status",
]
METRIC_COLUMNS = [
    "reference_tool", "comparator", "dataset", "database", "method",
    "reference_terms", "comparator_terms", "term_jaccard", "spearman",
    "max_abs_p_diff", "max_abs_q_diff", "max_abs_nes_diff",
    "sign_concordance", "significant_jaccard", "top20_jaccard", "status",
    "reason",
]
DETAIL_COLUMNS = [
    "reference_tool", "comparator", "dataset", "database", "method",
    "common_terms", "p_spearman", "median_abs_p_diff", "median_abs_q_diff",
    "valid_p_pairs", "valid_q_pairs", "valid_nes_pairs", "median_abs_nes_diff",
    "positive_significant_jaccard", "negative_significant_jaccard",
    "leading_edge_jaccard",
]


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def write_json(path: Path, value: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def open_text(path: Path):
    return gzip.open(path, "rt", encoding="utf-8") if path.suffix == ".gz" else path.open("r", encoding="utf-8")


def read_gene_list(path: Path) -> list[str]:
    return list(dict.fromkeys(line.strip() for line in path.read_text(encoding="utf-8").splitlines() if line.strip()))


def read_ranked(path: Path) -> list[tuple[str, float]]:
    frame = pd.read_csv(path, sep="\t")
    if frame.shape[1] < 2:
        raise ValueError(f"ranked input requires two columns: {path}")
    rows: list[tuple[str, float]] = []
    seen: set[str] = set()
    for gene, value in zip(frame.iloc[:, 0].astype(str), pd.to_numeric(frame.iloc[:, 1], errors="raise")):
        gene = gene.strip()
        if gene and gene not in seen and math.isfinite(float(value)):
            rows.append((gene, float(value)))
            seen.add(gene)
    rows.sort(key=lambda item: item[1], reverse=True)
    return rows


def read_gmt(path: Path) -> dict[str, dict[str, Any]]:
    terms: dict[str, dict[str, Any]] = {}
    with open_text(path) as stream:
        for line_number, line in enumerate(stream, start=1):
            fields = line.rstrip("\r\n").split("\t")
            if len(fields) < 3:
                continue
            term_id, term_name = fields[:2]
            if term_id in terms:
                raise ValueError(f"duplicate GMT term {term_id!r} at line {line_number}: {path}")
            terms[term_id] = {
                "name": term_name or term_id,
                "genes": list(dict.fromkeys(gene.strip() for gene in fields[2:] if gene.strip())),
            }
    if not terms:
        raise ValueError(f"GMT has no valid terms: {path}")
    return terms


def normalize_gmt(
    source: Path,
    target: Path,
    allowed_genes: set[str],
    minimum: int,
    maximum: int | None,
) -> dict[str, dict[str, Any]]:
    normalized: dict[str, dict[str, Any]] = {}
    for term_id, info in read_gmt(source).items():
        genes = sorted(set(map(str, info["genes"])) & allowed_genes)
        if len(genes) < minimum or (maximum is not None and len(genes) > maximum):
            continue
        normalized[term_id] = {"name": str(info["name"]), "genes": genes}
    target.parent.mkdir(parents=True, exist_ok=True)
    with target.open("w", encoding="utf-8", newline="\n") as stream:
        for term_id, info in normalized.items():
            stream.write("\t".join([term_id, info["name"], *info["genes"]]) + "\n")
    if not normalized:
        raise ValueError(f"no GMT terms remain after normalization: {source}")
    return normalized


def case_inputs(input_root: Path, dataset: str) -> dict[str, Path]:
    converted = input_root / dataset / "converted"
    return {
        "query": converted / "query_genes.txt",
        "background": converted / "background_genes.txt",
        "ranked": converted / "ranked_genes.tsv",
        "expression": converted / "expression_counts.tsv",
        "provenance": input_root / dataset / "provenance.json",
    }


def expected_result_sets(config: dict[str, Any]) -> list[dict[str, str]]:
    result: list[dict[str, str]] = []
    for dataset in config["datasets"]:
        for database in config["databases"]:
            for method in ("ORA", "GSEA"):
                for tool in config["methods"][method]["tools"]:
                    result.append({"dataset": dataset, "database": database, "method": method, "tool": tool})
        result.append({"dataset": dataset, "database": "GO", "method": "ORA", "tool": "g:Profiler"})
    return result


def case_id(spec: dict[str, str]) -> str:
    safe_tool = spec["tool"].replace(":", "").replace(" ", "_")
    return "__".join([spec["dataset"], spec["database"], spec["method"], safe_tool])


def find_rscript() -> str:
    configured = os.environ.get("ALLENRICHER_RSCRIPT")
    candidates = [configured, shutil.which("Rscript")]
    for candidate in candidates:
        if candidate and Path(candidate).is_file():
            return str(candidate)
    raise RuntimeError("Rscript was not found; set ALLENRICHER_RSCRIPT")


def package_version(package: str, wsl: bool = False) -> str:
    try:
        command = [find_rscript(), "-e", f"cat(as.character(packageVersion('{package}')))" ]
        if wsl:
            command = ["wsl.exe", "-d", "Ubuntu", "Rscript", "-e", f"cat(as.character(packageVersion('{package}')))" ]
        result = subprocess.run(
            command,
            capture_output=True, text=True, encoding="utf-8", errors="replace", timeout=30,
        )
        return result.stdout.strip() if result.returncode == 0 else "unavailable"
    except Exception:
        return "unavailable"


def resolved_tool_versions() -> dict[str, str]:
    return {
        "AllEnricher": allenricher_version,
        "clusterProfiler": package_version("clusterProfiler"),
        "WebGestaltR": package_version("WebGestaltR", wsl=True),
        "g:Profiler": "service-recorded",
        "getENRICH": "git-commit-recorded",
        "GSVA": package_version("GSVA"),
    }


def blank_normalized(spec: dict[str, str], version: str, status: str) -> pd.DataFrame:
    row = {column: "" for column in NORMALIZED_COLUMNS}
    row.update(spec)
    row.update({"tool_version": version, "species": spec.get("species", ""), "status": status})
    return pd.DataFrame([row], columns=NORMALIZED_COLUMNS)


def normalize_allenricher(
    frame: pd.DataFrame,
    spec: dict[str, str],
    version: str,
    term_data: dict[str, dict[str, Any]],
) -> pd.DataFrame:
    rows = []
    if spec["method"] == "ORA":
        for _, item in frame.iterrows():
            rows.append({
                "term_id": item["Term_ID"], "term_name": item["Term_Name"],
                "p_value": item["P_Value"], "adjusted_p_value": item["Adjusted_P_Value"],
                "es": np.nan, "nes": np.nan, "leading_edge": "",
                "overlap_count": item["Gene_Count"], "term_size": item["Background_Count"],
            })
    else:
        names = {term_id: info["name"] for term_id, info in term_data.items()}
        for _, item in frame.iterrows():
            term_id = str(item.get("Term_ID", item.get("pathway", "")))
            leading = item.get("leadingEdge", item.get("Lead_genes", ""))
            if isinstance(leading, (list, tuple, np.ndarray)):
                leading = ";".join(map(str, leading))
            rows.append({
                "term_id": term_id, "term_name": names.get(term_id, term_id),
                "p_value": item.get("pval", item.get("p_value", np.nan)),
                "adjusted_p_value": item.get("padj", item.get("FDR", np.nan)),
                "es": item.get("ES", np.nan), "nes": item.get("NES", np.nan),
                "leading_edge": leading, "overlap_count": np.nan,
                "term_size": item.get("size", item.get("setSize", np.nan)),
            })
    result = pd.DataFrame(rows)
    for key in ("tool", "dataset", "species", "database", "method"):
        result[key] = spec[key]
    result["tool_version"] = version
    result["status"] = "PASS"
    return result[NORMALIZED_COLUMNS]


def run_allenricher(
    spec: dict[str, str],
    terms: dict[str, dict[str, Any]],
    inputs: dict[str, Path],
) -> pd.DataFrame:
    method = "hypergeometric" if spec["method"] == "ORA" else "gsea"
    config = Config(
        species=spec["species"], databases=[spec["database"]], method=method,
        correction="BH", min_genes=3, output_all=True, n_jobs=1,
        gsea_min_size=15, gsea_max_size=500,
    )
    analyzer = EnrichmentAnalyzer(config)
    query = set(read_gene_list(inputs["query"])) if method == "hypergeometric" else set()
    background = set(read_gene_list(inputs["background"])) if method == "hypergeometric" else set()
    ranked = read_ranked(inputs["ranked"]) if method == "gsea" else None
    result = analyzer.run_analysis(
        query, background, {spec["database"]: terms}, parallel=False, ranked_gene_list=ranked
    )
    return result.get(spec["database"], pd.DataFrame())


def execute(command: list[str], case_dir: Path, timeout: int) -> subprocess.CompletedProcess[str]:
    case_dir.mkdir(parents=True, exist_ok=True)
    (case_dir / "command.txt").write_text(subprocess.list2cmdline(command) + "\n", encoding="utf-8")
    result = subprocess.run(
        command, cwd=PROJECT_ROOT, capture_output=True, text=True,
        encoding="utf-8", errors="replace", timeout=timeout,
    )
    (case_dir / "stdout.log").write_text(result.stdout, encoding="utf-8")
    (case_dir / "stderr.log").write_text(result.stderr, encoding="utf-8")
    (case_dir / "exit_code.txt").write_text(f"{result.returncode}\n", encoding="ascii")
    return result


def run_r_competitor(
    spec: dict[str, str],
    gmt_path: Path,
    inputs: dict[str, Path],
    raw_dir: Path,
    config: dict[str, Any],
) -> pd.DataFrame:
    output = raw_dir / "result.tsv"
    arguments = [
        "--tool", spec["tool"], "--method", spec["method"],
        "--gmt", str(gmt_path), "--query", str(inputs["query"]),
        "--background", str(inputs["background"]), "--ranked", str(inputs["ranked"]),
        "--output", str(output), "--session-info", str(raw_dir / "sessionInfo.txt"),
        "--min-size", str(config["methods"][spec["method"]]["min_gene_set_size"]),
        "--max-size", str(config["methods"].get("GSEA", {}).get("max_gene_set_size", 500)),
        "--seed", str(config["methods"].get("GSEA", {}).get("seed", 42)),
        "--permutations", str(config["methods"].get("GSEA", {}).get("webgestalt_permutations", 10000)),
    ]
    if spec["tool"] == "WebGestaltR":
        path_flags = {"--gmt", "--query", "--background", "--ranked", "--output", "--session-info"}
        translated = []
        translate_next = False
        for value in arguments:
            translated.append(windows_to_wsl_path(value) if translate_next else value)
            translate_next = value in path_flags
        command = ["wsl.exe", "-d", "Ubuntu", "Rscript", windows_to_wsl_path(R_ADAPTER), *translated]
    else:
        command = [find_rscript(), str(R_ADAPTER), *arguments]
    result = execute(command, raw_dir, timeout=3600)
    if result.returncode != 0:
        raise RuntimeError(f"{spec['tool']} exited with code {result.returncode}")
    return pd.read_csv(output, sep="\t")


def normalize_r_result(frame: pd.DataFrame, spec: dict[str, str], version: str) -> pd.DataFrame:
    result = frame.copy()
    detected_version = version
    if "tool_version" in result and not result["tool_version"].dropna().empty:
        detected_version = str(result["tool_version"].dropna().iloc[0])
    for column in ("p_value", "adjusted_p_value", "es", "nes", "overlap_count", "term_size"):
        if column not in result:
            result[column] = np.nan
    for column in ("term_id", "term_name", "leading_edge"):
        if column not in result:
            result[column] = ""
    for key in ("tool", "dataset", "species", "database", "method"):
        result[key] = spec[key]
    result["tool_version"] = detected_version
    result["status"] = "PASS"
    return result[NORMALIZED_COLUMNS]


def http_json(url: str, payload: dict[str, Any], timeout: int = 180) -> dict[str, Any]:
    request = urllib.request.Request(
        url, data=json.dumps(payload).encode("utf-8"), method="POST",
        headers={"Content-Type": "application/json", "User-Agent": "AllEnricher-benchmark/1"},
    )
    with urllib.request.urlopen(request, timeout=timeout) as response:
        return json.loads(response.read().decode("utf-8"))


def normalize_gprofiler_response(
    response: dict[str, Any], spec: dict[str, str], fallback_version: str
) -> pd.DataFrame:
    meta = response.get("meta", {})
    service_version = str(meta.get("version", meta.get("query_metadata", {}).get("data_version", fallback_version)))
    rows = []
    for item in response.get("result", []):
        rows.append({
            "term_id": item.get("native", ""), "term_name": item.get("name", ""),
            "p_value": np.nan, "adjusted_p_value": item.get("p_value", np.nan),
            "es": np.nan, "nes": np.nan, "leading_edge": "",
            "overlap_count": item.get("intersection_size", np.nan),
            "term_size": item.get("term_size", np.nan),
        })
    return normalize_r_result(pd.DataFrame(rows), spec, service_version)


def run_gprofiler(
    spec: dict[str, str],
    gmt_path: Path,
    inputs: dict[str, Path],
    raw_dir: Path,
    base_url: str,
) -> tuple[pd.DataFrame, str]:
    raw_dir.mkdir(parents=True, exist_ok=True)
    upload_payload = {"gmt": gmt_path.read_text(encoding="utf-8"), "name": gmt_path.name}
    write_json(raw_dir / "upload_request.json", upload_payload)
    upload = http_json(f"{base_url.rstrip('/')}/api/gost/custom/", upload_payload)
    write_json(raw_dir / "upload_response.json", upload)
    token = str(upload["organism"])
    query_payload = {
        "organism": token, "query": read_gene_list(inputs["query"]), "sources": None,
        "user_threshold": 1.0, "all_results": True, "ordered": False,
        "no_evidences": True, "combined": False, "measure_underrepresentation": False,
        "no_iea": False, "domain_scope": "custom", "numeric_ns": "",
        "significance_threshold_method": "fdr", "background": read_gene_list(inputs["background"]),
        "output": "json",
    }
    write_json(raw_dir / "profile_request.json", query_payload)
    response = http_json(f"{base_url.rstrip('/')}/api/gost/profile/", query_payload)
    write_json(raw_dir / "profile_response.json", response)
    return normalize_gprofiler_response(response, spec, "service-recorded"), token


def jaccard(left: Iterable[str], right: Iterable[str]) -> float:
    left_set, right_set = set(left), set(right)
    union = left_set | right_set
    return len(left_set & right_set) / len(union) if union else 1.0


def safe_spearman(left: pd.Series, right: pd.Series) -> float:
    valid = pd.to_numeric(left, errors="coerce").notna() & pd.to_numeric(right, errors="coerce").notna()
    if valid.sum() < 2:
        return np.nan
    return float(spearmanr(pd.to_numeric(left[valid]), pd.to_numeric(right[valid])).statistic)


def leading_sets(values: pd.Series) -> list[set[str]]:
    return [set(filter(None, str(value).replace(",", ";").split(";"))) for value in values]


def compare_pair(reference: pd.DataFrame, comparator: pd.DataFrame, spec: dict[str, str]) -> tuple[dict[str, Any], dict[str, Any]]:
    reference = reference[(reference["status"] == "PASS") & reference["term_id"].astype(str).ne("")].copy()
    comparator = comparator[(comparator["status"] == "PASS") & comparator["term_id"].astype(str).ne("")].copy()
    base = {
        "reference_tool": "AllEnricher", "comparator": spec["tool"],
        "dataset": spec["dataset"], "database": spec["database"], "method": spec["method"],
    }
    if reference.empty or comparator.empty:
        metric = {**base, **{column: np.nan for column in METRIC_COLUMNS if column not in base}}
        metric.update({
            "reference_terms": len(reference), "comparator_terms": len(comparator),
            "status": "UNAVAILABLE", "reason": "reference or comparator result is unavailable",
        })
        detail = {**base, **{column: np.nan for column in DETAIL_COLUMNS if column not in base}}
        detail["common_terms"] = 0
        return metric, detail
    merged = reference.merge(comparator, on="term_id", suffixes=("_ref", "_cmp"))
    significant_ref = set(reference.loc[pd.to_numeric(reference["adjusted_p_value"], errors="coerce") < 0.05, "term_id"])
    significant_cmp = set(comparator.loc[pd.to_numeric(comparator["adjusted_p_value"], errors="coerce") < 0.05, "term_id"])
    top_ref = set(reference.sort_values(["adjusted_p_value", "p_value", "term_id"], na_position="last").head(20)["term_id"])
    top_cmp = set(comparator.sort_values(["adjusted_p_value", "p_value", "term_id"], na_position="last").head(20)["term_id"])
    metric = {
        **base, "reference_terms": len(reference), "comparator_terms": len(comparator),
        "term_jaccard": jaccard(reference["term_id"], comparator["term_id"]),
        "spearman": np.nan, "max_abs_p_diff": np.nan, "max_abs_q_diff": np.nan,
        "max_abs_nes_diff": np.nan, "sign_concordance": np.nan,
        "significant_jaccard": jaccard(significant_ref, significant_cmp),
        "top20_jaccard": jaccard(top_ref, top_cmp), "status": "PASS", "reason": "",
    }
    detail = {**base, **{column: np.nan for column in DETAIL_COLUMNS if column not in base}}
    detail["common_terms"] = len(merged)
    p_ref = pd.to_numeric(merged["p_value_ref"], errors="coerce")
    p_cmp = pd.to_numeric(merged["p_value_cmp"], errors="coerce")
    q_ref = pd.to_numeric(merged["adjusted_p_value_ref"], errors="coerce")
    q_cmp = pd.to_numeric(merged["adjusted_p_value_cmp"], errors="coerce")
    valid_p, valid_q = p_ref.notna() & p_cmp.notna(), q_ref.notna() & q_cmp.notna()
    detail["valid_p_pairs"] = int(valid_p.sum())
    detail["valid_q_pairs"] = int(valid_q.sum())
    if valid_p.any():
        p_differences = np.abs(p_ref[valid_p] - p_cmp[valid_p])
        metric["max_abs_p_diff"] = float(p_differences.max())
        detail["median_abs_p_diff"] = float(p_differences.median())
        detail["p_spearman"] = safe_spearman(p_ref, p_cmp)
    if valid_q.any():
        q_differences = np.abs(q_ref[valid_q] - q_cmp[valid_q])
        metric["max_abs_q_diff"] = float(q_differences.max())
        detail["median_abs_q_diff"] = float(q_differences.median())
    if spec["method"] == "ORA":
        if valid_q.any():
            metric["spearman"] = safe_spearman(-np.log10(q_ref.clip(lower=np.finfo(float).tiny)), -np.log10(q_cmp.clip(lower=np.finfo(float).tiny)))
    else:
        nes_ref = pd.to_numeric(merged["nes_ref"], errors="coerce")
        nes_cmp = pd.to_numeric(merged["nes_cmp"], errors="coerce")
        valid = nes_ref.notna() & nes_cmp.notna()
        detail["valid_nes_pairs"] = int(valid.sum())
        if valid.any():
            differences = np.abs(nes_ref[valid] - nes_cmp[valid])
            metric["max_abs_nes_diff"] = float(differences.max())
            detail["median_abs_nes_diff"] = float(differences.median())
            metric["spearman"] = safe_spearman(nes_ref, nes_cmp)
            metric["sign_concordance"] = float((np.sign(nes_ref[valid]) == np.sign(nes_cmp[valid])).mean())
        for direction, key in ((1, "positive_significant_jaccard"), (-1, "negative_significant_jaccard")):
            left = set(reference.loc[(pd.to_numeric(reference["adjusted_p_value"], errors="coerce") < 0.05) & (np.sign(pd.to_numeric(reference["nes"], errors="coerce")) == direction), "term_id"])
            right = set(comparator.loc[(pd.to_numeric(comparator["adjusted_p_value"], errors="coerce") < 0.05) & (np.sign(pd.to_numeric(comparator["nes"], errors="coerce")) == direction), "term_id"])
            detail[key] = jaccard(left, right)
        common_sig = merged[(pd.to_numeric(merged["adjusted_p_value_ref"], errors="coerce") < 0.05) & (pd.to_numeric(merged["adjusted_p_value_cmp"], errors="coerce") < 0.05)]
        if not common_sig.empty:
            scores = [jaccard(left, right) for left, right in zip(leading_sets(common_sig["leading_edge_ref"]), leading_sets(common_sig["leading_edge_cmp"]))]
            detail["leading_edge_jaccard"] = float(np.median(scores))
    return metric, detail


def build_metrics(results: pd.DataFrame, expected: list[dict[str, str]], acceptance: dict[str, float]) -> tuple[pd.DataFrame, pd.DataFrame]:
    metrics, details = [], []
    for spec in expected:
        if spec["tool"] == "AllEnricher":
            continue
        selector = (
            (results["dataset"] == spec["dataset"]) & (results["database"] == spec["database"])
            & (results["method"] == spec["method"])
        )
        reference = results[selector & (results["tool"] == "AllEnricher")]
        comparator = results[selector & (results["tool"] == spec["tool"])]
        metric, detail = compare_pair(reference, comparator, spec)
        if metric["status"] == "PASS" and spec["tool"] in {"WebGestaltR", "g:Profiler"} and spec["method"] == "ORA":
            metric["status"] = "INCOMPARABLE"
            metric["reason"] = (
                "The comparator restricts the ORA universe to genes annotated in the custom GMT; "
                "AllEnricher used the complete detected background."
            )
        elif metric["status"] == "PASS" and spec["tool"] == "clusterProfiler" and spec["method"] == "ORA":
            metric["status"] = "FAIL_NUMERIC" if (
                pd.isna(metric["max_abs_p_diff"])
                or metric["max_abs_p_diff"] > acceptance["ora_max_abs_p_diff"]
                or pd.isna(metric["max_abs_q_diff"])
                or metric["max_abs_q_diff"] > acceptance["ora_max_abs_q_diff"]
            ) else "PASS"
            if metric["status"] != "PASS":
                metric["reason"] = (
                    "Raw ORA P values or BH-FDR values did not meet the preregistered threshold "
                    "under the shared positive-overlap testing family."
                )
        elif metric["status"] == "PASS" and spec["tool"] == "clusterProfiler" and spec["method"] == "GSEA":
            metric["status"] = "FAIL_NUMERIC" if (
                pd.isna(metric["max_abs_nes_diff"]) or metric["max_abs_nes_diff"] > acceptance["clusterprofiler_gsea_max_abs_nes_diff"]
                or pd.isna(metric["max_abs_q_diff"]) or metric["max_abs_q_diff"] > acceptance["clusterprofiler_gsea_max_abs_q_diff"]
            ) else "PASS"
            if metric["status"] != "PASS":
                metric["reason"] = (
                    "clusterProfiler 4.20.0 delegates GSEA to enrichit::gsea_gson; the preregistered "
                    "1e-8 equality threshold against fgseaMultilevel was not met."
                )
        elif metric["status"] == "PASS" and spec["tool"] == "WebGestaltR" and spec["method"] == "GSEA":
            metric["status"] = "DESCRIPTIVE"
            metric["reason"] = "WebGestaltR uses a distinct 10,000-permutation GSEA implementation; no post-hoc pass threshold was applied."
        metrics.append(metric)
        details.append(detail)
    return pd.DataFrame(metrics, columns=METRIC_COLUMNS), pd.DataFrame(details, columns=DETAIL_COLUMNS)


def input_statistics(config: dict[str, Any], input_root: Path, database_root: Path) -> pd.DataFrame:
    rows = []
    snapshot = config["annotation_snapshot"]
    for dataset, metadata in config["datasets"].items():
        inputs = case_inputs(input_root, dataset)
        query = set(read_gene_list(inputs["query"]))
        background = set(read_gene_list(inputs["background"]))
        ranked = {gene for gene, _ in read_ranked(inputs["ranked"])}
        for database in config["databases"]:
            source = database_root / "organism" / snapshot / metadata["species"] / f"{metadata['species']}.{database}.gmt.gz"
            terms = read_gmt(source)
            annotated = {gene for info in terms.values() for gene in info["genes"]}
            rows.append({
                "dataset": dataset, "accession": metadata["accession"], "species": metadata["species"],
                "database": database, "query_genes": len(query), "background_genes": len(background),
                "ranked_genes": len(ranked), "query_mapping_rate": len(query & annotated) / len(query) if query else np.nan,
                "background_mapping_rate": len(background & annotated) / len(background) if background else np.nan,
                "source_gmt": str(source.resolve()), "source_gmt_sha256": sha256(source),
            })
    return pd.DataFrame(rows)


def audit_getenrich(config: dict[str, Any], run_dir: Path, input_root: Path) -> pd.DataFrame:
    tool = config["tools"]["getENRICH"]
    rows = []
    audit_root = run_dir / "raw" / "getENRICH"
    audit_root.mkdir(parents=True, exist_ok=True)
    docker = shutil.which("docker")
    docker_check = None
    if docker:
        try:
            docker_check = subprocess.run([docker, "info"], capture_output=True, text=True, timeout=30)
        except Exception as exc:
            docker_check = exc
    available = isinstance(docker_check, subprocess.CompletedProcess) and docker_check.returncode == 0
    runtime_check = None
    if available:
        runtime_command = [docker, "run", "--rm", f"rocker/r-ver:{tool['r_version']}", "R", "--version"]
        (audit_root / "container_command.txt").write_text(subprocess.list2cmdline(runtime_command) + "\n", encoding="utf-8")
        try:
            runtime_check = subprocess.run(runtime_command, capture_output=True, text=True, timeout=600)
            (audit_root / "container_stdout.log").write_text(runtime_check.stdout, encoding="utf-8")
            (audit_root / "container_stderr.log").write_text(runtime_check.stderr, encoding="utf-8")
            (audit_root / "container_exit_code.txt").write_text(f"{runtime_check.returncode}\n", encoding="ascii")
        except Exception as exc:
            runtime_check = exc
            (audit_root / "container_stderr.log").write_text(f"{type(exc).__name__}: {exc}\n", encoding="utf-8")
    runtime_ready = (
        isinstance(runtime_check, subprocess.CompletedProcess) and runtime_check.returncode == 0
        and f"R version {tool['r_version']}" in runtime_check.stdout
    )
    evidence = {
        "repository": tool["repository"],
        "commit": tool["commit"],
        "required_runtime": f"R {tool['r_version']} isolated container",
        "container_image": f"rocker/r-ver:{tool['r_version']}",
        "official_statistical_call": "clusterProfiler::enrichKEGG",
        "official_supported_database": "KEGG",
        "official_identifier_contract": "NCBI Entrez Gene identifiers for a named organism",
        "forbidden_option": "-f",
    }
    write_json(audit_root / "workflow_contract.json", evidence)
    for dataset in config["workflow_audits"]["getENRICH"]["datasets"]:
        for database in config["workflow_audits"]["getENRICH"]["databases"]:
            case = audit_root / f"{dataset}__{database}"
            case.mkdir(parents=True, exist_ok=True)
            if database == "GO":
                status = "N/A"
                reason = "The pinned official getENRICH workflow implements KEGG enrichment only; no GO workflow is present."
            elif not docker:
                status = "UNAVAILABLE"
                reason = "The getENRICH workflow was not run because the Docker CLI is unavailable."
            elif isinstance(docker_check, Exception):
                status = "UNAVAILABLE"
                reason = f"The getENRICH workflow was not run: {type(docker_check).__name__}: {docker_check}"
            elif not available:
                status = "UNAVAILABLE"
                reason = "The getENRICH workflow was not run because the Docker engine is unavailable."
                (case / "docker_stdout.log").write_text(docker_check.stdout, encoding="utf-8")
                (case / "docker_stderr.log").write_text(docker_check.stderr, encoding="utf-8")
            elif not runtime_ready:
                status = "UNAVAILABLE"
                reason = "The isolated getENRICH R runtime preflight did not complete successfully."
            else:
                status = "INCOMPARABLE"
                reason = (
                    "The archived inputs use gene symbols, whereas the pinned getENRICH named-organism workflow requires "
                    "NCBI Entrez Gene identifiers. No unregistered identifier conversion was introduced."
                )
            record = {
                "tool": "getENRICH", "dataset": dataset, "database": database,
                "status": status, "reason": reason, "commit": tool["commit"],
                "required_runtime": evidence["required_runtime"], "forbidden_option": "-f",
                "container_runtime_verified": runtime_ready,
                "input_root": str(case_inputs(input_root, dataset)["query"].parent),
            }
            write_json(case / "status.json", record)
            rows.append({"tool": "getENRICH", "dataset": dataset, "database": database, "status": status, "reason": reason})
    return pd.DataFrame(rows)


def collect_manifest(run_dir: Path, config_path: Path, tool_versions: dict[str, str], expected: list[dict[str, str]]) -> None:
    files = []
    for path in sorted(run_dir.rglob("*")):
        if path.is_file() and path.name != "benchmark_manifest.json":
            files.append({
                "path": path.relative_to(run_dir).as_posix(), "size": path.stat().st_size,
                "sha256": sha256(path),
            })
    write_json(run_dir / "benchmark_manifest.json", {
        "created_at": utc_now(), "immutable_run_directory": str(run_dir.resolve()),
        "config": str(config_path.resolve()), "config_sha256": sha256(config_path),
        "python": sys.version, "platform": platform.platform(), "tool_versions": tool_versions,
        "expected_result_sets": len(expected), "files": files,
    })


def rebuild_from_raw(config_path: Path, source_run: Path, output: Path) -> Path:
    if output.exists():
        raise FileExistsError(f"run directory already exists and will not be overwritten: {output}")
    for required in (source_run / "raw", source_run / "inputs", source_run / "benchmark_manifest.json"):
        if not required.exists():
            raise FileNotFoundError(required)
    output.mkdir(parents=True)
    shutil.copytree(source_run / "raw", output / "raw", ignore=shutil.ignore_patterns("getENRICH"))
    shutil.copytree(source_run / "inputs", output / "inputs")
    if (source_run / "00_metadata").is_dir():
        shutil.copytree(source_run / "00_metadata", output / "00_metadata")
    else:
        (output / "00_metadata").mkdir()

    config = yaml.safe_load(config_path.read_text(encoding="utf-8"))
    expected = expected_result_sets(config)
    input_root = (config_path.parent / config["input_root"]).resolve()
    database_root = (config_path.parent / config["database_root"]).resolve()
    source_manifest_path = source_run / "benchmark_manifest.json"
    source_manifest = json.loads(source_manifest_path.read_text(encoding="utf-8"))
    versions = dict(source_manifest.get("tool_versions", resolved_tool_versions()))
    versions["getENRICH"] = config["tools"]["getENRICH"]["commit"]
    write_json(output / "00_metadata" / "offline_replay.json", {
        "created_at": utc_now(), "source_run": str(source_run.resolve()),
        "source_manifest_sha256": sha256(source_manifest_path), "network_access": "not used",
    })
    (output / "00_metadata" / "benchmark_matrix.yaml").write_text(
        config_path.read_text(encoding="utf-8"), encoding="utf-8"
    )

    all_results = []
    for spec_base in expected:
        metadata = config["datasets"][spec_base["dataset"]]
        spec = {**spec_base, "species": metadata["species"]}
        raw_dir = output / "raw" / case_id(spec)
        status_path = raw_dir / "status.json"
        status = json.loads(status_path.read_text(encoding="utf-8")) if status_path.is_file() else {"status": "FAIL"}
        try:
            if status.get("status") != "PASS":
                raise RuntimeError(status.get("reason", "archived case did not pass"))
            if spec["tool"] == "AllEnricher":
                frame = pd.read_csv(raw_dir / "result.tsv", sep="\t")
                terms = read_gmt(output / "inputs" / spec["dataset"] / spec["database"] / spec["method"] / "normalized.gmt")
                normalized = normalize_allenricher(frame, spec, versions["AllEnricher"], terms)
            elif spec["tool"] in {"clusterProfiler", "WebGestaltR"}:
                frame = pd.read_csv(raw_dir / "result.tsv", sep="\t")
                normalized = normalize_r_result(frame, spec, versions[spec["tool"]])
            else:
                response = json.loads((raw_dir / "profile_response.json").read_text(encoding="utf-8"))
                normalized = normalize_gprofiler_response(response, spec, versions["g:Profiler"])
            if normalized.empty:
                normalized = blank_normalized(spec, versions.get(spec["tool"], "unknown"), "EMPTY")
        except Exception as exc:
            write_json(raw_dir / "offline_replay_error.json", {
                "error_type": type(exc).__name__, "reason": str(exc), "created_at": utc_now(),
            })
            normalized = blank_normalized(spec, versions.get(spec["tool"], "unknown"), "FAIL")
        all_results.append(normalized)

    normalized_results = pd.concat(all_results, ignore_index=True)
    normalized_results.to_csv(output / "normalized_results.tsv", sep="\t", index=False, lineterminator="\n")
    metrics, details = build_metrics(normalized_results, expected, config["acceptance"])
    metrics.to_csv(output / "benchmark_metrics.tsv", sep="\t", index=False, lineterminator="\n")
    details.to_csv(output / "benchmark_metrics_detail.tsv", sep="\t", index=False, lineterminator="\n")
    input_statistics(config, input_root, database_root).to_csv(output / "input_statistics.tsv", sep="\t", index=False, lineterminator="\n")
    audit_getenrich(config, output, input_root).to_csv(output / "getenrich_workflow_audit.tsv", sep="\t", index=False, lineterminator="\n")
    collect_manifest(output, config_path, versions, expected)
    return output


def run(config_path: Path, output: Path, selected_case: str | None) -> Path:
    if output.exists():
        raise FileExistsError(f"run directory already exists and will not be overwritten: {output}")
    output.mkdir(parents=True)
    config = yaml.safe_load(config_path.read_text(encoding="utf-8"))
    input_root = (SCRIPT_ROOT / config["input_root"]).resolve()
    database_root = (SCRIPT_ROOT / config["database_root"]).resolve()
    expected = expected_result_sets(config)
    if len(expected) != int(config["expected_result_sets"]):
        raise ValueError(f"expected matrix size {config['expected_result_sets']}, generated {len(expected)}")
    if selected_case:
        expected = [spec for spec in expected if case_id(spec) == selected_case]
        if not expected:
            raise KeyError(f"unknown case: {selected_case}")
    versions = resolved_tool_versions()
    versions["getENRICH"] = config["tools"]["getENRICH"]["commit"]
    write_json(output / "00_metadata" / "environment.json", {
        "created_at": utc_now(), "python": sys.version, "platform": platform.platform(),
        "rscript": find_rscript(), "tool_versions": versions,
        "git_commit": subprocess.run(["git", "rev-parse", "HEAD"], cwd=PROJECT_ROOT, capture_output=True, text=True).stdout.strip(),
    })
    (output / "00_metadata" / "benchmark_matrix.yaml").write_text(config_path.read_text(encoding="utf-8"), encoding="utf-8")
    all_results = []
    tokens = []
    for spec_base in expected:
        metadata = config["datasets"][spec_base["dataset"]]
        spec = {**spec_base, "species": metadata["species"]}
        inputs = case_inputs(input_root, spec["dataset"])
        for name, path in inputs.items():
            if name != "expression" and not path.is_file():
                raise FileNotFoundError(path)
        allowed = set(read_gene_list(inputs["background"])) if spec["method"] == "ORA" else {gene for gene, _ in read_ranked(inputs["ranked"])}
        settings = config["methods"][spec["method"]]
        source_gmt = database_root / "organism" / config["annotation_snapshot"] / spec["species"] / f"{spec['species']}.{spec['database']}.gmt.gz"
        gmt_path = output / "inputs" / spec["dataset"] / spec["database"] / spec["method"] / "normalized.gmt"
        if gmt_path.exists():
            terms = read_gmt(gmt_path)
        else:
            terms = normalize_gmt(source_gmt, gmt_path, allowed, int(settings["min_gene_set_size"]), settings.get("max_gene_set_size"))
            write_json(gmt_path.parent / "input_manifest.json", {
                "dataset": spec["dataset"], "database": spec["database"], "method": spec["method"],
                "source_gmt": str(source_gmt), "source_gmt_sha256": sha256(source_gmt),
                "normalized_gmt_sha256": sha256(gmt_path), "normalized_terms": len(terms),
                "query_sha256": sha256(inputs["query"]), "background_sha256": sha256(inputs["background"]),
                "ranked_sha256": sha256(inputs["ranked"]),
            })
        raw_dir = output / "raw" / case_id(spec)
        raw_dir.mkdir(parents=True, exist_ok=True)
        try:
            if spec["tool"] == "AllEnricher":
                frame = run_allenricher(spec, terms, inputs)
                frame.to_csv(raw_dir / "result.tsv", sep="\t", index=False)
                normalized = normalize_allenricher(frame, spec, versions["AllEnricher"], terms)
            elif spec["tool"] in {"clusterProfiler", "WebGestaltR"}:
                frame = run_r_competitor(spec, gmt_path, inputs, raw_dir, config)
                normalized = normalize_r_result(frame, spec, versions[spec["tool"]])
            else:
                normalized, token = run_gprofiler(spec, gmt_path, inputs, raw_dir, config["tools"]["g:Profiler"]["base_url"])
                tokens.append({"dataset": spec["dataset"], "database": spec["database"], "token": token, "accessed_at": utc_now()})
            if normalized.empty:
                normalized = blank_normalized(spec, versions.get(spec["tool"], "unknown"), "EMPTY")
            write_json(raw_dir / "status.json", {"case_id": case_id(spec), "status": "PASS", "rows": len(normalized), "completed_at": utc_now()})
        except Exception as exc:
            write_json(raw_dir / "status.json", {
                "case_id": case_id(spec), "status": "FAIL", "completed_at": utc_now(),
                "error_type": type(exc).__name__, "reason": str(exc),
            })
            normalized = blank_normalized(spec, versions.get(spec["tool"], "unknown"), "FAIL")
        all_results.append(normalized)
    normalized_results = pd.concat(all_results, ignore_index=True) if all_results else pd.DataFrame(columns=NORMALIZED_COLUMNS)
    passed = normalized_results[(normalized_results["status"] == "PASS") & normalized_results["tool_version"].astype(str).ne("")]
    for tool, group in passed.groupby("tool"):
        detected = str(group["tool_version"].iloc[0])
        if detected != "unavailable":
            versions[tool] = detected
    normalized_results.to_csv(output / "normalized_results.tsv", sep="\t", index=False, lineterminator="\n")
    metrics, details = build_metrics(normalized_results, expected, config["acceptance"])
    metrics.to_csv(output / "benchmark_metrics.tsv", sep="\t", index=False, lineterminator="\n")
    details.to_csv(output / "benchmark_metrics_detail.tsv", sep="\t", index=False, lineterminator="\n")
    input_statistics(config, input_root, database_root).to_csv(output / "input_statistics.tsv", sep="\t", index=False, lineterminator="\n")
    if tokens:
        write_json(output / "00_metadata" / "gprofiler_tokens.json", tokens)
    if not selected_case:
        audit_getenrich(config, output, input_root).to_csv(output / "getenrich_workflow_audit.tsv", sep="\t", index=False, lineterminator="\n")
    collect_manifest(output, config_path, versions, expected)
    return output


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--case", help="Run one exact case ID")
    parser.add_argument("--from-raw-run", type=Path, help="Rebuild derived tables offline from an archived run")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    try:
        if args.from_raw_run:
            if args.case:
                raise ValueError("--case cannot be combined with --from-raw-run")
            rebuild_from_raw(args.config.resolve(), args.from_raw_run.resolve(), args.output.resolve())
        else:
            run(args.config.resolve(), args.output.resolve(), args.case)
    except Exception as exc:
        print(f"ERROR: {type(exc).__name__}: {exc}", file=sys.stderr)
        return 1
    print(args.output.resolve())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
