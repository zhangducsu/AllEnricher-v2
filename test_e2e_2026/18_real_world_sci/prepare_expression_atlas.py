#!/usr/bin/env python3
"""Download and convert fixed Expression Atlas experiments without re-running DE."""

from __future__ import annotations

import argparse
import hashlib
import json
import time
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd
import requests
import yaml


PROJECT_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_MATRIX = Path(__file__).with_name("case_matrix.yaml")
DEFAULT_INPUT_ROOT = PROJECT_ROOT / "test_e2e_2026" / "00_input_data" / "real_world_sci"
RESOURCE_PATHS = {
    "analytics.tsv": "resources/DifferentialSecondaryDataFiles.RnaSeq/analytics",
    "raw_counts.tsv": "resources/DifferentialSecondaryDataFiles.RnaSeq/raw-counts",
    "experiment_design.tsv": "resources/ExperimentDesignFile.RnaSeq/experiment-design",
}
QUERY_FILENAME = "query_results_padj0.05_abslog2fc1.tsv"
DOWNLOAD_ATTEMPTS = 4


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def manifest(root: Path) -> list[dict]:
    return [
        {
            "path": path.relative_to(root).as_posix(),
            "size": path.stat().st_size,
            "sha256": sha256(path),
        }
        for path in sorted(root.rglob("*"))
        if path.is_file() and path.name not in {"source_manifest.json", "input_manifest.json"}
    ]


def _matches_source_manifest(target: Path) -> bool:
    manifest_path = target.parent / "source_manifest.json"
    if not target.is_file() or not manifest_path.is_file():
        return False
    try:
        entries = json.loads(manifest_path.read_text(encoding="utf-8"))
        entry = next(item for item in entries if item.get("path") == target.name)
        return target.stat().st_size == int(entry["size"]) and sha256(target) == entry["sha256"]
    except (json.JSONDecodeError, KeyError, OSError, StopIteration, TypeError, ValueError):
        return False


def download(url: str, target: Path, offline: bool) -> None:
    if _matches_source_manifest(target):
        return
    if offline:
        raise FileNotFoundError(f"offline cache missing or unverified: {target}")
    target.parent.mkdir(parents=True, exist_ok=True)
    temporary = target.with_suffix(target.suffix + ".part")
    session = requests.Session()
    session.trust_env = False
    last_error: Exception | None = None
    try:
        for attempt in range(1, DOWNLOAD_ATTEMPTS + 1):
            temporary.unlink(missing_ok=True)
            try:
                with session.get(url, timeout=180, stream=True) as response:
                    response.raise_for_status()
                    written = 0
                    with temporary.open("wb") as stream:
                        for chunk in response.iter_content(1024 * 1024):
                            if chunk:
                                stream.write(chunk)
                                written += len(chunk)
                    if written == 0:
                        raise OSError("empty response")
                    content_length = response.headers.get("Content-Length")
                    if content_length and not response.headers.get("Content-Encoding"):
                        if written != int(content_length):
                            raise OSError(
                                f"incomplete response: expected {content_length} bytes, got {written}"
                            )
                temporary.replace(target)
                return
            except (OSError, ValueError, requests.RequestException) as exc:
                last_error = exc
                temporary.unlink(missing_ok=True)
                if attempt < DOWNLOAD_ATTEMPTS:
                    time.sleep(2 ** (attempt - 1))
    finally:
        session.close()
    raise RuntimeError(f"download failed after {DOWNLOAD_ATTEMPTS} attempts: {url}") from last_error


def stable_retrieved_at(case_root: Path, accession: str) -> str:
    provenance_path = case_root / "provenance.json"
    if provenance_path.is_file():
        try:
            existing = json.loads(provenance_path.read_text(encoding="utf-8"))
            if existing.get("accession") == accession and existing.get("retrieved_at"):
                return str(existing["retrieved_at"])
        except (OSError, UnicodeError, json.JSONDecodeError):
            pass
    return datetime.now(timezone.utc).isoformat()


def _clean_gene_names(frame: pd.DataFrame) -> pd.Series:
    names = frame["Gene Name"].astype("string").str.strip()
    return names.mask(names.eq(""))


def _analysis_gene_ids(frame: pd.DataFrame) -> pd.Series:
    """Priority is given to gene IDs; the name returns consistently to the stable Gene ID when missing."""
    names = _clean_gene_names(frame)
    gene_ids = frame["Gene ID"].astype("string").str.strip()
    gene_ids = gene_ids.mask(gene_ids.eq(""))
    return names.fillna(gene_ids)


def _map_query_to_analytics_ids(
    query_results: pd.DataFrame,
    analytics: pd.DataFrame,
) -> pd.Series:
    """Use the uniform analytical markings using Gene ID."""
    mapping = analytics.dropna(subset=["Analysis Gene"])[["Gene ID", "Analysis Gene"]].copy()
    mapping["Gene ID"] = mapping["Gene ID"].astype("string").str.strip()
    mapping = mapping.drop_duplicates("Gene ID").set_index("Gene ID")["Analysis Gene"]
    query_ids = query_results["Gene ID"].astype("string").str.strip()
    return query_ids.map(mapping)


def _pick_samples(
    design: pd.DataFrame,
    factor: str,
    reference: str,
    test: str,
    sample_filters: dict[str, str] | None = None,
) -> tuple[list[str], list[str]]:
    if factor not in design.columns:
        raise ValueError(f"missing design factor: {factor}")
    selected = (
        design["Analysed"].astype(str).str.casefold().eq("yes")
        if "Analysed" in design
        else pd.Series(True, index=design.index)
    )
    for column, value in (sample_filters or {}).items():
        if column not in design.columns:
            raise ValueError(f"missing sample filter column: {column}")
        selected &= design[column].astype(str).eq(str(value))
    reference_runs = design.loc[
        selected & design[factor].astype(str).eq(reference), "Run"
    ].astype(str).tolist()
    test_runs = design.loc[
        selected & design[factor].astype(str).eq(test), "Run"
    ].astype(str).tolist()
    if len(reference_runs) < 3 or len(test_runs) < 3:
        raise ValueError(
            f"expected at least three samples per group, got {len(reference_runs)} and {len(test_runs)}"
        )
    return reference_runs, test_runs


def _aggregate_technical_replicates(
    expression: pd.DataFrame,
    design: pd.DataFrame,
    reference_runs: list[str],
    test_runs: list[str],
    unit_column: str | None,
) -> tuple[pd.DataFrame, list[str], list[str]]:
    if not unit_column:
        return expression, reference_runs, test_runs
    if unit_column not in design.columns:
        raise ValueError(f"missing biological unit column: {unit_column}")

    run_to_unit = (
        design.loc[:, ["Run", unit_column]]
        .dropna()
        .astype(str)
        .drop_duplicates("Run")
        .set_index("Run")[unit_column]
    )

    def aggregate_group(runs: list[str], prefix: str) -> tuple[pd.DataFrame, list[str]]:
        missing = [run for run in runs if run not in run_to_unit]
        if missing:
            raise ValueError(f"biological unit missing for runs: {missing}")
        units: dict[str, list[str]] = {}
        for run in runs:
            unit = str(run_to_unit[run]).strip()
            if not unit:
                raise ValueError(f"empty biological unit for run: {run}")
            units.setdefault(unit, []).append(run)
        labels = [f"{prefix}_{unit}" for unit in units]
        aggregated = pd.DataFrame(
            {
                label: expression[runs_for_unit].sum(axis=1)
                for label, runs_for_unit in zip(labels, units.values())
            },
            index=expression.index,
        )
        return aggregated, labels

    reference, reference_labels = aggregate_group(reference_runs, "Control")
    test, test_labels = aggregate_group(test_runs, "Treatment")
    if len(reference_labels) < 3 or len(test_labels) < 3:
        raise ValueError(
            "expected at least three biological units per group after technical-replicate aggregation, "
            f"got {len(reference_labels)} and {len(test_labels)}"
        )
    return pd.concat([reference, test], axis=1), reference_labels, test_labels


def prepare_dataset(case_id: str, spec: dict, input_root: Path, base_url: str, offline: bool) -> dict:
    case_root = input_root / case_id
    source_root = case_root / "source"
    converted_root = case_root / "converted"
    converted_root.mkdir(parents=True, exist_ok=True)
    accession = spec["accession"]
    resource_base = f"{base_url.rstrip('/')}/experiments-content/{accession}/"
    for filename, resource in RESOURCE_PATHS.items():
        download(resource_base + resource, source_root / filename, offline)
    query_url = (
        resource_base
        + "download/RNASEQ_MRNA_DIFFERENTIAL"
        + "?cutoff=0.05&foldChangeCutoff=1&type=RNASEQ_MRNA_DIFFERENTIAL"
    )
    download(query_url, source_root / QUERY_FILENAME, offline)

    analytics = pd.read_csv(source_root / "analytics.tsv", sep="\t", low_memory=False)
    counts = pd.read_csv(source_root / "raw_counts.tsv", sep="\t", low_memory=False)
    design = pd.read_csv(source_root / "experiment_design.tsv", sep="\t", low_memory=False)
    query_results = pd.read_csv(source_root / QUERY_FILENAME, sep="\t", comment="#", low_memory=False)
    contrast = spec["contrast"]
    pvalue_column = f"{contrast}.p-value"
    rank_column = f"{contrast}.log2foldchange"
    missing = [column for column in ("Gene ID", "Gene Name", pvalue_column, rank_column) if column not in analytics]
    if missing:
        raise ValueError(f"{accession} analytics missing columns: {missing}")

    analytics = analytics.copy()
    analytics["Gene Name"] = _clean_gene_names(analytics)
    analytics["Analysis Gene"] = _analysis_gene_ids(analytics)
    analytics[pvalue_column] = pd.to_numeric(analytics[pvalue_column], errors="coerce")
    analytics[rank_column] = pd.to_numeric(analytics[rank_column], errors="coerce")
    mapped = analytics.dropna(subset=["Analysis Gene"]).copy()
    mapped["_abs_rank"] = mapped[rank_column].abs()
    mapped = mapped.sort_values("_abs_rank", ascending=False, kind="mergesort").drop_duplicates("Analysis Gene")
    ranked = mapped.dropna(subset=[rank_column]).sort_values(rank_column, ascending=False, kind="mergesort")

    query_fold_column = f"{contrast} .foldChange"
    query_pvalue_column = f"{contrast}.pValue"
    query_missing = [
        column
        for column in ("Gene ID", "Gene Name", query_fold_column, query_pvalue_column)
        if column not in query_results
    ]
    if query_missing:
        raise ValueError(f"{accession} filtered query result missing columns: {query_missing}")
    query_results = query_results.copy()
    query_results["Analysis Gene"] = _map_query_to_analytics_ids(query_results, analytics)
    query_results[query_fold_column] = pd.to_numeric(query_results[query_fold_column], errors="coerce")
    query_results[query_pvalue_column] = pd.to_numeric(query_results[query_pvalue_column], errors="coerce")
    target_query_rows = query_results.loc[
        query_results[query_fold_column].notna() & query_results[query_pvalue_column].notna()
    ].copy()
    query_mapping_rate = float(target_query_rows["Analysis Gene"].notna().mean()) if len(target_query_rows) else 0.0
    query = target_query_rows.dropna(subset=["Analysis Gene"])["Analysis Gene"].drop_duplicates()

    reference_runs, test_runs = _pick_samples(
        design,
        spec["factor_column"],
        str(spec["reference_value"]),
        str(spec["test_value"]),
        spec.get("sample_filters"),
    )
    samples = reference_runs + test_runs
    missing_samples = [sample for sample in samples if sample not in counts.columns]
    if missing_samples:
        raise ValueError(f"{accession} count matrix missing samples: {missing_samples}")
    counts = counts.copy()
    counts["Analysis Gene"] = _analysis_gene_ids(counts)
    expression = counts.dropna(subset=["Analysis Gene"])[["Analysis Gene", *samples]].copy()
    expression[samples] = expression[samples].apply(pd.to_numeric, errors="coerce").fillna(0)
    expression = expression.groupby("Analysis Gene", sort=True, as_index=True)[samples].sum()
    expression.index.name = "Gene"
    expression, reference_samples, test_samples = _aggregate_technical_replicates(
        expression,
        design,
        reference_runs,
        test_runs,
        spec.get("biological_unit_column"),
    )

    query_genes = set(query.astype(str))
    ranked_genes = set(ranked["Analysis Gene"].astype(str))
    expression_genes = set(expression.index.astype(str))
    query_in_rank = query_genes & ranked_genes
    query_in_background = query_genes & expression_genes

    query_path = converted_root / "query_genes.txt"
    rank_path = converted_root / "ranked_genes.tsv"
    expression_path = converted_root / "expression_counts.tsv"
    background_path = converted_root / "background_genes.txt"
    groups_path = converted_root / "groups.txt"
    mapping_path = converted_root / "id_mapping_audit.tsv"
    query_path.write_text("\n".join(query.astype(str)) + "\n", encoding="utf-8")
    ranked[["Analysis Gene", rank_column]].rename(columns={"Analysis Gene": "gene", rank_column: "weight"}).to_csv(
        rank_path, sep="\t", index=False, lineterminator="\n"
    )
    expression.to_csv(expression_path, sep="\t", lineterminator="\n")
    background_path.write_text("\n".join(expression.index.astype(str)) + "\n", encoding="utf-8")
    groups = f"Control:{','.join(reference_samples)};Treatment:{','.join(test_samples)}"
    groups_path.write_text(groups + "\n", encoding="utf-8")
    analytics[["Gene ID", "Gene Name", "Analysis Gene"]].assign(
        mapping_status=lambda frame: frame["Gene Name"].notna().map(
            {True: "GENE_NAME", False: "GENE_ID_FALLBACK"}
        )
    ).to_csv(mapping_path, sep="\t", index=False, lineterminator="\n")

    source_symbol_rate = float(analytics["Gene Name"].notna().mean())
    checks = {
        "query_source_rows": int(len(target_query_rows)),
        "query_genes": int(query.nunique()),
        "ranked_genes": int(ranked["Analysis Gene"].nunique()),
        "expression_genes": int(len(expression)),
        "reference_samples": len(reference_samples),
        "test_samples": len(test_samples),
        "reference_runs": len(reference_runs),
        "test_runs": len(test_runs),
        "query_mapping_rate": query_mapping_rate,
        "query_in_rank_rate": len(query_in_rank) / len(query_genes) if query_genes else 0.0,
        "query_in_background_rate": len(query_in_background) / len(query_genes) if query_genes else 0.0,
        "analytics_symbol_rate": source_symbol_rate,
        "rank_has_positive": bool((ranked[rank_column] > 0).any()),
        "rank_has_negative": bool((ranked[rank_column] < 0).any()),
    }
    errors = []
    if not checks["query_genes"]:
        errors.append("fixed Atlas cutoffs produced an empty query")
    elif checks["query_genes"] < 10:
        errors.append("fewer than 10 query genes at the fixed Atlas cutoffs")
    if checks["ranked_genes"] < 1000:
        errors.append("fewer than 1000 ranked genes")
    if query_mapping_rate < 0.70:
        errors.append("query gene-name mapping rate below 70%")
    if checks["query_in_rank_rate"] != 1.0:
        errors.append("mapped query is not a subset of the ranked genes")
    if checks["query_in_background_rate"] != 1.0:
        errors.append("mapped query is not a subset of the expression background")
    if not checks["rank_has_positive"] or not checks["rank_has_negative"]:
        errors.append("rank does not contain both signs")
    if errors:
        raise ValueError(f"{case_id} data gates failed: {errors}")

    provenance = {
        "case_id": case_id,
        "accession": accession,
        "retrieved_at": stable_retrieved_at(case_root, accession),
        "source": "EMBL-EBI Expression Atlas",
        "source_url": f"{base_url.rstrip('/')}/experiments/{accession}",
        "license": "CC BY 4.0",
        "contrast": contrast,
        "query_rule": "Expression Atlas filtered TSV: adjusted p-value < 0.05 and abs(log2 fold change) >= 1",
        "rank_rule": "Expression Atlas full analytics log2 fold change; no local DE model",
        "source_urls": {
            "filtered_query": query_url,
            **{filename: resource_base + resource for filename, resource in RESOURCE_PATHS.items()},
        },
        "groups": {"Control": reference_samples, "Treatment": test_samples},
        "source_runs": {"Control": reference_runs, "Treatment": test_runs},
        "biological_unit_column": spec.get("biological_unit_column"),
        "sample_filters": spec.get("sample_filters", {}),
        "checks": checks,
    }
    (case_root / "provenance.json").write_text(
        json.dumps(provenance, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    (source_root / "source_manifest.json").write_text(
        json.dumps(manifest(source_root), ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    (converted_root / "input_manifest.json").write_text(
        json.dumps(manifest(converted_root), ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    return provenance


def prepare_all(matrix_path: Path, input_root: Path, offline: bool = False) -> dict[str, dict]:
    matrix = yaml.safe_load(matrix_path.read_text(encoding="utf-8"))
    base_url = matrix["source"]["base_url"]
    return {
        case_id: prepare_dataset(case_id, spec, input_root, base_url, offline)
        for case_id, spec in matrix["datasets"].items()
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--matrix", type=Path, default=DEFAULT_MATRIX)
    parser.add_argument("--input-root", type=Path, default=DEFAULT_INPUT_ROOT)
    parser.add_argument("--offline", action="store_true")
    args = parser.parse_args()
    summaries = prepare_all(args.matrix.resolve(), args.input_root.resolve(), args.offline)
    print(json.dumps(summaries, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
