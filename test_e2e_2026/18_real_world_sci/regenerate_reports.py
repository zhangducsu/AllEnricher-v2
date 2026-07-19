#!/usr/bin/env python3
"""Rebuild the real data with the archived results sheet and graph E2E HTML report."""

from __future__ import annotations

import argparse
import hashlib
import json
import mimetypes
import re
from html import escape
from pathlib import Path

import pandas as pd

from allenricher.core.config import Config
from allenricher.report.generator import ReportGenerator


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def output_manifest(case_dir: Path) -> list[dict[str, object]]:
    records = []
    for path in sorted((case_dir / "output").rglob("*")):
        if not path.is_file():
            continue
        records.append({
            "path": path.relative_to(case_dir).as_posix(),
            "size": path.stat().st_size,
            "sha256": sha256(path),
            "suffix": path.suffix.lower(),
            "mime": mimetypes.guess_type(path.name)[0] or "application/octet-stream",
        })
    return records


def refresh_run_manifest(run_dir: Path) -> None:
    records = []
    for path in sorted(run_dir.rglob("*")):
        if not path.is_file() or path == run_dir / "run_manifest.json":
            continue
        records.append({
            "path": path.relative_to(run_dir).as_posix(),
            "size": path.stat().st_size,
            "sha256": sha256(path),
            "suffix": path.suffix.lower(),
            "mime": mimetypes.guess_type(path.name)[0] or "application/octet-stream",
        })
    (run_dir / "run_manifest.json").write_text(
        json.dumps(records, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )


def query_genes(case_dir: Path) -> list[str]:
    manifest_path = case_dir / "input_manifest.json"
    if not manifest_path.is_file():
        return []
    records = json.loads(manifest_path.read_text(encoding="utf-8"))
    query = next((record for record in records if record.get("role") == "query"), None)
    if not query:
        return []
    path = Path(query["absolute_path"])
    if not path.is_file():
        return []
    return [line.strip().split()[0] for line in path.read_text(encoding="utf-8").splitlines() if line.strip()]


def regenerate_case(case_dir: Path) -> None:
    output = case_dir / "output"
    metadata_path = output / "analysis_metadata.json"
    if not metadata_path.is_file():
        raise FileNotFoundError(f"Missing analytical metadata: {metadata_path}")
    metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
    method = metadata["analysis_method"]
    database = metadata["databases"][0]
    result_path = output / f"{database}_enrichment.tsv"
    if not result_path.is_file():
        raise FileNotFoundError(f"Lack of official results table: {result_path}")

    results = pd.read_csv(result_path, sep="\t")
    config = Config(method=method, pvalue_cutoff=1.0, qvalue_cutoff=1.0)
    ReportGenerator(str(output), config).generate(
        {database: results},
        str(output / "report.html"),
        gene_list=query_genes(case_dir),
        analysis_method=method,
        metadata=metadata,
    )
    (case_dir / "output_manifest.json").write_text(
        json.dumps(output_manifest(case_dir), indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )


def audit_case(case_dir: Path) -> dict[str, object]:
    output = case_dir / "output"
    errors: list[str] = []
    report_path = output / "report.html"
    metadata_path = output / "analysis_metadata.json"
    if not report_path.is_file() or not metadata_path.is_file():
        return {"case": case_dir.name, "status": "FAIL", "errors": ["Lack of metadata to report or analyse"]}

    html = report_path.read_text(encoding="utf-8", errors="strict")
    metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
    method = metadata["analysis_method"]
    database = metadata["databases"][0]
    result_path = output / f"{database}_enrichment.tsv"
    result = pd.read_csv(result_path, sep="\t")
    plots = [
        path for path in output.rglob("*")
        if path.is_file() and path.suffix.lower() in {".png", ".jpg", ".jpeg", ".svg", ".pdf"}
    ]
    plot_types = {
        path.relative_to(output).with_suffix("").as_posix()
        for path in plots
    }
    embedded_media = html.count("<img ") + html.count("<object ")

    expected_title = {
        "hypergeometric": "ORA Enrichment Analysis Report",
        "gsea": "GSEA Enrichment Analysis Report",
        "ssgsea": "ssGSEA Pathway Activity Report",
        "gsva": "GSVA Pathway Activity Report",
    }[method]
    if expected_title not in html:
        errors.append("Report title and method of analysis are inconsistent")
    if "\ufffd" in html:
        errors.append("Organisation")
    if "fonts.googleapis.com" in html:
        errors.append("The report still relies on online fonts")
    if embedded_media != len(plot_types):
        errors.append(f"Charts are not fully embedded: {embedded_media}P/{len(PLOT_types)}")
    if not result.empty:
        if 'class="no-results-box"' in html:
            errors.append("Non-empty results are incorrectly rendered as result pages")
        if "Term_Name" not in result.columns:
            errors.append("Formal results are missing")
        else:
            names = result["Term_Name"].dropna().astype(str).str.strip()
            names = names[names.ne("")]
            if names.empty or escape(names.iloc[0], quote=True) not in html:
                errors.append("The report does not present the name of the route in the official results")
    if "Term ID" not in html or "Term Name" not in html:
        errors.append("The results sheet is missing/TermName Header")
    targets = set(re.findall(r'id="([^"]+)"', html))
    missing_targets = sorted(set(re.findall(r'href="#([^"]+)"', html)) - targets)
    if missing_targets:
        errors.append(f"Unable to anchor the top navigation: {missing_targets}")
    return {
        "case": case_dir.name,
        "method": method,
        "database": database,
        "status": "PASS" if not errors else "FAIL",
        "plot_types": len(plot_types),
        "embedded_media": embedded_media,
        "report_size": report_path.stat().st_size,
        "errors": errors,
    }


def write_audit(run_dir: Path, records: list[dict[str, object]]) -> None:
    failures = [record for record in records if record["status"] == "FAIL"]
    payload = {
        "run": run_dir.name,
        "reports": len(records),
        "pass": len(records) - len(failures),
        "fail": len(failures),
        "records": records,
    }
    (run_dir / "REPORT_AUDIT.json").write_text(
        json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    lines = [
        "# HTML report completeness audit",
        "",
        f"- Number of reports: {len(records)}",
        f"- Through: {len(records) - len(failures)}",
        f"- Failed: {len(failures)}",
        "",
    ]
    if failures:
        lines.extend(["# Lose the item", ""])
        lines.extend(
            f"- `{record['case']}`: {'; '.join(record['errors'])}"
            for record in failures
        )
    else:
        lines.append("All reports are through methodological titles, Term ID/Name, chart embedding, UTF-8, offline font and navigation anchor check.")
    (run_dir / "REPORT_AUDIT.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("run_dirs", nargs="+", type=Path)
    parser.add_argument("--audit-only", action="store_true")
    args = parser.parse_args()
    total = 0
    total_failures = 0
    for run_dir in args.run_dirs:
        run_dir = run_dir.resolve()
        cases_dir = run_dir / "cases"
        case_dirs = sorted(path for path in cases_dir.iterdir() if path.is_dir())
        if not args.audit_only:
            for case_dir in case_dirs:
                regenerate_case(case_dir)
                total += 1
            print(f"{run_dir}: rebuilt {len(case_dirs)} reports")
        records = [audit_case(case_dir) for case_dir in case_dirs]
        write_audit(run_dir, records)
        refresh_run_manifest(run_dir)
        failures = sum(record["status"] == "FAIL" for record in records)
        total_failures += failures
        print(f"{run_dir}: audited {len(records)} reports, failures={failures}")
    if not args.audit_only:
        print(f"total rebuilt reports: {total}")
    return 1 if total_failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
