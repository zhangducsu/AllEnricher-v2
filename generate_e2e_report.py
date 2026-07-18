"""Generate a compact HTML summary from saved enrichment E2E outputs.

The script reports existing evidence only. It does not rerun an analysis or
infer biological conclusions from the result values.
"""

from __future__ import annotations

import html
import json
from pathlib import Path

import numpy as np
import pandas as pd


TEST_DATA_DIR = Path("test_data")
REPORT_JSON = TEST_DATA_DIR / "e2e_test_report.json"
OUTPUT_HTML = TEST_DATA_DIR / "e2e_test_report.html"


def _load_json(path: Path) -> dict:
    if not path.exists():
        raise FileNotFoundError(f"Missing E2E report metadata: {path}")
    return json.loads(path.read_text(encoding="utf-8"))


def _load_optional_matrix(name: str) -> pd.DataFrame | None:
    path = TEST_DATA_DIR / name
    return pd.read_csv(path, index_col=0) if path.exists() else None


def _format_value(value: object) -> str:
    if isinstance(value, (float, np.floating)):
        return "N/A" if np.isnan(value) else f"{value:.4g}"
    return html.escape(str(value))


def _mapping_table(rows: list[tuple[str, object]]) -> str:
    body = "".join(
        f"<tr><th>{html.escape(label)}</th><td>{_format_value(value)}</td></tr>"
        for label, value in rows
    )
    return f'<table class="summary-table"><tbody>{body}</tbody></table>'


def _matrix_table(matrix: pd.DataFrame | None, title: str) -> str:
    if matrix is None:
        return f'<section><h3>{html.escape(title)}</h3><p>No saved matrix was found.</p></section>'

    header = "".join(f"<th>{html.escape(str(column))}</th>" for column in matrix.columns)
    rows: list[str] = []
    for index, values in matrix.iterrows():
        cells = "".join(f"<td>{_format_value(value)}</td>" for value in values)
        rows.append(f"<tr><th>{html.escape(str(index))}</th>{cells}</tr>")

    return f"""
    <section>
      <h3>{html.escape(title)}</h3>
      <p>{matrix.shape[0]} gene sets by {matrix.shape[1]} samples</p>
      <div class="table-scroll">
        <table><thead><tr><th>Gene set</th>{header}</tr></thead><tbody>{''.join(rows)}</tbody></table>
      </div>
    </section>
    """


def _method_summary(report_data: dict) -> str:
    result_blocks: list[str] = []
    for method, details in report_data.get("results", {}).items():
        if not isinstance(details, dict):
            continue
        rows = [(key.replace("_", " ").title(), value) for key, value in details.items()]
        result_blocks.append(
            f"<section><h3>{html.escape(str(method))}</h3>{_mapping_table(rows)}</section>"
        )
    return "".join(result_blocks) or "<p>No method summary was recorded.</p>"


def build_html(report_data: dict) -> str:
    test_data = report_data.get("test_data", {})
    metadata_rows = [("Test date", report_data.get("test_date", "Not recorded"))]
    if isinstance(test_data, dict):
        metadata_rows.extend(
            (key.replace("_", " ").title(), value) for key, value in test_data.items()
        )

    matrices = [
        _matrix_table(_load_optional_matrix("ssgsea_results.csv"), "ssGSEA activity matrix"),
        _matrix_table(_load_optional_matrix("gsva_results.csv"), "GSVA activity matrix"),
        _matrix_table(_load_optional_matrix("gsva_plage_results.csv"), "PLAGE activity matrix"),
        _matrix_table(_load_optional_matrix("gsva_zscore_results.csv"), "Z-score activity matrix"),
    ]

    return f"""<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <title>AllEnricher enrichment E2E report</title>
  <style>
    :root {{ color-scheme: light; --ink: #1f2933; --muted: #617080; --line: #d8dee6; --accent: #2c5282; }}
    * {{ box-sizing: border-box; }}
    body {{ margin: 0; background: #f4f6f8; color: var(--ink); font: 15px/1.55 "Segoe UI", Arial, sans-serif; }}
    main {{ width: min(1180px, calc(100% - 32px)); margin: 24px auto; background: #fff; border: 1px solid var(--line); }}
    header, section {{ padding: 22px 26px; border-bottom: 1px solid var(--line); }}
    h1, h2, h3 {{ margin: 0 0 12px; font-family: Georgia, serif; }}
    h1 {{ color: var(--accent); font-size: 28px; }}
    h2 {{ font-size: 22px; }}
    h3 {{ font-size: 18px; }}
    p {{ margin: 6px 0; color: var(--muted); }}
    table {{ width: 100%; border-collapse: collapse; font-variant-numeric: tabular-nums; }}
    th, td {{ padding: 8px 10px; border: 1px solid var(--line); text-align: left; vertical-align: top; }}
    thead th, .summary-table th {{ background: #eef3f8; color: #243b53; }}
    .summary-table th {{ width: 220px; }}
    .table-scroll {{ overflow-x: auto; }}
    footer {{ padding: 18px 26px; color: var(--muted); }}
  </style>
</head>
<body>
  <main>
    <header>
      <h1>AllEnricher enrichment E2E report</h1>
      <p>Saved execution metadata and numerical outputs for manual review.</p>
    </header>
    <section><h2>Execution metadata</h2>{_mapping_table(metadata_rows)}</section>
    <section><h2>Method summaries</h2>{_method_summary(report_data)}</section>
    {''.join(matrices)}
    <footer>This report describes saved test evidence only; it does not interpret biological results.</footer>
  </main>
</body>
</html>
"""


def main() -> Path:
    report_data = _load_json(REPORT_JSON)
    OUTPUT_HTML.write_text(build_html(report_data), encoding="utf-8")
    print(f"HTML report generated: {OUTPUT_HTML}")
    return OUTPUT_HTML


if __name__ == "__main__":
    main()
