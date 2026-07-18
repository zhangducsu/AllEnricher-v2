"""Run a compact end-to-end check for GSEA, ssGSEA, and GSVA."""

from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd


sys.path.insert(0, str(Path(__file__).parent))

from allenricher.core.enrichment import GSEA, SSGSEA  # noqa: E402
from allenricher.core.gsva import GSVA  # noqa: E402


TEST_DATA_DIR = Path("test_data")


def load_test_data() -> tuple[list[str], dict[str, float], pd.DataFrame, dict[str, set[str]]]:
    """Load the ranked list, expression matrix, and GMT gene sets."""
    ranked_df = pd.read_csv(TEST_DATA_DIR / "ranked_genes.tsv", sep="\t")
    ranked_genes = ranked_df["gene"].astype(str).tolist()
    gene_weights = dict(zip(ranked_df["gene"].astype(str), ranked_df["weight"]))
    expression = pd.read_csv(
        TEST_DATA_DIR / "expression_matrix.tsv", sep="\t", index_col=0
    )

    gene_sets: dict[str, set[str]] = {}
    with (TEST_DATA_DIR / "gene_sets.gmt").open(encoding="utf-8") as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) >= 3:
                gene_sets[fields[0]] = set(fields[2:])

    print(f"Ranked gene list: {len(ranked_genes)} genes")
    print(f"Expression matrix: {expression.shape[0]} genes x {expression.shape[1]} samples")
    print(f"Gene-set collection: {len(gene_sets)} sets")
    return ranked_genes, gene_weights, expression, gene_sets


def test_gsea(
    ranked_genes: list[str],
    gene_weights: dict[str, float],
    gene_sets: dict[str, set[str]],
) -> pd.DataFrame:
    """Run a compact GSEA check on the first five gene sets."""
    started = time.perf_counter()
    gsea = GSEA(permutations=100)
    selected = dict(list(gene_sets.items())[:5])
    rows = []
    for name, genes in selected.items():
        es, nes, pvalue, leading_edge = gsea.calculate_normalized_es(
            ranked_genes, genes, gene_weights
        )
        rows.append(
            {
                "term_id": name,
                "term_name": name,
                "es": es,
                "nes": nes,
                "pvalue": pvalue,
                "leading_edge": leading_edge,
                "gene_count": len(genes),
            }
        )

    results = pd.DataFrame(rows)
    print(f"GSEA completed in {time.perf_counter() - started:.2f}s")
    if not results.empty:
        print(results[["term_name", "es", "nes", "pvalue", "gene_count"]].to_string(index=False))
    return results


def test_activity_method(
    method: str,
    expression: pd.DataFrame,
    gene_sets: dict[str, set[str]],
) -> pd.DataFrame:
    """Run one activity-scoring method on the first five gene sets."""
    selected = dict(list(gene_sets.items())[:5])
    started = time.perf_counter()
    if method == "ssgsea":
        analyzer = SSGSEA(min_size=10, max_size=500)
    else:
        analyzer = GSVA(method=method, kcdf="Gaussian", min_size=10, max_size=500)
    results = analyzer.analyze_matrix(expression, selected)
    print(
        f"{method} completed in {time.perf_counter() - started:.2f}s: "
        f"{results.shape[0]} gene sets x {results.shape[1]} samples"
    )
    return results


def generate_report(
    gsea_results: pd.DataFrame,
    activity_results: dict[str, pd.DataFrame],
    expression: pd.DataFrame,
    gene_set_count: int,
) -> dict:
    """Save a machine-readable summary of the executed checks."""
    report: dict[str, object] = {
        "test_date": time.strftime("%Y-%m-%d %H:%M:%S"),
        "test_data": {
            "expression_matrix": f"{expression.shape[0]} genes x {expression.shape[1]} samples",
            "gene_sets": gene_set_count,
            "samples": expression.columns.astype(str).tolist(),
        },
        "results": {},
    }
    results = report["results"]
    assert isinstance(results, dict)
    if not gsea_results.empty:
        results["GSEA"] = {
            "status": "completed",
            "gene_sets_analyzed": len(gsea_results),
            "positive_nes": int((gsea_results["nes"] > 0).sum()),
            "negative_nes": int((gsea_results["nes"] < 0).sum()),
            "nominal_p_below_0_05": int((gsea_results["pvalue"] < 0.05).sum()),
        }
    for method, matrix in activity_results.items():
        results[method] = {
            "status": "completed",
            "shape": list(matrix.shape),
            "score_range": [float(np.nanmin(matrix.values)), float(np.nanmax(matrix.values))],
        }

    path = TEST_DATA_DIR / "e2e_test_report.json"
    path.write_text(json.dumps(report, indent=2), encoding="utf-8")
    print(f"Report saved: {path}")
    return report


def main() -> dict:
    """Run all checks and return the saved summary."""
    print("=" * 60)
    print("GSEA, ssGSEA, and GSVA end-to-end check")
    print("=" * 60)
    ranked_genes, gene_weights, expression, gene_sets = load_test_data()
    gsea_results = test_gsea(ranked_genes, gene_weights, gene_sets)
    activity_results = {
        "ssGSEA": test_activity_method("ssgsea", expression, gene_sets),
        "GSVA": test_activity_method("gsva", expression, gene_sets),
        "PLAGE": test_activity_method("plage", expression, gene_sets),
        "Z-score": test_activity_method("zscore", expression, gene_sets),
    }
    return generate_report(gsea_results, activity_results, expression, len(gene_sets))


if __name__ == "__main__":
    main()
