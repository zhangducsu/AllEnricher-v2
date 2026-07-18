#!/usr/bin/env python3
"""Show the input contracts and configuration for GSEA, ssGSEA, and GSVA.

AllEnricher consumes an existing ranked-gene statistic for GSEA and an existing
gene-by-sample expression matrix for ssGSEA or GSVA. This example intentionally
does not perform upstream differential-expression analysis.
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd

from allenricher.core.config import Config


def validate_ranked_gene_table(path: str | Path) -> pd.DataFrame:
    """Load and validate a two-column GSEA ranked-gene table."""
    table = pd.read_csv(path, sep="\t")
    required = {"gene", "weight"}
    if not required.issubset(table.columns):
        raise ValueError("The ranked table must contain 'gene' and 'weight' columns")
    table = table.loc[:, ["gene", "weight"]].dropna()
    table["gene"] = table["gene"].astype(str).str.strip()
    table["weight"] = pd.to_numeric(table["weight"], errors="raise")
    if table["gene"].eq("").any():
        raise ValueError("Gene identifiers must not be blank")
    if table["gene"].duplicated().any():
        raise ValueError("Gene identifiers must be unique in a GSEA ranked table")
    return table.sort_values("weight", ascending=False, kind="mergesort")


def gsea_config() -> Config:
    """Create a GSEA configuration using the publication-oriented defaults."""
    return Config(
        method="gsea",
        species="hsa",
        databases=["GO", "KEGG"],
        gsea_permutations=1000,
        gsea_min_size=15,
        gsea_max_size=500,
        output_dir="./results/gsea",
    )


def ssgsea_config() -> Config:
    """Create an ssGSEA configuration for a gene-by-sample matrix."""
    return Config(
        method="ssgsea",
        species="hsa",
        databases=["GO", "KEGG"],
        gsea_min_size=1,
        gsea_max_size=None,
        output_dir="./results/ssgsea",
    )


def gsva_config() -> Config:
    """Create a GSVA configuration for count data."""
    return Config(
        method="gsva",
        species="hsa",
        databases=["GO", "KEGG"],
        gsva_method="gsva",
        gsva_kcdf="Poisson",
        gsea_min_size=1,
        gsea_max_size=None,
        output_dir="./results/gsva",
    )


def print_cli_examples() -> None:
    """Print equivalent CLI commands."""
    print("GSEA from an existing ranked statistic:")
    print(
        "  allenricher analyze -m gsea -r ranked_genes.tsv "
        "-s hsa -d GO,KEGG -o results/gsea"
    )
    print("\nssGSEA from an expression matrix:")
    print(
        "  allenricher analyze -m ssgsea -e expression_matrix.tsv "
        "-s hsa -d GO,KEGG -o results/ssgsea"
    )
    print("\nGSVA from a count matrix:")
    print(
        "  allenricher analyze -m gsva -e count_matrix.tsv "
        "-s hsa -d GO,KEGG -o results/gsva"
    )


if __name__ == "__main__":
    print("AllEnricher ranked-list and pathway-activity examples")
    print(f"GSEA method: {gsea_config().method}")
    print(f"ssGSEA method: {ssgsea_config().method}")
    print(f"GSVA kernel: {gsva_config().gsva_kcdf}")
    print_cli_examples()
