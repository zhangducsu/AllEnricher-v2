"""Generate deterministic synthetic fixtures for enrichment smoke tests.

The fixtures exercise GSEA with a ranked gene list and ssGSEA/GSVA with a
gene-by-sample expression matrix. They are not intended as biological evidence.
"""

from __future__ import annotations

import gzip
from pathlib import Path
from typing import Dict, Set

import numpy as np
import pandas as pd


CURATED_HUMAN_SYMBOLS = [
    "TP53", "BRCA1", "BRCA2", "EGFR", "KRAS", "MYC", "BRAF", "AKT1", "PTEN", "RB1",
    "APC", "CTNNB1", "PIK3CA", "ERBB2", "VEGFA", "CDK4", "CDKN2A", "MDM2", "BCL2", "BAX",
    "CCNA2", "CCNB1", "CCND1", "CCNE1", "CDC25A", "CDC25B", "CHEK1", "CHEK2", "AURKA", "AURKB",
    "ATM", "ATR", "RAD51", "PALB2", "MLH1", "MSH2", "PMS2", "PDCD1", "CD274", "CTLA4",
    "LAG3", "HAVCR2", "TIGIT", "IFNG", "IL6", "TNF", "CASP3", "CASP8", "CASP9", "MCL1",
    "BAK1", "PMAIP1", "BBC3", "WNT1", "WNT2", "AXIN1", "AXIN2", "GSK3B", "JUN", "FOS",
    "HK2", "LDHA", "PKM", "PGK1", "ENO1", "GAPDH", "TPI1", "PFKL", "G6PD", "TKT",
    "VEGFB", "VEGFC", "FLT1", "KDR", "FLT4", "ANGPT1", "ANGPT2", "TEK", "PDGFA", "HIF1A",
    "CDH1", "CDH2", "VIM", "SNAI1", "SNAI2", "ZEB1", "ZEB2", "TWIST1", "MMP2", "MMP9",
    "NANOG", "POU5F1", "SOX2", "KLF4", "LIN28A", "LIN28B", "ALDH1A1", "PROM1", "CD44", "ACTB",
    "B2M", "PPIA", "RPLP0", "RPS18", "EEF1A1", "TUBB", "MAPK1", "MAPK3", "MAPK8", "MAPK14",
    "AKT2", "AKT3", "MTOR", "SRC", "ABL1", "STAT1", "STAT3", "STAT5A", "STAT5B", "JAK1",
    "JAK2", "JAK3", "TYK2", "SOCS1", "SOCS3", "VHL", "EP300", "HDAC1", "PPARG", "VDR",
]
EXTRA_GENES = [f"GENE{i:04d}" for i in range(1, 2601)]
ALL_GENES = list(dict.fromkeys([*CURATED_HUMAN_SYMBOLS, *EXTRA_GENES]))


def _rng(seed: int) -> np.random.Generator:
    return np.random.default_rng(seed)


def generate_sorted_gene_list(n_genes: int = 2000, seed: int = 42) -> pd.DataFrame:
    """Generate a unique ranked gene list with synthetic signed weights."""
    if n_genes > len(ALL_GENES):
        raise ValueError(f"Requested {n_genes} genes, but only {len(ALL_GENES)} are available")
    rng = _rng(seed)
    genes = rng.choice(ALL_GENES, size=n_genes, replace=False)
    direction = rng.choice([-3.0, 0.0, 3.0], size=n_genes, p=[0.25, 0.5, 0.25])
    weights = rng.normal(size=n_genes) * 2 + direction
    frame = pd.DataFrame({"gene": genes, "weight": weights})
    frame = frame.sort_values("weight", ascending=False).reset_index(drop=True)
    frame["rank"] = np.arange(1, len(frame) + 1)
    return frame


def generate_expression_matrix(
    n_genes: int = 2000,
    n_samples: int = 6,
    n_samples_group: int = 3,
    seed: int = 42,
) -> pd.DataFrame:
    """Generate a two-group synthetic expression matrix."""
    if n_samples != 2 * n_samples_group:
        raise ValueError("n_samples must equal two times n_samples_group")
    if n_genes > len(ALL_GENES):
        raise ValueError(f"Requested {n_genes} genes, but only {len(ALL_GENES)} are available")

    rng = _rng(seed)
    genes = rng.choice(ALL_GENES, size=n_genes, replace=False)
    sample_names = [f"Control_{index + 1}" for index in range(n_samples_group)] + [
        f"Treatment_{index + 1}" for index in range(n_samples_group)
    ]
    control = rng.normal(8, 2, size=(n_genes, n_samples_group))
    treatment = rng.normal(8, 2, size=(n_genes, n_samples_group))
    changed = rng.choice(n_genes, size=min(100, n_genes), replace=False)
    midpoint = len(changed) // 2
    treatment[changed[:midpoint]] += rng.uniform(1, 3, size=(midpoint, n_samples_group))
    treatment[changed[midpoint:]] -= rng.uniform(
        1, 3, size=(len(changed) - midpoint, n_samples_group)
    )
    return pd.DataFrame(
        np.hstack([control, treatment]), index=genes, columns=sample_names
    )


def generate_gene_sets(n_sets: int = 10, genes_per_set: int = 50) -> dict[str, set[str]]:
    """Generate deterministic synthetic gene sets around curated pathway seeds."""
    pathway_seeds = [
        ("Cell_Cycle", ["CCNA2", "CCNB1", "CCND1", "CCNE1", "CDC25A", "CHEK1", "AURKA", "TP53", "RB1"]),
        ("DNA_Repair", ["BRCA1", "BRCA2", "ATM", "ATR", "RAD51", "PALB2", "CHEK1", "PMS2", "MLH1"]),
        ("PI3K_AKT", ["AKT1", "AKT2", "PTEN", "PIK3CA", "MTOR"]),
        ("MAPK_Signaling", ["MAPK1", "MAPK3", "EGFR", "KRAS", "BRAF", "MAPK8", "JUN"]),
        ("Apoptosis", ["BCL2", "BAX", "BAK1", "CASP3", "CASP8", "CASP9", "PMAIP1", "BBC3", "MCL1"]),
        ("Immune_Response", ["PDCD1", "CD274", "CTLA4", "IFNG", "IL6", "TNF", "LAG3", "HAVCR2", "TIGIT"]),
        ("EMT", ["CDH1", "CDH2", "VIM", "SNAI1", "SNAI2", "ZEB1", "ZEB2", "TWIST1", "MMP2", "MMP9"]),
        ("Angiogenesis", ["VEGFA", "VEGFC", "FLT1", "KDR", "ANGPT1", "ANGPT2", "TEK", "PDGFA", "HIF1A", "EP300"]),
        ("Metabolism", ["HK2", "LDHA", "PKM", "PGK1", "ENO1", "GAPDH", "G6PD", "TPI1", "PFKL"]),
        ("Stem_Cell", ["NANOG", "POU5F1", "SOX2", "KLF4", "MYC", "LIN28A", "ALDH1A1", "PROM1", "CD44"]),
    ]
    rng = _rng(42)
    result: dict[str, set[str]] = {}
    for index in range(n_sets):
        name, core = pathway_seeds[index] if index < len(pathway_seeds) else (f"Pathway_{index + 1}", [])
        candidates = sorted(set(core) | set(rng.choice(ALL_GENES, size=100, replace=False)))
        size = min(genes_per_set, len(candidates))
        result[f"HSA_{name}_{index + 1}"] = set(rng.choice(candidates, size=size, replace=False))
    return result


def save_test_data(output_dir: str = "test_data") -> tuple[pd.DataFrame, pd.DataFrame, dict[str, set[str]]]:
    """Write the default ranked list, expression matrix, GMT, and query list."""
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)
    ranked = generate_sorted_gene_list()
    expression = generate_expression_matrix()
    gene_sets = generate_gene_sets()

    ranked_path = output_path / "ranked_genes.tsv"
    expression_path = output_path / "expression_matrix.tsv"
    gmt_path = output_path / "gene_sets.gmt"
    query_path = output_path / "top_degs.txt"
    ranked.to_csv(ranked_path, sep="\t", index=False)
    expression.to_csv(expression_path, sep="\t")
    with gmt_path.open("w", encoding="utf-8") as handle:
        for name, genes in gene_sets.items():
            gene_fields = "\t".join(sorted(genes))
            handle.write(f"{name}\tsynthetic fixture\t{gene_fields}\n")
    query_path.write_text("\n".join(ranked.head(100)["gene"]), encoding="utf-8")

    print(f"Ranked gene list saved: {ranked_path} ({len(ranked)} genes)")
    print(f"Expression matrix saved: {expression_path} ({expression.shape[0]} genes x {expression.shape[1]} samples)")
    print(f"Gene sets saved: {gmt_path} ({len(gene_sets)} sets)")
    print(f"Query list saved: {query_path} (100 genes)")
    return ranked, expression, gene_sets


def extract_gene_pool_from_gmt(gmt_dir: str) -> Set[str]:
    """Collect unique genes from the available human GMT files."""
    genes: set[str] = set()
    for filename in ("hsa.GO.gmt.gz", "hsa.KEGG.gmt.gz", "hsa.Reactome.gmt.gz", "hsa.DO.gmt.gz"):
        path = Path(gmt_dir) / filename
        if not path.exists():
            continue
        with gzip.open(path, "rt", encoding="utf-8") as handle:
            for line in handle:
                fields = line.rstrip("\n").split("\t")
                if len(fields) >= 3:
                    genes.update(fields[2:])
    return genes


def select_test_pathways_from_gmt(
    gmt_dir: str,
    pathways_per_db: int = 5,
    min_genes: int = 10,
    max_genes: int = 200,
) -> Dict[str, Set[str]]:
    """Select medium-sized gene sets from each available human database."""
    database_files = {
        "GO": "hsa.GO.gmt.gz",
        "KEGG": "hsa.KEGG.gmt.gz",
        "Reactome": "hsa.Reactome.gmt.gz",
        "DO": "hsa.DO.gmt.gz",
    }
    selected_sets: dict[str, set[str]] = {}
    for database, filename in database_files.items():
        path = Path(gmt_dir) / filename
        if not path.exists():
            continue
        candidates: list[tuple[str, set[str]]] = []
        with gzip.open(path, "rt", encoding="utf-8") as handle:
            for line in handle:
                fields = line.rstrip("\n").split("\t")
                genes = set(fields[2:]) if len(fields) >= 3 else set()
                if min_genes <= len(genes) <= max_genes:
                    candidates.append((fields[0], genes))
        candidates.sort(key=lambda item: abs(len(item[1]) - 50))
        for name, genes in candidates[:pathways_per_db]:
            selected_sets[f"{database}_{name}"] = genes
    return selected_sets


def generate_sorted_gene_list_from_pool(
    gene_pool: Set[str], n_genes: int = 500, seed: int = 42
) -> pd.DataFrame:
    """Generate a ranked list from a supplied gene pool."""
    available = sorted(gene_pool)
    count = min(n_genes, len(available))
    if count < n_genes:
        print(f"Warning: gene pool contains {count} genes; requested {n_genes}")
    rng = _rng(seed)
    selected = rng.choice(available, size=count, replace=False)
    direction = rng.choice([-2.0, 0.0, 2.0], size=count, p=[0.3, 0.4, 0.3])
    frame = pd.DataFrame({"gene": selected, "weight": rng.normal(size=count) * 1.5 + direction})
    frame = frame.sort_values("weight", ascending=False).reset_index(drop=True)
    frame["rank"] = np.arange(1, len(frame) + 1)
    return frame


def generate_expression_matrix_from_pool(
    gene_pool: Set[str], n_genes: int = 6000, n_samples: int = 6, seed: int = 42
) -> pd.DataFrame:
    """Generate a two-group expression matrix from a supplied gene pool."""
    if n_samples < 2 or n_samples % 2:
        raise ValueError("n_samples must be an even integer of at least two")
    available = sorted(gene_pool)
    count = min(n_genes, len(available))
    if count < n_genes:
        print(f"Warning: gene pool contains {count} genes; requested {n_genes}")
    rng = _rng(seed)
    selected = rng.choice(available, size=count, replace=False)
    per_group = n_samples // 2
    columns = [f"Control_{i + 1}" for i in range(per_group)] + [
        f"Treatment_{i + 1}" for i in range(per_group)
    ]
    control = rng.normal(8, 1.5, size=(count, per_group))
    treatment = rng.normal(8, 1.5, size=(count, per_group))
    changed = rng.choice(count, size=min(100, count), replace=False)
    midpoint = len(changed) // 2
    treatment[changed[:midpoint]] += rng.uniform(1, 2, size=(midpoint, per_group))
    treatment[changed[midpoint:]] -= rng.uniform(1, 2, size=(len(changed) - midpoint, per_group))
    return pd.DataFrame(np.hstack([control, treatment]), index=selected, columns=columns)


if __name__ == "__main__":
    print("=" * 60)
    print("GSEA, ssGSEA, and GSVA synthetic fixture generator")
    print("=" * 60)
    ranked_genes, expression_matrix, generated_sets = save_test_data()
    print(f"Ranked genes: {len(ranked_genes)}")
    print(f"Expression matrix: {expression_matrix.shape}")
    print(f"Gene sets: {len(generated_sets)}")
