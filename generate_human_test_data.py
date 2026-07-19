#!/usr/bin/env python3
"""Generate human test data (from GMT files)"""

import sys
import json
from pathlib import Path
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
from test_data_generator import (
    extract_gene_pool_from_gmt,
    select_test_pathways_from_gmt,
    generate_sorted_gene_list_from_pool,
    generate_expression_matrix_from_pool
)

def main():
    print("=" * 60)
    print("Generate human test data (from GMT files)")
    print("=" * 60)

    gmt_dir = "database/organism/v20260515/hsa"
    output_dir = Path("test_data")
    output_dir.mkdir(exist_ok=True)

    # 1. Collect genes from the GMT file.
    print("\n[1/4] Collecting genes from the GMT file...")
    gene_pool = extract_gene_pool_from_gmt(gmt_dir)
    print(f"Gene pool size: {len(gene_pool):,} unique genes")

    # 2. Selection of test circuits
    print("\n[2/4Select test access from GMT...")
    test_pathways = select_test_pathways_from_gmt(gmt_dir, pathways_per_db=5)
    print(f"* Test number of routes: {len(test_pathways)}")
    for name, genes in list(test_pathways.items())[:5]:
        print(f"  - {name}: {len(genes)} genes")

    # 3. Generating 500 gene sequencing list
    print("\n[3/4_ Generate 500 gene ranked list...")
    n_sorted = min(500, len(gene_pool))
    ranked_genes = generate_sorted_gene_list_from_pool(gene_pool, n_genes=n_sorted)
    ranked_genes.to_csv(output_dir / "ranked_genes_500.tsv", sep='\t', index=False)
    print(f"*Save: ranked_genes_500.tsv ({len(ranked_genes)} genes)")

    # 4. Generation of 6, 000 x 6 matrix of expression
    print("\n [4/4] produces 6, 000 x 6 expression matrix...")
    n_expr = min(6000, len(gene_pool))
    expr_matrix = generate_expression_matrix_from_pool(gene_pool, n_genes=n_expr)
    expr_matrix.to_csv(output_dir / "expression_matrix_6000.tsv", sep='\t')
    print(f"*Save: conservation_matrix_6000.tsv ({expr_matrix.shape[0]} x {expr_matrix.shape[1]})")

    # 5. Preservation of test access
    with open(output_dir / "test_pathways_from_gmt.gmt", 'w') as f:
        for name, genes in test_pathways.items():
            genes_str = '\t'.join(genes)
            f.write(f"{name}\tfrom_gmt\t{genes_str}\n")
    print(f"*Save: test_pathways_frog_gmt.gmt ({len(test_pathways)} pathways)")

    # 6. Preservation of metadata
    metadata = {
        "timestamp": pd.Timestamp.now().isoformat(),
        "species": "hsa",
        "gene_pool_size": len(gene_pool),
        "test_data": {
            "sorted_genes": len(ranked_genes),
            "expression_matrix": f"{expr_matrix.shape[0]}×{expr_matrix.shape[1]}",
            "pathways": len(test_pathways)
        },
        "pathway_breakdown": {
            name: len(genes) for name, genes in test_pathways.items()
        }
    }

    with open(output_dir / "test_data_metadata.json", 'w') as f:
        json.dump(metadata, f, indent=2)
    print(f"*Save: test_data_metadata.json")

    print("\n" + "=" * 60)
    print("Test data generation complete!")
    print("=" * 60)

if __name__ == "__main__":
    main()
