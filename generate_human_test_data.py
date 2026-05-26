#!/usr/bin/env python3
"""生成人类测试数据（从GMT文件）"""

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
    print("生成人类测试数据（从GMT文件）")
    print("=" * 60)

    gmt_dir = "database/organism/v20260515/hsa"
    output_dir = Path("test_data")
    output_dir.mkdir(exist_ok=True)

    # 1. 提取基因池
    print("\n[1/4] 从GMT文件提取基因池...")
    gene_pool = extract_gene_pool_from_gmt(gmt_dir)
    print(f"✓ 基因池大小: {len(gene_pool):,} 个唯一基因")

    # 2. 选取测试通路
    print("\n[2/4] 从GMT选取测试通路...")
    test_pathways = select_test_pathways_from_gmt(gmt_dir, pathways_per_db=5)
    print(f"✓ 测试通路数: {len(test_pathways)}")
    for name, genes in list(test_pathways.items())[:5]:
        print(f"  - {name}: {len(genes)} genes")

    # 3. 生成500基因排序列表
    print("\n[3/4] 生成500基因排序列表...")
    n_sorted = min(500, len(gene_pool))
    ranked_genes = generate_sorted_gene_list_from_pool(gene_pool, n_genes=n_sorted)
    ranked_genes.to_csv(output_dir / "ranked_genes_500.tsv", sep='\t', index=False)
    print(f"✓ 保存: ranked_genes_500.tsv ({len(ranked_genes)} genes)")

    # 4. 生成6000×6表达矩阵
    print("\n[4/4] 生成6000×6表达矩阵...")
    n_expr = min(6000, len(gene_pool))
    expr_matrix = generate_expression_matrix_from_pool(gene_pool, n_genes=n_expr)
    expr_matrix.to_csv(output_dir / "expression_matrix_6000.tsv", sep='\t')
    print(f"✓ 保存: expression_matrix_6000.tsv ({expr_matrix.shape[0]}×{expr_matrix.shape[1]})")

    # 5. 保存测试通路
    with open(output_dir / "test_pathways_from_gmt.gmt", 'w') as f:
        for name, genes in test_pathways.items():
            genes_str = '\t'.join(genes)
            f.write(f"{name}\tfrom_gmt\t{genes_str}\n")
    print(f"✓ 保存: test_pathways_from_gmt.gmt ({len(test_pathways)} pathways)")

    # 6. 保存元数据
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
    print(f"✓ 保存: test_data_metadata.json")

    print("\n" + "=" * 60)
    print("测试数据生成完成!")
    print("=" * 60)

if __name__ == "__main__":
    main()
