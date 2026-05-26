#!/usr/bin/env python3
"""
调查 v1 和 v2 结果差异的脚本
"""

import sys
from pathlib import Path
import pandas as pd
import numpy as np

# 添加项目路径
project_root = Path(__file__).parent
sys.path.insert(0, str(project_root))

from allenricher.core.config import Config
from allenricher.core.enrichment import EnrichmentAnalyzer
from allenricher.database.manager import DatabaseManager


def load_genes():
    """加载基因列表"""
    gene_list_file = Path(__file__).parent / "example_genes.txt"
    with open(gene_list_file, 'r') as f:
        genes = [line.strip() for line in f if line.strip()]
    return genes


def check_background_genes():
    """检查背景基因数量"""
    print("=" * 80)
    print("检查背景基因")
    print("=" * 80)
    
    v1_db_dir = Path(__file__).parent.parent / "AllEnricher-v1" / "database" / "organism" / "v20190612" / "hsa"
    
    db_manager = DatabaseManager(
        database_dir=str(v1_db_dir),
        species="hsa"
    )
    db_manager.load_databases(["GO", "KEGG"])
    
    bg_genes = db_manager.get_background_genes()
    print(f"v2 加载的背景基因数: {len(bg_genes)}")
    print()
    
    # 检查各个数据库的基因总数
    db_data = db_manager.get_all_term_data()
    for db_name, data in db_data.items():
        all_genes = set()
        for term_id, genes in data.items():
            all_genes.update(genes)
        print(f"{db_name}: 共 {len(all_genes)} 个基因")
    
    print()
    return bg_genes, db_data


def check_single_term_calculation(db_data, gene_set, bg_genes):
    """检查单个 term 的计算"""
    print("=" * 80)
    print("检查单个 term 的富集分析计算")
    print("=" * 80)
    
    # 先检查 KEGG 的 hsa00010
    print("\n--- 检查 KEGG hsa00010 ---")
    if 'KEGG' in db_data and 'hsa00010' in db_data['KEGG']:
        term_genes = db_data['KEGG']['hsa00010']
        print(f"term 基因数: {len(term_genes)}")
        
        # 计算 2x2 列联表
        n_obs = len(gene_set & term_genes)
        n_term = len(term_genes)
        n_bg = len(bg_genes)
        n_gene = len(gene_set)
        
        a = n_obs
        b = n_gene - a
        c = n_term - a
        d = n_bg - n_gene - c
        
        print(f"2x2 表: a={a}, b={b}, c={c}, d={d}")
        
        # 手动计算 p 值
        from scipy.stats import fisher_exact
        oddsratio, pvalue = fisher_exact([[a, b], [c, d]], alternative='greater')
        print(f"手动计算 p 值: {pvalue}")
        
        # 使用 v2 代码计算
        config = Config(
            species="hsa",
            databases=["KEGG"],
            method="fisher",
            correction="BH",
            pvalue_cutoff=1.0,
            qvalue_cutoff=1.0
        )
        analyzer = EnrichmentAnalyzer(config)
        result = analyzer._analyze_single_term(
            term_genes=term_genes,
            gene_set=gene_set,
            background_set=bg_genes
        )
        print(f"v2 代码计算 p 值: {result['pvalue']}")
    
    # 再检查一个 GO term
    print("\n--- 检查 GO:0005576 ---")
    if 'GO' in db_data and 'GO:0005576' in db_data['GO']:
        term_genes = db_data['GO']['GO:0005576']
        print(f"term 基因数: {len(term_genes)}")
        
        n_obs = len(gene_set & term_genes)
        n_term = len(term_genes)
        n_bg = len(bg_genes)
        n_gene = len(gene_set)
        
        a = n_obs
        b = n_gene - a
        c = n_term - a
        d = n_bg - n_gene - c
        
        print(f"2x2 表: a={a}, b={b}, c={c}, d={d}")
        
        from scipy.stats import fisher_exact
        oddsratio, pvalue = fisher_exact([[a, b], [c, d]], alternative='greater')
        print(f"手动计算 p 值: {pvalue}")
        
        config = Config(
            species="hsa",
            databases=["GO"],
            method="fisher",
            correction="BH"
        )
        analyzer = EnrichmentAnalyzer(config)
        result = analyzer._analyze_single_term(
            term_genes=term_genes,
            gene_set=gene_set,
            background_set=bg_genes
        )
        print(f"v2 代码计算 p 值: {result['pvalue']}")


def main():
    genes = load_genes()
    print(f"输入基因数: {len(genes)}")
    
    bg_genes, db_data = check_background_genes()
    
    check_single_term_calculation(db_data, set(genes), bg_genes)


if __name__ == "__main__":
    main()
