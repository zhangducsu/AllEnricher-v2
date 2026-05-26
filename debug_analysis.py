#!/usr/bin/env python3
"""
调试v2的分析过程
"""

import sys
sys.path.insert(0, r'F:\OneDrive\Documents\TraeSOLO\AllEnricher\AllEnricher-v2')

from allenricher.database.manager import DatabaseManager
from allenricher.core.enrichment import EnrichmentAnalyzer
from allenricher.core.config import Config

# 创建配置
config = Config(
    input_file=r'F:\OneDrive\Documents\TraeSOLO\AllEnricher\AllEnricher-v1\example\allenricher\example.glist',
    output_dir=r'F:\OneDrive\Documents\TraeSOLO\AllEnricher\AllEnricher-v2\test_output\debug',
    species='hsa',
    databases=['GO'],
    method='fisher',
    correction='BH',
    pvalue_cutoff=0.05,
    qvalue_cutoff=0.05,
    min_genes=2,
)

# 加载数据库
db_manager = DatabaseManager(
    r'F:\OneDrive\Documents\TraeSOLO\AllEnricher\AllEnricher-v1\database\organism\v20190612\hsa',
    'hsa'
)
db_manager.load_databases(['GO'])

# 获取背景基因
background_set = db_manager.get_background_genes()
print(f"v2背景基因总数: {len(background_set)}")

# 加载输入基因
analyzer = EnrichmentAnalyzer(config)
gene_set = analyzer.load_gene_list(config.input_file)
print(f"输入基因数: {len(gene_set)}")

# 获取GO数据库数据
go_data = db_manager.databases['GO']
print(f"GO条目数: {len(go_data)}")

# 检查一个具体条目
term_id = 'GO:0051301'  # cell division
if term_id in go_data:
    term_info = go_data[term_id]
    term_genes = set(term_info['genes'])
    print(f"\n条目 {term_id} ({term_info['name']}):")
    print(f"  条目基因数: {len(term_genes)}")
    
    # 计算交集
    genes_in_term = gene_set & term_genes
    background_in_term = background_set & term_genes
    print(f"  输入基因中属于该条目: {len(genes_in_term)}")
    print(f"  背景基因中属于该条目: {len(background_in_term)}")
    
    # 计算期望值
    gene_total = len(gene_set)
    background_total = len(background_set)
    expected = (len(background_in_term) / background_total) * gene_total
    print(f"  期望值计算: ({len(background_in_term)} / {background_total}) * {gene_total} = {expected:.6f}")
    
    # v1的期望值
    v1_expected = 17.347928
    print(f"  v1的期望值: {v1_expected}")
    print(f"  差异: {abs(expected - v1_expected):.6f}")
