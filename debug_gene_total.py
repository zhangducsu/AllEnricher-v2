#!/usr/bin/env python3
"""
验证v1的gene_total计算方式
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

# 加载输入基因
analyzer = EnrichmentAnalyzer(config)
gene_set = analyzer.load_gene_list(config.input_file)
print(f"输入基因总数: {len(gene_set)}")

# 获取GO数据库数据
go_data = db_manager.databases['GO']

# 计算匹配到任何条目的基因（v1的gene_total）
matched_genes = set()
for term_id, term_info in go_data.items():
    term_genes = set(term_info['genes'])
    genes_in_term = gene_set & term_genes
    matched_genes.update(genes_in_term)

v1_gene_total = len(matched_genes)
print(f"v1的gene_total (匹配到条目的基因): {v1_gene_total}")

# 用v1的方式计算期望值
term_id = 'GO:0051301'
term_info = go_data[term_id]
term_genes = set(term_info['genes'])

background_in_term = background_set & term_genes
genes_in_term = gene_set & term_genes

v1_expected = (len(background_in_term) / len(background_set)) * v1_gene_total
v2_expected = (len(background_in_term) / len(background_set)) * len(gene_set)

print(f"\n条目 {term_id} ({term_info['name']}):")
print(f"  背景基因中属于该条目: {len(background_in_term)}")
print(f"  输入基因中属于该条目: {len(genes_in_term)}")
print(f"  背景基因总数: {len(background_set)}")
print(f"\n  v1期望值: ({len(background_in_term)} / {len(background_set)}) * {v1_gene_total} = {v1_expected:.6f}")
print(f"  v2期望值: ({len(background_in_term)} / {len(background_set)}) * {len(gene_set)} = {v2_expected:.6f}")
print(f"  v1实际值: 17.347928")
