#!/usr/bin/env python3
"""
精确模拟v1的计算逻辑
"""

import gzip
import sys
sys.path.insert(0, r'F:\OneDrive\Documents\TraeSOLO\AllEnricher\AllEnricher-v2')

from allenricher.database.manager import DatabaseManager
from allenricher.core.enrichment import EnrichmentAnalyzer
from allenricher.core.config import Config

# 读取gene_info作为background_list
gene_info_file = r'F:\OneDrive\Documents\TraeSOLO\AllEnricher\AllEnricher-v1\database\organism\v20190612\hsa\hsa.gene_info'
background_list = set()
with open(gene_info_file, 'r') as f:
    for line in f:
        if line.startswith('#'):
            continue
        parts = line.strip().split('\t')
        if len(parts) >= 3:
            gene_symbol = parts[2]
            if gene_symbol:
                background_list.add(gene_symbol)

print(f"background_list (gene_info) 基因数: {len(background_list)}")

# 加载数据库
db_manager = DatabaseManager(
    r'F:\OneDrive\Documents\TraeSOLO\AllEnricher\AllEnricher-v1\database\organism\v20190612\hsa',
    'hsa'
)
db_manager.load_databases(['GO'])

# 获取GO数据库数据
go_data = db_manager.databases['GO']

# 统计background_total（从GO文件读取的基因数，但只统计在background_list中的）
background_total = 0
for term_id, term_info in go_data.items():
    for gene in term_info['genes']:
        if gene in background_list:
            background_total += 1
            break  # 每个基因只统计一次

# 实际上background_total应该是GO文件中的唯一基因数
go_genes = set()
for term_id, term_info in go_data.items():
    for gene in term_info['genes']:
        go_genes.add(gene)

background_total = len(go_genes)
print(f"background_total (GO文件中的基因数): {background_total}")

# 加载输入基因
config = Config(
    input_file=r'F:\OneDrive\Documents\TraeSOLO\AllEnricher\AllEnricher-v1\example\allenricher\example.glist',
    output_dir=r'F:\OneDrive\Documents\TraeSOLO\AllEnricher\AllEnricher-v2\test_output\debug',
    species='hsa',
    databases=['GO'],
)
analyzer = EnrichmentAnalyzer(config)
gene_set = analyzer.load_gene_list(config.input_file)
gene_list = {g: 1 for g in gene_set}  # 模拟v1的%gene_list

# 模拟v1的分析逻辑
term_id = 'GO:0051301'
term_info = go_data[term_id]
term_genes = {g: 1 for g in term_info['genes']}  # 模拟v1的$ref{$term}

num_in_C = 0
num_in_O = 0
gene_list1 = {}  # 模拟v1的%gene_list1

for gene in term_genes:
    if gene in background_list:
        num_in_C += 1
    if gene in gene_list:
        num_in_O += 1
        gene_list1[gene] = 1

gene_total = len(gene_list1)
expected = num_in_C / background_total * gene_total if background_total > 0 else 0

print(f"\n条目 {term_id} ({term_info['name']}):")
print(f"  num_in_C (在background_list中的条目基因): {num_in_C}")
print(f"  num_in_O (在输入基因列表中的条目基因): {num_in_O}")
print(f"  gene_total (匹配的基因数): {gene_total}")
print(f"  background_total: {background_total}")
print(f"\n  计算期望值: ({num_in_C} / {background_total}) * {gene_total} = {expected:.6f}")
print(f"  v1实际期望值: 17.347928")
print(f"  差异: {abs(expected - 17.347928):.6f}")
