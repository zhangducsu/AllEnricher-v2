#!/usr/bin/env python3
"""
检查v2的背景基因计算
"""

import gzip
from pathlib import Path
import sys

# 检查GO数据库文件
go_file = Path(r"F:\OneDrive\Documents\TraeSOLO\AllEnricher\AllEnricher-v1\database\organism\v20190612\hsa\hsa.GO2gene.tab.gz")

# 统计GO文件中的基因数
with gzip.open(go_file, 'rt') as f:
    header = f.readline()
    genes_in_go = sum(1 for _ in f)

print(f"GO数据库文件中的基因数: {genes_in_go}")
print(f"表头中的条目数: {len(header.strip().split(chr(9))) - 1}")

# 读取gene_info中的基因
gene_info_file = Path(r"F:\OneDrive\Documents\TraeSOLO\AllEnricher\AllEnricher-v1\database\organism\v20190612\hsa\hsa.gene_info")
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

print(f"gene_info中的基因数: {len(background_list)}")

# 检查GO文件中的基因有多少在background_list中
with gzip.open(go_file, 'rt') as f:
    header = f.readline()
    genes_in_background = 0
    for line in f:
        gene = line.strip().split('\t')[0]
        if gene in background_list:
            genes_in_background += 1

print(f"GO文件中同时在background_list中的基因数: {genes_in_background}")
