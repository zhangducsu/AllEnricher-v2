#!/usr/bin/env python3
"""
Check the background genetic calculations of v2.
"""

import gzip
from pathlib import Path
import sys

# Check Go database files
go_file = Path(r"F:\OneDrive\Documents\TraeSOLO\AllEnricher\AllEnricher-v1\database\organism\v20190612\hsa\hsa.GO2gene.tab.gz")

# Count genes in Go files
with gzip.open(go_file, 'rt') as f:
    header = f.readline()
    genes_in_go = sum(1 for _ in f)

print(f"Gene count in GO database file: {genes_in_go}")
print(f"Number of entries in table header: {len(header.strip().split(chr(9))) - 1}")

# Read genes in gene_info
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

print(f"Genes in gene_info: {len(background_list)}")

# Checking how many genes in Go files are in the Background_list
with gzip.open(go_file, 'rt') as f:
    header = f.readline()
    genes_in_background = 0
    for line in f:
        gene = line.strip().split('\t')[0]
        if gene in background_list:
            genes_in_background += 1

print(f"The number of genes in the GO document is also in the Background_list: {genes_in_background}")
