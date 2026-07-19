#!/usr/bin/env python3
"""
Accurate simulation of the calculation logic of v1 - fixes the gene_tal calculation
"""

import gzip
import sys
sys.path.insert(0, r'F:\OneDrive\Documents\TraeSOLO\AllEnricher\AllEnricher-v2')

from allenricher.database.manager import DatabaseManager
from allenricher.core.enrichment import EnrichmentAnalyzer
from allenricher.core.config import Config

# Read gene_info as Background_list
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

print(f"Gene count: {len(background_list)}")

# Loading Database
db_manager = DatabaseManager(
    r'F:\OneDrive\Documents\TraeSOLO\AllEnricher\AllEnricher-v1\database\organism\v20190612\hsa',
    'hsa'
)
db_manager.load_databases(['GO'])

# Get GO database data
go_data = db_manager.databases['GO']

# Statistical Background_ttal (only gene in GO)
go_genes = set()
for term_id, term_info in go_data.items():
    for gene in term_info['genes']:
        go_genes.add(gene)

background_total = len(go_genes)
print(f"The genetics in the Go file: {background_total}")

# Loading input gene
config = Config(
    input_file=r'F:\OneDrive\Documents\TraeSOLO\AllEnricher\AllEnricher-v1\example\allenricher\example.glist',
    output_dir=r'F:\OneDrive\Documents\TraeSOLO\AllEnricher\AllEnricher-v2\test_output\debug',
    species='hsa',
    databases=['GO'],
)
analyzer = EnrichmentAnalyzer(config)
gene_set = analyzer.load_gene_list(config.input_file)
gene_list = {g: 1 for g in gene_set}  # Simulate v1%gene_list

# The complete analytical logic of simulation v1 - calculate the matching genes for all entries first
gene_list1 = {}  # Simulate v1%géne_list1 - To reassemble in all entries matching genes

for term_id, term_info in go_data.items():
    term_genes = {g: 1 for g in term_info['genes']}
    for gene in term_genes:
        if gene in gene_list:
            gene_list1[gene] = 1

gene_total = len(gene_list1)
print(f"Gene_ttal (de-weighting of all entries matching genes): {gene_total}")

# Calculating expectations for specific entries
term_id = 'GO:0051301'
term_info = go_data[term_id]
term_genes = {g: 1 for g in term_info['genes']}

num_in_C = 0
num_in_O = 0

for gene in term_genes:
    if gene in background_list:
        num_in_C += 1
    if gene in gene_list:
        num_in_O += 1

expected = num_in_C / background_total * gene_total if background_total > 0 else 0

print(f"\nEntry{term_id} ({term_info['name']}):")
print(f"Num_in_C (in blackground_list for entry): {num_in_C}")
print(f"Num_in_O (entry gene in the gene list): {num_in_O}")
print(f"  gene_total: {gene_total}")
print(f"  background_total: {background_total}")
print(f"\nCalculate expectations: ({num_in_C} / {background_total}) * {gene_total} = {expected: .6f}")
print(f"V1 Actual expected value: 17.347928")
print(f"Variance: {abs(expected - 17.347928): .6f}")
