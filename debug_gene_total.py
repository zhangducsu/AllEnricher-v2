#!/usr/bin/env python3
"""
Verify the gene_ttal calculation of v1
"""

import sys
sys.path.insert(0, r'F:\OneDrive\Documents\TraeSOLO\AllEnricher\AllEnricher-v2')

from allenricher.database.manager import DatabaseManager
from allenricher.core.enrichment import EnrichmentAnalyzer
from allenricher.core.config import Config

# Create Configuration
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

# Loading Database
db_manager = DatabaseManager(
    r'F:\OneDrive\Documents\TraeSOLO\AllEnricher\AllEnricher-v1\database\organism\v20190612\hsa',
    'hsa'
)
db_manager.load_databases(['GO'])

# Retrieving background genes
background_set = db_manager.get_background_genes()

# Loading input gene
analyzer = EnrichmentAnalyzer(config)
gene_set = analyzer.load_gene_list(config.input_file)
print(f"Total number of genes entered: {len(gene_set)}")

# Get GO database data
go_data = db_manager.databases['GO']

# Calculate genes matching any entry (gene_total of v1)
matched_genes = set()
for term_id, term_info in go_data.items():
    term_genes = set(term_info['genes'])
    genes_in_term = gene_set & term_genes
    matched_genes.update(genes_in_term)

v1_gene_total = len(matched_genes)
print(f"GENE_TOTAL OF V1 (genesis matching to entry): {v1_gene_total}")

# The calculation of the expected value in v1
term_id = 'GO:0051301'
term_info = go_data[term_id]
term_genes = set(term_info['genes'])

background_in_term = background_set & term_genes
genes_in_term = gene_set & term_genes

v1_expected = (len(background_in_term) / len(background_set)) * v1_gene_total
v2_expected = (len(background_in_term) / len(background_set)) * len(gene_set)

print(f"\nEntry{term_id} ({term_info['name']}):")
print(f"The background gene belongs to this entry: {len(background_in_term)}")
print(f"Enter the entry in the gene that belongs to: {len(genes_in_term)}")
print(f"Total number of background genes: {len(background_set)}")
print(f"\nv1 Expected value: (b){len(background_in_term)} / {len(background_set)}) * {v1_gene_total}  = {v1_expected:. 6f}")
print(f"v2Expectations: ({len(background_in_term)} / {len(background_set)}) * {len(gene_set)} = { v2_expected:. 6f}")
print(f"v1 Actual value: 17.347928")
