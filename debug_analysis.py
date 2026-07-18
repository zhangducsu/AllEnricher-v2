#!/usr/bin/env python3
"""Inspect one GO term while debugging the v1-to-v2 ORA comparison."""

import sys
sys.path.insert(0, r'F:\OneDrive\Documents\TraeSOLO\AllEnricher\AllEnricher-v2')

from allenricher.database.manager import DatabaseManager
from allenricher.core.enrichment import EnrichmentAnalyzer
from allenricher.core.config import Config

# Configure the historical comparison fixture.
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

# Load the v1 database snapshot used by this comparison.
db_manager = DatabaseManager(
    r'F:\OneDrive\Documents\TraeSOLO\AllEnricher\AllEnricher-v1\database\organism\v20190612\hsa',
    'hsa'
)
db_manager.load_databases(['GO'])

# Retrieve the annotated background.
background_set = db_manager.get_background_genes()
print(f"Total background genes: {len(background_set)}")

# Load the query genes.
analyzer = EnrichmentAnalyzer(config)
gene_set = analyzer.load_gene_list(config.input_file)
print(f"Query genes: {len(gene_set)}")

# Inspect the GO database.
go_data = db_manager.databases['GO']
print(f"GO terms: {len(go_data)}")

# Check one term used in the historical comparison.
term_id = 'GO:0051301'  # cell division
if term_id in go_data:
    term_info = go_data[term_id]
    term_genes = set(term_info['genes'])
    print(f"\nTerm {term_id} ({term_info['name']}):")
    print(f"Term genes: {len(term_genes)}")
    
    # Calculate query and background intersections.
    genes_in_term = gene_set & term_genes
    background_in_term = background_set & term_genes
    print(f"Query genes in term: {len(genes_in_term)}")
    print(f"Background genes in term: {len(background_in_term)}")
    
    # Calculate the expected overlap.
    gene_total = len(gene_set)
    background_total = len(background_set)
    expected = (len(background_in_term) / background_total) * gene_total
    print(
        f"Expected overlap: ({len(background_in_term)} / {background_total}) "
        f"* {gene_total} = {expected:.6f}"
    )
    
    # Historical v1 reference value.
    v1_expected = 17.347928
    print(f"V1 expected overlap: {v1_expected}")
    print(f"Absolute difference: {abs(expected - v1_expected):.6f}")
