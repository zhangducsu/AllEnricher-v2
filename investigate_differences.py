#!/usr/bin/env python3
"""
Scripts for differences in results of investigations v1 and v2
"""

import sys
from pathlib import Path
import pandas as pd
import numpy as np

# Add Item Path
project_root = Path(__file__).parent
sys.path.insert(0, str(project_root))

from allenricher.core.config import Config
from allenricher.core.enrichment import EnrichmentAnalyzer
from allenricher.database.manager import DatabaseManager


def load_genes():
    """Loading list of genes"""
    gene_list_file = Path(__file__).parent / "example_genes.txt"
    with open(gene_list_file, 'r') as f:
        genes = [line.strip() for line in f if line.strip()]
    return genes


def check_background_genes():
    """Checking the number of background genes"""
    print("=" * 80)
    print("Check background genes.")
    print("=" * 80)
    
    v1_db_dir = Path(__file__).parent.parent / "AllEnricher-v1" / "database" / "organism" / "v20190612" / "hsa"
    
    db_manager = DatabaseManager(
        database_dir=str(v1_db_dir),
        species="hsa"
    )
    db_manager.load_databases(["GO", "KEGG"])
    
    bg_genes = db_manager.get_background_genes()
    print(f"Number of background genes loaded in v2: {len(bg_genes)}")
    print()
    
    # Checking the total number of genes in each database
    db_data = db_manager.get_all_term_data()
    for db_name, data in db_data.items():
        all_genes = set()
        for term_id, genes in data.items():
            all_genes.update(genes)
        print(f"{db_name}Total: {len(all_genes)}Genome.")
    
    print()
    return bg_genes, db_data


def check_single_term_calculation(db_data, gene_set, bg_genes):
    """Check the calculation of the individual term"""
    print("=" * 80)
    print("Check for individual term enrichment analysis calculations")
    print("=" * 80)
    
    # Check KEG for hsa00010
    print("\n---Check KEG hsa00010---")
    if 'KEGG' in db_data and 'hsa00010' in db_data['KEGG']:
        term_genes = db_data['KEGG']['hsa00010']
        print(f"Genome count: {len(term_genes)}")
        
        # Calculate 2x2 Columns
        n_obs = len(gene_set & term_genes)
        n_term = len(term_genes)
        n_bg = len(bg_genes)
        n_gene = len(gene_set)
        
        a = n_obs
        b = n_gene - a
        c = n_term - a
        d = n_bg - n_gene - c
        
        print(f"2x2 Table: a={a}, b={b}, c={c}, d={d}")
        
        # Manually calculates the p value
        from scipy.stats import fisher_exact
        oddsratio, pvalue = fisher_exact([[a, b], [c, d]], alternative='greater')
        print(f"Manually calculate p: {pvalue}")
        
        # Calculating with v2 code
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
        print(f"v2 Code calculation p value: {result['pvalue']}")
    
    # Check another GO term
    print("\n---Check GO: 0005576---")
    if 'GO' in db_data and 'GO:0005576' in db_data['GO']:
        term_genes = db_data['GO']['GO:0005576']
        print(f"Genome count: {len(term_genes)}")
        
        n_obs = len(gene_set & term_genes)
        n_term = len(term_genes)
        n_bg = len(bg_genes)
        n_gene = len(gene_set)
        
        a = n_obs
        b = n_gene - a
        c = n_term - a
        d = n_bg - n_gene - c
        
        print(f"2x2 Table: a={a}, b={b}, c={c}, d={d}")
        
        from scipy.stats import fisher_exact
        oddsratio, pvalue = fisher_exact([[a, b], [c, d]], alternative='greater')
        print(f"Manually calculate p: {pvalue}")
        
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
        print(f"v2 Code calculation p value: {result['pvalue']}")


def main():
    genes = load_genes()
    print(f"Entering number of genes: {len(genes)}")
    
    bg_genes, db_data = check_background_genes()
    
    check_single_term_calculation(db_data, set(genes), bg_genes)


if __name__ == "__main__":
    main()
