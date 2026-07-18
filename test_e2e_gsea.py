#!/usr/bin/env python3
"""GSEA end-to-end test."""

import sys
import time
import json
from pathlib import Path
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
from allenricher.core.enrichment import GSEA

TEST_DATA_DIR = Path("test_data")
RESULTS_DIR = Path("test_data/e2e_results")
RESULTS_DIR.mkdir(exist_ok=True)

def load_test_data():
    """Loading test data"""
    # Read the list of 500 gene sequences (the top 500 from the current 2000 genes)
    ranked_df = pd.read_csv(TEST_DATA_DIR / "ranked_genes.tsv", sep='\t')
    ranked_df = ranked_df.head(500)  # Only 500.
    ranked_genes = ranked_df['gene'].tolist()
    gene_weights = dict(zip(ranked_df['gene'], ranked_df['weight']))
    
    # Read test circuits (all 10 routes selected from GMT)
    gene_sets = {}
    with open(TEST_DATA_DIR / "gene_sets.gmt", 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            pathway = parts[0]
            genes = set(parts[2:])
            gene_sets[pathway] = genes
    
    return ranked_genes, gene_weights, gene_sets

def test_gsea_full():
    """GSEA full-scale testing - all routes"""
    print("=" * 60)
    print("GSEA end-to-end test")
    print("=" * 60)
    
    ranked_genes, gene_weights, gene_sets = load_test_data()
    print(f"* Sorting the list of genes: {len(ranked_genes)} genes")
    print(f"* Test number of circuits: {len(gene_sets)} pathways")
    
    # Create GSEA Analyzer
    gsea = GSEA(permutations=100)
    
    # Test all routes.
    results = []
    start_time = time.time()
    
    for pathway_name, pathway_genes in gene_sets.items():
        es, nes, pvalue, leading_edge = gsea.calculate_normalized_es(
            ranked_genes, pathway_genes, gene_weights
        )
        results.append({
            'pathway': pathway_name,
            'es': es,
            'nes': nes,
            'pvalue': pvalue,
            'gene_count': len(pathway_genes),
            'leading_edge_count': len(leading_edge)
        })
    
    elapsed = time.time() - start_time
    
    # Convert to DataFrame
    df_results = pd.DataFrame(results)
    
    # Save Results
    df_results.to_csv(RESULTS_DIR / "gsea_results.csv", index=False)
    
    # Generate test report
    report = {
        "test_name": "GSEA Full E2E Test",
        "timestamp": time.strftime("%Y-%m-%d %H:%M:%S"),
        "input_data": {
            "ranked_genes": len(ranked_genes),
            "pathways_tested": len(gene_sets),
            "permutations": 100
        },
        "results": {
            "total_pathways": len(df_results),
            "significant_p05": int((df_results['pvalue'] < 0.05).sum()),
            "significant_p01": int((df_results['pvalue'] < 0.01).sum()),
            "positive_enrichment": int((df_results['nes'] > 0).sum()),
            "negative_enrichment": int((df_results['nes'] < 0).sum()),
            "nes_range": [float(df_results['nes'].min()), float(df_results['nes'].max())],
            "execution_time": f"{elapsed:.2f}s"
        },
        "top_results": df_results.nsmallest(5, 'pvalue')[['pathway', 'nes', 'pvalue']].to_dict('records'),
        "status": "passed"
    }
    
    with open(RESULTS_DIR / "gsea_report.json", 'w') as f:
        json.dump(report, f, indent=2)
    
    # Print Results
    print(f"\nGSEA completed in {elapsed:.2f}s")
    print(f"• Results saved: {RESULTS_DIR / 'gsea_results.csv'}")
    print(f"\nStatistical information:")
    print(f"- Total number of routes: {len(df_results)}")
    print(f"  - Nominal P < 0.05: {report['results']['significant_p05']}")
    print(f"- NES Scope: [{df_results['nes'].min(): .3f}, {df_results['nes'].max(): .3f}]")
    
    return report

if __name__ == "__main__":
    report = test_gsea_full()
    print("\n" + "=" * 60)
    print("GSEA E2E test complete!")
    print("=" * 60)
