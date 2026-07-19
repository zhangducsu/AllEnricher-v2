#!/usr/bin/env python3
"""
Compares allEnricher v1 and v2 analysis
"""

import pandas as pd
from pathlib import Path
import sys

def compare_database(db_name, v1_file, v2_file):
    """Comparison of analysis of individual databases"""
    print(f"\n{'='*60}")
    print(f"Database: {db_name}")
    print('='*60)
    
    # Read v1 result
    try:
        df_v1 = pd.read_csv(v1_file, sep='\t')
        v1_count = len(df_v1)
        print(f"v1 Entry: {v1_count}")
    except Exception as e:
        print(f"v1 file reading failed: {e}")
        v1_count = 0
        df_v1 = None
    
    # Read v2 results
    try:
        df_v2 = pd.read_csv(v2_file, sep='\t')
        v2_count = len(df_v2)
        print(f"Number of entries: {v2_count}")
    except Exception as e:
        print(f"v2 File Reading Failed: {e}")
        v2_count = 0
        df_v2 = None
    
    if df_v1 is None or df_v2 is None:
        return
    
    # Comparative entries
    diff = v2_count - v1_count
    if diff == 0:
        print(f"*Asymmetrical number of entries")
    else:
        print(f"⚠ Entry differences: {diff:+d} (v2 {'More' if diff > 0 else 'Less'} {abs(diff)} Article)")
    
    # Get Term ID List
    term_col = None
    for col in ['Term ID', 'term_id', 'ID', 'id']:
        if col in df_v1.columns:
            term_col = col
            break
    
    if term_col and term_col in df_v2.columns:
        v1_terms = set(df_v1[term_col].astype(str))
        v2_terms = set(df_v2[term_col].astype(str))
        
        common = v1_terms & v2_terms
        v1_only = v1_terms - v2_terms
        v2_only = v2_terms - v1_terms
        
        print(f"\nThe blog is a good example of how the country is a country where the world is not a land of its own")
        print(f"Common entry: {len(common)}")
        print(f"Only v1 has: {len(v1_only)}")
        print(f"Only v2 has: {len(v2_only)}")
        
        if len(v1_only) > 0:
            print(f"v1 unique example: {list(v1_only)[: 5]}")
        if len(v2_only) > 0:
            print(f"v2 is unique: {list(v2_only)[: 5]}")
    
    # Compare P values (for common entries)
    pval_col_v1 = None
    for col in ['P-Value', 'p_value', 'pvalue', 'PValue']:
        if col in df_v1.columns:
            pval_col_v1 = col
            break
    
    pval_col_v2 = None
    for col in ['P-Value', 'p_value', 'pvalue', 'PValue']:
        if col in df_v2.columns:
            pval_col_v2 = col
            break
    
    if pval_col_v1 and pval_col_v2 and term_col:
        # Merge comparison
        merged = df_v1[[term_col, pval_col_v1]].merge(
            df_v2[[term_col, pval_col_v2]], 
            on=term_col, 
            suffixes=('_v1', '_v2')
        )
        
        if len(merged) > 0:
            # Calculate P value variance
            pval_diff = abs(merged[pval_col_v1 + '_v1'] - merged[pval_col_v2 + '_v2'])
            max_diff = pval_diff.max()
            mean_diff = pval_diff.mean()
            
            print(f"\nP-value comparison (common entry){len(merged)}(a) The number of persons:")
            print(f"Max.: {max_diff: .2e}")
            print(f"Average difference: {mean_diff: .2e}")
            
            if max_diff < 1e-10:
                print(f"  ✓ PValue is equal (Variance < 1e-10)")
            elif max_diff < 1e-5:
                print(f"  ✓ PIt's basically the same value. (Variance < 1e-5)")
            else:
                print(f"Significant difference in P level")
                # Show the largest differences
                merged['diff'] = pval_diff
                top_diff = merged.nlargest(3, 'diff')[[term_col, pval_col_v1 + '_v1', pval_col_v2 + '_v2', 'diff']]
                print(f"The most significant difference:")
                print(top_diff.to_string(index=False))

def main():
    v1_dir = Path(r"F:\OneDrive\Documents\TraeSOLO\AllEnricher\AllEnricher-v1\example\allenricher\fisher\Q0.05")
    v2_dir = Path(r"F:\OneDrive\Documents\TraeSOLO\AllEnricher\AllEnricher-v2\test_output\v1v2_compare")
    
    # Database File Map
    databases = {
        'GO': ('example.glist.GO.xls', 'GO_enrichment.tsv'),
        'KEGG': ('example.glist.KEGG.xls', 'KEGG_enrichment.tsv'),
        'Reactome': ('example.glist.Reactome.xls', 'Reactome_enrichment.tsv'),
        'DO': ('example.glist.DO.xls', 'DO_enrichment.tsv'),
        'DisGeNET': ('example.glist.DisGeNET.xls', 'DisGeNET_enrichment.tsv'),
    }
    
    print("="*60)
    print("AllEnricher v1 vs v2")
    print("="*60)
    
    for db_name, (v1_file, v2_file) in databases.items():
        compare_database(
            db_name,
            v1_dir / v1_file,
            v2_dir / v2_file
        )
    
    print("\n" + "="*60)
    print("Comparison completed")
    print("="*60)

if __name__ == '__main__':
    main()
