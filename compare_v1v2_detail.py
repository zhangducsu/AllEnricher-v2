#!/usr/bin/env python3
"""
Details AllEnricher v1 and v2 Results of analysis
Deal with different forms of listing
"""

import pandas as pd
from pathlib import Path
import sys

def find_column(df, possible_names):
    """Find matching listings in DataFrame"""
    for name in possible_names:
        if name in df.columns:
            return name
    return None

def compare_database(db_name, v1_file, v2_file):
    """Comparison of analysis of individual databases"""
    print(f"\n{'='*70}")
    print(f"Database: {db_name}")
    print('='*70)
    
    # Read v1 result
    try:
        df_v1 = pd.read_csv(v1_file, sep='\t')
        v1_count = len(df_v1)
        print(f"v1 Entry: {v1_count}")
        print(f"v1 List: {list(df_v1.columns)}")
    except Exception as e:
        print(f"v1 file reading failed: {e}")
        return
    
    # Read v2 results
    try:
        df_v2 = pd.read_csv(v2_file, sep='\t')
        v2_count = len(df_v2)
        print(f"Number of entries: {v2_count}")
        print(f"v2 Listing: {list(df_v2.columns)}")
    except Exception as e:
        print(f"v2 File Reading Failed: {e}")
        return
    
    # Comparative entries
    diff = v2_count - v1_count
    if diff == 0:
        print(f"\n* Same number of entries: {v1_count}")
    else:
        print(f"\n⚠ Entry variance: {diff:+d} (v2 {'More' if diff > 0 else 'Less'} {abs(diff)} Article)")
    
    # Find Term Id Column
    term_col_v1 = find_column(df_v1, ['TermID', 'Term ID', 'term_id', 'ID'])
    term_col_v2 = find_column(df_v2, ['Term_ID', 'Term ID', 'term_id', 'ID'])
    
    if not term_col_v1:
        print(f"The  could not find the Term ID column in v1")
        return
    if not term_col_v2:
        print(f"The  cannot find the Term ID column in v2.")
        return
    
    print(f"\nTerm ID column: v1 ='{term_col_v1}', v2='{term_col_v2}'")
    
    # Compare Term ID
    v1_terms = set(df_v1[term_col_v1].astype(str))
    v2_terms = set(df_v2[term_col_v2].astype(str))
    
    common = v1_terms & v2_terms
    v1_only = v1_terms - v2_terms
    v2_only = v2_terms - v1_terms
    
    print(f"\nThe blog is a good example of how the country is a country where the world is not a land of its own")
    print(f"Common entry: {len(common)}")
    print(f"Only v1 has: {len(v1_only)}")
    print(f"Only v2 has: {len(v2_only)}")
    
    if len(v1_only) > 0:
        print(f"v1 unique example: {list(v1_only)[: 3]}")
    if len(v2_only) > 0:
        print(f"v2 is unique: {list(v2_only)[: 3]}")
    
    if len(common) == v1_count == v2_count:
        print("\nTerm IDs match exactly.")
    
    # Find P-value bar
    pval_col_v1 = find_column(df_v1, ['rawP', 'P-Value', 'p_value', 'pvalue', 'PValue'])
    pval_col_v2 = find_column(df_v2, ['P_Value', 'P-Value', 'p_value', 'pvalue', 'PValue'])
    
    if not pval_col_v1:
        print(f"\nCould not close temporary folder: %s")
        return
    if not pval_col_v2:
        print(f"\nP-value column not found in v2")
        return
    
    print(f"\nP column: v1={pval_col_v1}', v2='{pval_col_v2}'")
    
    # Compare P Value
    merged = df_v1[[term_col_v1, pval_col_v1]].merge(
        df_v2[[term_col_v2, pval_col_v2]], 
        left_on=term_col_v1,
        right_on=term_col_v2,
        how='inner'
    )
    
    if len(merged) > 0:
        # Calculate P value variance
        pval_v1 = merged[pval_col_v1].astype(float)
        pval_v2 = merged[pval_col_v2].astype(float)
        pval_diff = abs(pval_v1 - pval_v2)
        
        max_diff = pval_diff.max()
        mean_diff = pval_diff.mean()
        
        print(f"\nP-value comparison (common entry){len(merged)}(a) The number of persons:")
        print(f"Max.: {max_diff: .2e}")
        print(f"Average difference: {mean_diff: .2e}")
        
        if max_diff < 1e-10:
            print(f"  ✓ PValues are consistent (Variance < 1e-10)")
        elif max_diff < 1e-5:
            print(f"  ✓ PValues are broadly consistent. (Variance < 1e-5)")
        elif max_diff < 0.01:
            print(f"  ⚠ PThere's a slight difference in value. (Variance < 0.01)")
        else:
            print(f"  ⚠ PValues vary significantly (Variance >= 0.01)")
            # Show the largest differences
            merged['pval_diff'] = pval_diff
            top_diff = merged.nlargest(3, 'pval_diff')[[term_col_v1, pval_col_v1, pval_col_v2, 'pval_diff']]
            print(f"The most significant difference:")
            print(top_diff.to_string(index=False))
    
    # Compare Q/Correct P
    qval_col_v1 = find_column(df_v1, ['adjP', 'Q-Value', 'q_value', 'qvalue', 'Adjusted_P'])
    qval_col_v2 = find_column(df_v2, ['Adjusted_P_Value', 'Q-Value', 'q_value', 'qvalue', 'adjP'])
    
    if qval_col_v1 and qval_col_v2:
        print(f"\nQ-column: v1={qval_col_v1}', v2='{qval_col_v2}'")
        
        merged_q = df_v1[[term_col_v1, qval_col_v1]].merge(
            df_v2[[term_col_v2, qval_col_v2]], 
            left_on=term_col_v1,
            right_on=term_col_v2,
            how='inner'
        )
        
        if len(merged_q) > 0:
            qval_v1 = merged_q[qval_col_v1].astype(float)
            qval_v2 = merged_q[qval_col_v2].astype(float)
            qval_diff = abs(qval_v1 - qval_v2)
            
            max_diff_q = qval_diff.max()
            mean_diff_q = qval_diff.mean()
            
            print(f"\nQ Value Comparison (Common Entry){len(merged_q)}(a) The number of persons:")
            print(f"Max.: {max_diff_q: .2e}")
            print(f"Average difference: {mean_diff_q: .2e}")
            
            if max_diff_q < 1e-10:
                print(f"  ✓ QValues are consistent (Variance < 1e-10)")
            elif max_diff_q < 1e-5:
                print(f"  ✓ QValues are broadly consistent. (Variance < 1e-5)")
            elif max_diff_q < 0.01:
                print(f"  ⚠ QThere's a slight difference in value. (Variance < 0.01)")
            else:
                print(f"  ⚠ QValues vary significantly (Variance >= 0.01)")

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
    
    print("="*70)
    print("AllErricher v1 vs v2 Comparison of Detailed Results")
    print("="*70)
    
    for db_name, (v1_file, v2_file) in databases.items():
        compare_database(
            db_name,
            v1_dir / v1_file,
            v2_dir / v2_file
        )
    
    print("\n" + "="*70)
    print("Comparison completed")
    print("="*70)

if __name__ == '__main__':
    main()
