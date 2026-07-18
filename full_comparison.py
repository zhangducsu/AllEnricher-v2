#!/usr/bin/env python3
"""
Script for full comparison v1 and v2 results
"""

import sys
from pathlib import Path
import pandas as pd
import numpy as np


def load_and_normalize_v1():
    """Load and standardize v1 results"""
    v1_results_dir = Path(__file__).parent.parent / "AllEnricher-v1" / "example" / "allenricher" / "fisher" / "Q0.05"
    
    db_files = {
        'GO': 'example.glist.GO.xls',
        'KEGG': 'example.glist.KEGG.xls',
        'Reactome': 'example.glist.Reactome.xls',
        'DO': 'example.glist.DO.xls'
    }
    
    normalized_results = {}
    
    for db_name, filename in db_files.items():
        filepath = v1_results_dir / filename
        if not filepath.exists():
            continue
        
        df = pd.read_csv(filepath, sep='\t')
        df['TermID'] = df['TermID'].astype(str)
        
        # Standardized listing
        df_norm = pd.DataFrame()
        df_norm['Term_ID'] = df['TermID']
        df_norm['Term_Name'] = df['TermName']
        df_norm['P_Value'] = df['rawP']
        df_norm['Adjusted_P_Value'] = df['adjP']
        df_norm['Gene_Count'] = df['ObservedGeneNum']
        df_norm['Background_Count'] = df['TermGeneNum']
        
        normalized_results[db_name] = df_norm
    
    return normalized_results


def load_v2_results():
    """Loading v2 Results"""
    v2_dir = Path(__file__).parent / "comparison_output"
    
    results = {}
    for db_name in ['GO', 'KEGG', 'Reactome', 'DO']:
        filepath = v2_dir / f"v2_{db_name}.tsv"
        if filepath.exists():
            df = pd.read_csv(filepath, sep='\t')
            df['Term_ID'] = df['Term_ID'].astype(str)
            results[db_name] = df
    
    return results


def save_v1_normalized(v1_norm):
    """Keep standardized v1 results"""
    output_dir = Path(__file__).parent / "comparison_output"
    for db_name, df in v1_norm.items():
        out_file = output_dir / f"v1_{db_name}_normalized.tsv"
        df.to_csv(out_file, sep='\t', index=False)


def compare_databases(v1_norm, v2_results):
    """Compare all database results"""
    print("=" * 80)
    print("AllEnricher v1 vs v2 is a detailed comparison")
    print("=" * 80)
    print()
    
    all_match = True
    
    for db_name in ['GO', 'KEGG', 'Reactome', 'DO']:
        if db_name not in v1_norm or db_name not in v2_results:
            continue
        
        print(f"[ {db_name} ]")
        print("-" * 80)
        
        v1_df = v1_norm[db_name]
        v2_df = v2_results[db_name]
        
        print(f"v1 Entry: {len(v1_df)}")
        print(f"v2 Entry: {len(v2_df)}")
        
        # Press Term_ID Alignment
        v1_sorted = v1_df.sort_values('Term_ID').set_index('Term_ID')
        v2_sorted = v2_df.sort_values('Term_ID').set_index('Term_ID')
        
        # Common entry
        common_terms = v1_sorted.index.intersection(v2_sorted.index)
        v1_only = v1_sorted.index.difference(v2_sorted.index)
        v2_only = v2_sorted.index.difference(v1_sorted.index)
        
        print(f"Common entry: {len(common_terms)}")
        if len(v1_only) > 0:
            print(f"v1 Special: {len(v1_only)}")
        if len(v2_only) > 0:
            print(f"v2 Special: {len(v2_only)}")
        
        if len(common_terms) == 0:
            print()
            continue
        
        # Common entries only
        v1_common = v1_sorted.loc[common_terms]
        v2_common = v2_sorted.loc[common_terms]
        
        # Compare rows
        columns_to_compare = ['P_Value', 'Adjusted_P_Value', 'Gene_Count', 'Background_Count']
        db_match = True
        
        for col in columns_to_compare:
            if col not in v1_common.columns or col not in v2_common.columns:
                continue
            
            v1_vals = v1_common[col]
            v2_vals = v2_common[col]
            
            if col in ['P_Value', 'Adjusted_P_Value']:
                # Float Comparison
                if np.allclose(v1_vals, v2_vals, rtol=1e-5, atol=1e-8):
                    print(f"{col}: Perfectly matched *")
                else:
                    db_match = False
                    all_match = False
                    # Calculated variance
                    diff_mask = ~np.isclose(v1_vals, v2_vals, rtol=1e-5, atol=1e-8)
                    diff_count = np.sum(diff_mask)
                    print(f"{col}: {diff_count}Variance")
                    
                    # Show the first 3 differences
                    diff_terms = v1_common.index[diff_mask][:3]
                    for term in diff_terms:
                        print(f"    {term}: v1={v1_common.loc[term, col]:.6g}, v2={v2_common.loc[term, col]:.6g}")
            else:
                # Integer comparison
                if (v1_vals == v2_vals).all():
                    print(f"{col}: Perfectly matched *")
                else:
                    db_match = False
                    all_match = False
                    diff_count = np.sum(v1_vals != v2_vals)
                    print(f"{col}: {diff_count}Variance")
        
        if db_match:
            print(f"✓ {db_name}All values are exactly the same.")
        
        print()
    
    print("=" * 80)
    if all_match:
        print("* The results of all tests are identical!")
    else:
        print("Differences were found.")
    print("=" * 80)


def main():
    # Loading results
    v1_norm = load_and_normalize_v1()
    v2_results = load_v2_results()
    
    # Keep standardized v1 results
    save_v1_normalized(v1_norm)
    
    # Comparison
    compare_databases(v1_norm, v2_results)


if __name__ == "__main__":
    main()
