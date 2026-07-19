#!/usr/bin/env python3
"""Detailed script comparing v1 and v2 results"""

import pandas as pd
from pathlib import Path
import numpy as np


def compare_v1_v2_results():
    """Compare enrichment analysis of v1 and v2"""
    
    # Path Settings
    v1_dir = Path(__file__).parent.parent / "AllEnricher-v1" / "out"
    v2_dir = Path(__file__).parent / "comparison_output"
    output_dir = Path(__file__).parent / "comparison_output"
    output_dir.mkdir(exist_ok=True)
    
    databases = ["GO", "KEGG", "Reactome", "DO"]
    
    total_differences = 0
    
    print("=" * 80)
    print("AllEnricher v1 vs v2 is a detailed comparison")
    print("=" * 80)
    print()
    
    for db in databases:
        print(f"[ {db} ]")
        print("-" * 80)
        
        # Loading results
        v1_file = v1_dir / f"v1_{db}.tsv"
        v2_file = v2_dir / f"v2_{db}.tsv"
        
        if not v1_file.exists():
            print(f"V1 not found{db}Results")
            continue
        if not v2_file.exists():
            print(f"V2 not found{db}Results")
            continue
        
        # Read Results
        v1_df = pd.read_csv(v1_file, sep='\t')
        v2_df = pd.read_csv(v2_file, sep='\t')
        
        print(f"v1 Entry: {len(v1_df)}")
        print(f"v2 Entry: {len(v2_df)}")
        
        if len(v1_df) != len(v2_df):
            print(f"The number of entries is different!")
            total_differences += 1
        
        # Sort by Term_ID
        v1_sorted = v1_df.sort_values("Term_ID").reset_index(drop=True)
        v2_sorted = v2_df.sort_values("Term_ID").reset_index(drop=True)
        
        # Check for Term_ID to be fully consistent
        v1_terms = set(v1_sorted["Term_ID"])
        v2_terms = set(v2_sorted["Term_ID"])
        
        common_terms = v1_terms & v2_terms
        v1_only = v1_terms - v2_terms
        v2_only = v2_terms - v1_terms
        
        print(f"Common entry: {len(common_terms)}")
        if v1_only:
            print(f"v1 unique: {len(v1_only)}One.")
            total_differences += len(v1_only)
        if v2_only:
            print(f"v2 Unique: {len(v2_only)}One.")
            total_differences += len(v2_only)
        
        # Comparison of values with common entries
        if common_terms:
            v1_common = v1_sorted[v1_sorted["Term_ID"].isin(common_terms)].set_index("Term_ID")
            v2_common = v2_sorted[v2_sorted["Term_ID"].isin(common_terms)].set_index("Term_ID")
            
            # Compare Key Value Columns
            compare_columns = ["P_Value", "Adjusted_P_Value", "Gene_Count", 
                             "Background_Count", "Rich_Factor"]
            
            all_match = True
            
            for col in compare_columns:
                if col not in v1_common.columns or col not in v2_common.columns:
                    continue
                    
                diff_count = 0
                for term_id in common_terms:
                    v1_val = v1_common.loc[term_id, col]
                    v2_val = v2_common.loc[term_id, col]
                    
                    # Value comparison (consider floating point accuracy)
                    if isinstance(v1_val, (int, float)) and isinstance(v2_val, (int, float)):
                        if not np.isclose(v1_val, v2_val, rtol=1e-5, atol=1e-8):
                            diff_count += 1
                            if diff_count <= 5:  # Show only the top 5 differences
                                print(f"{col}Variance: {term_id}")
                                print(f"      v1: {v1_val}, v2: {v2_val}")
                    else:
                        if v1_val != v2_val:
                            diff_count += 1
                            if diff_count <= 5:
                                print(f"{col}Variance: {term_id}")
                                print(f"      v1: {v1_val}, v2: {v2_val}")
                
                if diff_count > 0:
                    all_match = False
                    total_differences += diff_count
                    print(f"{col}: {diff_count}Variance")
                else:
                    print(f"{col}: Perfectly matched *")
            
            if all_match:
                print(f"✓ {db}All values are exactly the same.")
        
        print()
    
    print("=" * 80)
    if total_differences == 0:
        print("* The results of all tests are identical!")
    else:
        print(f"Found{total_differences}Variance")
    print("=" * 80)


if __name__ == "__main__":
    compare_v1_v2_results()
