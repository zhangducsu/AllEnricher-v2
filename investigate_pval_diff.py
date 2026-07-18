#!/usr/bin/env python3
"""
Reasons for variance in P-values investigated
"""

import pandas as pd
from pathlib import Path

# Read Results
v1_file = Path(r"F:\OneDrive\Documents\TraeSOLO\AllEnricher\AllEnricher-v1\example\allenricher\fisher\Q0.05\example.glist.GO.xls")
v2_file = Path(r"F:\OneDrive\Documents\TraeSOLO\AllEnricher\AllEnricher-v2\test_output\v1v2_compare\GO_enrichment.tsv")

df_v1 = pd.read_csv(v1_file, sep='\t')
df_v2 = pd.read_csv(v2_file, sep='\t')

# Finds entries with large P-value differences
merged = df_v1[['TermID', 'rawP', 'TermGeneNum', 'ObservedGeneNum', 'ExpectedGeneNum']].merge(
    df_v2[['Term_ID', 'P_Value', 'Gene_Count', 'Background_Count', 'Expected_Count']], 
    left_on='TermID',
    right_on='Term_ID',
    how='inner'
)

merged['pval_diff'] = abs(merged['rawP'] - merged['P_Value'])

# Shows the top 10 differences
print("The top 10 entries for P-value differences:")
print(merged.nlargest(10, 'pval_diff')[['TermID', 'rawP', 'P_Value', 'pval_diff', 
                                         'TermGeneNum', 'Background_Count',
                                         'ObservedGeneNum', 'Gene_Count',
                                         'ExpectedGeneNum', 'Expected_Count']].to_string(index=False))

print("\n" + "="*80)
print("The 10 entries with the smallest P value (the most significant) are:")
print(merged.nsmallest(10, 'rawP')[['TermID', 'rawP', 'P_Value', 'pval_diff']].to_string(index=False))

# Check for the rich factor calculation
print("\n" + "="*80)
print("Check Rich Factor calculations:")
merged['rich_factor_v1'] = merged['ObservedGeneNum'] / merged['ExpectedGeneNum']
merged['rich_factor_v2'] = merged['Gene_Count'] / merged['Expected_Count']
merged['rf_diff'] = abs(merged['rich_factor_v1'] - merged['rich_factor_v2'])

print("\nRich Factor difference:")
print(f"Maximum difference: {merged['rf_diff'].max(): .6f}")
print(f"Average variance: {merged['rf_diff'].mean(): .6f}")

# Show top 10 entries
print("\n" + "="*80)
print("Details of the top 10 entries:")
top10 = merged.nsmallest(10, 'rawP')[['TermID', 'rawP', 'P_Value', 'pval_diff',
                                       'ObservedGeneNum', 'Gene_Count',
                                       'TermGeneNum', 'Background_Count',
                                       'ExpectedGeneNum', 'Expected_Count',
                                       'rich_factor_v1', 'rich_factor_v2']]
print(top10.to_string(index=False))
