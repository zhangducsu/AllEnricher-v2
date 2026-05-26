#!/usr/bin/env python3
"""
调查P值差异的原因
"""

import pandas as pd
from pathlib import Path

# 读取结果
v1_file = Path(r"F:\OneDrive\Documents\TraeSOLO\AllEnricher\AllEnricher-v1\example\allenricher\fisher\Q0.05\example.glist.GO.xls")
v2_file = Path(r"F:\OneDrive\Documents\TraeSOLO\AllEnricher\AllEnricher-v2\test_output\v1v2_compare\GO_enrichment.tsv")

df_v1 = pd.read_csv(v1_file, sep='\t')
df_v2 = pd.read_csv(v2_file, sep='\t')

# 查找P值差异大的条目
merged = df_v1[['TermID', 'rawP', 'TermGeneNum', 'ObservedGeneNum', 'ExpectedGeneNum']].merge(
    df_v2[['Term_ID', 'P_Value', 'Gene_Count', 'Background_Count', 'Expected_Count']], 
    left_on='TermID',
    right_on='Term_ID',
    how='inner'
)

merged['pval_diff'] = abs(merged['rawP'] - merged['P_Value'])

# 显示差异最大的10个
print("P值差异最大的10个条目:")
print(merged.nlargest(10, 'pval_diff')[['TermID', 'rawP', 'P_Value', 'pval_diff', 
                                         'TermGeneNum', 'Background_Count',
                                         'ObservedGeneNum', 'Gene_Count',
                                         'ExpectedGeneNum', 'Expected_Count']].to_string(index=False))

print("\n" + "="*80)
print("P值最小的10个条目（富集最显著的）:")
print(merged.nsmallest(10, 'rawP')[['TermID', 'rawP', 'P_Value', 'pval_diff']].to_string(index=False))

# 检查富集因子计算
print("\n" + "="*80)
print("检查Rich Factor计算:")
merged['rich_factor_v1'] = merged['ObservedGeneNum'] / merged['ExpectedGeneNum']
merged['rich_factor_v2'] = merged['Gene_Count'] / merged['Expected_Count']
merged['rf_diff'] = abs(merged['rich_factor_v1'] - merged['rich_factor_v2'])

print("\nRich Factor差异:")
print(f"  最大差异: {merged['rf_diff'].max():.6f}")
print(f"  平均差异: {merged['rf_diff'].mean():.6f}")

# 显示前10个条目对比
print("\n" + "="*80)
print("前10个条目的详细对比:")
top10 = merged.nsmallest(10, 'rawP')[['TermID', 'rawP', 'P_Value', 'pval_diff',
                                       'ObservedGeneNum', 'Gene_Count',
                                       'TermGeneNum', 'Background_Count',
                                       'ExpectedGeneNum', 'Expected_Count',
                                       'rich_factor_v1', 'rich_factor_v2']]
print(top10.to_string(index=False))
