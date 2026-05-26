#!/usr/bin/env python3
"""
完整比较 v1 和 v2 结果的脚本
"""

import sys
from pathlib import Path
import pandas as pd
import numpy as np


def load_and_normalize_v1():
    """加载并标准化 v1 结果"""
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
        
        # 标准化列名
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
    """加载 v2 结果"""
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
    """保存标准化的 v1 结果"""
    output_dir = Path(__file__).parent / "comparison_output"
    for db_name, df in v1_norm.items():
        out_file = output_dir / f"v1_{db_name}_normalized.tsv"
        df.to_csv(out_file, sep='\t', index=False)


def compare_databases(v1_norm, v2_results):
    """比较所有数据库的结果"""
    print("=" * 80)
    print("AllEnricher v1 vs v2 详细比较")
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
        
        print(f"  v1 条目数: {len(v1_df)}")
        print(f"  v2 条目数: {len(v2_df)}")
        
        # 按 Term_ID 对齐
        v1_sorted = v1_df.sort_values('Term_ID').set_index('Term_ID')
        v2_sorted = v2_df.sort_values('Term_ID').set_index('Term_ID')
        
        # 共同条目
        common_terms = v1_sorted.index.intersection(v2_sorted.index)
        v1_only = v1_sorted.index.difference(v2_sorted.index)
        v2_only = v2_sorted.index.difference(v1_sorted.index)
        
        print(f"  共同条目: {len(common_terms)}")
        if len(v1_only) > 0:
            print(f"  v1 特有: {len(v1_only)}")
        if len(v2_only) > 0:
            print(f"  v2 特有: {len(v2_only)}")
        
        if len(common_terms) == 0:
            print()
            continue
        
        # 只比较共同条目
        v1_common = v1_sorted.loc[common_terms]
        v2_common = v2_sorted.loc[common_terms]
        
        # 比较各个列
        columns_to_compare = ['P_Value', 'Adjusted_P_Value', 'Gene_Count', 'Background_Count']
        db_match = True
        
        for col in columns_to_compare:
            if col not in v1_common.columns or col not in v2_common.columns:
                continue
            
            v1_vals = v1_common[col]
            v2_vals = v2_common[col]
            
            if col in ['P_Value', 'Adjusted_P_Value']:
                # 浮点比较
                if np.allclose(v1_vals, v2_vals, rtol=1e-5, atol=1e-8):
                    print(f"  {col}: 完全匹配 ✓")
                else:
                    db_match = False
                    all_match = False
                    # 计算差异数
                    diff_mask = ~np.isclose(v1_vals, v2_vals, rtol=1e-5, atol=1e-8)
                    diff_count = np.sum(diff_mask)
                    print(f"  {col}: {diff_count} 个差异")
                    
                    # 显示前 3 个差异
                    diff_terms = v1_common.index[diff_mask][:3]
                    for term in diff_terms:
                        print(f"    {term}: v1={v1_common.loc[term, col]:.6g}, v2={v2_common.loc[term, col]:.6g}")
            else:
                # 整数比较
                if (v1_vals == v2_vals).all():
                    print(f"  {col}: 完全匹配 ✓")
                else:
                    db_match = False
                    all_match = False
                    diff_count = np.sum(v1_vals != v2_vals)
                    print(f"  {col}: {diff_count} 个差异")
        
        if db_match:
            print(f"  ✓ {db_name} 所有数值完全一致")
        
        print()
    
    print("=" * 80)
    if all_match:
        print("✓ 所有测试数据库的结果完全一致！")
    else:
        print("⚠ 发现一些差异")
    print("=" * 80)


def main():
    # 加载结果
    v1_norm = load_and_normalize_v1()
    v2_results = load_v2_results()
    
    # 保存标准化 v1 结果
    save_v1_normalized(v1_norm)
    
    # 比较
    compare_databases(v1_norm, v2_results)


if __name__ == "__main__":
    main()
