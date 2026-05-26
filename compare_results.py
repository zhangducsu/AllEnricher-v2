#!/usr/bin/env python3
"""比较 v1 和 v2 结果的详细脚本"""

import pandas as pd
from pathlib import Path
import numpy as np


def compare_v1_v2_results():
    """比较 v1 和 v2 的富集分析结果"""
    
    # 路径设置
    v1_dir = Path(__file__).parent.parent / "AllEnricher-v1" / "out"
    v2_dir = Path(__file__).parent / "comparison_output"
    output_dir = Path(__file__).parent / "comparison_output"
    output_dir.mkdir(exist_ok=True)
    
    databases = ["GO", "KEGG", "Reactome", "DO"]
    
    total_differences = 0
    
    print("=" * 80)
    print("AllEnricher v1 vs v2 详细比较")
    print("=" * 80)
    print()
    
    for db in databases:
        print(f"[ {db} ]")
        print("-" * 80)
        
        # 加载结果
        v1_file = v1_dir / f"v1_{db}.tsv"
        v2_file = v2_dir / f"v2_{db}.tsv"
        
        if not v1_file.exists():
            print(f"✗ 找不到 v1 {db} 结果")
            continue
        if not v2_file.exists():
            print(f"✗ 找不到 v2 {db} 结果")
            continue
        
        # 读取结果
        v1_df = pd.read_csv(v1_file, sep='\t')
        v2_df = pd.read_csv(v2_file, sep='\t')
        
        print(f"  v1 条目数: {len(v1_df)}")
        print(f"  v2 条目数: {len(v2_df)}")
        
        if len(v1_df) != len(v2_df):
            print(f"  ⚠ 条目数不同!")
            total_differences += 1
        
        # 按 Term_ID 排序
        v1_sorted = v1_df.sort_values("Term_ID").reset_index(drop=True)
        v2_sorted = v2_df.sort_values("Term_ID").reset_index(drop=True)
        
        # 检查 Term_ID 是否完全一致
        v1_terms = set(v1_sorted["Term_ID"])
        v2_terms = set(v2_sorted["Term_ID"])
        
        common_terms = v1_terms & v2_terms
        v1_only = v1_terms - v2_terms
        v2_only = v2_terms - v1_terms
        
        print(f"  共同条目数: {len(common_terms)}")
        if v1_only:
            print(f"  v1 独有的: {len(v1_only)} 个")
            total_differences += len(v1_only)
        if v2_only:
            print(f"  v2 独有的: {len(v2_only)} 个")
            total_differences += len(v2_only)
        
        # 对共同条目比较数值
        if common_terms:
            v1_common = v1_sorted[v1_sorted["Term_ID"].isin(common_terms)].set_index("Term_ID")
            v2_common = v2_sorted[v2_sorted["Term_ID"].isin(common_terms)].set_index("Term_ID")
            
            # 比较关键数值列
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
                    
                    # 数值比较 (考虑浮点精度)
                    if isinstance(v1_val, (int, float)) and isinstance(v2_val, (int, float)):
                        if not np.isclose(v1_val, v2_val, rtol=1e-5, atol=1e-8):
                            diff_count += 1
                            if diff_count <= 5:  # 只显示前5个差异
                                print(f"    {col} 差异: {term_id}")
                                print(f"      v1: {v1_val}, v2: {v2_val}")
                    else:
                        if v1_val != v2_val:
                            diff_count += 1
                            if diff_count <= 5:
                                print(f"    {col} 差异: {term_id}")
                                print(f"      v1: {v1_val}, v2: {v2_val}")
                
                if diff_count > 0:
                    all_match = False
                    total_differences += diff_count
                    print(f"    {col}: {diff_count} 个差异")
                else:
                    print(f"    {col}: 完全匹配 ✓")
            
            if all_match:
                print(f"  ✓ {db} 所有数值完全一致")
        
        print()
    
    print("=" * 80)
    if total_differences == 0:
        print("✓ 所有测试数据库的结果完全一致！")
    else:
        print(f"共发现 {total_differences} 个差异")
    print("=" * 80)


if __name__ == "__main__":
    compare_v1_v2_results()
