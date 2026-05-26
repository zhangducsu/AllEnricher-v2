#!/usr/bin/env python3
"""
详细对比 AllEnricher v1 和 v2 的分析结果
处理不同的列名格式
"""

import pandas as pd
from pathlib import Path
import sys

def find_column(df, possible_names):
    """在DataFrame中查找匹配的列名"""
    for name in possible_names:
        if name in df.columns:
            return name
    return None

def compare_database(db_name, v1_file, v2_file):
    """对比单个数据库的分析结果"""
    print(f"\n{'='*70}")
    print(f"数据库: {db_name}")
    print('='*70)
    
    # 读取v1结果
    try:
        df_v1 = pd.read_csv(v1_file, sep='\t')
        v1_count = len(df_v1)
        print(f"v1 条目数: {v1_count}")
        print(f"v1 列名: {list(df_v1.columns)}")
    except Exception as e:
        print(f"v1 文件读取失败: {e}")
        return
    
    # 读取v2结果
    try:
        df_v2 = pd.read_csv(v2_file, sep='\t')
        v2_count = len(df_v2)
        print(f"v2 条目数: {v2_count}")
        print(f"v2 列名: {list(df_v2.columns)}")
    except Exception as e:
        print(f"v2 文件读取失败: {e}")
        return
    
    # 对比条目数
    diff = v2_count - v1_count
    if diff == 0:
        print(f"\n✓ 条目数一致: {v1_count}")
    else:
        print(f"\n⚠ 条目数差异: {diff:+d} (v2 {'多' if diff > 0 else '少'} {abs(diff)} 条)")
    
    # 查找Term ID列
    term_col_v1 = find_column(df_v1, ['TermID', 'Term ID', 'term_id', 'ID'])
    term_col_v2 = find_column(df_v2, ['Term_ID', 'Term ID', 'term_id', 'ID'])
    
    if not term_col_v1:
        print(f"⚠ 无法在v1中找到Term ID列")
        return
    if not term_col_v2:
        print(f"⚠ 无法在v2中找到Term ID列")
        return
    
    print(f"\nTerm ID列: v1='{term_col_v1}', v2='{term_col_v2}'")
    
    # 对比Term ID
    v1_terms = set(df_v1[term_col_v1].astype(str))
    v2_terms = set(df_v2[term_col_v2].astype(str))
    
    common = v1_terms & v2_terms
    v1_only = v1_terms - v2_terms
    v2_only = v2_terms - v1_terms
    
    print(f"\nTerm ID对比:")
    print(f"  共同条目: {len(common)}")
    print(f"  仅v1有: {len(v1_only)}")
    print(f"  仅v2有: {len(v2_only)}")
    
    if len(v1_only) > 0:
        print(f"  v1独有示例: {list(v1_only)[:3]}")
    if len(v2_only) > 0:
        print(f"  v2独有示例: {list(v2_only)[:3]}")
    
    if len(common) == v1_count == v2_count:
        print(f"\n✓ Term ID完全一致！")
    
    # 查找P值列
    pval_col_v1 = find_column(df_v1, ['rawP', 'P-Value', 'p_value', 'pvalue', 'PValue'])
    pval_col_v2 = find_column(df_v2, ['P_Value', 'P-Value', 'p_value', 'pvalue', 'PValue'])
    
    if not pval_col_v1:
        print(f"\n⚠ 无法在v1中找到P值列")
        return
    if not pval_col_v2:
        print(f"\n⚠ 无法在v2中找到P值列")
        return
    
    print(f"\nP值列: v1='{pval_col_v1}', v2='{pval_col_v2}'")
    
    # 对比P值
    merged = df_v1[[term_col_v1, pval_col_v1]].merge(
        df_v2[[term_col_v2, pval_col_v2]], 
        left_on=term_col_v1,
        right_on=term_col_v2,
        how='inner'
    )
    
    if len(merged) > 0:
        # 计算P值差异
        pval_v1 = merged[pval_col_v1].astype(float)
        pval_v2 = merged[pval_col_v2].astype(float)
        pval_diff = abs(pval_v1 - pval_v2)
        
        max_diff = pval_diff.max()
        mean_diff = pval_diff.mean()
        
        print(f"\nP值对比 (共同条目 {len(merged)} 个):")
        print(f"  最大差异: {max_diff:.2e}")
        print(f"  平均差异: {mean_diff:.2e}")
        
        if max_diff < 1e-10:
            print(f"  ✓ P值高度一致 (差异 < 1e-10)")
        elif max_diff < 1e-5:
            print(f"  ✓ P值基本一致 (差异 < 1e-5)")
        elif max_diff < 0.01:
            print(f"  ⚠ P值存在轻微差异 (差异 < 0.01)")
        else:
            print(f"  ⚠ P值存在明显差异 (差异 >= 0.01)")
            # 显示差异最大的几个
            merged['pval_diff'] = pval_diff
            top_diff = merged.nlargest(3, 'pval_diff')[[term_col_v1, pval_col_v1, pval_col_v2, 'pval_diff']]
            print(f"  差异最大的条目:")
            print(top_diff.to_string(index=False))
    
    # 对比Q值/校正P值
    qval_col_v1 = find_column(df_v1, ['adjP', 'Q-Value', 'q_value', 'qvalue', 'Adjusted_P'])
    qval_col_v2 = find_column(df_v2, ['Adjusted_P_Value', 'Q-Value', 'q_value', 'qvalue', 'adjP'])
    
    if qval_col_v1 and qval_col_v2:
        print(f"\nQ值列: v1='{qval_col_v1}', v2='{qval_col_v2}'")
        
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
            
            print(f"\nQ值对比 (共同条目 {len(merged_q)} 个):")
            print(f"  最大差异: {max_diff_q:.2e}")
            print(f"  平均差异: {mean_diff_q:.2e}")
            
            if max_diff_q < 1e-10:
                print(f"  ✓ Q值高度一致 (差异 < 1e-10)")
            elif max_diff_q < 1e-5:
                print(f"  ✓ Q值基本一致 (差异 < 1e-5)")
            elif max_diff_q < 0.01:
                print(f"  ⚠ Q值存在轻微差异 (差异 < 0.01)")
            else:
                print(f"  ⚠ Q值存在明显差异 (差异 >= 0.01)")

def main():
    v1_dir = Path(r"F:\OneDrive\Documents\TraeSOLO\AllEnricher\AllEnricher-v1\example\allenricher\fisher\Q0.05")
    v2_dir = Path(r"F:\OneDrive\Documents\TraeSOLO\AllEnricher\AllEnricher-v2\test_output\v1v2_compare")
    
    # 数据库文件映射
    databases = {
        'GO': ('example.glist.GO.xls', 'GO_enrichment.tsv'),
        'KEGG': ('example.glist.KEGG.xls', 'KEGG_enrichment.tsv'),
        'Reactome': ('example.glist.Reactome.xls', 'Reactome_enrichment.tsv'),
        'DO': ('example.glist.DO.xls', 'DO_enrichment.tsv'),
        'DisGeNET': ('example.glist.DisGeNET.xls', 'DisGeNET_enrichment.tsv'),
    }
    
    print("="*70)
    print("AllEnricher v1 vs v2 详细结果对比")
    print("="*70)
    
    for db_name, (v1_file, v2_file) in databases.items():
        compare_database(
            db_name,
            v1_dir / v1_file,
            v2_dir / v2_file
        )
    
    print("\n" + "="*70)
    print("对比完成")
    print("="*70)

if __name__ == '__main__':
    main()
