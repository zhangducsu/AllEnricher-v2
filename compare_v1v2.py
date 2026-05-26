#!/usr/bin/env python3
"""
对比 AllEnricher v1 和 v2 的分析结果
"""

import pandas as pd
from pathlib import Path
import sys

def compare_database(db_name, v1_file, v2_file):
    """对比单个数据库的分析结果"""
    print(f"\n{'='*60}")
    print(f"数据库: {db_name}")
    print('='*60)
    
    # 读取v1结果
    try:
        df_v1 = pd.read_csv(v1_file, sep='\t')
        v1_count = len(df_v1)
        print(f"v1 条目数: {v1_count}")
    except Exception as e:
        print(f"v1 文件读取失败: {e}")
        v1_count = 0
        df_v1 = None
    
    # 读取v2结果
    try:
        df_v2 = pd.read_csv(v2_file, sep='\t')
        v2_count = len(df_v2)
        print(f"v2 条目数: {v2_count}")
    except Exception as e:
        print(f"v2 文件读取失败: {e}")
        v2_count = 0
        df_v2 = None
    
    if df_v1 is None or df_v2 is None:
        return
    
    # 对比条目数
    diff = v2_count - v1_count
    if diff == 0:
        print(f"✓ 条目数一致")
    else:
        print(f"⚠ 条目数差异: {diff:+d} (v2 {'多' if diff > 0 else '少'} {abs(diff)} 条)")
    
    # 获取Term ID列名
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
        
        print(f"\nTerm ID对比:")
        print(f"  共同条目: {len(common)}")
        print(f"  仅v1有: {len(v1_only)}")
        print(f"  仅v2有: {len(v2_only)}")
        
        if len(v1_only) > 0:
            print(f"  v1独有示例: {list(v1_only)[:5]}")
        if len(v2_only) > 0:
            print(f"  v2独有示例: {list(v2_only)[:5]}")
    
    # 对比P值（对于共同条目）
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
        # 合并对比
        merged = df_v1[[term_col, pval_col_v1]].merge(
            df_v2[[term_col, pval_col_v2]], 
            on=term_col, 
            suffixes=('_v1', '_v2')
        )
        
        if len(merged) > 0:
            # 计算P值差异
            pval_diff = abs(merged[pval_col_v1 + '_v1'] - merged[pval_col_v2 + '_v2'])
            max_diff = pval_diff.max()
            mean_diff = pval_diff.mean()
            
            print(f"\nP值对比 (共同条目 {len(merged)} 个):")
            print(f"  最大差异: {max_diff:.2e}")
            print(f"  平均差异: {mean_diff:.2e}")
            
            if max_diff < 1e-10:
                print(f"  ✓ P值高度一致 (差异 < 1e-10)")
            elif max_diff < 1e-5:
                print(f"  ✓ P值基本一致 (差异 < 1e-5)")
            else:
                print(f"  ⚠ P值存在明显差异")
                # 显示差异最大的几个
                merged['diff'] = pval_diff
                top_diff = merged.nlargest(3, 'diff')[[term_col, pval_col_v1 + '_v1', pval_col_v2 + '_v2', 'diff']]
                print(f"  差异最大的条目:")
                print(top_diff.to_string(index=False))

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
    
    print("="*60)
    print("AllEnricher v1 vs v2 结果对比")
    print("="*60)
    
    for db_name, (v1_file, v2_file) in databases.items():
        compare_database(
            db_name,
            v1_dir / v1_file,
            v2_dir / v2_file
        )
    
    print("\n" + "="*60)
    print("对比完成")
    print("="*60)

if __name__ == '__main__':
    main()
