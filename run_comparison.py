#!/usr/bin/env python3
"""
运行v2富集分析并与v1结果进行对比
"""

import sys
import os
from pathlib import Path
import pandas as pd
import numpy as np

# 添加项目路径
project_root = Path(__file__).parent
sys.path.insert(0, str(project_root))

from allenricher.core.config import Config
from allenricher.core.enrichment import EnrichmentAnalyzer
from allenricher.database.manager import DatabaseManager


def load_v1_results(v1_results_dir):
    """加载v1的结果文件"""
    results = {}
    v1_dir = Path(v1_results_dir)
    
    db_files = {
        'GO': 'example.glist.GO.xls',
        'KEGG': 'example.glist.KEGG.xls',
        'Reactome': 'example.glist.Reactome.xls',
        'DO': 'example.glist.DO.xls',
        'DisGeNET': 'example.glist.DisGeNET.xls'
    }
    
    for db_name, filename in db_files.items():
        filepath = v1_dir / filename
        if filepath.exists():
            try:
                df = pd.read_csv(filepath, sep='\t')
                # 确保TermID是字符串格式
                if 'TermID' in df.columns:
                    df['TermID'] = df['TermID'].astype(str)
                results[db_name] = df
                print(f"✓ 加载v1 {db_name} 结果: {len(df)} 个条目")
            except Exception as e:
                print(f"✗ 加载v1 {db_name} 失败: {e}")
    
    return results


def main():
    print("=" * 80)
    print("AllEnricher v1 vs v2 对比分析")
    print("=" * 80)
    
    # 路径设置
    v1_dir = Path(__file__).parent.parent / "AllEnricher-v1"
    v1_results_dir = v1_dir / "example" / "allenricher" / "fisher" / "Q0.05"
    v1_db_dir = v1_dir / "database" / "organism" / "v20190612" / "hsa"
    v2_output_dir = Path(__file__).parent / "comparison_output"
    
    # 创建输出目录
    v2_output_dir.mkdir(exist_ok=True)
    
    # 1. 加载基因列表
    print("\n1. 加载基因列表...")
    gene_list_file = Path(__file__).parent / "example_genes.txt"
    with open(gene_list_file, 'r') as f:
        genes = [line.strip() for line in f if line.strip()]
    print(f"   输入基因数: {len(genes)}")
    
    # 2. 配置分析
    print("\n2. 配置分析...")
    config = Config(
        species="hsa",
        databases=["GO", "KEGG", "Reactome", "DO"],
        method="fisher",
        correction="BH",
        pvalue_cutoff=1.0,  # 不做截断，获取完整结果以便对比
        qvalue_cutoff=1.0,
        output_dir=str(v2_output_dir),
        generate_report=False
    )
    
    # 3. 加载数据库
    print("\n3. 加载数据库...")
    try:
        db_manager = DatabaseManager(
            database_dir=str(v1_db_dir),
            species="hsa"
        )
        db_manager.load_databases(["GO", "KEGG", "Reactome", "DO"])
        print("✓ 数据库加载成功")
    except Exception as e:
        print(f"✗ 数据库加载失败: {e}")
        print(f"  尝试使用test_db...")
        test_db_dir = Path(__file__).parent / "test_db"
        if test_db_dir.exists():
            try:
                db_manager = DatabaseManager(
                    database_dir=str(test_db_dir),
                    species="hsa"
                )
                db_manager.load_databases(["GO", "KEGG"])
                print("✓ test_db 加载成功")
            except Exception as e2:
                print(f"✗ test_db 也加载失败: {e2}")
                return
        else:
            print(f"✗ test_db 不存在: {test_db_dir}")
            return
    
    # 4. 运行v2分析
    print("\n4. 运行v2富集分析...")
    try:
        analyzer = EnrichmentAnalyzer(config)
        database_data = db_manager.get_all_term_data()
        
        # 加载背景基因
        background_genes = db_manager.get_background_genes()
        print(f"   背景基因数: {len(background_genes)}")
        
        # 运行分析
        v2_results = analyzer.run_analysis(
            gene_set=set(genes),
            background_set=background_genes,
            database_data=database_data
        )
        
        print(f"✓ 分析完成，结果数据库数: {len(v2_results)}")
        
        # 保存v2结果
        for db_name, df in v2_results.items():
            output_file = v2_output_dir / f"v2_{db_name}.tsv"
            df.to_csv(output_file, sep='\t', index=False)
            print(f"   保存 {db_name}: {len(df)} 个条目 → {output_file.name}")
        
    except Exception as e:
        print(f"✗ 分析失败: {e}")
        import traceback
        traceback.print_exc()
        return
    
    # 5. 加载v1结果进行对比
    print("\n5. 加载v1结果...")
    v1_results = load_v1_results(v1_results_dir)
    
    # 6. 对比结果
    print("\n" + "=" * 80)
    print("结果对比")
    print("=" * 80)
    
    comparison_report = []
    common_dbs = set(v1_results.keys()) & set(v2_results.keys())
    
    for db_name in sorted(common_dbs):
        print(f"\n[{db_name}]")
        v1_df = v1_results[db_name]
        v2_df = v2_results[db_name]
        
        # 基本统计
        print(f"  v1 条目数: {len(v1_df)}")
        print(f"  v2 条目数: {len(v2_df)}")
        
        # 尝试找到共同的列名进行对比
        v1_cols = set(v1_df.columns)
        v2_cols = set(v2_df.columns)
        
        # 找出共同条目
        term_id_col_v1 = 'TermID' if 'TermID' in v1_cols else None
        term_id_col_v2 = None
        for col in ['TermID', 'term_id', 'ID']:
            if col in v2_cols:
                term_id_col_v2 = col
                break
        
        if term_id_col_v1 and term_id_col_v2:
            v1_terms = set(v1_df[term_id_col_v1].astype(str))
            v2_terms = set(v2_df[term_id_col_v2].astype(str))
            
            common = v1_terms & v2_terms
            v1_only = v1_terms - v2_terms
            v2_only = v2_terms - v1_terms
            
            print(f"  共同条目: {len(common)}")
            print(f"  v1特有: {len(v1_only)}")
            print(f"  v2特有: {len(v2_only)}")
            
            # 找出显著条目（Q<0.05）
            if 'adjP' in v1_df.columns and 'adjusted_pvalue' in v2_df.columns:
                v1_sig = set(v1_df[v1_df['adjP'] < 0.05][term_id_col_v1].astype(str))
                v2_sig = set(v2_df[v2_df['adjusted_pvalue'] < 0.05][term_id_col_v2].astype(str))
                
                sig_common = v1_sig & v2_sig
                sig_v1_only = v1_sig - v2_sig
                sig_v2_only = v2_sig - v1_sig
                
                print(f"\n  显著条目 (Q<0.05):")
                print(f"    v1: {len(v1_sig)}")
                print(f"    v2: {len(v2_sig)}")
                print(f"    共同显著: {len(sig_common)}")
                print(f"    v1特有显著: {len(sig_v1_only)}")
                print(f"    v2特有显著: {len(sig_v2_only)}")
                
                # 保存对比结果
                comparison_report.append({
                    'database': db_name,
                    'v1_total': len(v1_df),
                    'v2_total': len(v2_df),
                    'common': len(common),
                    'v1_sig': len(v1_sig),
                    'v2_sig': len(v2_sig),
                    'sig_common': len(sig_common),
                    'sig_v1_only': len(sig_v1_only),
                    'sig_v2_only': len(sig_v2_only)
                })
                
                # 详细对比前几个条目
                if len(sig_common) > 0:
                    print(f"\n  共同显著条目对比 (前5个):")
                    
                    # 获取共同条目的数据
                    common_sig_v1 = v1_df[v1_df[term_id_col_v1].astype(str).isin(sig_common)].head()
                    common_sig_v2 = v2_df[v2_df[term_id_col_v2].astype(str).isin(sig_common)].head()
                    
                    for idx, row in common_sig_v1.head(5).iterrows():
                        term_id = str(row[term_id_col_v1])
                        v2_row = v2_df[v2_df[term_id_col_v2].astype(str) == term_id]
                        
                        if not v2_row.empty:
                            v2_row = v2_row.iloc[0]
                            term_name = row.get('TermName', row.get('term_name', row.get('Description', 'N/A')))
                            v1_p = row.get('rawP', row.get('pvalue', 'N/A'))
                            v1_q = row.get('adjP', 'N/A')
                            v2_p = v2_row.get('pvalue', v2_row.get('rawP', 'N/A'))
                            v2_q = v2_row.get('adjusted_pvalue', v2_row.get('adjP', 'N/A'))
                            
                            print(f"    {term_id} ({term_name[:30]}...)")
                            print(f"      v1: p={v1_p:.6g}, q={v1_q:.6g}")
                            print(f"      v2: p={v2_p:.6g}, q={v2_q:.6g}")
    
    # 7. 保存对比报告
    print("\n" + "=" * 80)
    if comparison_report:
        comparison_df = pd.DataFrame(comparison_report)
        report_file = v2_output_dir / "comparison_report.tsv"
        comparison_df.to_csv(report_file, sep='\t', index=False)
        print(f"✓ 对比报告已保存至: {report_file}")
    
    print("\n完成!")
    print(f"输出目录: {v2_output_dir}")


if __name__ == "__main__":
    main()
