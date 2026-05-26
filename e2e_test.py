"""
GSEA/GSVA/ssGSEA 端对端测试

使用生成的测试数据运行三种分析方法的完整流程测试。
"""

import sys
import time
import json
from pathlib import Path
import pandas as pd
import numpy as np

# 添加项目路径
sys.path.insert(0, str(Path(__file__).parent))

from allenricher.core.enrichment import GSEA, SSGSEA
from allenricher.core.gsva import GSVA

# 测试数据路径
TEST_DATA_DIR = Path("test_data")

def load_test_data():
    """加载测试数据"""
    print("=" * 60)
    print("加载测试数据")
    print("=" * 60)
    
    # 1. 排序基因列表
    ranked_df = pd.read_csv(TEST_DATA_DIR / "ranked_genes.tsv", sep='\t')
    ranked_genes = ranked_df['gene'].tolist()
    gene_weights = dict(zip(ranked_df['gene'], ranked_df['weight']))
    print(f"✓ 排序基因列表: {len(ranked_genes)} genes")
    
    # 2. 表达矩阵
    expr_matrix = pd.read_csv(TEST_DATA_DIR / "expression_matrix.tsv", sep='\t', index_col=0)
    print(f"✓ 表达矩阵: {expr_matrix.shape[0]} genes x {expr_matrix.shape[1]} samples")
    
    # 3. 基因集
    gene_sets = {}
    with open(TEST_DATA_DIR / "gene_sets.gmt", 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            pathway = parts[0]
            genes = set(parts[2:])
            gene_sets[pathway] = genes
    print(f"✓ 基因集: {len(gene_sets)} sets")
    
    return ranked_genes, gene_weights, expr_matrix, gene_sets

def test_gsea(ranked_genes, gene_weights, gene_sets):
    """测试 GSEA"""
    print()
    print("=" * 60)
    print("测试 GSEA (基因集富集分析 - 排序基因列表)")
    print("=" * 60)
    
    start_time = time.time()
    
    # 创建 GSEA 分析器（减少置换次数加速测试）
    gsea = GSEA(permutations=100)
    
    # 只分析前5个通路
    subset_gene_sets = {k: v for i, (k, v) in enumerate(gene_sets.items()) if i < 5}
    
    # 对每个通路单独计算
    results = []
    for pathway_name, pathway_genes in subset_gene_sets.items():
        es, nes, pvalue, leading_edge = gsea.calculate_normalized_es(
            ranked_genes, pathway_genes, gene_weights
        )
        results.append({
            'term_id': pathway_name,
            'term_name': pathway_name,
            'es': es,
            'nes': nes,
            'pvalue': pvalue,
            'leading_edge': leading_edge,
            'gene_count': len(pathway_genes)
        })
    
    df_results = pd.DataFrame(results)
    elapsed = time.time() - start_time
    
    # 结果分析
    print(f"✓ GSEA 分析完成，耗时: {elapsed:.2f}s")
    print(f"✓ 结果条数: {len(df_results)}")
    
    if len(df_results) > 0:
        print()
        print("GSEA 结果预览:")
        print("-" * 80)
        display_cols = ["term_name", "es", "nes", "pvalue", "gene_count"]
        print(df_results[display_cols].to_string(index=False))
        
        # 统计
        print()
        print("统计信息:")
        print(f"  - 正向富集 (NES > 0): {(df_results['nes'] > 0).sum()}")
        print(f"  - 负向富集 (NES < 0): {(df_results['nes'] < 0).sum()}")
        print(f"  - 显著 (p < 0.05): {(df_results['pvalue'] < 0.05).sum()}")
        print(f"  - NES 范围: [{df_results['nes'].min():.3f}, {df_results['nes'].max():.3f}]")
    
    return df_results

def test_ssgsea(expr_matrix, gene_sets):
    """测试 ssGSEA"""
    print()
    print("=" * 60)
    print("测试 ssGSEA (单样本基因集富集分析 - 表达矩阵)")
    print("=" * 60)
    
    start_time = time.time()
    
    # 创建 ssGSEA 分析器
    ssgsea = SSGSEA(
        min_size=10,
        max_size=500
    )
    
    # 只分析前5个通路
    subset_gene_sets = {k: v for i, (k, v) in enumerate(gene_sets.items()) if i < 5}
    
    # 运行分析
    results = ssgsea.analyze_matrix(expr_matrix, subset_gene_sets)
    
    elapsed = time.time() - start_time
    
    # 结果分析
    print(f"✓ ssGSEA 分析完成，耗时: {elapsed:.2f}s")
    print(f"✓ 结果形状: {results.shape[0]} 通路 x {results.shape[1]} 样本")
    
    print()
    print("ssGSEA 结果预览:")
    print("-" * 80)
    print(results.round(3).to_string())
    
    # 统计
    print()
    print("统计信息:")
    print(f"  - 得分范围: [{results.values.min():.3f}, {results.values.max():.3f}]")
    print(f"  - 均值: {results.values.mean():.3f}")
    print(f"  - 标准差: {results.values.std():.3f}")
    
    # 样本间相关性
    sample_corr = results.corr()
    upper_tri = sample_corr.values[np.triu_indices_from(sample_corr.values, 1)]
    corr_range = f"[{upper_tri.min():.3f}, {upper_tri.max():.3f}]"
    print(f"  - 样本相关性范围: {corr_range}")
    
    return results

def test_gsva(expr_matrix, gene_sets):
    """测试 GSVA"""
    print()
    print("=" * 60)
    print("测试 GSVA (基因集变异分析 - 表达矩阵)")
    print("=" * 60)
    
    start_time = time.time()
    
    # 创建 GSVA 分析器
    gsva = GSVA(
        method="gsva",
        kcdf="Gaussian",
        tau=1.0,
        min_size=10,
        max_size=500
    )
    
    # 只分析前5个通路
    subset_gene_sets = {k: v for i, (k, v) in enumerate(gene_sets.items()) if i < 5}
    
    # 运行分析
    results = gsva.analyze_matrix(expr_matrix, subset_gene_sets)
    
    elapsed = time.time() - start_time
    
    # 结果分析
    print(f"✓ GSVA 分析完成，耗时: {elapsed:.2f}s")
    print(f"✓ 结果形状: {results.shape[0]} 通路 x {results.shape[1]} 样本")
    
    print()
    print("GSVA 结果预览:")
    print("-" * 80)
    print(results.round(3).to_string())
    
    # 统计
    print()
    print("统计信息:")
    print(f"  - 得分范围: [{results.values.min():.3f}, {results.values.max():.3f}]")
    print(f"  - 均值: {results.values.mean():.3f}")
    print(f"  - 标准差: {results.values.std():.3f}")
    
    # 样本间相关性
    sample_corr = results.corr()
    upper_tri = sample_corr.values[np.triu_indices_from(sample_corr.values, 1)]
    corr_range = f"[{upper_tri.min():.3f}, {upper_tri.max():.3f}]"
    print(f"  - 样本相关性范围: {corr_range}")
    
    return results

def test_gsva_methods(expr_matrix, gene_sets):
    """测试 GSVA 的三种方法变体"""
    print()
    print("=" * 60)
    print("测试 GSVA 方法变体 (PLAGE, Z-score)")
    print("=" * 60)
    
    # 只分析前3个通路
    subset_gene_sets = {k: v for i, (k, v) in enumerate(gene_sets.items()) if i < 3}
    
    # PLAGE 方法
    print()
    print("PLAGE 方法:")
    gsva_plage = GSVA(method="plage")
    results_plage = gsva_plage.analyze_matrix(expr_matrix, subset_gene_sets)
    print(f"  结果形状: {results_plage.shape}")
    print(f"  得分范围: [{results_plage.values.min():.3f}, {results_plage.values.max():.3f}]")
    
    # Z-score 方法
    print()
    print("Z-score 方法:")
    gsva_zscore = GSVA(method="zscore")
    results_zscore = gsva_zscore.analyze_matrix(expr_matrix, subset_gene_sets)
    print(f"  结果形状: {results_zscore.shape}")
    print(f"  得分范围: [{results_zscore.values.min():.3f}, {results_zscore.values.max():.3f}]")
    
    return results_plage, results_zscore

def generate_report(gsea_results, ssgsea_results, gsva_results):
    """生成测试报告"""
    print()
    print("=" * 60)
    print("生成测试报告")
    print("=" * 60)
    
    report = {
        "test_date": time.strftime("%Y-%m-%d %H:%M:%S"),
        "test_data": {
            "sorted_genes": 2000,
            "expression_matrix": "2000 genes x 6 samples",
            "gene_sets": 10,
            "samples": ["Normal_1", "Normal_2", "Normal_3", "Disease_1", "Disease_2", "Disease_3"]
        },
        "results": {}
    }
    
    # GSEA 结果
    if len(gsea_results) > 0:
        report["results"]["GSEA"] = {
            "status": "success",
            "pathways_analyzed": len(gsea_results),
            "positive_enrichment": int((gsea_results['nes'] > 0).sum()),
            "negative_enrichment": int((gsea_results['nes'] < 0).sum()),
            "significant_p05": int((gsea_results['pvalue'] < 0.05).sum()),
            "nes_range": [float(gsea_results['nes'].min()), float(gsea_results['nes'].max())],
            "top_results": gsea_results.head(5)[["term_name", "nes", "pvalue", "gene_count"]].to_dict('records')
        }
    
    # ssGSEA 结果
    report["results"]["ssGSEA"] = {
        "status": "success",
        "shape": list(ssgsea_results.shape),
        "score_range": [float(ssgsea_results.values.min()), float(ssgsea_results.values.max())],
        "mean_score": float(ssgsea_results.values.mean()),
        "std_score": float(ssgsea_results.values.std())
    }
    
    # GSVA 结果
    report["results"]["GSVA"] = {
        "status": "success",
        "shape": list(gsva_results.shape),
        "score_range": [float(gsva_results.values.min()), float(gsva_results.values.max())],
        "mean_score": float(gsva_results.values.mean()),
        "std_score": float(gsva_results.values.std())
    }
    
    # 保存报告
    report_path = TEST_DATA_DIR / "e2e_test_report.json"
    with open(report_path, 'w', encoding='utf-8') as f:
        json.dump(report, f, indent=2, ensure_ascii=False)
    
    print(f"✓ 报告已保存: {report_path}")
    
    return report

def main():
    """主函数"""
    print()
    print("╔" + "═" * 58 + "╗")
    print("║" + " " * 15 + "GSEA/GSVA/ssGSEA 端对端测试" + " " * 13 + "║")
    print("╚" + "═" * 58 + "╝")
    
    # 加载测试数据
    ranked_genes, gene_weights, expr_matrix, gene_sets = load_test_data()
    
    # 运行测试
    gsea_results = test_gsea(ranked_genes, gene_weights, gene_sets)
    ssgsea_results = test_ssgsea(expr_matrix, gene_sets)
    gsva_results = test_gsva(expr_matrix, gene_sets)
    results_plage, results_zscore = test_gsva_methods(expr_matrix, gene_sets)
    
    # 生成报告
    report = generate_report(gsea_results, ssgsea_results, gsva_results)
    
    # 打印总结
    print()
    print("╔" + "═" * 58 + "╗")
    print("║" + " " * 15 + "测试总结" + " " * 30 + "║")
    print("╚" + "═" * 58 + "╝")
    print()
    print(f"✓ GSEA:  {report['results']['GSEA']['pathways_analyzed']} 条通路, "
          f"{report['results']['GSEA']['significant_p05']} 条显著 (p<0.05)")
    print(f"✓ ssGSEA: {report['results']['ssGSEA']['shape'][0]} 通路 x "
          f"{report['results']['ssGSEA']['shape'][1]} 样本")
    print(f"✓ GSVA:  {report['results']['GSVA']['shape'][0]} 通路 x "
          f"{report['results']['GSVA']['shape'][1]} 样本")
    print()
    print("所有测试完成！✓")
    print()
    
    return report

if __name__ == "__main__":
    report = main()
