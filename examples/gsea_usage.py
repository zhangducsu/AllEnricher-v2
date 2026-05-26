#!/usr/bin/env python3
"""
GSEA / ssGSEA / GSVA 使用示例

本示例展示如何使用 AllEnricher v2.0 进行基于基因集富集分析的方法，
包括 GSEA、ssGSEA 和 GSVA 三种方法的使用方式。

使用前请确保已安装依赖：
    pip install allenricher[all]

目录:
    1. 从表达矩阵创建排序基因列表（GSEA输入）
    2. 使用 GSEA 进行基因集富集分析
    3. 使用 ssGSEA 进行单样本基因集富集分析
    4. 使用 GSVA 进行基因集变异分析
    5. 结果解读说明
"""

import numpy as np
import pandas as pd
from typing import Dict, List, Tuple

from allenricher.core.config import Config, EnrichmentMethod
from allenricher.core.enrichment import EnrichmentAnalyzer


# ===========================================================================
# 1. 从表达矩阵创建排序基因列表（GSEA 输入）
# ===========================================================================

def create_ranked_gene_list(
    expression_matrix: pd.DataFrame,
    sample: str,
    method: str = "log2fc"
) -> List[Tuple[str, float]]:
    """
    从表达矩阵中为指定样本创建排序基因列表

    GSEA 需要一个排序的基因列表作为输入，每个基因附带一个权重（如 log2FC、
    t-statistic 等）。本函数支持多种排序策略。

    参数:
        expression_matrix: 表达矩阵（行=基因，列=样本）
        sample: 目标样本名称（列名）
        method: 排序策略
            - "log2fc": 基于与所有样本均值的 log2 倍变化
            - "signal_to_noise": 信噪比（Signal-to-Noise Ratio）
            - "ttest": 基于与参考组（其他样本）的 t 检验统计量
            - "mean": 直接使用表达值

    返回值:
        List[Tuple[str, float]]: 排序后的 (基因名, 权重) 列表，按权重降序排列

    示例:
        >>> expr = pd.DataFrame(...)  # 表达矩阵
        >>> ranked = create_ranked_gene_list(expr, "tumor_01", method="log2fc")
        >>> print(ranked[:5])  # 打印前5个基因
    """
    if sample not in expression_matrix.columns:
        raise ValueError(f"样本 '{sample}' 不在表达矩阵中")

    sample_expr = expression_matrix[sample]

    if method == "log2fc":
        # log2 倍变化：样本表达值 vs 所有样本均值
        mean_expr = expression_matrix.mean(axis=1)
        # 避免除以零，添加伪计数
        mean_expr = mean_expr.replace(0, 1e-10)
        log2fc = np.log2(sample_expr / mean_expr)
        ranked = sorted(
            zip(expression_matrix.index, log2fc.values),
            key=lambda x: x[1],
            reverse=True
        )

    elif method == "signal_to_noise":
        # 信噪比：(mean_tumor - mean_normal) / (std_tumor + std_normal)
        # 这里将目标样本作为 "tumor"，其他样本作为 "normal"
        other_samples = [c for c in expression_matrix.columns if c != sample]
        mean_normal = expression_matrix[other_samples].mean(axis=1)
        std_sample = expression_matrix[sample].std()
        std_normal = expression_matrix[other_samples].std(axis=1)
        snr = (sample_expr - mean_normal) / (std_sample + std_normal + 1e-10)
        ranked = sorted(
            zip(expression_matrix.index, snr.values),
            key=lambda x: x[1],
            reverse=True
        )

    elif method == "mean":
        # 直接按表达值排序
        ranked = sorted(
            zip(expression_matrix.index, sample_expr.values),
            key=lambda x: x[1],
            reverse=True
        )

    else:
        raise ValueError(f"未知的排序策略: {method}，支持: log2fc/signal_to_noise/mean")

    return ranked


def save_ranked_gene_list(
    ranked_genes: List[Tuple[str, float]],
    output_file: str
) -> None:
    """
    将排序基因列表保存为文件（TSV 格式，两列: 基因名 权重）

    参数:
        ranked_genes: 排序后的 (基因名, 权重) 列表
        output_file: 输出文件路径

    示例:
        >>> save_ranked_gene_list(ranked, "ranked_genes.tsv")
    """
    with open(output_file, 'w', encoding='utf-8') as f:
        f.write("gene\tweight\n")
        for gene, weight in ranked_genes:
            f.write(f"{gene}\t{weight:.6f}\n")
    print(f"排序基因列表已保存到: {output_file}")


# ===========================================================================
# 2. 使用 GSEA 进行基因集富集分析
# ===========================================================================

def run_gsea_example():
    """
    GSEA 分析示例

    GSEA（Gene Set Enrichment Analysis）是最经典的基因集富集分析方法，
    适用于差异表达分析后对排序基因列表进行通路富集分析。

    使用场景：
    - RNA-seq 差异表达分析后，对按 log2FC 排序的基因列表进行通路分析
    - 不需要预先设定差异基因阈值（如 |log2FC| > 1）
    - 可以发现整体表达趋势一致的通路

    使用方式：
    1. 通过 CLI: allenricher analyze -i ranked_genes.tsv -m gsea -s hsa -d GO,KEGG
    2. 通过 Python API（如下所示）
    """
    print("=" * 60)
    print("GSEA 分析示例")
    print("=" * 60)

    # ---- 步骤1：准备排序基因列表 ----
    # 方式A：从表达矩阵生成
    # expression_matrix = pd.read_csv("expression_matrix.tsv", sep="\t", index_col=0)
    # ranked_genes = create_ranked_gene_list(expression_matrix, "tumor_01", method="log2fc")
    # save_ranked_gene_list(ranked_genes, "ranked_genes.tsv")

    # 方式B：直接提供排序基因列表文件（TSV格式，两列: 基因名 权重）
    # ranked_genes_file = "ranked_genes.tsv"

    # ---- 步骤2：创建配置 ----
    config = Config(
        method="gsea",           # 使用 GSEA 方法
        species="hsa",           # 人类
        databases=["GO", "KEGG"],  # 使用 GO 和 KEGG 数据库
        gsea_permutations=1000,  # 排列检验次数（默认1000，越大越精确但越慢）
        gsea_min_size=10,        # 基因集最小大小
        gsea_max_size=500,       # 基因集最大大小
        output_dir="./results/gsea"
    )

    # ---- 步骤3：创建分析器并运行 ----
    # analyzer = EnrichmentAnalyzer(config)
    # gene_list = analyzer.load_gene_list("ranked_genes.tsv")
    # results = analyzer.run_analysis(gene_list, background_set, database_data)
    # analyzer.save_results("./results/gsea")

    print("GSEA 配置已创建")
    print(f"  方法: {config.method}")
    print(f"  排列次数: {config.gsea_permutations}")
    print(f"  基因集大小范围: [{config.gsea_min_size}, {config.gsea_max_size}]")
    print()

    # ---- CLI 使用方式 ----
    print("CLI 使用方式:")
    print("  allenricher analyze -i ranked_genes.tsv -m gsea -s hsa -d GO,KEGG")
    print("  allenricher analyze -i ranked_genes.tsv -m gsea -s hsa -d GO,KEGG \\")
    print("      --gsea-permutations 1000 --gsea-min-size 10 --gsea-max-size 500")
    print()


# ===========================================================================
# 3. 使用 ssGSEA 进行单样本基因集富集分析
# ===========================================================================

def run_ssgsea_example():
    """
    ssGSEA 分析示例

    ssGSEA（Single Sample GSEA）是 GSEA 的单样本变体，每个样本独立计算
    归一化富集分数（NES），适用于：
    - 单细胞测序数据分析
    - 肿瘤样本分型
    - 免疫浸润评分
    - 少量样本或无重复样本的场景

    与 GSEA 的区别：
    - 不需要排列检验，计算速度快
    - 每个样本独立计算，输出 NES 分数
    - NES 范围在 [-1, 1] 之间
    """
    print("=" * 60)
    print("ssGSEA 分析示例")
    print("=" * 60)

    # ---- 步骤1：准备表达矩阵 ----
    # 表达矩阵格式：行=基因，列=样本（TSV 或 CSV）
    # expression_matrix = pd.read_csv("expression_matrix.tsv", sep="\t", index_col=0)

    # ---- 步骤2：创建配置 ----
    config = Config(
        method="ssgsea",         # 使用 ssGSEA 方法
        species="hsa",
        databases=["GO", "KEGG"],
        gsea_min_size=10,
        gsea_max_size=500,
        output_dir="./results/ssgsea"
    )

    print("ssGSEA 配置已创建")
    print(f"  方法: {config.method}")
    print(f"  基因集大小范围: [{config.gsea_min_size}, {config.gsea_max_size}]")
    print()

    # ---- CLI 使用方式 ----
    print("CLI 使用方式:")
    print("  allenricher analyze -i genes.txt -m ssgsea -s hsa -d GO,KEGG \\")
    print("      -e expression_matrix.tsv")
    print()


# ===========================================================================
# 4. 使用 GSVA 进行基因集变异分析
# ===========================================================================

def run_gsva_example():
    """
    GSVA 分析示例

    GSVA（Gene Set Variation Analysis）是一种非参数方法，将基因表达矩阵
    转换为基因集活性分数矩阵。每个样本的每个基因集都有一个活性分数。

    支持三种方法变体：
    - gsva:  标准 GSVA 方法（默认），基于累积密度函数和随机游走
    - plage: Pathway Level Analysis of Gene Expression，基于主成分分析
    - zscore: 简单的 Z-score 方法，计算基因集内基因表达均值的标准分数

    适用场景：
    - 样本聚类和分型
    - 通路活性热图可视化
    - 与临床表型的关联分析
    - 生物标志物发现
    """
    print("=" * 60)
    print("GSVA 分析示例")
    print("=" * 60)

    # ---- 步骤1：准备表达矩阵 ----
    # 表达矩阵格式：行=基因，列=样本（TSV 或 CSV）
    # expression_matrix = pd.read_csv("expression_matrix.tsv", sep="\t", index_col=0)

    # ---- 步骤2：创建配置（标准 GSVA 方法） ----
    config_gsva = Config(
        method="gsva",              # 使用 GSVA 方法
        species="hsa",
        databases=["GO", "KEGG"],
        gsva_method="gsva",         # GSVA 方法变体: gsva / plage / zscore
        gsva_kcdf="Gaussian",       # 核密度核函数: Gaussian / Poisson
        gsva_tau=1.0,               # 核密度带宽参数
        gsea_min_size=10,
        gsea_max_size=500,
        output_dir="./results/gsva"
    )

    print("GSVA 配置已创建（标准方法）")
    print(f"  方法: {config_gsva.method}")
    print(f"  GSVA 变体: {config_gsva.gsva_method}")
    print(f"  核密度函数: {config_gsva.gsva_kcdf}")
    print(f"  带宽参数 tau: {config_gsva.gsva_tau}")
    print()

    # ---- PLAGE 变体示例 ----
    config_plage = Config(
        method="gsva",
        gsva_method="plage",        # 使用 PLAGE 变体
        gsea_min_size=10,
        gsea_max_size=500,
        output_dir="./results/gsva_plage"
    )
    print("GSVA PLAGE 变体配置已创建")

    # ---- Z-score 变体示例 ----
    config_zscore = Config(
        method="gsva",
        gsva_method="zscore",       # 使用 Z-score 变体
        gsea_min_size=10,
        gsea_max_size=500,
        output_dir="./results/gsva_zscore"
    )
    print("GSVA Z-score 变体配置已创建")
    print()

    # ---- Python API 使用方式 ----
    print("Python API 使用方式:")
    print("  from allenricher.core.gsva import GSVA")
    print("  gsva_analyzer = GSVA(method='gsva', kcdf='Gaussian', tau=1.0)")
    print("  scores = gsva_analyzer.run_gsva(expression_matrix, gene_sets)")
    print("  scores.to_csv('gsva_scores.tsv', sep='\\t')")
    print()

    # ---- CLI 使用方式 ----
    print("CLI 使用方式:")
    print("  allenricher analyze -i genes.txt -m gsva -s hsa -d GO,KEGG \\")
    print("      -e expression_matrix.tsv")
    print()


# ===========================================================================
# 5. 结果解读说明
# ===========================================================================

def explain_results():
    """
    GSEA / ssGSEA / GSVA 结果解读指南
    """
    print("=" * 60)
    print("结果解读说明")
    print("=" * 60)

    print("""
GSEA 结果解读:
--------------
- NES (Normalized Enrichment Score): 归一化富集分数
  - NES > 0: 基因集在排序列表的前部（上调方向）富集
  - NES < 0: 基因集在排序列表的后部（下调方向）富集
  - |NES| 越大，富集程度越强
- FDR q-value: 假发现率，< 0.25 通常被认为有意义（GSEA推荐标准）
- Leading Edge: 前沿基因，对富集分数贡献最大的基因子集
  - tag: 前沿基因占基因集总基因数的比例
  - list: 前沿基因占排序列表总基因数的比例
  - signal: 前沿基因的富集信号强度

ssGSEA 结果解读:
----------------
- NES (Normalized Enrichment Score): 归一化富集分数，范围 [-1, 1]
  - NES 接近 1: 基因集在样本中高度活跃（上调）
  - NES 接近 -1: 基因集在样本中高度抑制（下调）
  - NES 接近 0: 基因集在样本中无明显活性
- ssGSEA 不提供 p 值/FDR，需要通过其他方法进行统计检验

GSVA 结果解读:
--------------
- GSVA 分数: 连续的活性分数，无固定范围
  - 正值: 基因集在样本中激活
  - 负值: 基因集在样本中抑制
  - 分数的绝对值越大，活性越强
- GSVA 不提供 p 值，通常需要后续统计分析（如 t 检验、ANOVA 等）
  来比较不同组别之间的基因集活性差异

三种方法的比较:
--------------
| 特性        | GSEA          | ssGSEA        | GSVA          |
|-------------|---------------|---------------|---------------|
| 输入        | 排序基因列表   | 表达矩阵      | 表达矩阵      |
| 输出        | NES + FDR     | NES/样本      | 活性分数矩阵  |
| 统计检验    | 排列检验      | 无            | 无            |
| 适用场景    | 差异表达后    | 单样本分析    | 样本分型/聚类 |
| 计算速度    | 慢            | 中等          | 中等          |
| 多样本支持  | 全局分析      | 每样本独立    | 每样本独立    |
""")


# ===========================================================================
# 主入口
# ===========================================================================

if __name__ == "__main__":
    print("AllEnricher v2.0 - GSEA/ssGSEA/GSVA 使用示例")
    print("=" * 60)
    print()

    # 运行各方法的示例
    run_gsea_example()
    run_ssgsea_example()
    run_gsva_example()
    explain_results()
