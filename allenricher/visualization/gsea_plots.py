"""
GSEA 发表级可视化模块
=====================

提供三种 GSEA 核心图表：
- plot_gsea_enrichment: GSEA 富集曲线图（三面板：ES曲线 / 基因集成员位置 / 基因权重）
- plot_gsea_nes_barplot: NES 条形图（多通路标准化富集分数比较）
- plot_gsea_dotplot: GSEA 气泡图（NES / 基因数量 / 显著性）

依赖：matplotlib, seaborn（不依赖 R 环境）
"""

import logging
from typing import Dict, List, Optional, Set

import matplotlib
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np
import pandas as pd
import seaborn as sns

from .plot_theme import PlotTheme, save_figure_dual

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# 工具函数
# ---------------------------------------------------------------------------

def _compute_running_es(
    ranked_genes: List[str],
    gene_set: Set[str],
    gene_weights: Optional[Dict[str, float]] = None,
) -> np.ndarray:
    """
    计算沿排序列表的 running enrichment score 曲线

    算法与 GSEA.calculate_enrichment_score 一致：
    - 命中增量 hit_inc = 1 / N_hit
    - 未命中增量 miss_inc = 1 / (N - nh)
    - 遍历排序列表，命中加分，未命中减分

    Args:
        ranked_genes: 排序基因列表
        gene_set: 基因集成员
        gene_weights: 可选基因权重

    Returns:
        np.ndarray: 长度为 len(ranked_genes) 的 running ES 数组
    """
    n = len(ranked_genes)
    nh = len(gene_set & set(ranked_genes))

    if nh == 0:
        return np.zeros(n)

    # 计算命中增量和未命中增量
    nr = sum(
        abs(gene_weights.get(g, 1.0)) for g in gene_set if g in gene_weights
    ) if gene_weights else nh
    hit_inc = 1.0 / nr if nr > 0 else 0
    miss_inc = 1.0 / (n - nh) if (n - nh) > 0 else 0

    running_es = np.zeros(n)
    running_sum = 0.0

    for i, gene in enumerate(ranked_genes):
        if gene in gene_set:
            weight = abs(gene_weights.get(gene, 1.0)) if gene_weights else 1.0
            running_sum += hit_inc * weight
        else:
            running_sum -= miss_inc
        running_es[i] = running_sum

    return running_es


def _get_gene_set_positions(
    ranked_genes: List[str],
    gene_set: Set[str],
) -> List[int]:
    """返回基因集成员在排序列表中的索引位置（0-based）"""
    return [i for i, g in enumerate(ranked_genes) if g in gene_set]


def _save_figure(fig: plt.Figure, output_file: Optional[str], dpi: int = 300):
    """保存图表到文件（如指定了 output_file）"""
    if output_file:
        save_figure_dual(fig, output_file, dpi=dpi)
        logger.info(f"图表已保存: {output_file}")


# ---------------------------------------------------------------------------
# plot_gsea_enrichment — GSEA 富集曲线图（三面板）
# ---------------------------------------------------------------------------

def plot_gsea_enrichment(
    ranked_genes: List[str],
    gene_weights: Dict[str, float],
    gene_set: Set[str],
    es: float,
    nes: float,
    pvalue: float,
    title: str = None,
    output_file: str = None,
    figsize: tuple = (10, 6),
    style: Optional[str] = None,
    palette: Optional[str] = None,
    dpi: int = 300,
) -> matplotlib.figure.Figure:
    """
    GSEA 富集曲线图（三面板）

    面板1（上）: Running Enrichment Score 曲线
    面板2（中）: 基因集成员位置标记（黑色竖线）
    面板3（下）: 基因排序度量值（权重柱状图）

    Args:
        ranked_genes: 排序基因列表
        gene_weights: 基因权重字典 {gene: weight}
        gene_set: 基因集成员集合
        es: 富集分数
        nes: 标准化富集分数
        pvalue: P 值
        title: 图表标题（默认自动生成）
        output_file: 输出文件路径（可选）
        figsize: 图表尺寸
        style: 图表风格（默认 nature）
        palette: 色板名（覆盖风格默认色板）

    Returns:
        matplotlib.figure.Figure
    """
    # 计算 running ES
    running_es = _compute_running_es(ranked_genes, gene_set, gene_weights)
    positions = _get_gene_set_positions(ranked_genes, gene_set)

    # 基因权重数组
    weights = np.array([gene_weights.get(g, 0.0) for g in ranked_genes])

    with PlotTheme.context(style or 'nature'):
        # 获取颜色
        color = PlotTheme.get_palette(palette, n=1)[0]
        color_pos, color_neg = PlotTheme.get_palette(palette, n=2)

        # 创建图形和 gridspec（高度比例 3:1:2）
        fig = plt.figure(figsize=figsize)
        gs = gridspec.GridSpec(
            3, 1, height_ratios=[3, 1, 2], hspace=0.05,
            left=0.12, right=0.95, top=0.90, bottom=0.08,
        )

        x = np.arange(len(ranked_genes))

        # ---- 面板1: Running Enrichment Score ----
        ax1 = fig.add_subplot(gs[0])
        ax1.fill_between(x, 0, running_es, where=(running_es >= 0),
                         color=color, alpha=0.3)
        ax1.fill_between(x, 0, running_es, where=(running_es < 0),
                         color=color_pos, alpha=0.3)
        ax1.plot(x, running_es, color=color, linewidth=1.0)
        ax1.axhline(y=0, color="gray", linestyle="--", linewidth=0.5)

        # 标注峰值 ES
        peak_idx = np.argmax(np.abs(running_es))
        peak_val = running_es[peak_idx]
        ax1.annotate(
            f"ES = {peak_val:.3f}",
            xy=(peak_idx, peak_val),
            xytext=(peak_idx + len(ranked_genes) * 0.05, peak_val),
            fontsize=8,
            arrowprops=dict(arrowstyle="->", color="black", lw=0.8),
        )

        ax1.set_ylabel("Running Enrichment Score", fontsize=9)
        ax1.set_xlim(0, len(ranked_genes) - 1)
        ax1.set_xticklabels([])
        ax1.tick_params(labelsize=8)

        # ---- 面板2: 基因集成员位置 ----
        ax2 = fig.add_subplot(gs[1], sharex=ax1)
        for pos in positions:
            ax2.axvline(x=pos, color="black", linewidth=0.5, alpha=0.7)
        ax2.set_xlim(0, len(ranked_genes) - 1)
        ax2.set_ylim(0, 1)
        ax2.set_yticks([])
        ax2.set_xticklabels([])
        ax2.set_ylabel("Gene Set", fontsize=9)
        ax2.tick_params(labelsize=8)

        # ---- 面板3: 基因权重柱状图 ----
        ax3 = fig.add_subplot(gs[2], sharex=ax1)
        bar_colors = [color_pos if w >= 0 else color_neg for w in weights]
        ax3.bar(x, weights, width=1.0, color=bar_colors, alpha=0.6)
        ax3.axhline(y=0, color="gray", linewidth=0.5)
        ax3.set_xlim(0, len(ranked_genes) - 1)
        ax3.set_ylabel("Weight (log2FC)", fontsize=9)
        ax3.set_xlabel("Rank in Ordered Gene List", fontsize=9)
        ax3.tick_params(labelsize=8)

        # ---- 标题 ----
        if title is None:
            title = f"NES = {nes:.2f}, P-value = {pvalue:.2e}"
        fig.suptitle(title, fontsize=11, fontweight="bold", y=0.96)

        _save_figure(fig, output_file, dpi=dpi)
    return fig


# ---------------------------------------------------------------------------
# plot_gsea_nes_barplot — NES 条形图
# ---------------------------------------------------------------------------

def plot_gsea_nes_barplot(
    results_df: pd.DataFrame,
    top_n: int = 20,
    title: str = "GSEA NES Ranking",
    output_file: str = None,
    figsize: tuple = (10, 8),
    style: Optional[str] = None,
    palette: Optional[str] = None,
    dpi: int = 300,
) -> matplotlib.figure.Figure:
    """
    GSEA NES 条形图

    按 NES 绝对值降序排列，展示多个通路的标准化富集分数。
    NES > 0 红色（正向富集），NES < 0 蓝色（负向富集）。
    添加显著性星号标注。

    Args:
        results_df: 包含 pathway, nes, pvalue 列的 DataFrame
        top_n: 展示前 N 个通路
        title: 图表标题
        output_file: 输出文件路径（可选）
        figsize: 图表尺寸
        style: 图表风格（默认 nature）
        palette: 色板名（覆盖风格默认色板）

    Returns:
        matplotlib.figure.Figure
    """
    df = results_df.copy()

    # 动态检测列名：支持 clusterProfiler 命名和旧命名
    nes_col = None
    for c in ['NES', 'nes', 'enrichmentScore']:
        if c in df.columns:
            nes_col = c
            break
    pval_col = None
    for c in ['p_value', 'NOM p-val', 'pvalue', 'P_Value']:
        if c in df.columns:
            pval_col = c
            break
    pathway_col = None
    for c in ['Description', 'pathway', 'Term_Name', 'term_name']:
        if c in df.columns:
            pathway_col = c
            break

    # 按 NES 绝对值降序排序
    df["_abs_nes"] = df[nes_col].abs() if nes_col else df.iloc[:, 0]
    df = df.sort_values("_abs_nes", ascending=True).tail(top_n).copy()
    df = df.drop(columns=["_abs_nes"])

    with PlotTheme.context(style or 'nature'):
        # 获取颜色
        colors = PlotTheme.get_palette(palette, n=3)
        color_pos = colors[2]
        color_neg = colors[0]

        # 颜色映射
        bar_colors = [color_pos if v >= 0 else color_neg for v in df[nes_col]]

        # 显著性标注
        def _sig_stars(p):
            if p < 0.001:
                return "***"
            elif p < 0.01:
                return "**"
            elif p < 0.05:
                return "*"
            return ""

        sig_labels = [_sig_stars(p) for p in df[pval_col]]

        # 绘图
        fig, ax = plt.subplots(figsize=figsize)
        bars = ax.barh(range(len(df)), df[nes_col], color=bar_colors, edgecolor="none", height=0.7)

        # 添加显著性标注
        for i, (nes_val, sig) in enumerate(zip(df[nes_col], sig_labels)):
            offset = 0.05 * (df[nes_col].max() - df[nes_col].min()) if len(df) > 1 else 0.1
            if sig:
                ax.text(
                    nes_val + (offset if nes_val >= 0 else -offset),
                    i,
                    sig,
                    va="center",
                    ha="left" if nes_val >= 0 else "right",
                    fontsize=8,
                )

        ax.set_yticks(range(len(df)))
        ax.set_yticklabels(df[pathway_col], fontsize=8)
        ax.axvline(x=0, color="gray", linewidth=0.8, linestyle="-")
        ax.set_xlabel("Normalized Enrichment Score (NES)", fontsize=10)
        ax.set_title(title, fontsize=12, fontweight="bold")
        ax.tick_params(labelsize=8)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

        plt.tight_layout()
        _save_figure(fig, output_file, dpi=dpi)
    return fig


# ---------------------------------------------------------------------------
# plot_gsea_dotplot — GSEA 气泡图
# ---------------------------------------------------------------------------

def plot_gsea_dotplot(
    results_df: pd.DataFrame,
    top_n: int = 20,
    title: str = "GSEA Dot Plot",
    output_file: str = None,
    figsize: tuple = (10, 8),
    style: Optional[str] = None,
    palette: Optional[str] = None,
    dpi: int = 300,
) -> matplotlib.figure.Figure:
    """
    GSEA 气泡图

    X 轴: NES 值
    Y 轴: 通路名称
    点大小: 基因数量
    点颜色: -log10(pvalue) 显著性

    Args:
        results_df: 包含 pathway, nes, pvalue, gene_count 列的 DataFrame
        top_n: 展示前 N 个通路
        title: 图表标题
        output_file: 输出文件路径（可选）
        figsize: 图表尺寸
        style: 图表风格（默认 nature）
        palette: 色板名（覆盖风格默认色板）

    Returns:
        matplotlib.figure.Figure
    """
    df = results_df.copy()

    # 动态检测列名：支持 clusterProfiler 命名和旧命名
    nes_col = None
    for c in ['NES', 'nes', 'enrichmentScore']:
        if c in df.columns:
            nes_col = c
            break
    pval_col = None
    for c in ['p_value', 'NOM p-val', 'pvalue', 'P_Value']:
        if c in df.columns:
            pval_col = c
            break
    pathway_col = None
    for c in ['Description', 'pathway', 'Term_Name', 'term_name']:
        if c in df.columns:
            pathway_col = c
            break
    gcount_col = None
    for c in ['setSize', 'gene_count', 'Gene_Count', 'GeneCount']:
        if c in df.columns:
            gcount_col = c
            break

    # 按 NES 绝对值降序排序
    df["_abs_nes"] = df[nes_col].abs() if nes_col else df.iloc[:, 0]
    df = df.sort_values("_abs_nes", ascending=True).tail(top_n).copy()
    df = df.drop(columns=["_abs_nes"])

    # 计算 -log10(pvalue)，避免 log(0)
    df["neg_log10_p"] = -np.log10(df[pval_col].clip(lower=1e-300))

    with PlotTheme.context(style or 'nature'):
        # 绘图
        fig, ax = plt.subplots(figsize=figsize)

        scatter = ax.scatter(
            df[nes_col],
            range(len(df)),
            s=df[gcount_col] * 3 if gcount_col else 30,
            c=df["neg_log10_p"],
            cmap=PlotTheme.get_diverging_cmap(palette),
            edgecolors="gray",
            linewidths=0.3,
            alpha=0.85,
        )

        ax.set_yticks(range(len(df)))
        ax.set_yticklabels(df[pathway_col], fontsize=8)
        ax.axvline(x=0, color="gray", linewidth=0.8, linestyle="-")
        ax.set_xlabel("Normalized Enrichment Score (NES)", fontsize=10)
        ax.set_title(title, fontsize=12, fontweight="bold")
        ax.tick_params(labelsize=8)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

        # 添加颜色条
        cbar = plt.colorbar(scatter, ax=ax, shrink=0.6, pad=0.02)
        cbar.set_label("-log10(P-value)", fontsize=9)
        cbar.ax.tick_params(labelsize=8)

        plt.tight_layout()
        _save_figure(fig, output_file, dpi=dpi)
    return fig
