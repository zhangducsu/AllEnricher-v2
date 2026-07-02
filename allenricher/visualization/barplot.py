"""
富集分析条形图模块
==================

提供发表级水平柱状图绘制功能，支持 GO/KEGG/Reactome/DO/DisGeNET 五种数据库的类别着色。

依赖：matplotlib, pandas, numpy
"""

import logging
from typing import Dict, List, Optional, Tuple, Union

import matplotlib
matplotlib.use("Agg", force=True)
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from .plot_theme import PlotTheme, save_figure_dual
from .plot_utils import clean_pathway_label, term_figure_size

logger = logging.getLogger(__name__)


# -----------------------------------------------------------------------------
# 数据库类别配色方案 - 从 ColorConfig 获取
# -----------------------------------------------------------------------------

def _get_category_colors(database: str, palette: Optional[str] = None) -> Dict[str, str]:
    """
    获取指定数据库的类别配色方案

    Args:
        database: 数据库类型 (GO, KEGG, Reactome, DO, DisGeNET)
        palette: 配色方案名称 (default, bright, vibrant, muted)

    Returns:
        类别到颜色的映射字典
    """
    from .color_config import ColorConfig
    
    config = ColorConfig()
    database = database.upper()
    
    if database == "GO":
        return config.get_categorical_colors('go', palette=palette)
    elif database == "KEGG":
        return config.get_categorical_colors('kegg', palette=palette)
    elif database == "REACTOME":
        colors = config.get_colors(palette or 'default', n=1)
        return {"default": colors[0]}
    elif database == "DO":
        colors = config.get_colors(palette or 'default', n=1)
        return {"default": colors[0]}
    elif database == "DISGENET":
        colors = config.get_colors(palette or 'default', n=1)
        return {"default": colors[0]}
    else:
        colors = config.get_colors(palette or 'default', n=1)
        return {"default": colors[0]}


def _parse_go_category(term_str: str) -> str:
    """
    从 GO term 字符串中解析类别

    Args:
        term_str: 格式如 "Cellular Component|Extracellular Matrix" 或 "biological_process:cell cycle"

    Returns:
        标准化类别名称 (biological_process / cellular_component / molecular_function)
    """
    # 实际数据格式: "Cellular Component|Extracellular Matrix" (| 分隔)
    # 兼容格式: "biological_process:cell cycle" (: 分隔)
    for sep in ("|", ":"):
        if sep in term_str:
            category = term_str.split(sep)[0].strip().lower()
            # 标准化类别名称（空格 → 下划线）
            category = category.replace(" ", "_")
            # 验证是否为有效 GO 类别
            valid_categories = {"biological_process", "cellular_component", "molecular_function"}
            if category in valid_categories:
                return category
    return "default"


def _parse_kegg_category(term_str: str) -> str:
    """
    从 KEGG term 字符串中解析类别

    Args:
        term_str: 格式如 "Metabolism|Carbohydrate metabolism|Glycolysis"

    Returns:
        类别名称
    """
    if "|" in term_str:
        category = term_str.split("|")[0].strip()
        # 将空格转为下划线，以匹配 color_config.py 中字典 key 的格式
        category = category.replace(" ", "_")
        return category
    return "default"


def _auto_figsize(n_terms: int, base_width: float = 10.0) -> Tuple[float, float]:
    """
    根据 term 数量自动计算图表尺寸

    Args:
        n_terms: term 数量
        base_width: 基础宽度

    Returns:
        (width, height) 元组
    """
    # 每个 term 约 0.35 英寸高度，最小 6 英寸，最大 16 英寸
    return term_figure_size(n_terms, width=base_width, min_height=2.8, row_height=0.42, max_height=16.0)


def _save_figure(fig: plt.Figure, output_file: Optional[str], dpi: int = 300):
    """保存图表到文件（如指定了 output_file），同时输出PNG和PDF"""
    if output_file:
        png_path, pdf_path = save_figure_dual(fig, output_file, dpi)
        logger.info(f"图表已保存: {png_path}, {pdf_path}")


# -----------------------------------------------------------------------------
# plot_barplot - 富集分析水平柱状图
# -----------------------------------------------------------------------------

def plot_barplot(
    data: pd.DataFrame,
    database: str = "GO",
    top_n: int = 20,
    style: str = "default",
    palette: str = "default",
    title: Optional[str] = None,
    output_file: Optional[str] = None,
    figsize: Optional[Tuple[float, float]] = None,
    dpi: int = 300,
    show_gene_count: bool = True,
    gene_count_col: str = "gene_count",
    qvalue_col: str = "qvalue",
    term_col: str = "term",
    rich_factor_col: Optional[str] = "rich_factor",
) -> matplotlib.figure.Figure:
    """
    绘制富集分析水平柱状图

    展示富集 term 的 -log10(Q-value) 水平柱状图，支持按类别着色，
    并在柱状图右侧标注基因数。

    Args:
        data: 富集结果 DataFrame，需包含 term、qvalue、gene_count 等列
        database: 数据库类型 (GO, KEGG, Reactome, DO, DisGeNET)
        top_n: 展示 top N 个 term（按 qvalue 排序）
        style: 图表风格 (default, nature, science, presentation)
        palette: 配色方案 (default, bright, vibrant, muted)
        title: 图表标题（默认自动生成）
        output_file: 输出文件路径（可选）
        figsize: 图表尺寸（默认自动计算）
        dpi: 输出分辨率
        show_gene_count: 是否在柱状图右侧显示基因数
        gene_count_col: 基因数量列名
        qvalue_col: Q-value 列名
        term_col: term 名称列名
        rich_factor_col: 富集因子列名（可选）

    Returns:
        matplotlib.figure.Figure
    """
    # 复制数据避免修改原始数据
    df = data.copy()

    # 确保必要的列存在
    required_cols = [term_col, qvalue_col, gene_count_col]
    for col in required_cols:
        if col not in df.columns:
            raise ValueError(f"缺少必要列: {col}")

    # 按 qvalue 排序并取 top_n
    df = df.sort_values(by=qvalue_col, ascending=True).head(top_n)

    if len(df) == 0:
        fig, ax = plt.subplots(figsize=figsize or (8, 6))
        ax.text(0.5, 0.5, "No enrichment terms to display", ha="center", va="center")
        ax.set_title(title or f"{database} Enrichment")
        _save_figure(fig, output_file, dpi)
        return fig

    # 计算 -log10(qvalue)
    df["neg_log10_q"] = -np.log10(df[qvalue_col].clip(lower=1e-300))

    # 获取类别配色
    category_colors = _get_category_colors(database, palette)

    # 解析每个 term 的类别并分配颜色
    colors = []
    categories = []
    for term in df[term_col]:
        if database.upper() == "GO":
            category = _parse_go_category(str(term))
        elif database.upper() == "KEGG":
            category = _parse_kegg_category(str(term))
        else:
            category = "default"

        categories.append(category)
        colors.append(category_colors.get(category, category_colors.get("default", "#4477AA")))

    df["category"] = categories
    df["color"] = colors

    # 准备基因数标签
    if show_gene_count:
        if rich_factor_col and rich_factor_col in df.columns:
            gene_labels = [
                f"{int(gc)}/{rf:.2f}" if pd.notna(rf) else str(int(gc))
                for gc, rf in zip(df[gene_count_col], df[rich_factor_col])
            ]
        else:
            gene_labels = [str(int(gc)) for gc in df[gene_count_col]]
    else:
        gene_labels = None

    # 自动计算 figsize
    if figsize is None:
        figsize = _auto_figsize(len(df))

    # 应用风格设置
    with PlotTheme.context(style or 'nature'):
        # 创建图形
        fig, ax = plt.subplots(figsize=figsize)

        # 反转顺序使最小的 qvalue 在顶部
        y_positions = np.arange(len(df))
        values = df["neg_log10_q"].values[::-1]
        bar_colors = df["color"].values[::-1]
        terms = df[term_col].values[::-1]

        if gene_labels:
            gene_labels = gene_labels[::-1]

        # 绘制水平柱状图
        bars = ax.barh(
            y_positions,
            values,
            color=bar_colors,
            edgecolor="none",
            height=0.7,
        )

        # 设置 y 轴标签（term 名称）
        # 解析并清理 term 名称
        clean_terms = []
        for term in terms:
            term_str = str(term)
            if database.upper() == "GO":
                # GO: 移除前面的分类前缀（Cellular Component| 或 biological_process:），只保留描述
                for sep in ("|", ":"):
                    if sep in term_str:
                        parts = term_str.split(sep, 1)
                        category = parts[0].strip().lower().replace(" ", "_")
                        if category in {"biological_process", "cellular_component", "molecular_function"}:
                            term_str = parts[1].strip()
                            break
                clean_terms.append(term_str)
            elif database.upper() == "KEGG" and "|" in term_str:
                # KEGG: 提取最后一部分
                parts = term_str.split("|")
                clean_terms.append(parts[-1].strip())
            else:
                clean_terms.append(term_str)

        clean_terms = [clean_pathway_label(term) for term in terms]
        ax.set_yticks(y_positions)
        ax.set_yticklabels(clean_terms, fontsize=9)

        # 在柱状图右侧标注基因数
        if gene_labels:
            for i, (bar, label) in enumerate(zip(bars, gene_labels)):
                width = bar.get_width()
                ax.text(
                    width + 0.05,
                    bar.get_y() + bar.get_height() / 2,
                    label,
                    ha="left",
                    va="center",
                    fontsize=8,
                    color="#333333",
                )

        # 设置 x 轴
        ax.set_xlabel("-log10(Q-value)", fontsize=10)
        max_value = max(float(np.nanmax(values)), 1e-6)
        ax.set_xlim(0, max_value * (1.35 if len(df) == 1 else 1.25))

        # 设置标题
        if title is None:
            title = f"{database} Enrichment"
            if show_gene_count and rich_factor_col and rich_factor_col in data.columns:
                title += " (Gene# / Rich Factor)"
            elif show_gene_count:
                title += " (Gene Count)"
        ax.set_title(title, fontsize=12, fontweight="bold")

        # 添加图例（仅 GO 和 KEGG 有多类别）
        if len(df) > 1 and database.upper() in ["GO", "KEGG"] and len(category_colors) > 1:
            legend_elements = [
                plt.Rectangle((0, 0), 1, 1, facecolor=color, edgecolor="none", label=cat.replace("_", " ").title())
                for cat, color in category_colors.items()
            ]
            ax.legend(
                handles=legend_elements,
                loc="lower right",
                fontsize=8,
                frameon=True,
                fancybox=False,
                edgecolor="#CCCCCC",
            )

        # 美化边框
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.tick_params(labelsize=9)

        # 添加网格线
        ax.xaxis.grid(True, linestyle="--", alpha=0.3)
        ax.set_axisbelow(True)

        plt.tight_layout()
        _save_figure(fig, output_file, dpi)

    return fig
