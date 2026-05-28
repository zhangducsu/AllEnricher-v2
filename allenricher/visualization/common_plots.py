"""
通用可视化组件模块
==================

提供所有富集分析方法共用的可视化图表：
- plot_enrichment_network: 通路关系网络图（基于基因重叠度）
- plot_upset: UpSet 图（基因集交集可视化）
- plot_volcano: 火山图（NES/log2FC vs -log10 pvalue）
- plot_method_comparison: 方法间比较散点图

依赖：matplotlib, seaborn, numpy, scipy, networkx（可选）
"""

import logging
from itertools import combinations
from typing import Dict, Optional, Set

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy import stats

from .plot_theme import PlotTheme, save_figure_dual
from .color_config import ColorConfig, VOLCANO_COLORS

logger = logging.getLogger(__name__)

# 尝试导入 networkx，不可用时回退到纯 matplotlib 实现
try:
    import networkx as nx
    HAS_NETWORKX = True
except ImportError:
    HAS_NETWORKX = False
    logger.warning("networkx 未安装，plot_enrichment_network 将使用简化布局")


# ---------------------------------------------------------------------------
# 工具函数
# ---------------------------------------------------------------------------

def _save_figure(fig: plt.Figure, output_file: Optional[str], dpi: int = 300):
    """保存图表到文件（如指定了 output_file）"""
    if output_file:
        save_figure_dual(fig, output_file, dpi=dpi)
        logger.info(f"图表已保存: {output_file}")


# ---------------------------------------------------------------------------
# plot_enrichment_network - 通路关系网络图
# ---------------------------------------------------------------------------

def plot_enrichment_network(
    gene_sets: Dict[str, Set[str]],
    results_df: pd.DataFrame = None,
    top_n: int = 30,
    min_overlap: int = 3,
    layout: str = "spring",
    title: str = "Pathway Network",
    output_file: str = None,
    figsize: tuple = (12, 10),
    dpi: int = 300,
    style: Optional[str] = None,
    palette: Optional[str] = None,
) -> matplotlib.figure.Figure:
    """
    通路关系网络图

    基于通路间基因重叠度（Jaccard 系数）构建网络图，节点大小映射通路基因数量，
    节点颜色映射 NES 或显著性（如有 results_df），边的粗细映射重叠度。

    Args:
        gene_sets: 通路名到基因集合的映射 {pathway_name: set(genes)}
        results_df: 可选，包含 pathway, nes/pvalue 列的 DataFrame
        top_n: 展示基因数量最多的前 N 个通路
        min_overlap: 最小基因重叠数，低于此值不画边
        layout: 布局算法，"spring", "circular", "kamada_kawai"
        title: 图表标题
        output_file: 输出文件路径（可选）
        figsize: 图表尺寸
        dpi: 输出分辨率
        style: 风格名称（可选，如 'nature', 'science', 'presentation' 等）
        palette: 色板名称（可选，覆盖风格默认色板）

    Returns:
        matplotlib.figure.Figure
    """
    # 按基因数量降序取前 top_n 个通路
    sorted_pathways = sorted(gene_sets.items(), key=lambda x: len(x[1]), reverse=True)
    top_pathways = dict(sorted_pathways[:top_n])
    pathway_names = list(top_pathways.keys())

    if len(pathway_names) < 2:
        fig, ax = plt.subplots(figsize=figsize)
        ax.text(0.5, 0.5, "Not enough pathways for network", ha="center", va="center")
        ax.set_title(title)
        _save_figure(fig, output_file, dpi)
        return fig

    # 计算通路间基因重叠度
    edges = []
    for i, j in combinations(range(len(pathway_names)), 2):
        name_i, name_j = pathway_names[i], pathway_names[j]
        set_i, set_j = top_pathways[name_i], top_pathways[name_j]
        overlap = len(set_i & set_j)
        if overlap >= min_overlap:
            jaccard = overlap / len(set_i | set_j)
            edges.append((i, j, overlap, jaccard))

    # 节点属性
    node_sizes = np.array([len(top_pathways[n]) for n in pathway_names])
    node_sizes_scaled = 200 + 800 * (node_sizes - node_sizes.min()) / max(
        node_sizes.max() - node_sizes.min(), 1
    )

    # 节点颜色
    if results_df is not None and not results_df.empty:
        pathway_col = None
        for col in ("pathway", "Pathway", "term", "Term", "name", "Name"):
            if col in results_df.columns:
                pathway_col = col
                break
        if pathway_col is None:
            pathway_col = results_df.columns[0]

        # 尝试映射 NES 值
        nes_col = None
        for col in ("nes", "NES", "score", "Score"):
            if col in results_df.columns:
                nes_col = col
                break

        if nes_col is not None:
            name_to_nes = dict(zip(results_df[pathway_col], results_df[nes_col]))
            node_colors = [name_to_nes.get(n, 0.0) for n in pathway_names]
            use_cmap = True
            cmap_label = "NES"
        else:
            node_colors = None
            use_cmap = False
            cmap_label = None
    else:
        node_colors = None
        use_cmap = False
        cmap_label = None

    # 使用风格上下文（style 为 None 时默认使用 nature）
    with PlotTheme.context(style or 'nature', palette):
        fig, ax = plt.subplots(figsize=figsize)

        # 获取发散色图（替代硬编码 RdBu_r）
        diverging_cmap = PlotTheme.get_diverging_cmap()
        # 获取色板颜色（替代硬编码 viridis）
        palette_colors = PlotTheme.get_palette(palette, n=len(pathway_names))
        # 获取边的颜色
        edge_color = PlotTheme.get_palette(palette, n=1)[0]

        if HAS_NETWORKX and len(edges) > 0:
            # 使用 networkx 构建图和计算布局
            G = nx.Graph()
            for i, name in enumerate(pathway_names):
                G.add_node(i, label=name)
            for i, j, overlap, jaccard in edges:
                G.add_edge(i, j, weight=jaccard, overlap=overlap)

            # 计算布局
            if layout == "circular":
                pos = nx.circular_layout(G)
            elif layout == "kamada_kawai":
                try:
                    pos = nx.kamada_kawai_layout(G)
                except nx.NetworkXError:
                    pos = nx.spring_layout(G, seed=42)
            else:
                pos = nx.spring_layout(G, seed=42, k=2.0 / np.sqrt(len(G.nodes())))

            # 绘制边 - 使用风格系统获取的灰色
            edge_weights = [G[u][v]["weight"] for u, v in G.edges()]
            edge_widths = [1 + 4 * w for w in edge_weights]
            nx.draw_networkx_edges(
                G, pos, ax=ax, width=edge_widths,
                edge_color=edge_color, alpha=0.6,
            )

            # 绘制节点
            if use_cmap and node_colors is not None:
                vmin = min(node_colors)
                vmax = max(node_colors)
                if vmin == vmax:
                    vmax = vmin + 1
                nc = diverging_cmap((np.array(node_colors) - vmin) / (vmax - vmin))
            else:
                nc = palette_colors

            nx.draw_networkx_nodes(
                G, pos, ax=ax, node_size=node_sizes_scaled,
                node_color=nc, edgecolors="white", linewidths=1.5, alpha=0.9,
            )

            # 绘制标签
            nx.draw_networkx_labels(
                G, pos, ax=ax,
                labels={i: pathway_names[i] for i in G.nodes()},
                font_size=7, font_weight="bold",
            )

            # 添加颜色条
            if use_cmap and node_colors is not None:
                sm = plt.cm.ScalarMappable(
                    cmap=diverging_cmap,
                    norm=plt.Normalize(vmin=vmin, vmax=vmax),
                )
                sm.set_array([])
                cbar = plt.colorbar(sm, ax=ax, shrink=0.5, pad=0.02)
                cbar.set_label(cmap_label, fontsize=10)
        else:
            # 纯 matplotlib 简化布局：圆形排列
            n = len(pathway_names)
            angles = np.linspace(0, 2 * np.pi, n, endpoint=False)
            radius = 1.0
            pos_x = radius * np.cos(angles)
            pos_y = radius * np.sin(angles)

            # 绘制边
            for i, j, overlap, jaccard in edges:
                lw = 1 + 4 * jaccard
                ax.plot(
                    [pos_x[i], pos_x[j]], [pos_y[i], pos_y[j]],
                    color=edge_color, alpha=0.6, linewidth=lw,
                )

            # 绘制节点
            if use_cmap and node_colors is not None:
                vmin = min(node_colors)
                vmax = max(node_colors)
                if vmin == vmax:
                    vmax = vmin + 1
                nc = diverging_cmap((np.array(node_colors) - vmin) / (vmax - vmin))
            else:
                nc = palette_colors

            ax.scatter(pos_x, pos_y, s=node_sizes_scaled, c=nc,
                       edgecolors="white", linewidths=1.5, alpha=0.9, zorder=5)

            # 绘制标签
            for i, name in enumerate(pathway_names):
                ax.annotate(
                    name, (pos_x[i], pos_y[i]),
                    textcoords="offset points", xytext=(0, 10),
                    ha="center", fontsize=7, fontweight="bold",
                )

            # 添加颜色条
            if use_cmap and node_colors is not None:
                sm = plt.cm.ScalarMappable(
                    cmap=diverging_cmap,
                    norm=plt.Normalize(vmin=vmin, vmax=vmax),
                )
                sm.set_array([])
                cbar = plt.colorbar(sm, ax=ax, shrink=0.5, pad=0.02)
                cbar.set_label(cmap_label, fontsize=10)

        ax.set_title(title, fontsize=14, fontweight="bold")
        ax.axis("off")
        plt.tight_layout()
        _save_figure(fig, output_file, dpi)
    return fig


# ---------------------------------------------------------------------------
# plot_upset - UpSet 图（基因集交集）
# ---------------------------------------------------------------------------

def plot_upset(
    gene_sets: Dict[str, Set[str]],
    top_n: int = 20,
    title: str = "Gene Set Overlap",
    output_file: str = None,
    figsize: tuple = (14, 8),
    dpi: int = 300,
    style: Optional[str] = None,
    palette: Optional[str] = None,
) -> matplotlib.figure.Figure:
    """
    UpSet 图（基因集交集可视化）

    取基因数量最多的前 top_n 个基因集，计算所有非空交集大小。
    上方为水平条形图展示交集大小，下方为竖线+圆点矩阵展示交集组成。

    Args:
        gene_sets: 通路名到基因集合的映射 {pathway_name: set(genes)}
        top_n: 展示基因数量最多的前 N 个基因集
        title: 图表标题
        output_file: 输出文件路径（可选）
        figsize: 图表尺寸
        dpi: 输出分辨率
        style: 风格名称（可选）
        palette: 色板名称（可选）

    Returns:
        matplotlib.figure.Figure
    """
    # 按基因数量降序取前 top_n
    sorted_pathways = sorted(gene_sets.items(), key=lambda x: len(x[1]), reverse=True)
    top_pathways = dict(sorted_pathways[:top_n])
    pathway_names = list(top_pathways.keys())
    n_sets = len(pathway_names)

    if n_sets == 0:
        fig, ax = plt.subplots(figsize=figsize)
        ax.text(0.5, 0.5, "No gene sets provided", ha="center", va="center")
        ax.set_title(title)
        _save_figure(fig, output_file, dpi)
        return fig

    # 计算所有非空交集的大小
    intersections = []
    for size in range(2, n_sets + 1):
        for combo in combinations(range(n_sets), size):
            combo_sets = [top_pathways[pathway_names[i]] for i in combo]
            intersect = set.intersection(*combo_sets)
            # 减去所有超集的交集（确保只计算"精确"交集）
            # 这里简化处理：直接计算交集大小
            if len(intersect) > 0:
                intersections.append((combo, len(intersect)))

    # 同时添加单个基因集的大小
    for i in range(n_sets):
        intersections.append(((i,), len(top_pathways[pathway_names[i]])))

    # 按交集大小降序排列
    intersections.sort(key=lambda x: x[1], reverse=True)

    # 限制展示数量（最多 25 个交集）
    max_intersections = 25
    intersections = intersections[:max_intersections]

    if not intersections:
        fig, ax = plt.subplots(figsize=figsize)
        ax.text(0.5, 0.5, "No intersections found", ha="center", va="center")
        ax.set_title(title)
        _save_figure(fig, output_file, dpi)
        return fig

    # 使用风格上下文（style 为 None 时默认使用 nature）
    with PlotTheme.context(style or 'nature', palette):
        # 获取主色
        main_color = PlotTheme.get_palette(palette, n=1)[0]
        # 获取灰色
        gray_color = PlotTheme.get_palette(palette, n=3)[2]
        # 获取浅色
        light_color = PlotTheme.get_palette(palette, n=5)[4]

        # 创建图形
        fig = plt.figure(figsize=figsize)
        # 上方条形图占 40%，下方矩阵占 60%
        gs = fig.add_gridspec(
            2, 1, height_ratios=[1, 1.2], hspace=0.05,
            left=0.15, right=0.95, top=0.92, bottom=0.08,
        )

        ax_bar = fig.add_subplot(gs[0])
        ax_matrix = fig.add_subplot(gs[1])

        n_intersections = len(intersections)
        combo_list = [item[0] for item in intersections]
        size_list = [item[1] for item in intersections]

        # 上方：水平条形图
        y_positions = range(n_intersections)
        ax_bar.barh(list(y_positions), size_list, color=main_color, edgecolor="none", height=0.7)
        ax_bar.set_yticks(list(y_positions))
        ax_bar.set_yticklabels([""] * n_intersections)
        ax_bar.invert_yaxis()
        ax_bar.set_xlabel("Intersection Size", fontsize=10)
        ax_bar.spines["top"].set_visible(False)
        ax_bar.spines["right"].set_visible(False)
        ax_bar.tick_params(labelsize=8)

        # 下方：圆点矩阵
        # 为每个基因集分配 x 位置
        set_x = np.arange(n_sets)

        for row_idx, (combo, size) in enumerate(intersections):
            for col_idx in range(n_sets):
                if col_idx in combo:
                    ax_matrix.plot(
                        set_x[col_idx], row_idx, "o",
                        color=main_color, markersize=8, markeredgecolor="white",
                        markeredgewidth=0.5,
                    )
                else:
                    ax_matrix.plot(
                        set_x[col_idx], row_idx, "o",
                        color=gray_color, markersize=4, markeredgecolor="none",
                    )

        ax_matrix.set_xlim(-0.5, n_sets - 0.5)
        ax_matrix.set_ylim(-0.5, n_intersections - 0.5)
        ax_matrix.invert_yaxis()
        ax_matrix.set_xticks(set_x)

        # 截断长名称
        short_names = []
        for name in pathway_names:
            if len(name) > 20:
                short_names.append(name[:18] + "..")
            else:
                short_names.append(name)
        ax_matrix.set_xticklabels(short_names, rotation=45, ha="right", fontsize=7)
        ax_matrix.set_yticks(range(n_intersections))
        ax_matrix.set_yticklabels([""] * n_intersections)
        ax_matrix.spines["top"].set_visible(False)
        ax_matrix.spines["right"].set_visible(False)
        ax_matrix.spines["left"].set_visible(False)
        ax_matrix.tick_params(axis="y", length=0)

        # 添加竖线连接上下两个面板
        for col_idx in range(n_sets):
            ax_matrix.axvline(x=set_x[col_idx], color=light_color, linewidth=0.5, zorder=0)

        fig.suptitle(title, fontsize=14, fontweight="bold")
        _save_figure(fig, output_file, dpi)
    return fig


# ---------------------------------------------------------------------------
# plot_volcano - 火山图
# ---------------------------------------------------------------------------

def plot_volcano(
    results_df: pd.DataFrame,
    nes_col: str = "nes",
    pvalue_col: str = "pvalue",
    nes_threshold: float = 1.0,
    pvalue_threshold: float = 0.05,
    title: str = "Volcano Plot",
    output_file: str = None,
    figsize: tuple = (10, 8),
    dpi: int = 300,
    style: Optional[str] = None,
    palette: Optional[str] = None,
) -> matplotlib.figure.Figure:
    """
    火山图

    X 轴为 NES（或 log2FC），Y 轴为 -log10(pvalue)。
    显著上调（红色）、显著下调（蓝色）、不显著（灰色）。
    标注 top 显著通路名称。

    Args:
        results_df: 包含 pathway, nes(或 log2fc), pvalue 列的 DataFrame
        nes_col: NES/log2FC 列名
        pvalue_col: pvalue 列名
        nes_threshold: NES 阈值（绝对值）
        pvalue_threshold: pvalue 阈值
        title: 图表标题
        output_file: 输出文件路径（可选）
        figsize: 图表尺寸
        dpi: 输出分辨率
        style: 风格名称（可选）
        palette: 色板名称（可选）

    Returns:
        matplotlib.figure.Figure
    """
    df = results_df.copy()

    # 查找通路名列
    pathway_col = None
    for col in ("pathway", "Pathway", "term", "Term", "name", "Name"):
        if col in df.columns:
            pathway_col = col
            break
    if pathway_col is None:
        pathway_col = df.columns[0]

    # 计算 -log10(pvalue)
    df["neg_log10_p"] = -np.log10(df[pvalue_col].clip(lower=1e-300))

    # 分类：显著上调、显著下调、不显著
    df["category"] = "not_sig"
    df.loc[
        (df[nes_col] >= nes_threshold) & (df[pvalue_col] <= pvalue_threshold),
        "category",
    ] = "up"
    df.loc[
        (df[nes_col] <= -nes_threshold) & (df[pvalue_col] <= pvalue_threshold),
        "category",
    ] = "down"

    # 使用风格上下文（style 为 None 时默认使用 nature）
    with PlotTheme.context(style or 'nature', palette):
        # 颜色映射 - 使用 VOLCANO_COLORS 常量
        color_map = {
            "up": VOLCANO_COLORS["up"],
            "down": VOLCANO_COLORS["down"],
            "not_sig": VOLCANO_COLORS["ns"],
        }
        df["color"] = df["category"].map(color_map)

        fig, ax = plt.subplots(figsize=figsize)

        # 绘制散点
        for cat in ("not_sig", "up", "down"):
            mask = df["category"] == cat
            if mask.sum() == 0:
                continue
            ax.scatter(
                df.loc[mask, nes_col],
                df.loc[mask, "neg_log10_p"],
                c=df.loc[mask, "color"],
                s=40,
                alpha=0.7,
                edgecolors="none",
                label=cat,
            )

        # 阈值线
        ax.axvline(x=nes_threshold, color="gray", linestyle="--", linewidth=0.8, alpha=0.7)
        ax.axvline(x=-nes_threshold, color="gray", linestyle="--", linewidth=0.8, alpha=0.7)
        ax.axhline(
            y=-np.log10(pvalue_threshold), color="gray",
            linestyle="--", linewidth=0.8, alpha=0.7,
        )

        # 标注 top 显著通路
        sig_df = df[df["category"] != "not_sig"]
        if len(sig_df) > 0:
            # 按 -log10(pvalue) * |NES| 排序，取前 5 个
            sig_df = sig_df.copy()
            sig_df["_importance"] = sig_df["neg_log10_p"] * sig_df[nes_col].abs()
            top_sig = sig_df.nlargest(min(5, len(sig_df)), "_importance")

            for _, row in top_sig.iterrows():
                ax.annotate(
                    row[pathway_col],
                    (row[nes_col], row["neg_log10_p"]),
                    textcoords="offset points",
                    xytext=(5, 5),
                    fontsize=7,
                    fontweight="bold",
                    color=row["color"],
                    arrowprops=dict(arrowstyle="-", color="gray", lw=0.5),
                )

        # 图例
        from matplotlib.lines import Line2D
        legend_elements = [
            Line2D([0], [0], marker="o", color="w", markerfacecolor=VOLCANO_COLORS["up"],
                   markersize=8, label=f"Up (NES >= {nes_threshold})"),
            Line2D([0], [0], marker="o", color="w", markerfacecolor=VOLCANO_COLORS["down"],
                   markersize=8, label=f"Down (NES <= -{nes_threshold})"),
            Line2D([0], [0], marker="o", color="w", markerfacecolor=VOLCANO_COLORS["ns"],
                   markersize=8, label="Not Significant"),
        ]
        ax.legend(handles=legend_elements, fontsize=8, loc="upper right")

        ax.set_xlabel(nes_col.upper(), fontsize=12)
        ax.set_ylabel("-log10(P-value)", fontsize=12)
        ax.set_title(title, fontsize=14, fontweight="bold")
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.tick_params(labelsize=9)

        plt.tight_layout()
        _save_figure(fig, output_file, dpi)
    return fig


# ---------------------------------------------------------------------------
# plot_method_comparison - 方法间比较散点图
# ---------------------------------------------------------------------------

def plot_method_comparison(
    results_a: pd.Series,
    results_b: pd.Series,
    method_a_name: str = "Method A",
    method_b_name: str = "Method B",
    title: str = "Method Comparison",
    output_file: str = None,
    figsize: tuple = (8, 8),
    dpi: int = 300,
    style: Optional[str] = None,
    palette: Optional[str] = None,
) -> matplotlib.figure.Figure:
    """
    方法间比较散点图

    取两种方法共有的通路，绘制散点图，添加对角线和回归线，
    计算并标注 Pearson 相关系数。

    Args:
        results_a: 方法 A 的结果 Series（通路名: 得分）
        results_b: 方法 B 的结果 Series（通路名: 得分）
        method_a_name: 方法 A 名称
        method_b_name: 方法 B 名称
        title: 图表标题
        output_file: 输出文件路径（可选）
        figsize: 图表尺寸
        dpi: 输出分辨率
        style: 风格名称（可选）
        palette: 色板名称（可选）

    Returns:
        matplotlib.figure.Figure
    """
    # 取共有通路
    common_pathways = results_a.index.intersection(results_b.index)
    if len(common_pathways) == 0:
        fig, ax = plt.subplots(figsize=figsize)
        ax.text(0.5, 0.5, "No common pathways between methods", ha="center", va="center")
        ax.set_title(title)
        _save_figure(fig, output_file, dpi)
        return fig

    vals_a = results_a.loc[common_pathways].values
    vals_b = results_b.loc[common_pathways].values

    # 计算 Pearson 相关系数（至少需要 2 个数据点）
    if len(common_pathways) >= 2:
        r, p_val = stats.pearsonr(vals_a, vals_b)
    else:
        r, p_val = np.nan, np.nan

    # 使用风格上下文（style 为 None 时默认使用 nature）
    with PlotTheme.context(style or 'nature', palette):
        # 获取颜色
        scatter_color = PlotTheme.get_palette(palette, n=1)[0]
        line_color = PlotTheme.get_palette(palette, n=2)[1]

        fig, ax = plt.subplots(figsize=figsize)

        # 绘制散点
        ax.scatter(vals_a, vals_b, s=50, alpha=0.6, edgecolors="white",
                   linewidths=0.5, color=scatter_color)

        # 对角线 y = x
        all_vals = np.concatenate([vals_a, vals_b])
        val_min, val_max = all_vals.min(), all_vals.max()
        margin = (val_max - val_min) * 0.05
        ax.plot(
            [val_min - margin, val_max + margin],
            [val_min - margin, val_max + margin],
            "k--", linewidth=0.8, alpha=0.5, label="y = x",
        )

        # 回归线
        if len(common_pathways) >= 3:
            slope, intercept, _, _, _ = stats.linregress(vals_a, vals_b)
            x_line = np.linspace(val_min, val_max, 100)
            y_line = slope * x_line + intercept
            ax.plot(x_line, y_line, color=line_color, linewidth=1.5, alpha=0.7,
                    label=f"y = {slope:.2f}x + {intercept:.2f}")

        # 标注 Pearson 相关系数
        if p_val < 0.001:
            p_text = f"p < 0.001"
        else:
            p_text = f"p = {p_val:.3e}"
        ax.text(
            0.05, 0.95,
            f"Pearson r = {r:.3f}\n{p_text}\nn = {len(common_pathways)}",
            transform=ax.transAxes,
            fontsize=10, verticalalignment="top",
            bbox=dict(boxstyle="round,pad=0.3", facecolor="wheat", alpha=0.5),
        )

        ax.set_xlabel(method_a_name, fontsize=12)
        ax.set_ylabel(method_b_name, fontsize=12)
        ax.set_title(title, fontsize=14, fontweight="bold")
        ax.legend(fontsize=9, loc="lower right")
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.tick_params(labelsize=9)

        plt.tight_layout()
        _save_figure(fig, output_file, dpi)
    return fig
