"""
ssGSEA/GSVA 发表级可视化模块

本模块为 ssGSEA 和 GSVA 分析结果提供发表级质量的可视化图表。
所有图表使用 matplotlib + seaborn 绘制，不依赖 R 环境。

主要图表：
- plot_pathway_heatmap: 样本-通路活性聚类热图
- plot_group_comparison: 组间通路活性比较（箱线图/小提琴图/柱状图）
- plot_pathway_dotplot: 通路活性气泡图
- plot_sample_correlation: 样本相关性热图
"""

import logging
import math
from typing import Dict, List, Optional

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy import stats

logger = logging.getLogger(__name__)


def plot_pathway_heatmap(
    scores_df: pd.DataFrame,
    annotation_col: pd.DataFrame = None,
    cluster_rows: bool = True,
    cluster_cols: bool = True,
    cmap: str = "RdBu_r",
    center: float = 0,
    title: str = "Pathway Activity Heatmap",
    output_file: str = None,
    figsize: tuple = None,
    dpi: int = 300,
) -> matplotlib.figure.Figure:
    """
    样本-通路活性聚类热图

    使用 seaborn.clustermap 绘制，支持行/列聚类和分组注释颜色条。

    参数:
        scores_df: 通路活性得分矩阵，行为通路名，列为样本名
        annotation_col: 样本分组注释 DataFrame，列为分组类别，行为样本名（可选）
        cluster_rows: 是否对通路（行）进行聚类
        cluster_cols: 是否对样本（列）进行聚类
        cmap: 颜色映射名称，默认 RdBu_r（红蓝发散色）
        center: 颜色映射中心值，默认 0
        title: 图表标题
        output_file: 输出文件路径（可选，支持 .png/.pdf/.svg）
        figsize: 图表尺寸 (宽, 高)，None 时自动根据数据维度计算
        dpi: 输出分辨率，默认 300

    返回值:
        matplotlib.figure.Figure: 生成的图表对象
    """
    # 自动计算 figsize
    if figsize is None:
        n_pathways, n_samples = scores_df.shape
        # 基础尺寸：每个通路 0.3 行高，每个样本 0.8 列宽
        height = max(6, n_pathways * 0.3 + 2)
        width = max(8, n_samples * 0.8 + 2)
        figsize = (width, height)

    # 构建颜色映射字典（用于分组注释）
    palette_colors = sns.color_palette("Set2", 8).as_hex()

    # 绘制聚类热图
    g = sns.clustermap(
        scores_df,
        method="average",
        metric="euclidean",
        cmap=cmap,
        center=center,
        figsize=figsize,
        row_cluster=cluster_rows,
        col_cluster=cluster_cols,
        linewidths=0.5,
        linecolor="white",
        cbar_kws={"label": "Activity Score", "shrink": 0.8},
        annot=False,
        dendrogram_ratio=(0.15, 0.1),
        colors_ratio=(0.03, 0.03),
        yticklabels=True,
        xticklabels=True,
    )

    # 设置标题
    g.fig.suptitle(title, fontsize=14, fontweight="bold", y=1.02)

    # 旋转 x 轴标签
    g.ax_heatmap.set_xticklabels(
        g.ax_heatmap.get_xticklabels(), rotation=45, ha="right", fontsize=9
    )
    g.ax_heatmap.set_yticklabels(
        g.ax_heatmap.get_yticklabels(), fontsize=9
    )

    # 添加分组注释颜色条
    if annotation_col is not None and not annotation_col.empty:
        # 确保注释 DataFrame 的行与 scores_df 的列对齐
        common_samples = scores_df.columns.intersection(annotation_col.index)
        if len(common_samples) > 0:
            annot_df = annotation_col.loc[common_samples]
            # 获取列聚类后的样本顺序
            if cluster_cols and g.dendrogram_col is not None:
                col_order = g.ax_heatmap.get_xticklabels()
                ordered_samples = [t.get_text() for t in col_order]
                annot_df = annot_df.reindex(ordered_samples).dropna(how="all")

            # 为每个分组类别绘制颜色条
            n_categories = annot_df.shape[1]
            for i, cat_col in enumerate(annot_df.columns):
                unique_groups = annot_df[cat_col].dropna().unique()
                color_map = {
                    g_name: palette_colors[j % len(palette_colors)]
                    for j, g_name in enumerate(unique_groups)
                }
                # 在热图上方添加颜色条
                y_pos = n_categories - i - 1
                for j, sample in enumerate(annot_df.index):
                    val = annot_df.loc[sample, cat_col]
                    if pd.notna(val) and val in color_map:
                        g.ax_col_dendrogram.barh(
                            y_pos, 1, left=j, color=color_map[val], edgecolor="none"
                        )
                # 添加类别标签
                g.ax_col_dendrogram.text(
                    -0.5,
                    y_pos,
                    cat_col,
                    ha="right",
                    va="center",
                    fontsize=8,
                    fontweight="bold",
                )

    # 保存文件
    if output_file is not None:
        g.savefig(output_file, dpi=dpi, bbox_inches="tight", facecolor="white")
        logger.info(f"热图已保存至: {output_file}")

    return g.fig


def plot_group_comparison(
    scores_df: pd.DataFrame,
    groups: Dict[str, List[str]],
    pathways: List[str] = None,
    plot_type: str = "box",
    title: str = "Group Comparison",
    output_file: str = None,
    figsize: tuple = None,
    dpi: int = 300,
) -> matplotlib.figure.Figure:
    """
    组间通路活性比较图

    将宽格式活性矩阵转为长格式，按分组绘制比较图。

    参数:
        scores_df: 通路活性得分矩阵，行为通路名，列为样本名
        groups: 分组字典，键为组名，值为该组样本名列表
        pathways: 要展示的通路列表，None 时取方差最大的前 10 个通路
        plot_type: 图表类型，"box"（箱线图+散点）、"violin"（小提琴图）、"bar"（柱状图）
        title: 图表标题
        output_file: 输出文件路径（可选）
        figsize: 图表尺寸 (宽, 高)，None 时自动计算
        dpi: 输出分辨率，默认 300

    返回值:
        matplotlib.figure.Figure: 生成的图表对象
    """
    if plot_type not in ("box", "violin", "bar"):
        raise ValueError(f"不支持的 plot_type: '{plot_type}'，有效值为 'box', 'violin', 'bar'")

    # 确定展示的通路
    if pathways is None:
        # 取方差最大的前 10 个通路
        pathway_vars = scores_df.var(axis=1)
        pathways = pathway_vars.nlargest(min(10, len(pathway_vars))).index.tolist()

    # 自动计算 figsize
    n_pathways = len(pathways)
    if figsize is None:
        n_cols = min(3, n_pathways)
        n_rows = math.ceil(n_pathways / n_cols)
        figsize = (5 * n_cols, 4 * n_rows)

    # 宽格式转长格式
    df_long = scores_df.loc[pathways].T.melt(
        var_name="Pathway", value_name="Score", ignore_index=False
    )
    df_long.index.name = "Sample"
    df_long = df_long.reset_index()

    # 添加分组列
    sample_to_group = {}
    for group_name, samples in groups.items():
        for s in samples:
            sample_to_group[s] = group_name
    df_long["Group"] = df_long["Sample"].map(sample_to_group)
    df_long = df_long.dropna(subset=["Group"])

    # 创建子图
    n_cols = min(3, n_pathways)
    n_rows = math.ceil(n_pathways / n_cols)
    fig, axes = plt.subplots(n_rows, n_cols, figsize=figsize)
    fig.suptitle(title, fontsize=14, fontweight="bold", y=1.02)

    # 确保 axes 为二维数组
    if n_pathways == 1:
        axes = np.array([[axes]])
    elif n_rows == 1:
        axes = axes.reshape(1, -1)

    # 配色
    group_colors = sns.color_palette("Set2", len(groups))

    for idx, pathway in enumerate(pathways):
        row, col = divmod(idx, n_cols)
        ax = axes[row, col]
        data = df_long[df_long["Pathway"] == pathway]

        if plot_type == "box":
            sns.boxplot(
                data=data, x="Group", y="Score", hue="Group", ax=ax,
                palette=group_colors, width=0.5, linewidth=1,
                fliersize=0, legend=False,
            )
            sns.stripplot(
                data=data, x="Group", y="Score", ax=ax,
                color="black", size=3, alpha=0.6, jitter=True, legend=False,
            )
        elif plot_type == "violin":
            sns.violinplot(
                data=data, x="Group", y="Score", hue="Group", ax=ax,
                palette=group_colors, width=0.7, linewidth=1,
                cut=0, legend=False,
            )
            # 叠加散点
            sns.stripplot(
                data=data, x="Group", y="Score", ax=ax,
                color="black", size=3, alpha=0.5, jitter=True, legend=False,
            )
        elif plot_type == "bar":
            sns.barplot(
                data=data, x="Group", y="Score", hue="Group", ax=ax,
                palette=group_colors, errorbar="sd", capsize=0.1,
                linewidth=1, edgecolor="black", legend=False,
            )

        # 统计显著性标注（t-test）
        group_names = list(groups.keys())
        if len(group_names) == 2:
            g1_data = data[data["Group"] == group_names[0]]["Score"]
            g2_data = data[data["Group"] == group_names[1]]["Score"]
            if len(g1_data) >= 2 and len(g2_data) >= 2:
                stat_val, p_val = stats.ttest_ind(g1_data, g2_data)
                y_max = data["Score"].max()
                y_range = data["Score"].max() - data["Score"].min()
                y_h = y_max + y_range * 0.1
                # 绘制显著性标注线
                ax.plot([0, 0, 1, 1], [y_h, y_h + y_range * 0.05, y_h + y_range * 0.05, y_h],
                        lw=1, c="black")
                if p_val < 0.001:
                    sig_text = "***"
                elif p_val < 0.01:
                    sig_text = "**"
                elif p_val < 0.05:
                    sig_text = "*"
                else:
                    sig_text = "ns"
                ax.text(0.5, y_h + y_range * 0.07, sig_text,
                        ha="center", va="bottom", fontsize=9, fontweight="bold")

        ax.set_title(pathway, fontsize=10, fontweight="bold")
        ax.set_xlabel("")
        ax.set_ylabel("Activity Score" if idx % n_cols == 0 else "")
        ax.tick_params(axis="x", labelsize=8)

    # 隐藏多余的子图
    for idx in range(n_pathways, n_rows * n_cols):
        row, col = divmod(idx, n_cols)
        axes[row, col].set_visible(False)

    plt.tight_layout()

    # 保存文件
    if output_file is not None:
        fig.savefig(output_file, dpi=dpi, bbox_inches="tight", facecolor="white")
        logger.info(f"组间比较图已保存至: {output_file}")

    return fig


def plot_pathway_dotplot(
    scores_df: pd.DataFrame,
    groups: Dict[str, List[str]] = None,
    top_n: int = 20,
    title: str = "Pathway Activity Dot Plot",
    output_file: str = None,
    figsize: tuple = (10, 8),
    dpi: int = 300,
) -> matplotlib.figure.Figure:
    """
    通路活性气泡图

    Y 轴为通路名称，X 轴为活性得分，点大小表示方差，点颜色表示分组。

    参数:
        scores_df: 通路活性得分矩阵，行为通路名，列为样本名
        groups: 分组字典，键为组名，值为该组样本名列表（可选）
        top_n: 展示的通路数量，按均值绝对值排序取前 N 个
        title: 图表标题
        output_file: 输出文件路径（可选）
        figsize: 图表尺寸，默认 (10, 8)
        dpi: 输出分辨率，默认 300

    返回值:
        matplotlib.figure.Figure: 生成的图表对象
    """
    # 计算每个通路的均值和方差
    pathway_stats = pd.DataFrame({
        "mean": scores_df.mean(axis=1),
        "var": scores_df.var(axis=1),
    })

    # 按均值绝对值排序，取前 top_n
    pathway_stats["abs_mean"] = pathway_stats["mean"].abs()
    top_pathways = pathway_stats.nlargest(min(top_n, len(pathway_stats)), "abs_mean").index.tolist()

    # 自动调整 figsize
    height = max(6, len(top_pathways) * 0.4 + 1)
    figsize = (figsize[0], height)

    fig, ax = plt.subplots(figsize=figsize)

    if groups is not None and len(groups) > 0:
        # 有分组：计算每组均值，按分组着色
        group_colors = sns.color_palette("Set2", len(groups))
        plot_data = []
        for group_name, samples in groups.items():
            valid_samples = [s for s in samples if s in scores_df.columns]
            if valid_samples:
                group_mean = scores_df.loc[top_pathways, valid_samples].mean(axis=1)
                for pw in top_pathways:
                    plot_data.append({
                        "Pathway": pw,
                        "Score": group_mean[pw],
                        "Group": group_name,
                        "Variance": pathway_stats.loc[pw, "var"],
                    })
        df_plot = pd.DataFrame(plot_data)

        # 归一化方差用于点大小
        var_min = df_plot["Variance"].min()
        var_max = df_plot["Variance"].max()
        if var_max > var_min:
            df_plot["size"] = 50 + 200 * (df_plot["Variance"] - var_min) / (var_max - var_min)
        else:
            df_plot["size"] = 100

        # 绘制气泡图
        for i, group_name in enumerate(groups.keys()):
            group_data = df_plot[df_plot["Group"] == group_name]
            ax.scatter(
                group_data["Score"],
                group_data["Pathway"],
                s=group_data["size"],
                c=[group_colors[i]],
                label=group_name,
                alpha=0.8,
                edgecolors="black",
                linewidths=0.5,
            )
        ax.legend(title="Group", fontsize=9, title_fontsize=10, loc="lower right")
    else:
        # 无分组：展示所有样本均值
        means = scores_df.loc[top_pathways].mean(axis=1)
        variances = pathway_stats.loc[top_pathways, "var"]

        # 归一化方差用于点大小
        var_min = variances.min()
        var_max = variances.max()
        if var_max > var_min:
            sizes = 50 + 200 * (variances - var_min) / (var_max - var_min)
        else:
            sizes = pd.Series(100, index=top_pathways)

        ax.scatter(
            means, top_pathways, s=sizes,
            c=sns.color_palette("viridis", len(top_pathways)),
            alpha=0.8, edgecolors="black", linewidths=0.5,
        )

    ax.set_xlabel("Activity Score (Mean)", fontsize=12)
    ax.set_ylabel("")
    ax.set_title(title, fontsize=14, fontweight="bold")
    ax.tick_params(axis="y", labelsize=9)
    ax.axvline(x=0, color="gray", linestyle="--", linewidth=0.8, alpha=0.5)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    plt.tight_layout()

    # 保存文件
    if output_file is not None:
        fig.savefig(output_file, dpi=dpi, bbox_inches="tight", facecolor="white")
        logger.info(f"气泡图已保存至: {output_file}")

    return fig


def plot_sample_correlation(
    scores_df: pd.DataFrame,
    method: str = "pearson",
    annotation_col: pd.DataFrame = None,
    title: str = "Sample Correlation",
    output_file: str = None,
    figsize: tuple = (8, 7),
    dpi: int = 300,
) -> matplotlib.figure.Figure:
    """
    样本相关性热图

    计算样本间相关性矩阵并绘制热图，支持分组注释。

    参数:
        scores_df: 通路活性得分矩阵，行为通路名，列为样本名
        method: 相关性计算方法，"pearson" 或 "spearman"
        annotation_col: 样本分组注释 DataFrame（可选）
        title: 图表标题
        output_file: 输出文件路径（可选）
        figsize: 图表尺寸，默认 (8, 7)
        dpi: 输出分辨率，默认 300

    返回值:
        matplotlib.figure.Figure: 生成的图表对象
    """
    if method not in ("pearson", "spearman"):
        raise ValueError(f"不支持的 method: '{method}'，有效值为 'pearson', 'spearman'")

    # 计算相关性矩阵（样本间）
    corr_matrix = scores_df.T.corr(method=method)

    # 自动调整 figsize
    n_samples = len(corr_matrix)
    if n_samples > 10:
        figsize = (max(10, n_samples * 0.6 + 2), max(9, n_samples * 0.6 + 2))

    fig, ax = plt.subplots(figsize=figsize)

    # 绘制热图
    mask = np.zeros_like(corr_matrix, dtype=bool)
    # 不遮掩任何格子（显示完整矩阵）

    hm = sns.heatmap(
        corr_matrix,
        annot=True,
        fmt=".2f",
        cmap="RdBu_r",
        center=1 if method == "pearson" else 1,
        vmin=0 if method == "spearman" else None,
        square=True,
        linewidths=0.5,
        linecolor="white",
        ax=ax,
        annot_kws={"size": 8},
        cbar_kws={"label": f"{method.capitalize()} Correlation", "shrink": 0.8},
    )

    ax.set_title(title, fontsize=14, fontweight="bold")
    ax.tick_params(axis="x", rotation=45, labelsize=9)
    ax.tick_params(axis="y", rotation=0, labelsize=9)

    # 如果有分组注释，在热图上方和左侧添加颜色条
    if annotation_col is not None and not annotation_col.empty:
        palette_colors = sns.color_palette("Set2", 8).as_hex()
        common_samples = corr_matrix.columns.intersection(annotation_col.index)

        if len(common_samples) > 0:
            annot_df = annotation_col.loc[common_samples]

            # 在顶部添加分组颜色条
            n_categories = annot_df.shape[1]
            for i, cat_col in enumerate(annot_df.columns):
                unique_groups = annot_df[cat_col].dropna().unique()
                color_map = {
                    g_name: palette_colors[j % len(palette_colors)]
                    for j, g_name in enumerate(unique_groups)
                }
                for j, sample in enumerate(corr_matrix.columns):
                    if sample in annot_df.index:
                        val = annot_df.loc[sample, cat_col]
                        if pd.notna(val) and val in color_map:
                            ax.add_patch(plt.Rectangle(
                                (j, -n_categories + i),
                                1, 1,
                                color=color_map[val],
                                clip_on=False,
                            ))
                # 添加类别标签
                ax.text(
                    -n_categories * 0.15,
                    -n_categories + i + 0.5,
                    cat_col,
                    ha="right", va="center",
                    fontsize=7, fontweight="bold",
                    transform=ax.get_yaxis_transform(),
                )

    plt.tight_layout()

    # 保存文件
    if output_file is not None:
        fig.savefig(output_file, dpi=dpi, bbox_inches="tight", facecolor="white")
        logger.info(f"相关性热图已保存至: {output_file}")

    return fig
