"""
AllEnricher 气泡图可视化模块
==============================

提供富集分析结果的气泡图可视化功能。

气泡图特点：
- X轴：RichFactor（富集因子）
- Y轴：Term_Name（通路/GO名称）
- 点大小：GeneCount（基因数量）
- 颜色：-log10(Qvalue)（显著性）

Usage:
    >>> from allenricher.visualization.bubble import plot_bubble
    >>> fig = plot_bubble(df, top_n=20, style='nature')
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg", force=True)
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from typing import Optional, Union, Tuple, Any

from allenricher.visualization.plot_theme import PlotTheme
from allenricher.visualization.color_config import ColorConfig
from allenricher.visualization.plot_utils import clean_pathway_label, term_figure_size


def plot_bubble(
    data: pd.DataFrame,
    top_n: Optional[int] = 20,
    style: Optional[str] = None,
    palette: Optional[str] = None,
    figsize: Optional[Tuple[float, float]] = None,
    title: Optional[str] = None,
    xlabel: str = "Rich Factor",
    ylabel: str = "",
    colorbar_label: str = "-log10(Qvalue)",
    ax: Optional[plt.Axes] = None,
    **kwargs
) -> Union[Figure, Tuple[Figure, plt.Axes]]:
    """
    绘制富集分析气泡图
    
    气泡图展示富集分析结果，通过气泡大小表示基因数量，
    颜色表示显著性水平（-log10(Qvalue)）。
    
    Args:
        data: 富集分析结果DataFrame，必须包含以下列：
            - Term_Name: 通路/GO名称
            - RichFactor: 富集因子
            - GeneCount: 基因数量
            - Qvalue: 校正后的p值
        top_n: 显示前N个条目，None表示显示全部
        style: 图表风格（nature, science, presentation等）
        palette: 色板名称，None则使用当前风格默认色板
        figsize: 图形尺寸 (width, height)，None则自动计算
        title: 图表标题
        xlabel: X轴标签
        ylabel: Y轴标签
        colorbar_label: 色条标签
        ax: 可选的Axes对象，None则创建新Figure
        **kwargs: 额外的scatter参数
        
    Returns:
        如果ax为None，返回Figure；否则返回(ax, scatter)
        
    Examples:
        >>> import pandas as pd
        >>> df = pd.DataFrame({
        ...     'Term_Name': ['GO:001', 'GO:002', 'GO:003'],
        ...     'RichFactor': [2.5, 3.0, 1.8],
        ...     'GeneCount': [10, 15, 8],
        ...     'Qvalue': [0.001, 0.0001, 0.01]
        ... })
        >>> fig = plot_bubble(df, top_n=20, style='nature')
        
        >>> # 使用自定义风格
        >>> fig = plot_bubble(df, style='presentation', palette='tol_sunset')
    """
    # 验证输入数据
    required_cols = ['Term_Name', 'RichFactor', 'GeneCount', 'Qvalue']
    missing_cols = [col for col in required_cols if col not in data.columns]
    if missing_cols:
        raise ValueError(f"数据缺少必需列: {missing_cols}")
    
    # 数据预处理
    df = data.copy()
    
    # 计算 -log10(Qvalue)
    df['neg_log10_qvalue'] = -np.log10(df['Qvalue'].replace(0, np.nan))
    df['neg_log10_qvalue'] = df['neg_log10_qvalue'].fillna(
        -np.log10(df['Qvalue'][df['Qvalue'] > 0].min() * 0.1)
    )
    
    # 按Qvalue排序，选择top_n
    df = df.sort_values('Qvalue', ascending=True).reset_index(drop=True)
    if top_n is not None:
        df = df.head(top_n)
    
    # 根据palette生成颜色映射
    n_bars = len(df)
    colors = PlotTheme.get_palette(palette, n=n_bars)
    
    # 记录是否为自动创建的图形（用于决定返回值）
    ax_was_none = ax is None
    
    # 应用风格
    with PlotTheme.context(style or 'nature'):
        # 创建图形
        if ax_was_none:
            if figsize is None:
                # 根据条目数量自动计算高度
                figsize = term_figure_size(len(df), width=8.0, min_height=2.8, row_height=0.44)
            fig, ax = plt.subplots(figsize=figsize)
        else:
            fig = ax.figure
        
        # 获取色图 - 三色渐变：取色板首尾两个不连续颜色实现最大差异
        from matplotlib.colors import LinearSegmentedColormap
        palette_colors = PlotTheme.get_palette(palette, n=6)
        cmap = LinearSegmentedColormap.from_list(
            'bubble_cmap',
            [palette_colors[0], '#FFFFFF', palette_colors[-1]]
        )
        
        # 计算气泡大小（基于GeneCount）
        # 使用面积映射，最小50，最大300
        min_count = df['GeneCount'].min()
        max_count = df['GeneCount'].max()
        if max_count > min_count:
            sizes = 50 + (df['GeneCount'] - min_count) / (max_count - min_count) * 250
        else:
            sizes = np.full(len(df), 100)
        
        # 绘制气泡
        scatter = ax.scatter(
            df['RichFactor'],
            range(len(df)),
            s=sizes,
            c=df['neg_log10_qvalue'],
            cmap=cmap,
            alpha=0.8,
            edgecolors='white',
            linewidths=0.5,
            **kwargs
        )
        
        # 设置Y轴标签
        ax.set_yticks(range(len(df)))
        ax.set_yticklabels([clean_pathway_label(term) for term in df['Term_Name']], fontsize=9)
        
        # 反转Y轴，使最显著的显示在顶部
        ax.set_ylim(len(df) - 0.5, -0.5)
        
        # 设置轴标签
        ax.set_xlabel(xlabel, fontsize=10)
        ax.set_ylabel(ylabel, fontsize=10)
        
        # 设置标题
        if title:
            ax.set_title(title, fontsize=12, fontweight='bold')
        
        # 添加网格线
        ax.grid(True, axis='x', alpha=0.3, linestyle='--')
        ax.set_axisbelow(True)
        
        # 添加色条
        color_has_range = df['neg_log10_qvalue'].nunique(dropna=True) > 1
        if color_has_range:
            norm = plt.Normalize(df['neg_log10_qvalue'].min(), df['neg_log10_qvalue'].max())
            sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
            sm.set_array([])
            cbar = plt.colorbar(sm, ax=ax, pad=0.02, shrink=0.45)
            cbar.set_label(colorbar_label, fontsize=9)
            cbar.ax.tick_params(labelsize=8)
        
        # 创建气泡大小数值的图例
        min_gc = df['GeneCount'].min()
        max_gc = df['GeneCount'].max()
        if max_gc > min_gc and len(df) > 1:
            legend_gc = np.linspace(min_gc, max_gc, 4).astype(int)
            legend_sizes = 50 + (legend_gc - min_gc) / (max_gc - min_gc) * 250
            for gc, sz in zip(legend_gc, legend_sizes):
                ax.scatter([], [], s=sz, c='#888888',
                           edgecolors='#888888', linewidth=0.5, label=str(gc))

            handles, labels = ax.get_legend_handles_labels()
            by_label = dict(zip(labels, handles))
            ax.legend(by_label.values(), by_label.keys(), scatterpoints=1,
                      frameon=False, labelspacing=1,
                      title='Gene Count', title_fontsize=8, fontsize=8)
        elif len(df) == 1:
            row = df.iloc[0]
            rich_factor = float(row['RichFactor'])
            x_range = max(abs(rich_factor), 1e-6)
            ax.set_xlim(0, rich_factor + x_range * 0.45)
            ax.text(
                rich_factor + x_range * 0.08,
                0,
                f"Gene count: {int(row['GeneCount'])}  FDR: {row['Qvalue']:.2g}  Rich factor: {rich_factor:.2f}",
                va='center',
                ha='left',
                fontsize=8.5,
                color='#333333',
            )
        
        # 调整布局
        plt.tight_layout()
    
    if ax_was_none:
        return fig
    return fig, ax


def plot_bubble_comparison(
    data_dict: dict,
    style: Optional[str] = None,
    palette: Optional[str] = None,
    figsize: Optional[Tuple[float, float]] = None,
    title: Optional[str] = None,
    sharey: bool = True,
    **kwargs
) -> Figure:
    """
    绘制多组富集分析气泡图对比
    
    Args:
        data_dict: 字典，键为分组名称，值为DataFrame
        style: 图表风格
        palette: 色板名称
        figsize: 图形尺寸
        title: 总标题
        sharey: 是否共享Y轴
        **kwargs: 传递给plot_bubble的额外参数
        
    Returns:
        Figure对象
    """
    n_groups = len(data_dict)
    
    if figsize is None:
        figsize = (6 * n_groups, 10)
    
    # 应用风格
    with PlotTheme.context(style or 'nature'):
        fig, axes = plt.subplots(1, n_groups, figsize=figsize, sharey=sharey)
        if n_groups == 1:
            axes = [axes]
        
        for ax, (group_name, data) in zip(axes, data_dict.items()):
            plot_bubble(
                data,
                style=None,  # 风格已在上面应用
                palette=palette,
                title=group_name,
                ax=ax,
                **kwargs
            )
        
        if title:
            fig.suptitle(title, fontsize=14, fontweight='bold', y=1.02)
        
        plt.tight_layout()
    return fig
