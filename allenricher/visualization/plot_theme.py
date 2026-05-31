"""
AllEnricher 可视化风格系统核心模块
=====================================

提供发表级图表的统一风格管理，包括：
- 19组专业配色方案（Paul Tol、Okabe-Ito、期刊风格等）
- 7组图表风格预设（Nature、Science、Cell、Colorblind、演示、OmicShare、Cute）
- 全局/上下文风格切换

Usage:
    # 全局应用风格
    PlotTheme.apply('nature')
    
    # 上下文临时切换
    with PlotTheme.context('presentation'):
        plot_volcano(data)
    
    # 同时保存PNG和PDF
    save_figure_dual(fig, 'output.png')
"""

from __future__ import annotations

import contextlib
import os
from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional, Tuple, Union

import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, ListedColormap


# =============================================================================
# Paul Tol 色板 - 色盲友好、高对比度
# https://personal.sron.nl/~pault/
# =============================================================================

TOL_BRIGHT = [
    "#4477AA", "#66CCEE", "#228833", "#CCBB44", "#EE6677",
    "#AA3377", "#BBBBBB"
]

TOL_HIGH_CONTRAST = [
    "#004488", "#DDAA33", "#BB5566"
]

TOL_VIBRANT = [
    "#0077BB", "#33BBEE", "#009988", "#EE7733", "#CC3311",
    "#EE3377", "#BBBBBB"
]

TOL_MUTED = [
    "#332288", "#88CCEE", "#44AA99", "#117733", "#999933",
    "#DDCC77", "#CC6677", "#882255", "#AA4499", "#DDDDDD"
]

TOL_MEDIUM_CONTRAST = [
    "#6699CC", "#004488", "#EECC66", "#994455", "#997700",
    "#EE99AA"
]

TOL_LIGHT = [
    "#77AADD", "#99DDFF", "#44BB99", "#BBCC77", "#AAAA00",
    "#EEDD88", "#EE8866", "#FFAABB", "#DDDDDD"
]

TOL_SUNSET = [
    "#364B9A", "#4A7BB7", "#6EA6CD", "#98CAE1", "#C2E4EF",
    "#EAECCC", "#FEDA8B", "#FDB366", "#F67E4B", "#DD3D2D",
    "#A50026"
]

TOL_BURGA = [
    "#F7F4F9", "#E7E1EF", "#D4B9DA", "#C994C7", "#DF65B0",
    "#E7298A", "#CE1256", "#980043", "#67001F"
]

TOL_PRGn = [
    "#762A83", "#9970AB", "#C2A5CF", "#E7D4E8", "#F7F7F7",
    "#D9F0D3", "#ACD39E", "#5AAE61", "#1B7837"
]


# =============================================================================
# Okabe-Ito 色板 - 色盲友好标准
# https://jfly.uni-koeln.de/color/
# =============================================================================

OKABE_ITO = [
    "#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442",
    "#0072B2", "#D55E00", "#CC79A7"
]


# =============================================================================
# GO/KEGG 类别色 - 生物信息学专用
# =============================================================================

GO_BP_COLORS = [
    "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",
    "#FFFF33", "#A65628", "#F781BF", "#999999"
]

GO_CC_COLORS = [
    "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854",
    "#FFD92F", "#E5C494", "#B3B3B3"
]

GO_MF_COLORS = [
    "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3",
    "#FDB462", "#B3DE69", "#FCCDE5", "#D9D9D9"
]

KEGG_PATHWAY_COLORS = [
    "#1F77B4", "#FF7F0E", "#2CA02C", "#D62728", "#9467BD",
    "#8C564B", "#E377C2", "#7F7F7F", "#BCBD22", "#17BECF"
]


# =============================================================================
# 火山图专用色
# =============================================================================

VOLCANO_COLORS = {
    "up": "#DC143C",      # 上调 - 深红
    "down": "#4169E1",    # 下调 - 皇家蓝
    "ns": "#808080",      # 不显著 - 灰色
}


# =============================================================================
# 科研期刊风格色板
# =============================================================================

NATURE_COLORS = [
    "#0C5DA5", "#FF9500", "#00B945", "#FF2C00", "#845B97",
    "#474747", "#9E9E9E"
]

SCIENCE_COLORS = [
    "#1F77B4", "#FF7F0E", "#2CA02C", "#D62728", "#9467BD",
    "#8C564B", "#E377C2", "#7F7F7F", "#BCBD22", "#17BECF"
]

CELL_COLORS = [
    "#0072B2", "#D55E00", "#CC79A7", "#F0E442", "#009E73",
    "#56B4E9", "#E69F00", "#000000"
]

LANCET_COLORS = [
    "#00468B", "#ED0000", "#42B540", "#0099B4", "#925E9F",
    "#FDAF91", "#AD002A", "#ADB6B6"
]

NEJM_COLORS = [
    "#BC3C29", "#0072B5", "#E18727", "#20854E", "#7876B1",
    "#6F99AD", "#FFDC91", "#EE4C97"
]

JAMA_COLORS = [
    "#374E55", "#DF8F44", "#00A1D5", "#B24745", "#79AF97",
    "#6A6599", "#80796B"
]


# =============================================================================
# 生物信息学工具风格色板
# =============================================================================

GSEA_COLORS = [
    "#58ACFA", "#BC8F8F", "#FF6347", "#4682B4", "#9ACD32",
    "#DDA0DD", "#F0E68C", "#FF69B4"
]

Cytoscape_COLORS = [
    "#FF9900", "#66CC00", "#0099FF", "#FF0066", "#9900CC",
    "#00CC99", "#FFCC00", "#CC3300"
]

IGV_COLORS = [
    "#0000FF", "#00FF00", "#FF0000", "#00FFFF", "#FF00FF",
    "#FFFF00", "#FFA500", "#800080"
]

TBTOOLS_COLORS = [
    "#FF6B6B", "#4ECDC4", "#45B7D1", "#96CEB4", "#FFEAA7",
    "#DDA0DD", "#98D8C8", "#F7DC6F"
]

OMICSHARE_COLORS = [
    "#FF6B9D", "#C44569", "#F8B500", "#4ECDC4", "#556270",
    "#36D1DC", "#5AB9EA", "#8860D0"
]


# =============================================================================
# 中国风格色板
# =============================================================================

CHINA_STYLE_COLORS = [
    "#C23531", "#2F4554", "#61A0A8", "#D48265", "#91C7AE",
    "#749F83", "#CA8622", "#BDA29A", "#6E7074", "#546570"
]


# =============================================================================
# PALETTES 注册表 - 19组配色方案
# =============================================================================

PALETTES: Dict[str, List[str]] = {
    # 默认 (1组)
    "default": TOL_BRIGHT,
    
    # Paul Tol 系列 (8组)
    "tol_bright": TOL_BRIGHT,
    "tol_high_contrast": TOL_HIGH_CONTRAST,
    "tol_vibrant": TOL_VIBRANT,
    "tol_muted": TOL_MUTED,
    "tol_medium_contrast": TOL_MEDIUM_CONTRAST,
    "tol_light": TOL_LIGHT,
    "tol_sunset": TOL_SUNSET,
    "tol_burga": TOL_BURGA,
    
    # Okabe-Ito - 色盲友好标准 (1组)
    "okabe_ito": OKABE_ITO,
    
    # 科研期刊 (6组)
    "nature": NATURE_COLORS,
    "science": SCIENCE_COLORS,
    "cell": CELL_COLORS,
    "lancet": LANCET_COLORS,
    "nejm": NEJM_COLORS,
    "jama": JAMA_COLORS,
    
    # 生物信息学工具 (2组)
    "gsea": GSEA_COLORS,
    "omicshare": OMICSHARE_COLORS,
    
    # 中国风格 (1组)
    "china_style": CHINA_STYLE_COLORS,

    # GO/KEGG 生物信息学专用 (4组)
    "go_bp": GO_BP_COLORS,
    "go_cc": GO_CC_COLORS,
    "go_mf": GO_MF_COLORS,
    "kegg_pathway": KEGG_PATHWAY_COLORS,
}


# =============================================================================
# 风格预设配置 - 7组风格
# =============================================================================

@dataclass
class StylePreset:
    """
    发表级图表风格预设 - 详细参数配置
    
    包含25+个详细参数，覆盖字体、线条、边框、刻度、网格、图形、背景等
    """
    
    # === 基本信息 ===
    name: str                           # 风格名称
    palette: str                        # 主色板名称
    
    # === 字体配置 (9个参数) ===
    font_family: str = "sans-serif"     # 字体家族 (Arial, Helvetica, Times New Roman, SimHei)
    font_weight: str = "normal"         # 字体粗细 (normal, bold)
    font_size: int = 10                 # 基础字号
    title_size: int = 12                # 图标题字号
    title_weight: str = "bold"          # 图标题粗细
    label_size: int = 10                # 坐标轴标签字号
    tick_label_size: int = 9            # 刻度标签字号
    legend_size: int = 9                # 图例字号
    legend_title_size: int = 10         # 图例标题字号
    
    # === 线条配置 (4个参数) ===
    line_width: float = 1.0             # 数据线宽
    line_style: str = "solid"           # 线条样式 (solid, dashed, dotted)
    marker_size: float = 6.0            # 标记点大小
    marker_edge_width: float = 0.5      # 标记点边缘宽度
    
    # === 边框配置 (4个参数) ===
    axes_line_width: float = 0.8        # 坐标轴线宽
    spine_top: bool = True              # 是否显示上边框
    spine_right: bool = True            # 是否显示右边框
    spine_color: str = "#000000"        # 边框颜色
    
    # === 刻度配置 (5个参数) ===
    tick_direction: str = "out"         # 刻度方向 (in, out, inout)
    tick_major_size: float = 5.0        # 主刻度长度
    tick_major_width: float = 0.8       # 主刻度宽度
    tick_minor_size: float = 3.0        # 次刻度长度
    tick_minor_visible: bool = False    # 是否显示次刻度
    
    # === 网格配置 (6个参数) ===
    grid: bool = False                  # 是否显示网格
    grid_axis: str = "both"             # 网格轴 (both, x, y)
    grid_alpha: float = 0.3             # 网格透明度
    grid_linewidth: float = 0.5         # 网格线宽
    grid_color: str = "#CCCCCC"         # 网格颜色
    grid_linestyle: str = "dashed"      # 网格样式 (solid, dashed, dotted)
    
    # === 图形配置 (6个参数) ===
    figure_dpi: int = 300               # 图形DPI
    figure_format: str = "png"          # 默认输出格式
    savefig_dpi: int = 300              # 保存DPI
    savefig_format: str = "png"         # 保存格式
    savefig_bbox: str = "tight"         # 保存边界框 (tight, standard)
    savefig_pad: float = 0.1            # 保存边距
    
    # === 背景配置 (2个参数) ===
    facecolor: str = "white"            # 图形背景色
    axes_facecolor: str = "white"       # 绘图区背景色
    
    # === 上下文 ===
    context: str = "paper"              # seaborn上下文 (paper, notebook, talk, poster)
    
    # === Matplotlib覆盖 ===
    rc_overrides: Dict[str, Any] = field(default_factory=dict)


# =============================================================================
# 7种风格详细配置
# =============================================================================

PRESETS: Dict[str, StylePreset] = {
    "nature": StylePreset(
        name="nature",
        palette="nature",
        # 字体配置
        font_family="Helvetica",
        font_weight="normal",
        font_size=8,
        title_size=10,
        title_weight="bold",
        label_size=8,
        tick_label_size=7,
        legend_size=7,
        legend_title_size=8,
        # 线条配置
        line_width=1.0,
        line_style="solid",
        marker_size=4.0,
        marker_edge_width=0.5,
        # 边框配置
        axes_line_width=0.5,
        spine_top=False,
        spine_right=False,
        spine_color="#333333",
        # 刻度配置
        tick_direction="in",
        tick_major_size=4.0,
        tick_major_width=0.5,
        tick_minor_size=2.0,
        tick_minor_visible=False,
        # 网格配置
        grid=False,
        grid_axis="both",
        grid_alpha=0.3,
        grid_linewidth=0.5,
        grid_color="#CCCCCC",
        grid_linestyle="dashed",
        # 图形配置
        figure_dpi=300,
        savefig_dpi=300,
        savefig_bbox="tight",
        savefig_pad=0.1,
        # 背景配置
        facecolor="white",
        axes_facecolor="white",
        # 上下文
        context="paper",
    ),
    
    "science": StylePreset(
        name="science",
        palette="science",
        # 字体配置
        font_family="Times New Roman",
        font_weight="normal",
        font_size=9,
        title_size=11,
        title_weight="bold",
        label_size=9,
        tick_label_size=8,
        legend_size=8,
        legend_title_size=9,
        # 线条配置
        line_width=1.5,
        line_style="solid",
        marker_size=6.0,
        marker_edge_width=0.8,
        # 边框配置
        axes_line_width=1.0,
        spine_top=True,
        spine_right=True,
        spine_color="#000000",
        # 刻度配置
        tick_direction="out",
        tick_major_size=5.0,
        tick_major_width=0.8,
        tick_minor_size=3.0,
        tick_minor_visible=True,
        # 网格配置
        grid=False,
        grid_axis="both",
        grid_alpha=0.3,
        grid_linewidth=0.5,
        grid_color="#CCCCCC",
        grid_linestyle="dashed",
        # 图形配置
        figure_dpi=300,
        savefig_dpi=300,
        savefig_bbox="tight",
        savefig_pad=0.1,
        # 背景配置
        facecolor="white",
        axes_facecolor="white",
        # 上下文
        context="paper",
    ),
    
    "cell": StylePreset(
        name="cell",
        palette="cell",
        # 字体配置
        font_family="Arial",
        font_weight="normal",
        font_size=9,
        title_size=11,
        title_weight="bold",
        label_size=9,
        tick_label_size=8,
        legend_size=8,
        legend_title_size=9,
        # 线条配置
        line_width=1.2,
        line_style="solid",
        marker_size=5.0,
        marker_edge_width=0.6,
        # 边框配置
        axes_line_width=0.8,
        spine_top=False,
        spine_right=False,
        spine_color="#333333",
        # 刻度配置
        tick_direction="in",
        tick_major_size=4.0,
        tick_major_width=0.6,
        tick_minor_size=2.0,
        tick_minor_visible=False,
        # 网格配置
        grid=False,
        grid_axis="both",
        grid_alpha=0.3,
        grid_linewidth=0.5,
        grid_color="#CCCCCC",
        grid_linestyle="dashed",
        # 图形配置
        figure_dpi=300,
        savefig_dpi=300,
        savefig_bbox="tight",
        savefig_pad=0.1,
        # 背景配置
        facecolor="white",
        axes_facecolor="white",
        # 上下文
        context="paper",
    ),
    
    "colorblind": StylePreset(
        name="colorblind",
        palette="okabe_ito",
        # 字体配置
        font_family="sans-serif",
        font_weight="bold",
        font_size=10,
        title_size=12,
        title_weight="bold",
        label_size=10,
        tick_label_size=9,
        legend_size=9,
        legend_title_size=10,
        # 线条配置
        line_width=1.0,
        line_style="solid",
        marker_size=6.0,
        marker_edge_width=0.8,
        # 边框配置
        axes_line_width=0.8,
        spine_top=False,
        spine_right=False,
        spine_color="#000000",
        # 刻度配置
        tick_direction="in",
        tick_major_size=5.0,
        tick_major_width=0.6,
        tick_minor_size=3.0,
        tick_minor_visible=False,
        # 网格配置
        grid=False,
        grid_axis="both",
        grid_alpha=0.3,
        grid_linewidth=0.5,
        grid_color="#CCCCCC",
        grid_linestyle="dashed",
        # 图形配置
        figure_dpi=300,
        savefig_dpi=300,
        savefig_bbox="tight",
        savefig_pad=0.1,
        # 背景配置
        facecolor="white",
        axes_facecolor="white",
        # 上下文
        context="paper",
    ),
    
    "presentation": StylePreset(
        name="presentation",
        palette="tol_vibrant",
        # 字体配置
        font_family="sans-serif",
        font_weight="bold",
        font_size=14,
        title_size=16,
        title_weight="bold",
        label_size=14,
        tick_label_size=12,
        legend_size=12,
        legend_title_size=13,
        # 线条配置
        line_width=2.0,
        line_style="solid",
        marker_size=10.0,
        marker_edge_width=1.2,
        # 边框配置
        axes_line_width=1.5,
        spine_top=True,
        spine_right=True,
        spine_color="#000000",
        # 刻度配置
        tick_direction="out",
        tick_major_size=8.0,
        tick_major_width=1.2,
        tick_minor_size=4.0,
        tick_minor_visible=True,
        # 网格配置
        grid=True,
        grid_axis="both",
        grid_alpha=0.5,
        grid_linewidth=1.0,
        grid_color="#EEEEEE",
        grid_linestyle="solid",
        # 图形配置
        figure_dpi=300,
        savefig_dpi=300,
        savefig_bbox="tight",
        savefig_pad=0.2,
        # 背景配置
        facecolor="white",
        axes_facecolor="white",
        # 上下文
        context="talk",
    ),
    
    "omicshare": StylePreset(
        name="omicshare",
        palette="omicshare",
        # 字体配置
        font_family="Microsoft YaHei",
        font_weight="normal",
        font_size=10,
        title_size=12,
        title_weight="normal",
        label_size=10,
        tick_label_size=9,
        legend_size=9,
        legend_title_size=10,
        # 线条配置
        line_width=1.2,
        line_style="solid",
        marker_size=6.0,
        marker_edge_width=0.6,
        # 边框配置
        axes_line_width=0.8,
        spine_top=True,
        spine_right=True,
        spine_color="#333333",
        # 刻度配置
        tick_direction="in",
        tick_major_size=5.0,
        tick_major_width=0.6,
        tick_minor_size=3.0,
        tick_minor_visible=False,
        # 网格配置
        grid=False,
        grid_axis="both",
        grid_alpha=0.3,
        grid_linewidth=0.5,
        grid_color="#CCCCCC",
        grid_linestyle="dashed",
        # 图形配置
        figure_dpi=300,
        savefig_dpi=300,
        savefig_bbox="tight",
        savefig_pad=0.1,
        # 背景配置
        facecolor="white",
        axes_facecolor="white",
        # 上下文
        context="paper",
    ),
}


# =============================================================================
# 便捷函数：同时保存PNG和PDF
# =============================================================================

def save_figure_dual(
    fig,
    output_path: str,
    dpi: int = 300,
    bbox_inches: str = "tight",
    facecolor: Optional[str] = None,
    **kwargs
) -> Tuple[str, str]:
    """
    同时保存PNG和PDF格式
    
    自动同时生成PNG（位图）和PDF（矢量图）两种格式，
    满足发表级图表的输出需求。
    
    Args:
        fig: matplotlib Figure对象
        output_path: 输出路径（可以是.png、.pdf或不含扩展名）
        dpi: PNG图像的DPI，默认300
        bbox_inches: 边界框设置，默认'tight'
        facecolor: 背景色，None则使用fig.get_facecolor()
        **kwargs: 其他传递给savefig的参数
        
    Returns:
        Tuple[str, str]: (png_path, pdf_path) 两个输出文件路径
        
    Example:
        >>> fig, ax = plt.subplots()
        >>> ax.plot([1, 2, 3], [1, 4, 2])
        >>> png_path, pdf_path = save_figure_dual(fig, 'output.png')
        >>> print(f"Saved: {png_path}, {pdf_path}")
    """
    # 获取基础路径（去除扩展名）
    base_path = output_path
    for ext in ['.png', '.pdf', '.jpg', '.jpeg', '.svg', '.eps']:
        if base_path.lower().endswith(ext):
            base_path = base_path[:-len(ext)]
            break
    
    png_path = base_path + '.png'
    pdf_path = base_path + '.pdf'
    
    # 获取背景色
    if facecolor is None:
        facecolor = fig.get_facecolor()
    
    # 保存PNG（位图）
    fig.savefig(
        png_path,
        format='png',
        dpi=dpi,
        bbox_inches=bbox_inches,
        facecolor=facecolor,
        **kwargs
    )
    
    # 保存PDF（矢量图）
    fig.savefig(
        pdf_path,
        format='pdf',
        bbox_inches=bbox_inches,
        facecolor=facecolor,
        **kwargs
    )
    
    return png_path, pdf_path


# =============================================================================
# PlotTheme 主类
# =============================================================================

class PlotTheme:
    """
    AllEnricher 可视化风格管理器
    
    提供全局风格应用、上下文临时切换、色板获取等功能。
    
    Examples:
        >>> # 查看可用资源
        >>> PlotTheme.available_styles()
        ['nature', 'science', 'cell', 'colorblind', 'presentation', 'omicshare']
        >>> PlotTheme.available_palettes()
        ['default', 'tol_bright', 'tol_high_contrast', ...]
        
        >>> # 全局应用风格
        >>> PlotTheme.apply('nature')
        
        >>> # 临时切换风格
        >>> with PlotTheme.context('presentation'):
        ...     plot_volcano(data)
        
        >>> # 获取当前激活配置
        >>> theme = PlotTheme.get_active()
        >>> colors = PlotTheme.get_palette('tol_vibrant', n=5)
    """
    
    _active_preset: Optional[str] = None
    _active_config: Optional[StylePreset] = None
    _original_rc: Optional[Dict[str, Any]] = None
    
    @classmethod
    def available_styles(cls) -> List[str]:
        """返回所有可用的风格预设名称"""
        return list(PRESETS.keys())
    
    @classmethod
    def available_palettes(cls) -> List[str]:
        """返回所有可用的色板名称"""
        return list(PALETTES.keys())
    
    @classmethod
    def apply(cls, style: str) -> None:
        """
        全局应用指定风格
        
        Args:
            style: 风格名称，必须是 available_styles() 返回的名称之一
            
        Raises:
            ValueError: 如果风格名称不存在
        """
        if style not in PRESETS:
            available = ", ".join(cls.available_styles())
            raise ValueError(f"Unknown style '{style}'. Available: {available}")
        
        preset = PRESETS[style]
        cls._active_preset = style
        cls._active_config = preset
        
        # 保存原始配置（首次应用时）
        if cls._original_rc is None:
            cls._original_rc = plt.rcParams.copy()
        
        # 构建 rcParams - 完整映射所有25+参数
        rc_params = {
            # === 字体配置 ===
            "font.family": preset.font_family,
            "font.sans-serif": ["Helvetica", "Arial", "DejaVu Sans", "Liberation Sans",
                                 "Tahoma", "Verdana", "sans-serif"],
            "font.weight": preset.font_weight,
            "font.size": preset.font_size,
            "axes.titlesize": preset.title_size,
            "axes.titleweight": preset.title_weight,
            "axes.labelsize": preset.label_size,
            "xtick.labelsize": preset.tick_label_size,
            "ytick.labelsize": preset.tick_label_size,
            "legend.fontsize": preset.legend_size,
            "legend.title_fontsize": preset.legend_title_size,
            
            # === 线条配置 ===
            "lines.linewidth": preset.line_width,
            "lines.linestyle": preset.line_style,
            "lines.markersize": preset.marker_size,
            "lines.markeredgewidth": preset.marker_edge_width,
            
            # === 边框配置 ===
            "axes.linewidth": preset.axes_line_width,
            "axes.spines.top": preset.spine_top,
            "axes.spines.right": preset.spine_right,
            "axes.edgecolor": preset.spine_color,
            
            # === 刻度配置 ===
            "xtick.direction": preset.tick_direction,
            "ytick.direction": preset.tick_direction,
            "xtick.major.size": preset.tick_major_size,
            "ytick.major.size": preset.tick_major_size,
            "xtick.major.width": preset.tick_major_width,
            "ytick.major.width": preset.tick_major_width,
            "xtick.minor.size": preset.tick_minor_size,
            "ytick.minor.size": preset.tick_minor_size,
            "xtick.minor.visible": preset.tick_minor_visible,
            "ytick.minor.visible": preset.tick_minor_visible,
            
            # === 网格配置 ===
            "axes.grid": preset.grid,
            "axes.grid.axis": preset.grid_axis,
            "grid.alpha": preset.grid_alpha,
            "grid.linewidth": preset.grid_linewidth,
            "grid.color": preset.grid_color,
            "grid.linestyle": preset.grid_linestyle,
            
            # === 图形配置 ===
            "figure.dpi": preset.figure_dpi,
            "savefig.dpi": preset.savefig_dpi,
            "savefig.bbox": preset.savefig_bbox,
            "savefig.pad_inches": preset.savefig_pad,
            
            # === 背景配置 ===
            "figure.facecolor": preset.facecolor,
            "axes.facecolor": preset.axes_facecolor,
            
            # === 颜色循环 ===
            "axes.prop_cycle": plt.cycler(color=PALETTES[preset.palette]),
        }
        
        # 应用覆盖参数
        rc_params.update(preset.rc_overrides)
        
        # 应用到 matplotlib
        plt.rcParams.update(rc_params)
    
    @classmethod
    @contextlib.contextmanager
    def context(cls, style: str, palette: Optional[str] = None):
        """
        上下文管理器，临时切换风格
        
        Args:
            style: 风格名称
            palette: 可选色板名称（覆盖风格默认色板）
            
        Example:
            >>> with PlotTheme.context('presentation'):
            ...     fig, ax = plt.subplots()
            ...     ax.plot([1, 2, 3], [1, 4, 2])
            ...     plt.show()
        """
        # 保存当前状态
        previous_preset = cls._active_preset
        previous_rc = plt.rcParams.copy()
        
        try:
            cls.apply(style)
            # 如果指定了palette，覆盖颜色循环
            if palette is not None and palette in PALETTES:
                plt.rcParams["axes.prop_cycle"] = plt.cycler(color=PALETTES[palette])
            yield cls
        finally:
            # 恢复之前的状态
            cls._active_preset = previous_preset
            if previous_preset:
                cls._active_config = PRESETS.get(previous_preset)
            else:
                cls._active_config = None
            plt.rcParams.update(previous_rc)
    
    @classmethod
    def get_active(cls) -> Optional[StylePreset]:
        """获取当前激活的风格配置"""
        return cls._active_config
    
    @classmethod
    def get_palette(cls, name: Optional[str] = None, n: Optional[int] = None, palette: Optional[str] = None) -> List[str]:
        """
        获取指定色板的颜色列表
        
        Args:
            name: 色板名称，None则使用当前激活风格的色板（向后兼容）
            n: 需要的颜色数量，None则返回全部
            palette: 色板名称（新参数，与name功能相同，优先级更高）
            
        Returns:
            颜色列表（十六进制字符串）
        """
        # 优先使用palette参数
        if palette is not None:
            name = palette
        
        if name is None:
            if cls._active_config is None:
                name = "default"
            else:
                name = cls._active_config.palette
        
        if name not in PALETTES:
            available = ", ".join(cls.available_palettes())
            raise ValueError(f"Unknown palette '{name}'. Available: {available}")
        
        colors = PALETTES[name]
        
        if n is None:
            return colors
        
        # 循环扩展或截断
        if n <= len(colors):
            return colors[:n]
        else:
            # 循环使用颜色
            result = []
            for i in range(n):
                result.append(colors[i % len(colors)])
            return result
    
    @classmethod
    def get_sequential_cmap(
        cls,
        name: Optional[str] = None,
        colors: Optional[List[str]] = None
    ) -> LinearSegmentedColormap:
        """
        获取连续色图（Sequential Colormap）
        
        Args:
            name: 使用预定义色板名称生成连续色图
            colors: 直接指定颜色列表（优先级高于name）
            
        Returns:
            matplotlib LinearSegmentedColormap
        """
        if colors is None:
            if name is None:
                name = cls._active_config.palette if cls._active_config else "tol_sunset"
            colors = PALETTES.get(name, TOL_SUNSET)
        
        return LinearSegmentedColormap.from_list("sequential", colors)
    
    @staticmethod
    def get_continuous_cmap(palette: str = "tol_sunset") -> LinearSegmentedColormap:
        """
        获取连续色图（Continuous Colormap）
        
        使用PlotTheme.get_palette获取256色，创建LinearSegmentedColormap。
        适用于需要连续渐变色的场景，如热图、密度图等。
        
        Args:
            palette: 色板名称，默认为"tol_sunset"
            
        Returns:
            matplotlib LinearSegmentedColormap，包含256个连续颜色
            
        Example:
            >>> cmap = PlotTheme.get_continuous_cmap("tol_sunset")
            >>> plt.imshow(data, cmap=cmap)
        """
        colors = PlotTheme.get_palette(palette, n=256)
        return LinearSegmentedColormap.from_list(f"{palette}_continuous", colors)
    
    @classmethod
    def get_diverging_cmap(
        cls,
        name: Optional[str] = None,
        colors: Optional[List[str]] = None
    ) -> LinearSegmentedColormap:
        """
        获取发散色图（Diverging Colormap）
        
        Args:
            name: 使用预定义发散色板名称
            colors: 直接指定颜色列表（优先级高于name）
            
        Returns:
            matplotlib LinearSegmentedColormap
        """
        diverging_palettes = {
            "tol_prgn": TOL_PRGn,
            "tol_sunset": TOL_SUNSET,
            "rdbu": ["#2166AC", "#F7F7F7", "#B2182B"],
            "brbg": ["#543005", "#F7F7F7", "#003C30"],
        }
        
        if colors is None:
            if name and name in diverging_palettes:
                colors = diverging_palettes[name]
            else:
                colors = TOL_PRGn
        
        return LinearSegmentedColormap.from_list("diverging", colors)
    
    @classmethod
    def get_category_colors(
        cls,
        categories: List[str],
        palette: Optional[str] = None
    ) -> Dict[str, str]:
        """
        为类别分配颜色
        
        Args:
            categories: 类别列表
            palette: 使用的色板名称
            
        Returns:
            类别到颜色的映射字典
        """
        colors = cls.get_palette(palette, n=len(categories))
        return dict(zip(categories, colors))
    
    @classmethod
    def reset(cls) -> None:
        """重置为 matplotlib 默认配置"""
        cls._active_preset = None
        cls._active_config = None
        if cls._original_rc is not None:
            plt.rcParams.update(cls._original_rc)
            cls._original_rc = None
        else:
            plt.rcdefaults()


# =============================================================================
# 便捷函数
# =============================================================================

def set_theme(style: str) -> None:
    """便捷函数：全局应用风格"""
    PlotTheme.apply(style)


def get_colors(palette: Optional[str] = None, n: Optional[int] = None) -> List[str]:
    """便捷函数：获取色板颜色"""
    return PlotTheme.get_palette(palette, n)


# =============================================================================
# 模块初始化
# =============================================================================

# 默认应用 'default' 色板（不改变其他设置）
plt.rcParams["axes.prop_cycle"] = plt.cycler(color=PALETTES["default"])
