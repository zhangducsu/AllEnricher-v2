"""
GSEA 可视化配置模块
===================

定义 GSEA 发表级图表的全局配置参数，包括输出格式、分辨率、
配色方案、字体设置等。
"""

from dataclasses import dataclass


@dataclass
class PlotConfig:
    """GSEA 可视化全局配置"""

    figure_format: str = "png"       # 输出格式: png, pdf, svg
    figure_dpi: int = 300            # 输出分辨率（发表级 300dpi）
    color_palette: str = "RdBu_r"    # 配色方案名称
    font_family: str = "sans-serif"  # 字体族
    font_size: int = 10              # 基础字号
    top_n_pathways: int = 20         # 默认展示的通路数量
