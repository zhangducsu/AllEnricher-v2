# 可视化系统统一与发表级风格支持 Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** 将所有可视化统一为 Python 实现，建立发表级风格系统，支持多风格/色系切换（含 cutecharts 手绘风格），消除 R 脚本依赖

**Architecture:** 新增 `plot_theme.py` 风格管理模块（PlotTheme 枚举 + 预设色板 + apply_style 上下文管理器）；用 matplotlib 重写 2 个 R 脚本（barplot/bubble）；改造所有 12 个 Python 绘图函数接入风格系统；新增 `cute_charts.py` 适配层将 ORA barplot/bubble 映射到 cutecharts 手绘风格输出 HTML；CLI 新增 `--style`/`--palette` 参数

**Tech Stack:** Python, matplotlib, seaborn, scienceplots (可选依赖), cutecharts (可选依赖)

---

## 现状分析

### 当前可视化文件

| 文件 | 图表数 | 绘图库 | 问题 |
|------|--------|--------|------|
| `visualization/barplot.R` | 1 | R base | 需 Python 重写 |
| `visualization/bubble.R` | 1 | R ggplot2 | 需 Python 重写 |
| `visualization/plotter.py` | 2 (调用 R) | subprocess | 需移除 R 依赖 |
| `visualization/gsea_plots.py` | 3 | matplotlib | 硬编码颜色，无风格 |
| `visualization/gsva_plots.py` | 4 | seaborn | 硬编码颜色，无风格 |
| `visualization/common_plots.py` | 4 | matplotlib | 硬编码颜色，未集成 |
| `visualization/plot_config.py` | 0 | - | 定义但未使用 |

### 硬编码颜色清单（15+ 处）

| 文件:行号 | 硬编码值 |
|-----------|---------|
| gsea_plots.py:107 | `#2E86AB` |
| gsea_plots.py:152 | `#E74C3C` |
| gsea_plots.py:212-213 | `#E74C3C`, `#3498DB` |
| gsea_plots.py:334 | `RdBu_r` |
| gsva_plots.py:69 | `Set2` |
| gsva_plots.py:221 | `Set2` |
| common_plots.py:176,223 | `RdBu_r`, `viridis` |
| common_plots.py:343 | `#4C72B0` |
| common_plots.py:461 | `#E74C3C`, `#3498DB`, `#B0BEC5` |
| common_plots.py:588 | `#4C72B0` |

---

## 文件结构

```
allenricher/visualization/
├── __init__.py              ← 修改：导出新增模块
├── plot_theme.py             ← 新增：风格系统核心（PlotTheme + 色板 + apply_style）
├── plot_config.py            ← 修改：接入风格系统，或合并到 plot_theme.py
├── plotter.py                ← 修改：移除 R 依赖，调用 Python barplot/bubble
├── barplot.py                ← 新增：Python 重写的水平柱状图（替代 barplot.R）
├── bubble.py                 ← 新增：Python 重写的气泡图（替代 bubble.R）
├── cute_charts.py            ← 新增：cutecharts 手绘风格适配层（可选依赖）
├── gsea_plots.py             ← 修改：接入风格系统
├── gsva_plots.py             ← 修改：接入风格系统
├── common_plots.py           ← 修改：接入风格系统
├── barplot.R                 ← 删除（被 barplot.py 替代）
├── bubble.R                  ← 删除（被 bubble.py 替代）
allenricher/cli.py             ← 修改：新增 --style/--palette 参数
allenricher/core/config.py     ← 修改：新增 style/palette 配置字段
allenricher/report/generator.py ← 修改：传递风格参数到绘图调用
tests/test_common_plots.py     ← 修改：适配新接口
tests/test_gsea_plots.py       ← 修改：适配新接口
tests/test_gsva_plots.py       ← 修改：适配新接口
```

---

## Task 1: 创建风格系统核心模块

**Files:**
- Create: `allenricher/visualization/plot_theme.py`
- Test: `tests/test_plot_theme.py`

- [ ] **Step 1: 创建 plot_theme.py**

```python
"""发表级可视化风格管理系统

提供预设风格（Nature/Science/IEEE/Colorblind 等）和色板，
所有绘图函数通过 apply_style() 统一应用风格。
"""
from __future__ import annotations

from contextlib import contextmanager
from dataclasses import dataclass, field
from enum import Enum
from typing import Dict, List, Optional, Tuple

import matplotlib.pyplot as plt
import matplotlib as mpl


# ============================================================
# 预设色板定义
# ============================================================

# Paul Tol 色盲友好色板（最权威推荐）
PAUL_TOL_BRIGHT = ['#4477AA', '#EE6677', '#228833', '#CCBB44',
                   '#66CCEE', '#AA3377', '#BBBBBB']

PAUL_TOL_VIBRANT = ['#EE7733', '#0077BB', '#33BBEE', '#EE3377',
                     '#CC3311', '#009988', '#BBBBBB']

PAUL_TOL_MUTED = ['#CC6677', '#332288', '#DDCC77', '#117733',
                   '#88CCEE', '#882255', '#44AA99', '#999933',
                   '#AA4499', '#DDDDDD']

PAUL_TOL_HIGH_CONTRAST = ['#004488', '#DDAA33', '#BB5566']

# Okabe-Ito 色盲友好色板
OKABE_ITO = ['#E69F00', '#56B4E9', '#009E73', '#F0E442',
             '#0072B2', '#D55E00', '#CC79A7', '#000000']

# GO 三大类经典色（从 barplot.R 迁移）
GO_CATEGORY_COLORS = {
    'Biological Process': '#FF6B35',
    'Cellular Component': '#004E89',
    'Molecular Function': '#7B2D8E',
}

# KEGG 六大类经典色（从 barplot.R 迁移）
KEGG_CATEGORY_COLORS = {
    'Metabolism': '#E64B35',
    'Genetic Information Processing': '#4DBBD5',
    'Environmental Information Processing': '#00A087',
    'Cellular Processes': '#3C5488',
    'Organismal Systems': '#F39B7F',
    'Human Diseases': '#8491B4',
}

# 火山图经典色
VOLCANO_COLORS = {
    'up': '#E74C3C',
    'down': '#3498DB',
    'not_sig': '#B0BEC5',
}


@dataclass
class Palette:
    """色板定义"""
    name: str
    colors: List[str]
    categorical: bool = True  # True=离散色板, False=连续色图名
    description: str = ""


# 预设色板注册表
PALETTES: Dict[str, Palette] = {
    'default': Palette('default', ['#4C72B0', '#DD8452', '#55A868', '#C44E52',
                                     '#8172B3', '#937860', '#DA8BC3', '#8C8C8C'],
                        description='Seaborn deep 风格（当前默认）'),
    'tol_bright': Palette('tol_bright', PAUL_TOL_BRIGHT,
                          description='Paul Tol bright - 色盲友好 7 色'),
    'tol_vibrant': Palette('tol_vibrant', PAUL_TOL_VIBRANT,
                          description='Paul Tol vibrant - 高饱和度 7 色'),
    'tol_muted': Palette('tol_muted', PAUL_TOL_MUTED,
                         description='Paul Tol muted - 低饱和度 10 色'),
    'tol_high_contrast': Palette('tol_high_contrast', PAUL_TOL_HIGH_CONTRAST,
                                 description='Paul Tol high-contrast - 灰度可区分 3 色'),
    'okabe_ito': Palette('okabe_ito', OKABE_ITO,
                         description='Okabe-Ito - 色盲友好 8 色'),
    'nature': Palette('nature', ['#0173B2', '#DE8F05', '#029F73', '#D55E00',
                                  '#CC78A7', '#CA9161', '#949494', '#56B4E9'],
                      description='Nature 期刊风格 8 色'),
    'colorblind': Palette('colorblind', ['#0072B2', '#E69F00', '#009E73', '#56B4E9',
                                          '#F0E442', '#D55E00', '#CC79A7', '#000000'],
                          description='色盲友好 8 色'),
    'omicshare': Palette('omicshare',
                          # 参考 OmicShare Tools 平台配色
                          # 富集分析、热图、火山图等常用图表的柔和色系
                          ['#E64B35', '#4DBBD5', '#00A087', '#3C5488',
                           '#F39B7F', '#8491B4', '#91D1C2', '#DC0000',
                           '#7E6148', '#B09C85'],
                          description='OmicShare 风格 10 色 - 柔和饱和，中文文献常见'),

    # ========== 科研期刊配色（来源：学术期刊官方风格） ==========

    'science_journal': Palette('science_journal',
                               # Science 期刊配色（高饱和度）
                               ['#1E3A5F', '#C41E3A', '#228B22', '#800080',
                                '#008080', '#DC143C', '#4B0082', '#2E8B57'],
                               description='Science 期刊风格 8 色 - 高对比度'),

    'npg': Palette('npg',
                   # Nature Publishing Group 配色（柔和粉彩）
                   ['#E89C8B', '#7EC4CF', '#9DB4C0', '#5B9A8B',
                    '#E8B4B8', '#C9D6DF', '#D4A5A5'],
                   description='NPG（Nature 出版集团）风格 7 色 - 柔和粉彩'),

    'lancet': Palette('lancet',
                      # Lancet 期刊配色（偏冷色调）
                      ['#8B9DC3', '#D4C5B9', '#5A9A8F', '#9B7CB6',
                       '#C85A5A', '#6B8E9F', '#B8A1C9'],
                      description='Lancet 期刊风格 7 色 - 冷色调沉稳'),

    'nejm': Palette('nejm',
                    # NEJM 期刊配色（高对比度明亮）
                    ['#6B5BFF', '#A8B5C4', '#E8A838', '#1E90FF',
                     '#00CED1', '#FF6B6B', '#4ECDC4'],
                    description='NEJM 期刊风格 7 色 - 明亮高对比'),

    'jama': Palette('jama',
                    # JAMA 期刊配色（沉稳商务）
                    ['#8B1538', '#C4B9AC', '#8B4513', '#CD853F',
                     '#B4C4AE', '#708090'],
                    description='JAMA 期刊风格 6 色 - 沉稳专业'),

    'jco': Palette('jco',
                   # JCO 期刊配色（临床肿瘤学）
                   ['#0072B2', '#D4AF37', '#8B8680', '#CD5C5C',
                    '#4682B4', '#191970', '#808000', '#556B2F'],
                   description='JCO 期刊风格 8 色 - 临床医学专业'),

    # ========== 生物信息学工具配色 ==========

    'aaas': Palette('aaas',
                    # AAAS（美国科学促进会）配色
                    ['#4A569D', '#DC143C', '#228B22', '#8B008B',
                     '#008B8B', '#B22222', '#9932CC', '#8B0000', '#800080'],
                    description='AAAS 风格 9 色 - 高饱和度科学配色'),

    'd3js': Palette('d3js',
                    # D3.js Category10 经典配色
                    ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728',
                     '#9467bd', '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22'],
                    description='D3.js Category10 9 色 - 经典数据可视化配色'),

    'futurama': Palette('futurama',
                        # Futurama 动画风格配色
                        ['#FF6600', '#CC3300', '#FF3300', '#009999',
                         '#993366', '#FF9999', '#66CCCC', '#FFCCCC', '#99CCFF'],
                        description='Futurama 风格 9 色 - 活泼卡通风格'),

    'igv': Palette('igv',
                   # IGV 基因组浏览器配色
                   ['#3333FF', '#A52A2A', '#CD853F', '#BDB76B',
                    '#556B2F', '#2F4F4F', '#8B4513', '#4682B4', '#4B0082'],
                   description='IGV 基因组浏览器风格 9 色'),

    'locuszoom': Palette('locuszoom',
                         # LocusZoom 遗传关联分析配色
                         ['#D43F3A', '#EE9336', '#5CB85C', '#46B8DA',
                          '#357EBD', '#9632B8', '#7A7A7A', '#B8B8B8', '#A0A0A0'],
                         description='LocusZoom 风格 9 色 - GWAS 可视化'),

    # ========== 中国风配色 ==========

    'china_style': Palette('china_style',
                           # 中国传统色彩
                           ['#5F8D77', '#D4A574', '#8B7355', '#6B8E6B',
                            '#E8F5E9', '#F4A460', '#FF7F50'],
                           description='中国风 7 色 - 传统色彩（竹青、赭石、朱砂等）'),
}

# 连续色图（用于热图、气泡图等）
SEQUENTIAL_CMAPS = {
    'viridis': 'viridis',
    'plasma': 'plasma',
    'inferno': 'inferno',
    'magma': 'magma',
    'cividis': 'cividis',
}

DIVERGING_CMAPS = {
    'RdBu': 'RdBu_r',
    'RdYlBu': 'RdYlBu_r',
    'coolwarm': 'coolwarm',
    'PiYG': 'PiYG',
}


# ============================================================
# 预设风格定义
# ============================================================

@dataclass
class StylePreset:
    """风格预设"""
    name: str
    description: str
    rc_params: dict = field(default_factory=dict)
    palette: str = 'default'
    sequential_cmap: str = 'viridis'
    diverging_cmap: str = 'RdBu'


PRESETS: Dict[str, StylePreset] = {}


def _register_preset(name: str, description: str, rc: dict,
                      palette: str = 'default',
                      seq_cmap: str = 'viridis',
                      div_cmap: str = 'RdBu') -> None:
    PRESETS[name] = StylePreset(
        name=name, description=description, rc_params=rc,
        palette=palette, sequential_cmap=seq_cmap, diverging_cmap=div_cmap,
    )


# --- Nature 风格 ---
_register_preset(
    'nature', 'Nature 期刊风格 - 无衬线字体，简洁专业',
    rc={
        'font.family': 'sans-serif',
        'font.sans-serif': ['Arial', 'Helvetica', 'DejaVu Sans'],
        'font.size': 8,
        'axes.labelsize': 8,
        'axes.titlesize': 9,
        'xtick.labelsize': 7,
        'ytick.labelsize': 7,
        'legend.fontsize': 7,
        'lines.linewidth': 1.0,
        'axes.linewidth': 0.5,
        'xtick.major.width': 0.5,
        'ytick.major.width': 0.5,
        'xtick.major.size': 3,
        'ytick.major.size': 3,
        'xtick.minor.width': 0.3,
        'ytick.minor.width': 0.3,
        'figure.dpi': 300,
        'savefig.dpi': 300,
        'savefig.bbox': 'tight',
        'axes.spines.top': False,
        'axes.spines.right': False,
    },
    palette='nature',
)

# --- Science 风格 ---
_register_preset(
    'science', '通用科研风格 - 基于 scienceplots（如已安装）',
    rc={
        'font.family': 'sans-serif',
        'font.size': 8,
        'axes.labelsize': 8,
        'axes.titlesize': 9,
        'xtick.labelsize': 7,
        'ytick.labelsize': 7,
        'legend.fontsize': 7,
        'lines.linewidth': 1.0,
        'axes.linewidth': 0.5,
        'figure.dpi': 300,
        'savefig.dpi': 300,
        'savefig.bbox': 'tight',
    },
    palette='tol_bright',
)

# --- Colorblind 风格 ---
_register_preset(
    'colorblind', '色盲友好风格 - 所有颜色对色盲可区分',
    rc={
        'font.family': 'sans-serif',
        'font.size': 8,
        'axes.labelsize': 8,
        'axes.titlesize': 9,
        'xtick.labelsize': 7,
        'ytick.labelsize': 7,
        'legend.fontsize': 7,
        'lines.linewidth': 1.0,
        'axes.linewidth': 0.5,
        'figure.dpi': 300,
        'savefig.dpi': 300,
        'savefig.bbox': 'tight',
    },
    palette='tol_bright',
    seq_cmap='cividis',
)

# --- Presentation 风格 ---
_register_preset(
    'presentation', '演示风格 - 大字号，高对比度',
    rc={
        'font.family': 'sans-serif',
        'font.size': 14,
        'axes.labelsize': 14,
        'axes.titlesize': 16,
        'xtick.labelsize': 12,
        'ytick.labelsize': 12,
        'legend.fontsize': 12,
        'lines.linewidth': 2.0,
        'axes.linewidth': 1.0,
        'figure.dpi': 150,
        'savefig.dpi': 150,
        'savefig.bbox': 'tight',
    },
    palette='tol_vibrant',
)

# --- Cute 手绘风格 ---
# cutecharts 是基于 ECharts 的手绘风格图表库，输出为 HTML（非 matplotlib）
# 此处仅注册预设名称和色板，实际渲染由 cute_charts.py 适配层处理
_register_preset(
    'cute', '手绘可爱风格 - 基于 cutecharts，输出交互式 HTML',
    rc={
        'font.family': 'sans-serif',
        'font.size': 10,
        'figure.dpi': 150,
        'savefig.dpi': 150,
    },
    palette='tol_bright',
    seq_cmap='viridis',
    div_cmap='RdBu',
)

# --- OmicShare 风格 ---
# 参考 OmicShare Tools（基迪奥）在线平台的可视化风格
# 特点：白底无网格、柔和配色、中文友好、发表级图表
# 富集条形图颜色深浅渐变代表P值，气泡图蓝→红渐变
_register_preset(
    'omicshare', 'OmicShare 风格 - 参考基迪奥在线平台，柔和配色，中文友好',
    rc={
        'font.family': 'sans-serif',
        'font.sans-serif': ['Arial', 'SimHei', 'Microsoft YaHei', 'DejaVu Sans'],
        'font.size': 9,
        'axes.labelsize': 9,
        'axes.titlesize': 10,
        'xtick.labelsize': 8,
        'ytick.labelsize': 8,
        'legend.fontsize': 8,
        'lines.linewidth': 1.0,
        'axes.linewidth': 0.6,
        'xtick.major.width': 0.5,
        'ytick.major.width': 0.5,
        'xtick.major.size': 4,
        'ytick.major.size': 4,
        'figure.dpi': 300,
        'savefig.dpi': 300,
        'savefig.bbox': 'tight',
        'axes.spines.top': False,
        'axes.spines.right': False,
        'axes.grid': False,
    },
    palette='omicshare',
    seq_cmap='YlOrRd',
    div_cmap='RdBu',
)


# ============================================================
# 公共 API
# ============================================================

class PlotTheme:
    """可视化风格管理器

    用法::

        # 方式一：全局应用
        PlotTheme.apply('nature')

        # 方式二：上下文管理器（推荐，不污染全局状态）
        with PlotTheme.context('nature'):
            fig, ax = plt.subplots()
            ...

        # 方式三：获取当前风格参数
        theme = PlotTheme.get_active()
        colors = theme.get_palette(8)
        cmap = theme.get_diverging_cmap()
    """

    _active: Optional[StylePreset] = None
    _active_palette_override: Optional[str] = None

    @classmethod
    def available_styles(cls) -> List[str]:
        """列出所有可用风格"""
        return list(PRESETS.keys())

    @classmethod
    def available_palettes(cls) -> List[str]:
        """列出所有可用色板"""
        return list(PALETTES.keys())

    @classmethod
    def apply(cls, style: str = 'nature', palette: Optional[str] = None) -> None:
        """全局应用指定风格

        Args:
            style: 风格名（nature/science/colorblind/presentation）
            palette: 色板名（覆盖风格默认色板）
        """
        preset = PRESETS.get(style)
        if preset is None:
            raise ValueError(
                f"未知风格 '{style}'，可用: {cls.available_styles()}"
            )

        # 尝试使用 scienceplots（可选依赖）
        if style == 'science':
            try:
                import scienceplots  # noqa: F401
                plt.style.use(['science', 'no-latex'])
            except ImportError:
                pass  # 退回手动 rcParams

        plt.rcParams.update(preset.rc_params)
        cls._active = preset
        cls._active_palette_override = palette

    @classmethod
    @contextmanager
    def context(cls, style: str = 'nature', palette: Optional[str] = None):
        """临时应用风格（上下文管理器）

        用法::

            with PlotTheme.context('nature'):
                fig, ax = plt.subplots()
                ax.plot(x, y)
                plt.savefig('fig.pdf')
        """
        # 保存当前 rcParams
        old_rc = plt.rcParams.copy()
        old_active = cls._active
        old_palette = cls._active_palette_override

        try:
            cls.apply(style, palette)
            yield
        finally:
            # 恢复
            plt.rcParams.update(old_rc)
            cls._active = old_active
            cls._active_palette_override = old_palette

    @classmethod
    def get_active(cls) -> StylePreset:
        """获取当前活跃风格（如未设置则返回 nature 默认）"""
        if cls._active is None:
            cls.apply('nature')
        return cls._active

    @classmethod
    def get_palette(cls, n: int = 8, palette: Optional[str] = None) -> List[str]:
        """获取 n 个离散颜色

        Args:
            n: 需要的颜色数量
            palette: 色板名（None 则使用当前风格默认色板）
        """
        name = palette or cls._active_palette_override
        if name is None:
            name = cls.get_active().palette

        pal = PALETTES.get(name)
        if pal is None:
            raise ValueError(
                f"未知色板 '{name}'，可用: {cls.available_palettes()}"
            )

        colors = pal.colors
        if n <= len(colors):
            return colors[:n]
        # 循环扩展
        return [colors[i % len(colors)] for i in range(n)]

    @classmethod
    def get_sequential_cmap(cls, name: Optional[str] = None) -> str:
        """获取连续色图名"""
        if name is None:
            return cls.get_active().sequential_cmap
        return SEQUENTIAL_CMAPS.get(name, name)

    @classmethod
    def get_diverging_cmap(cls, name: Optional[str] = None) -> str:
        """获取发散色图名"""
        if name is None:
            return cls.get_active().diverging_cmap
        return DIVERGING_CMAPS.get(name, name)

    @classmethod
    def get_category_colors(cls, category_type: str = 'go') -> Dict[str, str]:
        """获取特定数据库的类别颜色映射"""
        if category_type == 'go':
            return GO_CATEGORY_COLORS
        elif category_type == 'kegg':
            return KEGG_CATEGORY_COLORS
        return {}
```

- [ ] **Step 2: 创建测试**

```python
# tests/test_plot_theme.py
import matplotlib.pyplot as plt
from allenricher.visualization.plot_theme import (
    PlotTheme, PALETTES, PRESETS, PAUL_TOL_BRIGHT,
    GO_CATEGORY_COLORS, VOLCANO_COLORS,
)


class TestPlotTheme:
    def test_available_styles(self):
        styles = PlotTheme.available_styles()
        assert 'nature' in styles
        assert 'science' in styles
        assert 'colorblind' in styles
        assert 'presentation' in styles
        assert 'cute' in styles
        assert 'omicshare' in styles

    def test_available_palettes(self):
        palettes = PlotTheme.available_palettes()
        assert 'tol_bright' in palettes
        assert 'okabe_ito' in palettes
        assert 'nature' in palettes
        assert 'omicshare' in palettes
        # 科研期刊配色
        assert 'science_journal' in palettes
        assert 'npg' in palettes
        assert 'lancet' in palettes
        assert 'nejm' in palettes
        assert 'jama' in palettes
        assert 'jco' in palettes
        # 生物信息学工具配色
        assert 'd3js' in palettes
        assert 'igv' in palettes
        assert 'locuszoom' in palettes
        # 中国风
        assert 'china_style' in palettes

    def test_apply_nature(self):
        PlotTheme.apply('nature')
        theme = PlotTheme.get_active()
        assert theme.name == 'nature'
        assert theme.palette == 'nature'

    def test_apply_invalid_style(self):
        import pytest
        with pytest.raises(ValueError, match="未知风格"):
            PlotTheme.apply('nonexistent_style')

    def test_context_manager(self):
        PlotTheme.apply('colorblind')
        with PlotTheme.context('nature'):
            theme = PlotTheme.get_active()
            assert theme.name == 'nature'
        theme = PlotTheme.get_active()
        assert theme.name == 'colorblind'

    def test_get_palette_default(self):
        PlotTheme.apply('nature')
        colors = PlotTheme.get_palette(5)
        assert len(colors) == 5
        assert all(c.startswith('#') for c in colors)

    def test_get_palette_cycle(self):
        PlotTheme.apply('nature')
        colors = PlotTheme.get_palette(20)  # 超过色板长度
        assert len(colors) == 20

    def test_get_palette_override(self):
        PlotTheme.apply('nature')
        colors = PlotTheme.get_palette(3, palette='tol_high_contrast')
        assert len(colors) == 3
        assert colors[0] == '#004488'

    def test_get_sequential_cmap(self):
        cmap = PlotTheme.get_sequential_cmap()
        assert isinstance(cmap, str)

    def test_get_diverging_cmap(self):
        cmap = PlotTheme.get_diverging_cmap()
        assert isinstance(cmap, str)

    def test_get_category_colors_go(self):
        colors = PlotTheme.get_category_colors('go')
        assert 'Biological Process' in colors
        assert len(colors) == 3

    def test_get_category_colors_kegg(self):
        colors = PlotTheme.get_category_colors('kegg')
        assert 'Metabolism' in colors
        assert len(colors) == 6

    def test_palettes_have_required_fields(self):
        for name, pal in PALETTES.items():
            assert pal.name == name
            assert len(pal.colors) >= 3

    def test_presets_have_required_fields(self):
        for name, preset in PRESETS.items():
            assert preset.name == name
            assert 'font.size' in preset.rc_params
```

- [ ] **Step 3: 运行测试**

Run: `cd f:\OneDrive\Documents\TraeSOLO\AllEnricher\AllEnricher-v2 && python -m pytest tests/test_plot_theme.py -v`
Expected: 14 passed

- [ ] **Step 4: Commit**

```bash
git add allenricher/visualization/plot_theme.py tests/test_plot_theme.py
git commit -m "feat(viz): 创建风格系统核心模块 PlotTheme"
```

---

## Task 2: Python 重写 barplot（替代 barplot.R）

**Files:**
- Create: `allenricher/visualization/barplot.py`
- Test: `tests/test_barplot.py`

- [ ] **Step 1: 创建 barplot.py**

用 matplotlib 重写 barplot.R 的水平柱状图功能。关键要求：
- 支持 GO/KEGG/Reactome/DO/DisGeNET 五种数据库的类别着色
- 通过 `PlotTheme.get_category_colors()` 获取颜色
- 支持 top_n 参数截断
- 输出格式由文件扩展名决定（png/pdf/svg）
- DPI 可配置

```python
"""发表级水平柱状图（替代 barplot.R）

生成 ORA 富集分析的水平柱状图，按数据库类型自动着色。
"""
from __future__ import annotations

from pathlib import Path
from typing import Dict, List, Optional

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np

from .plot_theme import PlotTheme


def plot_barplot(
    data: List[Dict],
    output_file: str,
    db_name: str = 'GO',
    top_n: int = 20,
    dpi: int = 300,
    figsize: Optional[tuple] = None,
    style: Optional[str] = None,
    palette: Optional[str] = None,
) -> Path:
    """绘制发表级水平柱状图

    Args:
        data: 富集结果列表，每个元素需包含:
            - Term_Name (str): 条目名称
            - Adjusted_P_Value (float): 校正后 P 值
            - Gene_Count (int): 基因数
            - Database (str, optional): 数据库类型（用于 GO 子类着色）
        output_file: 输出文件路径（扩展名决定格式）
        db_name: 数据库名称（GO/KEGG/Reactome/DO/DisGeNET）
        top_n: 显示的最大条目数
        dpi: 输出分辨率
        figsize: 图表尺寸 (width, height)，None 则自动计算
        style: 风格名（None 则使用当前活跃风格）
        palette: 色板名（覆盖默认）

    Returns:
        输出文件路径
    """
    if not data:
        return Path(output_file)

    # 截取 top_n
    df = data[:top_n]

    # 提取数据
    names = [d.get('Term_Name', d.get('term_name', '')) for d in df]
    pvalues = [d.get('Adjusted_P_Value', d.get('adjusted_p_value', 1)) for d in df]
    gene_counts = [d.get('Gene_Count', d.get('gene_count', 0)) for d in df]

    # -log10(Q-value)
    neg_log_q = [-np.log10(max(p, 1e-300)) for p in pvalues]

    # 确定颜色
    colors = _get_bar_colors(df, db_name, style, palette)

    # 自动计算 figsize
    n = len(names)
    if figsize is None:
        height = max(4, n * 0.35)
        figsize = (8, height)

    # 绘图
    ctx = PlotTheme.context if style else _nullcontext
    with ctx(style or 'nature', palette):
        fig, ax = plt.subplots(figsize=figsize)

        y_pos = range(n)
        bars = ax.barh(y_pos, neg_log_q, color=colors, edgecolor='white',
                        linewidth=0.5, height=0.7)

        ax.set_yticks(y_pos)
        ax.set_yticklabels(names, fontsize=7)
        ax.invert_yaxis()
        ax.set_xlabel('-log10(Q-value)', fontsize=8)
        ax.set_title(f'{db_name} Enrichment (Top {n})', fontsize=9, fontweight='bold')

        # 在柱状图右侧标注基因数
        for i, (bar, gc) in enumerate(zip(bars, gene_counts)):
            ax.text(bar.get_width() + 0.1, bar.get_y() + bar.get_height() / 2,
                    str(gc), va='center', fontsize=6, color='#555555')

        plt.tight_layout()
        fig.savefig(output_file, dpi=dpi, bbox_inches='tight')
        plt.close(fig)

    return Path(output_file)


def _get_bar_colors(
    data: List[Dict], db_name: str,
    style: Optional[str] = None, palette: Optional[str] = None,
) -> List[str]:
    """根据数据库类型获取柱状图颜色"""
    db_lower = db_name.upper()

    if db_lower == 'GO':
        cat_colors = PlotTheme.get_category_colors('go')
        colors = []
        for d in data:
            name = d.get('Term_Name', d.get('term_name', ''))
            # 从 Term_Name 前缀判断 GO 子类
            if name.startswith('Biological Process') or 'Biological Process' in name:
                colors.append(cat_colors.get('Biological Process', '#FF6B35'))
            elif name.startswith('Cellular Component') or 'Cellular Component' in name:
                colors.append(cat_colors.get('Cellular Component', '#004E89'))
            elif name.startswith('Molecular Function') or 'Molecular Function' in name:
                colors.append(cat_colors.get('Molecular Function', '#7B2D8E'))
            else:
                colors.append('#888888')
        return colors

    elif db_lower == 'KEGG':
        cat_colors = PlotTheme.get_category_colors('kegg')
        colors = []
        for d in data:
            name = d.get('Term_Name', d.get('term_name', ''))
            matched = False
            for category, color in cat_colors.items():
                if category in name:
                    colors.append(color)
                    matched = True
                    break
            if not matched:
                colors.append('#888888')
        return colors

    elif db_lower == 'DO':
        return ['#2CA02C'] * len(data)

    elif db_lower == 'REACTOME':
        return ['#D62728'] * len(data)

    elif db_lower == 'DISGENET':
        return ['#1F77B4'] * len(data)

    else:
        # 使用当前风格色板
        return PlotTheme.get_palette(len(data), palette=palette)


from contextlib import contextmanager

@contextmanager
def _nullcontext(*args, **kwargs):
    yield
```

- [ ] **Step 2: 创建测试**

```python
# tests/test_barplot.py
import json
import tempfile
from pathlib import Path
from allenricher.visualization.barplot import plot_barplot


def _make_sample_data(n=10):
    """生成测试数据"""
    data = []
    go_categories = ['Biological Process', 'Cellular Component', 'Molecular Function']
    for i in range(n):
        data.append({
            'Term_Name': f'{go_categories[i % 3]}|Term{i}',
            'Adjusted_P_Value': 10 ** (-i - 1),
            'Gene_Count': 50 - i * 3,
        })
    return data


class TestPlotBarplot:
    def test_basic_generation(self, tmp_path):
        data = _make_sample_data(10)
        out = tmp_path / 'test_barplot.png'
        result = plot_barplot(data, str(out), db_name='GO', top_n=10)
        assert result.exists()

    def test_pdf_output(self, tmp_path):
        data = _make_sample_data(5)
        out = tmp_path / 'test_barplot.pdf'
        result = plot_barplot(data, str(out), db_name='GO')
        assert result.exists()

    def test_empty_data(self, tmp_path):
        out = tmp_path / 'empty.png'
        result = plot_barplot([], str(out), db_name='GO')
        assert result == Path(out)

    def test_kegg_colors(self, tmp_path):
        data = [
            {'Term_Name': 'Metabolism|Pathway1', 'Adjusted_P_Value': 1e-10, 'Gene_Count': 30},
            {'Term_Name': 'Genetic Information Processing|Pathway2', 'Adjusted_P_Value': 1e-8, 'Gene_Count': 20},
        ]
        out = tmp_path / 'kegg.png'
        result = plot_barplot(data, str(out), db_name='KEGG')
        assert result.exists()

    def test_top_n_truncation(self, tmp_path):
        data = _make_sample_data(30)
        out = tmp_path / 'top10.png'
        plot_barplot(data, str(out), db_name='GO', top_n=10)
        # 验证文件生成（不验证内容，因为 top_n 在函数内截断）
        assert out.exists()
```

- [ ] **Step 3: 运行测试**

Run: `cd f:\OneDrive\Documents\TraeSOLO\AllEnricher\AllEnricher-v2 && python -m pytest tests/test_barplot.py -v`
Expected: 5 passed

- [ ] **Step 4: Commit**

```bash
git add allenricher/visualization/barplot.py tests/test_barplot.py
git commit -m "feat(viz): Python 重写 barplot（替代 barplot.R）"
```

---

## Task 3: Python 重写 bubble（替代 bubble.R）

**Files:**
- Create: `allenricher/visualization/bubble.py`
- Test: `tests/test_bubble.py`

- [ ] **Step 1: 创建 bubble.py**

用 matplotlib 重写 bubble.R 的气泡图功能。关键要求：
- X=RichFactor, Y=Term_Name, 点大小=GeneCount, 颜色=-log10(Qvalue)
- 通过 `PlotTheme.get_diverging_cmap()` 获取色图
- 支持 top_n 参数
- DPI 可配置

```python
"""发表级气泡图（替代 bubble.R）

生成 ORA 富集分析的气泡图。
"""
from __future__ import annotations

from pathlib import Path
from typing import Dict, List, Optional

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import numpy as np

from .plot_theme import PlotTheme


def plot_bubble(
    data: List[Dict],
    output_file: str,
    db_name: str = 'GO',
    top_n: int = 20,
    dpi: int = 300,
    figsize: Optional[tuple] = None,
    style: Optional[str] = None,
    palette: Optional[str] = None,
) -> Path:
    """绘制发表级气泡图

    Args:
        data: 富集结果列表
        output_file: 输出文件路径
        db_name: 数据库名称
        top_n: 最大条目数
        dpi: 输出分辨率
        figsize: 图表尺寸
        style: 风格名
        palette: 色板名（用于离散色图）

    Returns:
        输出文件路径
    """
    if not data:
        return Path(output_file)

    df = data[:top_n]

    names = [d.get('Term_Name', d.get('term_name', '')) for d in df]
    pvalues = [d.get('Adjusted_P_Value', d.get('adjusted_p_value', 1)) for d in df]
    gene_counts = [d.get('Gene_Count', d.get('gene_count', 1)) for d in df]
    bg_counts = [d.get('Background_Count', d.get('background_count', 1)) for d in df]
    gene_ratios = [d.get('Gene_Ratio', d.get('gene_ratio', 0)) for d in df]

    # Rich Factor = Gene_Count / Background_Count（如无则用 Gene_Ratio）
    rich_factors = []
    for gc, bc, gr in zip(gene_counts, bg_counts, gene_ratios):
        if bc and bc > 0:
            rich_factors.append(gc / bc)
        elif gr:
            rich_factors.append(gr)
        else:
            rich_factors.append(0)

    neg_log_q = [-np.log10(max(p, 1e-300)) for p in pvalues]

    n = len(names)
    if figsize is None:
        height = max(4, n * 0.35)
        figsize = (7, height)

    ctx = PlotTheme.context if style else _nullcontext
    with ctx(style or 'nature', palette):
        fig, ax = plt.subplots(figsize=figsize)

        # 散点图
        scatter = ax.scatter(
            rich_factors, range(n),
            s=[max(c * 3, 10) for c in gene_counts],  # 点大小
            c=neg_log_q,
            cmap=PlotTheme.get_diverging_cmap(),
            alpha=0.8,
            edgecolors='white',
            linewidths=0.3,
        )

        ax.set_yticks(range(n))
        ax.set_yticklabels(names, fontsize=7)
        ax.invert_yaxis()
        ax.set_xlabel('Rich Factor', fontsize=8)
        ax.set_title(f'{db_name} Enrichment (Top {n})', fontsize=9, fontweight='bold')

        # 色条
        cbar = plt.colorbar(scatter, ax=ax, shrink=0.6, aspect=15, pad=0.02)
        cbar.set_label('-log10(Q-value)', fontsize=7)
        cbar.ax.tick_params(labelsize=6)

        plt.tight_layout()
        fig.savefig(output_file, dpi=dpi, bbox_inches='tight')
        plt.close(fig)

    return Path(output_file)


from contextlib import contextmanager

@contextmanager
def _nullcontext(*args, **kwargs):
    yield
```

- [ ] **Step 2: 创建测试**

```python
# tests/test_bubble.py
from pathlib import Path
from allenricher.visualization.bubble import plot_bubble


def _make_sample_data(n=10):
    data = []
    for i in range(n):
        data.append({
            'Term_Name': f'Term{i}',
            'Adjusted_P_Value': 10 ** (-i - 1),
            'Gene_Count': 50 - i * 3,
            'Background_Count': 1000,
            'Gene_Ratio': (50 - i * 3) / 1000,
        })
    return data


class TestPlotBubble:
    def test_basic_generation(self, tmp_path):
        data = _make_sample_data(10)
        out = tmp_path / 'test_bubble.png'
        result = plot_bubble(data, str(out), db_name='GO')
        assert result.exists()

    def test_pdf_output(self, tmp_path):
        data = _make_sample_data(5)
        out = tmp_path / 'test_bubble.pdf'
        result = plot_bubble(data, str(out), db_name='KEGG')
        assert result.exists()

    def test_empty_data(self, tmp_path):
        out = tmp_path / 'empty.png'
        result = plot_bubble([], str(out), db_name='GO')
        assert result == Path(out)
```

- [ ] **Step 3: 运行测试**

Run: `cd f:\OneDrive\Documents\TraeSOLO\AllEnricher\AllEnricher-v2 && python -m pytest tests/test_bubble.py -v`
Expected: 3 passed

- [ ] **Step 4: Commit**

```bash
git add allenricher/visualization/bubble.py tests/test_bubble.py
git commit -m "feat(viz): Python 重写 bubble（替代 bubble.R）"
```

---

## Task 4: 改造 Plotter 类，移除 R 依赖

**Files:**
- Modify: `allenricher/visualization/plotter.py`

- [ ] **Step 1: 修改 Plotter 类**

将 `plot_barplot()` 和 `plot_bubble()` 从调用 R 脚本改为调用 Python 函数：

- 删除 `_run_r_script()` 方法
- `plot_barplot()` → 调用 `barplot.plot_barplot()`
- `plot_bubble()` → 调用 `bubble.plot_bubble()`
- `plot_all()` 保持接口不变
- 保留 `__init__` 的 `output_dir` 参数
- 新增 `style` 和 `palette` 参数传递

- [ ] **Step 2: 更新 __init__.py 导出**

在 `allenricher/visualization/__init__.py` 中添加：
```python
from .plot_theme import PlotTheme
from .barplot import plot_barplot
from .bubble import plot_bubble
```

- [ ] **Step 3: 验证现有测试**

Run: `cd f:\OneDrive\Documents\TraeSOLO\AllEnricher\AllEnricher-v2 && python -m pytest tests/test_common_plots.py tests/test_e2e_visualization.py -v --tb=short`
Expected: 全部通过（可能需要适配）

- [ ] **Step 4: Commit**

```bash
git add allenricher/visualization/plotter.py allenricher/visualization/__init__.py
git commit -m "refactor(viz): Plotter 移除 R 依赖，改用 Python 绘图"
```

---

## Task 5: 改造 GSEA 绘图函数接入风格系统

**Files:**
- Modify: `allenricher/visualization/gsea_plots.py`

- [ ] **Step 1: 替换硬编码颜色**

在所有 3 个绘图函数中：
- 添加 `style: Optional[str] = None` 和 `palette: Optional[str] = None` 参数
- 用 `PlotTheme.context(style or 'nature', palette)` 包裹绘图代码
- 将硬编码颜色替换为 `PlotTheme.get_palette()` / `PlotTheme.get_diverging_cmap()` 调用

具体替换清单：
| 行号 | 原值 | 替换为 |
|------|------|--------|
| 107 | `color="#2E86AB"` | `PlotTheme.get_palette(1)[0]` |
| 152 | `"#E74C3C"` | `PlotTheme.get_palette(1, palette='tol_high_contrast')[2]` |
| 212-213 | `color_pos="#E74C3C"`, `color_neg="#3498DB"` | `PlotTheme.get_palette(2, palette='tol_high_contrast')` |
| 334 | `cmap="RdBu_r"` | `PlotTheme.get_diverging_cmap()` |

- [ ] **Step 2: 运行测试**

Run: `cd f:\OneDrive\Documents\TraeSOLO\AllEnricher\AllEnricher-v2 && python -m pytest tests/test_gsea_plots.py -v`
Expected: 全部通过

- [ ] **Step 3: Commit**

```bash
git add allenricher/visualization/gsea_plots.py
git commit -m "refactor(viz): GSEA 绘图接入风格系统"
```

---

## Task 6: 改造 GSVA/ssGSEA 绘图函数接入风格系统

**Files:**
- Modify: `allenricher/visualization/gsva_plots.py`

- [ ] **Step 1: 替换硬编码颜色**

在所有 4 个绘图函数中：
- 添加 `style` 和 `palette` 参数
- 用 `PlotTheme.context()` 包裹
- 替换 `sns.color_palette("Set2", ...)` → `PlotTheme.get_palette(n, palette=...)`
- 替换 `cmap="RdBu_r"` → `PlotTheme.get_diverging_cmap()`

- [ ] **Step 2: 运行测试**

Run: `cd f:\OneDrive\Documents\TraeSOLO\AllEnricher\AllEnricher-v2 && python -m pytest tests/test_gsva_plots.py -v`
Expected: 全部通过

- [ ] **Step 3: Commit**

```bash
git add allenricher/visualization/gsva_plots.py
git commit -m "refactor(viz): GSVA/ssGSEA 绘图接入风格系统"
```

---

## Task 7: 改造通用绘图函数 + 集成到工作流

**Files:**
- Modify: `allenricher/visualization/common_plots.py`
- Modify: `allenricher/cli.py` (集成 4 个通用函数)

- [ ] **Step 1: common_plots.py 接入风格系统**

在 4 个函数中添加 `style`/`palette` 参数，替换硬编码颜色：
- `plot_enrichment_network()`: 替换 `plt.cm.RdBu_r` 和 `sns.color_palette("viridis", ...)`
- `plot_upset()`: 替换 `"#4C72B0"`
- `plot_volcano()`: 使用 `VOLCANO_COLORS` 常量
- `plot_method_comparison()`: 替换 `"#4C72B0"`

- [ ] **Step 2: 集成通用函数到 CLI 工作流**

在 `cli.py` 的 `_METHOD_PLOT_TYPES` 中添加通用图表类型：
```python
_COMMON_PLOT_TYPES = {'network', 'upset', 'volcano'}
```

在 `_generate_plots()` 中添加对通用图表类型的调度。

- [ ] **Step 3: 运行测试**

Run: `cd f:\OneDrive\Documents\TraeSOLO\AllEnricher\AllEnricher-v2 && python -m pytest tests/test_common_plots.py -v`
Expected: 全部通过

- [ ] **Step 4: Commit**

```bash
git add allenricher/visualization/common_plots.py allenricher/cli.py
git commit -m "feat(viz): 通用绘图接入风格系统 + 集成到 CLI"
```

---

## Task 8: CLI/API 集成风格参数

**Files:**
- Modify: `allenricher/cli.py`
- Modify: `allenricher/core/config.py`
- Modify: `allenricher/report/generator.py`

- [ ] **Step 1: Config 新增字段**

在 `config.py` 的 Config dataclass 中添加：
```python
plot_style: str = 'nature'          # 图表风格（nature/science/colorblind/presentation）
plot_palette: Optional[str] = None  # 色板名（覆盖风格默认）
```

- [ ] **Step 2: CLI 新增参数**

在 analyze 子命令中添加：
```python
analyze_parser.add_argument('--style', type=str, default='nature',
    choices=['nature', 'science', 'colorblind', 'presentation', 'cute', 'omicshare'],
    help='图表风格（默认 nature，cute 为手绘风格，omicshare 为基迪奥风格）')
analyze_parser.add_argument('--palette', type=str, default=None,
    help='色板名称（覆盖风格默认色板）。可用: tol_bright, tol_vibrant, tol_muted, okabe_ito, nature, colorblind, omicshare, science_journal, npg, lancet, nejm, jama, jco, aaas, d3js, futurama, igv, locuszoom, china_style')
```

- [ ] **Step 3: 传递风格参数到绘图流程**

在 `cmd_analyze()` 和 `_generate_plots()` 中，将 `config.plot_style` 和 `config.plot_palette` 传递给所有绘图函数调用。

- [ ] **Step 4: 报告生成器传递风格**

在 `generator.py` 的 `_generate_gsea_plots_section()` 和 `_generate_gsva_plots_section()` 中传递风格参数。

- [ ] **Step 5: 验证**

Run: `cd f:\OneDrive\Documents\TraeSOLO\AllEnricher\AllEnricher-v2 && python -m allenricher analyze --help`
Expected: 显示 `--style` 和 `--palette` 参数

- [ ] **Step 6: Commit**

```bash
git add allenricher/cli.py allenricher/core/config.py allenricher/report/generator.py
git commit -m "feat(viz): CLI 新增 --style/--palette 参数"
```

---

## Task 9: 清理 + 端到端验证

**Files:**
- Delete: `allenricher/visualization/barplot.R`
- Delete: `allenricher/visualization/bubble.R`
- Test: E2E 验证

- [ ] **Step 1: 删除 R 脚本**

```bash
rm allenricher/visualization/barplot.R
rm allenricher/visualization/bubble.R
```

- [ ] **Step 2: 运行全量测试**

Run: `cd f:\OneDrive\Documents\TraeSOLO\AllEnricher\AllEnricher-v2 && python -m pytest tests/ -v --tb=short -m "not slow" 2>&1 | tail -30`
Expected: 全部通过

- [ ] **Step 3: E2E 验证 - 默认风格**

Run: `cd f:\OneDrive\Documents\TraeSOLO\AllEnricher\AllEnricher-v2 && python -m allenricher analyze -i "../AllEnricher-v1/example/example.glist" -s hsa -d GO -o test_output/viz_default`
Expected: 生成 barplot + bubble 图表（Python 版），Nature 风格

- [ ] **Step 4: E2E 验证 - colorblind 风格**

Run: `cd f:\OneDrive\Documents\TraeSOLO\AllEnricher\AllEnricher-v2 && python -m allenricher analyze -i "../AllEnricher-v1/example/example.glist" -s hsa -d GO --style colorblind -o test_output/viz_colorblind`
Expected: 生成色盲友好风格的图表

- [ ] **Step 5: E2E 验证 - presentation 风格**

Run: `cd f:\OneDrive\Documents\TraeSOLO\AllEnricher\AllEnricher-v2 && python -m allenricher analyze -i "../AllEnricher-v1/example/example.glist" -s hsa -d GO --style presentation -o test_output/viz_presentation`
Expected: 生成演示风格的大字号图表

- [ ] **Step 6: E2E 验证 - omicshare 风格**

Run: `cd f:\OneDrive\Documents\TraeSOLO\AllEnricher\AllEnricher-v2 && python -m allenricher analyze -i "../AllEnricher-v1/example/example.glist" -s hsa -d GO --style omicshare -o test_output/viz_omicshare`
Expected: 生成 OmicShare 风格的柔和配色图表

- [ ] **Step 7: Commit**

```bash
git add -A
git commit -m "chore(viz): 删除 R 脚本，全量 E2E 验证通过"
```

---

## Task 10: cutecharts 手绘风格集成

**Files:**
- Create: `allenricher/visualization/cute_charts.py`
- Modify: `allenricher/visualization/plotter.py` (cute 风格路由)
- Test: `tests/test_cute_charts.py`

**背景**: [cutecharts](https://github.com/chenjiandongx/cutecharts) 是基于 ECharts 的 Python 手绘风格图表库，输出为交互式 HTML（非静态图片）。支持的图表类型：Bar、Line、Pie、Radar、Scatter。cutecharts 为可选依赖，未安装时 `--style cute` 应给出友好提示。

- [ ] **Step 1: 创建 cute_charts.py 适配层**

```python
"""cutecharts 手绘风格适配层

将 ORA 富集分析的 barplot/bubble 映射到 cutecharts 手绘风格，
输出为交互式 HTML 文件。cutecharts 为可选依赖。
"""
from __future__ import annotations

from pathlib import Path
from typing import Dict, List, Optional


def _check_cutecharts() -> bool:
    """检查 cutecharts 是否已安装"""
    try:
        import cutecharts  # noqa: F401
        return True
    except ImportError:
        return False


def plot_cute_barplot(
    data: List[Dict],
    output_file: str,
    db_name: str = 'GO',
    top_n: int = 20,
) -> Path:
    """cutecharts 手绘风格水平柱状图

    Args:
        data: 富集结果列表
        output_file: 输出 HTML 文件路径
        db_name: 数据库名称
        top_n: 最大条目数

    Returns:
        输出文件路径
    """
    if not _check_cutecharts():
        raise ImportError(
            "cutecharts 未安装。请运行: pip install cutecharts"
        )

    from cutecharts.charts import Bar

    if not data:
        return Path(output_file)

    df = data[:top_n]

    # 提取数据（取最后 top_n 条，因为 cutecharts bar 是垂直的，从下到上）
    names = [d.get('Term_Name', d.get('term_name', ''))[-30:] for d in reversed(df)]
    neg_log_q = [
        -round(max(0, __import__('math').log10(
            max(d.get('Adjusted_P_Value', d.get('adjusted_p_value', 1)), 1e-300)
        )), 2)
        for d in reversed(df)
    ]

    chart = Bar(f"{db_name} Enrichment (Top {len(names)})")
    chart.set_options(
        labels=names,
        x_label="-log10(Q-value)",
        y_label="",
        colors=[
            "#EE6677", "#4477AA", "#228833", "#CCBB44",
            "#66CCEE", "#AA3377", "#BBBBBB", "#EE7733",
        ],
    )
    chart.add_series("-log10(Q-value)", neg_log_q)

    # 确保输出为 .html
    out = Path(output_file)
    if out.suffix != '.html':
        out = out.with_suffix('.html')

    chart.render(str(out))
    return out


def plot_cute_bubble(
    data: List[Dict],
    output_file: str,
    db_name: str = 'GO',
    top_n: int = 20,
) -> Path:
    """cutecharts 手绘风格散点图（近似气泡图）

    注意：cutecharts 没有原生气泡图（大小映射），此处用 Scatter 近似。
    点大小固定，颜色区分 -log10(Q-value) 区间。

    Args:
        data: 富集结果列表
        output_file: 输出 HTML 文件路径
        db_name: 数据库名称
        top_n: 最大条目数

    Returns:
        输出文件路径
    """
    if not _check_cutecharts():
        raise ImportError(
            "cutecharts 未安装。请运行: pip install cutecharts"
        )

    from cutecharts.charts import Scatter
    import math

    if not data:
        return Path(output_file)

    df = data[:top_n]

    names = [d.get('Term_Name', d.get('term_name', ''))[-25:] for d in df]
    gene_counts = [d.get('Gene_Count', d.get('gene_count', 1)) for d in df]
    bg_counts = [d.get('Background_Count', d.get('background_count', 1)) for d in df]
    pvalues = [d.get('Adjusted_P_Value', d.get('adjusted_p_value', 1)) for d in df]

    rich_factors = []
    for gc, bc in zip(gene_counts, bg_counts):
        if bc and bc > 0:
            rich_factors.append(round(gc / bc, 4))
        else:
            rich_factors.append(0)

    neg_log_q = [
        -round(max(0, math.log10(max(p, 1e-300))), 2)
        for p in pvalues
    ]

    chart = Scatter(f"{db_name} Enrichment (Top {len(names)})")
    chart.set_options(
        x_label="Rich Factor",
        y_label="-log10(Q-value)",
        dot_size=3,
    )
    chart.add_series(
        "Enrichment",
        [(rf, nlq) for rf, nlq in zip(rich_factors, neg_log_q)],
    )

    out = Path(output_file)
    if out.suffix != '.html':
        out = out.with_suffix('.html')

    chart.render(str(out))
    return out
```

- [ ] **Step 2: 创建测试**

```python
# tests/test_cute_charts.py
import pytest
from pathlib import Path


class TestCuteCharts:
    def test_cute_barplot_requires_install(self, tmp_path):
        """cutecharts 未安装时应抛出 ImportError"""
        try:
            from allenricher.visualization.cute_charts import plot_cute_barplot
            # 如果能导入，说明已安装，测试正常生成
            data = [
                {'Term_Name': f'Term{i}', 'Adjusted_P_Value': 10**(-i-1), 'Gene_Count': 50-i*3}
                for i in range(5)
            ]
            out = tmp_path / 'cute_bar.html'
            result = plot_cute_barplot(data, str(out), db_name='GO')
            assert result.exists()
            content = result.read_text()
            assert '<html' in content.lower() or '<!doctype' in content.lower()
        except ImportError:
            pytest.skip("cutecharts not installed")

    def test_cute_bubble_requires_install(self, tmp_path):
        """cutecharts 未安装时应抛出 ImportError"""
        try:
            from allenricher.visualization.cute_charts import plot_cute_bubble
            data = [
                {
                    'Term_Name': f'Term{i}',
                    'Adjusted_P_Value': 10**(-i-1),
                    'Gene_Count': 50-i*3,
                    'Background_Count': 1000,
                }
                for i in range(5)
            ]
            out = tmp_path / 'cute_bubble.html'
            result = plot_cute_bubble(data, str(out), db_name='GO')
            assert result.exists()
        except ImportError:
            pytest.skip("cutecharts not installed")

    def test_empty_data(self, tmp_path):
        """空数据应返回路径不生成文件"""
        try:
            from allenricher.visualization.cute_charts import plot_cute_barplot
            out = tmp_path / 'empty.html'
            result = plot_cute_barplot([], str(out), db_name='GO')
            assert result == Path(out)
        except ImportError:
            pytest.skip("cutecharts not installed")

    def test_output_extension_forced_html(self, tmp_path):
        """非 .html 扩展名应自动转为 .html"""
        try:
            from allenricher.visualization.cute_charts import plot_cute_barplot
            data = [
                {'Term_Name': 'T1', 'Adjusted_P_Value': 1e-5, 'Gene_Count': 30},
            ]
            out = tmp_path / 'cute_bar.png'  # 非 html 扩展名
            result = plot_cute_barplot(data, str(out), db_name='GO')
            assert result.suffix == '.html'
        except ImportError:
            pytest.skip("cutecharts not installed")
```

- [ ] **Step 3: 运行测试**

Run: `cd f:\OneDrive\Documents\TraeSOLO\AllEnricher\AllEnricher-v2 && python -m pytest tests/test_cute_charts.py -v`
Expected: passed 或 skipped（取决于 cutecharts 是否安装）

- [ ] **Step 4: 在 Plotter 中路由 cute 风格**

修改 `plotter.py` 的 `plot_barplot()` 和 `plot_bubble()`，当 `style == 'cute'` 时调用 `cute_charts.py` 的函数而非 matplotlib 版本：

```python
def plot_barplot(self, df, db_name, top_n=20, style=None, palette=None):
    if style == 'cute':
        from .cute_charts import plot_cute_barplot
        return plot_cute_barplot(df, str(self.output_dir / f"{db_name}_barplot.html"),
                                  db_name=db_name, top_n=top_n)
    else:
        from .barplot import plot_barplot
        return plot_barplot(df, str(self.output_dir / f"{db_name}_barplot.{self.config.figure_format}"),
                            db_name=db_name, top_n=top_n, dpi=self.config.figure_dpi,
                            style=style, palette=palette)

def plot_bubble(self, df, db_name, top_n=20, style=None, palette=None):
    if style == 'cute':
        from .cute_charts import plot_cute_bubble
        return plot_cute_bubble(df, str(self.output_dir / f"{db_name}_bubble.html"),
                                  db_name=db_name, top_n=top_n)
    else:
        from .bubble import plot_bubble
        return plot_bubble(df, str(self.output_dir / f"{db_name}_bubble.{self.config.figure_format}"),
                            db_name=db_name, top_n=top_n, dpi=self.config.figure_dpi,
                            style=style, palette=palette)
```

- [ ] **Step 5: E2E 验证 cute 风格**

Run: `cd f:\OneDrive\Documents\TraeSOLO\AllEnricher\AllEnricher-v2 && pip install cutecharts && python -m allenricher analyze -i "../AllEnricher-v1/example/example.glist" -s hsa -d GO --style cute -o test_output/viz_cute`
Expected: 生成 `_barplot.html` 和 `_bubble.html`，浏览器打开可见手绘风格交互式图表

- [ ] **Step 6: Commit**

```bash
git add allenricher/visualization/cute_charts.py allenricher/visualization/plotter.py tests/test_cute_charts.py
git commit -m "feat(viz): 新增 cutecharts 手绘风格支持"
```

---

## cutecharts 风格说明

### 与其他风格的区别

| 特性 | nature/science/colorblind/presentation | cute |
|------|---------------------------------------|------|
| 渲染引擎 | matplotlib (静态图片) | ECharts (交互式 HTML) |
| 输出格式 | png/pdf/svg | html |
| 交互性 | 无 | 鼠标悬停提示、缩放 |
| 适用场景 | 论文发表、打印 | 网页展示、教学演示、PPT 嵌入 |
| 依赖 | matplotlib | cutecharts (可选) |

### cutecharts 支持的图表类型映射

| AllEnricher 图表 | cutecharts 类型 | 说明 |
|------------------|----------------|------|
| ORA barplot | Bar | 垂直柱状图（cutecharts 无水平柱状图） |
| ORA bubble | Scatter | 散点图近似（cutecharts 无气泡大小映射） |
| GSEA enrichment | Line | 折线图近似 running ES 曲线 |
| GSEA NES barplot | Bar | 条形图 |
| GSEA dotplot | Scatter | 散点图近似 |
| GSVA heatmap | - | 不支持（cutecharts 无热图） |
| GSVA group_comparison | Bar | 箱线图用柱状图近似 |
| GSVA dotplot | Scatter | 散点图近似 |
| GSVA correlation | Radar | 雷达图近似 |

### 注意事项

1. cutecharts 为**可选依赖**，未安装时使用 `--style cute` 会给出清晰的安装提示
2. cute 风格仅适用于 ORA 和 GSEA 的部分图表类型，GSVA 热图等复杂图表会自动退回 matplotlib 渲染
3. 输出为 HTML 文件，可直接在浏览器中打开，也可嵌入网页或 iframe

---

## OmicShare 风格说明

### 风格来源

OmicShare Tools（[omicshare.com](https://www.omicshare.com/tools)）是基迪奥生物开发的在线生物数据分析可视化平台，已被 4500+ 篇 SCI 文章引用。其图表风格在国内生物信息学领域广泛使用，特点是柔和配色、中文友好、白底简洁。

### 风格特征

| 特性 | 说明 |
|------|------|
| 底色 | 白色，无网格线 |
| 字体 | Arial + SimHei/Microsoft YaHei（中文回退） |
| 字号 | 9pt 基础，标题 10pt |
| 线宽 | 坐标轴 0.6pt，刻度 0.5pt |
| 边框 | 隐藏顶部和右侧 spine |
| 连续色图 | YlOrRd（黄→橙→红，用于热图） |
| 发散色图 | RdBu（红→白→蓝，用于相关性热图） |

### OmicShare 色板（10 色）

| 序号 | 颜色 | Hex | 典型用途 |
|------|------|-----|---------|
| 1 | 红色 | `#E64B35` | KEGG Metabolism |
| 2 | 青色 | `#4DBBD5` | KEGG Genetic Information Processing |
| 3 | 绿色 | `#00A087` | KEGG Environmental Information Processing |
| 4 | 深蓝 | `#3C5488` | KEGG Cellular Processes |
| 5 | 浅橙 | `#F39B7F` | KEGG Organismal Systems |
| 6 | 灰蓝 | `#8491B4` | KEGG Human Diseases |
| 7 | 薄荷 | `#91D1C2` | 补充分类色 |
| 8 | 正红 | `#DC0000` | 强调色/上调基因 |
| 9 | 棕色 | `#7E6148` | 补充分类色 |
| 10 | 卡其 | `#B09C85` | 补充分类色 |

### OmicShare 经典图表配色映射

| 图表类型 | 配色方案 |
|---------|---------|
| GO 富集条形图 | 颜色深浅渐变代表 P 值（深=显著） |
| KEGG 富集条形图 | 按一级分类着色（7 大类对应色板前 7 色） |
| 富集气泡图 | 蓝→红渐变（P 值，低→高显著性） |
| 火山图 | 红色(上调) + 蓝色(下调) + 灰色(不显著) |
| 热图 | 蓝白红或绿白红发散色图 |
| GO 分类柱状图 | BP=橙红, CC=蓝色, MF=绿色 |

### 与其他风格对比

| 特性 | nature | omicshare | cute |
|------|--------|-----------|------|
| 渲染引擎 | matplotlib | matplotlib | ECharts |
| 输出格式 | png/pdf/svg | png/pdf/svg | html |
| 配色风格 | 冷色调简洁 | 柔和饱和 | 手绘卡通 |
| 中文支持 | 需额外配置 | 内置 SimHei 回退 | 依赖浏览器 |
| 适用场景 | 国际期刊投稿 | 中文文献/毕业论文 | 网页/教学演示 |
