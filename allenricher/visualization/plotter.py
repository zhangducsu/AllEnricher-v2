"""
可视化模块 - AllEnricher v2.0
=============================

本模块负责富集分析结果的可视化展示。
- barplot.py: Python matplotlib 绘图
- bubble.py: Python matplotlib 绘图

依赖：
- matplotlib, pandas, numpy
"""

import logging
from pathlib import Path
from typing import Dict, List, Optional

import pandas as pd

logger = logging.getLogger(__name__)


class Plotter:
    """可视化绘图器类"""

    def __init__(self, output_dir: str, config=None):
        """初始化绘图器"""
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.config = config

    def _prepare_barplot_data(self, data: pd.DataFrame) -> pd.DataFrame:
        """准备 barplot.py 所需的数据格式

        Args:
            data: 原始富集结果 DataFrame

        Returns:
            格式化后的 DataFrame，包含以下列：
            - term: Term 名称
            - qvalue: 校正后 P 值
            - gene_count: 基因数量
            - rich_factor: 富集因子（可选，GSEA 时使用 NES 替代）
        """
        df = data.copy()
        result = pd.DataFrame()

        # 检测是否为 GSEA 结果
        is_gsea = 'NES' in df.columns or 'nes' in df.columns

        # 映射常见列名到标准列名
        column_mappings = {
            'term': ['Term_Name', 'term_name', 'Term', 'term', 'Description', 'description'],
            'qvalue': ['Adjusted_P_Value', 'adjusted_p_value', 'adjP', 'qvalue', 'Qvalue', 'FDR', 'fdr'],
            'gene_count': ['Gene_Count', 'gene_count', 'setSize', 'ObservedGeneNum', 'Count', 'count'],
        }
        if not is_gsea:
            column_mappings['rich_factor'] = ['Rich_Factor', 'rich_factor', 'RichFactor', 'richfactor']

        for std_col, possible_cols in column_mappings.items():
            for col in possible_cols:
                if col in df.columns:
                    result[std_col] = df[col]
                    break

        if is_gsea:
            # GSEA：用 NES 作为数值轴，替代 Rich_Factor
            for col in ['NES', 'nes', 'ES', 'es']:
                if col in df.columns:
                    result['rich_factor'] = df[col]
                    break
        else:
            # 如果没有找到 rich_factor，尝试计算
            if 'rich_factor' not in result.columns and 'gene_count' in result.columns:
                bg_col = None
                for col in ['Background_Count', 'background_count', 'TermGeneNum']:
                    if col in df.columns:
                        bg_col = df[col]
                        break
                if bg_col is not None:
                    result['rich_factor'] = result['gene_count'] / bg_col

        return result

    def _prepare_bubble_data(self, data: pd.DataFrame) -> pd.DataFrame:
        """准备 bubble.py 所需的数据格式

        Args:
            data: 原始富集结果 DataFrame

        Returns:
            格式化后的 DataFrame，包含以下列：
            - Term_Name: Term 名称
            - RichFactor: 富集因子（GSEA 时使用 NES 替代）
            - GeneCount: 基因数量
            - Qvalue: 校正后 P 值
        """
        df = data.copy()
        result = pd.DataFrame()

        is_gsea = 'NES' in df.columns or 'nes' in df.columns

        # 映射常见列名到标准列名
        column_mappings = {
            'Term_Name': ['Term_Name', 'term_name', 'Term', 'term', 'Description', 'description'],
            'GeneCount': ['Gene_Count', 'gene_count', 'setSize', 'ObservedGeneNum', 'Count', 'count'],
            'Qvalue': ['Adjusted_P_Value', 'adjusted_p_value', 'adjP', 'qvalue', 'Qvalue', 'FDR', 'fdr'],
        }
        if not is_gsea:
            column_mappings['RichFactor'] = ['Rich_Factor', 'rich_factor', 'RichFactor', 'richfactor', 'Rich Factor']

        for std_col, possible_cols in column_mappings.items():
            for col in possible_cols:
                if col in df.columns:
                    result[std_col] = df[col]
                    break

        if is_gsea:
            # GSEA：用 NES 作为气泡 X 轴，替代 RichFactor
            for col in ['NES', 'nes', 'ES', 'es']:
                if col in df.columns:
                    result['RichFactor'] = df[col]
                    break
        else:
            # 如果没有 RichFactor，尝试计算
            if 'RichFactor' not in result.columns and 'GeneCount' in result.columns:
                bg_col = None
                for col in ['Background_Count', 'background_count', 'TermGeneNum']:
                    if col in df.columns:
                        bg_col = df[col]
                        break
                if bg_col is not None:
                    result['RichFactor'] = result['GeneCount'] / bg_col

        return result

    def plot_barplot(
        self,
        data: pd.DataFrame,
        database: str,
        output_file: str,
        top_n: int = 20,
        style: Optional[str] = None,
        palette: Optional[str] = None,
    ) -> str:
        """生成富集分析柱状图

        Args:
            data: 富集结果 DataFrame
            database: 数据库名称 (GO, KEGG, Reactome, DO, DisGeNET)
            output_file: 输出文件名
            top_n: 显示前 N 条数据
            style: 图表风格 (nature, science, presentation, colorblind, omicshare)
            palette: 色板名称

        Returns:
            输出文件路径
        """
        # 构建输出路径
        output_path = self.output_dir / output_file

        from .barplot import plot_barplot as _plot_barplot

        # 准备数据
        plot_data = self._prepare_barplot_data(data)

        # 获取 DPI 设置
        dpi = 300
        if self.config and hasattr(self.config, 'figure_dpi'):
            dpi = self.config.figure_dpi

        try:
            _plot_barplot(
                data=plot_data,
                database=database,
                top_n=top_n,
                style=style or 'default',
                palette=palette or 'default',
                output_file=str(output_path),
                dpi=dpi,
            )
        except Exception as e:
            logger.warning(f"柱状图生成失败: {e}")

        return str(output_path)

    def plot_bubble(
        self,
        data: pd.DataFrame,
        output_file: str,
        database: str = 'GO',
        top_n: int = 20,
        style: Optional[str] = None,
        palette: Optional[str] = None,
    ) -> str:
        """生成富集分析气泡图

        Args:
            data: 富集结果 DataFrame
            output_file: 输出文件名
            database: 数据库名称 (GO, KEGG, Reactome, DO, DisGeNET)
            top_n: 显示前 N 条数据
            style: 图表风格
            palette: 色板名称

        Returns:
            输出文件路径
        """
        # 构建输出路径
        output_path = self.output_dir / output_file

        from .bubble import plot_bubble as _plot_bubble

        # 准备数据
        plot_data = self._prepare_bubble_data(data)

        # 获取 DPI 设置
        dpi = 300
        if self.config and hasattr(self.config, 'figure_dpi'):
            dpi = self.config.figure_dpi

        try:
            fig = _plot_bubble(
                data=plot_data,
                top_n=top_n,
                style=style,
                palette=palette,
                title=f"{database} Enrichment",
            )
            if fig:
                # 使用 save_figure_dual 同时生成 PNG 和 PDF
                from .plot_theme import save_figure_dual
                save_figure_dual(fig, str(output_path), dpi=dpi)
                import matplotlib.pyplot as plt
                plt.close(fig)
        except Exception as e:
            logger.warning(f"气泡图生成失败: {e}")

        return str(output_path)

    def plot_all(
        self,
        data: pd.DataFrame,
        database: str,
        top_n: int = 20,
        style: Optional[str] = None,
        palette: Optional[str] = None,
    ) -> Dict[str, str]:
        """生成所有标准图表

        Args:
            data: 富集结果 DataFrame
            database: 数据库名称
            top_n: 显示前 N 条数据
            style: 图表风格
            palette: 色板名称

        Returns:
            包含图表路径的字典
        """
        plots = {}

        # 根据配置确定输出格式
        fmt = 'png'
        if self.config and hasattr(self.config, 'figure_format'):
            fmt = self.config.figure_format

        # 生成柱状图
        bar_file = f"{database}_barplot.{fmt}"
        bar_path = self.plot_barplot(
            data, database, bar_file, top_n,
            style=style, palette=palette
        )
        plots["barplot"] = bar_path

        # 生成气泡图
        bubble_file = f"{database}_bubble.{fmt}"
        bubble_path = self.plot_bubble(
            data, bubble_file, database, top_n,
            style=style, palette=palette
        )
        plots["bubble"] = bubble_path

        return plots
