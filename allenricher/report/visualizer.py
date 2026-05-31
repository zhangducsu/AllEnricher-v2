"""
TF富集分析可视化模块 - 基于 plotly.graph_objects 的交互式图表

提供转录因子(Transcription Factor)富集分析结果的交互式可视化图表：
- 水平条形图：TF富集显著性排名
- 饼图：TF调控模式分布
- 热图：TF-Target重叠度(Jaccard相似系数)
"""

import logging
from typing import Dict, Optional

import numpy as np
import pandas as pd
import plotly.graph_objects as go

logger = logging.getLogger(__name__)

# TF调控模式配色方案
MODE_COLORS = {
    "activator": "#2ecc71",   # 绿色 - 激活因子
    "repressor": "#e74c3c",    # 红色 - 抑制因子
    "mixed": "#f39c12",        # 橙色 - 混合调控
    "unknown": "#3498db",      # 蓝色 - 未知模式
}


class Visualizer:
    """TF富集分析可视化器 - 基于 plotly 的交互式图表"""

    def plot_tf_enrichment_bar(
        self,
        result_df: pd.DataFrame,
        top_n: int = 20,
        title: str = "Transcription Factor Enrichment",
        color_by_mode: bool = True,
    ) -> go.Figure:
        """生成TF富集分析水平条形图

        Args:
            result_df: TF富集分析结果DataFrame，需包含 TF名称、Pvalue/FDR、Overlap 等列
            top_n: 显示前N个最显著的TF
            title: 图表标题
            color_by_mode: 是否按调控模式(Mode列)着色

        Returns:
            go.Figure: plotly交互式水平条形图
        """
        df = result_df.copy()

        # 标准化列名映射
        pval_col = self._find_column(df, ["Pvalue", "P_Value", "pvalue", "P-value", "p_value"])
        overlap_col = self._find_column(df, ["Overlap", "overlap", "Gene_Count", "gene_count"])
        tf_col = self._find_column(df, ["TF", "tf", "Term_Name", "term_name", "Name", "name"])

        if tf_col is None or pval_col is None:
            raise ValueError(
                "result_df must contain TF name and P-value columns. "
                "Expected columns: TF/Name/Term_Name and Pvalue/P_Value"
            )

        # 按 Pvalue 升序排列，取 top_n
        df = df.sort_values(by=pval_col, ascending=True).head(top_n).copy()

        # 计算 -log10(Pvalue)
        df["_neg_log10_pval"] = -np.log10(df[pval_col].astype(float).clip(lower=1e-300))

        # 确定条形颜色
        if color_by_mode and "Mode" in df.columns:
            bar_colors = df["Mode"].map(MODE_COLORS).fillna(MODE_COLORS["unknown"]).tolist()
        else:
            bar_colors = ["#3498db"] * len(df)

        # 构建条形图
        fig = go.Figure()

        fig.add_trace(
            go.Bar(
                y=df[tf_col].tolist(),
                x=df["_neg_log10_pval"].tolist(),
                orientation="h",
                marker_color=bar_colors,
                textposition="outside",
                textfont=dict(size=10, color="#333333"),
                hovertemplate=(
                    "<b>%{y}</b><br>"
                    "-log10(Pvalue): %{x:.2f}<br>"
                    "<extra></extra>"
                ),
            )
        )

        # 条形末端显示 Overlap 数值
        if overlap_col is not None:
            overlap_vals = df[overlap_col].tolist()
            fig.add_trace(
                go.Scatter(
                    x=df["_neg_log10_pval"].tolist(),
                    y=df[tf_col].tolist(),
                    mode="text",
                    text=[str(v) for v in overlap_vals],
                    textposition="middle right",
                    textfont=dict(size=9, color="#555555"),
                    hoverinfo="skip",
                    showlegend=False,
                )
            )

        # 布局设置
        fig.update_layout(
            title=dict(text=title, font=dict(size=16)),
            xaxis_title="-log10(Pvalue)",
            yaxis=dict(
                title="Transcription Factor",
                categoryorder="total ascending",
                automargin=True,
            ),
            height=400 + top_n * 15,
            margin=dict(l=180, r=80, t=50, b=50),
            showlegend=False,
            template="plotly_white",
        )

        return fig

    def plot_tf_mode_pie(
        self,
        result_df: pd.DataFrame,
        title: str = "TF Regulation Mode Distribution",
    ) -> go.Figure:
        """生成TF调控模式饼图

        仅统计显著TF（FDR < 0.05）的调控模式分布。

        Args:
            result_df: TF富集分析结果DataFrame，需包含 Mode 和 FDR 列
            title: 图表标题

        Returns:
            go.Figure: plotly交互式饼图
        """
        df = result_df.copy()

        # 查找FDR列
        fdr_col = self._find_column(df, ["FDR", "fdr", "Adjusted_P_Value", "adjusted_p_value", "Qvalue"])

        # 筛选显著TF
        if fdr_col is not None:
            df = df[df[fdr_col].astype(float) < 0.05].copy()

        if "Mode" not in df.columns or len(df) == 0:
            # 无显著结果时返回空饼图
            fig = go.Figure(
                data=[
                    go.Pie(
                        labels=["No significant TFs (FDR < 0.05)"],
                        values=[1],
                        hole=0.4,
                    )
                ]
            )
            fig.update_layout(title=dict(text=title, font=dict(size=16)))
            return fig

        # 统计各模式数量
        mode_counts = df["Mode"].value_counts()

        # 确保所有已知模式都有对应颜色
        labels = []
        values = []
        colors = []
        for mode in ["activator", "repressor", "mixed", "unknown"]:
            if mode in mode_counts.index:
                labels.append(mode.capitalize())
                values.append(int(mode_counts[mode]))
                colors.append(MODE_COLORS[mode])

        # 处理未知模式名
        for mode_name, count in mode_counts.items():
            if mode_name not in MODE_COLORS:
                labels.append(str(mode_name))
                values.append(int(count))
                colors.append("#95a5a6")

        fig = go.Figure(
            data=[
                go.Pie(
                    labels=labels,
                    values=values,
                    marker_colors=colors,
                    hole=0.4,
                    textinfo="label+percent",
                    textposition="outside",
                    textfont=dict(size=12),
                    hovertemplate="<b>%{label}</b><br>Count: %{value}<br>Percent: %{percent}<extra></extra>",
                )
            ]
        )

        fig.update_layout(
            title=dict(text=title, font=dict(size=16)),
            showlegend=True,
            legend=dict(orientation="h", yanchor="bottom", y=-0.1),
            margin=dict(t=50, b=80),
            template="plotly_white",
        )

        return fig

    def plot_tf_overlap_heatmap(
        self,
        result_df: pd.DataFrame,
        tf_to_targets: Dict[str, set],
        top_n: int = 15,
        title: str = "TF-Target Overlap Heatmap",
    ) -> go.Figure:
        """生成TF-Target重叠度热图(Jaccard相似系数)

        Args:
            result_df: TF富集分析结果DataFrame（用于确定top TF）
            tf_to_targets: TF到其靶基因集合的映射 {TF_name: set_of_genes}
            top_n: 显示前N个最显著的TF
            title: 图表标题

        Returns:
            go.Figure: plotly交互式热图
        """
        # 确定top TF列表
        pval_col = self._find_column(result_df, ["Pvalue", "P_Value", "pvalue", "P-value", "p_value"])
        tf_col = self._find_column(result_df, ["TF", "tf", "Term_Name", "term_name", "Name", "name"])

        if tf_col is None:
            raise ValueError("result_df must contain a TF name column (TF/Name/Term_Name)")

        if pval_col is not None:
            top_tfs = result_df.sort_values(by=pval_col, ascending=True).head(top_n)[tf_col].tolist()
        else:
            top_tfs = result_df.head(top_n)[tf_col].tolist()

        # 过滤出有靶基因数据的TF
        top_tfs = [tf for tf in top_tfs if tf in tf_to_targets]

        if len(top_tfs) < 2:
            raise ValueError(
                f"Need at least 2 TFs with target data for heatmap, got {len(top_tfs)}. "
                "Ensure tf_to_targets contains entries for the top TFs."
            )

        # 计算 Jaccard 相似系数矩阵
        n = len(top_tfs)
        jaccard_matrix = np.zeros((n, n))

        for i in range(n):
            for j in range(n):
                set_i = tf_to_targets[top_tfs[i]]
                set_j = tf_to_targets[top_tfs[j]]
                intersection = len(set_i & set_j)
                union = len(set_i | set_j)
                jaccard_matrix[i, j] = intersection / union if union > 0 else 0.0

        # 构建热图
        fig = go.Figure(
            data=go.Heatmap(
                z=jaccard_matrix,
                x=top_tfs,
                y=top_tfs,
                colorscale="YlOrRd",
                zmin=0,
                zmax=1,
                text=jaccard_matrix,
                texttemplate="%{text:.2f}",
                textfont=dict(size=8),
                hovertemplate=(
                    "<b>%{y} vs %{x}</b><br>"
                    "Jaccard: %{z:.3f}<br>"
                    "<extra></extra>"
                ),
                colorbar=dict(title="Jaccard Similarity"),
            )
        )

        fig.update_layout(
            title=dict(text=title, font=dict(size=16)),
            xaxis=dict(tickangle=45, side="bottom"),
            yaxis=dict(autorange="reversed"),
            height=max(500, 300 + n * 25),
            width=max(600, 400 + n * 30),
            margin=dict(l=150, r=80, t=50, b=150),
            template="plotly_white",
        )

        return fig

    @staticmethod
    def _find_column(df: pd.DataFrame, candidates: list) -> Optional[str]:
        """在DataFrame中查找第一个匹配的列名

        Args:
            df: 目标DataFrame
            candidates: 候选列名列表

        Returns:
            匹配到的列名，未找到返回None
        """
        for col in candidates:
            if col in df.columns:
                return col
        return None
