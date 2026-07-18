"""Interactive Plotly figures for transcription factor enrichment results."""

import logging
from typing import Dict, Optional

import numpy as np
import pandas as pd
import plotly.graph_objects as go

logger = logging.getLogger(__name__)

# Colors used for transcription-factor regulatory modes.
MODE_COLORS = {
    "activator": "#2ecc71",
    "repressor": "#e74c3c",
    "mixed": "#f39c12",
    "unknown": "#3498db",
}


class Visualizer:
    """Build interactive figures for transcription factor analyses."""

    def plot_tf_enrichment_bar(
        self,
        result_df: pd.DataFrame,
        top_n: int = 20,
        title: str = "Transcription Factor Enrichment",
        color_by_mode: bool = True,
    ) -> go.Figure:
        """Plot the most significant transcription factors as horizontal bars.

        Args:
            result_df: Enrichment results containing TF names and P values.
            top_n: Maximum number of transcription factors to display.
            title: Figure title.
            color_by_mode: Color bars by the ``Mode`` column when available.

        Returns:
            An interactive Plotly figure.
        """
        df = result_df.copy()

        # Resolve supported column aliases without changing the input table.
        pval_col = self._find_column(df, ["Pvalue", "P_Value", "pvalue", "P-value", "p_value"])
        overlap_col = self._find_column(df, ["Overlap", "overlap", "Gene_Count", "gene_count"])
        tf_col = self._find_column(
            df, ["Term_Name", "term_name", "TF", "tf", "Name", "name"]
        )

        if tf_col is None or pval_col is None:
            raise ValueError(
                "result_df must contain TF name and P-value columns. "
                "Expected columns: TF/Name/Term_Name and Pvalue/P_Value"
            )

        # Select the most significant rows.
        df = df.sort_values(by=pval_col, ascending=True).head(top_n).copy()

        # Transform P values for display while protecting against log10(0).
        df["_neg_log10_pval"] = -np.log10(df[pval_col].astype(float).clip(lower=1e-300))

        # Use regulatory-mode colors when that annotation is available.
        if color_by_mode and "Mode" in df.columns:
            bar_colors = df["Mode"].map(MODE_COLORS).fillna(MODE_COLORS["unknown"]).tolist()
        else:
            bar_colors = ["#3498db"] * len(df)

        # Draw the primary bar layer.
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

        # Label bars with overlap counts when provided.
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

        # Keep long TF names readable without crowding the plot.
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
        """Plot the regulatory-mode distribution among significant TFs.

        When an FDR column is available, only rows with FDR below 0.05 are
        included.

        Args:
            result_df: Enrichment results containing ``Mode`` and optionally FDR.
            title: Figure title.

        Returns:
            An interactive Plotly figure.
        """
        df = result_df.copy()

        # Resolve the adjusted-P-value column when present.
        fdr_col = self._find_column(df, ["FDR", "fdr", "Adjusted_P_Value", "adjusted_p_value", "Qvalue"])

        # Restrict the summary to statistically significant TFs.
        if fdr_col is not None:
            df = df[df[fdr_col].astype(float) < 0.05].copy()

        if "Mode" not in df.columns or len(df) == 0:
            # Render an explicit empty state instead of a misleading distribution.
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

        # Count TFs in each regulatory mode.
        mode_counts = df["Mode"].value_counts()

        # Apply the established color mapping to known modes.
        labels = []
        values = []
        colors = []
        for mode in ["activator", "repressor", "mixed", "unknown"]:
            if mode in mode_counts.index:
                labels.append(mode.capitalize())
                values.append(int(mode_counts[mode]))
                colors.append(MODE_COLORS[mode])

        # Preserve unexpected mode labels with a neutral fallback color.
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
        """Plot pairwise Jaccard similarity between top TF target sets.

        Args:
            result_df: TF enrichment results used to select the top TFs.
            tf_to_targets: Mapping from each TF name to its target-gene set.
            top_n: Maximum number of TFs to include.
            title: Figure title.

        Returns:
            An interactive Plotly heatmap.
        """
        # Select top TFs using P values when available.
        pval_col = self._find_column(result_df, ["Pvalue", "P_Value", "pvalue", "P-value", "p_value"])
        tf_col = self._find_column(result_df, ["TF", "tf", "Term_Name", "term_name", "Name", "name"])

        if tf_col is None:
            raise ValueError("result_df must contain a TF name column (TF/Name/Term_Name)")

        if pval_col is not None:
            top_tfs = result_df.sort_values(by=pval_col, ascending=True).head(top_n)[tf_col].tolist()
        else:
            top_tfs = result_df.head(top_n)[tf_col].tolist()

        # Keep only TFs for which target-gene data are available.
        top_tfs = [tf for tf in top_tfs if tf in tf_to_targets]

        if len(top_tfs) < 2:
            raise ValueError(
                f"Need at least 2 TFs with target data for heatmap, got {len(top_tfs)}. "
                "Ensure tf_to_targets contains entries for the top TFs."
            )

        # Calculate the pairwise Jaccard similarity matrix.
        n = len(top_tfs)
        jaccard_matrix = np.zeros((n, n))

        for i in range(n):
            for j in range(n):
                set_i = tf_to_targets[top_tfs[i]]
                set_j = tf_to_targets[top_tfs[j]]
                intersection = len(set_i & set_j)
                union = len(set_i | set_j)
                jaccard_matrix[i, j] = intersection / union if union > 0 else 0.0

        # Build the annotated similarity heatmap.
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
        """Return the first candidate column present in a DataFrame.

        Args:
            df: DataFrame to inspect.
            candidates: Candidate column names in priority order.

        Returns:
            The matching column name, or ``None`` when no candidate exists.
        """
        for col in candidates:
            if col in df.columns:
                return col
        return None
