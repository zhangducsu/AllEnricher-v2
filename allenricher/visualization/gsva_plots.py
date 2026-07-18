"""Create activity heatmaps, group comparisons, and sample-correlation figures for ssGSEA and GSVA."""

import itertools
import logging
import math
import textwrap
from pathlib import Path
from typing import Dict, List, Optional

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.cm import ScalarMappable
from matplotlib.colors import LinearSegmentedColormap, ListedColormap, Normalize
from matplotlib.patches import FancyBboxPatch, Polygon, Rectangle
from scipy.cluster.hierarchy import dendrogram, leaves_list, linkage
from scipy.spatial.distance import pdist, squareform
from scipy import stats
from statsmodels.stats.multitest import multipletests

from allenricher.visualization.plot_theme import (
    PlotTheme,
    apply_figure_style,
    resolve_style,
    save_figure,
)

logger = logging.getLogger(__name__)

_ACTIVITY_COLORS = [
    "#3F0D5E", "#8A5AA8", "#D8CBE2", "#F7F7F7",
    "#D7ECD0", "#7DBB73", "#1B5E2A",
]

_ANNOTATION_COLORS = [
    "#4E79A7", "#F28E2B", "#59A14F", "#E15759",
    "#B07AA1", "#76B7B2", "#EDC948", "#FF9DA7",
]

DEFAULT_ACTIVITY_HEATMAP_TOP_N = 40


def _clean_pathway_label(value: object) -> str:
    label = str(value).rsplit("|", 1)[-1].strip().replace("_", " ")
    return " ".join(label.split())


def _wrap_panel_label(value: object, width: int = 32) -> str:
    label = _clean_pathway_label(value)
    lines = textwrap.wrap(label, width=width, break_long_words=False, break_on_hyphens=False)
    if len(lines) <= 2:
        return "\n".join(lines)
    remainder = textwrap.shorten(" ".join(lines[1:]), width=width, placeholder="...")
    return f"{lines[0]}\n{remainder}"


def select_activity_heatmap_scores(
    scores_df: pd.DataFrame,
    annotation_col: pd.DataFrame = None,
    top_n: int = DEFAULT_ACTIVITY_HEATMAP_TOP_N,
) -> pd.DataFrame:
    """Select informative display rows without altering the saved activity matrix."""
    if top_n <= 0 or len(scores_df) <= top_n:
        return scores_df

    numeric = scores_df.apply(pd.to_numeric, errors="coerce")
    variance = numeric.var(axis=1, skipna=True).fillna(0.0)
    display_score = variance
    if annotation_col is not None and not annotation_col.empty:
        group_column = "Group" if "Group" in annotation_col.columns else annotation_col.columns[0]
        labels = annotation_col[group_column].reindex(numeric.columns)
        groups = [value for value in dict.fromkeys(labels.dropna().astype(str))]
        if len(groups) >= 2:
            means = pd.concat(
                [numeric.loc[:, labels.astype(str).eq(group)].mean(axis=1) for group in groups],
                axis=1,
            )
            display_score = (means.max(axis=1) - means.min(axis=1)).fillna(0.0)

    ranking = pd.DataFrame({
        "score": display_score.to_numpy(),
        "variance": variance.to_numpy(),
        "label": scores_df.index.map(str),
        "position": np.arange(len(scores_df)),
    }).sort_values(
        ["score", "variance", "label"],
        ascending=[False, False, True],
        kind="mergesort",
    )
    return scores_df.iloc[ranking.head(top_n)["position"].to_numpy()]


def _activity_heatmap_layout(
    n_pathways: int,
    n_samples: int,
    max_label_chars: int,
    n_annotations: int,
) -> dict:
    """Return bounded dimensions and label density for an activity heatmap."""
    heat_width = float(np.clip(0.23 * n_samples, 3.2, 9.2))
    heat_height = float(np.clip(0.28 * n_pathways, 3.0, 9.0))
    label_width = float(np.clip(0.055 * max_label_chars, 1.8, 4.2))
    annotation_height = 0.14 * max(1, n_annotations)
    figure_width = float(np.clip(1.2 + 0.70 + heat_width + label_width + 1.25, 8.4, 15.0))
    figure_height = float(np.clip(1.55 + 0.55 + annotation_height + heat_height, 5.2, 12.0))

    sample_capacity = max(1, int(heat_width * 72 / 8.0))
    pathway_capacity = max(1, int(heat_height * 72 / 8.0))
    sample_stride = max(1, math.ceil(n_samples / sample_capacity))
    pathway_stride = max(1, math.ceil(n_pathways / pathway_capacity))
    sample_fontsize = float(np.clip(heat_width * 72 / max(1, n_samples) * 0.55, 4.5, 8.5))
    pathway_fontsize = float(np.clip(heat_height * 72 / max(1, n_pathways) * 0.55, 4.5, 9.2))
    return {
        "figsize": (figure_width, figure_height),
        "heat_width": heat_width,
        "heat_height": heat_height,
        "label_width": label_width,
        "annotation_height": annotation_height,
        "sample_stride": sample_stride,
        "pathway_stride": pathway_stride,
        "sample_fontsize": sample_fontsize,
        "pathway_fontsize": pathway_fontsize,
    }

def _activity_linkage(values: np.ndarray) -> np.ndarray:
    distances = pdist(values, metric="correlation")
    finite = distances[np.isfinite(distances)]
    distances[~np.isfinite(distances)] = finite.max() if finite.size else 1.0
    return linkage(np.maximum(distances, 0), method="average", optimal_ordering=True)


def plot_pathway_heatmap(
    scores_df: pd.DataFrame,
    annotation_col: pd.DataFrame = None,
    cluster_rows: bool = True,
    cluster_cols: bool = True,
    cmap: str = None,
    center: float = 0,
    title: str = "Pathway Activity Heatmap",
    output_file: str = None,
    figsize: tuple = None,
    dpi: int = 300,
    style: str = "nature",
    palette: str = None,
    scale: str = "row",
    top_n: int = DEFAULT_ACTIVITY_HEATMAP_TOP_N,
) -> matplotlib.figure.Figure:
    """Plot a clustered pathway-by-sample activity matrix with rounded cells."""
    if scale not in {"row", "column", "none"}:
        raise ValueError("Scale must be a row, a column or none")

    scores_df = select_activity_heatmap_scores(scores_df, annotation_col, top_n)
    scores = scores_df.apply(pd.to_numeric, errors="coerce").copy()
    values = scores.to_numpy(dtype=float, copy=True)
    for row in range(values.shape[0]):
        finite = np.isfinite(values[row])
        values[row, ~finite] = np.median(values[row, finite]) if finite.any() else 0

    if scale != "none":
        axis = 1 if scale == "row" else 0
        means = values.mean(axis=axis, keepdims=True)
        stds = values.std(axis=axis, ddof=1, keepdims=True)
        stds[~np.isfinite(stds) | (stds == 0)] = 1
        values = (values - means) / stds
    scores = pd.DataFrame(values, index=scores.index, columns=scores.columns)

    row_link = _activity_linkage(values) if cluster_rows and len(scores) > 1 else None
    col_link = _activity_linkage(values.T) if cluster_cols and scores.shape[1] > 1 else None
    row_order = leaves_list(row_link) if row_link is not None else np.arange(len(scores))
    col_order = leaves_list(col_link) if col_link is not None else np.arange(scores.shape[1])
    scores = scores.iloc[row_order, col_order]

    metadata = pd.DataFrame(index=scores.columns)
    if annotation_col is not None and not annotation_col.empty:
        metadata = annotation_col.reindex(scores.columns).dropna(axis=1, how="all").astype(str)

    limit = float(np.quantile(np.abs(scores.to_numpy()), 0.98))
    if not np.isfinite(limit) or limit <= 0:
        limit = 1.0
    if cmap is not None:
        heat_cmap = plt.get_cmap(cmap)
    else:
        heat_colors = PlotTheme.get_plot_colors(
            style, palette, default=_ACTIVITY_COLORS, divergent=True
        )
        heat_cmap = LinearSegmentedColormap.from_list("activity", heat_colors)
    norm = Normalize(-limit, limit)
    category_count = sum(metadata[column].nunique() for column in metadata.columns)
    group_colors = PlotTheme.get_plot_colors(
        style, palette, default=_ANNOTATION_COLORS, n=max(1, category_count)
    )

    mappings = {}
    offset = 0
    for column in metadata.columns:
        categories = list(dict.fromkeys(metadata[column]))
        mappings[column] = {
            category: group_colors[offset + index]
            for index, category in enumerate(categories)
        }
        offset += len(categories)

    n_pathways, n_samples = scores.shape
    n_annotations = max(1, metadata.shape[1])
    max_label_chars = max((len(_clean_pathway_label(value)) for value in scores.index), default=1)
    layout = _activity_heatmap_layout(
        n_pathways, n_samples, max_label_chars, n_annotations
    )
    if figsize is None:
        figsize = layout["figsize"]

    with PlotTheme.context(style or "nature", palette):
        plt.rcParams["font.family"] = "sans-serif"
        plt.rcParams["font.sans-serif"] = ["Arial", "DejaVu Sans"]
        fig = plt.figure(figsize=figsize, dpi=dpi, facecolor="white")
        grid = fig.add_gridspec(
            3, 4,
            width_ratios=[0.70, layout["heat_width"], layout["label_width"], 1.25],
            height_ratios=[0.55, layout["annotation_height"], layout["heat_height"]],
            left=0.025, right=0.985, top=0.915, bottom=0.17,
            wspace=0.025, hspace=0.02,
        )
        ax_col = fig.add_subplot(grid[0, 1])
        ax_annotation = fig.add_subplot(grid[1, 1])
        ax_row = fig.add_subplot(grid[2, 0])
        ax_heat = fig.add_subplot(grid[2, 1])
        ax_labels = fig.add_subplot(grid[2, 2])
        ax_legend = fig.add_subplot(grid[:, 3], label="activity_legend")

        if col_link is not None:
            dendrogram(col_link, ax=ax_col, no_labels=True, color_threshold=0,
                       above_threshold_color="black", link_color_func=lambda _: "black")
        ax_col.axis("off")
        if row_link is not None:
            dendrogram(row_link, ax=ax_row, orientation="left", no_labels=True,
                       color_threshold=0, above_threshold_color="black",
                       link_color_func=lambda _: "black")
        ax_row.axis("off")

        if metadata.shape[1]:
            ax_annotation.set_xlim(-0.5, n_samples - 0.5)
            ax_annotation.set_ylim(metadata.shape[1] - 0.5, -0.5)
            for index, column in enumerate(metadata.columns):
                categories = list(mappings[column])
                codes = np.array([categories.index(value) for value in metadata[column]])[None, :]
                ax_annotation.imshow(
                    codes, aspect="auto", interpolation="nearest",
                    cmap=ListedColormap([mappings[column][value] for value in categories]),
                    vmin=-0.5, vmax=max(0.5, len(categories) - 0.5),
                    extent=[-0.5, n_samples - 0.5, index + 0.5, index - 0.5],
                )
        ax_annotation.set_xticks([])
        ax_annotation.set_yticks([])
        for spine in ax_annotation.spines.values():
            spine.set_visible(False)

        image = ScalarMappable(norm=norm, cmap=heat_cmap)
        image.set_array(scores.to_numpy())
        for row, values in enumerate(scores.to_numpy()):
            for column, value in enumerate(values):
                ax_heat.add_patch(FancyBboxPatch(
                    (column - 0.46, row - 0.44), 0.92, 0.88,
                    boxstyle="round,pad=0,rounding_size=0.10",
                    facecolor=heat_cmap(norm(value)), edgecolor="none", zorder=3,
                ))
        ax_heat.set_xlim(-0.5, n_samples - 0.5)
        ax_heat.set_ylim(n_pathways - 0.5, -0.5)
        sample_ticks = np.arange(0, n_samples, layout["sample_stride"])
        ax_heat.set_xticks(
            sample_ticks,
            scores.columns[sample_ticks],
            rotation=90,
            fontsize=layout["sample_fontsize"],
        )
        ax_heat.set_yticks([])
        ax_heat.set_xticks(np.arange(-0.5, n_samples, 1), minor=True)
        ax_heat.set_yticks(np.arange(-0.5, n_pathways, 1), minor=True)
        ax_heat.grid(which="minor", color="#F7F7F7", linewidth=0.5, zorder=1)
        ax_heat.tick_params(which="both", left=False, bottom=False)
        for spine in ax_heat.spines.values():
            spine.set_visible(True)
            spine.set_color("#333333")
            spine.set_linewidth(0.7)

        ax_labels.set_xlim(0, 1)
        ax_labels.set_ylim(n_pathways - 0.5, -0.5)
        ax_labels.axis("off")
        for index in range(0, n_pathways, layout["pathway_stride"]):
            pathway = scores.index[index]
            ax_labels.text(
                0.02, index, _clean_pathway_label(pathway),
                ha="left", va="center", fontsize=layout["pathway_fontsize"],
            )

        ax_legend.axis("off")
        legend_rows = sum(1 + len(mappings[column]) for column in metadata.columns)
        figure_height = fig.get_size_inches()[1]
        legend_step = min(0.045, 0.24 / figure_height, 0.54 / max(1, legend_rows))
        heading_gap = min(0.05, 0.28 / figure_height)
        section_gap = min(0.06, 0.34 / figure_height)
        legend_y = 0.94
        for column in metadata.columns:
            ax_legend.text(0, legend_y, column, ha="left", va="top",
                           fontsize=8.8, fontweight="bold")
            categories = list(mappings[column].items())
            legend_y -= heading_gap
            for index, (label, color) in enumerate(categories):
                y_pos = legend_y - legend_step * index
                ax_legend.add_patch(Rectangle(
                    (0.02, y_pos - 0.009), 0.075, 0.018,
                    transform=ax_legend.transAxes, facecolor=color, edgecolor="none",
                ))
                ax_legend.text(0.12, y_pos, label, ha="left", va="center", fontsize=8.1)
            legend_y -= legend_step * max(0, len(categories) - 1) + section_gap

        color_title = "Activity score" + ("\n(row Z-score)" if scale == "row" else "")
        ax_legend.text(0, legend_y, color_title, ha="left", va="top",
                       fontsize=8.6, fontweight="bold")
        color_top = legend_y - ((0.42 if scale == "row" else 0.24) / figure_height)
        color_height = min(0.18, 1.35 / figure_height)
        ax_color = ax_legend.inset_axes(
            [0.02, max(0.04, color_top - color_height), 0.10, color_height],
            label="activity_colorbar",
        )
        colorbar = fig.colorbar(image, cax=ax_color, orientation="vertical")
        colorbar.set_ticks([-limit, 0, limit])
        colorbar.set_ticklabels([f"{-limit:.2g}", "0", f"{limit:.2g}"])
        colorbar.ax.tick_params(labelsize=8.2, length=2, pad=2)
        colorbar.ax.get_yticklabels()[0].set_verticalalignment("bottom")
        colorbar.ax.get_yticklabels()[-1].set_verticalalignment("top")
        colorbar.outline.set_linewidth(0.6)
        fig.suptitle(title, fontsize=12.5, fontweight="bold", y=0.972)
        apply_figure_style(fig, style, axes=[ax_heat], border="full")

        if output_file is not None:
            save_figure(fig, output_file, dpi=dpi, facecolor="white")
            logger.info("Activity heatmap saved: %s", output_file)
        return fig


def plot_group_comparison(
    scores_df: pd.DataFrame,
    groups: Dict[str, List[str]],
    pathways: List[str] = None,
    plot_type: str = "box",
    title: str = "Pathway activity between groups",
    output_file: str = None,
    figsize: tuple = None,
    dpi: int = 300,
    style: str = "nature",
    palette: str = None,
    statistics_file: str = None,
    top_n: int = 6,
    global_test: str = "kruskal",
    pairwise_test: str = "mannwhitney",
    p_adjust: str = "BH",
    comparison_mode: str = "all",
    reference_group: str = None,
    ncols: int = 3,
    random_seed: int = 123,
) -> matplotlib.figure.Figure:
    """Plot group activity distributions and write the corresponding statistical table."""
    if plot_type not in ("box", "violin", "bar"):
        raise ValueError(f"Unsupported plot_type: '{plot_type}'")
    if global_test not in ("kruskal", "anova"):
        raise ValueError("global_test must be kruskal or anova")
    if pairwise_test not in ("mannwhitney", "ttest"):
        raise ValueError("pairwise_test must be 'mannwhitney' or 'ttest'")
    if p_adjust.lower() not in ("bh", "holm", "bonferroni", "none"):
        raise ValueError("p_adjust must be 'BH', 'holm', 'bonferroni', or 'none'")
    if comparison_mode not in ("all", "reference", "none"):
        raise ValueError("comparison_mode must be 'all', 'reference', or 'none'")

    scores = scores_df.apply(pd.to_numeric, errors="coerce")
    group_names = list(groups)
    if len(group_names) < 2:
        raise ValueError("Group comparison requires at least two sample groups")
    two_groups = len(group_names) == 2
    grouped_samples = {
        group: [sample for sample in groups[group] if sample in scores.columns]
        for group in group_names
    }
    empty_groups = [group for group, samples in grouped_samples.items() if not samples]
    if empty_groups:
        raise ValueError("These groups contain no matching sample columns: " + ", ".join(empty_groups))

    def grouped_values(pathway):
        row = scores.loc[pathway]
        return {group: row[samples].to_numpy(float) for group, samples in grouped_samples.items()}

    def global_p(values):
        values = [value[np.isfinite(value)] for value in values if np.isfinite(value).any()]
        if len(values) < 2:
            return np.nan
        try:
            test = stats.kruskal(*values) if global_test == "kruskal" else stats.f_oneway(*values)
            return float(test.pvalue)
        except ValueError:
            return np.nan

    def pair_p(left, right):
        left, right = left[np.isfinite(left)], right[np.isfinite(right)]
        if not len(left) or not len(right):
            return np.nan
        test = (stats.mannwhitneyu(left, right, alternative="two-sided")
                if pairwise_test == "mannwhitney"
                else stats.ttest_ind(left, right, equal_var=False, nan_policy="omit"))
        return float(test.pvalue)

    if pathways is None:
        ranking = []
        for pathway in scores.index:
            values = grouped_values(pathway)
            pvalue = (pair_p(values[group_names[0]], values[group_names[1]]) if two_groups
                      else global_p(list(values.values())))
            ranking.append((pathway, pvalue))
        ranking.sort(key=lambda item: (np.inf if not np.isfinite(item[1]) else item[1], str(item[0])))
        pathways = [item[0] for item in ranking[:max(1, min(top_n, len(ranking)))]]
    else:
        missing = [pathway for pathway in pathways if pathway not in scores.index]
        if missing:
            raise ValueError("The activity matrix does not contain these pathways: " + ", ".join(map(str, missing)))
    if not pathways:
        raise ValueError("No pathways are available for group comparison")

    if comparison_mode == "none":
        comparisons = []
    elif comparison_mode == "reference":
        reference = reference_group or group_names[0]
        if reference not in group_names:
            raise ValueError(f"Reference group not found: {reference}")
        comparisons = [(reference, group) for group in group_names if group != reference]
    else:
        comparisons = list(itertools.combinations(group_names, 2))

    ncols = max(1, min(ncols, len(pathways)))
    nrows = math.ceil(len(pathways) / ncols)
    if figsize is None:
        figsize = (max(8.4, 3.55 * ncols), max(5.4, 4.25 * nrows))
    default_colors = [
        "#5B8DB8", "#80B99A", "#C489A6", "#D99B65",
        "#8E79B9", "#63A6A0", "#D3B75E", "#D88883",
    ]
    colors = PlotTheme.get_plot_colors(
        style, palette, default=default_colors, n=len(group_names)
    )
    rng = np.random.default_rng(random_seed)
    statistic_rows = []
    pairwise_label = "Wilcoxon" if pairwise_test == "mannwhitney" else "Welch t-test"
    global_label = "Kruskal-Wallis" if global_test == "kruskal" else "ANOVA"

    def format_p(value):
        if not np.isfinite(value):
            return "NA"
        return f"{value:.1e}" if value < 0.001 else f"{value:.3f}"

    def adjusted_pvalues(values):
        result = np.full(len(values), np.nan)
        finite = np.isfinite(values)
        if finite.any():
            if p_adjust.lower() == "none":
                result[finite] = np.asarray(values)[finite]
            else:
                methods = {"bh": "fdr_bh", "holm": "holm", "bonferroni": "bonferroni"}
                result[finite] = multipletests(
                    np.asarray(values)[finite], method=methods[p_adjust.lower()]
                )[1]
        return result

    def half_violin(ax, values, x, color, grid):
        values = values[np.isfinite(values)]
        if not len(values):
            return
        span = max(float(np.ptp(grid)), 0.1)
        if len(values) < 2 or np.std(values) == 0:
            bandwidth = max(span / 30, 0.05)
        else:
            sd = float(np.std(values, ddof=1))
            q75, q25 = np.percentile(values, [75, 25])
            scale = min(sd, (q75 - q25) / 1.34) if q75 > q25 else sd
            bandwidth = max(0.9 * max(scale, 0.1) * len(values) ** (-0.2), span / 250)
        density = np.exp(-0.5 * ((grid[:, None] - values) / bandwidth) ** 2).sum(axis=1)
        density = density / density.max() * 0.28
        visible = density >= density.max() * 0.02
        ax.fill_betweenx(
            grid, x - density, x, where=visible,
            color=color, alpha=0.18, linewidth=0, zorder=1,
        )
        ax.plot(np.where(visible, x - density, np.nan), grid,
                color="#696969", linewidth=0.85, zorder=2)

    with PlotTheme.context(style or "nature", palette):
        plt.rcParams["font.family"] = "sans-serif"
        plt.rcParams["font.sans-serif"] = ["Arial", "DejaVu Sans"]
        fig, axes = plt.subplots(nrows, ncols, figsize=figsize, dpi=dpi,
                                 squeeze=False, facecolor="white")
        for index, pathway in enumerate(pathways):
            ax = axes.flat[index]
            grouped = grouped_values(pathway)
            finite_values = [value[np.isfinite(value)] for value in grouped.values() if np.isfinite(value).any()]
            if not finite_values:
                ax.axis("off")
                continue
            all_values = np.concatenate(finite_values)
            ymin, ymax = float(all_values.min()), float(all_values.max())
            value_range = max(ymax - ymin, 1.0)
            pathway_global_p = np.nan if two_groups else global_p(list(grouped.values()))
            raw_pvalues = [pair_p(grouped[left], grouped[right]) for left, right in comparisons]
            adjusted = np.asarray(raw_pvalues) if two_groups else adjusted_pvalues(raw_pvalues)
            step, bracket_height = value_range * 0.13, value_range * 0.025
            start = ymax + value_range * 0.11
            upper = start + max(0, len(comparisons) - 1) * step + bracket_height + value_range * 0.15
            grid = np.linspace(ymin - value_range * 0.04, ymax + value_range * 0.04, 256)

            for group_index, group in enumerate(group_names):
                x = group_index + 1
                finite = grouped[group][np.isfinite(grouped[group])]
                half_violin(ax, finite, x - 0.05, colors[group_index], grid)
                ax.scatter(x + rng.uniform(-0.12, 0.05, len(finite)), finite, s=15,
                           facecolor=colors[group_index], edgecolor="white", linewidth=0.35,
                           alpha=0.72, zorder=3)
                ax.boxplot([finite], positions=[x + 0.10], widths=0.22, patch_artist=True,
                           showfliers=False, medianprops={"color": "white", "linewidth": 1.4},
                           boxprops={"facecolor": colors[group_index], "edgecolor": "#555555",
                                     "linewidth": 0.9, "alpha": 0.82},
                           whiskerprops={"color": "#555555", "linewidth": 0.9},
                           capprops={"color": "#555555", "linewidth": 0.9})

            for pair_index, ((left, right), raw_p, adjusted_p) in enumerate(
                zip(comparisons, raw_pvalues, adjusted)
            ):
                y = start + pair_index * step
                height = bracket_height
                x1, x2 = group_names.index(left) + 1, group_names.index(right) + 1
                ax.plot([x1, x1, x2, x2], [y, y + height, y + height, y],
                        color="#333333", linewidth=0.85, clip_on=False)
                prefix = "" if two_groups or p_adjust.lower() == "none" else f"{p_adjust} "
                ax.text((x1 + x2) / 2, y + height,
                        f"{pairwise_label}, {prefix}P={format_p(adjusted_p)}",
                        ha="center", va="bottom", fontsize=7.7)
                statistic_rows.append({
                    "Pathway": pathway, "test_type": "pairwise", "group1": left,
                    "group2": right, "raw_pvalue": raw_p,
                    "adjusted_pvalue": np.nan if two_groups else adjusted_p,
                    "p_adjust_method": "" if two_groups else p_adjust,
                })
            if not two_groups:
                statistic_rows.append({
                    "Pathway": pathway, "test_type": "global", "group1": "", "group2": "",
                    "raw_pvalue": pathway_global_p, "adjusted_pvalue": np.nan,
                    "p_adjust_method": "",
                })

            ax.set_xlim(0.45, len(group_names) + 0.55)
            ax.set_ylim(ymin - value_range * 0.10, upper)
            ax.set_xticks(np.arange(1, len(group_names) + 1), group_names, fontsize=9)
            ax.tick_params(axis="x", length=0)
            ax.tick_params(axis="y", labelsize=8.7)
            ax.grid(axis="y", color="#E7E7E7", linewidth=0.7)
            ax.set_axisbelow(True)
            for spine in ax.spines.values():
                spine.set_visible(True)
                spine.set_color("#404040")
                spine.set_linewidth(0.85)
            panel_label = _wrap_panel_label(pathway)
            label_lines = panel_label.count("\n") + 1
            strip_height = 0.075 if label_lines == 1 else 0.125
            ax.add_patch(Rectangle((0, 1), 1, strip_height, transform=ax.transAxes,
                                   facecolor="#E8E8E8", edgecolor="#707070",
                                   linewidth=0.75, clip_on=False))
            ax.text(0.5, 1 + strip_height / 2, panel_label,
                    transform=ax.transAxes, ha="center", va="center",
                    fontsize=9.4 if label_lines == 1 else 8.8,
                    fontweight="bold", linespacing=0.95)
            if not two_groups:
                ax.text(0.02, 0.985,
                        f"{global_label}, P={format_p(pathway_global_p)}",
                        transform=ax.transAxes, ha="left", va="top", fontsize=8.1)
            ax.set_xlabel("")
            ax.set_ylabel("Activity score" if index % ncols == 0 else "")

        for empty in range(len(pathways), nrows * ncols):
            axes.flat[empty].axis("off")
        fig.suptitle(title, fontsize=13, y=0.995)
        fig.subplots_adjust(left=0.08, right=0.985, top=0.89, bottom=0.10,
                            wspace=0.26, hspace=0.46)
        apply_figure_style(
            fig,
            style,
            axes=[axes.flat[index] for index in range(len(pathways))],
            grid_axis="y",
            border="full",
        )

        if output_file is not None:
            save_figure(fig, output_file, dpi=dpi, facecolor="white")
            logger.info("Group-comparison figure saved: %s", output_file)
        if statistics_file is None and output_file is not None:
            statistics_file = str(Path(output_file).with_suffix(".statistics.tsv"))
        if statistics_file is not None:
            Path(statistics_file).parent.mkdir(parents=True, exist_ok=True)
            pd.DataFrame(statistic_rows).to_csv(statistics_file, sep="\t", index=False)
            logger.info("Group-comparison statistics saved: %s", statistics_file)
        return fig


def plot_sample_correlation(
    scores_df: pd.DataFrame,
    method: str = "pearson",
    annotation_col: pd.DataFrame = None,
    title: str = "Sample Correlation",
    output_file: str = None,
    figsize: tuple = None,
    dpi: int = 300,
    style: str = "nature",
    palette: str = None,
) -> matplotlib.figure.Figure:
    """Plot clustered sample correlations with ellipses below and values above the diagonal."""
    if method not in ("pearson", "spearman"):
        raise ValueError(f"Unsupported correlation method: '{method}'. Expected 'pearson' or 'spearman'.")

    corr = scores_df.apply(pd.to_numeric, errors="coerce").corr(method=method)
    values = corr.to_numpy(dtype=float, copy=True)
    values[~np.isfinite(values)] = 0
    np.fill_diagonal(values, 1)
    corr = pd.DataFrame(values, index=corr.index, columns=corr.columns)
    n = len(corr)
    if n < 2:
        raise ValueError("Sample correlation requires at least two samples")

    distance = np.clip(1 - values, 0, None)
    link = linkage(squareform(distance, checks=False), method="average", optimal_ordering=True)
    order = leaves_list(link)
    corr = corr.iloc[order, order]
    samples = corr.columns.tolist()

    metadata = pd.DataFrame(index=samples)
    if annotation_col is not None and not annotation_col.empty:
        metadata = annotation_col.reindex(samples).dropna(axis=1, how="all").astype(str)

    category_count = sum(metadata[column].nunique() for column in metadata.columns)
    annotation_colors = PlotTheme.get_plot_colors(
        style, palette, default=_ANNOTATION_COLORS, n=max(1, category_count)
    )
    mappings = {}
    offset = 0
    for column in metadata.columns:
        categories = list(dict.fromkeys(metadata[column]))
        mapping = {
            category: annotation_colors[offset + index]
            for index, category in enumerate(categories)
        }
        mappings[column] = mapping
        offset += len(categories)

    def ellipse(correlation: float, x: float, y: float) -> np.ndarray:
        theta = np.linspace(0, 2 * np.pi, 180)
        points = np.vstack([np.cos(theta), np.sin(theta)])
        scales = np.diag([np.sqrt(max(1 + correlation, 0)), np.sqrt(max(1 - correlation, 0))])
        angle = np.pi / 4
        rotation = np.array([[np.cos(angle), -np.sin(angle)], [np.sin(angle), np.cos(angle)]])
        result = rotation @ scales @ points * 0.40
        result[0] += x
        result[1] += y
        return result.T

    correlation_colors = PlotTheme.get_plot_colors(
        style,
        palette,
        default=["#3B4CC0", "#F7F7F7", "#CA0020"],
        divergent=True,
    )
    cmap = LinearSegmentedColormap.from_list("sample_correlation", correlation_colors)
    norm = Normalize(-1, 1)
    n_annotations = metadata.shape[1]
    heat_units = max(4.9, n * 0.295)
    label_units = max(0.78, min(1.60, max(map(len, samples)) * 0.070))
    legend_units = 1.25
    if figsize is None:
        width = max(9.2, 1.9 + heat_units + label_units + legend_units)
        height = width * (0.72 + max(0.06, 0.11 * n_annotations) + heat_units) / (
            0.65 + heat_units + label_units + legend_units
        )
        figsize = (width, height)
    style_name = resolve_style(style)
    font_scale = 1.0

    with PlotTheme.context(style_name, palette):
        plt.rcParams["font.family"] = "sans-serif"
        plt.rcParams["font.sans-serif"] = ["Arial", "DejaVu Sans"]
        fig = plt.figure(figsize=figsize, dpi=dpi, facecolor="white")
        grid = fig.add_gridspec(
            3, 4,
            width_ratios=[0.65, heat_units, label_units, legend_units],
            height_ratios=[0.72, max(0.06, 0.11 * n_annotations), heat_units],
            left=0.02, right=0.99, top=0.90, bottom=0.11,
            wspace=0.01, hspace=0.01,
        )
        ax_col = fig.add_subplot(grid[0, 1])
        ax_annotation = fig.add_subplot(grid[1, 1])
        ax_row = fig.add_subplot(grid[2, 0])
        ax_heat = fig.add_subplot(grid[2, 1])
        ax_labels = fig.add_subplot(grid[2, 2], sharey=ax_heat)
        ax_legend = fig.add_subplot(grid[:, 3], label="sample_correlation_legend")

        dendrogram(link, ax=ax_col, no_labels=True, color_threshold=0,
                   above_threshold_color="black", link_color_func=lambda _: "black")
        dendrogram(link, ax=ax_row, orientation="left", no_labels=True, color_threshold=0,
                   above_threshold_color="black", link_color_func=lambda _: "black")
        for axis in (ax_col, ax_row):
            axis.axis("off")
            for collection in axis.collections:
                collection.set_linewidth(0.8)

        if n_annotations:
            ax_annotation.set_xlim(-0.5, n - 0.5)
            ax_annotation.set_ylim(n_annotations - 0.5, -0.5)
            for row, column in enumerate(metadata.columns):
                categories = list(mappings[column])
                codes = np.array([categories.index(value) for value in metadata[column]])[None, :]
                ax_annotation.imshow(
                    codes, aspect="auto", interpolation="nearest",
                    cmap=ListedColormap([mappings[column][value] for value in categories]),
                    vmin=-0.5, vmax=max(0.5, len(categories) - 0.5),
                    extent=[-0.5, n - 0.5, row + 0.5, row - 0.5],
                )
        ax_annotation.axis("off")

        ax_heat.set_xlim(-0.5, n - 0.5)
        ax_heat.set_ylim(n - 0.5, -0.5)
        ax_heat.set_aspect("equal")
        ax_heat.set_xticks(np.arange(-0.5, n, 1), minor=True)
        ax_heat.set_yticks(np.arange(-0.5, n, 1), minor=True)
        ax_heat.grid(which="minor", color="#E6E6E6", linewidth=0.5)
        ordered = corr.to_numpy()
        for row in range(n):
            for column in range(n):
                if row == column:
                    ax_heat.text(column, row, "1", ha="center", va="center",
                                 color="#404040", fontsize=8.0 * font_scale,
                                 fontweight="bold", zorder=4)
                elif row > column:
                    ax_heat.add_patch(Polygon(
                        ellipse(ordered[row, column], column, row), closed=True,
                        facecolor=cmap(norm(ordered[row, column])), edgecolor="none", zorder=3,
                    ))
                else:
                    value = ordered[row, column]
                    ax_heat.text(
                        column, row, f"{value:.2f}", ha="center", va="center",
                        color="#404040",
                        fontsize=7.5 * font_scale, fontweight="bold", zorder=4,
                    )
        ax_heat.set_xticks(np.arange(n), samples, rotation=90, fontsize=8.2 * font_scale)
        ax_heat.set_yticks([])
        ax_heat.tick_params(which="both", length=0, pad=2)
        for spine in ax_heat.spines.values():
            spine.set_visible(True)
            spine.set_color("#404040")
            spine.set_linewidth(0.8)

        ax_labels.set_xlim(0, 1)
        ax_labels.set_ylim(n - 0.5, -0.5)
        ax_labels.axis("off")
        for row, sample in enumerate(samples):
            ax_labels.text(
                0, row, sample, ha="left", va="center", fontsize=8.4 * font_scale
            )

        ax_legend.axis("off")
        legend_rows = sum(1 + len(mappings[column]) for column in metadata.columns)
        figure_height = fig.get_size_inches()[1]
        legend_step = min(0.045, 0.24 / figure_height, 0.54 / max(1, legend_rows))
        heading_gap = min(0.05, 0.28 / figure_height)
        section_gap = min(0.06, 0.34 / figure_height)
        legend_y = 0.94
        for column in metadata.columns:
            ax_legend.text(
                0, legend_y, column, ha="left", va="top",
                fontsize=9.0 * font_scale, fontweight="bold",
            )
            categories = list(mappings[column].items())
            first_y = legend_y - heading_gap
            for index, (category, color) in enumerate(categories):
                y_pos = first_y - legend_step * index
                ax_legend.add_patch(Rectangle((0.02, y_pos - 0.009), 0.075, 0.018,
                                              transform=ax_legend.transAxes, facecolor=color, edgecolor="none"))
                ax_legend.text(
                    0.12, y_pos, category, ha="left", va="center",
                    fontsize=8.0 * font_scale,
                )
            legend_y = first_y - legend_step * max(0, len(categories) - 1) - section_gap
        ax_legend.text(0, legend_y, "Correlation", ha="left", va="top",
                       fontsize=9.0 * font_scale, fontweight="bold")
        color_top = legend_y - 0.28 / figure_height
        color_height = min(0.18, 1.35 / figure_height)
        color_axis = ax_legend.inset_axes(
            [0.02, max(0.04, color_top - color_height), 0.10, color_height],
            label="sample_correlation_colorbar",
        )
        colorbar = fig.colorbar(ScalarMappable(norm=norm, cmap=cmap), cax=color_axis)
        colorbar.set_ticks([-1, 0, 1])
        colorbar.ax.tick_params(labelsize=7.4 * font_scale, length=2, pad=1)
        colorbar.outline.set_linewidth(0.6)

        fig.canvas.draw()
        heat_position = ax_heat.get_position()
        column_position = ax_col.get_position()
        ax_col.set_position([heat_position.x0, column_position.y0, heat_position.width, column_position.height])
        annotation_position = ax_annotation.get_position()
        ax_annotation.set_position([heat_position.x0, annotation_position.y0,
                                    heat_position.width, annotation_position.height])
        row_position = ax_row.get_position()
        ax_row.set_position([row_position.x0, heat_position.y0,
                             max(0.02, heat_position.x0 - 0.003 - row_position.x0), heat_position.height])
        label_position = ax_labels.get_position()
        ax_labels.set_position([heat_position.x1 + 0.003, heat_position.y0,
                                max(0.03, label_position.x1 - heat_position.x1 - 0.003), heat_position.height])

        title_x = (heat_position.x0 + heat_position.x1) / 2
        title_y = min(0.97, ax_col.get_position().y1 + 0.028)
        fig.suptitle(
            title, fontsize=12.5 * font_scale, fontweight="bold", x=title_x, y=title_y
        )
        apply_figure_style(fig, style_name, axes=[ax_heat], border="full")

        if output_file is not None:
            save_figure(fig, output_file, dpi=dpi, facecolor="white")
            logger.info("Sample-correlation figure saved: %s", output_file)
        return fig
