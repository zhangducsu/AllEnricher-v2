"""Create publication-oriented single- and multi-pathway GSEA figures."""

import logging
import re
import textwrap
from typing import Dict, List, Optional, Set

import matplotlib
matplotlib.use("Agg")
from matplotlib import colors as mcolors
from matplotlib.patches import Rectangle
from matplotlib.ticker import MaxNLocator
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np
import pandas as pd
import seaborn as sns

from .plot_theme import PlotTheme, apply_figure_style, save_figure

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Shared calculations and file output.
# ---------------------------------------------------------------------------

def _compute_running_es(
    ranked_genes: List[str],
    gene_set: Set[str],
    gene_weights: Optional[Dict[str, float]] = None,
) -> np.ndarray:
    """Calculate the running enrichment-score trajectory along a ranked gene list."""
    n = len(ranked_genes)
    nh = len(gene_set & set(ranked_genes))

    if nh == 0:
        return np.zeros(n)

    # Weighted hits increase the score; misses contribute a uniform decrement.
    nr = sum(
        abs(gene_weights.get(g, 1.0)) for g in gene_set if g in gene_weights
    ) if gene_weights else nh
    hit_inc = 1.0 / nr if nr > 0 else 0
    miss_inc = 1.0 / (n - nh) if (n - nh) > 0 else 0

    running_es = np.zeros(n)
    running_sum = 0.0

    for i, gene in enumerate(ranked_genes):
        if gene in gene_set:
            weight = abs(gene_weights.get(gene, 1.0)) if gene_weights else 1.0
            running_sum += hit_inc * weight
        else:
            running_sum -= miss_inc
        running_es[i] = running_sum

    return running_es


def _get_gene_set_positions(
    ranked_genes: List[str],
    gene_set: Set[str],
) -> List[int]:
    """Return zero-based ranks occupied by members of one gene set."""
    return [i for i, g in enumerate(ranked_genes) if g in gene_set]


def _save_figure(
    fig: plt.Figure,
    output_file: Optional[str],
    dpi: int = 300,
    pad_inches: Optional[float] = None,
):
    """Save the figure using the format implied by the output filename."""
    if output_file:
        kwargs = {"pad_inches": pad_inches} if pad_inches is not None else {}
        save_figure(fig, output_file, dpi=dpi, **kwargs)
        plt.close(fig)
        logger.info("GSEA figure saved: %s", output_file)


# ---------------------------------------------------------------------------
# Standard three-panel GSEA running-score figure.
# ---------------------------------------------------------------------------

def _first_existing_column(df: pd.DataFrame, candidates: List[str]) -> Optional[str]:
    lookup = {str(col).lower(): col for col in df.columns}
    for candidate in candidates:
        if candidate in df.columns:
            return candidate
        matched = lookup.get(candidate.lower())
        if matched is not None:
            return matched
    return None


def _parse_ratio_value(value) -> float:
    if pd.isna(value):
        return np.nan
    if isinstance(value, (int, float, np.number)):
        return float(value)

    text = str(value).strip()
    if "/" in text:
        left, right = text.split("/", 1)
        try:
            denominator = float(right)
            return float(left) / denominator if denominator else np.nan
        except ValueError:
            return np.nan

    try:
        return float(text)
    except ValueError:
        return np.nan


def _clean_lollipop_label(value: str, width: int = 34) -> str:
    label = str(value).strip()
    if "|" in label:
        label = label.rsplit("|", 1)[-1].strip()
    label = " ".join(label.replace("_", " ").split())
    return "\n".join(textwrap.wrap(label, width=width, break_long_words=False))


def _size_area(values: pd.Series, max_diameter_mm: float = 8.8) -> np.ndarray:
    arr = np.asarray(values, dtype=float)
    if arr.size == 0:
        return np.array([])

    max_value = np.nanmax(arr)
    if not np.isfinite(max_value) or max_value <= 0:
        return np.full(arr.shape, 120.0)

    max_diameter_pt = max_diameter_mm / 25.4 * 72.0
    max_area = max_diameter_pt ** 2
    return np.maximum(arr / max_value * max_area, 10.0)


def _nice_count_breaks(values: pd.Series, n_breaks: int = 3) -> list[int]:
    arr = np.asarray(values, dtype=float)
    arr = arr[np.isfinite(arr)]
    if arr.size == 0:
        return [1]
    vmin = int(np.floor(arr.min()))
    vmax = int(np.ceil(arr.max()))
    if vmin == vmax:
        return [vmax]
    ticks = sorted(set(int(round(x)) for x in np.linspace(vmin, vmax, n_breaks)))
    if vmax not in ticks:
        ticks[-1] = vmax
    return ticks


def _fdr_tick_labels(score_breaks: np.ndarray) -> list[str]:
    return [
        f"{10 ** (-float(score)):.1e}".replace("e-0", "e-").replace("e+0", "e+")
        for score in score_breaks
    ]


def _compute_lollipop_x_layout(values: pd.Series, signed: bool) -> tuple[float, float, float, float]:
    arr = np.asarray(values, dtype=float)
    arr = arr[np.isfinite(arr)]
    if arr.size == 0:
        return 0.0, 1.0, 0.0, 1.0
    xmin = float(arr.min())
    xmax = float(arr.max())
    data_range = xmax - xmin
    if data_range <= 0:
        data_range = max(abs(xmax), 1.0) * 0.20
    if signed and xmin < 0:
        limit = max(abs(xmin), abs(xmax), 1.0) * 1.12
        return 0.0, limit, -limit, limit
    xlim_high = max(xmax + 0.14 * data_range, xmax * 1.08, 1.0)
    return 0.0, xlim_high * 0.98, 0.0, xlim_high


def _prepare_lollipop_df(results_df: pd.DataFrame, top_n: int) -> tuple[pd.DataFrame, str]:
    df = results_df.copy()

    term_col = _first_existing_column(
        df,
        [
            "Term_Name", "term_name", "Description", "description", "name", "Name",
            "term", "Term", "pathway", "Pathway", "Term_ID", "ID",
        ],
    )
    if term_col is None:
        raise ValueError("lollipop plot requires a pathway/term column")

    enrich_col = _first_existing_column(
        df,
        [
            "EnrichFactor", "enrichfactor", "enrich_factor", "EnrichmentFactor",
            "FoldEnrichment", "fold_enrichment", "RichFactor", "rich_factor", "Rich_Factor",
        ],
    )

    metric = None
    metric_label = "EnrichFactor"
    if enrich_col is not None:
        metric = pd.to_numeric(df[enrich_col], errors="coerce")
    else:
        gene_ratio_col = _first_existing_column(df, ["GeneRatio", "gene_ratio", "geneRatio"])
        bg_ratio_col = _first_existing_column(df, ["BgRatio", "bg_ratio", "BackgroundRatio", "background_ratio"])
        if gene_ratio_col and bg_ratio_col:
            gene_ratio = df[gene_ratio_col].map(_parse_ratio_value)
            bg_ratio = df[bg_ratio_col].map(_parse_ratio_value)
            metric = gene_ratio / bg_ratio.replace(0, np.nan)

    if metric is None or not np.isfinite(metric).any():
        nes_col = _first_existing_column(df, ["NES", "nes", "enrichmentScore"])
        if nes_col is None:
            raise ValueError("lollipop plot requires EnrichFactor or NES columns")
        metric = pd.to_numeric(df[nes_col], errors="coerce")
        metric_label = "NES"

    fdr_col = _first_existing_column(
        df,
        [
            "Adjusted_P_Value", "FDR", "FDR q-val", "p.adjust", "padj",
            "qvalue", "q_value", "adj_p", "pvalue", "P_Value", "NOM p-val",
        ],
    )
    if fdr_col is None:
        raise ValueError("lollipop plot requires an FDR/p-value column")

    gene_count_col = _first_existing_column(
        df,
        ["size", "Gene_Count", "gene_count", "GeneCount", "setSize", "Count", "count"],
    )
    if gene_count_col is None:
        raise ValueError("lollipop plot requires a gene count column")

    out = pd.DataFrame({
        "term": df[term_col].astype(str),
        "metric": pd.to_numeric(metric, errors="coerce"),
        "fdr": pd.to_numeric(df[fdr_col], errors="coerce"),
        "gene_count": pd.to_numeric(df[gene_count_col], errors="coerce"),
    })
    out = out.dropna(subset=["metric", "fdr", "gene_count"]).copy()
    out = out[(out["fdr"] > 0) & np.isfinite(out["metric"]) & np.isfinite(out["gene_count"])]
    if out.empty:
        raise ValueError("no valid rows available for lollipop plot")

    out["_metric_abs"] = out["metric"].abs()
    out["term_wrap"] = out["term"].map(_clean_lollipop_label)

    if metric_label == "NES" and (out["metric"] > 0).any() and (out["metric"] < 0).any():
        sig_df = out[out["fdr"] < 0.05].copy()
        if sig_df.empty:
            sig_df = out.sort_values("fdr", ascending=True).head(top_n).copy() if top_n else out.copy()
        top_each = max(1, top_n // 2) if top_n else len(sig_df)
        up = sig_df[sig_df["metric"] > 0].sort_values("metric", ascending=False).head(top_each)
        down = sig_df[sig_df["metric"] < 0].sort_values("metric", ascending=True).head(top_each)
        down = down.sort_values("metric", ascending=False)
        out = pd.concat([up, down], ignore_index=True)
    else:
        if metric_label == "EnrichFactor" and top_n:
            out = out.sort_values(["fdr", "_metric_abs"], ascending=[True, False]).head(top_n)
        out = out.sort_values("metric", ascending=False)
        if top_n and metric_label != "EnrichFactor":
            out = out.head(top_n)

    return out.reset_index(drop=True), metric_label


def _count_core_genes(value) -> float:
    if pd.isna(value):
        return np.nan
    if isinstance(value, (int, float, np.number)):
        return float(value)
    text = str(value).strip()
    if not text or text.upper() == "NA":
        return np.nan
    genes = [g for g in re.split(r"[/,;|\s]+", text) if g]
    return float(len(set(genes)))


def _format_sig_short(value, digits: int = 1) -> str:
    if pd.isna(value) or not np.isfinite(value):
        return "NA"
    text = f"{float(value):.{digits}e}"
    text = re.sub(r"e([+-])0+", r"e\1", text)
    return text.replace("e+", "e").replace(".0e", "e")


def plot_gsea_barplot(
    results_df: pd.DataFrame,
    database: str = "",
    top_n: int = 20,
    title: Optional[str] = None,
    output_file: Optional[str] = None,
    figsize: Optional[tuple] = None,
    style: Optional[str] = None,
    palette: Optional[str] = None,
    dpi: int = 300,
) -> matplotlib.figure.Figure:
    """Plot positive and negative GSEA terms ordered by normalized enrichment score."""
    df, metric_label = _prepare_lollipop_df(results_df, top_n)
    if metric_label != "NES":
        raise ValueError("GSEA barplot requires an NES column")

    df = df.sort_values("metric", ascending=True).reset_index(drop=True)
    df["fdr_score"] = np.minimum(
        -np.log10(np.maximum(df["fdr"].astype(float), np.finfo(float).tiny)),
        10.0,
    )
    score_min = float(df["fdr_score"].min())
    score_max = float(df["fdr_score"].max())
    if np.isclose(score_min, score_max):
        score_min -= 0.5
        score_max += 0.5

    colors = PlotTheme.get_plot_colors(
        style,
        palette,
        sequential=True,
        n=6,
    )
    cmap = mcolors.LinearSegmentedColormap.from_list("gsea_barplot_fdr", colors)
    norm = mcolors.Normalize(vmin=score_min, vmax=score_max)
    max_abs = max(float(df["metric"].abs().max()), 1.0)
    offset = max_abs * 0.03
    y = np.arange(len(df))
    value_labels = [
        f"{int(round(count))} / {_format_sig_short(fdr)}"
        for count, fdr in zip(df["gene_count"], df["fdr"])
    ]

    if figsize is None:
        figsize = (7.45, max(4.25, len(df) * 0.205 + 1.35))

    with PlotTheme.context(style or "nature", palette):
        fig = plt.figure(figsize=figsize, dpi=dpi)
        grid = fig.add_gridspec(1, 2, width_ratios=[1.0, 0.055], wspace=0.07)
        ax = fig.add_subplot(grid[0, 0])
        cbar_ax = fig.add_subplot(grid[0, 1])

        ax.barh(
            y,
            df["metric"],
            height=0.62,
            color=[cmap(norm(value)) for value in df["fdr_score"]],
            edgecolor="#4A4A4A",
            linewidth=0.35,
        )
        ax.axvline(0, color="#333333", linewidth=0.65)

        for index, row in df.iterrows():
            positive = row["metric"] >= 0
            ax.text(
                -offset if positive else offset,
                index,
                row["term_wrap"],
                ha="right" if positive else "left",
                va="center",
                fontsize=8.7,
                linespacing=0.92,
            )
            ax.text(
                row["metric"] + offset if positive else row["metric"] - offset,
                index,
                value_labels[index],
                ha="left" if positive else "right",
                va="center",
                fontsize=7.8,
                color="#333333",
            )

        ax.text(0.76, 1.006, "UP", transform=ax.transAxes, ha="center", va="bottom",
                fontsize=9.5, fontweight="bold")
        ax.text(0.24, 1.006, "DOWN", transform=ax.transAxes, ha="center", va="bottom",
                fontsize=9.5, fontweight="bold")
        ax.set_xlim(-max_abs * 1.48, max_abs * 1.56)
        ax.set_ylim(-0.65, len(df) - 0.15)
        ax.set_yticks([])
        ax.set_xlabel("Normalized enrichment score (NES)", fontsize=10)
        chart_title = title or f"{database + ' ' if database else ''}GSEA NES Ranking"
        fig.suptitle(f"{chart_title} (core# / FDR)", fontsize=11.4, y=0.985)
        ax.xaxis.set_major_locator(MaxNLocator(nbins=7))
        ax.grid(axis="x", color="#D9D9D9", linewidth=0.35, linestyle="--")
        ax.grid(axis="y", visible=False)
        for spine in ("top", "right", "left"):
            ax.spines[spine].set_visible(False)

        scalar = plt.cm.ScalarMappable(norm=norm, cmap=cmap)
        scalar.set_array([])
        colorbar = fig.colorbar(scalar, cax=cbar_ax)
        colorbar.ax.set_title(r"$-log_{10}$(FDR)", fontsize=8.2, pad=6)
        colorbar.ax.tick_params(labelsize=7.5, length=2)
        colorbar.outline.set_linewidth(0.45)

        apply_figure_style(fig, style, axes=[ax], grid_axis="x", border=None)
        fig.subplots_adjust(left=0.045, right=0.94, top=0.88, bottom=0.13)
        _save_figure(fig, output_file, dpi=dpi, pad_inches=0.06)
    return fig


def plot_gsea_enrichment(
    ranked_genes: List[str],
    gene_weights: Dict[str, float],
    gene_set: Set[str],
    es: float,
    nes: float,
    pvalue: float,
    padj: Optional[float] = None,
    title: str = None,
    output_file: str = None,
    figsize: tuple = (6.60, 4.40),
    style: Optional[str] = None,
    palette: Optional[str] = None,
    dpi: int = 300,
) -> matplotlib.figure.Figure:
    """Plot one pathway as the standard three-panel GSEA running-score figure."""
    running_es = _compute_running_es(ranked_genes, gene_set, gene_weights)
    positions = _get_gene_set_positions(ranked_genes, gene_set)
    weights = np.array([gene_weights.get(g, 0.0) for g in ranked_genes], dtype=float)
    n = len(ranked_genes)
    x = np.arange(1, n + 1)

    direction_colors = PlotTheme.get_plot_colors(
        style, palette, default=["#FF9999", "#99CC00"], divergent=True
    )
    heat_colors = PlotTheme.get_plot_colors(
        style, palette, default=["#FF9999", "#FFFFFF", "#99CC00"], divergent=True
    )
    rank_colors = heat_colors

    def _fmt_p(value):
        if value is None or pd.isna(value) or not np.isfinite(float(value)):
            return "NA"
        value = float(value)
        if value < 0.001:
            return "< 0.001"
        if value < 0.01:
            return f"{value:.2e}"
        return f"{value:.3f}"

    display_title = _clean_lollipop_label(title, width=58) if title else f"NES = {float(nes):.2f}"

    with PlotTheme.context(style or "nature", palette):
        plt.rcParams["font.family"] = "sans-serif"
        plt.rcParams["font.sans-serif"] = ["Arial", "DejaVu Sans"]
        plt.rcParams["svg.fonttype"] = "none"
        plt.rcParams["pdf.fonttype"] = 42

        fig = plt.figure(figsize=figsize, dpi=dpi)
        fig.suptitle(display_title, y=0.995, fontsize=1, alpha=0)
        gs = gridspec.GridSpec(3, 1, height_ratios=[2.85, 0.60, 1.22], hspace=0.012, figure=fig)
        ax1 = fig.add_subplot(gs[0])

        curve_color = direction_colors[-1] if float(nes) >= 0 else direction_colors[0]
        ax1.plot(x, running_es, color=curve_color, linewidth=2.35)
        ax1.axhline(0, linestyle="--", linewidth=0.9, color="black")
        y_pad = max(0.05, 0.06 * (float(np.nanmax(running_es)) - float(np.nanmin(running_es))))
        ax1.set_xlim(1, n)
        ax1.set_ylim(float(np.nanmin(running_es)) - y_pad, float(np.nanmax(running_es)) + y_pad)
        ax1.set_ylabel("Running Enrichment Score", fontsize=13.5)
        ax1.set_title(display_title, fontsize=15, fontweight="bold", pad=3)
        ax1.tick_params(axis="x", which="both", bottom=False, labelbottom=False)
        ax1.tick_params(axis="y", labelsize=11.5)
        ax1.text(
            0.985,
            0.955,
            f"NES: {float(nes):.2f}\nP value: {_fmt_p(pvalue)}\nAdjusted P value: {_fmt_p(padj)}",
            transform=ax1.transAxes,
            ha="right",
            va="top",
            fontsize=10.8,
            color="0.20",
            fontstyle="italic",
            linespacing=1.05,
        )

        ax2 = fig.add_subplot(gs[1], sharex=ax1)
        heat_cmap = mcolors.LinearSegmentedColormap.from_list("gsea_heat", heat_colors)
        heat_values = np.linspace(-1, 1, n, dtype=float)[None, :]
        ax2.imshow(
            heat_values,
            aspect="auto",
            cmap=heat_cmap,
            extent=[1, n, 0.00, 0.28],
            origin="lower",
            interpolation="nearest",
            alpha=0.90,
        )
        hit_positions = np.array(positions, dtype=int) + 1
        ax2.vlines(hit_positions, 0.28, 1.0, color="black", linewidth=0.62)
        ax2.set_ylim(0, 1)
        ax2.set_yticks([])
        ax2.tick_params(axis="x", which="both", bottom=False, labelbottom=False)

        ax3 = fig.add_subplot(gs[2], sharex=ax1)
        score_limit = float(np.nanmax(np.abs(weights)))
        if not np.isfinite(score_limit) or score_limit <= 0:
            score_limit = 1.0
        rank_cmap = mcolors.LinearSegmentedColormap.from_list("gsea_rank", rank_colors)
        rank_norm = mcolors.Normalize(vmin=-score_limit, vmax=score_limit)
        step = max(1, int(np.ceil(n / 3000)))
        ax3.bar(
            x[::step],
            weights[::step],
            width=step + 0.35,
            color=rank_cmap(rank_norm(weights[::step])),
            edgecolor="none",
            align="center",
        )
        ax3.axhline(0, linestyle="--", linewidth=0.9, color="black")
        ax3.set_ylabel("Ranked List", fontsize=13.5)
        ax3.set_xlabel("Rank in Ordered Dataset", fontsize=13.5)
        ax3.tick_params(axis="both", labelsize=11.5)

        for axis in (ax1, ax2, ax3):
            for spine in axis.spines.values():
                spine.set_linewidth(0.8)
                spine.set_color("black")
            axis.set_facecolor("white")
            axis.grid(False)

        apply_figure_style(
            fig, style, axes=[ax1, ax3], grid_axis="y", border="style"
        )
        fig.subplots_adjust(left=0.140, right=0.992, top=0.935, bottom=0.115)
        _save_figure(fig, output_file, dpi=dpi, pad_inches=0.02)
    return fig


def plot_gsea_multi_enrichment(
    results_df: pd.DataFrame,
    selected_ids: List[str],
    ranked_genes: List[str],
    gene_weights: Dict[str, float],
    gene_sets: Dict[str, Set[str]],
    output_file: str = None,
    figsize: tuple = None,
    style: Optional[str] = None,
    palette: Optional[str] = None,
    dpi: int = 300,
) -> matplotlib.figure.Figure:
    """Plot multiple running-score trajectories with one hit strip per pathway."""
    if not selected_ids:
        raise ValueError("at least one pathway is required")
    if len(selected_ids) > 8:
        raise ValueError("multi-pathway GSEA supports at most 8 pathways")

    id_col = _first_existing_column(results_df, ["pathway", "Term_ID", "term_id", "ID", "id"])
    name_col = _first_existing_column(
        results_df, ["Term_Name", "Description", "pathway", "term_name"]
    )
    nes_col = _first_existing_column(results_df, ["NES", "nes"])
    p_col = _first_existing_column(results_df, ["pval", "p_value", "pvalue", "P_Value", "NOM p-val"])
    padj_col = _first_existing_column(
        results_df, ["Adjusted_P_Value", "FDR", "p.adjust", "padj", "qvalue"]
    )
    if id_col is None:
        raise ValueError("multi-pathway GSEA requires a pathway ID column")

    def _fmt(value) -> str:
        if value is None or pd.isna(value) or not np.isfinite(float(value)):
            return "NA"
        value = float(value)
        if value < 0.001:
            return "<0.001"
        if value < 0.01:
            return f"{value:.2e}"
        return f"{value:.3f}"

    default_colors = [
        "#99CC00", "#FF9999", "#9B6BA6", "#4DBBD5",
        "#E6AB02", "#7E6148", "#00A087", "#F39B7F",
    ]
    curves = []
    hit_arrays = []
    descriptions = []
    legend_labels = []
    nes_values = []

    for set_id in selected_ids:
        if set_id not in gene_sets:
            raise ValueError(f"gene set not found: {set_id}")
        rows = results_df[results_df[id_col].astype(str) == str(set_id)]
        if rows.empty:
            raise ValueError(f"GSEA result not found: {set_id}")
        row = rows.iloc[0]
        description = str(row[name_col]) if name_col else str(set_id)
        description = description.rsplit("|", 1)[-1].strip().replace("_", " ")
        curve = _compute_running_es(ranked_genes, gene_sets[set_id], gene_weights)
        hits = np.array([gene in gene_sets[set_id] for gene in ranked_genes], dtype=bool)
        nes = float(row[nes_col]) if nes_col and pd.notna(row[nes_col]) else float(curve[np.argmax(np.abs(curve))])
        pvalue = row[p_col] if p_col else None
        padj = row[padj_col] if padj_col else None

        curves.append(curve)
        hit_arrays.append(hits)
        descriptions.append(description)
        nes_values.append(nes)
        legend_labels.append(f"{description} | NES {nes:.2f} | P {_fmt(pvalue)} | FDR {_fmt(padj)}")

    nsets = len(selected_ids)
    n = len(ranked_genes)
    x = np.arange(1, n + 1)
    colors = PlotTheme.get_plot_colors(
        style, palette, default=default_colors, n=nsets
    )
    heat_colors = PlotTheme.get_plot_colors(
        style, palette, default=["#FF9999", "#FFFFFF", "#99CC00"], divergent=True
    )
    if figsize is None:
        figsize = (7.2, 2.55 + 0.46 * nsets)

    with PlotTheme.context(style or "nature", palette):
        plt.rcParams["font.family"] = "sans-serif"
        plt.rcParams["font.sans-serif"] = ["Arial", "DejaVu Sans"]
        plt.rcParams["svg.fonttype"] = "none"
        plt.rcParams["pdf.fonttype"] = 42

        fig = plt.figure(figsize=figsize, dpi=dpi)
        gs = gridspec.GridSpec(
            nsets + 1,
            1,
            height_ratios=[2.7] + [0.43] * nsets,
            hspace=0.012,
            figure=fig,
        )
        ax_top = fig.add_subplot(gs[0])
        ymin = min(float(np.min(curve)) for curve in curves)
        ymax = max(float(np.max(curve)) for curve in curves)
        data_span = max(ymax - ymin, 0.10)
        ypad = max(0.04, 0.06 * data_span)
        legend_pad = max(ypad, data_span * (0.08 + 0.055 * nsets))

        for curve, color, label in zip(curves, colors, legend_labels):
            ax_top.plot(x, curve, color=color, linewidth=2.15, label=label)
        ax_top.axhline(0, linestyle="--", color="black", linewidth=0.9)
        ax_top.set_xlim(1, n)
        ax_top.set_ylim(ymin - ypad, ymax + legend_pad)
        ax_top.set_ylabel("Running ES", fontsize=14)
        ax_top.tick_params(axis="x", bottom=False, labelbottom=False)
        ax_top.tick_params(axis="y", labelsize=11.5)
        ax_top.legend(
            loc="upper right" if np.nanmean(nes_values) >= 0 else "upper left",
            frameon=True,
            facecolor="white",
            framealpha=0.92,
            edgecolor="none",
            fontsize=9.2,
            handlelength=1.8,
            labelspacing=0.45,
            borderaxespad=0.25,
        )

        heat_cmap = mcolors.LinearSegmentedColormap.from_list(
            "gsea_multi_heat", heat_colors
        )
        heat_values = np.linspace(-1, 1, n)[None, :]
        for index, (hits, color, description) in enumerate(
            zip(hit_arrays, colors, descriptions), start=1
        ):
            ax = fig.add_subplot(gs[index], sharex=ax_top)
            ax.imshow(
                heat_values,
                aspect="auto",
                cmap=heat_cmap,
                extent=[1, n, 0.0, 0.30],
                origin="lower",
                interpolation="nearest",
                alpha=0.90,
            )
            ax.vlines(np.flatnonzero(hits) + 1, 0.30, 1.0, color=color, linewidth=0.65)
            ax.set_ylim(0, 1)
            ax.set_yticks([])
            short_label = "\n".join(textwrap.wrap(description, width=24, break_long_words=False))
            ax.text(
                -0.012, 0.53, short_label,
                transform=ax.transAxes,
                ha="right", va="center",
                fontsize=9.5, color=color, fontweight="bold", clip_on=False,
            )
            if index < nsets:
                ax.tick_params(axis="x", bottom=False, labelbottom=False)
            else:
                ax.set_xlabel("Rank in Ordered Dataset", fontsize=14)
                ax.tick_params(axis="x", labelsize=11.5)

        for axis in fig.axes:
            for spine in axis.spines.values():
                spine.set_linewidth(0.8)
                spine.set_color("black")
            axis.set_facecolor("white")
            axis.grid(False)

        apply_figure_style(
            fig,
            style,
            axes=[ax_top],
            grid_axis="y",
            border="style",
            max_text_scale=1.15,
        )
        fig.subplots_adjust(left=0.195, right=0.985, top=0.985, bottom=0.115)
        if output_file:
            fig.savefig(
                output_file,
                dpi=dpi,
                bbox_inches="tight",
                pad_inches=0.01,
                facecolor="white",
            )
            plt.close(fig)
    return fig


def plot_gsea_ridgeplot(
    results_df: pd.DataFrame,
    ranked_genes: List[str],
    gene_weights: Dict[str, float],
    gene_sets: Dict[str, Set[str]],
    top_n: int = 10,
    min_genes: int = 5,
    density_adjust: float = 1.0,
    output_file: str = None,
    figsize: tuple = None,
    style: Optional[str] = None,
    palette: Optional[str] = None,
    dpi: int = 300,
) -> matplotlib.figure.Figure:
    """Plot the distribution of ranking statistics for selected GSEA pathways."""
    id_col = _first_existing_column(results_df, ["pathway", "Term_ID", "term_id", "ID", "id"])
    name_col = _first_existing_column(
        results_df, ["Term_Name", "Description", "pathway", "term_name"]
    )
    sig_col = _first_existing_column(
        results_df,
        ["padj", "pval", "p_value", "pvalue", "P_Value", "NOM p-val", "Adjusted_P_Value", "FDR", "p.adjust", "qvalue"],
    )
    nes_col = _first_existing_column(results_df, ["NES", "nes"])
    if id_col is None or name_col is None or sig_col is None:
        raise ValueError("ridgeplot requires pathway ID, name, and significance columns")

    rows = []
    for _, row in results_df.iterrows():
        term_id = str(row[id_col])
        significance = pd.to_numeric(pd.Series([row[sig_col]]), errors="coerce").iloc[0]
        if term_id not in gene_sets or not np.isfinite(significance) or significance <= 0:
            continue
        values = np.asarray(
            [gene_weights[g] for g in gene_sets[term_id] if g in gene_weights],
            dtype=float,
        )
        values = values[np.isfinite(values)]
        if values.size < min_genes:
            continue
        nes = pd.to_numeric(pd.Series([row[nes_col]]), errors="coerce").iloc[0] if nes_col else np.nan
        rows.append({
            "term_id": term_id,
            "description": str(row[name_col]).rsplit("|", 1)[-1].strip().replace("_", " "),
            "significance": float(significance),
            "nes": float(nes) if np.isfinite(nes) else np.nan,
            "values": values,
        })

    if not rows:
        raise ValueError("no pathways contain enough ranked genes for ridgeplot")
    rows.sort(key=lambda row: (row["significance"], -abs(row["nes"]) if np.isfinite(row["nes"]) else 0))
    rows = rows[:top_n]

    all_scores = np.asarray([gene_weights[g] for g in ranked_genes if g in gene_weights], dtype=float)
    all_scores = all_scores[np.isfinite(all_scores)]
    if all_scores.size < 2:
        raise ValueError("ranked gene weights contain fewer than two finite values")
    score_min, score_max = float(all_scores.min()), float(all_scores.max())
    padding = max((score_max - score_min) * 0.025, 0.05)
    x_limits = (score_min - padding, score_max + padding)
    x_grid = np.linspace(*x_limits, 512)

    def density(values: np.ndarray) -> np.ndarray:
        span = x_limits[1] - x_limits[0]
        if np.unique(values).size < 2:
            bandwidth = max(span / 45.0, 0.05)
        else:
            sd = float(np.std(values, ddof=1))
            q75, q25 = np.percentile(values, [75, 25])
            scale = min(sd, (q75 - q25) / 1.34) if q75 > q25 else sd
            bandwidth = max(0.9 * max(scale, 0.1) * values.size ** (-0.2) * density_adjust, span / 300)
        z = (x_grid[:, None] - values[None, :]) / bandwidth
        return np.exp(-0.5 * z * z).sum(axis=1) / (values.size * bandwidth * np.sqrt(2 * np.pi))

    color_values = np.asarray([-np.log10(row["significance"]) for row in rows])
    color_min = float(np.floor(color_values.min()))
    color_max = float(np.ceil(color_values.max()))
    if color_max <= color_min:
        color_max = color_min + 1
    ridge_colors = PlotTheme.get_plot_colors(
        style,
        palette,
        default=["#08060D", "#35105F", "#741287", "#B7277B", "#E04B68", "#F58B63", "#FFF2A6"],
        sequential=True,
    )
    cmap = mcolors.LinearSegmentedColormap.from_list("gsea_ridge", ridge_colors)
    hit_color = ridge_colors[-1]
    norm = mcolors.Normalize(color_min, color_max)

    with PlotTheme.context(style or "nature", palette):
        plt.rcParams["font.family"] = "sans-serif"
        plt.rcParams["font.sans-serif"] = ["Arial", "DejaVu Sans"]
        plt.rcParams["svg.fonttype"] = "none"
        plt.rcParams["pdf.fonttype"] = 42
        if figsize is None:
            figsize = (8.8, max(4.8, 1.4 + 0.48 * len(rows)))
        fig, ax = plt.subplots(figsize=figsize, dpi=dpi)
        for index, row in enumerate(rows):
            baseline = len(rows) - index
            ridge = density(row["values"])
            ridge = ridge / ridge.max() * 0.78 if ridge.max() > 0 else ridge
            color = cmap(norm(-np.log10(row["significance"])))
            ax.hlines(baseline, *x_limits, color="black", linewidth=0.7)
            ax.fill_between(x_grid, baseline, baseline + ridge, color=color, linewidth=0)
            ax.plot(x_grid, baseline + ridge, color="black", linewidth=0.72)
            ax.vlines(row["values"], baseline - 0.105, baseline + 0.025, color=hit_color, linewidth=0.72)
            label = "\n".join(textwrap.wrap(row["description"], width=47, break_long_words=False))
            ax.text(x_limits[0] - 0.02 * (x_limits[1] - x_limits[0]), baseline, label,
                    ha="right", va="center", fontsize=10.2)

        ax.set_xlim(x_limits)
        ax.set_ylim(0.52, len(rows) + 0.88)
        ax.set_yticks([])
        ax.set_xlabel("Ranking statistic", fontsize=12.2)
        ax.tick_params(axis="x", labelsize=10.2)
        for spine in ax.spines.values():
            spine.set_linewidth(0.72)

        sm = plt.cm.ScalarMappable(norm=norm, cmap=cmap)
        sm.set_array([])
        cax = fig.add_axes([0.84, 0.27, 0.022, 0.48])
        cbar = fig.colorbar(sm, cax=cax)
        cbar.set_label(f"-log10({sig_col})", fontsize=11.2, labelpad=8)
        cbar.ax.tick_params(labelsize=9.2)
        apply_figure_style(fig, style, axes=[ax], grid_axis="x", border="style")
        fig.subplots_adjust(left=0.34, right=0.79, top=0.98, bottom=0.12)
        if output_file:
            fig.savefig(output_file, dpi=dpi, facecolor="white")
            plt.close(fig)
    return fig


def plot_gsea_lollipop(
    results_df: pd.DataFrame,
    top_n: int = 20,
    title: str = "GSEA Lollipop Plot",
    output_file: str = None,
    figsize: tuple = None,
    style: Optional[str] = None,
    palette: Optional[str] = None,
    dpi: int = 300,
) -> matplotlib.figure.Figure:
    """Plot GSEA terms as signed normalized-enrichment-score lollipops."""
    df, metric_label = _prepare_lollipop_df(results_df, top_n)
    df = df.copy()
    df["fdr_score"] = -np.log10(np.maximum(df["fdr"].astype(float), np.finfo(float).tiny))
    df["y"] = np.arange(len(df), 0, -1)

    score_min = float(df["fdr_score"].min())
    score_max = float(df["fdr_score"].max())
    constant_fdr = np.isclose(score_min, score_max)
    fdr_breaks = np.array([score_min]) if constant_fdr else np.linspace(score_min, score_max, 6)
    if constant_fdr:
        score_min -= 0.5
        score_max += 0.5

    if figsize is None:
        figsize = (9.2, max(4.6, 1.35 + 0.40 * len(df)))

    fdr_colors = PlotTheme.get_plot_colors(
        style,
        palette,
        default=["#234DA0", "#1F7AB5", "#2FA4C2", "#5FC1C0", "#A0DAB8", "#D6EFB3"],
        sequential=True,
    )
    cmap = mcolors.LinearSegmentedColormap.from_list("gsea_lollipop_fdr", fdr_colors)
    norm = mcolors.Normalize(vmin=score_min, vmax=score_max)
    fdr_labels = _fdr_tick_labels(fdr_breaks)
    signed_metric = metric_label == "NES" and (df["metric"] < 0).any()
    x_start, x_bg_right, xlim_low, xlim_high = _compute_lollipop_x_layout(df["metric"], signed_metric)
    point_sizes = _size_area(df["gene_count"])

    if figsize is None:
        figsize = (9.2, max(4.2, 1.25 + 0.38 * len(df)))

    with PlotTheme.context(style or "nature", palette):
        # ponytail: local fallback avoids noisy Helvetica warnings on machines without that font.
        plt.rcParams["font.family"] = "sans-serif"
        plt.rcParams["font.sans-serif"] = ["Arial", "DejaVu Sans"]
        plt.rcParams["svg.fonttype"] = "none"
        plt.rcParams["pdf.fonttype"] = 42
        fig = plt.figure(figsize=figsize, constrained_layout=True)
        gs = fig.add_gridspec(1, 2, width_ratios=[1.0, 0.22], wspace=0.08)
        ax = fig.add_subplot(gs[0, 0])
        side_gs = gs[0, 1].subgridspec(4, 1, height_ratios=[0.22, 0.30, 0.30, 0.18], hspace=0.30)
        cbar_gs = side_gs[1, 0].subgridspec(
            1, 2, width_ratios=[0.40, 0.60], wspace=0.10
        )
        cbar_ax = fig.add_subplot(cbar_gs[0, 0])
        fig.add_subplot(cbar_gs[0, 1]).axis("off")
        legend_ax = fig.add_subplot(side_gs[2, 0])

        for _, row in df.iterrows():
            left = xlim_low if signed_metric else x_start
            right = xlim_high if signed_metric else x_bg_right
            ax.add_patch(
                Rectangle(
                    (left, row["y"] - 0.45),
                    right - left,
                    0.90,
                    facecolor=cmap(norm(row["fdr_score"])),
                    edgecolor="none",
                    alpha=0.12,
                    zorder=0,
                )
            )
            color = cmap(norm(row["fdr_score"]))
            ax.hlines(row["y"], x_start, row["metric"], color=color, linewidth=1.8, zorder=1)
            ax.hlines(
                row["y"], x_start, row["metric"], color="white",
                linewidth=1.8, alpha=0.18, zorder=2,
            )
            ax.vlines(
                row["metric"],
                row["y"] - 0.26,
                row["y"] + 0.26,
                color="#8C8C8C",
                linewidth=0.45,
                linestyle=(0, (2, 2)),
                alpha=0.80,
                zorder=3,
            )

        ax.scatter(
            df["metric"],
            df["y"],
            s=point_sizes,
            c=[cmap(norm(v)) for v in df["fdr_score"]],
            edgecolors="#8C8C8C",
            linewidths=1.1,
            alpha=1.0,
            zorder=4,
        )

        if signed_metric:
            ax.axvline(0, color="#8C8C8C", linewidth=0.8, zorder=0)
        ax.set_xlim(xlim_low, xlim_high)
        ax.xaxis.set_major_locator(MaxNLocator(nbins=7))
        ax.set_ylim(0.5, len(df) + 0.5)
        ax.set_yticks(df["y"])
        ax.set_yticklabels(df["term_wrap"], ha="right", color="black")
        ax.set_xlabel(metric_label, labelpad=8)
        ax.set_ylabel("")
        ax.set_title(title, fontsize=11, pad=6)
        ax.grid(axis="x", color="#D3D3D3", linewidth=0.35, linestyle="--")
        ax.grid(axis="y", visible=False)
        ax.tick_params(axis="both", length=0, colors="black")
        for spine in ax.spines.values():
            spine.set_visible(False)

        sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
        sm.set_array([])
        cbar = fig.colorbar(sm, cax=cbar_ax, ticks=fdr_breaks)
        cbar.ax.set_title("FDR", pad=6, fontsize=9)
        cbar.ax.set_yticklabels(fdr_labels)
        cbar.outline.set_visible(False)
        cbar.ax.tick_params(length=0, labelsize=8.5)

        legend_counts = _nice_count_breaks(df["gene_count"], n_breaks=3)
        max_count_for_legend = max(float(df["gene_count"].max()), max(legend_counts))
        max_diameter_pt = 8.8 / 25.4 * 72.0
        max_area = max_diameter_pt ** 2
        legend_ax.set_xlim(0, 1)
        legend_ax.set_ylim(0, 1)
        legend_ax.axis("off")
        legend_ax.text(0.0, 0.98, "Gene count", ha="left", va="top", fontsize=9)
        y_positions = [0.50] if len(legend_counts) == 1 else np.linspace(0.68, 0.18, len(legend_counts))
        for count, y_pos in zip(legend_counts, y_positions):
            size = max(count / max_count_for_legend * max_area, 10.0)
            legend_ax.scatter(
                [0.25], [y_pos], s=size,
                facecolors=[cmap(norm(score_max))], edgecolors="#8C8C8C", linewidths=1.1,
            )
            legend_ax.text(0.50, y_pos, str(count), ha="left", va="center", fontsize=8.5)

        apply_figure_style(fig, style, axes=[ax], grid_axis="x", border="style")
        _save_figure(fig, output_file, dpi=dpi)
    return fig
