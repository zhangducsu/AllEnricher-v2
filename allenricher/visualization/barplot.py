"""Create publication-oriented ORA and GSEA bar plots with optional hierarchy categories."""

import logging
from typing import Dict, Iterable, Mapping, Optional, Tuple

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import numpy as np
import pandas as pd

from .plot_theme import PlotTheme, apply_figure_style, save_figure
from .color_config import PaletteLike, categorical_colors, palette_name_for_role

logger = logging.getLogger(__name__)


# -----------------------------------------------------------------------------
# Database hierarchy colors resolved through ColorConfig.
# -----------------------------------------------------------------------------

def _get_category_colors(database: str, palette: PaletteLike = None) -> Dict[str, str]:
    """Return category colors for a supported hierarchical database."""
    from .color_config import ColorConfig
    
    config = ColorConfig()
    database = database.upper()
    
    if database == "GO":
        return config.get_categorical_colors('go', palette=palette)
    elif database == "KEGG":
        return config.get_categorical_colors('kegg', palette=palette)
    elif database == "REACTOME":
        colors = categorical_colors(palette, n=1)
        return {"default": colors[0]}
    elif database == "DO":
        colors = categorical_colors(palette, n=1)
        return {"default": colors[0]}
    elif database == "DISGENET":
        colors = categorical_colors(palette, n=1)
        return {"default": colors[0]}
    else:
        colors = categorical_colors(palette, n=1)
        return {"default": colors[0]}


def _parse_go_category(term_str: str) -> str:
    """Extract the displayed Gene Ontology namespace from term metadata."""
    # Accept both project hierarchy paths and legacy "namespace: term" labels.
    for sep in ("|", ":"):
        if sep in term_str:
            category = term_str.split(sep)[0].strip().lower()
            # Normalize legacy spaces to the canonical underscore form.
            category = category.replace(" ", "_")
            # Only the three GO namespaces are meaningful category labels.
            valid_categories = {"biological_process", "cellular_component", "molecular_function"}
            if category in valid_categories:
                return category
    return "default"


def _parse_kegg_category(term_str: str) -> str:
    """Extract the displayed KEGG top-level category from term metadata."""
    if "|" in term_str:
        category = term_str.split("|")[0].strip()
        # Normalize category labels before color lookup.
        category = category.replace(" ", "_")
        return category
    return "default"


def _split_hierarchy(value: object) -> list[str]:
    if value is None or pd.isna(value):
        return []
    return [part.strip() for part in str(value).split("|") if part.strip()]


def _hierarchy_level_counts(values: Iterable[object]) -> Dict[int, int]:
    """Count distinct non-terminal hierarchy labels at each level."""
    paths = [_split_hierarchy(value) for value in values]
    paths = [path for path in paths if len(path) >= 2]
    if not paths:
        return {}
    return {
        level: len({path[level] for path in paths if len(path) > level + 1})
        for level in range(max(len(path) - 1 for path in paths))
        if any(len(path) > level + 1 for path in paths)
    }


def _select_hierarchy_level(
    values: Iterable[object], max_categories: int = 6,
) -> Tuple[Optional[int], Dict[int, int]]:
    counts = _hierarchy_level_counts(values)
    selected = next(
        (level for level, count in counts.items() if 2 <= count <= max_categories),
        None,
    )
    return selected, counts


def _category_at_level(value: object, level: int) -> Optional[str]:
    path = _split_hierarchy(value)
    return path[level] if len(path) > level + 1 else None


def _resolve_hierarchy(
    data: pd.DataFrame,
    term_col: str,
    term_id_col: str,
    hierarchy_col: str,
    hierarchy_map: Optional[Mapping[str, str]],
) -> Tuple[pd.Series, list[object]]:
    if hierarchy_map and term_id_col in data.columns:
        normalized = {str(key): value for key, value in hierarchy_map.items()}
        mapped = data[term_id_col].astype(str).map(normalized)
        if hierarchy_col in data.columns:
            mapped = mapped.fillna(data[hierarchy_col])
        return mapped.fillna(""), list(normalized.values())
    if hierarchy_col in data.columns:
        return data[hierarchy_col].fillna(""), data[hierarchy_col].tolist()
    inferred = data[term_col].where(data[term_col].astype(str).str.contains("|", regex=False), "")
    return inferred, inferred.tolist()


def _auto_figsize(n_terms: int, base_width: float = 10.0) -> Tuple[float, float]:
    """Return compact figure dimensions based on term count and label length."""
    # Allocate a bounded row height so large result sets remain reviewable.
    height = max(6.0, min(16.0, n_terms * 0.35))
    return (base_width, height)


def _save_figure(fig: plt.Figure, output_file: Optional[str], dpi: int = 300):
    """Save the figure using the format implied by the output filename."""
    if output_file:
        saved_path = save_figure(fig, output_file, dpi)
        logger.info("Bar plot saved: %s", saved_path)


# -----------------------------------------------------------------------------
# Public enrichment bar-plot function.
# -----------------------------------------------------------------------------

def plot_barplot(
    data: pd.DataFrame,
    database: str = "GO",
    top_n: int = 20,
    style: str = "nature",
    palette: PaletteLike = None,
    title: Optional[str] = None,
    output_file: Optional[str] = None,
    figsize: Optional[Tuple[float, float]] = None,
    dpi: int = 300,
    show_gene_count: bool = True,
    gene_count_col: str = "gene_count",
    qvalue_col: str = "qvalue",
    term_col: str = "term",
    term_id_col: str = "term_id",
    hierarchy_col: str = "hierarchy",
    hierarchy_map: Optional[Mapping[str, str]] = None,
    rich_factor_col: Optional[str] = "rich_factor",
) -> matplotlib.figure.Figure:
    """Plot significant terms as compact horizontal bars."""
    # Never mutate the canonical result table.
    df = data.copy()

    # Validate the normalized plotting schema.
    required_cols = [term_col, qvalue_col, gene_count_col]
    for col in required_cols:
        if col not in df.columns:
            raise ValueError(f"Missing column: {col}")

    hierarchy, hierarchy_population = _resolve_hierarchy(
        df, term_col, term_id_col, hierarchy_col, hierarchy_map,
    )
    df["_hierarchy"] = hierarchy
    hierarchy_level, hierarchy_counts = _select_hierarchy_level(hierarchy_population)
    if hierarchy_counts:
        logger.info(
            "ORA Level Classification Statistics: %s",
            ", ".join(
                f"level {level + 1}={count}" for level, count in hierarchy_counts.items()
            ),
        )

    # Select the most significant terms before plotting.
    df = df.sort_values(by=qvalue_col, ascending=True).head(top_n)

    if len(df) == 0:
        fig, ax = plt.subplots(figsize=figsize or (8, 6))
        ax.text(0.5, 0.5, "No enrichment terms to display", ha="center", va="center")
        ax.set_title(title or f"{database} Enrichment")
        apply_figure_style(fig, style, axes=[ax], border="style")
        _save_figure(fig, output_file, dpi)
        return fig

    # Use capped -log10(q) for a stable horizontal scale.
    df["neg_log10_q"] = -np.log10(df[qvalue_col].clip(lower=1e-300))

    use_hierarchy = hierarchy_level is not None
    if use_hierarchy:
        df["category"] = df["_hierarchy"].map(
            lambda value: _category_at_level(value, hierarchy_level)
        )
        use_hierarchy = df["category"].notna().all()

    # Show both hit count and enrichment factor at the end of each bar.
    if show_gene_count:
        if rich_factor_col and rich_factor_col in df.columns:
            gene_labels = [
                f"{int(gc)}/{rf:.2f}" if pd.notna(rf) else str(int(gc))
                for gc, rf in zip(df[gene_count_col], df[rich_factor_col])
            ]
        else:
            gene_labels = [str(int(gc)) for gc in df[gene_count_col]]
    else:
        gene_labels = None

    # Adapt the canvas to row count and label length.
    if figsize is None:
        figsize = _auto_figsize(len(df))

    # Figure style controls typography and spacing, not data semantics.
    with PlotTheme.context(style or 'nature', palette):
        color_norm = None
        color_cmap = None
        category_colors: Dict[str, str] = {}
        if use_hierarchy:
            category_order = list(dict.fromkeys(
                category for value in hierarchy_population
                if (category := _category_at_level(value, hierarchy_level)) is not None
            ))
            category_colors = dict(zip(
                category_order, categorical_colors(palette, n=len(category_order))
            ))
            df["color"] = df["category"].map(category_colors)
        else:
            color_cmap = PlotTheme.get_sequential_cmap(
                palette_name_for_role(palette, "sequential")
            )
            score_min = float(df["neg_log10_q"].min())
            score_max = float(df["neg_log10_q"].max())
            if score_min == score_max:
                score_min = 0.0
            color_norm = mcolors.Normalize(vmin=score_min, vmax=score_max)
            df["color"] = [color_cmap(color_norm(value)) for value in df["neg_log10_q"]]

        # Build one axes so labels and category legend share a consistent layout.
        fig, ax = plt.subplots(figsize=figsize)

        # Display the most significant term at the top.
        y_positions = np.arange(len(df))
        values = df["neg_log10_q"].values[::-1]
        bar_colors = df["color"].values[::-1]
        terms = df[term_col].values[::-1]

        if gene_labels:
            gene_labels = gene_labels[::-1]

        # Draw horizontal bars with a thin publication-style outline.
        bars = ax.barh(
            y_positions,
            values,
            color=bar_colors,
            edgecolor="none",
            height=0.7,
        )

        # Display descriptive term names without hierarchy prefixes.
        clean_terms = []
        for term in terms:
            term_str = str(term)
            if database.upper() == "GO":
                # GO namespace remains available through color and legend.
                for sep in ("|", ":"):
                    if sep in term_str:
                        parts = term_str.split(sep, 1)
                        category = parts[0].strip().lower().replace(" ", "_")
                        if category in {"biological_process", "cellular_component", "molecular_function"}:
                            term_str = parts[1].strip()
                            break
                clean_terms.append(term_str)
            elif database.upper() == "KEGG" and "|" in term_str:
                # KEGG hierarchy remains available through color and legend.
                parts = term_str.split("|")
                clean_terms.append(parts[-1].strip())
            else:
                clean_terms.append(term_str)

        ax.set_yticks(y_positions)
        ax.set_yticklabels(clean_terms, fontsize=9)

        # Keep numeric annotations just beyond the bar end.
        if gene_labels:
            for i, (bar, label) in enumerate(zip(bars, gene_labels)):
                width = bar.get_width()
                ax.text(
                    width + 0.05,
                    bar.get_y() + bar.get_height() / 2,
                    label,
                    ha="left",
                    va="center",
                    fontsize=8,
                    color="#333333",
                )

        # Label the significance axis explicitly.
        ax.set_xlabel("-log10(Q-value)", fontsize=10)
        max_value = float(np.nanmax(values)) if len(values) else 0.0
        ax.set_xlim(0, max(max_value * 1.25, 1.0))  # Leave room for count labels.

        # Title identifies the database and displayed metrics.
        if title is None:
            title = f"{database} Enrichment"
            if show_gene_count and rich_factor_col and rich_factor_col in data.columns:
                title += " (Gene# / Rich Factor)"
            elif show_gene_count:
                title += " (Gene Count)"
        ax.set_title(title, fontsize=12, fontweight="bold")

        if use_hierarchy:
            present_categories = set(df["category"])
            legend_elements = [
                plt.Rectangle((0, 0), 1, 1, facecolor=color, edgecolor="none", label=cat)
                for cat, color in category_colors.items()
                if cat in present_categories
            ]
            ax.legend(
                handles=legend_elements,
                title=f"Hierarchy level {hierarchy_level + 1}",
                loc="lower right",
                fontsize=8,
                frameon=True,
                fancybox=False,
                edgecolor="#CCCCCC",
            )
        else:
            colorbar = fig.colorbar(
                plt.cm.ScalarMappable(norm=color_norm, cmap=color_cmap),
                ax=ax,
                fraction=0.025,
                pad=0.02,
                shrink=0.28,
                aspect=14,
            )
            colorbar.set_label("-log10(Q-value)", fontsize=9)
            colorbar.ax.tick_params(labelsize=8)

        # Keep only the axes needed to read values.
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.tick_params(labelsize=9)

        # A light x-grid supports value comparison without dominating the figure.
        ax.xaxis.grid(True, linestyle="--", alpha=0.3)
        ax.set_axisbelow(True)

        apply_figure_style(fig, style, axes=[ax], grid_axis="x", border="style")
        plt.tight_layout()
        _save_figure(fig, output_file, dpi)

    return fig
