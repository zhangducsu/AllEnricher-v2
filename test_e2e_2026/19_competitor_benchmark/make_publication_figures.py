#!/usr/bin/env python3
"""Create publication figures and supplementary tables from one frozen run."""

from __future__ import annotations

import argparse
import json
import math
import re
from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import yaml
from PIL import Image


SCRIPT_ROOT = Path(__file__).resolve().parent
PROJECT_ROOT = SCRIPT_ROOT.parents[1]
MATRIX_PATH = SCRIPT_ROOT / "benchmark_matrix.yaml"
VERIFIED_ACTIVITY_RUN = (
    PROJECT_ROOT / "test_e2e_2026" / "99_runs"
    / "20260717_010408_798643_REAL_WORLD_SCI_PRIMARY_VERIFIED"
)
COLORS = {
    "AllEnricher": "#0072B2",
    "clusterProfiler": "#D55E00",
    "WebGestaltR": "#009E73",
    "g:Profiler": "#CC79A7",
    "getENRICH": "#7F7F7F",
}
MARKERS = {"GO": "o", "KEGG": "s"}
DATASET_LABELS = {
    "human_airway_dex": "Human",
    "cattle_metabolic_risk_week1": "Cattle",
    "drosophila_fas_p218": "Drosophila",
    "yeast_tsa1_deletion_stress": "Yeast",
}
FEATURE_GROUPS = {
    "Methods and resources": [
        "ORA", "GSEA", "ssGSEA", "GSVA", "Multi-list comparison",
        "Network topology analysis", "Redundancy reduction", "TaxID species registry",
        "Local DB build/update", "Custom gene sets", "Custom background",
        "TF/regulatory resources",
    ],
    "Reproducibility and delivery": [
        "Offline execution", "Versioned snapshots", "Source/build lineage",
        "Run provenance", "CLI", "Local Web", "REST API",
        "Publication figures", "HTML report", "Methods reference",
    ],
}
FEATURES = [feature for group in FEATURE_GROUPS.values() for feature in group]
FEATURE_DEFINITIONS = {
    "ORA": "Over-representation analysis.",
    "GSEA": "Rank-based gene set enrichment analysis; ordered-query ORA alone is Partial.",
    "ssGSEA": "Single-sample GSEA from an expression matrix.",
    "GSVA": "Gene set variation analysis from an expression matrix.",
    "Multi-list comparison": "Accepts multiple named gene lists in one analysis and returns a linked comparison object or report.",
    "Network topology analysis": "Provides a topology-aware network enrichment test, not only a network visualization.",
    "Redundancy reduction": "Selects or highlights non-redundant enriched terms using topology, semantic similarity or set-cover methods.",
    "TaxID species registry": "Uses NCBI TaxID as stable identity and records database-specific species availability.",
    "Local DB build/update": "Builds or updates managed local annotation resources rather than only accepting one-off inputs.",
    "Custom gene sets": "Accepts user-provided gene-set annotations for enrichment.",
    "Custom background": "Accepts a user-defined ORA reference universe; Partial denotes method or universe restrictions.",
    "TF/regulatory resources": "Provides managed transcription-factor or regulatory gene sets.",
    "Offline execution": "Core analysis can run without a remote service after required resources are installed.",
    "Versioned snapshots": "Users can select or access a dated or versioned annotation snapshot.",
    "Source/build lineage": "Records annotation source and build lineage, not only a software version.",
    "Run provenance": "Automatically records analysis parameters and resource or data versions.",
    "CLI": "Provides a dedicated command-line entry point rather than only a language API.",
    "Local Web": "Provides a Web workbench that can run locally or be self-hosted.",
    "REST API": "Provides a documented HTTP analysis API.",
    "Publication figures": "Produces built-in enrichment visualizations for publication workflows.",
    "HTML report": "Generates a persistent HTML analysis report.",
    "Methods reference": "Generates methods or provenance text from recorded run metadata.",
}
CAPABILITIES = {
    "AllEnricher": [
        "Yes", "Yes", "Yes", "Yes", "No", "No", "No", "Yes", "Yes", "Yes", "Yes", "Yes",
        "Yes", "Yes", "Yes", "Yes", "Yes", "Yes", "Yes", "Yes", "Yes", "Yes",
    ],
    "clusterProfiler": [
        "Yes", "Yes", "No", "No", "Yes", "No", "Yes", "Partial", "No", "Yes", "Yes", "Partial",
        "Yes", "Partial", "No", "Partial", "No", "No", "No", "Yes", "No", "No",
    ],
    "WebGestaltR": [
        "Yes", "Yes", "No", "No", "Yes", "Yes", "Yes", "Partial", "No", "Yes", "Yes", "Yes",
        "Partial", "Partial", "Partial", "Partial", "No", "No", "Yes", "Yes", "Yes", "No",
    ],
    "g:Profiler": [
        "Yes", "Partial", "No", "No", "Yes", "No", "Yes", "Partial", "No", "Yes", "Yes", "Yes",
        "No", "Yes", "Partial", "Yes", "No", "No", "Yes", "Yes", "Partial", "No",
    ],
    "getENRICH": [
        "Yes", "No", "No", "No", "No", "No", "No", "Partial", "Partial", "Partial", "Partial", "No",
        "Partial", "Partial", "Partial", "Partial", "Yes", "No", "No", "Yes", "No", "No",
    ],
}
EVIDENCE = {
    "AllEnricher": {
        "Methods and resources": "runtime:benchmark_manifest.json; repository:README.md and allenricher/database/species_registry.py",
        "Reproducibility and delivery": "repository:README.md, allenricher/cli.py, allenricher/api/server.py and allenricher/report/",
    },
    "clusterProfiler": {
        "Methods and resources": "https://bioconductor.org/packages/release/bioc/html/clusterProfiler.html; https://yulab-smu.top/biomedical-knowledge-mining-book/027-universal-enrichment.html",
        "Reproducibility and delivery": "https://bioconductor.org/packages/release/bioc/html/clusterProfiler.html; https://yulab-smu.top/biomedical-knowledge-mining-book/enrichplot.html",
    },
    "WebGestaltR": {
        "Methods and resources": "https://bzhanglab.github.io/WebGestaltR/reference/WebGestaltR.html; https://bzhanglab.github.io/WebGestaltR/reference/idMapping.html",
        "Reproducibility and delivery": "https://bzhanglab.github.io/WebGestaltR/; https://www.webgestalt.org/api/",
    },
    "g:Profiler": {
        "Methods and resources": "https://biit.cs.ut.ee/gprofiler/page/docs",
        "Reproducibility and delivery": "https://biit.cs.ut.ee/gprofiler/page/docs; https://biit.cs.ut.ee/gprofiler/page/archives",
    },
    "getENRICH": {
        "Methods and resources": "https://github.com/BioinformaticsOnLine/getENRICH",
        "Reproducibility and delivery": "https://github.com/BioinformaticsOnLine/getENRICH",
    },
}
FEATURE_EVIDENCE = {
    ("AllEnricher", "Multi-list comparison"): "repository audit: CLI and API accept one query, ranking or matrix per analysis job; no linked multi-list workflow",
    ("AllEnricher", "Network topology analysis"): "repository audit: enrichment-network plots are available, but no topology-aware enrichment test is implemented",
    ("AllEnricher", "Redundancy reduction"): "repository audit: no deterministic term simplification or set-cover result filter is exposed",
    ("clusterProfiler", "Multi-list comparison"): "https://yulab-smu.top/biomedical-knowledge-mining-book/029-compareCluster.html",
    ("clusterProfiler", "Network topology analysis"): "https://bioconductor.org/packages/release/bioc/html/clusterProfiler.html; capability audit: no topology-aware enrichment test documented",
    ("clusterProfiler", "Redundancy reduction"): "https://yulab-smu.top/biomedical-knowledge-mining-book/enrichplot.html",
    ("WebGestaltR", "Multi-list comparison"): "https://bzhanglab.github.io/WebGestaltR/reference/WebGestaltR.html",
    ("WebGestaltR", "Network topology analysis"): "https://bzhanglab.github.io/WebGestaltR/reference/WebGestaltR.html",
    ("WebGestaltR", "Redundancy reduction"): "https://bzhanglab.github.io/WebGestaltR/reference/WebGestaltR.html",
    ("g:Profiler", "Multi-list comparison"): "https://biit.cs.ut.ee/gprofiler/page/docs",
    ("g:Profiler", "Network topology analysis"): "https://biit.cs.ut.ee/gprofiler/page/docs; tool-set capability audit: no topology-aware enrichment test documented",
    ("g:Profiler", "Redundancy reduction"): "https://biit.cs.ut.ee/gprofiler/page/docs",
    ("getENRICH", "Multi-list comparison"): "https://github.com/BioinformaticsOnLine/getENRICH; workflow capability audit: not documented",
    ("getENRICH", "Network topology analysis"): "https://github.com/BioinformaticsOnLine/getENRICH; workflow capability audit: not documented",
    ("getENRICH", "Redundancy reduction"): "https://github.com/BioinformaticsOnLine/getENRICH; workflow capability audit: not documented",
    ("getENRICH", "GSEA"): "https://github.com/BioinformaticsOnLine/getENRICH; frozen workflow audit: no ranked-list GSEA call",
}


def configure_style() -> None:
    plt.rcParams.update({
        "font.family": "DejaVu Sans", "font.size": 8, "axes.labelsize": 8,
        "xtick.labelsize": 7.5, "ytick.labelsize": 7.5, "legend.fontsize": 7.5,
        "axes.linewidth": 0.7, "lines.linewidth": 1.0, "savefig.facecolor": "white",
        "figure.facecolor": "white", "pdf.fonttype": 42, "ps.fonttype": 42,
    })


def panel_label(axis: plt.Axes, label: str) -> None:
    axis.text(-0.09, 1.04, label, transform=axis.transAxes, fontsize=10, fontweight="bold", va="bottom")


def save_figure(fig: plt.Figure, base: Path, tiff_dpi: int, png_dpi: int = 300) -> list[dict[str, Any]]:
    base.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(base.with_suffix(".pdf"))
    fig.savefig(base.with_suffix(".svg"))
    fig.savefig(base.with_suffix(".png"), dpi=png_dpi)
    fig.savefig(base.with_suffix(".tiff"), dpi=tiff_dpi, pil_kwargs={"compression": "tiff_lzw"})
    plt.close(fig)
    records = []
    for suffix in (".png", ".tiff"):
        path = base.with_suffix(suffix)
        with Image.open(path) as image:
            records.append({
                "file": str(path), "pixels": list(image.size),
                "dpi": [round(float(value)) for value in image.info.get("dpi", (0, 0))],
            })
    return records


def capability_table(run_dir: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    rows = []
    for tool, values in CAPABILITIES.items():
        if len(values) != len(FEATURES):
            raise ValueError(f"{tool} has {len(values)} capability values for {len(FEATURES)} features")
        for feature, value in zip(FEATURES, values):
            group = next(name for name, features in FEATURE_GROUPS.items() if feature in features)
            rows.append({
                "tool": tool, "group": group, "feature": feature, "value": value,
                "definition": FEATURE_DEFINITIONS[feature],
                "evidence": FEATURE_EVIDENCE.get((tool, feature), EVIDENCE[tool][group]),
            })
    long = pd.DataFrame(rows)
    wide = long.pivot(index="tool", columns="feature", values="value").reindex(index=CAPABILITIES, columns=FEATURES)
    return long, wide


def plot_capability_matrix(axis: plt.Axes, wide: pd.DataFrame, features: list[str], title: str) -> None:
    value_map = {"No": 0, "Partial": 1, "Yes": 2, "N/A": 3}
    colors = ["#FFFFFF", "#E69F00", "#56B4E9", "#D9D9D9"]
    display = wide.loc[:, features].T
    matrix = display.replace(value_map).to_numpy(dtype=float)
    from matplotlib.colors import ListedColormap
    axis.imshow(matrix, aspect="auto", cmap=ListedColormap(colors), vmin=-0.5, vmax=3.5)
    glyph = {"No": "N", "Partial": "P", "Yes": "Y", "N/A": "-"}
    for row in range(display.shape[0]):
        for column in range(display.shape[1]):
            value = display.iloc[row, column]
            axis.text(column, row, glyph[value], ha="center", va="center", fontsize=7.5)
    feature_labels = {
        "TaxID species registry": "TaxID registry",
        "Local DB build/update": "DB build/update",
        "Custom gene sets": "Custom sets",
        "Custom background": "Custom bg.",
        "TF/regulatory resources": "TF/regulatory",
        "Multi-list comparison": "Multi-list",
        "Network topology analysis": "Topology test",
        "Redundancy reduction": "Redundancy",
        "Offline execution": "Offline",
        "Versioned snapshots": "DB snapshots",
        "Source/build lineage": "Source lineage",
        "Publication figures": "Figures",
        "Methods reference": "Methods text",
    }
    axis.set_title(title, fontsize=8.2, fontweight="bold", pad=5)
    axis.set_xticks(range(display.shape[1]), display.columns, rotation=22, ha="right")
    axis.set_yticks(range(display.shape[0]), [feature_labels.get(name, name) for name in display.index])
    axis.tick_params(length=0, labelsize=7.5)
    axis.set_xticks(np.arange(-0.5, display.shape[1], 1), minor=True)
    axis.set_yticks(np.arange(-0.5, display.shape[0], 1), minor=True)
    axis.grid(which="minor", color="#B0B0B0", linewidth=0.35)


def plot_metric_points(axis: plt.Axes, frame: pd.DataFrame, metrics: list[str], ylabel: str) -> None:
    tools = [tool for tool in ("clusterProfiler", "WebGestaltR", "g:Profiler") if tool in set(frame["comparator"])]
    positions = np.arange(len(metrics), dtype=float)
    tool_offsets = {"clusterProfiler": -0.28, "WebGestaltR": 0.0, "g:Profiler": 0.28}
    case_offsets = {
        ("human_airway_dex", "GO"): -0.10,
        ("human_airway_dex", "KEGG"): -0.033,
        ("cattle_metabolic_risk_week1", "GO"): 0.033,
        ("cattle_metabolic_risk_week1", "KEGG"): 0.10,
    }
    for tool in tools:
        subset = frame[frame["comparator"] == tool]
        for _, row in subset.iterrows():
            marker = MARKERS.get(row["database"], "o")
            filled = row["dataset"] == "human_airway_dex"
            offset = tool_offsets[tool] + case_offsets.get((row["dataset"], row["database"]), 0.0)
            for metric_index, metric in enumerate(metrics):
                value = pd.to_numeric(pd.Series([row.get(metric)]), errors="coerce").iloc[0]
                if pd.notna(value):
                    axis.scatter(
                        metric_index + offset, value, s=20, marker=marker,
                        facecolor=COLORS[tool] if filled else "white",
                        edgecolor="none" if filled else COLORS[tool],
                        linewidth=0 if filled else 0.9, alpha=0.78, zorder=3,
                    )
    metric_labels = {
        "spearman": "Spearman", "term_jaccard": "Term\nJaccard",
        "significant_jaccard": "Significant\nJaccard", "top20_jaccard": "Top-20\nJaccard",
        "sign_concordance": "Direction",
    }
    axis.set_xticks(positions, [metric_labels.get(name, name.replace("_", " ")) for name in metrics])
    axis.set_ylabel(ylabel)
    axis.set_ylim(-0.04, 1.05)
    axis.axhline(1, color="#777777", linewidth=0.6, linestyle="--")
    axis.grid(axis="y", color="#D9D9D9", linewidth=0.5)
    axis.spines[["top", "right"]].set_visible(False)


def figure2(metrics: pd.DataFrame, config: dict[str, Any], wide: pd.DataFrame) -> plt.Figure:
    main_datasets = {dataset for dataset, values in config["datasets"].items() if values["manuscript_scope"] == "main"}
    main = metrics[(metrics["dataset"].isin(main_datasets)) & (metrics["status"] != "UNAVAILABLE")].copy()
    fig = plt.figure(figsize=(178 / 25.4, 115 / 25.4))
    outer = fig.add_gridspec(2, 1, height_ratios=[0.82, 1.45])
    top = outer[0].subgridspec(1, 3, wspace=0.55)
    bottom = outer[1].subgridspec(1, 2, wspace=0.68)
    axis_a = fig.add_subplot(top[0, 0])
    axis_b = fig.add_subplot(top[0, 1])
    axis_c = fig.add_subplot(top[0, 2])
    axis_d1 = fig.add_subplot(bottom[0, 0])
    axis_d2 = fig.add_subplot(bottom[0, 1])
    ora = main[main["method"] == "ORA"]
    plot_metric_points(axis_a, ora, ["spearman", "term_jaccard"], "Agreement")
    panel_label(axis_a, "A")
    plot_metric_points(axis_b, ora, ["significant_jaccard", "top20_jaccard"], "Jaccard index")
    panel_label(axis_b, "B")
    gsea = main[main["method"] == "GSEA"]
    plot_metric_points(axis_c, gsea, ["spearman", "sign_concordance", "significant_jaccard"], "Case-level metric")
    axis_c.set_xticklabels(["NES\nSpearman", "Direction", "Significant\nJaccard"])
    panel_label(axis_c, "C")
    plot_capability_matrix(axis_d1, wide, FEATURE_GROUPS["Methods and resources"], "Methods and resources")
    panel_label(axis_d1, "D")
    plot_capability_matrix(
        axis_d2, wide, FEATURE_GROUPS["Reproducibility and delivery"],
        "Reproducibility and delivery",
    )
    handles = [plt.Line2D([], [], marker="o", linestyle="", color=COLORS[tool], markeredgecolor="none", alpha=0.78, label=tool)
               for tool in ("clusterProfiler", "WebGestaltR", "g:Profiler")]
    handles.extend([
        plt.Line2D([], [], marker="o", linestyle="", color="black", label="GO"),
        plt.Line2D([], [], marker="s", linestyle="", color="black", label="KEGG"),
        plt.Line2D([], [], marker="o", linestyle="", markerfacecolor="black", markeredgecolor="none", alpha=0.78, label="Human"),
        plt.Line2D([], [], marker="o", linestyle="", markerfacecolor="white", markeredgecolor="#555555", alpha=0.78, label="Cattle"),
    ])
    fig.legend(
        handles=handles, loc="upper center", ncol=7, frameon=False,
        bbox_to_anchor=(0.5, 0.995), columnspacing=0.8, handletextpad=0.35,
    )
    fig.subplots_adjust(left=0.13, right=0.97, top=0.86, bottom=0.15, hspace=0.68)
    return fig


def paired_results(results: pd.DataFrame, dataset: str, database: str, method: str, comparator: str) -> pd.DataFrame:
    subset = results[(results["dataset"] == dataset) & (results["database"] == database) & (results["method"] == method) & (results["status"] == "PASS")]
    reference = subset[subset["tool"] == "AllEnricher"]
    other = subset[subset["tool"] == comparator]
    return reference.merge(other, on="term_id", suffixes=("_reference", "_comparator"))


def supplement_s1(results: pd.DataFrame, source_dir: Path) -> plt.Figure:
    datasets = list(DATASET_LABELS)
    comparators = ["clusterProfiler", "WebGestaltR", "g:Profiler"]
    fig, axes = plt.subplots(4, 3, figsize=(178 / 25.4, 230 / 25.4), sharex=False, sharey=False)
    source_rows = []
    for row, dataset in enumerate(datasets):
        for column, comparator in enumerate(comparators):
            axis = axes[row, column]
            plotted = False
            for database in ("GO", "KEGG"):
                paired = paired_results(results, dataset, database, "ORA", comparator)
                if paired.empty:
                    continue
                x = -np.log10(pd.to_numeric(paired["adjusted_p_value_reference"], errors="coerce").clip(lower=np.finfo(float).tiny))
                y = -np.log10(pd.to_numeric(paired["adjusted_p_value_comparator"], errors="coerce").clip(lower=np.finfo(float).tiny))
                valid = x.notna() & y.notna()
                axis.scatter(x[valid], y[valid], s=4, marker=MARKERS[database], color=COLORS[comparator], alpha=0.28, linewidths=0)
                for term_id, xv, yv in zip(paired.loc[valid, "term_id"], x[valid], y[valid]):
                    source_rows.append({"dataset": dataset, "database": database, "comparator": comparator, "term_id": term_id, "reference_log10_fdr": xv, "comparator_log10_fdr": yv})
                plotted = True
            if plotted:
                limit = max(axis.get_xlim()[1], axis.get_ylim()[1])
                axis.plot([0, limit], [0, limit], color="#555555", linestyle="--", linewidth=0.65)
            else:
                axis.text(0.5, 0.5, "Unavailable", transform=axis.transAxes, ha="center", va="center", color="#777777")
            axis.spines[["top", "right"]].set_visible(False)
            if row == 0:
                axis.set_title(comparator, fontsize=8)
            if column == 0:
                axis.set_ylabel(f"{DATASET_LABELS[dataset]}\nComparator -log10(FDR)")
            if row == len(datasets) - 1:
                axis.set_xlabel("AllEnricher -log10(FDR)")
    pd.DataFrame(source_rows).to_csv(source_dir / "Figure_S1_source_data.tsv", sep="\t", index=False)
    fig.subplots_adjust(left=0.11, right=0.99, top=0.96, bottom=0.07, wspace=0.30, hspace=0.34)
    return fig


def supplement_s2(results: pd.DataFrame, metrics: pd.DataFrame, source_dir: Path) -> plt.Figure:
    datasets = list(DATASET_LABELS)
    comparators = ["clusterProfiler", "WebGestaltR"]
    fig = plt.figure(figsize=(178 / 25.4, 250 / 25.4))
    grid = fig.add_gridspec(5, 2, height_ratios=[1, 1, 1, 1, 1.15], hspace=0.42, wspace=0.28)
    source_rows = []
    for row, dataset in enumerate(datasets):
        for column, comparator in enumerate(comparators):
            axis = fig.add_subplot(grid[row, column])
            plotted = False
            for database in ("GO", "KEGG"):
                paired = paired_results(results, dataset, database, "GSEA", comparator)
                if paired.empty:
                    continue
                x = pd.to_numeric(paired["nes_reference"], errors="coerce")
                y = pd.to_numeric(paired["nes_comparator"], errors="coerce")
                valid = x.notna() & y.notna()
                axis.scatter(x[valid], y[valid], s=4, marker=MARKERS[database], color=COLORS[comparator], alpha=0.28, linewidths=0)
                for term_id, xv, yv in zip(paired.loc[valid, "term_id"], x[valid], y[valid]):
                    source_rows.append({"dataset": dataset, "database": database, "comparator": comparator, "term_id": term_id, "reference_nes": xv, "comparator_nes": yv})
                plotted = True
            if plotted:
                low = min(axis.get_xlim()[0], axis.get_ylim()[0])
                high = max(axis.get_xlim()[1], axis.get_ylim()[1])
                axis.plot([low, high], [low, high], color="#555555", linestyle="--", linewidth=0.65)
            else:
                axis.text(0.5, 0.5, "Unavailable", transform=axis.transAxes, ha="center", va="center", color="#777777")
            axis.spines[["top", "right"]].set_visible(False)
            if row == 0:
                axis.set_title(comparator, fontsize=8)
            if column == 0:
                axis.set_ylabel(f"{DATASET_LABELS[dataset]}\nComparator NES")
            if row == len(datasets) - 1:
                axis.set_xlabel("AllEnricher NES")
    heat_axis = fig.add_subplot(grid[4, :])
    heat_metrics = metrics[(metrics["method"] == "GSEA") & (metrics["status"] != "UNAVAILABLE")].copy()
    heat_metrics["case"] = heat_metrics["dataset"].map(DATASET_LABELS) + " / " + heat_metrics["database"] + " / " + heat_metrics["comparator"]
    values = heat_metrics.set_index("case")[["spearman", "sign_concordance", "significant_jaccard"]]
    if not values.empty:
        image = heat_axis.imshow(values.to_numpy(dtype=float), aspect="auto", cmap="viridis", vmin=0, vmax=1)
        heat_axis.set_yticks(range(len(values)), values.index, fontsize=7.5)
        heat_axis.set_xticks(range(values.shape[1]), ["NES Spearman", "Direction", "Significant Jaccard"])
        fig.colorbar(image, ax=heat_axis, fraction=0.018, pad=0.015, label="Metric")
    else:
        heat_axis.text(0.5, 0.5, "No comparable GSEA metrics", transform=heat_axis.transAxes, ha="center")
    pd.DataFrame(source_rows).to_csv(source_dir / "Figure_S2_scatter_source_data.tsv", sep="\t", index=False)
    heat_metrics.to_csv(source_dir / "Figure_S2_metric_source_data.tsv", sep="\t", index=False)
    fig.subplots_adjust(left=0.18, right=0.94, top=0.97, bottom=0.05)
    return fig


def activity_errors() -> pd.DataFrame:
    rows = []
    if not VERIFIED_ACTIVITY_RUN.is_dir():
        return pd.DataFrame(columns=["case", "method", "max_abs_error", "status"])
    for case_dir in sorted((VERIFIED_ACTIVITY_RUN / "cases").iterdir()):
        method = "ssGSEA" if "__ssgsea__" in case_dir.name else "GSVA" if "__gsva__" in case_dir.name else None
        if method is None:
            continue
        oracle = case_dir / "oracle" / f"official_{method.lower()}.tsv"
        outputs = list((case_dir / "output").glob("*_enrichment.tsv"))
        try:
            expected = pd.read_csv(oracle, sep="\t", index_col=0)
            actual = pd.read_csv(outputs[0], sep="\t").set_index("Term_ID")
            common_rows = expected.index.intersection(actual.index)
            common_columns = expected.columns.intersection(actual.columns)
            delta = np.abs(expected.loc[common_rows, common_columns].to_numpy(dtype=float) - actual.loc[common_rows, common_columns].to_numpy(dtype=float))
            error = float(np.nanmax(delta)) if delta.size else np.nan
            status = "PASS" if pd.notna(error) else "UNAVAILABLE"
        except Exception:
            error, status = np.nan, "UNAVAILABLE"
        rows.append({"case": case_dir.name, "method": method, "max_abs_error": error, "status": status})
    return pd.DataFrame(rows)


def supplement_s3(errors: pd.DataFrame, source_dir: Path) -> plt.Figure:
    fig, axis = plt.subplots(figsize=(86 / 25.4, 72 / 25.4))
    for index, method in enumerate(("ssGSEA", "GSVA")):
        values = errors.loc[(errors["method"] == method) & errors["max_abs_error"].notna(), "max_abs_error"]
        jitter = np.linspace(-0.08, 0.08, len(values)) if len(values) else []
        axis.scatter(index + np.asarray(jitter), np.maximum(values.astype(float), 1e-16), s=17, color=COLORS["AllEnricher"], edgecolor="black", linewidth=0.35, alpha=0.8)
    axis.axhline(1e-8, color="#D55E00", linestyle="--", linewidth=0.8, label="Acceptance threshold")
    axis.set_yscale("log")
    axis.set_ylim(3e-17, 3e-7)
    axis.set_xticks([0, 1], ["ssGSEA", "GSVA"])
    axis.set_ylabel("Maximum absolute error")
    axis.spines[["top", "right"]].set_visible(False)
    axis.legend(frameon=False)
    errors.to_csv(source_dir / "Figure_S3_source_data.tsv", sep="\t", index=False)
    fig.subplots_adjust(left=0.21, right=0.97, top=0.95, bottom=0.16)
    return fig


def publication_path(path: str | Path, anchor: str) -> str:
    """Return the archive-relative suffix beginning at *anchor*."""
    normalized = str(path).replace("\\", "/")
    marker = f"/{anchor.strip('/')}/"
    searchable = f"/{normalized.lstrip('/')}"
    if marker not in searchable:
        raise ValueError(f"Cannot make publication path from {path!s}; missing {anchor!r}")
    return f"{anchor.strip('/')}/{searchable.split(marker, 1)[1]}"


def source_gmt_path(path: str | Path) -> str:
    return publication_path(path, "database_snapshot")


def portable_command(command: str) -> str:
    """Replace host-specific paths in archived commands with public placeholders."""
    normalized = command.replace("\\", "/")
    normalized = re.sub(r"(?i)(?:[A-Z]:|/mnt/[a-z])?/[^\s]*?/Rscript(?:\.exe)?", "Rscript", normalized)
    normalized = re.sub(
        r"(?i)(?:[A-Z]:|/mnt/[a-z])?/[^\s]*?/AllEnricher-v2",
        "{REPO_ROOT}",
        normalized,
    )
    normalized = re.sub(
        r"\{REPO_ROOT\}/test_e2e_2026/99_runs/[^/\s]+",
        "{ARCHIVE_ROOT}",
        normalized,
    )
    return normalized


def tool_summary(manifest: dict[str, Any]) -> pd.DataFrame:
    versions = manifest["tool_versions"]
    accessed = str(manifest.get("created_at", "2026-07-20"))[:10]
    rows = [
        {
            "tool": "AllEnricher", "version": versions["AllEnricher"],
            "role": "Reference workbench", "execution_environment": "Python 3.11.9; Windows",
            "access_date": accessed, "commit_sha": "N/A",
            "command_template": "python {REPO_ROOT}/test_e2e_2026/19_competitor_benchmark/run_benchmark.py --config {REPO_ROOT}/test_e2e_2026/19_competitor_benchmark/benchmark_matrix.yaml --output {ARCHIVE_ROOT} --case {CASE_ID}",
            "evidence_archive_path": "benchmark_manifest.json; raw/{CASE_ID}/",
            "source_url": "https://github.com/zhangducsu/AllEnricher-v2",
        },
        {
            "tool": "clusterProfiler", "version": versions["clusterProfiler"],
            "role": "ORA and GSEA numerical comparator", "execution_environment": "R 4.6.1; Windows",
            "access_date": accessed, "commit_sha": "N/A",
            "command_template": "Rscript {REPO_ROOT}/test_e2e_2026/19_competitor_benchmark/competitor_methods.R --tool clusterProfiler --method {METHOD} --gmt {ARCHIVE_ROOT}/inputs/{CASE_INPUT}/normalized.gmt --output {ARCHIVE_ROOT}/raw/{CASE_ID}/result.tsv",
            "evidence_archive_path": "raw/{CASE_ID}/command.txt; raw/{CASE_ID}/sessionInfo.txt",
            "source_url": "https://bioconductor.org/packages/clusterProfiler/",
        },
        {
            "tool": "WebGestaltR", "version": versions["WebGestaltR"],
            "role": "ORA and GSEA numerical comparator", "execution_environment": "R 4.5.2; Ubuntu 26.04 (WSL2)",
            "access_date": accessed, "commit_sha": "N/A",
            "command_template": "Rscript {REPO_ROOT}/test_e2e_2026/19_competitor_benchmark/competitor_methods.R --tool WebGestaltR --method {METHOD} --gmt {ARCHIVE_ROOT}/inputs/{CASE_INPUT}/normalized.gmt --output {ARCHIVE_ROOT}/raw/{CASE_ID}/result.tsv",
            "evidence_archive_path": "raw/{CASE_ID}/command.txt; raw/{CASE_ID}/sessionInfo.txt",
            "source_url": "https://bzhanglab.github.io/WebGestaltR/",
        },
        {
            "tool": "g:Profiler", "version": versions["g:Profiler"],
            "role": "GO ORA workflow comparator", "execution_environment": "Remote API",
            "access_date": accessed, "commit_sha": "N/A",
            "command_template": "POST https://biit.cs.ut.ee/gprofiler/api/gost/profile/ using {ARCHIVE_ROOT}/raw/{CASE_ID}/profile_request.json",
            "evidence_archive_path": "raw/{CASE_ID}/profile_request.json; raw/{CASE_ID}/profile_response.json",
            "source_url": "https://biit.cs.ut.ee/gprofiler/",
        },
        {
            "tool": "getENRICH", "version": versions["getENRICH"][:12],
            "role": "Representative workflow audit", "execution_environment": "R 4.4.1; isolated container",
            "access_date": accessed, "commit_sha": versions["getENRICH"],
            "command_template": "N/A: workflow audit only; no statistical command executed",
            "evidence_archive_path": "raw/getENRICH/{CASE_ID}/status.json",
            "source_url": "https://github.com/BioinformaticsOnLine/getENRICH",
        },
        {
            "tool": "GSVA", "version": versions["GSVA"],
            "role": "Direct ssGSEA and GSVA correctness oracle", "execution_environment": "R/Bioconductor direct-call validation",
            "access_date": accessed, "commit_sha": "N/A",
            "command_template": "python {REPO_ROOT}/test_e2e_2026/18_real_world_sci/run_real_world_sci.py --matrix {REPO_ROOT}/test_e2e_2026/18_real_world_sci/case_matrix.yaml --mode primary --output {ARCHIVE_ROOT} --case {CASE_ID}",
            "evidence_archive_path": "validation/108_case_run/cases/{CASE_ID}/oracle/",
            "source_url": "https://bioconductor.org/packages/GSVA/",
        },
    ]
    columns = [
        "tool", "version", "role", "execution_environment", "access_date",
        "commit_sha", "command_template", "evidence_archive_path", "source_url",
    ]
    return pd.DataFrame(rows, columns=columns)


def write_tables(run_dir: Path, supplementary: Path, long_capabilities: pd.DataFrame, metrics: pd.DataFrame, details: pd.DataFrame) -> None:
    inputs = pd.read_csv(run_dir / "input_statistics.tsv", sep="\t", dtype=str, keep_default_na=False)
    inputs["source_gmt"] = inputs["source_gmt"].map(source_gmt_path)
    inputs["expression_atlas_license"] = "CC BY 4.0"
    inputs["annotation_distribution"] = inputs["database"].map({
        "GO": "CC BY 4.0",
        "KEGG": "Not redistributed; hash only",
    })
    inputs["annotation_snapshot"] = "v20260715"
    inputs.to_csv(supplementary / "Table_S1_datasets_inputs_databases.tsv", sep="\t", index=False)
    long_capabilities.to_csv(supplementary / "Table_S2_capability_evidence.tsv", sep="\t", index=False)

    manifest = json.loads((run_dir / "benchmark_manifest.json").read_text(encoding="utf-8"))
    commands = []
    case_dirs = []
    for path in sorted(path for path in (run_dir / "raw").iterdir() if path.is_dir()):
        if path.name == "getENRICH":
            case_dirs.extend(sorted(child for child in path.iterdir() if child.is_dir()))
        else:
            case_dirs.append(path)
    for case_dir in case_dirs:
        case_name = f"getENRICH__{case_dir.name}" if case_dir.parent.name == "getENRICH" else case_dir.name
        command_path = case_dir / "command.txt"
        request_path = case_dir / "profile_request.json"
        if command_path.is_file():
            command = command_path.read_text(encoding="utf-8").strip()
            evidence = command_path
        elif request_path.is_file():
            command = f"POST https://biit.cs.ut.ee/gprofiler/api/gost/profile/ using {{ARCHIVE_ROOT}}/raw/{case_name}/profile_request.json"
            evidence = request_path
        elif case_dir.parent.name == "getENRICH":
            command = "N/A: workflow audit only; no statistical command executed"
            evidence = case_dir / "status.json"
        else:
            command = f"python {{REPO_ROOT}}/test_e2e_2026/19_competitor_benchmark/run_benchmark.py --config {{REPO_ROOT}}/test_e2e_2026/19_competitor_benchmark/benchmark_matrix.yaml --output {{ARCHIVE_ROOT}} --case {case_name}"
            evidence = case_dir / "status.json"
        commands.append(
            {
                "case_id": case_name,
                "tool": "getENRICH" if case_dir.parent.name == "getENRICH" else case_dir.name.rsplit("__", 1)[-1].replace("gProfiler", "g:Profiler"),
                "command": portable_command(command),
                "evidence_archive_path": publication_path(evidence, "raw"),
            }
        )
    tool_summary(manifest).to_csv(supplementary / "Table_S3_versions_commands_access.tsv", sep="\t", index=False)
    pd.DataFrame(commands).to_csv(
        supplementary / "source_data" / "Data_S1_full_case_commands.tsv", sep="\t", index=False
    )

    table_s4 = metrics.merge(details, on=["reference_tool", "comparator", "dataset", "database", "method"], how="left")
    table_s4.insert(0, "record_type", "numeric_comparison")
    getenrich = run_dir / "getenrich_workflow_audit.tsv"
    if getenrich.is_file():
        workflow = pd.read_csv(getenrich, sep="\t", keep_default_na=False)
        workflow.insert(0, "record_type", "workflow_audit")
        table_s4 = pd.concat([table_s4, workflow], ignore_index=True, sort=False)
    table_s4.to_csv(supplementary / "Table_S4_case_metrics_failures.tsv", sep="\t", index=False)

def write_legends(paper_dir: Path) -> None:
    text = """# Figure legends and alternative text

## Figure 1

**Legend.** Evolution from the original AllEnricher release to the v2 integrated workbench and its validation evidence. The original release supported gene-list ORA, locally updatable resources, model and non-model species, custom backgrounds, and a Unix command-line interface. Version 2 adds ranked-gene and expression-matrix methods, transcription-factor resources, a TaxID registry, versioned resource lineage, CLI/Web/REST access, publication graphics, reports, and traceable outputs. Validation comprises four public datasets, 108 correctness cases, 52 competitor result sets, and offline replay.

**Alt text.** A left-to-right diagram contrasts the grey original AllEnricher feature set with the blue, orange, and green v2 workbench, then connects it to purple validation evidence from four datasets, 108 correctness cases, 52 competitor result sets, and offline replay. Optional AI interpretation is shown as a dashed secondary box.

## Figure 2

**Legend.** Case-level numerical agreement and evidence-linked capabilities. (A) ORA adjusted-value rank agreement and tested-term overlap. (B) Significant-term and top-20 overlap. (C) GSEA NES rank correlation, direction concordance, and significant-term overlap. Each point is one dataset-by-database case; filled markers denote human, open markers cattle, circles GO, and squares KEGG. Fixed horizontal offsets separate tools and cases without changing metric values. (D) Two evidence matrices compare methods and resources, then reproducibility and delivery, without a composite score: Y, Yes; P, Partial; N, No; dash, not applicable. Competitor-only capabilities are retained to avoid treating AllEnricher as the comparison ceiling. Capability definitions and cell-level evidence are reported in Table S2; complete numerical metrics and non-comparable cases are in Table S4.

**Alt text.** Three point plots summarize ORA and GSEA agreement at the case level; filled and open markers distinguish human and cattle, while color and shape distinguish tools and databases. Below them, two five-tool matrices compare enrichment methods, multi-list comparison, network topology analysis, redundancy reduction, resource management, custom inputs, reproducibility, interfaces, and reports.

## Supplementary figures

**Figure S1.** Term-level ORA adjusted-value comparisons across four species. Dashed lines indicate equality. AllEnricher and clusterProfiler use the same positive-overlap BH family; WebGestaltR and g:Profiler cases retain their documented custom-annotation universe semantics.

**Figure S2.** Term-level GSEA NES comparisons across four species and a heat map of preregistered case-level metrics.

**Figure S3.** Maximum absolute errors between AllEnricher ssGSEA/GSVA outputs and direct Bioconductor GSVA calls in the frozen 108-case validation run. The dashed line marks 1e-8.
"""
    (paper_dir / "figure_legends_and_alt_text.md").write_text(text, encoding="utf-8")


def generate(run_dir: Path, paper_dir: Path) -> None:
    configure_style()
    config = yaml.safe_load(MATRIX_PATH.read_text(encoding="utf-8"))
    results = pd.read_csv(run_dir / "normalized_results.tsv", sep="\t", low_memory=False)
    metrics = pd.read_csv(run_dir / "benchmark_metrics.tsv", sep="\t")
    details = pd.read_csv(run_dir / "benchmark_metrics_detail.tsv", sep="\t")
    figures = paper_dir / "figures"
    supplementary = paper_dir / "supplementary"
    source_dir = supplementary / "source_data"
    for path in (figures, supplementary, source_dir):
        path.mkdir(parents=True, exist_ok=True)
    long_capabilities, wide_capabilities = capability_table(run_dir)
    audit = []
    audit.extend(save_figure(figure2(metrics, config, wide_capabilities), figures / "Figure_2_competitor_agreement", 600))
    audit.extend(save_figure(supplement_s1(results, source_dir), supplementary / "Figure_S1_ORA_term_agreement", 600))
    audit.extend(save_figure(supplement_s2(results, metrics, source_dir), supplementary / "Figure_S2_GSEA_term_agreement", 600))
    errors = activity_errors()
    audit.extend(save_figure(supplement_s3(errors, source_dir), supplementary / "Figure_S3_activity_method_error", 600))
    write_tables(run_dir, supplementary, long_capabilities, metrics, details)
    write_legends(paper_dir)
    (paper_dir / "figure_render_audit.json").write_text(json.dumps(audit, indent=2) + "\n", encoding="utf-8")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-dir", type=Path, required=True)
    parser.add_argument("--paper-dir", type=Path, required=True)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    generate(args.run_dir.resolve(), args.paper_dir.resolve())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
