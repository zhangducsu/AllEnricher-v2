#!/usr/bin/env python3
"""Command-line interface for AllEnricher.

The CLI exposes enrichment analysis, database management, a local Web service,
species discovery, configuration, and transcription-factor workflows.
"""

import argparse
import json
import logging
import re
import sys
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple

import pandas as pd
from jinja2 import Template

from allenricher import __version__
from allenricher.core.config import Config, SPECIES_CONFIGS, database_catalog_entry
from allenricher.core.enrichment import EnrichmentAnalyzer, add_result_term_metadata
from allenricher.database.manager import DatabaseManager, validate_tf_database_species
from allenricher.visualization.color_config import (
    PUBLIC_CATEGORICAL_PALETTES,
    PUBLIC_DIVERGING_PALETTES,
    PUBLIC_PALETTES,
    PUBLIC_SEQUENTIAL_PALETTES,
    PaletteLike,
    resolve_palette_selection,
)
from allenricher.visualization.plotter import Plotter
from allenricher.report.generator import ReportGenerator
from allenricher.ai.interpreter import create_interpreter, create_interpreter_from_config


_TF_DATABASE_NAMES = {'TRRUST', 'CHEA3', 'ANIMALTFDB', 'HTFTARGET'}


def _analysis_species_metadata(species: str, database_dir: str) -> Dict[str, object]:
    """Resolve a recorded species code without inventing unavailable labels."""
    species_config = SPECIES_CONFIGS.get(species)
    if species_config is not None:
        return {
            "species_name": species_config.name,
            "species_taxonomy_id": species_config.taxonomy_id,
        }
    try:
        from allenricher.database.species_registry import SpeciesRegistry

        entry = SpeciesRegistry.load_default(database_dir).query_by_kegg_code(species)
        if entry is not None:
            return {
                "species_name": entry.latin_name,
                "species_taxonomy_id": entry.taxid,
            }
    except Exception:
        logger.debug("Species metadata is unavailable for %s", species, exc_info=True)
    return {}


def _recorded_analysis_parameters(
    config: Config,
    background_mode: str,
    databases: Optional[List[str]] = None,
) -> Dict[str, object]:
    """Record only method parameters that actually affect the analysis."""
    from allenricher.analysis.tf_enrichment import resolve_tf_size_limits

    method = config.method
    limits: Dict[str, Dict[str, object]] = {}
    for database in databases or config.databases:
        is_tf = str(database).upper() in _TF_DATABASE_NAMES
        if is_tf:
            minimum, maximum = resolve_tf_size_limits(
                method,
                getattr(config, 'tf_min_size', None),
                getattr(config, 'tf_max_size', None),
            )
        elif method == 'hypergeometric':
            minimum = config.min_genes
            maximum = config.max_genes
        elif method == 'gsea':
            minimum = 15 if config.gsea_min_size is None else config.gsea_min_size
            maximum = 500 if config.gsea_max_size is None else config.gsea_max_size
        else:
            minimum = 1 if config.gsea_min_size is None else config.gsea_min_size
            maximum = config.gsea_max_size
        limits[str(database)] = {
            "min": minimum,
            "max": "Inf" if maximum is None or maximum == float('inf') else maximum,
        }

    parameters: Dict[str, object] = {"gene_set_size_by_database": limits}
    if method == 'hypergeometric':
        parameters.update({
            "background_mode": "custom" if config.background_file else background_mode,
            "correction": config.correction,
            "pvalue_cutoff": config.pvalue_cutoff,
            "qvalue_cutoff": config.qvalue_cutoff,
        })
    elif method == 'gsea':
        parameters.update({
            "pvalue_cutoff": config.pvalue_cutoff,
            "qvalue_cutoff": config.qvalue_cutoff,
        })
    elif method == 'gsva':
        parameters.update({
            "gsva_method": config.gsva_method,
            "gsva_kcdf": config.gsva_kcdf,
            "gsva_tau": config.gsva_tau,
        })
    elif method == 'ssgsea':
        parameters["ssgsea_tau"] = 0.25
    return parameters


def _recorded_analysis_databases(
    configured_databases: List[str],
    tf_database: Optional[str],
    tf_only: bool,
    has_tf_results: bool,
) -> List[str]:
    """Return databases that actually contributed output to this run."""
    databases = [] if tf_only else list(configured_databases)
    if has_tf_results and tf_database:
        tf_names = ['TRRUST', 'ChEA3'] if tf_database == 'both' else [tf_database]
        canonical = {
            'trrust': 'TRRUST',
            'chea3': 'ChEA3',
            'animaltfdb': 'AnimalTFDB',
            'htftarget': 'hTFtarget',
        }
        databases.extend(canonical.get(str(name).lower(), str(name)) for name in tf_names)
    return list(dict.fromkeys(databases))

# Configure the command-line log format without changing library loggers.
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logging.getLogger('fontTools').setLevel(logging.WARNING)
logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# GSEA, ssGSEA, and GSVA plotting helpers
# ---------------------------------------------------------------------------

# Plot types supported by each analysis method.
_METHOD_PLOT_TYPES = {
    'gsea': {'enrichment', 'enrichment2', 'barplot', 'lollipop', 'ridgeplot', 'emapplot'},
    'ssgsea': {'heatmap', 'group_comparison', 'correlation'},
    'gsva': {'heatmap', 'group_comparison', 'correlation'},
}

_DEFAULT_METHOD_PLOT_TYPES = {
    # emapplot requires R and aPEAR, so it is generated only when requested.
    'gsea': _METHOD_PLOT_TYPES['gsea'] - {'emapplot'},
    'ssgsea': _METHOD_PLOT_TYPES['ssgsea'],
    'gsva': _METHOD_PLOT_TYPES['gsva'],
}

def _parse_gmt_term_data(gmt_file: str) -> Dict[str, Dict[str, object]]:
    """Read gene sets and display metadata from a GMT file.

    Args:
        gmt_file: Path to a ``.gmt`` or ``.gmt.gz`` file. Each row contains a
            gene-set ID, a description, and one or more gene IDs.

    Returns:
        A mapping from gene-set ID to its name, optional hierarchy, and genes.

    Raises:
        FileNotFoundError: If the GMT file does not exist.
    """
    import gzip

    gmt_path = Path(gmt_file)
    if not gmt_path.exists():
        raise FileNotFoundError(f"The GMT file does not exist: {gmt_file}")

    # Open compressed and uncompressed GMT files through the same interface.
    if gmt_path.suffix.lower() == '.gz':
        opener = gzip.open
        mode = 'rt'
    else:
        opener = open
        mode = 'r'

    term_data: Dict[str, Dict[str, object]] = {}
    with opener(gmt_path, mode, encoding='utf-8') as f:
        for line_num, line in enumerate(f, 1):
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            parts = line.split('\t')
            if len(parts) < 3:
                logger.warning(
                    "Skipping malformed GMT row %d: expected at least three tab-separated columns",
                    line_num,
                )
                continue
            set_name = parts[0].strip()
            description = parts[1].strip()
            # Columns from the third onward contain gene identifiers.
            genes = {g.strip() for g in parts[2:] if g.strip()}
            if genes:
                name = description if description.lower() not in {'', 'na', 'n/a', 'null'} else set_name
                term_data[set_name] = {'name': name, 'genes': genes}
                if '|' in name:
                    term_data[set_name]['hierarchy'] = name

    logger.info("Loaded %d gene sets from %s", len(term_data), gmt_file)
    return term_data


def _parse_gmt_file(gmt_file: str) -> Dict[str, Set[str]]:
    """Return a gene-set mapping compatible with legacy callers."""
    return {
        term_id: set(info['genes'])
        for term_id, info in _parse_gmt_term_data(gmt_file).items()
    }


def _parse_groups(groups_str: str) -> Dict[str, List[str]]:
    """Parse a compact sample-group specification.

    Args:
        groups_str: Text such as ``Control:S1,S2;Treatment:S3,S4``.

    Returns:
        A mapping from each group name to its sample names.

    Raises:
        ValueError: If the specification is malformed or assigns a sample to
            more than one group.
    """
    if not groups_str:
        return {}

    groups: Dict[str, List[str]] = {}
    assigned_samples: Dict[str, str] = {}
    for group_def in groups_str.split(';'):
        group_def = group_def.strip()
        if not group_def:
            continue
        if ':' not in group_def:
            raise ValueError(
                f"Grouping definition error: '{group_def}',"
                f"Expected format is 'GroupName: sample1, sample2'"
            )
        group_name, samples_str = group_def.split(':', 1)
        group_name = group_name.strip()
        samples = [s.strip() for s in samples_str.split(',') if s.strip()]
        if not group_name:
            raise ValueError("Group name cannot be empty")
        if not samples:
            raise ValueError(f"Group '{group_name}' does not contain any samples")
        if group_name in groups:
            raise ValueError(f"Duplicate group name: '{group_name}'")
        duplicate_samples = [sample for sample in samples if sample in assigned_samples]
        if duplicate_samples:
            sample = duplicate_samples[0]
            raise ValueError(
                f"Sample '{sample}' is assigned to both "
                f"'{assigned_samples[sample]}' and '{group_name}'"
            )
        groups[group_name] = samples
        assigned_samples.update({sample: group_name for sample in samples})

    logger.info("Parsed %d sample groups: %s", len(groups), list(groups))
    return groups


def _safe_plot_stem(name: str) -> str:
    """Convert a term ID into a cross-platform filename stem."""
    stem = re.sub(r"[^\w.-]+", "_", str(name)).strip("._")
    return stem or "term"


def _with_plot_term_names(
    data: pd.DataFrame,
    term_names: Optional[Dict[str, str]],
    database: Optional[str] = None,
) -> pd.DataFrame:
    """Add display names to a plotting copy without altering result tables."""
    plot_data = data.copy()
    id_col = next(
        (column for column in ("pathway", "Term_ID", "term_id", "ID", "id")
         if column in plot_data.columns),
        None,
    )
    if id_col is None:
        return plot_data

    term_ids = plot_data[id_col].astype(str)
    existing_name_col = next(
        (column for column in ("Term_Name", "Description", "term_name")
         if column in plot_data.columns),
        None,
    )
    fallback = (
        plot_data[existing_name_col].astype(str)
        if existing_name_col is not None else term_ids
    )
    name_map = {str(key): str(value) for key, value in (term_names or {}).items()}
    mapped = term_ids.map(name_map)
    valid = mapped.notna() & mapped.astype(str).str.strip().ne("")

    plot_data["Term_ID"] = term_ids
    plot_data["Term_Name"] = mapped.where(valid, fallback)
    if database:
        plot_data["Database"] = str(database)
    return plot_data


def _normalize_ranked_genes(
    ranked_genes: Optional[list],
    gene_weights: Optional[dict],
) -> Tuple[List[str], Dict[str, float]]:
    """Normalize ranked genes and weights for deterministic GSEA plots."""
    if not ranked_genes:
        return [], gene_weights or {}

    normalized_genes: List[str] = []
    normalized_weights: Dict[str, float] = dict(gene_weights or {})
    for item in ranked_genes:
        if isinstance(item, (tuple, list)) and len(item) >= 2:
            gene = str(item[0])
            normalized_genes.append(gene)
            normalized_weights[gene] = float(item[1])
        else:
            gene = str(item)
            normalized_genes.append(gene)
            normalized_weights.setdefault(gene, 1.0)
    return normalized_genes, normalized_weights


def _calculate_running_es_rows(
    term_id: str,
    term_name: str,
    ranked_genes: List[str],
    gene_weights: Dict[str, float],
    gene_set: Set[str],
) -> List[dict]:
    """Calculate the running enrichment-score trajectory for one gene set."""
    n = len(ranked_genes)
    hits = gene_set & set(ranked_genes)
    nh = len(hits)
    if n == 0 or nh == 0:
        return []

    nr = sum(abs(gene_weights.get(gene, 1.0)) for gene in hits)
    hit_inc = 1.0 / nr if nr > 0 else 0.0
    miss_inc = 1.0 / (n - nh) if (n - nh) > 0 else 0.0
    running_sum = 0.0
    rows = []

    for idx, gene in enumerate(ranked_genes, start=1):
        weight = float(gene_weights.get(gene, 1.0))
        is_hit = gene in gene_set
        if is_hit:
            running_sum += hit_inc * abs(weight)
        else:
            running_sum -= miss_inc
        rows.append({
            "Term_ID": term_id,
            "Term_Name": term_name,
            "Rank": idx,
            "Gene": gene,
            "Weight": weight,
            "Hit": is_hit,
            "Running_ES": running_sum,
        })
    return rows


def _write_running_es_file(
    output_file: Path,
    pathways: pd.DataFrame,
    ranked_genes: Optional[list],
    gene_weights: Optional[dict],
    gene_sets: Optional[Dict[str, Set[str]]],
) -> Optional[str]:
    """Write the running-ES table consumed by the R plotting scripts."""
    normalized_genes, normalized_weights = _normalize_ranked_genes(ranked_genes, gene_weights)
    if not normalized_genes or not gene_sets:
        logger.warning("Ranked genes or gene sets are unavailable; skipping R enrichment plots")
        return None

    term_id_col = next((c for c in ["pathway", "Term_ID", "term_id", "ID", "id"] if c in pathways.columns), None)
    term_name_col = next((c for c in ["Term_Name", "Description", "pathway", "term_name"] if c in pathways.columns), None)
    if not term_id_col:
        logger.warning("GSEA results do not contain a pathway or Term_ID column; skipping R enrichment plots")
        return None

    rows = []
    for _, row in pathways.iterrows():
        term_id = str(row[term_id_col])
        gene_set = gene_sets.get(term_id)
        if not gene_set:
            continue
        term_name = str(row[term_name_col]) if term_name_col else term_id
        rows.extend(_calculate_running_es_rows(
            term_id, term_name, normalized_genes, normalized_weights, gene_set
        ))

    if not rows:
        logger.warning("Selected pathways do not match the supplied gene sets; skipping R enrichment plots")
        return None

    pd.DataFrame(rows).to_csv(output_file, sep="\t", index=False)
    return str(output_file)


def _gene_sets_for_database(
    database: str,
    gene_sets: Optional[Dict[str, Set[str]]],
    gene_sets_by_database: Optional[Dict[str, Dict[str, Set[str]]]],
) -> Optional[Dict[str, Set[str]]]:
    """Return the gene-set namespace that belongs to one database."""
    if gene_sets_by_database:
        target = str(database)
        for name, sets in gene_sets_by_database.items():
            if str(name).casefold() == target.casefold():
                return sets
    return gene_sets


def _select_gsea_enrichment_terms(
    df: pd.DataFrame,
    top_up: int = 5,
    top_down: int = 5,
) -> pd.DataFrame:
    """Select default single-pathway enrichment plots as top NES up/down terms."""
    if df is None or df.empty:
        return pd.DataFrame()
    if top_up <= 0 and top_down <= 0:
        return pd.DataFrame(columns=df.columns)

    nes_col = next((c for c in ['NES', 'nes'] if c in df.columns), None)
    if nes_col is None:
        return df.head(max(0, top_up) + max(0, top_down)).copy()

    sig_col = next((c for c in ['FDR', 'Adjusted_P_Value', 'FDR q-val', 'qvalue', 'padj', 'p.adjust'] if c in df.columns), None)
    pval_col = next((c for c in ['pval', 'p_value', 'P_Value', 'NOM p-val', 'pvalue', 'P-value'] if c in df.columns), None)

    work = df.copy()
    work[nes_col] = pd.to_numeric(work[nes_col], errors='coerce')
    work = work.dropna(subset=[nes_col])
    if work.empty:
        return pd.DataFrame()

    if sig_col:
        q = pd.to_numeric(work[sig_col], errors='coerce')
        filtered = work[q < 0.05].copy()
        if not filtered.empty:
            work = filtered
        elif pval_col:
            p = pd.to_numeric(work[pval_col], errors='coerce')
            filtered = work[p < 0.05].copy()
            if not filtered.empty:
                work = filtered
    elif pval_col:
        p = pd.to_numeric(work[pval_col], errors='coerce')
        filtered = work[p < 0.05].copy()
        if not filtered.empty:
            work = filtered

    up = work[work[nes_col] > 0].sort_values(nes_col, ascending=False).head(max(0, top_up))
    down = work[work[nes_col] < 0].sort_values(nes_col, ascending=True).head(max(0, top_down))
    selected = pd.concat([up, down], ignore_index=False)

    if selected.empty:
        selected = work.reindex(work[nes_col].abs().sort_values(ascending=False).index).head(max(1, top_up + top_down))

    return selected.drop_duplicates().copy()


def _generate_plots_one_format(
    method: str,
    results: dict,
    ranked_genes: Optional[list],
    gene_weights: Optional[dict],
    gene_sets: Optional[Dict[str, Set[str]]],
    gene_sets_by_database: Optional[Dict[str, Dict[str, Set[str]]]],
    expr_matrix,
    groups: Optional[Dict[str, List[str]]],
    plot_types: List[str],
    output_dir: str,
    plot_format: str = 'png',
    plot_dpi: int = 300,
    plot_style: str = 'nature',
    plot_palette: Optional[str] = None,
    use_r_plots: bool = False,
    emapplot_qvalue: float = 0.05,
    emapplot_min_count: int = 3,
    emapplot_top_n: int = 30,
    gsea_enrichment_top_up: int = 5,
    gsea_enrichment_top_down: int = 5,
    gsea_multi_top_up: int = 3,
    gsea_multi_top_down: int = 3,
    activity_heatmap_top_n: int = 40,
    plot_figsize: Optional[Tuple[float, float]] = None,
    term_name_maps: Optional[Dict[str, Dict[str, str]]] = None,
) -> List[str]:
    """Generate the plots requested for an enrichment method.

    Args:
        method: Analysis method: ``gsea``, ``ssgsea``, or ``gsva``.
        results: Mapping from database name to result table.
        ranked_genes: Ranked genes used by GSEA.
        gene_weights: Ranking statistic keyed by gene ID.
        gene_sets: Gene sets used to calculate GSEA trajectories.
        expr_matrix: Expression matrix used by ssGSEA and GSVA.
        groups: Optional mapping from group names to sample names.
        plot_types: Plot types requested by the user.
        output_dir: Analysis output directory.
        plot_format: Output format: ``png``, ``pdf``, or ``svg``.
        plot_dpi: Raster resolution in dots per inch.
        plot_style: Figure style preset.
        plot_palette: Optional palette override.

    Returns:
        Paths to the generated plot files.
    """
    generated_files: List[str] = []

    if not plot_types or not results:
        return generated_files

    # Ignore plot types that do not apply to the selected method.
    supported = _METHOD_PLOT_TYPES.get(method, set())
    for pt in plot_types:
        if pt not in supported:
            logger.warning(
                "Plot type '%s' is not supported for %s and will be skipped. "
                "Supported types: %s",
                pt,
                method,
                sorted(supported),
            )

    valid_types = [pt for pt in plot_types if pt in supported]
    if not valid_types:
        return generated_files

    if method == 'gsea':
        ranked_genes, gene_weights = _normalize_ranked_genes(ranked_genes, gene_weights)

    # Keep all method-specific plots in a predictable subdirectory.
    plot_dir = Path(output_dir) / "gsea_plots"
    plot_dir.mkdir(parents=True, exist_ok=True)

    # GSEA plots
    if method == 'gsea':
        if use_r_plots:
            # R plotting mode
            from allenricher.visualization.r_plotter import (
                check_r_environment,
                plot_gsea_barplot_r, plot_gsea_lollipop_r,
                plot_gsea_ridgeplot_r, plot_gsea_emapplot_r,
                plot_gsea_enrichment_r, plot_gsea_enrichment2_r,
            )

            if not check_r_environment():
                logger.warning("R environment not found, falling back to Python plotting")
                use_r_plots = False  # fallback
            else:
                # Write the standardized table consumed by R scripts.
                for db_name, df in results.items():
                    if df is None or len(df) == 0:
                        continue
                    db_gene_sets = _gene_sets_for_database(
                        db_name, gene_sets, gene_sets_by_database
                    )
                    name_map = (term_name_maps or {}).get(str(db_name).upper(), {})
                    df = _with_plot_term_names(df, name_map, db_name)
                    tsv_path = str(plot_dir / f"{db_name}_enrichment.tsv")
                    df.to_csv(tsv_path, sep='\t', index=False)

                    running_es_path: Optional[str] = None
                    top_pathways_for_enrichment: Optional[pd.DataFrame] = None
                    top_pathways_for_enrichment2: Optional[pd.DataFrame] = None
                    nes_col = next((c for c in ['NES', 'nes'] if c in df.columns), None)
                    needs_running_es = any(
                        pt in valid_types for pt in ['ridgeplot', 'enrichment', 'enrichment2']
                    )
                    if needs_running_es and nes_col:
                        _abs_col = f'_{nes_col}_abs'
                        df_for_top = df.copy()
                        df_for_top[_abs_col] = df_for_top[nes_col].abs()
                        selected_for_es = []
                        if 'ridgeplot' in valid_types:
                            sig_col = next(
                                (c for c in ['padj', 'pval', 'p_value', 'pvalue', 'P_Value', 'NOM p-val',
                                             'Adjusted_P_Value', 'FDR', 'p.adjust', 'qvalue']
                                 if c in df_for_top.columns),
                                None,
                            )
                            ridge_terms = (
                                df_for_top.sort_values([sig_col, _abs_col], ascending=[True, False]).head(10)
                                if sig_col else df_for_top.nlargest(10, _abs_col)
                            )
                            selected_for_es.append(ridge_terms.drop(columns=[_abs_col]))
                        if 'enrichment' in valid_types:
                            top_pathways_for_enrichment = _select_gsea_enrichment_terms(
                                df,
                                top_up=gsea_enrichment_top_up,
                                top_down=gsea_enrichment_top_down,
                            )
                            selected_for_es.append(top_pathways_for_enrichment)
                        if 'enrichment2' in valid_types:
                            top_pathways_for_enrichment2 = _select_gsea_enrichment_terms(
                                df,
                                top_up=gsea_multi_top_up,
                                top_down=gsea_multi_top_down,
                            )
                            selected_for_es.append(top_pathways_for_enrichment2)
                        top_pathways_for_es = (
                            pd.concat(selected_for_es, ignore_index=True)
                            .drop_duplicates()
                            if selected_for_es else pd.DataFrame()
                        )
                        running_es_path = _write_running_es_file(
                            plot_dir / f"{db_name}_running_es.tsv",
                            top_pathways_for_es,
                            ranked_genes,
                            gene_weights,
                            db_gene_sets,
                        )

                    # barplot
                    if 'barplot' in valid_types:
                        out_file = str(plot_dir / f"{db_name}_barplot.{plot_format}")
                        if plot_gsea_barplot_r(
                            tsv_path, out_file, top_n=20, dpi=plot_dpi,
                            style=plot_style, palette=plot_palette,
                        ):
                            generated_files.append(out_file)

                    # lollipop
                    if 'lollipop' in valid_types:
                        out_file = str(plot_dir / f"{db_name}_lollipop.{plot_format}")
                        if plot_gsea_lollipop_r(
                            tsv_path, out_file, top_n=20, dpi=plot_dpi,
                            style=plot_style, palette=plot_palette,
                        ):
                            generated_files.append(out_file)

                    # ridgeplot
                    if 'ridgeplot' in valid_types:
                        out_file = str(plot_dir / f"{db_name}_ridgeplot.{plot_format}")
                        if running_es_path and plot_gsea_ridgeplot_r(
                            tsv_path,
                            out_file,
                            top_n=10,
                            running_es_path=running_es_path or "",
                            dpi=plot_dpi,
                            style=plot_style,
                            palette=plot_palette,
                        ):
                            generated_files.append(out_file)
                        elif not running_es_path:
                            logger.warning(
                                "GSEA ridgesplent needs a match of rankedgens%s",
                                db_name,
                            )

                    # emapplot
                    if 'emapplot' in valid_types:
                        out_file = str(plot_dir / f"{db_name}_emapplot.{plot_format}")
                        if plot_gsea_emapplot_r(
                            tsv_path,
                            out_file,
                            top_n=emapplot_top_n,
                            qvalue=emapplot_qvalue,
                            min_count=emapplot_min_count,
                            dpi=plot_dpi,
                            style=plot_style,
                            palette=plot_palette,
                        ):
                            generated_files.append(out_file)

                    # Running-ES plots for selected single and multi-pathway views.
                    if 'enrichment' in valid_types or 'enrichment2' in valid_types:
                        frames = [
                            frame for frame in (top_pathways_for_enrichment, top_pathways_for_enrichment2)
                            if frame is not None and not frame.empty
                        ]
                        top_pathways = pd.concat(frames, ignore_index=True).drop_duplicates() if frames else None
                        if top_pathways is not None and not top_pathways.empty:
                            running_es_path = _write_running_es_file(
                                plot_dir / f"{db_name}_running_es.tsv",
                                top_pathways,
                                ranked_genes,
                                gene_weights,
                                db_gene_sets,
                            )
                            if running_es_path:
                                single_term_ids = {
                                    str(row.get('pathway', row.get('Term_ID', row.get('term_id', ''))))
                                    for _, row in (top_pathways_for_enrichment.iterrows()
                                                   if top_pathways_for_enrichment is not None
                                                   else [])
                                }
                                for _, row in top_pathways.iterrows():
                                    term_id = row.get('pathway', row.get('Term_ID', row.get('term_id', '')))
                                    if not term_id:
                                        continue
                                    term_id = str(term_id)
                                    if 'enrichment' in valid_types and term_id in single_term_ids:
                                        safe_name = _safe_plot_stem(term_id)
                                        out_file = str(plot_dir / f"{_safe_plot_stem(db_name)}_{safe_name}_enrichment.{plot_format}")
                                        if plot_gsea_enrichment_r(
                                            tsv_path, term_id, out_file, running_es_path, dpi=plot_dpi,
                                            style=plot_style, palette=plot_palette,
                                        ):
                                            generated_files.append(out_file)
                                if 'enrichment2' in valid_types and top_pathways_for_enrichment2 is not None:
                                    id_col = next(
                                        (c for c in ['pathway', 'Term_ID', 'term_id', 'ID', 'id']
                                         if c in top_pathways_for_enrichment2.columns),
                                        None,
                                    )
                                    nes_col = next(
                                        (c for c in ['NES', 'nes'] if c in top_pathways_for_enrichment2.columns),
                                        None,
                                    )
                                    if id_col and nes_col:
                                        for direction, mask in (
                                            ('up', top_pathways_for_enrichment2[nes_col] > 0),
                                            ('down', top_pathways_for_enrichment2[nes_col] < 0),
                                        ):
                                            term_ids = top_pathways_for_enrichment2.loc[mask, id_col].astype(str).tolist()
                                            if not term_ids:
                                                continue
                                            out_file = str(plot_dir / f"{db_name}_enrichment2_{direction}.{plot_format}")
                                            if plot_gsea_enrichment2_r(
                                                tsv_path, term_ids, out_file, running_es_path, dpi=plot_dpi,
                                                style=plot_style, palette=plot_palette,
                                            ):
                                                generated_files.append(out_file)
                        else:
                            logger.info("No up- or down-regulated pathways met the selection criteria; skipping the plot")

        if not use_r_plots:
            # Python plotting mode
            from allenricher.visualization.gsea_plots import (
                plot_gsea_enrichment,
                plot_gsea_lollipop,
                plot_gsea_multi_enrichment,
                plot_gsea_ridgeplot,
            )

            for db_name, df in results.items():
                if df is None or len(df) == 0:
                    continue
                db_gene_sets = _gene_sets_for_database(
                    db_name, gene_sets, gene_sets_by_database
                )
                name_map = (term_name_maps or {}).get(str(db_name).upper(), {})
                df = _with_plot_term_names(df, name_map, db_name)

                if 'lollipop' in valid_types:
                    out_file = str(plot_dir / f"{db_name}_lollipop.{plot_format}")
                    try:
                        plot_gsea_lollipop(df, output_file=out_file, dpi=plot_dpi,
                                           title=f"{db_name} Gene Set Enrichment Analysis (GSEA)",
                                           style=plot_style, palette=plot_palette,
                                           figsize=plot_figsize)
                        generated_files.append(out_file)
                        logger.info("Generated GSEA lollipop plot: %s", out_file)
                    except Exception as e:
                        logger.error("Failed to generate GSEA lollipop plot for %s: %s", db_name, e)

                if 'barplot' in valid_types:
                    out_file = str(plot_dir / f"{db_name}_barplot.{plot_format}")
                    try:
                        from types import SimpleNamespace
                        plotter = Plotter(str(plot_dir), SimpleNamespace(plot_dpi=plot_dpi))
                        plotter.plot_barplot(
                            df, db_name, Path(out_file).name, top_n=20,
                            style=plot_style, palette=plot_palette, figsize=plot_figsize,
                        )
                        generated_files.append(out_file)
                        logger.info("Generated GSEA bar plot: %s", out_file)
                    except Exception as e:
                        logger.error("Failed to generate GSEA bar plot for %s: %s", db_name, e)

                if 'ridgeplot' in valid_types and ranked_genes and db_gene_sets:
                    out_file = str(plot_dir / f"{db_name}_ridgeplot.{plot_format}")
                    try:
                        plot_gsea_ridgeplot(
                            results_df=df,
                            ranked_genes=ranked_genes,
                            gene_weights=gene_weights or {},
                            gene_sets=db_gene_sets,
                            top_n=10,
                            output_file=out_file,
                            dpi=plot_dpi,
                            figsize=plot_figsize,
                            style=plot_style,
                            palette=plot_palette,
                        )
                        generated_files.append(out_file)
                    except Exception as e:
                        logger.error(f"GSEA ridgeplot failed ({db_name}): {e}")

                # Single-pathway running enrichment-score plots.
                if 'enrichment' in valid_types and ranked_genes and db_gene_sets:
                    top_pathways = _select_gsea_enrichment_terms(
                        df,
                        top_up=gsea_enrichment_top_up,
                        top_down=gsea_enrichment_top_down,
                    )
                    if len(top_pathways) == 0:
                        logger.info("No up- or down-regulated pathways met the selection criteria; skipping the plot")
                    else:
                        term_id_col = next((c for c in ['pathway', 'Term_ID', 'term_id', 'ID', 'id'] if c in df.columns), None)
                        pathway_col = next((c for c in ['Description', 'pathway', 'Term_Name', 'term_name'] if c in df.columns), None)

                        if term_id_col:
                            top_names = set(top_pathways[term_id_col].astype(str).tolist())
                        elif pathway_col:
                            top_names = set(top_pathways[pathway_col].astype(str).tolist())
                        else:
                            top_names = set()

                        for set_name, gene_set in db_gene_sets.items():
                            if str(set_name) not in top_names:
                                continue

                            match = None
                            if term_id_col:
                                match_rows = df[df[term_id_col].astype(str) == str(set_name)]
                                if len(match_rows) > 0:
                                    match = match_rows.iloc[0]
                            if match is None and pathway_col:
                                match_rows = df[df[pathway_col].astype(str) == str(set_name)]
                                if len(match_rows) > 0:
                                    match = match_rows.iloc[0]
                            if match is None:
                                continue

                            es_val = match.get('ES', match.get('enrichmentScore', match.get('es', 0.0)))
                            nes_val = match.get('NES', match.get('nes', 0.0))
                            pval = match.get('pval', match.get('p_value', match.get('NOM p-val', match.get('pvalue', match.get('P_Value', 1.0)))))
                            padj = match.get('padj', match.get('Adjusted_P_Value', match.get('FDR', match.get('p.adjust', match.get('qvalue', None)))))

                            out_file = str(plot_dir / f"{_safe_plot_stem(db_name)}_{_safe_plot_stem(set_name)}_enrichment.{plot_format}")
                            try:
                                plot_gsea_enrichment(
                                    ranked_genes=ranked_genes,
                                    gene_weights=gene_weights or {},
                                    gene_set=gene_set,
                                    es=es_val,
                                    nes=nes_val,
                                    pvalue=pval,
                                    padj=padj,
                                    title=str(match.get('Term_Name', set_name)),
                                    output_file=out_file,
                                    dpi=plot_dpi,
                                    style=plot_style,
                                    palette=plot_palette,
                                    figsize=plot_figsize or (6.60, 4.40),
                                )
                                generated_files.append(out_file)
                                logger.info(f"  GSEA enrichment plot generated: {set_name} (NES={float(nes_val):.2f})")
                            except Exception as e:
                                logger.error(f"GSEA enrichment plot failed ({set_name}): {e}")

                if 'enrichment2' in valid_types and ranked_genes and db_gene_sets:
                    top_pathways = _select_gsea_enrichment_terms(
                        df,
                        top_up=gsea_multi_top_up,
                        top_down=gsea_multi_top_down,
                    )
                    id_col = next(
                        (c for c in ['pathway', 'Term_ID', 'term_id', 'ID', 'id'] if c in top_pathways.columns),
                        None,
                    )
                    nes_col = next((c for c in ['NES', 'nes'] if c in top_pathways.columns), None)
                    if id_col and nes_col:
                        for direction, mask in (
                            ('up', top_pathways[nes_col] > 0),
                            ('down', top_pathways[nes_col] < 0),
                        ):
                            selected_ids = top_pathways.loc[mask, id_col].astype(str).tolist()
                            selected_ids = [term_id for term_id in selected_ids if term_id in db_gene_sets]
                            if not selected_ids:
                                continue
                            out_file = str(plot_dir / f"{db_name}_enrichment2_{direction}.{plot_format}")
                            try:
                                plot_gsea_multi_enrichment(
                                    results_df=df,
                                    selected_ids=selected_ids,
                                    ranked_genes=ranked_genes,
                                    gene_weights=gene_weights or {},
                                    gene_sets=db_gene_sets,
                                    output_file=out_file,
                                    dpi=plot_dpi,
                                    figsize=plot_figsize,
                                    style=plot_style,
                                    palette=plot_palette,
                                )
                                generated_files.append(out_file)
                            except Exception as e:
                                logger.error(
                                    f"GSEA multi-pathway {direction} plot failed ({db_name}): {e}"
                                )

            if 'emapplot' in valid_types:
                logger.warning("emapplot requires the R/aPEAR backend; rerun with --use-r-plots")

    # ssGSEA and GSVA plots
    elif method in ('ssgsea', 'gsva'):
        from allenricher.visualization.gsva_plots import (
            plot_pathway_heatmap,
            plot_group_comparison,
            plot_sample_correlation,
            select_activity_heatmap_scores,
        )
        method_label = 'ssGSEA' if method == 'ssgsea' else 'GSVA'
        activity_frames = []
        for db_name, df in results.items():
            if df is None or len(df) == 0 or not isinstance(df, pd.DataFrame):
                continue
            numeric_cols = df.select_dtypes(include='number').columns
            non_metric_cols = {
                'p_value', 'FDR', 'NOM p-val', 'FDR q-val', 'FWER p-val',
                'pvalue', 'P_Value', 'Adjusted_P_Value', 'p.adjust', 'qvalues',
                'nes', 'es', 'fdr', 'gene_count', 'Gene_Count', 'NES',
                'enrichmentScore', 'setSize',
            }
            sample_cols = [column for column in numeric_cols if column not in non_metric_cols]
            if sample_cols:
                name_col = next(
                    (
                        column for column in
                        ['Description', 'pathway', 'Term_Name', 'Term_ID', df.index.name]
                        if column and column in df.columns
                    ),
                    None,
                )
                scores_df = df.set_index(name_col)[sample_cols] if name_col else df[sample_cols]
                scores_df.index.name = scores_df.index.name or 'pathway'
                activity_frames.append((str(db_name), scores_df))
            elif expr_matrix is not None:
                activity_frames.append((str(db_name), expr_matrix))

        if not activity_frames:
            logger.warning("No pathway-activity matrix was available; skipping ssGSEA/GSVA plots")
            return generated_files

        if 'group_comparison' in valid_types and not groups:
            logger.warning("Group comparison requires --groups and will be skipped")

        multi_database = len(activity_frames) > 1
        for db_name, scores_df in activity_frames:
            prefix = f"{_safe_plot_stem(db_name)}_" if multi_database else ""
            activity_title = (
                f"{db_name} {method_label} Pathway Activity"
                if multi_database else f"{method_label} Pathway Activity"
            )
            annotation_df = None
            if groups:
                sample_to_group = {
                    sample: group_name
                    for group_name, samples in groups.items()
                    for sample in samples
                }
                common_samples = [sample for sample in scores_df.columns if sample in sample_to_group]
                if common_samples:
                    annotation_df = pd.DataFrame(
                        {'Group': [sample_to_group[sample] for sample in common_samples]},
                        index=common_samples,
                    )

            scores_path = plot_dir / f"{prefix}activity_scores.tsv"
            metadata_path = plot_dir / f"{prefix}sample_metadata.tsv"

            if 'heatmap' in valid_types:
                out_file = str(plot_dir / f"{prefix}activity_heatmap.{plot_format}")
                heatmap_scores = select_activity_heatmap_scores(
                    scores_df, annotation_df, activity_heatmap_top_n
                )
                try:
                    generated = False
                    if use_r_plots and annotation_df is not None:
                        from allenricher.visualization.r_plotter import (
                            check_r_environment,
                            plot_activity_heatmap_r,
                        )
                        if check_r_environment():
                            heatmap_scores.to_csv(scores_path, sep='\t', index_label='Pathway')
                            annotation_df.to_csv(metadata_path, sep='\t', index_label='Sample')
                            generated = plot_activity_heatmap_r(
                                str(scores_path), str(metadata_path), out_file,
                                analysis_method=method, top_n=activity_heatmap_top_n,
                                dpi=plot_dpi,
                                style=plot_style, palette=plot_palette,
                            )
                    elif use_r_plots:
                        logger.info("R activity heatmaps require sample-group metadata; using the Python renderer")
                    if not generated:
                        plot_pathway_heatmap(
                            heatmap_scores,
                            annotation_col=annotation_df,
                            title=activity_title,
                            output_file=out_file,
                            dpi=plot_dpi,
                            style=plot_style,
                            palette=plot_palette,
                            figsize=plot_figsize,
                        )
                        generated = True
                    if generated:
                        generated_files.append(out_file)
                        logger.info(f"Active heatmap generated ({db_name}): {out_file}")
                except Exception as e:
                    logger.error("Failed to generate the activity heatmap for %s: %s", db_name, e)

            if 'group_comparison' in valid_types and groups:
                out_file = str(plot_dir / f"{prefix}group_comparison.{plot_format}")
                statistics_file = str(plot_dir / f"{prefix}group_comparison.statistics.tsv")
                try:
                    generated = False
                    if use_r_plots and annotation_df is not None:
                        from allenricher.visualization.r_plotter import (
                            check_r_environment,
                            plot_group_comparison_r,
                        )
                        if check_r_environment():
                            scores_df.to_csv(scores_path, sep='\t', index_label='Pathway')
                            annotation_df.to_csv(metadata_path, sep='\t', index_label='Sample')
                            generated = plot_group_comparison_r(
                                str(scores_path), str(metadata_path), out_file, statistics_file,
                                dpi=plot_dpi,
                                style=plot_style, palette=plot_palette,
                            )
                    if not generated:
                        plot_group_comparison(
                            scores_df,
                            groups=groups,
                            output_file=out_file,
                            statistics_file=statistics_file,
                            dpi=plot_dpi,
                            style=plot_style,
                            palette=plot_palette,
                            figsize=plot_figsize,
                        )
                        generated = True
                    if generated:
                        generated_files.append(out_file)
                        logger.info(
                            f"The inter-group comparison map and statistical table has been generated ({db_name}):"
                            f"{out_file}, {statistics_file}"
                        )
                except Exception as e:
                    logger.error("Failed to generate the group comparison figure for %s: %s", db_name, e)

            if 'correlation' in valid_types:
                if scores_df.shape[1] < 2:
                    logger.warning("Sample correlation for %s requires at least two samples and will be skipped", db_name)
                    continue
                out_file = str(plot_dir / f"{prefix}sample_correlation.{plot_format}")
                try:
                    generated = False
                    if use_r_plots and annotation_df is not None:
                        from allenricher.visualization.r_plotter import (
                            check_r_environment,
                            plot_sample_correlation_r,
                        )
                        if check_r_environment():
                            scores_df.to_csv(scores_path, sep='\t', index_label='Pathway')
                            annotation_df.to_csv(metadata_path, sep='\t', index_label='Sample')
                            generated = plot_sample_correlation_r(
                                str(scores_path), str(metadata_path), out_file, dpi=plot_dpi,
                                style=plot_style, palette=plot_palette,
                            )
                    elif use_r_plots:
                        logger.info("R sample correlation requires sample-group metadata; using the Python renderer")
                    if not generated:
                        plot_sample_correlation(
                            scores_df,
                            annotation_col=annotation_df,
                            output_file=out_file,
                            dpi=plot_dpi,
                            style=plot_style,
                            palette=plot_palette,
                            figsize=plot_figsize,
                        )
                        generated = True
                    if generated:
                        generated_files.append(out_file)
                        logger.info(f"Sample-correlation heatmap generated ({db_name}): {out_file}")
                except Exception as e:
                    logger.error("Failed to generate the sample-correlation figure for %s: %s", db_name, e)

    return generated_files


def _generate_plots(
    method: str,
    results: dict,
    ranked_genes: Optional[list],
    gene_weights: Optional[dict],
    gene_sets: Optional[Dict[str, Set[str]]],
    expr_matrix,
    groups: Optional[Dict[str, List[str]]],
    plot_types: List[str],
    output_dir: str,
    gene_sets_by_database: Optional[Dict[str, Dict[str, Set[str]]]] = None,
    plot_format: str = 'png',
    plot_formats: Optional[List[str]] = None,
    plot_dpi: int = 300,
    plot_style: str = 'nature',
    plot_palette: PaletteLike = None,
    use_r_plots: bool = False,
    emapplot_qvalue: float = 0.05,
    emapplot_min_count: int = 3,
    emapplot_top_n: int = 30,
    gsea_enrichment_top_up: int = 5,
    gsea_enrichment_top_down: int = 5,
    gsea_multi_top_up: int = 3,
    gsea_multi_top_down: int = 3,
    activity_heatmap_top_n: int = 40,
    plot_figsize: Optional[Tuple[float, float]] = None,
    term_name_maps: Optional[Dict[str, Dict[str, str]]] = None,
) -> List[str]:
    """Generate each requested plot in one or more output formats."""
    formats = plot_formats or [plot_format]
    formats = list(dict.fromkeys(str(fmt).lower() for fmt in formats))
    invalid = [fmt for fmt in formats if fmt not in {'png', 'pdf', 'svg'}]
    if invalid:
        raise ValueError(f"Unsupported plot format: {invalid}")

    generated: List[str] = []
    for fmt in formats:
        generated.extend(_generate_plots_one_format(
            method=method,
            results=results,
            ranked_genes=ranked_genes,
            gene_weights=gene_weights,
            gene_sets=gene_sets,
            gene_sets_by_database=gene_sets_by_database,
            expr_matrix=expr_matrix,
            groups=groups,
            plot_types=plot_types,
            output_dir=output_dir,
            plot_format=fmt,
            plot_dpi=plot_dpi,
            plot_style=plot_style,
            plot_palette=plot_palette,
            use_r_plots=use_r_plots,
            emapplot_qvalue=emapplot_qvalue,
            emapplot_min_count=emapplot_min_count,
            emapplot_top_n=emapplot_top_n,
            gsea_enrichment_top_up=gsea_enrichment_top_up,
            gsea_enrichment_top_down=gsea_enrichment_top_down,
            gsea_multi_top_up=gsea_multi_top_up,
            gsea_multi_top_down=gsea_multi_top_down,
            activity_heatmap_top_n=activity_heatmap_top_n,
            plot_figsize=plot_figsize,
            term_name_maps=term_name_maps,
        ))
    return generated


def _positive_int(value: str) -> int:
    parsed = int(value)
    if parsed <= 0:
        raise argparse.ArgumentTypeError("Value must be a positive integer")
    return parsed


def _nonnegative_int(value: str) -> int:
    parsed = int(value)
    if parsed < 0:
        raise argparse.ArgumentTypeError("Value must be a non-negative integer")
    return parsed


def _probability(value: str) -> float:
    parsed = float(value)
    if not 0 <= parsed <= 1:
        raise argparse.ArgumentTypeError("Must be between 0 and 1")
    return parsed


def create_parser() -> argparse.ArgumentParser:
    """Create the top-level command-line parser and all subcommands."""
    # Top-level command and shared version information.
    parser = argparse.ArgumentParser(
        prog='allenricher',
        description='AllEnricher v2.0 - Gene Set Enrichment Analysis Tool',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog='''
Examples:
  # Basic analysis
  allenricher analyze -i genes.txt -s hsa -d GO,KEGG -o results/

  # With AI interpretation
  allenricher analyze -i genes.txt -s hsa --ai openai --ai-key YOUR_KEY

  # Download databases
  allenricher download -d GO,KEGG -s hsa

  # Start API server
  allenricher serve --port 8000
        '''
    )

    # Expose the installed package version.
    parser.add_argument('-v', '--version', action='version', version=f'%(prog)s {__version__}')

    # Create Subcommand Solver Group
    subparsers = parser.add_subparsers(dest='command', help='Available commands')

    # Canonical analysis workflow shared by the CLI, API, and Web application.
    # Configure statistics, plots, reports, and optional AI interpretation.
    analyze_parser = subparsers.add_parser('analyze', help='Run enrichment analysis')
    analyze_parser.add_argument('-i', '--input', required=True, help='Input gene list file')           # Enter the path of the gene list file (required)
    analyze_parser.add_argument('-s', '--species', default='hsa', help='Species code (default: hsa)')  # species code, default human (hsa)
    analyze_parser.add_argument('-d', '--databases', default='GO,KEGG', help='Comma-separated databases')  # Comma-separated list of database names
    analyze_parser.add_argument('-o', '--output', default='./results', help='Output directory')        # Output directory, default to./results
    analyze_parser.add_argument('-b', '--background', help='Background gene list file')
    analyze_parser.add_argument('--background-mode', dest='background_mode',
                                choices=['annotated', 'genome', 'custom'], default='annotated',
                                help='Background gene set mode: annotated (default), genome, custom')
    analyze_parser.add_argument('-m', '--method', default='hypergeometric', choices=['hypergeometric', 'gsea', 'ssgsea', 'gsva'], help='Enrichment method')
    analyze_parser.add_argument('-c', '--correction', default='BH', choices=['BH', 'BY', 'bonferroni', 'holm', 'none'], help='Multiple testing correction')  # Multiple Test Correction Method
    analyze_parser.add_argument('-p', '--pvalue', type=float, default=0.05, help='P-value cutoff')    # P value threshold, default 0.05
    analyze_parser.add_argument('-q', '--qvalue', type=float, default=0.05, help='Q-value cutoff')    # Q Value (after correction P) threshold, default 0.05
    analyze_parser.add_argument('-n', '--min-genes', type=int, default=3, help='Minimum query-gene hits per ORA term')
    analyze_parser.add_argument('-j', '--jobs', type=int, default=1, help='Number of parallel jobs')   # Number of parallel tasks
    analyze_parser.add_argument('--no-plot', action='store_true', help='Skip figure generation')
    analyze_parser.add_argument('--no-report', action='store_true', help='Skip HTML report generation')
    analyze_parser.add_argument(
        '--methods-language',
        choices=['zh', 'en'],
        default='en',
        help='Language for the Materials and Methods writing reference (English only)',
    )
    analyze_parser.add_argument('--only-significant', action='store_true', help='Only output significant terms (filter by p/q cutoff)')  # Only output significant entries (at p)/qThreshold filter, default not enabled, output all entries)
    analyze_parser.add_argument('--ai', choices=['openai', 'claude', 'deepseek', 'glm', 'minimax', 'ollama', 'mock'],
                                help='AI backend for interpretation (override YAML config)')  # AI Read Backend Selection
    analyze_parser.add_argument('--ai-mode', choices=['summary', 'reviewer', 'caption'], default='summary',
                                help='AI interpretation profile (default: summary)')
    analyze_parser.add_argument('--ai-top-n', type=int, default=None, metavar='N',
                                help='AI evidence count override; default ORA=15, GSEA=10 per NES direction, ssGSEA/GSVA=10')
    analyze_parser.add_argument('--ai-key', help='AI API key (override YAML config, optional if set in YAML)')  # AI Service API Key
    analyze_parser.add_argument('--ai-model', help='AI model name (override YAML config)')  # AI model name
    analyze_parser.add_argument('--config', help='Configuration file (YAML/JSON)')                      # Path to external profile
    analyze_parser.add_argument('--database-dir', help='Database directory')                       # Database Directory Path
    analyze_parser.add_argument('--use-version', type=str, default=None,
                                help='Specifies the version of the database used (e. g. v2020515), using the latest version by default')
    analyze_parser.add_argument('-e', '--expression-matrix', default=None, help='Expression matrix (TSV/CSV; genes in rows and samples in columns) for GSEA/ssGSEA/GSVA')
    analyze_parser.add_argument('-r', '--ranked-genes', default=None, help="Ranked gene table for GSEA with 'gene' and numeric 'weight' columns")
    analyze_parser.add_argument('-g', '--gmt', default=None, help='GMT gene-set file (supports .gmt and .gmt.gz)')
    analyze_parser.add_argument('-pt', '--plot-types', default=None, help='Comma-separated figure types. GSEA: enrichment,enrichment2,barplot,lollipop,ridgeplot,emapplot. GSVA/ssGSEA: heatmap,group_comparison,correlation')
    analyze_parser.add_argument('--groups', default=None, help='Sample group definition, format: Group1:sample1,sample2;Group2:sample3,sample4')  # Sample group definition for intergroup comparison
    analyze_parser.add_argument('--plot-format', default=None, choices=['png', 'pdf', 'svg'], help='Plot output format (default: config or png/pdf)')  # Chart Output Format
    analyze_parser.add_argument('--plot-dpi', type=_positive_int, default=None, help='Plot resolution/DPI (default: config or 300)')  # Chart Resolution (DPI)
    analyze_parser.add_argument('--style', default=None, choices=['nature', 'science', 'presentation', 'cell', 'omicshare'], help='Plot theme: nature, science, or presentation; cell/omicshare are legacy aliases')  # Chart-style themes
    analyze_parser.add_argument(
        '--palette', default=None, choices=sorted(PUBLIC_PALETTES),
        help='Legacy color override; automatically applied only to its compatible color role',
    )
    analyze_parser.add_argument(
        '--categorical-palette', default=None,
        choices=sorted(PUBLIC_CATEGORICAL_PALETTES),
        help='Palette for categories, groups and pathway curves',
    )
    analyze_parser.add_argument(
        '--sequential-palette', default=None,
        choices=sorted(PUBLIC_SEQUENTIAL_PALETTES),
        help='Palette for one-direction continuous values such as FDR significance',
    )
    analyze_parser.add_argument(
        '--diverging-palette', default=None,
        choices=sorted(PUBLIC_DIVERGING_PALETTES),
        help='Palette for values centered at zero such as NES, activity and correlation',
    )
    analyze_parser.add_argument('--verbose', action='store_true', help='Enable verbose (DEBUG) logging')  # Enable detailed log output
    analyze_parser.add_argument(
        '--tf-database',
        choices=['trrust', 'chea3', 'animaltfdb', 'htftarget', 'both'],
        help='Include an additional TF ORA; for other methods pass the TF database to -d',
    )
    analyze_parser.add_argument('--tf-library', default=None,
                                help='Comma-separated ChEA3 libraries to include')
    analyze_parser.add_argument('--tf-tissue', default=None,
                                help='Comma-separated hTFtarget tissue contexts to include')
    analyze_parser.add_argument('--tf-regulation', default='all',
                                choices=['all', 'activation', 'repression', 'unknown'],
                                help='TRRUST edge regulation filter (default: all)')
    analyze_parser.add_argument('--tf-min-size', type=_positive_int, default=None,
                                help='Minimum TF target-set size (defaults: ORA 3, GSEA 15, ssGSEA/GSVA 1)')
    analyze_parser.add_argument('--tf-max-size', type=_positive_int, default=None,
                                help='Maximum TF target-set size (defaults: ORA/ssGSEA/GSVA unbounded, GSEA 5000)')
    analyze_parser.add_argument('--tf-combine', default='none',
                                choices=['none', 'meanrank', 'toprank'],
                                help='Optional rank consensus; source-level results are always retained')
    analyze_parser.add_argument('--tf-only', action='store_true',
                                help='Only perform TF enrichment, skip standard databases (GO/KEGG/Reactome etc.)')
    analyze_parser.add_argument('--use-r-plots', action='store_true', help='Use R scripts for GSEA plotting (requires R environment)')
    analyze_parser.add_argument('--emapplot-qvalue', type=_probability, default=0.05,
                                help='R emapplot p.adjust/FDR cutoff (default: 0.05)')
    analyze_parser.add_argument('--emapplot-min-count', type=_positive_int, default=3,
                                help='R emapplot minimum Count/Gene_Count filter (default: 3)')
    analyze_parser.add_argument('--emapplot-top-n', type=_positive_int, default=30,
                                help='R emapplot top pathway count after filtering (default: 30)')
    analyze_parser.add_argument('--gsea-enrichment-top-up', type=_nonnegative_int, default=5,
                                help='Number of positive NES pathways for single-pathway GSEA enrichment plots (default: 5)')
    analyze_parser.add_argument('--gsea-enrichment-top-down', type=_nonnegative_int, default=5,
                                help='Number of negative NES pathways for single-pathway GSEA enrichment plots (default: 5)')
    analyze_parser.add_argument('--gsea-multi-top-up', type=_nonnegative_int, default=3,
                                help='Number of positive NES pathways for the multi-pathway GSEA plot (default: 3)')
    analyze_parser.add_argument('--gsea-multi-top-down', type=_nonnegative_int, default=3,
                                help='Number of negative NES pathways for the multi-pathway GSEA plot (default: 3)')

    # Download shared source data and database-specific files.
    download_parser = subparsers.add_parser('download', help='Download databases')
    download_parser.add_argument('-d', '--databases', required=True, help='Comma-separated databases to download')  # Database name to download (required)
    download_parser.add_argument('-s', '--species', default='hsa', help='Species code')                # Species Code
    download_parser.add_argument('--database-dir', default='./database', help='Database directory')     # Database Storage Directory
    download_parser.add_argument('--workers', type=int, default=4, help='Multi-thread download workers (default: 4)')  # Multiline downloads
    download_parser.add_argument('--no-multi-thread', action='store_true', help='Disable multi-thread download')       # Disable Multiline
    download_parser.add_argument('--no-verify', action='store_true', help='Skip post-download integrity checks')
    download_parser.add_argument('--force', action='store_true', help='Force re-download even if local is the latest version')
    download_parser.add_argument('--trrust', action='store_true', help='Download the TRRUST TF-target database')
    download_parser.add_argument('--chea3', action='store_true', help='Download ChEA3 TF-target database')  # Download ChEA3 transcription factor-target gene database

    # Build analysis-ready database artifacts for one species.
    build_parser = subparsers.add_parser('build', help='Build species database')
    build_parser.add_argument('-s', '--species', required=True, help='Species code')                    # Species Code (necessary)
    build_parser.add_argument('-t', '--taxonomy', required=True, type=int, help='Taxonomy ID')          # NCBI Classification ID (Necessary)
    build_parser.add_argument('-d', '--databases', default='GO,KEGG,Reactome', help='Comma-separated databases to build')  # List of databases to be constructed
    build_parser.add_argument('--database-dir', default='./database', help='Database directory')        # Database Storage Directory
    build_parser.add_argument('--gene-info', help='Path to NCBI gene_info.gz file')                    # NCBI Gene_info.gz file path (GO and Reactome build needs)

    # Optional custom annotation input.
    build_parser.add_argument('--go-annot',
        help='Path to GO annotation file (TSV: gene<TAB>go_id<TAB>go_name[<TAB>hierarchy])')
    build_parser.add_argument('--kegg-annot',
        help='Path to KEGG annotation file (TSV: gene<TAB>pathway_id<TAB>pathway_name[<TAB>hierarchy])')
    build_parser.add_argument('--custom-annot',
        help='Path to custom annotation file (TSV format with hierarchy support)')
    build_parser.add_argument('--custom-db-name', default='CUSTOM',
        help='Database name for custom annotation (default: CUSTOM)')
    build_parser.add_argument('--annot-format',
        choices=['three_column', 'four_column', 'two_column', 'auto'],
        default='auto', help='Annotation file format (default: auto-detect)')
    build_parser.add_argument('--hierarchy-sep', default='|',
        help='Hierarchy level separator (default: |)')
    build_parser.add_argument('--latin-name', type=str, default='',
                              help='Latin species (underline format, for example, Bos_taurus)')

    # Start the local REST API and Web application.
    serve_parser = subparsers.add_parser('serve', help='Start API server')
    serve_parser.add_argument('--host', default='127.0.0.1', help='Server host')                        # Only listen to the current system by default
    serve_parser.add_argument('--port', type=int, default=8000, help='Server port')                     # Server listener port, default 8000
    serve_parser.add_argument('--reload', action='store_true', help='Enable development auto-reload')
    serve_parser.add_argument('--config', help='Server-side YAML/JSON configuration, including optional AI settings')

    # List of supported species or available database resources
    list_parser = subparsers.add_parser('list', help='List available resources')
    list_parser.add_argument('resource', choices=['species', 'databases'], help='Resource to list')      # Resource type to view: species (species) or datas (database)

    # Generate default YAML profile, which the user can modify on
    config_parser = subparsers.add_parser('config', help='Generate configuration file')
    config_parser.add_argument('-o', '--output', default='allenricher.yaml', help='Output configuration file')

    # Check if remote data sources are updated
    check_update_parser = subparsers.add_parser('check-update', help='Check if remote data sources are updated')
    check_update_parser.add_argument('--database-dir', default=None, help='Database Directory Path')
    check_update_parser.add_argument('--json', action='store_true', help='Output results in JSON format')

    # Clear old version of database files
    cleanup_parser = subparsers.add_parser('cleanup', help='Clear old version of database files')
    cleanup_parser.add_argument('--keep', type=int, default=2, help='Number of most recent versions retained (default: 2)')
    cleanup_parser.add_argument('--dry-run', action='store_true', help='Preview only, not physically deleted')
    cleanup_parser.add_argument('--database-dir', default=None, help='Database Directory Path')

    # View locally installed database version
    list_versions_parser = subparsers.add_parser('list-versions', help='View locally installed database version')
    list_versions_parser.add_argument('--database-dir', default=None, help='Database Directory Path')
    list_versions_parser.add_argument('--json', action='store_true', help='Output in JSON format')
    list_versions_parser.add_argument('--lineage', action='store_true', help='Show database build provenance')

    # List species recorded in the unified species registry.
    list_species_parser = subparsers.add_parser('list-species', help='List supported species from registry')
    list_species_parser.add_argument('--go', action='store_true', default=False, help='Filter by GO support')
    list_species_parser.add_argument('--kegg', action='store_true', default=False, help='Filter by KEGG support')
    list_species_parser.add_argument('--reactome', action='store_true', default=False, help='Filter by Reactome support')
    list_species_parser.add_argument('--do', action='store_true', default=False, help='Filter by DO support')
    list_species_parser.add_argument('--disgenet', action='store_true', default=False, help='Filter by DisGeNET support')
    list_species_parser.add_argument('--wikipathways', action='store_true', default=False, help='Filter by WikiPathways support')
    list_species_parser.add_argument('--trrust', action='store_true', default=False, help='Filter by TRRUST support')
    list_species_parser.add_argument('--chea3', action='store_true', default=False, help='Filter by ChEA3 support')
    list_species_parser.add_argument('--animaltfdb', action='store_true', default=False, help='Filter by AnimalTFDB support')
    list_species_parser.add_argument('--htftarget', action='store_true', default=False, help='Filter by hTFtarget support')
    list_species_parser.add_argument('--format', choices=['table', 'tsv', 'json'], default='table', help='Output format (default: table)')
    list_species_parser.add_argument('--summary', action='store_true', default=False, help='Show summary statistics')

    # Search for details of specific species
    query_species_parser = subparsers.add_parser('query-species', help='Query species detail from registry')
    query_species_parser.add_argument('--taxid', type=int, default=None, help='Query by NCBI Taxonomy ID')
    query_species_parser.add_argument('--name', type=str, default=None, help='Query by Latin name')
    query_species_parser.add_argument('--kegg', type=str, default=None, help='Query by KEGG organism code')

    # Dedicated TF-target enrichment entry point retained for compatibility.
    tf_enrich_parser = subparsers.add_parser('tf-enrich', help='Transcription factor enrichment analysis')
    tf_enrich_parser.add_argument(
        '-i', '--input', required=True,
        help="ORA: one gene per line; GSEA: TSV with 'gene' and 'weight' columns",
    )
    tf_enrich_parser.add_argument('-s', '--species', default='hsa', help='Species code (default: hsa)')
    tf_enrich_parser.add_argument(
        '-d', '--database', default='trrust',
        choices=['trrust', 'chea3', 'animaltfdb', 'htftarget'],
        help='TF database (default: trrust)',
    )
    tf_enrich_parser.add_argument('-o', '--output', default='./results', help='Output directory')
    tf_enrich_parser.add_argument(
        '-b', '--background', default=None,
        help='ORA background: one measurable gene per line (not valid for GSEA)',
    )
    tf_enrich_parser.add_argument('--report', action='store_true', help='Generate an HTML report')
    tf_enrich_parser.add_argument('--top-n', type=int, default=20, help='Display the top N TFs (default: 20)')
    tf_enrich_parser.add_argument('--database-dir', default=None, help='Database directory')
    tf_enrich_parser.add_argument('--method', default='ora', choices=['ora', 'gsea'], help='Enrichment method (default: ora)')
    tf_enrich_parser.add_argument('--tf-library', default=None,
                                  help='Comma-separated ChEA3 libraries to include')
    tf_enrich_parser.add_argument('--tf-tissue', default=None,
                                  help='Comma-separated hTFtarget tissue contexts to include')
    tf_enrich_parser.add_argument('--tf-regulation', default='all',
                                  choices=['all', 'activation', 'repression', 'unknown'],
                                  help='TRRUST edge regulation filter (default: all)')
    tf_enrich_parser.add_argument('--tf-min-size', type=_positive_int, default=None,
                                  help='Minimum target-set size (ORA default 3; GSEA default 15)')
    tf_enrich_parser.add_argument('--tf-max-size', type=_positive_int, default=None,
                                  help='Maximum target-set size (ORA unbounded; GSEA default 5000)')
    tf_enrich_parser.add_argument('--tf-combine', default='none',
                                  choices=['none', 'meanrank', 'toprank'],
                                  help='Optional ChEA3 library rank consensus')
    tf_enrich_parser.add_argument('--online', action='store_true',
                                  help='Use ChEA3 API for online analysis (requires internet)')  # Online analytical mode

    return parser


def _resolve_db_dir(args) -> str:
    """Resolve the database directory shared by CLI subcommands."""
    if hasattr(args, 'database_dir') and args.database_dir:
        return args.database_dir
    return "./database"


def _cli_option_was_provided(args, *option_names: str, fallback: bool = False) -> bool:
    """Return whether an option was explicitly present on the command line."""
    provided = getattr(args, "_provided_options", None)
    if provided is None:
        return fallback
    return any(option in provided for option in option_names)


def cmd_analyze(args) -> int:
    """Run the complete enrichment-analysis workflow.

    The command loads configuration and inputs, resolves local databases and
    the analysis background, executes the selected method, and writes result
    tables, plots, and the optional HTML report.

    Args:
        args: Parsed arguments for the ``analyze`` subcommand.

    Returns:
        ``0`` on success and ``1`` when the analysis cannot be completed.
    """
    # Enable detailed logs only when explicitly requested.
    if args.verbose:
        logging.getLogger('allenricher').setLevel(logging.DEBUG)
        logger.debug("Verbose logging enabled at DEBUG level")

    try:
        logger.info(f"AllEnricher v{__version__} - Starting analysis")

        # Load a configuration file when supplied; explicit CLI options override it.
        if args.config:
            config = Config.from_file(args.config)
            if args.input:
                config.input_file = args.input
            if _cli_option_was_provided(
                args, '-s', '--species', fallback=args.species != 'hsa'
            ):
                config.species = args.species
            if _cli_option_was_provided(
                args, '-d', '--databases', fallback=args.databases != 'GO,KEGG'
            ):
                config.databases = args.databases.split(',')
            if _cli_option_was_provided(
                args, '-m', '--method', fallback=args.method != 'hypergeometric'
            ):
                config.method = args.method
            if _cli_option_was_provided(
                args, '-c', '--correction', fallback=args.correction != 'BH'
            ):
                config.correction = args.correction
            if _cli_option_was_provided(
                args, '-p', '--pvalue', fallback=args.pvalue != 0.05
            ):
                config.pvalue_cutoff = args.pvalue
            if _cli_option_was_provided(
                args, '-q', '--qvalue', fallback=args.qvalue != 0.05
            ):
                config.qvalue_cutoff = args.qvalue
            if _cli_option_was_provided(
                args, '-n', '--min-genes', fallback=args.min_genes != 3
            ):
                config.min_genes = args.min_genes
            if _cli_option_was_provided(
                args, '-j', '--jobs', fallback=args.jobs != 1
            ):
                config.n_jobs = args.jobs
            if _cli_option_was_provided(
                args, '-o', '--output', fallback=args.output != './results'
            ):
                config.output_dir = args.output
            if args.background:
                config.background_file = args.background
            # CLI mark overoutput_all (False when default True, --only-significant)
            if args.only_significant:
                config.output_all = False
            # Chart Style Parameters
            if hasattr(args, 'style') and args.style:
                config.plot_style = args.style
            if hasattr(args, 'palette') and args.palette:
                config.plot_palette = args.palette
            if getattr(args, 'categorical_palette', None):
                config.categorical_palette = args.categorical_palette
            if getattr(args, 'sequential_palette', None):
                config.sequential_palette = args.sequential_palette
            if getattr(args, 'diverging_palette', None):
                config.diverging_palette = args.diverging_palette
            if getattr(args, 'plot_format', None):
                config.plot_formats = [args.plot_format]
            if getattr(args, 'plot_dpi', None) is not None:
                config.plot_dpi = args.plot_dpi
        else:
            config = Config(
                input_file=args.input,
                output_dir=args.output,
                species=args.species,
                databases=args.databases.split(','),  # Split Comma-separated strings into a list
                method=args.method,
                correction=args.correction,
                pvalue_cutoff=args.pvalue,
                qvalue_cutoff=args.qvalue,
                min_genes=args.min_genes,
                n_jobs=args.jobs,
                background_file=args.background,
                output_all=not args.only_significant,  # Default output of all entries (same as v1)
                plot_formats=[args.plot_format] if getattr(args, 'plot_format', None) else ["pdf", "png"],
                plot_dpi=getattr(args, 'plot_dpi', None) or 300,
                plot_style=getattr(args, 'style', None) or 'nature',
                plot_palette=getattr(args, 'palette', None),
                categorical_palette=getattr(args, 'categorical_palette', None),
                sequential_palette=getattr(args, 'sequential_palette', None),
                diverging_palette=getattr(args, 'diverging_palette', None),
            )
        # TF term-size defaults depend on the analysis method and are resolved
        # only when a TF database is analyzed.
        config.tf_min_size = getattr(args, 'tf_min_size', None)
        config.tf_max_size = getattr(args, 'tf_max_size', None)

        # ----Step 2: Validation configuration----
        # Validate paths, method settings, and species identifiers before analysis.
        errors = config.validate()
        if errors:
            for error in errors:
                logger.error(f"Configuration error: {error}")
            return 1

        palette_selection = resolve_palette_selection(
            legacy_palette=config.plot_palette,
            categorical_palette=config.categorical_palette,
            sequential_palette=config.sequential_palette,
            diverging_palette=config.diverging_palette,
        )

        # ----Step 3: Create an output directory----
        # Autocreate (including parent) if the output directory does not exist
        output_dir = Path(config.output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

        # Step 4: load the query gene list.
        # Prefer the configured input file, then fall back to the CLI value.
        input_path = config.input_file or args.input
        if not input_path:
            logger.error("No query gene file was provided. Use -i/--input or set input_file in the configuration.")
            return 1
        logger.info(f"Loading gene list from {input_path}")
        analyzer = EnrichmentAnalyzer(config)
        gene_set = analyzer.load_gene_list(input_path)

        # Step 4.5: load the ranked list or expression matrix required by the method.
        expression_matrix = None
        ranked_gene_list = None
        if hasattr(args, 'expression_matrix') and args.expression_matrix:
            logger.info(f"Loading expression matrix from {args.expression_matrix}")
            import pandas as pd
            expr_path = Path(args.expression_matrix)
            if expr_path.suffix.lower() in ('.csv',):
                expression_matrix = pd.read_csv(expr_path, index_col=0)
            else:
                # Read non-CSV matrices as tab-delimited files.
                expression_matrix = pd.read_csv(expr_path, sep='\t', index_col=0)
            logger.info(f"Expression matrix loaded: {expression_matrix.shape[0]} genes x {expression_matrix.shape[1]} samples")
        if hasattr(args, 'ranked_genes') and args.ranked_genes:
            logger.info(f"Loading ranked gene list from {args.ranked_genes}")
            ranked_gene_list = analyzer.load_ranked_gene_list(args.ranked_genes)
            logger.info(f"Ranked gene list loaded: {len(ranked_gene_list)} genes")
        if config.method == 'gsea' and not ranked_gene_list:
            logger.error("GSEA requires --ranked-genes with 'gene' and 'weight' columns")
            return 1

        # Load GMT gene sets used by GSEA, ssGSEA, and GSVA.
        gene_sets = None
        gmt_term_data = None
        if hasattr(args, 'gmt') and args.gmt:
            gmt_term_data = _parse_gmt_term_data(args.gmt)
            gene_sets = {
                term_id: set(term_info['genes'])
                for term_id, term_info in gmt_term_data.items()
            }

        # Step 4.7: parse visualization parameters.
        plot_types_list = None
        if hasattr(args, 'plot_types') and args.plot_types:
            plot_types_list = [pt.strip() for pt in args.plot_types.split(',') if pt.strip()]
            logger.info("Requested plot types: %s", plot_types_list)

        # Parse the sample-group specification.
        groups_dict = None
        if hasattr(args, 'groups') and args.groups:
            groups_dict = _parse_groups(args.groups)
            if expression_matrix is not None:
                matrix_samples = set(map(str, expression_matrix.columns))
                missing_samples = sorted({
                    sample for samples in groups_dict.values() for sample in samples
                    if sample not in matrix_samples
                })
                if missing_samples:
                    raise ValueError(
                        "--groupsInclude samples that do not exist in the expression matrix:" + ", ".join(missing_samples)
                    )

        # Fetch Chart Output Format and DPI
        plot_formats = list(dict.fromkeys(str(fmt).lower() for fmt in config.plot_formats))
        plot_format = plot_formats[0]
        plot_dpi = config.plot_dpi
        plot_figsize = (
            (float(config.plot_width), float(config.plot_height))
            if config.plot_width is not None and config.plot_height is not None
            else None
        )
        if args.use_r_plots:
            if plot_figsize is not None:
                logger.warning(
                    "plot_width and plot_height apply to Python figures only; R figures use their publication layouts"
                )

        # Load the selected species databases.
        # Resolve every configured database through the standard loader.
        # Legacy TF-only mode bypasses the standard database list.
        tf_only_mode = getattr(args, 'tf_only', False)
        tf_database = getattr(args, 'tf_database', None)
        db_dir = args.database_dir if args.database_dir else config.database_dir

        if tf_database and config.method != 'hypergeometric':
            logger.error(
                "--tf-database and --tf-only are legacy TF ORA options. "
                "For GSEA, ssGSEA, or GSVA, select the TF database with --databases."
            )
            return 1

        if tf_database:
            for database_name in ('trrust', 'chea3') if tf_database == 'both' else (tf_database,):
                validate_tf_database_species(database_name, config.species)

        if tf_only_mode:
            logger.info("TF-only mode: skipping standard database loading")
            results = {}
            db_manager = DatabaseManager(db_dir, config.species)
        else:
            logger.info(f"Loading databases: {config.databases} from {db_dir}")
            db_manager = DatabaseManager(db_dir, config.species)
            # Version locked: CLI--use-versionOver YAML config.use_version
            use_ver = getattr(args, 'use_version', None) or getattr(config, 'use_version', None)
            if use_ver:
                logger.info(f"Use specified database version: {use_ver}")
            db_manager.load_databases(config.databases, version=use_ver)

        # Resolve the ORA background from explicit input or the selected mode.
        background_mode = getattr(args, 'background_mode', 'annotated')

        if config.background_file:
            # An explicit background file takes precedence over background-mode.
            logger.info(f"Loading background genes from {config.background_file}")
            background_set = analyzer.load_gene_list(config.background_file)
        elif background_mode == 'custom':
            logger.error("Background mode 'custom' requires --background")
            return 1
        elif tf_only_mode and background_mode == 'genome':
            background_set = db_manager.get_genome_genes(species_code=config.species)
            if not background_set:
                logger.error("TF-only genome background is unavailable for this species")
                return 1
        elif tf_only_mode:
            # TF ORA derives its annotated background from the selected TF library.
            background_set = set()
        elif background_mode == 'annotated':
            logger.info("Using annotated genes as background (background_mode='annotated')")
            background_set = db_manager.get_background_genes()
        elif background_mode == 'genome':
            # Genome background is derived from NCBI gene_info for the selected TaxID.
            logger.info("Using genome genes as background (background_mode='genome')")
            try:
                # Directly using species code to get the whole genome.
                background_set = db_manager.get_genome_genes(species_code=config.species)
                if background_set:
                    logger.info(f"Loaded {len(background_set)} genome genes from gene_info.gz")
                else:
                    logger.warning("No genome genes found in gene_info.gz, falling back to annotated genes")
                    background_set = db_manager.get_background_genes()
            except Exception as e:
                logger.warning(f"Failed to load genome genes: {e}, falling back to annotated genes")
                background_set = db_manager.get_background_genes()
        else:
            # Default retreat
            logger.info("Using all database genes as background")
            background_set = db_manager.get_background_genes()

        # Execute the selected analysis method.
        # Analysis of the abundance of each database running, supporting parallel calculations (when n_jobs > 1)
        term_name_maps: Dict[str, Dict[str, str]] = {}
        gsea_gene_sets_by_database: Optional[Dict[str, Dict[str, Set[str]]]] = None
        if not tf_only_mode:
            logger.info("Running enrichment analysis...")
            database_data = db_manager.get_all_term_data()
            term_name_maps = {
                str(database).upper(): {
                    str(term_id): str(term_info.get("name") or term_id)
                    for term_id, term_info in terms.items()
                }
                for database, terms in database_data.items()
            }
            if config.method == 'gsea':
                if gene_sets is None:
                    gsea_gene_sets_by_database = {
                        str(database): {
                            str(term_id): set(term_info.get("genes", []))
                            for term_id, term_info in terms.items()
                        }
                        for database, terms in database_data.items()
                    }
                    gene_sets = {
                        term_id: genes
                        for sets in gsea_gene_sets_by_database.values()
                        for term_id, genes in sets.items()
                    }
                else:
                    gsea_gene_sets_by_database = {str(config.databases[0]): gene_sets}
            if config.method in ('ssgsea', 'gsva'):
                if expression_matrix is None:
                    logger.error(f"{config.method} requires --expression-matrix")
                    return 1
                if gene_sets:
                    activity_gene_sets = {config.databases[0]: gene_sets}
                    activity_term_data = {config.databases[0]: gmt_term_data or {}}
                else:
                    activity_gene_sets = {}
                    activity_term_data = database_data
                    for database, terms in database_data.items():
                        activity_gene_sets[database] = {
                            str(term_id): set(term_info.get("genes", []))
                            for term_id, term_info in terms.items()
                        }
                results = {
                    database: add_result_term_metadata(
                        analyzer.analyze_activity_database(expression_matrix, sets, database),
                        activity_term_data.get(database, {}),
                    )
                    for database, sets in activity_gene_sets.items()
                }
                results = {database: frame for database, frame in results.items() if not frame.empty}
            else:
                results = analyzer.run_analysis(
                    gene_set, background_set, database_data,
                    parallel=config.n_jobs > 1,
                    ranked_gene_list=ranked_gene_list
                )

        # Run the legacy additional TF analysis when explicitly requested.
        tf_results = None
        if tf_database:
            logger.info(f"Running TF enrichment analysis with database: {tf_database}")
            try:
                tf_background = (
                    background_set
                    if config.background_file or background_mode == 'genome'
                    else None
                )
                tf_results = _run_tf_analysis(
                    args,
                    list(gene_set),
                    config.species,
                    background_genes=tf_background,
                )
                if tf_results is not None and not tf_results.empty:
                    logger.info(f"TF enrichment analysis completed: {len(tf_results)} TFs found")
                else:
                    logger.warning("TF enrichment analysis found no significant TFs")
            except Exception as e:
                logger.error(f"TF enrichment analysis failed: {e}", exc_info=True)
                if tf_only_mode:
                    return 1

        # An empty result is a valid analysis outcome, not a command failure.
        has_standard_results = results and len(results) > 0
        has_tf_results = tf_results is not None and not tf_results.empty

        if not has_standard_results and not has_tf_results:
            logger.warning("=" * 60)
            logger.warning("No terms passed the recorded analysis filters.")
            logger.warning("Check gene identifier compatibility, input overlap, significance thresholds, and the ORA background where applicable.")
            logger.warning("=" * 60)
            return 0  # You're out of the ordinary. No mistakes.

        # ----Step 8: Save the results of the analysis----
        # Save complete result tables before applying display-only significance filters.
        logger.info("Saving results...")

        # Construct metadata dictionary, record version and analysis information
        from datetime import datetime, timezone
        recorded_databases = _recorded_analysis_databases(
            config.databases,
            tf_database,
            tf_only_mode,
            has_tf_results,
        )
        metadata = {
            "allenricher_version": __version__,
            "analysis_date": datetime.now(timezone.utc).isoformat(),
            "database_version": db_manager.active_version or "unknown",
            "database_versions": db_manager.database_versions,
            "species": config.species,
            "databases": recorded_databases,
            "methods_language": getattr(args, 'methods_language', 'en'),
            "parameters": _recorded_analysis_parameters(
                config,
                background_mode,
                recorded_databases,
            ),
        }
        metadata.update(_analysis_species_metadata(config.species, db_dir))
        build_meta = db_manager.get_build_metadata()
        source_versions = dict((build_meta or {}).get("source_versions", {}))
        for database in recorded_databases:
            catalog_version = database_catalog_entry(database).get("source_version")
            if catalog_version:
                source_versions[database] = catalog_version
        if source_versions:
            metadata["source_versions"] = source_versions
        if build_meta:
            metadata["built_at"] = build_meta.get("built_at", "")
        metadata["analysis_method"] = config.method
        with open(output_dir / "analysis_metadata.json", "w", encoding="utf-8") as handle:
            json.dump(metadata, handle, ensure_ascii=False, indent=2)
            handle.write("\n")

        # Save one canonical result table per database.
        if not tf_only_mode:
            if config.method in ('ssgsea', 'gsva'):
                for database, activity_scores in results.items():
                    output_file = output_dir / f"{database}_enrichment.tsv"
                    activity_scores.to_csv(
                        output_file,
                        sep='\t',
                        index=True,
                        index_label='Term_ID',
                        lineterminator='\n',
                    )
                    logger.info(f"Saved {database} activity scores to {output_file}")
            else:
                analyzer.save_results(str(output_dir), metadata=metadata)

        # Save legacy additional TF results separately when present.
        if tf_results is not None and not tf_results.empty:
            tf_output_file = output_dir / "TF_enrichment_results.csv"
            tf_results.to_csv(tf_output_file, index=False)
            logger.info(f"TF enrichment results saved to: {tf_output_file}")
            if getattr(args, 'tf_combine', 'none') != 'none':
                consensus = _tf_rank_consensus(
                    tf_results, getattr(args, 'tf_combine', 'none')
                )
                consensus.to_csv(output_dir / "TF_enrichment_consensus.csv", index=False)

        # ----Step 8.5: Screening of significant results for mapping and reporting----
        # Apply significance thresholds only to downstream figures and report views.
        significant_results = {}
        if config.method in ('ssgsea', 'gsva'):
            significant_results = results.copy()
        else:
            for db_name, df in results.items():
                if len(df) == 0:
                    continue
                # Dynamicly detect listing (compatible OLA and GSEA)
                pval_col = next((column for column in ('pval', 'p_value', 'NOM p-val', 'pvalue', 'P_Value') if column in df.columns), None)
                adj_pval_col = next((column for column in ('padj', 'FDR', 'FDR q-val', 'p.adjust', 'Adjusted_P_Value') if column in df.columns), None)
                if not pval_col or not adj_pval_col:
                    logger.warning("%s has no recognized P-value or adjusted-P-value column; skipping significance filtering", db_name)
                    continue
                filtered = df[
                    (df[pval_col] <= config.pvalue_cutoff) &
                    (df[adj_pval_col] <= config.qvalue_cutoff)
                ].copy()
                if len(filtered) > 0:
                    significant_results[db_name] = filtered

        # Add TF results to significant outcomes (if present and significant)
        if tf_results is not None and not tf_results.empty:
            tf_pval_col = 'Pvalue' if 'Pvalue' in tf_results.columns else 'pvalue'
            tf_adj_col = 'FDR' if 'FDR' in tf_results.columns else 'p.adjust'
            if tf_pval_col in tf_results.columns and tf_adj_col in tf_results.columns:
                tf_significant = tf_results[
                    (tf_results[tf_pval_col] <= config.pvalue_cutoff) &
                    (tf_results[tf_adj_col] <= config.qvalue_cutoff)
                ].copy()
                if len(tf_significant) > 0:
                    significant_results['TF'] = tf_significant
                    logger.info(f"TF significant results: {len(tf_significant)} TFs")

        sig_total = sum(len(df) for df in significant_results.values())
        all_total = sum(len(df) for df in results.values())
        if not tf_only_mode and config.method not in ('ssgsea', 'gsva'):
            logger.info(f"Significant results: {sig_total}/{all_total} terms (p<={config.pvalue_cutoff}, q<={config.qvalue_cutoff})")
        elif config.method in ('ssgsea', 'gsva'):
            logger.info(f"Activity score matrix: {all_total} pathways")

        # Generate method-appropriate figures.
        # GSEA, ssGSEA, and GSVA use their method-specific figure sets.
        if config.method not in _METHOD_PLOT_TYPES:
            if not args.no_plot and significant_results:
                logger.info("Generating ORA-style plots (significant results only)...")
                plotter = Plotter(str(output_dir / "plots"), config)
                for db_name, df in significant_results.items():
                    if len(df) > 0:
                        hierarchy_map = {
                            str(term_id): str(term_info["hierarchy"])
                            for term_id, term_info in database_data.get(db_name, {}).items()
                            if term_info.get("hierarchy")
                        }
                        plotter.plot_all(df, db_name, top_n=config.top_terms,
                                         style=config.plot_style, palette=palette_selection,
                                         hierarchy_map=hierarchy_map or None)

        # Step 9.5: generate method-specific GSEA, ssGSEA, or GSVA figures.
        # Apply method-specific plot defaults unless the user supplied --plot-types.
        gsea_plot_files = []
        _effective_plot_types = plot_types_list
        if not _effective_plot_types and config.method in _METHOD_PLOT_TYPES:
            _effective_plot_types = sorted(_DEFAULT_METHOD_PLOT_TYPES[config.method])
        if (_effective_plot_types and config.method in _METHOD_PLOT_TYPES
                and not args.no_plot and results):
            logger.info(f"Generating {config.method} specific plots...")
            # Build the gene-right dictionary (GSEA curve needs)
            gene_weights = None
            if ranked_gene_list and hasattr(ranked_gene_list, 'items'):
                gene_weights = dict(ranked_gene_list)
            elif ranked_gene_list and isinstance(ranked_gene_list, list):
                # ranked_gene_list is [(gene, weight), _] list
                gene_weights = {g: w for g, w in ranked_gene_list}

            gsea_plot_files = _generate_plots(
                method=config.method,
                results=results,
                ranked_genes=list(ranked_gene_list) if ranked_gene_list else None,
                gene_weights=gene_weights,
                gene_sets=gene_sets,
                gene_sets_by_database=gsea_gene_sets_by_database,
                expr_matrix=expression_matrix,
                groups=groups_dict,
                plot_types=_effective_plot_types,
                output_dir=str(output_dir),
                plot_format=plot_format,
                plot_formats=plot_formats,
                plot_dpi=plot_dpi,
                plot_style=config.plot_style,
                plot_palette=palette_selection,
                use_r_plots=args.use_r_plots,
                emapplot_qvalue=getattr(args, 'emapplot_qvalue', 0.05),
                emapplot_min_count=getattr(args, 'emapplot_min_count', 3),
                emapplot_top_n=getattr(args, 'emapplot_top_n', 30),
                gsea_enrichment_top_up=getattr(args, 'gsea_enrichment_top_up', 5),
                gsea_enrichment_top_down=getattr(args, 'gsea_enrichment_top_down', 5),
                gsea_multi_top_up=getattr(args, 'gsea_multi_top_up', 3),
                gsea_multi_top_down=getattr(args, 'gsea_multi_top_down', 3),
                activity_heatmap_top_n=config.activity_heatmap_top_n,
                plot_figsize=plot_figsize,
                term_name_maps=term_name_maps,
            )
            if gsea_plot_files:
                logger.info("Generated %d method-specific figures", len(gsea_plot_files))

        # Step 10: generate the optional AI interpretation.
        ai_interpretation = None

        # Determines whether AI interpretation is enabled: command line--ai > YAML ai_interpretation
        ai_enabled = args.ai or (config.ai_interpretation and config.ai_backend)

        ai_interpretation_error = None
        if ai_enabled:
            # Determines the backend used: Command Line--ai > YAML ai_backend
            ai_backend = args.ai or config.ai_backend

            logger.info(f"Generating AI interpretation using {ai_backend}...")
            try:
                # If the command line provides--ai-keyor--ai-model, use traditional method (command line parameter first)
                if args.ai_key or (args.ai_model and not config.ai_backends):
                    interpreter = create_interpreter(
                        backend=ai_backend,
                        api_key=args.ai_key,
                        model=args.ai_model
                    )
                else:
                    # Create an interpreter from a Config object (support YAML ai_backends configuration)
                    interpreter = create_interpreter_from_config(config, backend=ai_backend)

                ai_context = {
                    "species": metadata.get("species_name") or metadata.get("species"),
                    "taxid": metadata.get("taxid"),
                    "databases": list(significant_results),
                    "groups": {name: len(samples) for name, samples in (groups_dict or {}).items()},
                }
                if config.method not in ("ssgsea", "gsva"):
                    ai_context["significance_thresholds"] = {
                        "p_value": config.pvalue_cutoff,
                        "adjusted_p_value": config.qvalue_cutoff,
                    }
                ai_interpretation = interpreter.interpret_structured_results(
                    significant_results,
                    method=config.method,
                    mode=getattr(args, 'ai_mode', 'summary'),
                    groups=groups_dict,
                    context="Recorded analysis context: " + json.dumps(ai_context, ensure_ascii=False),
                    top_n=getattr(args, 'ai_top_n', None),
                )
                ai_interpretation["backend"] = ai_backend
                ai_interpretation["model"] = getattr(interpreter.interpreter, "model", args.ai_model)

                # Save AI's reading as a JSON file.
                with open(output_dir / "ai_interpretation.json", 'w', encoding='utf-8') as f:
                    json.dump(ai_interpretation, f, ensure_ascii=False, indent=2)
                    f.write("\n")
            except Exception as exc:
                # AI is an optional interpretation layer; preserve the completed analysis.
                ai_interpretation = None
                ai_interpretation_error = {
                    "error_code": "AI_INTERPRETATION_FAILED",
                    "backend": ai_backend,
                    "mode": getattr(args, 'ai_mode', 'summary'),
                    "message": str(exc),
                }
                logger.error(
                    "AI interpretation failed; analysis results remain available [%s]: %s",
                    ai_interpretation_error["error_code"],
                    exc,
                )
                (output_dir / "ai_interpretation.json").unlink(missing_ok=True)
                with open(output_dir / "ai_interpretation_error.json", 'w', encoding='utf-8') as f:
                    json.dump(ai_interpretation_error, f, ensure_ascii=False, indent=2)
                    f.write("\n")

        # Generate the HTML analysis report.
        # Unless specified by the user--no-report, otherwise the HTML report is generated based on significant results
        if not args.no_report:
            logger.info("Generating HTML report...")
            report_gen = ReportGenerator(str(output_dir), config)

            # extract ssGSEA/GSVAActive score matrix (for visualization in the report)
            _gsva_scores_df = None
            if config.method in ('ssgsea', 'gsva') and results:
                import pandas as _pd
                for _db_name, _df in results.items():
                    if _df is not None and len(_df) > 0 and isinstance(_df, _pd.DataFrame):
                        _numeric_cols = _df.select_dtypes(include='number').columns
                        _non_metric_cols = {'p_value', 'FDR',
                                            'NOM p-val', 'FDR q-val', 'FWER p-val',
                                            'pvalue', 'P_Value', 'Adjusted_P_Value',
                                            'p.adjust', 'qvalues',
                                            'nes', 'es', 'fdr', 'gene_count', 'Gene_Count',
                                            'NES', 'enrichmentScore', 'setSize'}
                        _sample_cols = [c for c in _numeric_cols if c not in _non_metric_cols]
                        if _sample_cols:
                            _name_col = None
                            for col in ['Description', 'pathway', 'Term_Name', 'Term_ID', _df.index.name]:
                                if col and col in _df.columns:
                                    _name_col = col
                                    break
                            if _name_col:
                                _gsva_scores_df = _df.set_index(_name_col)[_sample_cols]
                            else:
                                _gsva_scores_df = _df[_sample_cols]
                                _gsva_scores_df.index.name = 'pathway'
                            break

            report_gen.generate(
                significant_results,
                str(output_dir / "report.html"),
                gene_list=list(gene_set),
                ai_interpretation=ai_interpretation,
                ai_interpretation_error=ai_interpretation_error,
                metadata=metadata,
                gsva_results=_gsva_scores_df,
                gsva_groups=groups_dict,
                analysis_method=config.method,
                plot_types=_effective_plot_types
            )

        # ----Print Analysis Summary----
        logger.info("=" * 50)
        logger.info("Analysis Complete!")
        logger.info("=" * 50)
        for db_name, df in results.items():
            unit = "activity pathways" if config.method in ('ssgsea', 'gsva') else "enriched terms"
            logger.info(f"  {db_name}: {len(df)} {unit}")
        logger.info(f"Results saved to: {output_dir}")

        return 0

    except FileNotFoundError as e:
        # Report missing inputs and database files as user-facing errors.
        logger.error(f"File not found. Check path: {e}")
        return 1

    except ValueError as e:
        # Parameter error: Invalid configuration parameter, empty file content, etc.
        logger.error(f"parameter error, check input: {e}")
        return 1

    except KeyboardInterrupt:
        # User break (Ctrl+C)
        logger.warning("Analysis interrupted by the user")
        return 130  # Unix Practice: 128 + SIGNT (2) = 130

    except ImportError as e:
        # Reliance library missing error
        logger.error(f"Lack of necessary library: {e}")
        logger.error("Install the optional analysis dependencies with `pip install 'allenricher[all]'`")
        return 1

    except Exception as e:
        # Generic anomaly capture: full error information recorded (including stack tracking)
        logger.error("Analysis failed unexpectedly: %s", e, exc_info=True)
        return 1


def cmd_download(args) -> int:
    """Download shared source data used to build species databases.

    GO, Reactome, and Disease Ontology files are stored under the database
    ``basic`` directory. DisGeNET is not downloaded because current releases are
    not freely redistributable; AllEnricher can reuse the frozen v20190612 snapshot
    built by v1.

    Args:
        args: Parsed arguments for the ``download`` subcommand.

    Returns:
        ``0`` on success and ``1`` on failure."""
    from allenricher.database.downloader import DataDownloader

    # TF resources have dedicated download and registry workflows.
    if getattr(args, 'trrust', False):
        return _cmd_download_trrust(args)
    if getattr(args, 'chea3', False):
        return _cmd_download_chea3(args)
    if getattr(args, 'animaltfdb', False):
        return _cmd_download_animaltfdb(args)

    databases = [d.strip().lower() for d in args.databases.split(',')]
    download_dir = args.database_dir or "./database"

    # Avoid replacing current snapshots unless an update or --force requires it.
    if not getattr(args, 'force', False):
        try:
            from allenricher.database.version import RemoteVersionChecker, DatabaseVersionManager
            checker = RemoteVersionChecker()
            vm = DatabaseVersionManager(download_dir)
            update_status = checker.check_updates(vm)

            # Download only sources with a newer remote snapshot.
            sources_with_update = [s for s, info in update_status.items() if info["has_update"]]
            sources_checked = list(update_status.keys())

            if sources_checked and not sources_with_update:
                logger.info("All checked data sources are up to date")
                logger.info("Use --force to replace an existing local database snapshot")
                return 0
            elif sources_with_update:
                logger.info("Updates are available for: %s", ", ".join(sources_with_update))
                logger.info("Starting download")
        except Exception as e:
            logger.warning("Could not check remote database versions: %s", e)

    logger.info("Downloading shared source data to %s", download_dir)
    logger.info("Databases: %s", ", ".join(databases))

    downloader = DataDownloader(
        root_dir=download_dir,
        overwrite=getattr(args, 'force', False),
        max_workers=getattr(args, 'workers', 4),
        use_multi_thread=not getattr(args, 'no_multi_thread', False),
        verify_integrity=not getattr(args, 'no_verify', False),
    )

    try:
        downloaded = downloader.download_all(databases)
        for db_type, path in downloaded.items():
            logger.info("Downloaded %s: %s", db_type, path)
        logger.info("Download completed")
        logger.info("")
        logger.info("Next: build a database of designated species")
        logger.info("  allenricher build -s hsa -t 9606 -d GO,Reactome")
        return 0
    except Exception as e:
        logger.error("Database download failed: %s", e)
        return 1


def _cmd_download_trrust(args) -> int:
    """Download TRRUST TF-target gene sets and update the species registry.

    Args:
        args: Parsed download arguments.

    Returns:
        ``0`` on success and ``1`` on failure."""
    from allenricher.database.downloader import DataDownloader
    from allenricher.database.trrust_fetcher import TRRUSTFetcher

    download_dir = args.database_dir or "./database"
    database_dir = download_dir.rstrip('/')

    logger.info("Downloading TRRUST to %s/basic/trrust/", database_dir)
    fetcher = TRRUSTFetcher(basic_dir=database_dir + "/basic")

    try:
        results = fetcher.download_all(overwrite=args.force)
        for name, path in results.items():
            logger.info(f"  {name}: {path}")
        registry_downloader = DataDownloader(root_dir=download_dir)
        registry_downloader.record_database_species(
            "TRRUST", fetcher.get_supported_species_records()
        )
        registry_downloader.refresh_supported_species_registry()
        logger.info("TRRUST download completed")
        return 0
    except Exception as e:
        logger.error("TRRUST download failed: %s", e)
        return 1


def _cmd_download_chea3(args) -> int:
    """Download local ChEA3 TF-target libraries and update the registry.

    Args:
        args: Parsed download arguments.

    Returns:
        ``0`` on success and ``1`` on failure."""
    from allenricher.database.chea3_fetcher import ChEA3Fetcher
    from allenricher.database.downloader import DataDownloader

    download_dir = args.database_dir or "./database"
    database_dir = download_dir.rstrip('/')

    logger.info("Downloading ChEA3 libraries to %s/basic/chea3/", database_dir)
    fetcher = ChEA3Fetcher(basic_dir=database_dir + "/basic")

    try:
        results = fetcher.download_all_gmt_libraries(overwrite=args.force)
        for name, path in results.items():
            logger.info(f"  {name}: {path}")
        registry_downloader = DataDownloader(root_dir=download_dir)
        registry_downloader.record_database_species(
            "ChEA3", fetcher.get_supported_species_records()
        )
        registry_downloader.refresh_supported_species_registry()
        logger.info("ChEA3 download completed")
        return 0
    except Exception as e:
        logger.error(f"ChEA3 download failed: {e}")
        return 1


def _cmd_download_animaltfdb(args) -> int:
    """Download AnimalTFDB resources and the human hTFtarget library."""
    from allenricher.database.animaltfdb_fetcher import AnimalTFDBFetcher
    from allenricher.database.downloader import DataDownloader
    from allenricher.database.htftarget_fetcher import HTFtargetFetcher

    download_dir = args.database_dir or "./database"
    fetcher = AnimalTFDBFetcher(basic_dir=download_dir.rstrip('/') + "/basic")
    registry_downloader = DataDownloader(root_dir=download_dir)

    print("Downloading the human hTFtarget TF-target library...")
    fetcher.download_htftarget(overwrite=args.force)
    registry_downloader.record_database_species(
        "hTFtarget", HTFtargetFetcher.get_supported_species_records()
    )
    registry_downloader.record_database_species(
        "AnimalTFDB", fetcher.get_supported_species_records()
    )
    registry_downloader.refresh_supported_species_registry()

    species_list = args.species.split(',') if args.species else []

    if not species_list:
        print("No AnimalTFDB species were requested; only hTFtarget was downloaded.")
        print("Use --species Bos_taurus,Sus_scrofa to request AnimalTFDB species data.")
        return 0

    print(f"Downloading AnimalTFDB data for {len(species_list)} species...")
    for sp in species_list:
        print(f"\n--- Downloading {sp} ---")
        fetcher.download_species_data(sp, overwrite=args.force)

    print("\nAnimalTFDB download completed")
    return 0


def cmd_build(args) -> int:
    """Build analysis-ready databases for one species.

    Shared source data are filtered by species and written under
    ``database/organism``. User annotations can be built as GO, KEGG, or custom
    gene-set databases through the same command.

    Args:
        args: Parsed arguments for the ``build`` subcommand.

    Returns:
        ``0`` on success and ``1`` on failure."""
    from allenricher.database.builder import DatabaseBuilder

    databases = [d.strip().upper() for d in args.databases.split(',')]
    build_dir = args.database_dir or "./database"

    species = args.species
    taxid = args.taxonomy

    # Build user-supplied annotations before processing standard databases.
    has_custom_annotation = getattr(args, 'go_annot', None) or \
                           getattr(args, 'kegg_annot', None) or \
                           getattr(args, 'custom_annot', None)
    custom_built_databases = set()

    if has_custom_annotation:
        try:
            from allenricher.database.custom_builder import CustomDatabaseBuilder
            custom_builder = CustomDatabaseBuilder(root_dir=build_dir)

            if getattr(args, 'go_annot', None):
                fmt = None if getattr(args, 'annot_format', 'auto') == 'auto' else args.annot_format
                custom_builder.build_from_annotation(
                    annotation_file=args.go_annot,
                    species=species,
                    taxid=taxid,
                    db_name='GO',
                    format_type=fmt,
                    hierarchy_separator=getattr(args, 'hierarchy_sep', '|')
                )
                custom_built_databases.add('GO')

            if getattr(args, 'kegg_annot', None):
                fmt = None if getattr(args, 'annot_format', 'auto') == 'auto' else args.annot_format
                custom_builder.build_from_annotation(
                    annotation_file=args.kegg_annot,
                    species=species,
                    taxid=taxid,
                    db_name='KEGG',
                    format_type=fmt,
                    hierarchy_separator=getattr(args, 'hierarchy_sep', '|')
                )
                custom_built_databases.add('KEGG')

            if getattr(args, 'custom_annot', None):
                db_name = getattr(args, 'custom_db_name', 'CUSTOM')
                fmt = None if getattr(args, 'annot_format', 'auto') == 'auto' else args.annot_format
                custom_builder.build_from_annotation(
                    annotation_file=args.custom_annot,
                    species=species,
                    taxid=taxid,
                    db_name=db_name,
                    format_type=fmt,
                    hierarchy_separator=getattr(args, 'hierarchy_sep', '|')
                )
                custom_built_databases.update({'CUSTOM', db_name.upper()})
        except ImportError:
            print("Warning: CustomDatabaseBuilder not available. Skipping custom annotation build.")

    # ----Standard build process----
    databases = [database for database in databases if database not in custom_built_databases]
    if not databases:
        logger.info("Custom Quetz database construction completed and no standard database builder is required")
        return 0
    # Try to retrieve species information from SpeciesRegistry
    try:
        from .database.species_registry import SpeciesRegistry
        registry = SpeciesRegistry.load_default()
        entry = registry.query_by_taxid(taxid)
        if entry:
            logger.info(f"Species found in registry: {entry.latin_name} (TaxID: {entry.taxid})")
            logger.info(f"  Database support - GO: {'Yes' if entry.has_go else 'No'}, "
                       f"KEGG: {'Yes' if entry.has_kegg else 'No'}, "
                       f"DO: {'Yes' if entry.has_do else 'No'}, "
                       f"Reactome: {'Yes' if entry.has_reactome else 'No'}")
            # Prefer the scientific name recorded in the species registry.
            species_display_name = entry.latin_name
        else:
            species_display_name = species
    except Exception:
        # Stand back quietly and use the species code from the command line.
        species_display_name = species

    logger.info(f"Build a unique database of species: {species_display_name} (TaxID: {taxid})")
    logger.info(f"Database root directory: {build_dir}")
    logger.info(f"Database to be constructed: {', '.join(databases)}")

    builder = DatabaseBuilder(root_dir=build_dir)

    try:
        build_kwargs = dict(
            species=species,
            taxid=taxid,
            databases=databases
        )
        # If specified, pass to Builder
        latin_name = getattr(args, 'latin_name', '')
        if latin_name:
            build_kwargs['latin_name'] = latin_name
        outdir = builder.build_species_db(**build_kwargs)
        logger.info(f"Build complete! Output directory: {outdir}")
        logger.info("")
        logger.info("Next: Run enrichment analysis")
        logger.info(f"  allenricher analyze -i genes.txt -s {args.species} --database-dir {outdir}")
        return 0
    except FileNotFoundError as e:
        logger.error(f"Basic data not found: {e}")
        logger.info("Download the required source data first:")
        logger.info("  allenricher download -d go,reactome")
        return 1
    except Exception as e:
        logger.error("Database build failed: %s", e)
        return 1


def cmd_serve(args) -> int:
    """Start the local AllEnricher API and Web application.

    Args:
        args: Parsed arguments for the ``serve`` subcommand.

    Returns:
        ``0`` when the server exits normally and ``1`` on startup failure."""
    if getattr(args, 'config', None):
        os.environ['ALLENRICHER_CONFIG'] = args.config
    logger.info(f"Starting API server on {args.host}:{args.port}")

    # Delay import API server module to avoid loading dependency when not required
    from allenricher.api.server import start_api
    start_api(host=args.host, port=args.port, reload=args.reload)

    return 0



def _record_value(record, key: str, default=None):
    """Return a field from a registry dataclass or detail dictionary."""
    if isinstance(record, dict):
        return record.get(key, default)
    return getattr(record, key, default)


def _visible_database_support(record, database: str, support_key: str) -> bool:
    """Mirror runtime database/species validation before displaying support."""
    if not _record_value(record, support_key, False):
        return False
    db_name = str(database).split()[0]
    species_code = str(_record_value(record, 'kegg_code', '') or '').strip().lower()
    taxid = _record_value(record, 'taxid', None)
    if not species_code:
        species_code = {9606: 'hsa', 10090: 'mmu'}.get(taxid, str(taxid or '').lower())
    try:
        validate_tf_database_species(db_name, species_code)
    except ValueError:
        return False
    return True

def cmd_list(args) -> int:
    """List supported species or database resources.

    Args:
        args: Parsed arguments for the ``list`` subcommand.

    Returns:
        ``0`` on success and ``1`` for an unsupported request."""
    if args.resource == 'species':
        # List of supported species
        # Select to load from SpeciesRegistry, and then back to SPECIES_CONFIGS
        registry_loaded = False
        try:
            from .database.species_registry import SpeciesRegistry
            registry = SpeciesRegistry.load_default()
            if registry and registry.entries:
                registry_loaded = True
                print("\nSupported Species (from SpeciesRegistry):")
                print("-" * 165)
                columns = [
                    ('GO', 'has_go'), ('KEGG', 'has_kegg'), ('Reactome', 'has_reactome'),
                    ('DO', 'has_do'), ('DisGeNET', 'has_disgenet'), ('WikiPathways', 'has_wikipathways'),
                    ('TRRUST', 'has_trrust'), ('ChEA3', 'has_chea3'),
                    ('AnimalTFDB', 'has_animaltfdb'), ('hTFtarget', 'has_htftarget'),
                ]
                print(f"{'Code':<10} {'Name':<30} {'TaxID':<10} " + ' '.join(f"{label:<12}" for label, _ in columns))
                print("-" * 165)
                for entry in registry.entries.values():
                    # Use the KEGG organism code when available; otherwise display the TaxID.
                    code = entry.kegg_code if entry.has_kegg else str(entry.taxid)
                    name = entry.latin_name if entry.latin_name else f"TaxID {entry.taxid}"
                    statuses = ' '.join(f"{('Yes' if _visible_database_support(entry, label, key) else 'No'):<12}" for label, key in columns)
                    print(f"{code:<10} {name:<30} {entry.taxid:<10} {statuses}")
        except Exception:
            pass

        if not registry_loaded:
            from allenricher.core.config import SPECIES_CONFIGS

            print("\nSupported Species:")
            print("-" * 50)
            print(f"{'Code':<10} {'Name':<25} {'Taxonomy ID':<12}")
            print("-" * 50)
            for code, config in SPECIES_CONFIGS.items():
                print(f"{code:<10} {config.display_name:<25} {config.taxonomy_id:<12}")

    elif args.resource == 'databases':
        # List the types of databases supported
        from allenricher.core.config import DATABASE_CATALOG

        print("\nSupported Databases:")
        print("-" * 50)
        for database in DATABASE_CATALOG:
            print(f"  - {database['name']}")

    return 0


def cmd_config(args) -> int:
    """Write an example AllEnricher configuration file.

    Args:
        args: Parsed arguments for the ``config`` subcommand.

    Returns:
        ``0`` on success and ``1`` on failure."""
    # Create default profile and write to file
    config = Config()
    config.to_file(args.output)  # Convert Configuration Sequence into YAML File
    logger.info(f"Configuration file generated: {args.output}")
    return 0


def cmd_check_update(args) -> int:
    """Check whether newer remote database snapshots are available.

    Args:
        args: Parsed arguments for the ``check-update`` subcommand.

    Returns:
        ``0`` on success and ``1`` when the check cannot be completed."""
    import json
    from allenricher.database.version import RemoteVersionChecker, DatabaseVersionManager

    db_dir = _resolve_db_dir(args)
    checker = RemoteVersionChecker()
    vm = DatabaseVersionManager(db_dir)

    logger.info(f"Checking remote data source updates (Database Directory: {db_dir})...")
    update_status = checker.check_updates(vm)

    # Print a concise human-readable update table.
    print("\nRemote data source update check")
    print("=" * 80)
    print(f"  {'Data Sources':<20} {'Status':<10} {'Local version':<25} {'Remote version'}")
    print(f"  {'-'*20} {'-'*10} {'-'*25} {'-'*25}")

    has_any_update = False
    for source, info in sorted(update_status.items()):
        if info["has_update"]:
            status = "Update available"
            has_any_update = True
        else:
            status = "Latest"
        local_ver = info["local"].get("remote_version") or info["local"].get("version") or "-"
        remote_ver = info["remote"].get("remote_version") or "-"
        print(f"  {source:<20} {status:<10} {local_ver:<25} {remote_ver}")

    print("=" * 80)

    if has_any_update:
        print("Update available. Run `allenricher download --force` to refresh local data.")
    else:
        print("All local data sources are current.")

    # Preserve the readable table and append JSON when requested.
    if args.json:
        print("\n--- JSON Output ---")
        print(json.dumps(update_status, indent=2, ensure_ascii=False))

    return 0


def cmd_cleanup(args) -> int:
    """Remove old database versions according to the retention policy.

    Args:
        args: Parsed arguments for the ``cleanup`` subcommand.

    Returns:
        ``0`` on success and ``1`` on failure."""
    from allenricher.database.version import DatabaseVersionManager

    db_dir = _resolve_db_dir(args)
    vm = DatabaseVersionManager(db_dir)

    if args.dry_run:
        logger.info(f"[dry-run] Preview cleanup operation (Database directory: {db_dir}, with the latest version of {args.keep})")
    else:
        logger.info(f"Performing cleanup operations (database directory: {db_dir}Keep the update.{args.keep}Version)")

    removed = vm.remove_stale_versions(keep_count=args.keep, dry_run=args.dry_run)

    if not removed:
        print("No stale database versions are eligible for cleanup.")
        return 0

    print("\nClean-up results:")
    print("-" * 60)
    total_count = 0
    for source, versions in removed.items():
        if versions:
            action = "Delete" if args.dry_run else "Deleted"
            print(f"  [{source}] {action}: {', '.join(versions)}")
            total_count += len(versions)

    print("-" * 60)
    action = "Delete" if args.dry_run else "Deleted"
    print(f"Total{action} {total_count}An old version of the directory.")

    if args.dry_run:
        print("\nTip: This is preview mode. Remove it.--dry-runThe parameters are cleared for actual implementation.")

    return 0


def cmd_list_versions(args) -> int:
    """List locally installed database versions.

    Args:
        args: Parsed arguments for the ``list-versions`` subcommand.

    Returns:
        ``0`` on success and ``1`` on failure."""
    import json
    from allenricher.database.version import DatabaseVersionManager

    db_dir = _resolve_db_dir(args)
    vm = DatabaseVersionManager(db_dir)

    if args.json:
        print(json.dumps(vm.get_summary_json(), indent=2, ensure_ascii=False))
    elif args.lineage:
        print(vm.get_full_lineage_report())
    else:
        print(vm.get_summary_table())

    return 0


def main():
    """Parse command-line arguments and dispatch the selected subcommand.

    Args:
        argv: Optional argument vector. ``sys.argv`` is used when omitted.

    Returns:
        The subcommand exit code."""
    # Create Parameter parser and parsing command line parameters
    parser = create_parser()
    argv = sys.argv[1:]
    args = parser.parse_args(argv)
    args._provided_options = {
        token.split("=", 1)[0] for token in argv if token.startswith("-")
    }

    # If the user does not specify a sub-command, print help messages and exit
    if args.command is None:
        parser.print_help()
        return 0

    # Build sub-command name to map of the processing function
    commands = {
        'analyze': cmd_analyze,
        'download': cmd_download,
        'build': cmd_build,
        'serve': cmd_serve,
        'list': cmd_list,
        'config': cmd_config,
        'list-species': _cmd_list_species,
        'query-species': _cmd_query_species,
        'tf-enrich': _cmd_tf_enrich,
        'check-update': cmd_check_update,
        'cleanup': cmd_cleanup,
        'list-versions': cmd_list_versions,
    }

    # Dispatch to the selected subcommand handler.
    handler = commands.get(args.command)
    if handler:
        return handler(args)  # Call the processing function and return its exit code
    else:
        parser.print_help()   # Unknown command, print help
        return 1


def _cmd_list_species(args) -> int:
    """List species recorded in the unified species registry."""
    from .database.species_registry import SpeciesRegistry
    import json

    registry = SpeciesRegistry.load_default()

    if args.summary:
        summary = registry.get_summary()
        print(f"\n{'='*50}")
        print("Species Registry Summary")
        print(f"{'='*50}")
        print(f"Total species: {summary['total_species']:,}")
        print(f"\nBy Database:")
        for db, stats in summary.items():
            if db != 'total_species':
                print(f"  - {db.upper()}: {stats['count']:,} species")
        print(f"{'='*50}\n")
        return 0

    entries = registry.filter_by_databases(
        go=args.go or None,
        kegg=args.kegg or None,
        reactome=args.reactome or None,
        do=args.do or None,
        disgenet=args.disgenet or None,
        wikipathways=args.wikipathways or None,
        trrust=args.trrust or None,
        chea3=args.chea3 or None,
        animaltfdb=args.animaltfdb or None,
        htftarget=args.htftarget or None,
    )

    if args.format == "table":
        # When used--trrustor--chea3Show Species, Code, TaxID rows when filtering
        if args.disgenet or args.trrust or args.chea3 or args.animaltfdb or args.htftarget:
            print(f"{'Species':<30} {'Code':<8} {'TaxID':<10}")
            print("-" * 50)
            for e in entries[:100]:
                code = e.kegg_code or "N/A"
                print(f"{e.latin_name:<30} {code:<8} {e.taxid:<10}")
            if len(entries) > 100:
                print(f"... and {len(entries) - 100} more species")
        # Show Data Type column when filtering with --wikipathways
        elif args.wikipathways:
            print(f"{'Species':<30} {'Code':<8} {'Data Type':<10} {'TaxID':<10}")
            print("-" * 62)
            for e in entries[:100]:
                data_type = e.wikipathways_data_type or "N/A"
                code = e.kegg_code or "N/A"
                print(f"{e.latin_name:<30} {code:<8} {data_type:<10} {e.taxid:<10}")
            if len(entries) > 100:
                print(f"... and {len(entries) - 100} more species")
        else:
            columns = [
                ('GO', 'has_go'), ('KEGG', 'has_kegg'), ('Reactome', 'has_reactome'),
                ('DO', 'has_do'), ('DisGeNET', 'has_disgenet'), ('WikiPathways', 'has_wikipathways'),
                ('TRRUST', 'has_trrust'), ('ChEA3', 'has_chea3'),
                ('AnimalTFDB', 'has_animaltfdb'), ('hTFtarget', 'has_htftarget'),
            ]
            print(f"{'TaxID':<10} {'Latin Name':<30} " + ' '.join(f"{label:<12}" for label, _ in columns))
            print("-" * 165)
            for e in entries[:100]:
                statuses = ' '.join(f"{('Y' if _visible_database_support(e, label, key) else 'N'):<12}" for label, key in columns)
                print(f"{e.taxid:<10} {e.latin_name:<30} {statuses}")
            if len(entries) > 100:
                print(f"... and {len(entries) - 100} more species")
    elif args.format == "tsv":
        if args.disgenet or args.trrust or args.chea3 or args.animaltfdb or args.htftarget:
            print("species\tcode\ttaxid")
            for e in entries:
                code = e.kegg_code or "N/A"
                print(f"{e.latin_name}\t{code}\t{e.taxid}")
        elif args.wikipathways:
            print("species\tcode\tdata_type\ttaxid")
            for e in entries:
                data_type = e.wikipathways_data_type or "N/A"
                code = e.kegg_code or "N/A"
                print(f"{e.latin_name}\t{code}\t{data_type}\t{e.taxid}")
        else:
            print("taxid\tlatin_name\thas_go\thas_kegg\thas_reactome\thas_do\thas_disgenet\thas_wikipathways\thas_trrust\thas_chea3\thas_animaltfdb\thas_htftarget")
            for e in entries:
                print(f"{e.taxid}\t{e.latin_name}\t{e.has_go}\t{e.has_kegg}\t{e.has_reactome}\t{e.has_do}\t{e.has_disgenet}\t{e.has_wikipathways}\t{e.has_trrust}\t{e.has_chea3}\t{e.has_animaltfdb}\t{e.has_htftarget}")
    elif args.format == "json":
        if args.disgenet or args.trrust or args.chea3 or args.animaltfdb or args.htftarget:
            data = [{"species": e.latin_name, "code": e.kegg_code or "N/A", "taxid": e.taxid}
                    for e in entries]
        elif args.wikipathways:
            data = [{"species": e.latin_name, "code": e.kegg_code or "N/A",
                     "data_type": e.wikipathways_data_type or "N/A", "taxid": e.taxid}
                    for e in entries]
        else:
            data = [{"taxid": e.taxid, "latin_name": e.latin_name, "has_go": e.has_go,
                     "has_kegg": e.has_kegg, "has_reactome": e.has_reactome, "has_do": e.has_do,
                     "has_disgenet": e.has_disgenet,
                     "has_wikipathways": e.has_wikipathways, "has_trrust": e.has_trrust,
                     "has_chea3": e.has_chea3, "has_animaltfdb": e.has_animaltfdb,
                     "has_htftarget": e.has_htftarget}
                    for e in entries]
        print(json.dumps(data, indent=2))

    return 0


def _print_local_build_status(taxid: int, kegg_code: str = None) -> None:
    """Print locally built databases available for a species."""
    from .database.version import DatabaseVersionManager
    from pathlib import Path

    species_code = kegg_code

    # Determine Database Directory
    db_dir = Path("database")
    if not db_dir.exists():
        for candidate in [Path.cwd() / "database", Path(__file__).parent.parent.parent / "database"]:
            if candidate.exists():
                db_dir = candidate
                break

    manager = DatabaseVersionManager(database_dir=str(db_dir))
    org_versions = manager.list_installed_organism_versions()

    print(f"\nLocal Database Status:")
    print(f"-" * 60)

    if not org_versions:
        print(f"  Not built (no organism database found)")
        return

    # Find the species in all oranism versions
    found_builds = []
    for version in org_versions:
        build_info = manager.get_organism_build_info(version)
        species_list = build_info.get(version, [])

        matched = False
        if species_code and species_code in species_list:
            matched = True
        elif species_code is None:
            for sp in species_list:
                lineage = manager.get_build_lineage(version, sp)
                if lineage and lineage.get("taxid") == taxid:
                    species_code = sp
                    matched = True
                    break

        if matched:
            lineage = manager.get_build_lineage(version, species_code)
            found_builds.append((version, species_code, lineage))

    if not found_builds:
        print(f"  Not built")
    else:
        for version, sp, lineage in found_builds:
            print(f"  Version: {version}")
            print(f"  Species Code: {sp}")
            if lineage:
                built_at = lineage.get("built_at", "-")
                if built_at and built_at != "-":
                    built_at = built_at[:16]
                print(f"  Built at: {built_at}")
                databases = lineage.get("databases", [])
                if databases:
                    print(f"  Databases: {', '.join(databases)}")
                deps = lineage.get("dependencies", {})
                if deps:
                    for db_name, dep_info in deps.items():
                        basic_dir = dep_info.get("basic_dir", "-")
                        print(f"    {db_name} <- {basic_dir}")


def _print_species_detail(registry, entry, match_kind: str = "") -> None:
    """Print a concise species-registry record."""
    detail = registry.get_species_detail(entry.taxid)

    match_info = f" [{match_kind}]" if match_kind and match_kind != 'exact' else ""
    print(f"\n{'='*60}")
    print(f"Species Information{match_info}")
    print(f"{'='*60}")
    print(f"Taxonomy ID: {detail['taxid']}")
    print(f"Latin Name:  {detail['latin_name']}")
    if detail.get('common_name') and detail['common_name'] != '-':
        print(f"Common Name: {detail['common_name']}")
    if detail.get('synonyms') and detail['synonyms'] != '-':
        syn_list = detail['synonyms'].split(';')
        syn_display = [s for s in syn_list if s.strip() != detail['latin_name']]
        if syn_display:
            print(f"Other Names: {', '.join(syn_display[:5])}")
            if len(syn_display) > 5:
                print(f"  ... and {len(syn_display) - 5} more")
    print(f"\nDatabase Support:")
    print(f"-" * 60)

    for db_name, db_key in [
        ("GO", "has_go"),
        ("KEGG", "has_kegg"),
        ("Reactome", "has_reactome"),
        ("DO", "has_do"),
        ("DisGeNET (v20190612)", "has_disgenet"),
        ("WikiPathways", "has_wikipathways"),
        ("TRRUST", "has_trrust"),
        ("ChEA3", "has_chea3"),
        ("AnimalTFDB", "has_animaltfdb"),
        ("hTFtarget", "has_htftarget"),
    ]:
        if _visible_database_support(detail, db_name, db_key):
            print(f"\n{db_name}: Supported")
            if db_name == "GO" and detail.get('go_source') and detail['go_source'] != '-':
                print(f"  Source: {detail['go_source']}")
            if db_name == "KEGG" and detail.get('kegg_code') and detail['kegg_code'] != '-':
                print(f"  Code: {detail['kegg_code']} (source: {detail.get('kegg_code_source', '-')})")

    # Local database build state
    _print_local_build_status(detail['taxid'], detail.get('kegg_code'))

    print(f"\nBuild Command:")
    print(f"  allenricher build --taxonomy {detail['taxid']}")
    print(f"{'='*60}\n")


def _cmd_query_species(args) -> int:
    """Query the species registry by Taxonomy ID, name, or KEGG code."""
    from .database.species_registry import SpeciesRegistry

    registry = SpeciesRegistry.load_default()

    entries = []
    match_type = ""

    if args.taxid:
        entry = registry.query_by_taxid(args.taxid)
        if entry:
            entries = [(entry, 1.0, 'exact')]
        match_type = f"TaxID={args.taxid}"
    elif args.name:
        query_name = args.name.strip()
        entries = registry.fuzzy_search(query_name, cutoff=0.5)
        match_type = f"Name='{query_name}'"
    elif args.kegg:
        entry = registry.query_by_kegg_code(args.kegg)
        if entry:
            entries = [(entry, 1.0, 'exact')]
        match_type = f"KEGG={args.kegg}"

    if not entries:
        print(f"Species not found in registry ({match_type}).")
        return 1

    # If there is only one match, show the details.
    if len(entries) == 1:
        entry, score, match_kind = entries[0]
        _print_species_detail(registry, entry, match_kind)
        return 0

    # Multiple matches, display list for user selection
    print(f"\n{'='*60}")
    print(f"Found {len(entries)} matching species for {match_type}:")
    print(f"{'='*60}")

    for i, (entry, score, match_kind) in enumerate(entries[:10], 1):
        match_info = f"[{match_kind}]" if match_kind != 'exact' else ""
        print(f"{i}. {entry.latin_name} (TaxID: {entry.taxid}) {match_info}")
        if entry.common_name and entry.common_name != '-':
            print(f"   Common: {entry.common_name}")
        dbs = []
        for db_name, db_key in [
            ("GO", "has_go"), ("KEGG", "has_kegg"), ("Reactome", "has_reactome"),
            ("DO", "has_do"), ("DisGeNET (v20190612)", "has_disgenet"),
            ("WikiPathways", "has_wikipathways"), ("TRRUST", "has_trrust"),
            ("ChEA3", "has_chea3"), ("AnimalTFDB", "has_animaltfdb"),
            ("hTFtarget", "has_htftarget"),
        ]:
            if _visible_database_support(entry, db_name, db_key):
                dbs.append(db_name)
        print(f"   Databases: {', '.join(dbs) if dbs else 'None'}")
        print()

    if len(entries) > 10:
        print(f"... and {len(entries) - 10} more matches.")

    print(f"Showing details for best match: {entries[0][0].latin_name}")
    print(f"{'='*60}\n")
    _print_species_detail(registry, entries[0][0], entries[0][2])

    return 0


def _convert_api_result_to_df(api_result: Dict, gene_set_size: int) -> "pd.DataFrame":
    """Convert a ChEA3 API response into the standard TF result table.

    The returned columns match the local TF enrichment workflow so downstream
    reporting and plotting do not need a separate code path."""
    import pandas as pd
    import numpy as np
    from statsmodels.stats.multitest import multipletests

    all_entries = []

    for lib_name, entries in api_result.items():
        for entry in entries:
            # Parsing Overlap Fields (Form: "3"/100"or "3"")
            overlap_str = str(entry.get('Overlap', '0'))
            if '/' in overlap_str:
                overlap = int(overlap_str.split('/')[0])
            else:
                overlap = int(overlap_str) if overlap_str.isdigit() else 0

            # Parsing TargetCount
            target_count_str = str(entry.get('TargetCount', '0'))
            if '/' in target_count_str:
                target_count = int(target_count_str.split('/')[1]) if '/' in target_count_str else int(target_count_str.split('/')[0])
            else:
                target_count = int(target_count_str) if target_count_str.isdigit() else 0

            # Pvalue
            pvalue_str = str(entry.get('Pvalue', '1.0'))
            try:
                pvalue = float(pvalue_str)
            except (ValueError, TypeError):
                pvalue = 1.0

            all_entries.append({
                'Term_ID': f"{lib_name}|{entry.get('TF', '')}",
                'Term_Name': f"{entry.get('TF', '')} [{lib_name}]",
                'TF': str(entry.get('TF', '')),
                'Overlap': overlap,
                'TF_Targets': target_count,
                'GeneSet_Size': gene_set_size,
                'Overlap_Genes': '',  # API does not return specific overlapping genes
                'Pvalue': pvalue,
                'Mode': 'unknown',  # API does not return the mode of control
                'Library': lib_name,
                'Context': '',
                'Evidence_Type': '',
                'Inference_Type': 'direct',
                '_Rank': int(entry.get('Rank', 999)) if str(entry.get('Rank', '')).isdigit() else 999,
            })

    if not all_entries:
        return pd.DataFrame(columns=['TF', 'Overlap', 'TF_Targets', 'GeneSet_Size',
                                     'Overlap_Genes', 'Pvalue', 'FDR', 'Mode', 'Library'])

    result_df = pd.DataFrame(all_entries)

    # ChEA3 Independent correction of the libraries; not cross-bank integration of the same name TF.
    result_df['FDR'] = result_df.groupby('Library')['Pvalue'].transform(
        lambda values: multipletests(values.to_numpy(), method='fdr_bh')[1]
        if len(values) > 1 else values.to_numpy()
    )
    result_df['Target_Set_Size'] = result_df['TF_Targets']
    return result_df.sort_values(['Pvalue', 'Library', 'TF']).drop(
        columns=['_Rank']
    ).reset_index(drop=True)


def _tf_size_limits(args, method: str):
    from allenricher.analysis.tf_enrichment import resolve_tf_size_limits

    return resolve_tf_size_limits(
        method,
        getattr(args, 'tf_min_size', None),
        getattr(args, 'tf_max_size', None),
    )


def _load_tf_database(manager, database: str, species: str):
    loaders = {
        'trrust': manager.load_trrust,
        'chea3': manager.load_chea3,
        'animaltfdb': manager.load_animaltfdb,
        'htftarget': manager.load_htftarget,
    }
    return loaders[database](species=species)


def _tf_rank_consensus(results_df: pd.DataFrame, method: str) -> pd.DataFrame:
    from allenricher.analysis.tf_meta_analyzer import TFMetaAnalyzer

    source_columns = [column for column in ('TF_Database', 'Library') if column in results_df]
    if not source_columns:
        return pd.DataFrame()
    source = results_df[source_columns].fillna('').astype(str).agg('|'.join, axis=1)
    frames = {
        name.strip('|'): results_df.loc[source == name]
        for name in sorted(source.unique())
    }
    return TFMetaAnalyzer.rank_consensus(frames, method=method)


def _run_tf_analysis(
    args,
    gene_list: List[str],
    species: str,
    background_genes: Optional[Set[str]] = None,
) -> Optional[pd.DataFrame]:
    """Run local TF enrichment against one or more TF-target libraries.

    Args:
        gene_list: Query genes.
        databases: TF databases to analyze.
        species: Species identifier used to validate database support.
        method: ``ora`` or ``gsea``.

    Returns:
        Combined TF enrichment results."""
    import pandas as pd
    from allenricher.analysis.tf_enrichment import TFEnrichmentAnalyzer
    from allenricher.database.manager import DatabaseManager

    tf_database_choice = args.tf_database
    database_dir = args.database_dir or "./database"

    # Determines the list of databases to be analysed
    if tf_database_choice == 'both':
        tf_databases = ['trrust', 'chea3']
    else:
        tf_databases = [tf_database_choice]

    all_results = []
    for tf_db in tf_databases:
        logger.info(f"Loading TF database: {tf_db}")

        # Load the TF database.
        db_manager = DatabaseManager(database_dir, species)
        tf_database = _load_tf_database(db_manager, tf_db, species)

        if tf_database is None:
            logger.warning(
                "Cannot load %s for species %s; skipping this database. "
                "Build or download the database before retrying.",
                tf_db,
                species,
            )
            continue

        # Execute ORA analysis
        analyzer = TFEnrichmentAnalyzer(tf_database=tf_database)
        min_size, max_size = _tf_size_limits(args, 'ora')
        results_df = analyzer.ora(
            gene_set=gene_list,
            library=getattr(args, 'tf_library', None),
            tissue=getattr(args, 'tf_tissue', None),
            regulation=getattr(args, 'tf_regulation', 'all'),
            min_size=min_size,
            max_size=max_size,
            background_genes=background_genes,
        )
        if results_df is not None and not results_df.empty:
            # Record the source database.
            results_df['TF_Database'] = tf_db.upper()
            all_results.append(results_df)
            logger.info("%s analysis completed: %d TFs", tf_db.upper(), len(results_df))
        else:
            logger.warning("No significant transcription factors were found in %s", tf_db.upper())

    if not all_results:
        return None

    # Merge all results
    combined_df = pd.concat(all_results, ignore_index=True)
    # Sort by Pvalue
    if 'Pvalue' in combined_df.columns:
        combined_df = combined_df.sort_values('Pvalue').reset_index(drop=True)

    return combined_df


def _load_tf_gene_file(path: Path) -> List[str]:
    """Load an ORA gene list while preserving first-occurrence order."""
    genes: List[str] = []
    seen = set()
    with path.open('r', encoding='utf-8-sig') as handle:
        for line_number, raw_line in enumerate(handle, start=1):
            gene = raw_line.strip()
            if not gene or gene.startswith('#'):
                continue
            if any(separator in gene for separator in ('\t', ',', ';')):
                raise ValueError(
                    f"Line {line_number} contains multiple fields; TF ORA requires one gene ID per line"
                )
            if gene not in seen:
                seen.add(gene)
                genes.append(gene)
    if not genes:
        raise ValueError(f"Gene input file is empty: {path}")
    return genes


def _cmd_tf_enrich(args) -> int:
    """Run the dedicated transcription-factor enrichment workflow.

    Args:
        args: Parsed arguments for the ``tf-enrich`` subcommand.

    Returns:
        ``0`` on success and ``1`` on failure."""
    import pandas as pd
    from allenricher.analysis.tf_enrichment import TFEnrichmentAnalyzer
    from allenricher.core.enrichment import EnrichmentAnalyzer
    from allenricher.database.manager import DatabaseManager
    from allenricher.report.visualizer import Visualizer

    # Parsing parameters
    input_file = args.input
    species = args.species
    database = args.database
    output_dir = Path(args.output)
    top_n = args.top_n
    database_dir = args.database_dir or "./database"
    online_mode = getattr(args, 'online', False)
    method = args.method.lower()

    try:
        validate_tf_database_species(database, species)
    except ValueError as exc:
        logger.error(str(exc))
        return 1

    # Ensure that the output directory exists
    output_dir.mkdir(parents=True, exist_ok=True)

    # Reads the input for the method. GSEA must use the true sorting statistics and cannot be forged by the line number.
    input_path = Path(input_file)
    if not input_path.exists():
        logger.error(f"The input file does not exist: {input_file}")
        return 1

    if method == 'gsea' and online_mode:
        logger.error("--online uses the ChEA3 ORA endpoint and does not support --method gsea")
        return 1

    background_file = getattr(args, 'background', None)
    if method == 'gsea' and background_file:
        logger.error("--background applies to TF ORA only; TF GSEA uses the complete ranked gene list")
        return 1
    if online_mode and background_file:
        logger.error("The ChEA3 online API does not accept a custom background; use local TF ORA instead")
        return 1

    ranked_genes = None
    if method == 'gsea':
        try:
            ranked_genes = EnrichmentAnalyzer.load_ranked_gene_list(input_file)
        except (OSError, ValueError) as exc:
            logger.error(f"Could not load GSEA ranking file: {exc}")
            return 1
        gene_list = [gene for gene, _ in ranked_genes]
    else:
        try:
            gene_list = _load_tf_gene_file(input_path)
        except (OSError, ValueError) as exc:
            logger.error("Could not load the TF ORA gene file: %s", exc)
            return 1

    if not gene_list:
        logger.error(f"Enter a list of genes empty: {input_file}")
        return 1

    logger.info(f"Number of genes entered: {len(gene_list)}")
    logger.info(f"Species: {species}, TF database: {database}")

    # Online Applied Mode (ChEA3 only)
    if online_mode and database == 'chea3':
        logger.info("Use ChEA3 API for online analysis...")
        from allenricher.database.chea3_fetcher import ChEA3Fetcher
        from allenricher.database.parsers.chea3 import ChEA3Parser

        try:
            fetcher = ChEA3Fetcher(basic_dir=database_dir)
            query_name = Path(input_file).stem
            api_result = fetcher.enrich_api(gene_list, query_name=query_name)

            # Parsing API Results
            parsed_result = ChEA3Parser.parse_api_result(api_result)

            # Convert to Standard DataFrame
            results_df = _convert_api_result_to_df(parsed_result, len(gene_list))
            tf_library = getattr(args, 'tf_library', None)
            if tf_library:
                libraries = {item.strip().lower() for item in tf_library.split(',')}
                results_df = results_df[
                    results_df['Library'].str.lower().isin(libraries)
                ].copy()
            min_size, max_size = _tf_size_limits(args, 'ora')
            before_size_filter = len(results_df)
            below_min = int((results_df['TF_Targets'] < min_size).sum())
            above_max = int((results_df['TF_Targets'] > max_size).sum()) if max_size is not None else 0
            results_df = results_df[results_df['TF_Targets'] >= min_size]
            if max_size is not None:
                results_df = results_df[results_df['TF_Targets'] <= max_size]
            logger.info(
                "ChEA3 online TF ORA size filter: min=%d, max=%s, before=%d, "
                "after=%d, below_min=%d, above_max=%d",
                min_size, 'Inf' if max_size is None else max_size,
                before_size_filter, len(results_df), below_min, above_max,
            )

            logger.info("ChEA3 API analysis completed: %d TFs", len(results_df))
        except Exception as e:
            logger.error(f"ChEA3 API analysis failed: {e}")
            return 1
    else:
        # Local database analysis mode
        if online_mode and database != 'chea3':
            logger.warning("--online supports ChEA3 only; other TF databases will run locally")

        # Load the TF database.
        db_manager = DatabaseManager(database_dir, species)
        try:
            tf_database = _load_tf_database(db_manager, database, species)
        except ValueError as exc:
            logger.error(str(exc))
            return 1

        if tf_database is None:
            logger.error(
                "Cannot load %s for species %s. Build or download the database before retrying.",
                database,
                species,
            )
            return 1

        # Run the selected analysis method.
        analyzer = TFEnrichmentAnalyzer(tf_database=tf_database)
        min_size, max_size = _tf_size_limits(args, method)
        common_options = {
            'library': getattr(args, 'tf_library', None),
            'tissue': getattr(args, 'tf_tissue', None),
            'regulation': getattr(args, 'tf_regulation', 'all'),
            'min_size': min_size,
            'max_size': max_size,
        }
        if method == 'ora' and background_file:
            try:
                common_options['background_genes'] = set(
                    _load_tf_gene_file(Path(background_file))
                )
            except (OSError, ValueError) as exc:
                logger.error("Could not load the TF ORA background file: %s", exc)
                return 1
        try:
            if method == 'ora':
                results_df = analyzer.ora(gene_set=gene_list, **common_options)
            elif method == 'gsea':
                results_df = analyzer.gsea(ranked_genes=ranked_genes, **common_options)
            else:
                logger.error(f"Unsupported analytical methods: {method}")
                return 1
        except ValueError as exc:
            logger.error(str(exc))
            return 1

    if results_df.empty:
        logger.warning("No transcription factors passed the recorded significance filters.")
        return 0

    # Preserve the official fgsea columns while adding stable TF identifiers and names.
    is_fgsea_output = set(('pathway', 'pval', 'padj', 'size', 'leadingEdge')).issubset(results_df.columns)
    if is_fgsea_output:
        metadata_frame = analyzer.metadata_frame() if 'analyzer' in locals() else pd.DataFrame()
        metadata = (
            metadata_frame.set_index('Term_ID')
            if isinstance(metadata_frame, pd.DataFrame) and 'Term_ID' in metadata_frame.columns
            else pd.DataFrame()
        )
        term_ids = results_df['pathway'].astype(str)
        term_data = {
            term_id: {
                'name': metadata.at[term_id, 'Term_Name']
                if term_id in metadata.index and 'Term_Name' in metadata.columns
                else term_id,
            }
            for term_id in term_ids
        }
        results_df = add_result_term_metadata(results_df, term_data)
        metadata_position = results_df.columns.get_loc('Term_Name') + 1
        for column, fallback in (
            ('TF', results_df['Term_ID']),
            ('Library', database),
            ('Context', ''),
            ('Evidence_Type', ''),
            ('Inference_Type', 'direct'),
        ):
            values = (
                results_df['Term_ID'].map(metadata[column])
                if column in metadata.columns else fallback
            )
            if column not in results_df.columns:
                results_df.insert(metadata_position, column, values)
                metadata_position += 1

        display_results = pd.DataFrame({
            'Term_ID': results_df['Term_ID'],
            'Term_Name': results_df['Term_Name'],
            'TF': results_df['TF'],
            'Library': results_df['Library'],
            'Context': results_df['Context'],
            'Overlap': results_df['size'],
            'Pvalue': results_df['pval'],
            'FDR': results_df['padj'],
            'Mode': 'unknown',
            'Overlap_Genes': results_df['leadingEdge'],
            'NES': results_df['NES'],
        })
        display_results = display_results.sort_values(
            ['Pvalue', 'Term_ID']
        ).reset_index(drop=True)
    else:
        display_results = results_df

    # Save CSV results
    csv_path = output_dir / f"tf_enrichment_{database}_{species}.csv"
    results_df.to_csv(csv_path, index=False)
    logger.info(f"Results saved: {csv_path}")
    if 'analyzer' in locals():
        analyzer.metadata_frame().to_csv(
            output_dir / f"tf_term_metadata_{database}_{species}.tsv",
            sep='\t', index=False,
        )
    tf_combine = getattr(args, 'tf_combine', 'none')
    if tf_combine != 'none':
        consensus = _tf_rank_consensus(display_results, tf_combine)
        consensus.to_csv(
            output_dir / f"tf_enrichment_consensus_{database}_{species}.csv",
            index=False,
        )

    # Print Top N Results
    top_results = display_results.head(top_n)
    print(f"\n{'='*60}")
    print(f"TF Enrichment Analysis Results (Top {top_n})")
    print(f"Database: {database.upper()}, Species: {species}")
    print(f"{'='*60}")
    print(f"{'Term':<44} {'Overlap':<10} {'Pvalue':<12} {'FDR':<12} {'Mode':<10}")
    print("-" * 90)
    for _, row in top_results.iterrows():
        mode = row.get('Mode', 'unknown')
        term_label = str(row.get('Term_Name') or row.get('TF', ''))
        term_label = term_label if len(term_label) <= 43 else term_label[:40] + '...'
        print(
            f"{term_label:<44} {row['Overlap']:<10} {row['Pvalue']:<12.2e} "
            f"{row['FDR']:<12.2e} {mode:<10}"
        )
    print(f"{'='*60}\n")

    # Generate interactive TF figures.
    try:
        viz = Visualizer()

        # Create the TF enrichment bar chart.
        fig_bar = viz.plot_tf_enrichment_bar(display_results, top_n=args.top_n)
        fig_bar.write_html(str(output_dir / "tf_enrichment_bar.html"))
        logger.info(f"Bar chart saved: {output_dir / 'tf_enrichment_bar.html'}")
        try:
            fig_bar.write_image(str(output_dir / "tf_enrichment_bar.png"))
        except Exception as e:
            if "kaleido" in str(e).lower():
                logger.info("Kaleido is not installed; retained the interactive TF bar chart and skipped PNG export")
            else:
                logger.warning(f"TF bar chart PNG export failed, HTML retained: {e}")

        # The regulatory-mode pie chart is available for ORA-style TF results.
        if not is_fgsea_output:
            fig_pie = viz.plot_tf_mode_pie(display_results)
            fig_pie.write_html(str(output_dir / "tf_mode_distribution.html"))
            logger.info("TF regulatory-mode chart saved: %s", output_dir / 'tf_mode_distribution.html')
            try:
                fig_pie.write_image(str(output_dir / "tf_mode_distribution.png"))
            except Exception as e:
                if "kaleido" in str(e).lower():
                    logger.info("Kaleido is not installed; retained the interactive TF pie chart and skipped PNG export")
                else:
                    logger.warning("TF pie chart PNG export failed; retained HTML output: %s", e)
    except Exception as e:
        logger.warning("TF figure generation failed; result tables remain valid: %s", e)

    # Generate HTML reports (optional)
    if args.report:
        try:
            _generate_tf_enrichment_report(
                display_results,
                output_dir,
                database,
                species,
                top_n,
                "gsea" if is_fgsea_output else "ora",
            )
            logger.info(f"HTML report generated: {output_dir / 'tf_enrichment_report.html'}")
        except Exception as e:
            logger.warning("Failed to generate the TF HTML report: %s", e)

    return 0


def _generate_tf_enrichment_report(
    results_df, output_dir: Path, database: str, species: str, top_n: int, method: str = "ora"
) -> None:
    """Generate the standalone TF enrichment HTML report."""
    # Read Jinja2 template
    template_path = Path(__file__).parent / "report" / "templates" / "tf_report.html"
    with open(template_path, 'r', encoding='utf-8') as f:
        template = Template(f.read())

    # Read interactive chart HTML content
    bar_chart_html = ""
    pie_chart_html = ""

    bar_html_path = output_dir / "tf_enrichment_bar.html"
    pie_html_path = output_dir / "tf_mode_distribution.html"

    if bar_html_path.exists():
        with open(bar_html_path, 'r', encoding='utf-8') as f:
            bar_chart_html = f.read()

    if pie_html_path.exists():
        with open(pie_html_path, 'r', encoding='utf-8') as f:
            pie_chart_html = f.read()

    # Prepare template variables
    significant_count = len(results_df[results_df['FDR'] < 0.05]) if 'FDR' in results_df.columns else 0
    gene_count = len(results_df) if not results_df.empty else 0

    # Convert DataFrame as a dictionary list
    top_results = results_df.head(top_n)
    results_list = []
    for _, row in top_results.iterrows():
        overlap_genes = str(row.get('Overlap_Genes', ''))
        results_list.append({
            'TF': row['TF'],
            'Mode': row.get('Mode', 'unknown'),
            'Overlap': row['Overlap'],
            'Pvalue': row['Pvalue'],
            'FDR': row['FDR'],
            'Overlap_Genes': overlap_genes[:100] + ('...' if len(overlap_genes) > 100 else '')
        })

    # Render Template
    html_content = template.render(
        db_name=database.upper(),
        species=species,
        gene_count=gene_count,
        significant_count=significant_count,
        method="fgseaMultilevel" if method == "gsea" else "Hypergeometric Test",
        pie_chart_html=pie_chart_html if pie_chart_html else '<p>Pie chart not available.</p>',
        bar_chart_html=bar_chart_html if bar_chart_html else '<p>Bar chart not available.</p>',
        results=results_list,
        timestamp=datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    )

    # Save HTML file
    report_path = output_dir / "tf_enrichment_report.html"
    with open(report_path, 'w', encoding='utf-8') as f:
        f.write(html_content)


if __name__ == '__main__':
    sys.exit(main())
