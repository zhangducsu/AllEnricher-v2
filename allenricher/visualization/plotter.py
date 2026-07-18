"""Coordinate the supported Python enrichment figures."""

import logging
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import pandas as pd

logger = logging.getLogger(__name__)


class Plotter:
    """Generate the supported Python figures for enrichment result tables."""

    def __init__(self, output_dir: str, config=None):
        """Initialize the plotting facade with an optional analysis configuration."""
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.config = config

    def _prepare_barplot_data(self, data: pd.DataFrame) -> pd.DataFrame:
        """Normalize an enrichment result table for the shared bar-plot function."""
        df = data.copy()
        result = pd.DataFrame()

        # Detect the canonical fgsea result schema.
        is_gsea = 'NES' in df.columns or 'nes' in df.columns

        # Normalize ORA and GSEA column names for the shared plotting functions.
        column_mappings = {
            'term_id': ['Term_ID', 'term_id', 'ID', 'id', 'pathway'],
            'term': ['Term_Name', 'term_name', 'Term', 'term', 'Description', 'description', 'pathway'],
            'qvalue': ['Adjusted_P_Value', 'adjusted_p_value', 'adjP', 'qvalue', 'Qvalue', 'FDR', 'fdr', 'padj'],
            'gene_count': ['Gene_Count', 'gene_count', 'setSize', 'ObservedGeneNum', 'Count', 'count', 'size'],
            'hierarchy': ['Hierarchy', 'hierarchy'],
        }
        if not is_gsea:
            column_mappings['rich_factor'] = ['Rich_Factor', 'rich_factor', 'RichFactor', 'richfactor']

        for std_col, possible_cols in column_mappings.items():
            for col in possible_cols:
                if col in df.columns:
                    result[std_col] = df[col]
                    break

        if is_gsea:
            # GSEA uses signed NES rather than an ORA enrichment factor.
            for col in ['NES', 'nes', 'ES', 'es']:
                if col in df.columns:
                    result['rich_factor'] = df[col]
                    break
        else:
            # Derive the ORA enrichment factor from recorded ratios when needed.
            if 'rich_factor' not in result.columns and 'gene_count' in result.columns:
                bg_col = None
                for col in ['Background_Count', 'background_count', 'TermGeneNum']:
                    if col in df.columns:
                        bg_col = df[col]
                        break
                if bg_col is not None:
                    result['rich_factor'] = result['gene_count'] / bg_col

        return result

    def plot_barplot(
        self,
        data: pd.DataFrame,
        database: str,
        output_file: str,
        top_n: int = 20,
        style: Optional[str] = None,
        palette: Optional[str] = None,
        figsize: Optional[Tuple[float, float]] = None,
        hierarchy_map: Optional[Dict[str, str]] = None,
    ) -> str:
        """Plot significant terms as compact horizontal bars."""
        # Resolve the requested output format from configuration.
        output_path = self.output_dir / output_file

        from .barplot import plot_barplot as _plot_barplot

        # Normalize the table without mutating analysis results.
        plot_data = self._prepare_barplot_data(data)

        # Use the run-level DPI setting for raster output.
        dpi = 300
        if self.config and hasattr(self.config, 'plot_dpi'):
            dpi = self.config.plot_dpi

        try:
            if 'NES' in data.columns or 'nes' in data.columns:
                from .gsea_plots import plot_gsea_barplot

                plot_gsea_barplot(
                    results_df=data,
                    database=database,
                    top_n=top_n,
                    style=style or 'nature',
                    palette=palette,
                    output_file=str(output_path),
                    dpi=dpi,
                    figsize=figsize,
                )
            else:
                _plot_barplot(
                    data=plot_data,
                    database=database,
                    top_n=top_n,
                    style=style or 'nature',
                    palette=palette,
                    output_file=str(output_path),
                    dpi=dpi,
                    figsize=figsize,
                    hierarchy_map=hierarchy_map,
                )
        except Exception as e:
            logger.error("Failed to generate the enrichment bar plot: %s", e)
            raise

        return str(output_path)

    def plot_lollipop(
        self,
        data: pd.DataFrame,
        database: str,
        output_file: str,
        top_n: int = 20,
        style: Optional[str] = None,
        palette: Optional[str] = None,
        figsize: Optional[Tuple[float, float]] = None,
    ) -> str:
        """Plot ORA terms as enrichment-factor lollipops."""
        output_path = self.output_dir / output_file

        from .gsea_plots import plot_gsea_lollipop

        plot_data = self._prepare_barplot_data(data)
        dpi = 300
        if self.config and hasattr(self.config, 'plot_dpi'):
            dpi = self.config.plot_dpi

        try:
            import matplotlib.pyplot as plt
            fig = plot_gsea_lollipop(
                plot_data,
                top_n=top_n,
                title=f"{database} Over-Representation Analysis (ORA)",
                output_file=str(output_path),
                style=style or 'nature',
                palette=palette,
                dpi=dpi,
                figsize=figsize,
            )
            plt.close(fig)
        except Exception as e:
            logger.error("Failed to generate the enrichment lollipop plot: %s", e)
            raise

        return str(output_path)

    def plot_all(
        self,
        data: pd.DataFrame,
        database: str,
        top_n: int = 20,
        style: Optional[str] = None,
        palette: Optional[str] = None,
        hierarchy_map: Optional[Dict[str, str]] = None,
    ) -> Dict[str, str]:
        """Generate every supported standard figure for one result table."""
        plots = {}

        formats = ['png']
        if self.config and getattr(self.config, 'plot_formats', None):
            formats = list(dict.fromkeys(fmt.lower() for fmt in self.config.plot_formats))
        invalid = [fmt for fmt in formats if fmt not in {'png', 'pdf', 'svg'}]
        if invalid:
            raise ValueError(f"Unsupported figure format: {invalid}")

        figsize = None
        if (self.config and getattr(self.config, 'plot_width', None) is not None
                and getattr(self.config, 'plot_height', None) is not None):
            figsize = (float(self.config.plot_width), float(self.config.plot_height))

        for fmt in formats:
            bar_path = self.plot_barplot(
                data, database, f"{database}_barplot.{fmt}", top_n,
                style=style, palette=palette, figsize=figsize,
                hierarchy_map=hierarchy_map,
            )
            plots["barplot"] = bar_path

            lollipop_path = self.plot_lollipop(
                data, database, f"{database}_lollipop.{fmt}", top_n,
                style=style, palette=palette, figsize=figsize,
            )
            plots["lollipop"] = lollipop_path

        return plots
