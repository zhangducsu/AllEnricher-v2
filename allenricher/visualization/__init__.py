"""Public visualization interfaces for AllEnricher."""

from allenricher.visualization.barplot import plot_barplot
from allenricher.visualization.gsea_plots import (
    plot_gsea_barplot,
    plot_gsea_enrichment,
    plot_gsea_lollipop,
    plot_gsea_multi_enrichment,
    plot_gsea_ridgeplot,
)
from allenricher.visualization.plot_theme import PlotTheme
from allenricher.visualization.plotter import Plotter

__all__ = [
    "Plotter",
    "PlotTheme",
    "plot_gsea_barplot",
    "plot_gsea_enrichment",
    "plot_gsea_lollipop",
    "plot_gsea_multi_enrichment",
    "plot_gsea_ridgeplot",
    "plot_barplot",
]
