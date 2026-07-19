"""AllEnricher: gene-set enrichment analysis across species and databases.

The package provides over-representation analysis, GSEA, ssGSEA, and GSVA;
standard pathway, disease, transcription-factor, and custom gene-set libraries;
HTML reports; a REST API; and optional AI-assisted result interpretation.
"""

__version__ = "2.1.0"
__author__ = "AllEnricher Team"
__license__ = "MIT"

# Export the primary public API at the package root.
from allenricher.core.enrichment import EnrichmentAnalyzer    # Fuzzy Analysis Engine
from allenricher.core.config import Config                    # Configuration management
from allenricher.database.manager import DatabaseManager      # Database Manager
from allenricher.visualization.plotter import Plotter          # Visual Drawing
from allenricher.report.generator import ReportGenerator      # HTML Report Generator

__all__ = [
    "EnrichmentAnalyzer",   # Fuzzy Analysis Engine
    "Config",               # Configure Classes
    "DatabaseManager",      # Database Manager
    "Plotter",              # Visualise Drawing
    "ReportGenerator",      # Report Generator
]
