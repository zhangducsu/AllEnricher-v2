"""Core module - contains configuration management and enrichment analysis engines"""
from allenricher.core.config import Config                    # Configure Management Category
from allenricher.core.enrichment import EnrichmentAnalyzer    # Fuzzy Analysis Engine

__all__ = ["Config", "EnrichmentAnalyzer"]
