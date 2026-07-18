"""Reporting module - Generate interactive HTML-enrichment analysis, including statistical summaries, data tables and graphs"""
from allenricher.report.generator import ReportGenerator  # HTML Report Generator
from allenricher.report.visualizer import Visualizer  # TF Fusion Analytic Visualizer

__all__ = ["ReportGenerator", "Visualizer"]
