"""核心模块 - 包含配置管理和富集分析引擎"""
from allenricher.core.config import Config                    # 配置管理类
from allenricher.core.enrichment import EnrichmentAnalyzer    # 富集分析引擎

__all__ = ["Config", "EnrichmentAnalyzer"]
