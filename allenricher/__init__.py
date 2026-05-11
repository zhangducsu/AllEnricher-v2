"""
AllEnricher v2.0 - 综合性基因集功能富集分析工具

一个现代化、高性能的基因集富集分析工具，支持多种物种、数据库和算法，
并提供AI驱动的结果解读功能。

主要功能：
- 多种富集算法（Fisher精确检验、超几何检验、GSEA、ssGSEA）
- 多种数据库（GO、KEGG、Reactome、WikiPathways、MSigDB、DO、DisGeNET）
- 交互式HTML报告生成
- REST API接口
- AI驱动的结果解读
- 并行处理支持
"""

__version__ = "2.0.0"       # 版本号
__author__ = "AllEnricher Team"  # 作者
__license__ = "MIT"              # 开源许可证

# 导出核心类，供外部直接使用：from allenricher import EnrichmentAnalyzer, Config, ...
from allenricher.core.enrichment import EnrichmentAnalyzer    # 富集分析引擎
from allenricher.core.config import Config                    # 配置管理
from allenricher.database.manager import DatabaseManager      # 数据库管理器
from allenricher.visualization.plotter import Plotter          # 可视化绘图
from allenricher.report.generator import ReportGenerator      # HTML报告生成器

__all__ = [
    "EnrichmentAnalyzer",   # 富集分析引擎
    "Config",               # 配置类
    "DatabaseManager",      # 数据库管理器
    "Plotter",              # 可视化绘图器
    "ReportGenerator",      # 报告生成器
]
