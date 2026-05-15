"""
数据库构建模块

提供数据下载、解析和数据库构建功能，用于生成 AllEnricher 所需的各种
功能注释数据库文件（GO、KEGG、Reactome、DO、DisGeNET）。
"""

from .downloader import DataDownloader
from .builder import DatabaseBuilder

__all__ = ['DataDownloader', 'DatabaseBuilder']
