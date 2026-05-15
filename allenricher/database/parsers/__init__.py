"""
数据库解析器模块

包含各种数据库的解析器，用于将原始数据文件转换为 AllEnricher 标准格式。
"""

from .go import GOParser
from .kegg import KEGGParser
from .reactome import ReactomeParser
from .do import DOParser
from .disgenet import DisGeNETParser

__all__ = ['GOParser', 'KEGGParser', 'ReactomeParser', 'DOParser', 'DisGeNETParser']
