"""
层级注释文件解析器

解析用户提供的条目注释文件（非GMT文件），支持层级结构。

支持三种TSV格式：
1. 四列（带层级）: gene<TAB>term_id<TAB>term_name<TAB>hierarchy
2. 三列: gene<TAB>term_id<TAB>term_name
3. 两列: gene<TAB>term（term_name 同时作为 term_id）
"""

import gzip
from pathlib import Path
from typing import Dict, List, Optional, Set


class AnnotationRecord:
    """单个注释记录"""

    def __init__(self, gene: str, term_id: str, term_name: str,
                 hierarchy: Optional[str] = None):
        self.gene = gene
        self.term_id = term_id
        self.term_name = term_name
        self.hierarchy = hierarchy

    @property
    def hierarchy_levels(self) -> List[str]:
        """返回层级列表

        Returns:
            层级路径拆分后的列表，如 ["Biological Process", "Cellular Process"]
            如果没有层级信息则返回空列表
        """
        if not self.hierarchy:
            return []
        return self.hierarchy.split('|')

    def __repr__(self) -> str:
        return (
            f"AnnotationRecord(gene={self.gene!r}, term_id={self.term_id!r}, "
            f"term_name={self.term_name!r}, hierarchy={self.hierarchy!r})"
        )


class AnnotationParser:
    """注释文件解析器

    解析用户提供的条目注释文件，支持自动格式检测和层级结构提取。

    支持的文件格式（tab分隔）：
    - 2列: gene<TAB>term（term_name 同时作为 term_id）
    - 3列: gene<TAB>term_id<TAB>term_name
    - 4列: gene<TAB>term_id<TAB>term_name<TAB>hierarchy

    Args:
        file_path: 注释文件路径
        format_type: 文件格式类型，可选值: 'four_column', 'three_column', 'two_column'
            如果为 None 则自动检测
        hierarchy_separator: 层级路径分隔符，默认为 '|'
    """

    def __init__(self, file_path: Optional[str] = None, format_type: Optional[str] = None,
                 hierarchy_separator: str = '|', filepath: Optional[str] = None):
        # 兼容 file_path 和 filepath 两种参数名
        _path = filepath or file_path
        if _path is None:
            raise ValueError("必须提供 file_path 或 filepath 参数")
        self.file_path = Path(_path)
        self.format_type = format_type
        self.hierarchy_separator = hierarchy_separator
        self._records: Optional[List[AnnotationRecord]] = None

    def _detect_format(self) -> str:
        """自动检测文件格式

        根据第一个有效数据行的列数判断格式：
        - 4列及以上 → four_column
        - 3列 → three_column
        - 2列 → two_column

        Returns:
            格式类型字符串

        Raises:
            FileNotFoundError: 文件不存在时抛出
            ValueError: 无法检测格式时抛出
        """
        if not self.file_path.exists():
            raise FileNotFoundError(f"注释文件不存在: {self.file_path}")

        with open(self.file_path, 'r', encoding='utf-8') as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith('#'):
                    continue
                col_count = len(line.split('\t'))
                if col_count >= 4:
                    return 'four_column'
                elif col_count == 3:
                    return 'three_column'
                elif col_count == 2:
                    return 'two_column'
                else:
                    raise ValueError(
                        f"无法识别的文件格式：第一行有 {col_count} 列，"
                        f"期望 2、3 或 4 列"
                    )

        raise ValueError("注释文件为空或仅包含注释行")

    def parse(self) -> List[AnnotationRecord]:
        """解析注释文件

        Returns:
            AnnotationRecord 列表

        Raises:
            FileNotFoundError: 文件不存在时抛出
            ValueError: 文件格式无法识别时抛出
        """
        if self._records is not None:
            return self._records

        if not self.file_path.exists():
            raise FileNotFoundError(f"注释文件不存在: {self.file_path}")

        fmt = self.format_type or self._detect_format()
        self._records = []

        opener = gzip.open if str(self.file_path).endswith('.gz') else open
        with opener(self.file_path, 'rt', encoding='utf-8') as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith('#'):
                    continue

                parts = line.split('\t')

                if fmt == 'four_column' and len(parts) >= 4:
                    gene = parts[0]
                    term_id = parts[1]
                    term_name = parts[2]
                    hierarchy = parts[3]
                    self._records.append(
                        AnnotationRecord(gene, term_id, term_name, hierarchy)
                    )
                elif fmt == 'three_column' and len(parts) >= 3:
                    gene = parts[0]
                    term_id = parts[1]
                    term_name = parts[2]
                    self._records.append(
                        AnnotationRecord(gene, term_id, term_name)
                    )
                elif fmt == 'two_column' and len(parts) >= 2:
                    gene = parts[0]
                    term = parts[1]
                    self._records.append(
                        AnnotationRecord(gene, term, term)
                    )

        return self._records

    def get_term_genes(self) -> Dict[str, Set[str]]:
        """获取 term_id 到基因集合的映射

        Returns:
            字典，key 为 term_id，value 为关联的基因符号集合
        """
        records = self.parse()
        term_genes: Dict[str, Set[str]] = {}
        for rec in records:
            if rec.term_id not in term_genes:
                term_genes[rec.term_id] = set()
            term_genes[rec.term_id].add(rec.gene)
        return term_genes

    def get_term_names(self) -> Dict[str, str]:
        """获取 term_id 到 term_name 的映射

        Returns:
            字典，key 为 term_id，value 为 term_name
        """
        records = self.parse()
        term_names: Dict[str, str] = {}
        for rec in records:
            if rec.term_id not in term_names:
                term_names[rec.term_id] = rec.term_name
        return term_names

    def get_term_hierarchies(self) -> Dict[str, str]:
        """获取 term_id 到层级字符串的映射

        Returns:
            字典，key 为 term_id，value 为层级字符串；
            如果没有层级信息则不包含该 term_id
        """
        records = self.parse()
        term_hierarchies: Dict[str, str] = {}
        for rec in records:
            if rec.hierarchy and rec.term_id not in term_hierarchies:
                term_hierarchies[rec.term_id] = rec.hierarchy
        return term_hierarchies

    def get_hierarchy_tree(self) -> Dict:
        """获取层级树结构

        根据所有记录的层级信息构建树形结构。
        对于没有层级信息的记录，term_id 作为顶层节点。

        Returns:
            嵌套字典表示的层级树，格式为:
            {
                "level1_name": {
                    "level2_name": {
                        "term_id": {"genes": set, "term_name": str},
                        ...
                    },
                    ...
                },
                ...
            }
            对于没有层级的 term，直接放在顶层:
            {
                "term_id": {"genes": set, "term_name": str},
                ...
            }
        """
        records = self.parse()
        tree: Dict = {}

        # 先按 term_id 聚合基因
        term_data: Dict[str, Dict] = {}
        for rec in records:
            if rec.term_id not in term_data:
                term_data[rec.term_id] = {
                    'genes': set(),
                    'term_name': rec.term_name,
                }
            term_data[rec.term_id]['genes'].add(rec.gene)

        # 构建层级树
        for rec in records:
            levels = rec.hierarchy_levels
            if not levels:
                # 没有层级信息，直接放在顶层
                if rec.term_id not in tree:
                    tree[rec.term_id] = term_data[rec.term_id]
                continue

            # 逐层构建树
            current = tree
            for level in levels:
                if level not in current:
                    current[level] = {}
                current = current[level]

            # 在最底层放置 term 数据
            if rec.term_id not in current:
                current[rec.term_id] = term_data[rec.term_id]

        return tree
