"""
统一物种注册表管理模块

提供物种信息的注册、查询、持久化等功能，支持 GO、KEGG、Reactome、DO
等多个数据库的物种覆盖状态追踪。数据以 TSV 格式存储于 supported_species.tsv。
"""

from __future__ import annotations

import csv
from dataclasses import dataclass, asdict
from pathlib import Path
from typing import Optional, Dict, List, Any


# TSV 列定义，保持与文件格式一致
_FIELD_NAMES: List[str] = [
    "taxid", "latin_name", "common_name",
    "has_go", "go_source", "go_filename", "go_file_size", "go_gene_count", "go_term_count",
    "has_kegg", "kegg_code", "kegg_code_source", "kegg_gene_count", "kegg_pathway_count",
    "has_reactome", "reactome_code", "reactome_gene_count", "reactome_pathway_count",
    "has_do", "do_gene_count", "do_term_count",
]


@dataclass
class SpeciesEntry:
    """物种注册条目数据类

    存储单个物种的完整信息，包括分类学标识、常用名，
    以及各数据库（GO、KEGG、Reactome、DO）的覆盖状态与统计信息。

    Attributes:
        taxid: NCBI 分类学 ID（主键）
        latin_name: 物种拉丁学名
        common_name: 物种常用名
        has_go: 是否有 GO 注释数据
        go_source: GO 数据来源
        go_filename: GO 数据文件名
        go_file_size: GO 数据文件大小（字节）
        go_gene_count: GO 注释基因数量
        go_term_count: GO 术语数量
        has_kegg: 是否有 KEGG 数据
        kegg_code: KEGG 物种缩写代码
        kegg_code_source: KEGG 代码来源
        kegg_gene_count: KEGG 基因数量
        kegg_pathway_count: KEGG 通路数量
        has_reactome: 是否有 Reactome 数据
        reactome_code: Reactome 物种代码
        reactome_gene_count: Reactome 基因数量
        reactome_pathway_count: Reactome 通路数量
        has_do: 是否有 DO（疾病本体）数据
        do_gene_count: DO 关联基因数量
        do_term_count: DO 术语数量
    """
    taxid: int
    latin_name: str
    common_name: Optional[str] = None
    # GO 相关字段
    has_go: bool = False
    go_source: Optional[str] = None
    go_filename: Optional[str] = None
    go_file_size: Optional[int] = None
    go_gene_count: Optional[int] = None
    go_term_count: Optional[int] = None
    # KEGG 相关字段
    has_kegg: bool = False
    kegg_code: Optional[str] = None
    kegg_code_source: Optional[str] = None
    kegg_gene_count: Optional[int] = None
    kegg_pathway_count: Optional[int] = None
    # Reactome 相关字段
    has_reactome: bool = False
    reactome_code: Optional[str] = None
    reactome_gene_count: Optional[int] = None
    reactome_pathway_count: Optional[int] = None
    # DO 相关字段
    has_do: bool = False
    do_gene_count: Optional[int] = None
    do_term_count: Optional[int] = None


class SpeciesRegistry:
    """统一物种注册表

    管理所有已支持物种的信息，提供按 TaxID、拉丁名、KEGG 代码等多种方式
    的查询功能，并支持按数据库覆盖状态进行过滤。数据以 TSV 格式持久化。

    Attributes:
        registry_path: 注册表文件路径（supported_species.tsv）
        entries: 以 taxid 为键的物种条目字典
    """

    def __init__(self, registry_path: Path) -> None:
        """初始化物种注册表

        Args:
            registry_path: 注册表 TSV 文件路径
        """
        self.registry_path = Path(registry_path)
        self.entries: Dict[int, SpeciesEntry] = {}

    def load(self) -> None:
        """从 TSV 文件加载物种注册表数据

        读取 registry_path 指向的 supported_species.tsv 文件，
        解析每一行为 SpeciesEntry 对象。布尔值识别 True/False，
        缺失值（- 或空）映射为 None。
        """
        self.entries.clear()

        if not self.registry_path.exists():
            return

        with open(self.registry_path, "r", encoding="utf-8") as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            for row in reader:
                entry = self._parse_row(row)
                if entry is not None:
                    self.entries[entry.taxid] = entry

    def save(self) -> None:
        """将当前注册表数据保存到 TSV 文件

        将 entries 中所有条目按 taxid 升序写入 registry_path，
        布尔值输出为 True/False，缺失值输出为 -。
        """
        self.registry_path.parent.mkdir(parents=True, exist_ok=True)

        field_names = _FIELD_NAMES

        with open(self.registry_path, "w", encoding="utf-8", newline="") as fh:
            writer = csv.DictWriter(fh, fieldnames=field_names, delimiter="\t")
            writer.writeheader()

            for taxid in sorted(self.entries):
                writer.writerow(self._format_row(self.entries[taxid]))

    def add_entry(self, entry: SpeciesEntry) -> None:
        """添加或更新物种条目

        Args:
            entry: 要添加的物种条目，若 taxid 已存在则覆盖
        """
        self.entries[entry.taxid] = entry

    def query_by_taxid(self, taxid: int) -> Optional[SpeciesEntry]:
        """通过 NCBI TaxID 精确查询物种

        Args:
            taxid: NCBI 分类学 ID

        Returns:
            对应的 SpeciesEntry，未找到则返回 None
        """
        return self.entries.get(taxid)

    def query_by_latin_name(self, name: str) -> List[SpeciesEntry]:
        """通过拉丁学名模糊查询物种（大小写不敏感）

        Args:
            name: 查询关键词，与 latin_name 进行子串匹配

        Returns:
            所有匹配的 SpeciesEntry 列表
        """
        keyword = name.strip().lower()
        if not keyword:
            return []
        return [
            entry for entry in self.entries.values()
            if keyword in entry.latin_name.lower()
        ]

    def query_by_kegg_code(self, code: str) -> Optional[SpeciesEntry]:
        """通过 KEGG 物种代码精确查询

        Args:
            code: KEGG 物种缩写代码（如 'hsa'）

        Returns:
            对应的 SpeciesEntry，未找到则返回 None
        """
        code_normalized = code.strip().lower()
        for entry in self.entries.values():
            if entry.has_kegg and entry.kegg_code is not None and entry.kegg_code.lower() == code_normalized:
                return entry
        return None

    def filter_by_databases(
        self,
        go: Optional[bool] = None,
        kegg: Optional[bool] = None,
        reactome: Optional[bool] = None,
        do: Optional[bool] = None,
    ) -> List[SpeciesEntry]:
        """按数据库覆盖状态过滤物种

        各参数为 None 时表示不限制该条件，为 True/False 时要求对应字段严格匹配。

        Args:
            go: 是否要求有 GO 数据
            kegg: 是否要求有 KEGG 数据
            reactome: 是否要求有 Reactome 数据
            do: 是否要求有 DO 数据

        Returns:
            满足所有过滤条件的 SpeciesEntry 列表
        """
        results: List[SpeciesEntry] = []
        for entry in self.entries.values():
            if go is not None and entry.has_go != go:
                continue
            if kegg is not None and entry.has_kegg != kegg:
                continue
            if reactome is not None and entry.has_reactome != reactome:
                continue
            if do is not None and entry.has_do != do:
                continue
            results.append(entry)
        return results

    def get_summary(self) -> Dict[str, Any]:
        """获取各数据库覆盖统计汇总

        Returns:
            包含以下键的字典：
            - total_species: 总物种数
            - go: {count, with_gene_count, with_term_count}
            - kegg: {count, with_gene_count, with_pathway_count}
            - reactome: {count, with_gene_count, with_pathway_count}
            - do: {count, with_gene_count, with_term_count}
        """
        go_count = 0
        go_with_genes = 0
        go_with_terms = 0
        kegg_count = 0
        kegg_with_genes = 0
        kegg_with_pathways = 0
        reactome_count = 0
        reactome_with_genes = 0
        reactome_with_pathways = 0
        do_count = 0
        do_with_genes = 0
        do_with_terms = 0

        for entry in self.entries.values():
            if entry.has_go:
                go_count += 1
                if entry.go_gene_count is not None:
                    go_with_genes += 1
                if entry.go_term_count is not None:
                    go_with_terms += 1
            if entry.has_kegg:
                kegg_count += 1
                if entry.kegg_gene_count is not None:
                    kegg_with_genes += 1
                if entry.kegg_pathway_count is not None:
                    kegg_with_pathways += 1
            if entry.has_reactome:
                reactome_count += 1
                if entry.reactome_gene_count is not None:
                    reactome_with_genes += 1
                if entry.reactome_pathway_count is not None:
                    reactome_with_pathways += 1
            if entry.has_do:
                do_count += 1
                if entry.do_gene_count is not None:
                    do_with_genes += 1
                if entry.do_term_count is not None:
                    do_with_terms += 1

        return {
            "total_species": len(self.entries),
            "go": {
                "count": go_count,
                "with_gene_count": go_with_genes,
                "with_term_count": go_with_terms,
            },
            "kegg": {
                "count": kegg_count,
                "with_gene_count": kegg_with_genes,
                "with_pathway_count": kegg_with_pathways,
            },
            "reactome": {
                "count": reactome_count,
                "with_gene_count": reactome_with_genes,
                "with_pathway_count": reactome_with_pathways,
            },
            "do": {
                "count": do_count,
                "with_gene_count": do_with_genes,
                "with_term_count": do_with_terms,
            },
        }

    def get_species_detail(self, taxid: int) -> Optional[Dict[str, Any]]:
        """获取物种详细信息字典

        Args:
            taxid: NCBI 分类学 ID

        Returns:
            包含物种所有字段的字典，未找到则返回 None
        """
        entry = self.entries.get(taxid)
        if entry is None:
            return None
        return asdict(entry)

    @staticmethod
    def generate_kegg_abbreviation(latin_name: str) -> str:
        """根据拉丁学名生成 KEGG 物种缩写代码

        规则：属名首字母（小写）+ 种名前 2 个字母（小写）。

        Examples:
            "Homo sapiens" -> "hsa"
            "Mus musculus" -> "mmu"
            "Arabidopsis thaliana" -> "ath"

        Args:
            latin_name: 物种拉丁学名

        Returns:
            3 字符的 KEGG 缩写代码
        """
        parts = latin_name.strip().split()
        if len(parts) < 2:
            # 仅有属名的情况，取前 3 个字母
            genus = parts[0].lower()[:3]
            return genus.ljust(3, "x")[:3]
        genus_initial = parts[0].lower()[0]
        species_prefix = parts[1].lower()[:2]
        return genus_initial + species_prefix

    @classmethod
    def load_default(cls, root_dir: str = "./database") -> "SpeciesRegistry":
        """从默认路径加载物种注册表

        在 root_dir 目录下查找 supported_species.tsv 文件并加载。

        Args:
            root_dir: 数据库根目录路径，默认为 "./database"

        Returns:
            已加载数据的 SpeciesRegistry 实例
        """
        registry = cls(registry_path=Path(root_dir) / "supported_species.tsv")
        registry.load()
        return registry

    # ------------------------------------------------------------------
    # 内部辅助方法
    # ------------------------------------------------------------------

    @staticmethod
    def _parse_bool(value: str) -> bool:
        """将字符串解析为布尔值

        Args:
            value: 字符串值，期望 "True" 或 "False"

        Returns:
            对应的布尔值
        """
        return value.strip() == "True"

    @staticmethod
    def _parse_optional_int(value: str) -> Optional[int]:
        """将字符串解析为可选整数

        Args:
            value: 字符串值，"-" 或空字符串映射为 None

        Returns:
            整数值或 None
        """
        stripped = value.strip()
        if stripped in ("", "-"):
            return None
        return int(stripped)

    @staticmethod
    def _parse_optional_str(value: str) -> Optional[str]:
        """将字符串解析为可选字符串

        Args:
            value: 字符串值，"-" 或空字符串映射为 None

        Returns:
            字符串值或 None
        """
        stripped = value.strip()
        if stripped in ("", "-"):
            return None
        return stripped

    def _parse_row(self, row: Dict[str, str]) -> Optional[SpeciesEntry]:
        """将 TSV 行字典解析为 SpeciesEntry

        Args:
            row: 从 csv.DictReader 读取的行字典

        Returns:
            解析后的 SpeciesEntry，解析失败返回 None
        """
        try:
            taxid = int(row["taxid"].strip())
        except (ValueError, KeyError):
            return None

        return SpeciesEntry(
            taxid=taxid,
            latin_name=row["latin_name"].strip(),
            common_name=self._parse_optional_str(row.get("common_name", "")),
            has_go=self._parse_bool(row.get("has_go", "False")),
            go_source=self._parse_optional_str(row.get("go_source", "")),
            go_filename=self._parse_optional_str(row.get("go_filename", "")),
            go_file_size=self._parse_optional_int(row.get("go_file_size", "")),
            go_gene_count=self._parse_optional_int(row.get("go_gene_count", "")),
            go_term_count=self._parse_optional_int(row.get("go_term_count", "")),
            has_kegg=self._parse_bool(row.get("has_kegg", "False")),
            kegg_code=self._parse_optional_str(row.get("kegg_code", "")),
            kegg_code_source=self._parse_optional_str(row.get("kegg_code_source", "")),
            kegg_gene_count=self._parse_optional_int(row.get("kegg_gene_count", "")),
            kegg_pathway_count=self._parse_optional_int(row.get("kegg_pathway_count", "")),
            has_reactome=self._parse_bool(row.get("has_reactome", "False")),
            reactome_code=self._parse_optional_str(row.get("reactome_code", "")),
            reactome_gene_count=self._parse_optional_int(row.get("reactome_gene_count", "")),
            reactome_pathway_count=self._parse_optional_int(row.get("reactome_pathway_count", "")),
            has_do=self._parse_bool(row.get("has_do", "False")),
            do_gene_count=self._parse_optional_int(row.get("do_gene_count", "")),
            do_term_count=self._parse_optional_int(row.get("do_term_count", "")),
        )

    @staticmethod
    def _format_optional(value: Optional[Any]) -> str:
        """将可选值格式化为 TSV 字符串

        Args:
            value: 布尔值、整数、字符串或 None

        Returns:
            格式化后的字符串，None 输出为 "-"
        """
        if value is None:
            return "-"
        if isinstance(value, bool):
            return "True" if value else "False"
        return str(value)

    def _format_row(self, entry: SpeciesEntry) -> Dict[str, str]:
        """将 SpeciesEntry 格式化为 TSV 行字典

        Args:
            entry: 物种条目

        Returns:
            键为列名、值为格式化字符串的字典
        """
        entry_dict = asdict(entry)
        return {
            key: self._format_optional(entry_dict[key])
            for key in _FIELD_NAMES
        }
