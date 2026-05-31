"""
Database management for AllEnricher v2.0

数据库管理模块，负责加载和管理各种富集分析数据库。
支持从 v1 的数据库文件格式（.tab.gz）加载。
"""

import gzip
import csv
import json
import logging
from pathlib import Path
from typing import Dict, List, Optional, Set

import pandas as pd

logger = logging.getLogger(__name__)


# KEGG 物种代码到 NCBI Gene taxid 的映射
# 注意：gene_info.gz 使用 NCBI Gene taxid，而非 GOA 的 taxid
KEGG_CODE_TO_TAXID: Dict[str, int] = {
    # 常见模式生物
    'hsa': 9606,   # Homo sapiens (Human)
    'hpy': 835,    # Helicobacter pylori 26695
    'mta': 4530,   # Oryza sativa (Rice)
    'ath': 3702,   # Arabidopsis thaliana
    'bta': 9913,   # Bos taurus (Bovine)
    'cel': 6239,   # Caenorhabditis elegans
    'cfa': 9615,   # Canis familiaris (Dog)
    'dre': 7955,   # Danio rerio (Zebrafish)
    'dme': 7227,   # Drosophila melanogaster (Fruit fly)
    'gga': 9031,   # Gallus gallus (Chicken)
    'mcf': 594,    # Mycobacterium tuberculosis CDC1551
    'mmu': 10090,  # Mus musculus (Mouse)
    'rno': 10116,  # Rattus norvegicus (Rat)
    'sce': 4932,   # Saccharomyces cerevisiae (Yeast)
    'spo': 4896,   # Schizosaccharomyces pombe (Fission yeast)
    'xla': 8355,   # Xenopus laevis (African clawed frog)
    'xtr': 8364,   # Xenopus tropicalis
    # 其他常见物种
    'eco': 562,    # Escherichia coli K-12
    'bsu': 224308, # Bacillus subtilis 168
    'pae': 208964, # Pseudomonas aeruginosa PAO1
    'syf': 1148,   # Synechococcus elongatus PCC 7942
    'syn': 1140,   # Synechocystis sp. PCC 6803
    'mtu': 83332,  # Mycobacterium tuberculosis H37Rv
}


class DatabaseManager:
    """数据库管理器

    负责加载和管理富集分析所需的各种数据库（GO、KEGG、Reactome等）。
    兼容 v1 的数据库文件格式。

    Attributes:
        database_dir: 数据库文件所在目录
        species: 物种代码（如 hsa、mmu）
        databases: 已加载的数据库字典
        term_names: Term ID 到名称的映射字典
    """

    def __init__(self, database_dir: str, species: str):
        """初始化数据库管理器

        Args:
            database_dir: 数据库文件目录路径
            species: 物种代码（如 hsa）
        """
        self.database_dir = Path(database_dir)
        self.species = species
        self.databases: Dict[str, Dict] = {}
        self.term_names: Dict[str, Dict[str, str]] = {}  # {db_name: {term_id: term_name}}
        self._active_version: Optional[str] = None

    def _find_species_dir(self, database_dir: Path, species: str, version: Optional[str] = None) -> Path:
        """自动查找物种数据库目录

        支持两种目录结构:
        1. v2 格式: database/organism/v{date}/{species}/  (如 database/organism/2024-01-01/hsa/)
        2. v1 格式: database/  (直接在该目录下查找文件)

        Args:
            database_dir: 基础数据库目录
            species: 物种代码
            version: 指定使用的数据库版本（如 v20260515），为 None 时自动使用最新版本

        Returns:
            实际的物种数据库目录路径
        """
        # 模式 1: database/organism/v{date}/{species}/ (v2 结构)
        organism_dir = database_dir / "organism"

        # 如果指定了版本，直接使用
        if version and organism_dir.exists():
            species_dir = organism_dir / version / species
            if species_dir.exists():
                self._active_version = version
                return species_dir
            # 版本不存在时列出可用版本
            available = sorted(
                d.name for d in organism_dir.iterdir()
                if d.is_dir() and (d / species).exists()
            )
            if available:
                logger.error("版本 '%s' 的物种 '%s' 不存在。可用版本: %s", version, species, ", ".join(available))
            else:
                logger.error("物种 '%s' 没有任何已构建的版本。", species)

        # 自动查找最新版本
        if organism_dir.exists():
            for version_dir in sorted(organism_dir.iterdir(), reverse=True):
                if version_dir.is_dir():
                    species_dir = version_dir / species
                    if species_dir.exists():
                        self._active_version = version_dir.name
                        return species_dir

        # v1 兼容
        if (database_dir / f"{species}.GO2gene.tab.gz").exists():
            self._active_version = "v1-legacy"
            return database_dir

        self._active_version = None
        return database_dir

    @property
    def active_version(self) -> Optional[str]:
        """当前活跃的数据库版本号

        Returns:
            版本号字符串（如 'v20260515'、'v1-legacy'），未加载时为 None
        """
        return self._active_version

    def get_build_metadata(self) -> Optional[Dict]:
        """获取当前活跃版本的构建元数据

        从 build_manifest.json 读取构建时记录的元信息，
        包含构建时间、依赖版本等。

        Returns:
            元数据字典，文件不存在或未加载版本时返回 None
        """
        if not self._active_version:
            return None
        manifest_path = Path(self.database_dir) / "organism" / self._active_version / self.species / "build_manifest.json"
        if not manifest_path.exists():
            return None
        try:
            with open(manifest_path, "r", encoding="utf-8") as f:
                return json.load(f)
        except Exception as e:
            logger.warning("读取 build_manifest.json 失败: %s", e)
            return None

    def load_databases(self, database_names: List[str], version: Optional[str] = None) -> None:
        """加载指定的数据库

        Args:
            database_names: 数据库名称列表，如 ["GO", "KEGG"]
            version: 指定使用的数据库版本，为 None 时自动使用最新版本
        """
        for name in database_names:
            self.load_database(name, version=version)

    def load_database(self, name: str, version: Optional[str] = None) -> None:
        """加载单个数据库

        从 v1 格式的 .tab.gz 文件加载数据库。
        文件格式：第一行是表头（Gene\tTermID1\tTermID2...），
        后续行是基因和所属条目的 0/1 矩阵。

        Args:
            name: 数据库名称（如 "GO"、"KEGG"）
            version: 指定使用的数据库版本，为 None 时自动使用最新版本
        """
        # 自动查找物种目录（支持 v2 的 organism/v{date}/{species}/ 结构）
        self.database_dir = self._find_species_dir(self.database_dir, self.species, version=version)

        # 先加载 Term 名称映射（从 .tab.id.gz 或 .2disc.gz）
        self._load_term_names(name)

        # 构建文件路径（v1 格式：{species}.{name}2gene.tab.gz）
        # 注意部分数据库有特殊文件名映射
        name_to_prefix = {
            'GO': 'GO',
            'KEGG': 'kegg',
            'REACTOME': 'Reactome',
            'DO': 'DO',
            'DISGENET': 'CUI',  # DisGeNET 使用 CUI 前缀
            'WIKIPATHWAYS': 'WikiPathways',
            'TRRUST': 'TF2target',  # TRRUST 使用 TF2target 前缀
            'CHEA3': 'ChEA3_2gene',  # ChEA3 使用 ChEA3_2gene 前缀
            'ANIMALTFDB': 'AnimalTFDB_2gene',
            'HTFTARGET': 'hTF_2gene',
        }
        prefix = name_to_prefix.get(name.upper(), name)

        filename = f"{self.species}.{prefix}2gene.tab.gz"
        filepath = self.database_dir / filename

        if not filepath.exists():
            # 尝试小写
            filename = f"{self.species}.{prefix.lower()}2gene.tab.gz"
            filepath = self.database_dir / filename

        if not filepath.exists():
            # 也尝试数据库名称本身
            filename = f"{self.species}.{name}2gene.tab.gz"
            filepath = self.database_dir / filename

        if not filepath.exists():
            filename = f"{self.species}.{name.lower()}2gene.tab.gz"
            filepath = self.database_dir / filename

        if not filepath.exists():
            raise FileNotFoundError(f"数据库文件不存在: {self.database_dir} 下未找到 {self.species}.{prefix}2gene.tab.gz")

        # 解析数据库文件（此时名称会使用 term_names 映射）
        term_data = self._parse_tab_file(filepath, name)
        self.databases[name] = term_data

    @staticmethod
    def _capitalize(text: str) -> str:
        """首字母大写，同时保留全大写词（如 DNA、RNA、mRNA、ATP）

        Args:
            text: 输入文本

        Returns:
            首字母大写的文本
        """
        # 保留全大写词的映射（如 DNA -> DNA）
        upper_map = {
            'DNA': 'DNA', 'RNA': 'RNA', 'mRNA': 'mRNA', 'tRNA': 'tRNA',
            'rRNA': 'rRNA', 'ATP': 'ATP', 'ADP': 'ADP', 'GTP': 'GTP',
            'NAD': 'NAD', 'NADH': 'NADH', 'FAD': 'FAD', 'CoA': 'CoA',
            'AMP': 'AMP', 'GMP': 'GMP', 'UMP': 'UMP', 'CMP': 'CMP',
            'MAPK': 'MAPK', 'PI3K': 'PI3K', 'AKT': 'AKT', 'EGF': 'EGF',
            'TNF': 'TNF', 'IL': 'IL', 'IFN': 'IFN', 'TGF': 'TGF',
            'VEGF': 'VEGF', 'PDGF': 'PDGF', 'FGF': 'FGF', 'IGF': 'IGF',
            'JAK': 'JAK', 'STAT': 'STAT', 'NF': 'NF', 'AP': 'AP',
            'HIF': 'HIF', 'PPAR': 'PPAR', 'RXR': 'RXR', 'LXR': 'LXR',
            'FXR': 'FXR', 'CAR': 'CAR', 'PXR': 'PXR', 'SHP': 'SHP',
            'SREBP': 'SREBP', 'PGC': 'PGC', 'cAMP': 'cAMP', 'cGMP': 'cGMP',
            'PKA': 'PKA', 'PKC': 'PKC', 'PKG': 'PKG', 'AMPk': 'AMPK',
        }
        words = text.split()
        result = []
        for word in words:
            upper_word = word.upper()
            if upper_word in upper_map:
                result.append(upper_map[upper_word])
            else:
                result.append(word.capitalize())
        return ' '.join(result)

    def _format_term_name(self, db_name: str, raw_name: str) -> str:
        """格式化 Term 名称，统一为 "Category|Name" 格式，首字母大写

        Args:
            db_name: 数据库名称 (GO, KEGG, Reactome, DO)
            raw_name: 原始名称字符串

        Returns:
            格式化后的名称字符串
        """
        # 将下划线替换为空格（用于 KEGG 的 pathway_name_with_underscores）
        name = raw_name.replace('_', ' ')

        if db_name.upper() == 'GO':
            # GO 格式: "biological_process:mitochondrion inheritance"
            # 转为: "Biological Process|Mitochondrion Inheritance"
            if ':' in name:
                namespace, term = name.split(':', 1)
                return f"{namespace.title()}|{self._capitalize(term)}"
            return self._capitalize(name)

        elif db_name.upper() == 'KEGG':
            # KEGG 格式可能是:
            #   "Category|SubCategory|PathwayName" (有分类，三层)
            #   "Uncategorized|Uncategorized|PathwayName" (无分类，三层但无用)
            #   "PathwayName" (只有名称)
            # 根据实际层级数决定输出格式
            if '|' in name:
                parts = name.split('|')
                pathway_name = self._capitalize(parts[-1])

                # 如果一级分类是 Uncategorized，只显示通路名（1级）
                if parts[0].lower() == 'uncategorized':
                    return pathway_name

                # 有有效分类，显示 Category|SubCategory|PathwayName（三层）
                if len(parts) >= 3:
                    cat = self._capitalize(parts[0])
                    subcat = self._capitalize(parts[1])
                    return f"{cat}|{subcat}|{pathway_name}"
                elif len(parts) == 2:
                    return f"{parts[0].title()}|{pathway_name}"
            return self._capitalize(name)

        # Reactome, DO 等无层级结构，首字母大写即可
        return self._capitalize(name)

    def _load_term_names(self, db_name: str) -> None:
        """加载 Term ID 到名称的映射

        v1 数据库中，*.tab.id.gz 或 *.2disc.gz 文件包含 Term 名称。
        例如 GO 有 hsa.GO.tab.id.gz，KEGG 有 hsa.KEGG2disc.gz

        Args:
            db_name: 数据库名称
        """
        self.term_names[db_name] = {}

        # 数据库名到文件前缀的映射（与 load_database 一致）
        name_to_prefix = {
            'GO': 'GO',
            'KEGG': 'kegg',
            'REACTOME': 'Reactome',
            'DO': 'DO',
            'DISGENET': 'CUI',
            'WIKIPATHWAYS': 'WikiPathways',
        }
        prefix = name_to_prefix.get(db_name.upper(), db_name)

        # 方法1: 尝试加载 {species}.{db}.tab.id.gz 文件
        id_file = self.database_dir / f"{self.species}.{db_name}.tab.id.gz"
        if not id_file.exists():
            id_file = self.database_dir / f"{self.species}.{db_name.lower()}.tab.id.gz"

        if id_file.exists():
            try:
                with gzip.open(id_file, 'rt', encoding='utf-8') as f:
                    for line in f:
                        line = line.strip()
                        if not line:
                            continue
                        # 格式: TermID\tTermName\tParentTerms
                        parts = line.split('\t')
                        if len(parts) >= 2:
                            self.term_names[db_name][parts[0]] = parts[1]
                return  # 如果成功加载，直接返回
            except Exception as e:
                print(f"警告: 加载 {id_file} 失败: {e}")

        # 方法2: 尝试从 *.2disc.gz 加载（使用 prefix 构建路径，与 load_database 一致）
        disc_file = self.database_dir / f"{self.species}.{prefix}2disc.gz"
        if not disc_file.exists():
            disc_file = self.database_dir / f"{self.species}.{prefix.lower()}2disc.gz"
        if not disc_file.exists():
            # 回退：尝试使用 db_name 本身
            disc_file = self.database_dir / f"{self.species}.{db_name}2disc.gz"
        if not disc_file.exists():
            disc_file = self.database_dir / f"{self.species}.{db_name.lower()}2disc.gz"
        if not disc_file.exists():
            # 回退：尝试不带物种前缀的文件（如 GO2disc.gz）
            disc_file = self.database_dir / f"{prefix}2disc.gz"
        if not disc_file.exists():
            disc_file = self.database_dir / f"{prefix.lower()}2disc.gz"

        if disc_file.exists():
            try:
                count = 0
                with gzip.open(disc_file, 'rt', encoding='utf-8') as f:
                    for line in f:
                        line = line.strip()
                        if not line:
                            continue
                        # 格式: TermID\tTermName (TermName 使用 | 分隔层级)
                        parts = line.split('\t')
                        if len(parts) >= 2:
                            term_id = parts[0]
                            raw_name = parts[1]
                            term_name = self._format_term_name(db_name, raw_name)
                            self.term_names[db_name][term_id] = term_name
                            count += 1
                if count > 0:
                    print(f"    从 {disc_file.name} 加载了 {count} 个 Term 名称")
            except Exception as e:
                print(f"警告: 加载 {disc_file} 失败: {e}")

    def _parse_tab_file(self, filepath: Path, db_name: str) -> Dict[str, Dict]:
        """解析 .tab.gz 文件

        Args:
            filepath: 数据库文件路径
            db_name: 数据库名称

        Returns:
            Dict: {term_id: {"name": str, "genes": List[str]}}
        """
        term_data: Dict[str, Dict] = {}

        # 获取该数据库的名称映射
        name_map = self.term_names.get(db_name, {})

        with gzip.open(filepath, 'rt', encoding='utf-8') as f:
            reader = csv.reader(f, delimiter='\t')

            # 读取表头
            header = next(reader)
            # 表头格式：Gene\tTermID1\tTermID2\t...
            term_ids = header[1:]  # 第一列是 Gene，后面是 Term ID

            # 初始化 term_data，使用 term_names 中的名称或 Term ID 作为名称
            for term_id in term_ids:
                # 如果有名称映射，使用映射的名称；否则使用 Term ID
                term_name = name_map.get(term_id, term_id)
                term_data[term_id] = {"name": term_name, "genes": []}

            # 读取数据行
            for row in reader:
                if len(row) < 2:
                    continue

                gene = row[0]
                # 后面的列表示该基因是否属于对应条目（1=属于，0=不属于）
                for i, value in enumerate(row[1:]):
                    if i < len(term_ids) and value == '1':
                        term_id = term_ids[i]
                        term_data[term_id]["genes"].append(gene)

        return term_data

    def get_all_term_data(self) -> Dict[str, Dict]:
        """获取所有已加载数据库的条目数据

        Returns:
            Dict: {db_name: {term_id: {"name": str, "genes": List[str]}}}
        """
        return self.databases

    def get_background_genes(self) -> Set[str]:
        """获取背景基因集（所有数据库中基因的并集）

        Returns:
            Set[str]: 背景基因集合
        """
        background = set()
        for db_data in self.databases.values():
            for term_info in db_data.values():
                background.update(term_info["genes"])
        return background

    def get_database_genes(self, db_name: str) -> Set[str]:
        """获取指定数据库的所有基因

        Args:
            db_name: 数据库名称

        Returns:
            Set[str]: 该数据库中所有基因的集合
        """
        if db_name not in self.databases:
            return set()

        genes = set()
        for term_info in self.databases[db_name].values():
            genes.update(term_info["genes"])
        return genes

    def get_genome_genes(self, taxid: Optional[int] = None, species_code: Optional[str] = None) -> Set[str]:
        """获取该物种的全基因组基因（来自 gene_info.gz）

        gene_info.gz 包含该物种的所有基因（无论是否有 GO 注释），
        与 gene2go.gz（只有 GO 注释的基因）不同。

        Args:
            taxid: NCBI Gene 分类学 ID，如果为 None 则从 species_code 推断
            species_code: KEGG 物种代码（如 'hsa'），用于从 KEGG_CODE_TO_TAXID 推断 taxid

        Returns:
            Set[str]: 全基因组 Gene Symbol 集合（第3列，与 GO/KEGG 注释标识符一致）
        """
        # 如果提供了 species_code，尝试从中推断 taxid
        if taxid is None and species_code is not None:
            taxid = KEGG_CODE_TO_TAXID.get(species_code.lower())

        if taxid is None:
            # 如果仍然没有 taxid，返回空集合
            return set()

        genome_file = self.database_dir / "gene_info.gz"
        if not genome_file.exists():
            # 尝试在 basic/go 目录下查找
            # 向上查找 database 根目录
            db_root = self.database_dir
            for _ in range(5):
                if db_root.parent:
                    db_root = db_root.parent
                if db_root.name in ('database', 'AllEnricher-v2', 'AllEnricher'):
                    break

            basic_go_dir = db_root / "basic" / "go"
            if basic_go_dir.exists():
                # 使用最新版本
                versions = sorted([d for d in basic_go_dir.iterdir() if d.is_dir()], reverse=True)
                for version_dir in versions:
                    genome_file = version_dir / "gene_info.gz"
                    if genome_file.exists():
                        break
                    genome_file = self.database_dir / "gene_info.gz"  # 重置
                if not genome_file.exists():
                    return set()
            else:
                return set()

        genes = set()
        # 只包含真正的基因类型，排除 biological-region（基因组区域，非基因）、pseudo（假基因）
        _VALID_GENE_TYPES = frozenset({
            "protein-coding", "ncRNA", "snoRNA", "rRNA", "tRNA", "snRNA",
            "scRNA", "other",
        })
        with gzip.open(genome_file, "rt", encoding="utf-8") as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                parts = line.split("\t")
                if len(parts) >= 10:
                    file_taxid = parts[0]
                    gene_symbol = parts[2]  # 第3列：Gene Symbol（如 A1BG）
                    gene_type = parts[9]    # 第10列：基因类型
                    # 过滤只保留指定 taxid 且是真正基因类型的条目
                    try:
                        if int(file_taxid) == taxid and gene_type in _VALID_GENE_TYPES:
                            if gene_symbol and gene_symbol != "-" and not gene_symbol.startswith("NEW|"):
                                genes.add(gene_symbol)
                    except ValueError:
                        continue
        return genes

    def load_trrust(self, species: Optional[str] = None) -> Optional[Dict[str, 'pd.DataFrame']]:
        """加载 TRRUST 转录因子-靶基因数据库

        查找并加载以下文件:
        - {species}.TF2target.tab.gz: TF 到靶基因的映射
        - {species}.gene2TF.tab.gz: 基因到 TF 的映射
        - {species}.TF2disc.gz: TF 描述信息

        Args:
            species: 物种代码，为 None 时使用实例的 self.species

        Returns:
            包含三个 DataFrame 的字典 {'tf2target': DataFrame, 'gene2tf': DataFrame, 'tf_info': DataFrame}，
            如果任一必要文件不存在则返回 None
        """
        sp = species or self.species
        base_dir = self._find_species_dir(self.database_dir, sp)

        tf2target_file = base_dir / f"{sp}.TF2target.tab.gz"
        gene2tf_file = base_dir / f"{sp}.gene2TF.tab.gz"
        tf2disc_file = base_dir / f"{sp}.TF2disc.gz"

        # 检查必要文件是否存在
        if not tf2target_file.exists() and not gene2tf_file.exists():
            logger.warning("TRRUST 数据库文件不存在: %s", base_dir)
            return None

        result: Dict[str, 'pd.DataFrame'] = {}

        if tf2target_file.exists():
            result['tf2target'] = pd.read_csv(tf2target_file, sep='\t', compression='gzip')
        else:
            result['tf2target'] = pd.DataFrame()

        if gene2tf_file.exists():
            result['gene2tf'] = pd.read_csv(gene2tf_file, sep='\t', compression='gzip')
        else:
            result['gene2tf'] = pd.DataFrame()

        if tf2disc_file.exists():
            result['tf_info'] = pd.read_csv(tf2disc_file, sep='\t', compression='gzip')
        else:
            result['tf_info'] = pd.DataFrame()

        return result

    def load_chea3(self, species: Optional[str] = None) -> Optional[Dict[str, 'pd.DataFrame']]:
        """加载 ChEA3 转录因子-靶基因数据库

        查找并加载以下文件:
        - {species}.ChEA3_2gene.tab.gz: 基因到 TF 的映射
        - {species}.ChEA3_2disc.gz: TF 描述信息

        Args:
            species: 物种代码，为 None 时使用实例的 self.species

        Returns:
            包含两个 DataFrame 的字典 {'gene2tf': DataFrame, 'tf_info': DataFrame}，
            如果必要文件不存在则返回 None
        """
        sp = species or self.species
        base_dir = self._find_species_dir(self.database_dir, sp)

        gene2tf_file = base_dir / f"{sp}.ChEA3_2gene.tab.gz"
        tf2disc_file = base_dir / f"{sp}.ChEA3_2disc.gz"

        # 检查必要文件是否存在
        if not gene2tf_file.exists():
            logger.warning("ChEA3 数据库文件不存在: %s", base_dir)
            return None

        result: Dict[str, 'pd.DataFrame'] = {}

        result['gene2tf'] = pd.read_csv(gene2tf_file, sep='\t', compression='gzip')

        if tf2disc_file.exists():
            result['tf_info'] = pd.read_csv(tf2disc_file, sep='\t', compression='gzip')
        else:
            result['tf_info'] = pd.DataFrame()

        return result

    def load_htftarget(self, species: Optional[str] = None) -> Optional[Dict[str, pd.DataFrame]]:
        """加载 hTFtarget 数据库

        Returns:
            {'gene2tf': DataFrame, 'tf_info': DataFrame}
        """
        sp = species or self.species
        db_dir = self._find_species_db_dir(sp)
        if db_dir is None:
            return None

        gene2tf_file = db_dir / f"{sp}.hTF_2gene.tab.gz"
        disc_file = db_dir / f"{sp}.hTF_2disc.gz"

        if not gene2tf_file.exists():
            return None

        result = {}
        result['gene2tf'] = pd.read_csv(gene2tf_file, sep='\t', compression='gzip', low_memory=False)

        if disc_file.exists():
            result['tf_info'] = pd.read_csv(disc_file, sep='\t', compression='gzip')

        return result

    def load_animaltfdb(self, species: Optional[str] = None) -> Optional[Dict[str, pd.DataFrame]]:
        """加载 AnimalTFDB 数据库（同源映射结果）

        Returns:
            {'gene2tf': DataFrame, 'tf_info': DataFrame}
        """
        sp = species or self.species
        db_dir = self._find_species_db_dir(sp)
        if db_dir is None:
            return None

        gene2tf_file = db_dir / f"{sp}.AnimalTFDB_2gene.tab.gz"
        disc_file = db_dir / f"{sp}.AnimalTFDB_mapped_2disc.gz"

        if not gene2tf_file.exists():
            return self.load_htftarget(sp)

        result = {}
        result['gene2tf'] = pd.read_csv(gene2tf_file, sep='\t', compression='gzip', low_memory=False)

        if disc_file.exists():
            result['tf_info'] = pd.read_csv(disc_file, sep='\t', compression='gzip')

        return result
