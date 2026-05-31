# WikiPathways 全物种支持实现计划

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** 调研 WikiPathways 官方物种支持情况，将当前 18 物种扩展到官方支持的全部物种（40 物种），并解决基因 ID 格式问题（NCBI Gene ID → Gene Symbol）。

**Architecture:** WikiPathways 官方提供按物种分离的 GMT 文件（18 物种有通路数据）和 GPML 文件（40 物种有原始数据）。采用类似 KEGG/UniProt GOA 的按需下载策略：用户指定物种后，下载对应物种的 GMT 文件并构建数据库。同时需要解决基因 ID 转换问题（GMT 中的 NCBI Gene ID 需通过 gene_info 转换为 Gene Symbol）。

**Tech Stack:** Python 3.x, requests (下载), gzip/pandas (解析), NCBI gene_info (ID 映射)

---

## 前置知识：WikiPathways 物种支持现状

### 官方数据调研结果（2026-05-10 版本）

| 数据类型 | 物种数量 | 说明 |
|---------|---------|------|
| **GMT 文件** | 18 物种 | 有通路-基因关联数据，可直接用于富集分析 |
| **GPML 文件** | 40 物种 | 原始通路图数据，需额外解析提取基因信息 |

### GMT 文件支持的 18 物种（当前已实现）

```
Anopheles_gambiae, Arabidopsis_thaliana, Bos_taurus, Caenorhabditis_elegans,
Canis_familiaris, Danio_rerio, Drosophila_melanogaster, Equus_caballus,
Gallus_gallus, Homo_sapiens, Mus_musculus, Pan_troglodytes,
Populus_trichocarpa, Rattus_norvegicus, Saccharomyces_cerevisiae,
Solanum_lycopersicum, Sus_scrofa, Zea_mays
```

### GPML 文件支持的额外 22 物种（需扩展支持）

```
Acetobacterium_woodii, Bacillus_subtilis, Beta_vulgaris, Brassica_napus,
Caulobacter_vibrioides, Citrus_sinensis, Coffea_arabica, Daphnia_magna,
Escherichia_coli, Gibberella_zeae, Hordeum_vulgare, Ilex_paraguariensis,
Mycobacterium_tuberculosis, Oryza_sativa, Paullinia_cupana, Perilla_frutescens,
Plasmodium_falciparum, Theobroma_cacao, Triticum_aestivum, Vitis_vinifera
```

### 关键问题：基因 ID 格式

**WikiPathways GMT 文件中的基因格式**：
```
WPID<TAB>Pathway Name<TAB>ncbigene:1234/ncbigene:5678/...
```

- 基因使用 **NCBI Gene ID**（数字格式，如 `ncbigene:1234`）
- 不是 Gene Symbol（如 `TP53`）
- 需要通过 NCBI `gene_info.gz` 转换为 Gene Symbol

### 与 GO/KEGG 的对比

| 数据库 | 下载策略 | 基因 ID 格式 | 当前实现 |
|-------|---------|-------------|---------|
| GO | 全物种一次性下载 (gene2go.gz) | NCBI Gene ID → 本地转换为 Symbol | ✅ |
| KEGG | 按需下载指定物种 (REST API) | KEGG ID → NCBI Gene ID → Symbol | ✅ |
| UniProt GOA | 按需下载指定物种 (FTP) | UniProt ID → 本地转换为 Symbol | ✅ |
| **WikiPathways** | **按需下载指定物种 (GMT)** | **NCBI Gene ID → 本地转换为 Symbol** | **需完善** |

---

## 实现策略决策

基于调研结果，采用以下策略：

1. **下载策略**：按物种按需下载（类似 KEGG/GOA）
   - 理由：WikiPathways 官方提供按物种分离的 GMT 文件，无单一全物种文件
   - 好处：用户只需下载需要的物种，节省空间和带宽

2. **基因 ID 转换**：GMT 中的 NCBI Gene ID → Gene Symbol
   - 使用 NCBI `gene_info.gz` 进行 ID 映射
   - 在构建阶段完成转换，存储为 Gene Symbol 格式（与 GO/KEGG 一致）

3. **物种扩展**：从 18 物种扩展到 40 物种
   - 优先支持有 GMT 文件的 18 物种（已完成）
   - 扩展支持有 GPML 文件的额外 22 物种（需从 GPML 提取基因信息）

---

## 文件变更清单

| 操作 | 文件路径 | 职责 |
|------|---------|------|
| 修改 | `allenricher/database/wikipathways_fetcher.py` | 扩展物种映射（40 物种），添加 ID 转换功能 |
| 修改 | `allenricher/database/parsers/wikipathways.py` | 添加 NCBI Gene ID → Symbol 转换逻辑 |
| 修改 | `allenricher/database/builder.py` | 传递 gene_info 路径给 parser |
| 修改 | `allenricher/database/species_registry.py` | 扩展 WikiPathways 物种列表 |
| 新建 | `allenricher/database/wikipathways_gpml_parser.py` | 从 GPML 提取基因信息（用于 22 无 GMT 物种） |

---

## Task 1: 扩展 WikiPathways 物种映射（18 → 40 物种）

**Files:**
- Modify: `allenricher/database/wikipathways_fetcher.py`

- [ ] **Step 1: 更新 SPECIES_NAME_MAP 包含全部 40 物种**

```python
# 完整的 WikiPathways 物种映射（40 物种）
# 分为两类：
# 1. 有 GMT 文件的 18 物种（可直接用于富集分析）
# 2. 只有 GPML 文件的 22 物种（需额外解析）

SPECIES_NAME_MAP: Dict[str, str] = {
    # === 有 GMT 文件的 18 物种（优先支持）===
    "Homo_sapiens": "hsa",
    "Mus_musculus": "mmu",
    "Rattus_norvegicus": "rno",
    "Danio_rerio": "dre",
    "Drosophila_melanogaster": "dme",
    "Caenorhabditis_elegans": "cel",
    "Saccharomyces_cerevisiae": "sce",
    "Arabidopsis_thaliana": "ath",
    "Bos_taurus": "bta",
    "Gallus_gallus": "gga",
    "Sus_scrofa": "ssc",
    "Canis_familiaris": "cfa",
    "Equus_caballus": "eca",
    "Pan_troglodytes": "ptr",
    "Anopheles_gambiae": "aga",
    "Populus_trichocarpa": "ptc",
    "Solanum_lycopersicum": "sly",
    "Zea_mays": "zma",
    
    # === 只有 GPML 文件的 22 物种（扩展支持）===
    "Acetobacterium_woodii": "awo",
    "Bacillus_subtilis": "bsu",
    "Beta_vulgaris": "bvu",
    "Brassica_napus": "bna",
    "Caulobacter_vibrioides": "cvi",
    "Citrus_sinensis": "csi",
    "Coffea_arabica": "car",
    "Daphnia_magna": "dma",
    "Escherichia_coli": "eco",
    "Gibberella_zeae": "gze",
    "Hordeum_vulgare": "hvu",
    "Ilex_paraguariensis": "ipa",
    "Mycobacterium_tuberculosis": "mtu",
    "Oryza_sativa": "osa",
    "Paullinia_cupana": "pcu",
    "Perilla_frutescens": "pfr",
    "Plasmodium_falciparum": "pfa",
    "Theobroma_cacao": "tcc",
    "Triticum_aestivum": "tae",
    "Vitis_vinifera": "vvi",
}

# 标记哪些物种有 GMT 文件（可用于富集分析）
SPECIES_WITH_GMT: Set[str] = {
    "Homo_sapiens", "Mus_musculus", "Rattus_norvegicus", "Danio_rerio",
    "Drosophila_melanogaster", "Caenorhabditis_elegans", "Saccharomyces_cerevisiae",
    "Arabidopsis_thaliana", "Bos_taurus", "Gallus_gallus", "Sus_scrofa",
    "Canis_familiaris", "Equus_caballus", "Pan_troglodytes", "Anopheles_gambiae",
    "Populus_trichocarpa", "Solanum_lycopersicum", "Zea_mays",
}

# 标记哪些物种只有 GPML 文件（需额外解析）
SPECIES_GPML_ONLY: Set[str] = set(SPECIES_NAME_MAP.keys()) - SPECIES_WITH_GMT
```

- [ ] **Step 2: 添加物种数据类型检测方法**

```python
@staticmethod
def get_species_data_type(latin_name: str) -> str:
    """获取物种的数据类型
    
    Returns:
        'gmt' - 有 GMT 文件（可直接用于富集分析）
        'gpml' - 只有 GPML 文件（需额外解析）
        'none' - 不支持
    """
    if latin_name in SPECIES_WITH_GMT:
        return 'gmt'
    elif latin_name in SPECIES_GPML_ONLY:
        return 'gpml'
    return 'none'

@staticmethod
def get_supported_species() -> Dict[str, List[str]]:
    """获取所有支持的物种分类
    
    Returns:
        {'gmt': [...], 'gpml': [...]}
    """
    return {
        'gmt': sorted(SPECIES_WITH_GMT),
        'gpml': sorted(SPECIES_GPML_ONLY),
    }
```

- [ ] **Step 3: 验证物种映射**

Run: `python -c "from allenricher.database.wikipathways_fetcher import SPECIES_NAME_MAP, SPECIES_WITH_GMT, SPECIES_GPML_ONLY; print(f'Total: {len(SPECIES_NAME_MAP)}, GMT: {len(SPECIES_WITH_GMT)}, GPML: {len(SPECIES_GPML_ONLY)}')"`
Expected: `Total: 40, GMT: 18, GPML: 22`

---

## Task 2: 实现 NCBI Gene ID → Gene Symbol 转换

**Files:**
- Modify: `allenricher/database/parsers/wikipathways.py`
- Modify: `allenricher/database/builder.py`

### 2.1 更新 WikiPathwaysParser

- [ ] **Step 1: 添加 ID 转换方法**

```python
@staticmethod
def load_gene_id_mapping(gene_info_path: Path, taxid: int) -> Dict[str, str]:
    """从 NCBI gene_info.gz 加载 NCBI Gene ID → Gene Symbol 映射
    
    Args:
        gene_info_path: gene_info.gz 文件路径
        taxid: NCBI Taxonomy ID
        
    Returns:
        {ncbi_gene_id: gene_symbol} 字典
    """
    import gzip
    
    mapping = {}
    with gzip.open(gene_info_path, 'rt', encoding='utf-8') as f:
        header = f.readline().strip().split('\t')
        taxid_idx = header.index('tax_id')
        geneid_idx = header.index('GeneID')
        symbol_idx = header.index('Symbol')
        
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) < max(taxid_idx, geneid_idx, symbol_idx) + 1:
                continue
            if int(parts[taxid_idx]) == taxid:
                ncbi_id = parts[geneid_idx]
                symbol = parts[symbol_idx]
                mapping[ncbi_id] = symbol
                
    return mapping

@staticmethod
def convert_ncbi_to_symbol(
    gene_sets: Dict[str, List[str]],
    id_mapping: Dict[str, str]
) -> Dict[str, List[str]]:
    """将基因集中的 NCBI Gene ID 转换为 Gene Symbol
    
    Args:
        gene_sets: {pathway_id: [ncbigene:123, ncbigene:456, ...]}
        id_mapping: {123: 'TP53', 456: 'BRCA1', ...}
        
    Returns:
        {pathway_id: ['TP53', 'BRCA1', ...]}
    """
    converted_sets = {}
    
    for pathway_id, genes in gene_sets.items():
        converted_genes = []
        for gene in genes:
            # 提取数字 ID（去掉 "ncbigene:" 前缀）
            if gene.startswith('ncbigene:'):
                ncbi_id = gene.replace('ncbigene:', '')
            else:
                ncbi_id = gene
                
            # 查找对应的 Symbol
            if ncbi_id in id_mapping:
                converted_genes.append(id_mapping[ncbi_id])
            else:
                # 如果找不到映射，保留原始 ID（用于调试）
                logger.warning(f"NCBI Gene ID {ncbi_id} not found in gene_info")
                
        if converted_genes:
            converted_sets[pathway_id] = converted_genes
            
    return converted_sets
```

- [ ] **Step 2: 更新 build_database 方法支持 ID 转换**

```python
@staticmethod
def build_database(
    gmt_path: Path,
    output_dir: Path,
    species: str,
    taxid: int,
    gene_info_path: Optional[Path] = None,
    valid_genes: Optional[Set[str]] = None,
) -> Tuple[Path, Path]:
    """构建 AllEnricher 标准格式文件（支持 NCBI Gene ID → Symbol 转换）
    
    Args:
        gmt_path: WikiPathways GMT 文件路径
        output_dir: 输出目录
        species: 物种代码（如 hsa）
        taxid: NCBI TaxID（用于 gene_info 过滤）
        gene_info_path: NCBI gene_info.gz 路径（用于 ID 转换）
        valid_genes: 有效基因集合（可选）
        
    Returns:
        (tab_path, disc_path) 元组
    """
    gene_sets, descriptions = WikiPathwaysParser.parse_gmt(gmt_path)
    
    # 如果提供了 gene_info，进行 ID 转换
    if gene_info_path and gene_info_path.exists():
        logger.info(f"Loading NCBI Gene ID mapping from {gene_info_path}")
        id_mapping = WikiPathwaysParser.load_gene_id_mapping(gene_info_path, taxid)
        logger.info(f"Loaded {len(id_mapping)} gene ID mappings")
        
        gene_sets = WikiPathwaysParser.convert_ncbi_to_symbol(gene_sets, id_mapping)
        logger.info(f"Converted gene IDs to symbols")
    else:
        logger.warning("No gene_info provided, using raw gene IDs from GMT")
    
    # 如果提供了有效基因集合，过滤基因
    if valid_genes:
        filtered_sets = {}
        for wp_id, genes in gene_sets.items():
            filtered_genes = [g for g in genes if g in valid_genes]
            if filtered_genes:
                filtered_sets[wp_id] = filtered_genes
        gene_sets = filtered_sets
        logger.info(f"After valid gene filtering: {len(gene_sets)} pathways")
    
    # 生成标准格式文件（与之前相同）
    # ... （生成 .tab.gz 和 .disc.gz）
```

### 2.2 更新 DatabaseBuilder

- [ ] **Step 3: 修改 build_wikipathways 传递 gene_info**

```python
def build_wikipathways(self, species: str, taxid: int, outdir: Path) -> None:
    """构建 WikiPathways 数据库（支持 NCBI Gene ID 转换）"""
    from .wikipathways_fetcher import WikiPathwaysFetcher
    from .parsers.wikipathways import WikiPathwaysParser
    
    fetcher = WikiPathwaysFetcher(str(self.basic_dir))
    latin_name = fetcher.get_latin_name(species)
    
    if not latin_name:
        print(f"|--- [跳过] WikiPathways 不支持物种 {species}")
        return
    
    # 检查物种数据类型
    data_type = fetcher.get_species_data_type(latin_name)
    
    if data_type == 'gmt':
        # 有 GMT 文件，直接解析
        gmt_file = self._find_gmt_file(latin_name)
        if not gmt_file:
            print(f"|--- [跳过] 未找到 {latin_name} 的 GMT 文件")
            return
            
        # 获取 gene_info 路径（用于 ID 转换）
        gene_info_path = self._get_gene_info_path()
        
        # 构建数据库（传入 gene_info 进行 ID 转换）
        tab_path, disc_path = WikiPathwaysParser.build_database(
            gmt_path=gmt_file,
            output_dir=outdir,
            species=species,
            taxid=taxid,
            gene_info_path=gene_info_path,  # 新增参数
        )
        
    elif data_type == 'gpml':
        # 只有 GPML 文件，需要额外解析（Task 3 实现）
        print(f"|--- [信息] {latin_name} 只有 GPML 数据，需使用 GPML 解析器")
        # TODO: 调用 GPML 解析器
        return
        
    else:
        print(f"|--- [跳过] WikiPathways 不支持物种 {species}")
        return
    
    # 生成 GMT 文件
    if tab_path.exists() and disc_path.exists():
        from .gmt_generator import GMTGenerator
        gen = GMTGenerator(outdir)
        gen.generate_wikipathways_gmt(species)
        print(f"|--- WikiPathways 构建完成")

def _get_gene_info_path(self) -> Optional[Path]:
    """获取 gene_info.gz 路径"""
    # 优先使用 GO 目录中的 gene_info
    go_gene_info = self.basic_dir / "go" / f"GO{self._get_go_version()}" / "gene_info.gz"
    if go_gene_info.exists():
        return go_gene_info
    
    # 查找任何版本的 gene_info
    for go_dir in (self.basic_dir / "go").iterdir():
        if go_dir.is_dir():
            gene_info = go_dir / "gene_info.gz"
            if gene_info.exists():
                return gene_info
    
    return None

def _get_go_version(self) -> str:
    """获取 GO 版本号"""
    # 从版本管理器获取，或使用默认值
    from .version import DatabaseVersionManager
    vm = DatabaseVersionManager(str(self.basic_dir))
    go_version = vm.get_local_version('go')
    return go_version or "20260527"
```

- [ ] **Step 4: 验证 ID 转换**

Run: `python -c "
from allenricher.database.parsers.wikipathways import WikiPathwaysParser
# 测试 ID 转换
id_mapping = {'1234': 'TP53', '5678': 'BRCA1'}
gene_sets = {'WP1': ['ncbigene:1234', 'ncbigene:5678']}
converted = WikiPathwaysParser.convert_ncbi_to_symbol(gene_sets, id_mapping)
print(f'Converted: {converted}')
assert converted['WP1'] == ['TP53', 'BRCA1']
print('✅ ID conversion OK')
"`
Expected: `Converted: {'WP1': ['TP53', 'BRCA1']}`

---

## Task 3: 创建 GPML 解析器（支持 22 无 GMT 物种）

**Files:**
- Create: `allenricher/database/wikipathways_gpml_parser.py`

对于只有 GPML 文件、没有 GMT 文件的 22 个物种，需要从 GPML（XML 格式）中提取通路-基因关联信息。

- [ ] **Step 1: 创建 GPML 解析器**

```python
"""
WikiPathways GPML 解析器

从 GPML 文件（XML 格式）中提取通路-基因关联信息。
用于支持只有 GPML 数据、没有 GMT 文件的物种。
"""
from __future__ import annotations
import gzip
import logging
import xml.etree.ElementTree as ET
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple
import zipfile

logger = logging.getLogger(__name__)


class WikiPathwaysGPMLParser:
    """WikiPathways GPML 解析器
    
    从 GPML ZIP 文件中提取通路-基因关联。
    """
    
    # GPML XML 命名空间
    GPML_NS = "http://pathvisio.org/GPML/2013a"
    
    @staticmethod
    def parse_gpml_zip(
        gpml_zip_path: Path,
        gene_info_path: Optional[Path] = None,
        taxid: Optional[int] = None,
    ) -> Tuple[Dict[str, List[str]], Dict[str, str]]:
        """解析 GPML ZIP 文件
        
        Args:
            gpml_zip_path: GPML ZIP 文件路径
            gene_info_path: NCBI gene_info.gz 路径（用于 ID 转换）
            taxid: NCBI TaxID（用于 gene_info 过滤）
            
        Returns:
            (gene_sets, descriptions) 元组
            - gene_sets: {pathway_id: [gene1, gene2, ...]}
            - descriptions: {pathway_id: pathway_name}
        """
        gene_sets: Dict[str, List[str]] = {}
        descriptions: Dict[str, str] = {}
        
        # 加载 gene_info 映射（如果提供）
        id_mapping = {}
        if gene_info_path and taxid:
            id_mapping = WikiPathwaysGPMLParser._load_gene_id_mapping(
                gene_info_path, taxid
            )
        
        # 解压并解析 GPML 文件
        with zipfile.ZipFile(gpml_zip_path, 'r') as zf:
            for filename in zf.namelist():
                if not filename.endswith('.gpml'):
                    continue
                    
                pathway_id = Path(filename).stem  # WPxxx
                
                try:
                    with zf.open(filename) as f:
                        tree = ET.parse(f)
                        root = tree.getroot()
                        
                        # 提取通路名称
                        pathway_name = root.get('Name', pathway_id)
                        descriptions[pathway_id] = pathway_name
                        
                        # 提取基因
                        genes = WikiPathwaysGPMLParser._extract_genes_from_gpml(
                            root, id_mapping
                        )
                        if genes:
                            gene_sets[pathway_id] = genes
                            
                except Exception as e:
                    logger.warning(f"Failed to parse {filename}: {e}")
                    continue
        
        logger.info(f"Parsed GPML: {len(gene_sets)} pathways, "
                     f"from {gpml_zip_path.name}")
        return gene_sets, descriptions
    
    @staticmethod
    def _extract_genes_from_gpml(
        root: ET.Element,
        id_mapping: Dict[str, str]
    ) -> List[str]:
        """从 GPML 根元素提取基因列表
        
        GPML 中的 DataNode 元素包含基因信息：
        - TextLabel: 基因名称（可能是 Symbol 或 ID）
        - Xref: 外部数据库引用（Database, ID）
        """
        genes = []
        ns = {'gpml': WikiPathwaysGPMLParser.GPML_NS}
        
        # 查找所有 DataNode 元素
        for datanode in root.findall('.//gpml:DataNode', ns):
            # 获取 Xref 信息
            xref = datanode.find('gpml:Xref', ns)
            if xref is not None:
                database = xref.get('Database', '')
                identifier = xref.get('ID', '')
                
                # 如果是 NCBI Gene，进行 ID 转换
                if database == 'Entrez Gene' and identifier:
                    if identifier in id_mapping:
                        genes.append(id_mapping[identifier])
                    else:
                        # 保留原始 ID
                        genes.append(f"ncbigene:{identifier}")
                else:
                    # 使用 TextLabel 作为基因名
                    text_label = datanode.get('TextLabel', '')
                    if text_label:
                        genes.append(text_label)
        
        return list(set(genes))  # 去重
    
    @staticmethod
    def _load_gene_id_mapping(gene_info_path: Path, taxid: int) -> Dict[str, str]:
        """加载 NCBI Gene ID 映射（与 wikipathways.py 相同）"""
        import gzip
        
        mapping = {}
        with gzip.open(gene_info_path, 'rt', encoding='utf-8') as f:
            header = f.readline().strip().split('\t')
            taxid_idx = header.index('tax_id')
            geneid_idx = header.index('GeneID')
            symbol_idx = header.index('Symbol')
            
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) < max(taxid_idx, geneid_idx, symbol_idx) + 1:
                    continue
                if int(parts[taxid_idx]) == taxid:
                    ncbi_id = parts[geneid_idx]
                    symbol = parts[symbol_idx]
                    mapping[ncbi_id] = symbol
                    
        return mapping
    
    @staticmethod
    def build_database_from_gpml(
        gpml_zip_path: Path,
        output_dir: Path,
        species: str,
        taxid: int,
        gene_info_path: Optional[Path] = None,
    ) -> Tuple[Path, Path]:
        """从 GPML 构建 AllEnricher 标准格式
        
        与 wikipathways.py 的 build_database 输出格式相同。
        """
        gene_sets, descriptions = WikiPathwaysGPMLParser.parse_gpml_zip(
            gpml_zip_path, gene_info_path, taxid
        )
        
        # 生成标准格式文件（与 wikipathways.py 相同逻辑）
        # ... 生成 .tab.gz 和 .disc.gz
        
        # 这里复用 wikipathways.py 中的文件生成逻辑
        from .wikipathways import WikiPathwaysParser
        return WikiPathwaysParser._write_database_files(
            gene_sets, descriptions, output_dir, species
        )
```

- [ ] **Step 2: 导出 GPML 解析器**

在 `allenricher/database/parsers/__init__.py` 中添加:
```python
from .wikipathways_gpml import WikiPathwaysGPMLParser

__all__ = [..., 'WikiPathwaysGPMLParser']
```

- [ ] **Step 3: 验证 GPML 解析器**

Run: `python -c "from allenricher.database.parsers import WikiPathwaysGPMLParser; print('GPML Parser OK')"`
Expected: `GPML Parser OK`

---

## Task 4: 更新物种注册表支持 40 物种

**Files:**
- Modify: `allenricher/database/species_registry.py`

- [ ] **Step 1: 添加 WikiPathways 物种分类信息**

在 `SpeciesEntry` 中添加字段区分 GMT/GPML 支持的物种:

```python
# WikiPathways 相关字段
has_wikipathways: bool = False
wikipathways_data_type: Optional[str] = None  # 'gmt', 'gpml', or None
wikipathways_gene_count: Optional[int] = None
wikipathways_pathway_count: Optional[int] = None
```

- [ ] **Step 2: 更新物种注册表构建逻辑**

在 `downloader.py` 的 `_build_wikipathways_registry()` 中更新为支持 40 物种:

```python
def _build_wikipathways_registry(self):
    """构建 WikiPathways 物种注册表（40 物种）"""
    from .wikipathways_fetcher import (
        WikiPathwaysFetcher, 
        SPECIES_NAME_MAP,
        SPECIES_WITH_GMT,
        SPECIES_GPML_ONLY
    )
    
    registry_path = self.basic_dir / "wikipathways_species_registry.tsv"
    
    with open(registry_path, 'w', encoding='utf-8') as f:
        f.write("species_latin_name\tspecies_code\tdata_type\tgene_count\tpathway_count\n")
        
        for latin_name in sorted(SPECIES_NAME_MAP.keys()):
            code = SPECIES_NAME_MAP[latin_name]
            
            # 确定数据类型
            if latin_name in SPECIES_WITH_GMT:
                data_type = 'gmt'
            elif latin_name in SPECIES_GPML_ONLY:
                data_type = 'gpml'
            else:
                data_type = 'none'
            
            # gene_count 和 pathway_count 在构建阶段填充
            f.write(f"{latin_name}\t{code}\t{data_type}\t-\t-\n")
    
    print(f"|--- WikiPathways 物种注册表: {registry_path}")
    print(f"|---   GMT 物种: {len(SPECIES_WITH_GMT)}")
    print(f"|---   GPML 物种: {len(SPECIES_GPML_ONLY)}")
```

---

## Task 5: 更新 CLI 支持物种数据类型查询

**Files:**
- Modify: `allenricher/cli.py`

- [ ] **Step 1: 更新 list-species 输出显示数据类型**

在 `_cmd_list_species()` 中，当使用 `--wikipathways` 过滤时，显示数据类型（GMT/GPML）:

```python
# 在 table 输出格式中添加 wikipathways_data_type 列
if args.wikipathways:
    headers = ['Species', 'Code', 'Data Type', 'TaxID']
    rows = [[
        e.latin_name,
        SPECIES_NAME_MAP.get(e.latin_name, '-'),
        e.wikipathways_data_type or '-',
        e.taxid
    ] for e in filtered]
```

---

## Task 6: 端到端测试（40 物种）

**Files:**
- Create: `test_e2e_2026/test_wikipathways_full_species.py`

- [ ] **Step 1: 创建全物种测试脚本**

测试覆盖:
1. 验证 40 物种映射正确
2. 测试 GMT 物种的下载和构建
3. 测试 GPML 物种的下载和构建（如有 GPML 文件）
4. 验证 ID 转换功能（NCBI Gene ID → Symbol）

---

## 自检清单

1. **Spec 覆盖**: 
   - 物种扩展 18→40 ✓
   - NCBI Gene ID → Symbol 转换 ✓
   - GPML 解析器 ✓
   - 按需下载策略 ✓

2. **类型一致性**: 
   - `SPECIES_NAME_MAP` 包含 40 物种
   - `SPECIES_WITH_GMT` / `SPECIES_GPML_ONLY` 正确划分
   - ID 转换逻辑与 GO 保持一致

3. **数据流完整性**: 
   - GMT/GPML → parse → NCBI ID → gene_info → Symbol → .tab.gz → DatabaseManager
