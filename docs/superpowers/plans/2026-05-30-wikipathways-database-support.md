# WikiPathways 数据库支持 — 全流程实现计划

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** 为 AllEnricher v2 添加 WikiPathways 数据库的完整支持，包括数据下载、物种注册表扩展、物种数据库构建、富集分析集成。

**Architecture:** WikiPathways 数据从 `data.wikipathways.org` 获取每月发布的 GMT 文件（按物种分离），通过新建 `WikiPathwaysParser` 解析为 AllEnricher 标准格式（`.tab.gz` + `.disc.gz`），与现有 GO/KEGG/Reactome 管道完全对齐。运行时 `EnrichmentAnalyzer` 无需修改，因为它只消费 `{term_id: {name, genes}}` 字典。

**Tech Stack:** Python 3.x, requests (下载), gzip/pandas (解析), dataclasses (数据结构)

---

## 前置知识：WikiPathways 数据源

### 数据下载地址
- **GMT 文件（按物种）**: `https://data.wikipathways.org/{YYYYMMDD}/gmt/wikipathways-{YYYYMMDD}-gmt-{Species_Latin_Name}.gmt`
  - 例: `https://data.wikipathways.org/20260510/gmt/wikipathways-20260510-gmt-Homo_sapiens.gmt`
  - 每月更新，保留最近 12 个月
  - 18 个物种有 GMT 文件（2026-05-10 版本）
- **RDF 数据（全物种）**: `https://data.wikipathways.org/{YYYYMMDD}/rdf/wikipathways-{YYYYMMDD}-rdf-wp.zip`
  - 包含所有物种的通路-基因关联（TTL 格式）
  - 可作为 GMT 文件的补充/替代数据源
- **GPML 数据（原始格式）**: `https://data.wikipathways.org/{YYYYMMDD}/gpml/`
  - 原始通路图数据，解析复杂，不建议直接使用

### 支持的物种（GMT 文件）
Anopheles_gambiae, Arabidopsis_thaliana, Bos_taurus, Caenorhabditis_elegans, Canis_familiaris, Danio_rerio, Drosophila_melanogaster, Equus_caballus, Gallus_gallus, Homo_sapiens, Mus_musculus, Pan_troglodytes, Populus_trichocarpa, Rattus_norvegicus, Saccharomyces_cerevisiae, Solanum_lycopersicum, Sus_scrofa, Zea_mays

### GMT 文件格式
```
WPID (e.g. WP254)    Pathway Name    Gene1/Gene2/...    (基因用 / 分隔，非 Tab)
```
**注意**: WikiPathways GMT 格式与标准 GMT 不同——基因用 `/` 分隔而非 Tab。

### 物种名到拉丁名的映射
GMT 文件名使用 `Species_Latin_Name`（下划线分隔），需要建立与 AllEnricher 物种代码（如 hsa）的映射。

---

## 文件变更清单

| 操作 | 文件路径 | 职责 |
|------|---------|------|
| 新建 | `allenricher/database/parsers/wikipathways.py` | WikiPathways 解析器：GMT→标准格式 |
| 新建 | `allenricher/database/wikipathways_fetcher.py` | WikiPathways 数据获取器（下载+缓存） |
| 修改 | `allenricher/database/parsers/__init__.py` | 导出 WikiPathwaysParser |
| 修改 | `allenricher/database/mirrors.py` | 添加 WikiPathways 镜像源 |
| 修改 | `allenricher/database/downloader.py` | 添加 download_wikipathways_basic() |
| 修改 | `allenricher/database/version.py` | 添加 WikiPathways 版本管理 |
| 修改 | `allenricher/database/species_registry.py` | SpeciesEntry 添加 has_wikipathways 字段 |
| 修改 | `allenricher/database/builder.py` | build_species_db 添加 WIKIPATHWAYS 分支 |
| 修改 | `allenricher/database/manager.py` | name_to_prefix 添加 WIKIPATHWAYS 映射 |
| 修改 | `allenricher/database/gmt_generator.py` | 添加 generate_wikipathways_gmt() |
| 修改 | `allenricher/cli.py` | list-species/query-species/build 支持 WikiPathways |

---

## Task 1: 创建 WikiPathways 数据获取器

**Files:**
- Create: `allenricher/database/wikipathways_fetcher.py`

- [ ] **Step 1: 创建 WikiPathwaysFetcher 类**

```python
"""
WikiPathways 数据获取器

从 data.wikipathways.org 下载按物种分离的 GMT 文件，
并缓存到 database/basic/wikipathways/ 目录。
"""
from __future__ import annotations
import gzip
import logging
from pathlib import Path
from typing import Dict, List, Optional, Tuple
from datetime import datetime

logger = logging.getLogger(__name__)

# WikiPathways 物种拉丁名 → AllEnricher 物种代码映射
# 键: GMT 文件名中的物种名（下划线分隔）
# 值: AllEnricher 物种代码（3-4字母小写）
SPECIES_NAME_MAP: Dict[str, str] = {
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
}

# 反向映射: AllEnricher 物种代码 → WikiPathways 拉丁名
SPECIES_CODE_TO_NAME: Dict[str, str] = {v: k for k, v in SPECIES_NAME_MAP.items()}

# 基础数据下载 URL
BASE_URL = "https://data.wikipathways.org"


class WikiPathwaysFetcher:
    """WikiPathways 数据获取器

    从 data.wikipathways.org 下载 GMT 文件并缓存。
    """

    def __init__(self, basic_dir: str = "./database/basic"):
        self.basic_dir = Path(basic_dir)
        self.wp_dir = self.basic_dir / "wikipathways"

    def get_available_species(self, version: Optional[str] = None) -> List[str]:
        """获取指定版本中有 GMT 文件的物种列表

        Args:
            version: 版本号（如 '20260510'），默认自动获取最新

        Returns:
            物种拉丁名列表（下划线分隔格式）
        """
        if version is None:
            version = self._detect_latest_version()
        if version is None:
            logger.warning("无法检测 WikiPathways 最新版本")
            return []

        url = f"{BASE_URL}/{version}/gmt/"
        # 通过 HTTP 请求获取目录列表，解析物种文件名
        import requests
        try:
            resp = requests.get(url, timeout=30)
            resp.raise_for_status()
            # 解析 HTML 表格中的文件名
            species_list = []
            for line in resp.text.split('\n'):
                if 'gmt' in line and '.gmt' in line:
                    # 提取物种名: wikipathways-{version}-gmt-{Species_Name}.gmt
                    import re
                    match = re.search(r'wikipathways-\d+-gmt-([A-Za-z_]+)\.gmt', line)
                    if match:
                        species_list.append(match.group(1))
            return species_list
        except Exception as e:
            logger.error(f"获取 WikiPathways 物种列表失败: {e}")
            return []

    def _detect_latest_version(self) -> Optional[str]:
        """从 data.wikipathways.org 检测最新版本号

        Returns:
            版本字符串（如 '20260510'），或 None
        """
        import requests
        try:
            resp = requests.get(f"{BASE_URL}/", timeout=30)
            resp.raise_for_status()
            # 解析 HTML，找最新的日期目录名
            import re
            versions = re.findall(r'>(\d{8})<', resp.text)
            if versions:
                return sorted(versions)[-1]
        except Exception as e:
            logger.error(f"检测 WikiPathways 最新版本失败: {e}")
        return None

    def download_gmt(self, species_latin_name: str, version: Optional[str] = None,
                     overwrite: bool = False) -> Optional[Path]:
        """下载单个物种的 GMT 文件

        Args:
            species_latin_name: 物种拉丁名（下划线分隔，如 Homo_sapiens）
            version: 版本号，默认自动获取最新
            overwrite: 是否覆盖已有文件

        Returns:
            下载文件的路径，或 None（失败时）
        """
        if version is None:
            version = self._detect_latest_version()
        if version is None:
            return None

        self.wp_dir.mkdir(parents=True, exist_ok=True)
        version_dir = self.wp_dir / f"WP{version}"
        version_dir.mkdir(parents=True, exist_ok=True)

        filename = f"wikipathways-{version}-gmt-{species_latin_name}.gmt"
        dest_path = version_dir / filename

        if dest_path.exists() and not overwrite:
            logger.info(f"WikiPathways GMT 已存在: {dest_path}")
            return dest_path

        url = f"{BASE_URL}/{version}/gmt/{filename}"
        logger.info(f"下载 WikiPathways GMT: {url}")

        import requests
        try:
            resp = requests.get(url, timeout=120, stream=True)
            resp.raise_for_status()
            with open(dest_path, 'wb') as f:
                for chunk in resp.iter_content(chunk_size=8192):
                    f.write(chunk)
            logger.info(f"WikiPathways GMT 下载完成: {dest_path} ({dest_path.stat().st_size} bytes)")
            return dest_path
        except Exception as e:
            logger.error(f"下载 WikiPathways GMT 失败 ({species_latin_name}): {e}")
            return None

    def download_all_gmt(self, version: Optional[str] = None,
                         overwrite: bool = False) -> Dict[str, Path]:
        """下载所有物种的 GMT 文件

        Args:
            version: 版本号，默认自动获取最新
            overwrite: 是否覆盖已有文件

        Returns:
            {物种拉丁名: 文件路径} 字典
        """
        species_list = self.get_available_species(version)
        results = {}
        for species_name in species_list:
            path = self.download_gmt(species_name, version, overwrite)
            if path:
                results[species_name] = path
        return results

    def get_species_code(self, latin_name: str) -> Optional[str]:
        """将 WikiPathways 拉丁名转换为 AllEnricher 物种代码"""
        return SPECIES_NAME_MAP.get(latin_name)

    def get_latin_name(self, species_code: str) -> Optional[str]:
        """将 AllEnricher 物种代码转换为 WikiPathways 拉丁名"""
        return SPECIES_CODE_TO_NAME.get(species_code)
```

- [ ] **Step 2: 验证文件可导入且无语法错误**

Run: `python -c "from allenricher.database.wikipathways_fetcher import WikiPathwaysFetcher, SPECIES_NAME_MAP; print(f'Species count: {len(SPECIES_NAME_MAP)}'); print(SPECIES_NAME_MAP)"`
Expected: 打印物种映射字典，18 个物种

---

## Task 2: 创建 WikiPathways 解析器

**Files:**
- Create: `allenricher/database/parsers/wikipathways.py`
- Modify: `allenricher/database/parsers/__init__.py`

- [ ] **Step 1: 创建 WikiPathwaysParser 类**

```python
"""
WikiPathways 解析器

将 WikiPathways GMT 文件解析为 AllEnricher 标准格式：
- {species}.WikiPathways2gene.tab.gz  (基因-通路 0/1 矩阵)
- {species}.WikiPathways2disc.gz     (通路描述)
"""
from __future__ import annotations
import gzip
import logging
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple

logger = logging.getLogger(__name__)


class WikiPathwaysParser:
    """WikiPathways 数据解析器

    将 WikiPathways GMT 文件转换为 AllEnricher 标准格式。
    """

    @staticmethod
    def parse_gmt(gmt_path: Path) -> Tuple[Dict[str, List[str]], Dict[str, str]]:
        """解析 WikiPathways GMT 文件

        WikiPathways GMT 格式:
            WPID<TAB>Pathway Name<TAB>Gene1/Gene2/Gene3/...

        注意: 基因用 '/' 分隔（非标准 GMT 的 Tab 分隔）

        Args:
            gmt_path: GMT 文件路径

        Returns:
            (gene_sets, descriptions) 元组
            - gene_sets: {pathway_id: [gene1, gene2, ...]}
            - descriptions: {pathway_id: pathway_name}
        """
        gene_sets: Dict[str, List[str]] = {}
        descriptions: Dict[str, str] = {}

        with open(gmt_path, 'r', encoding='utf-8') as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith('#'):
                    continue
                parts = line.split('\t')
                if len(parts) < 3:
                    continue

                pathway_id = parts[0].strip()
                pathway_name = parts[1].strip()
                # WikiPathways GMT: 基因用 '/' 分隔
                genes_str = parts[2].strip()
                genes = [g.strip() for g in genes_str.split('/') if g.strip()]

                if pathway_id and genes:
                    gene_sets[pathway_id] = genes
                    descriptions[pathway_id] = pathway_name

        logger.info(f"解析 WikiPathways GMT: {len(gene_sets)} 条通路, "
                     f"来自 {gmt_path.name}")
        return gene_sets, descriptions

    @staticmethod
    def build_database(
        gmt_path: Path,
        output_dir: Path,
        species: str,
        gene_info_path: Optional[Path] = None,
        valid_genes: Optional[Set[str]] = None,
    ) -> Tuple[Path, Path]:
        """构建 AllEnricher 标准格式文件

        Args:
            gmt_path: WikiPathways GMT 文件路径
            output_dir: 输出目录
            species: 物种代码（如 hsa）
            gene_info_path: NCBI gene_info.gz 路径（用于基因符号验证，可选）
            valid_genes: 有效基因集合（如果提供，只保留其中的基因）

        Returns:
            (tab_path, disc_path) 元组
            - tab_path: {species}.WikiPathways2gene.tab.gz
            - disc_path: {species}.WikiPathways2disc.gz
        """
        gene_sets, descriptions = WikiPathwaysParser.parse_gmt(gmt_path)

        # 如果提供了有效基因集合，过滤基因
        if valid_genes:
            filtered_sets = {}
            for wp_id, genes in gene_sets.items():
                filtered_genes = [g for g in genes if g in valid_genes]
                if filtered_genes:
                    filtered_sets[wp_id] = filtered_genes
            gene_sets = filtered_sets
            logger.info(f"基因过滤后: {len(gene_sets)} 条通路 (有效基因集大小: {len(valid_genes)})")

        # 收集所有基因和通路 ID
        all_genes = set()
        for genes in gene_sets.values():
            all_genes.update(genes)
        all_pathways = sorted(gene_sets.keys())

        if not all_genes or not all_pathways:
            logger.warning("WikiPathways 解析结果为空，跳过构建")
            return Path(""), Path("")

        # 生成 .tab.gz 文件（0/1 矩阵）
        tab_filename = f"{species}.WikiPathways2gene.tab.gz"
        tab_path = output_dir / tab_filename
        with gzip.open(tab_path, 'wt', encoding='utf-8') as f:
            # 表头: Gene\tWPID1\tWPID2\t...
            header = ['Gene'] + all_pathways
            f.write('\t'.join(header) + '\n')
            # 数据行
            for gene in sorted(all_genes):
                row = [gene]
                for wp_id in all_pathways:
                    row.append('1' if gene in gene_sets[wp_id] else '0')
                f.write('\t'.join(row) + '\n')
        logger.info(f"生成 {tab_path} ({len(all_genes)} 基因 x {len(all_pathways)} 通路)")

        # 生成 .disc.gz 文件（通路描述）
        disc_filename = f"{species}.WikiPathways2disc.gz"
        disc_path = output_dir / disc_filename
        with gzip.open(disc_path, 'wt', encoding='utf-8') as f:
            for wp_id in all_pathways:
                desc = descriptions.get(wp_id, wp_id)
                f.write(f"{wp_id}\t{desc}\n")
        logger.info(f"生成 {disc_path} ({len(all_pathways)} 条描述)")

        return tab_path, disc_path
```

- [ ] **Step 2: 在 parsers/__init__.py 中导出**

修改 `allenricher/database/parsers/__init__.py`，添加:
```python
from .wikipathways import WikiPathwaysParser
```
并更新 `__all__`:
```python
__all__ = ['GOParser', 'KEGGParser', 'ReactomeParser', 'DOParser', 'DisGeNETParser', 'WikiPathwaysParser']
```

- [ ] **Step 3: 验证解析器可导入**

Run: `python -c "from allenricher.database.parsers import WikiPathwaysParser; print('OK')"`
Expected: OK

---

## Task 3: 扩展物种注册表

**Files:**
- Modify: `allenricher/database/species_registry.py`

- [ ] **Step 1: 在 _FIELD_NAMES 中添加 WikiPathways 列**

在 `_FIELD_NAMES` 列表中，`"synonyms"` 之前添加:
```python
"has_wikipathways", "wikipathways_gene_count", "wikipathways_pathway_count",
```

- [ ] **Step 2: 在 SpeciesEntry dataclass 中添加字段**

在 `do_term_count` 字段之后、`synonyms` 之前添加:
```python
# WikiPathways 相关字段
has_wikipathways: bool = False
wikipathways_gene_count: Optional[int] = None
wikipathways_pathway_count: Optional[int] = None
```

- [ ] **Step 3: 验证 SpeciesEntry 兼容性**

Run: `python -c "from allenricher.database.species_registry import SpeciesEntry, _FIELD_NAMES; e = SpeciesEntry(taxid=9606, latin_name='Homo sapiens'); print(e); print(_FIELD_NAMES)"`
Expected: 正常打印，新字段有默认值

---

## Task 4: 扩展数据库下载器

**Files:**
- Modify: `allenricher/database/mirrors.py`
- Modify: `allenricher/database/downloader.py`
- Modify: `allenricher/database/version.py`

- [ ] **Step 1: 在 mirrors.py 中添加 WikiPathways 镜像源**

在 `MIRROR_REGISTRY` 中添加:
```python
# WikiPathways 数据镜像
MirrorSource(
    name="wikipathways-official",
    base_url="https://data.wikipathways.org/",
    priority=1,
    region="US",
),
```

并在 `get_mirrors()` 函数的 `db_type` 分支中添加 `'wikipathways'` case:
```python
elif db_type == 'wikipathways':
    return [m for m in MIRROR_REGISTRY if m.name == "wikipathways-official" and m.enabled]
```

- [ ] **Step 2: 在 downloader.py 中添加 download_wikipathways_basic() 方法**

在 `DataDownloader` 类中添加新方法:
```python
def download_wikipathways_basic(self, version=None) -> str:
    """下载 WikiPathways 基础数据（所有物种的 GMT 文件）

    Args:
        version: 版本号（如 '20260510'），默认自动检测最新

    Returns:
        版本目录路径
    """
    from .wikipathways_fetcher import WikiPathwaysFetcher

    fetcher = WikiPathwaysFetcher(str(self.basic_dir))
    if version is None:
        version = fetcher._detect_latest_version()

    results = fetcher.download_all_gmt(version=version, overwrite=self.overwrite)

    version_str = version or "unknown"
    version_dir = str(self.basic_dir / "wikipathways" / f"WP{version_str}")

    # 记录版本
    from .version import DatabaseVersionManager
    vm = DatabaseVersionManager(str(self.basic_dir))
    vm.record_download("wikipathways", version_str, version_dir)

    print(f"|--- WikiPathways: 下载了 {len(results)} 个物种的 GMT 文件")
    print(f"|--- 版本: WP{version_str}, 目录: {version_dir}")

    return version_dir
```

并在 `download_all()` 方法的数据库分支中添加:
```python
elif db_type in ('wikipathways',):
    result['wikipathways'] = self.download_wikipathways_basic()
```

- [ ] **Step 3: 在 version.py 的 SOURCE_URLS 中添加 WikiPathways**

在 `SOURCE_URLS` 字典中添加:
```python
"wikipathways": "https://data.wikipathways.org/",
```

- [ ] **Step 4: 在 downloader.py 的注册表构建流水线中添加 WikiPathways**

在 `download_all()` 方法末尾的注册表构建部分，添加 `_build_wikipathways_registry()` 调用:
```python
self._build_wikipathways_registry()
```

并在 `DataDownloader` 类中添加:
```python
def _build_wikipathways_registry(self):
    """构建 WikiPathways 物种注册表"""
    from .wikipathways_fetcher import WikiPathwaysFetcher, SPECIES_NAME_MAP

    fetcher = WikiPathwaysFetcher(str(self.basic_dir))
    registry_path = self.basic_dir / "wikipathways_species_registry.tsv"

    # 获取最新版本中的物种列表
    species_list = fetcher.get_available_species()

    with open(registry_path, 'w', encoding='utf-8') as f:
        f.write("species_latin_name\tspecies_code\tgene_count\tpathway_count\n")
        for latin_name in sorted(species_list):
            code = fetcher.get_species_code(latin_name)
            code_str = code or "-"
            # gene_count 和 pathway_count 在构建阶段填充，下载阶段留空
            f.write(f"{latin_name}\t{code_str}\t-\t-\n")

    print(f"|--- WikiPathways 物种注册表: {registry_path} ({len(species_list)} 物种)")
```

- [ ] **Step 5: 验证下载命令**

Run: `python -m allenricher download -d wikipathways --force 2>&1 | head -20`
Expected: 显示下载进度，最终输出 "WikiPathways: 下载了 N 个物种的 GMT 文件"

---

## Task 5: 扩展物种数据库构建器

**Files:**
- Modify: `allenricher/database/builder.py`
- Modify: `allenricher/database/gmt_generator.py`

- [ ] **Step 1: 在 builder.py 的 build_species_db() 中添加 WIKIPATHWAYS 分支**

在 `build_species_db()` 方法的 `for db_name in databases:` 循环中，`elif db_upper == 'DISGENET':` 之后添加:
```python
elif db_upper == 'WIKIPATHWAYS':
    self.build_wikipathways(species, taxid, outdir)
```

并在 `DatabaseBuilder` 类中添加 `build_wikipathways()` 方法:
```python
def build_wikipathways(self, species: str, taxid: int, outdir: Path) -> None:
    """构建 WikiPathways 数据库

    从 database/basic/wikipathways/ 中找到对应物种的 GMT 文件，
    解析为 AllEnricher 标准格式。

    Args:
        species: 物种代码（如 hsa）
        taxid: NCBI TaxID
        outdir: 输出目录
    """
    from .wikipathways_fetcher import WikiPathwaysFetcher
    from .parsers.wikipathways import WikiPathwaysParser

    fetcher = WikiPathwaysFetcher(str(self.basic_dir))
    latin_name = fetcher.get_latin_name(species)

    if not latin_name:
        print(f"|--- [跳过] WikiPathways 不支持物种 {species}（无拉丁名映射）")
        return

    # 查找已下载的 GMT 文件
    wp_dir = self.basic_dir / "wikipathways"
    gmt_file = None
    if wp_dir.exists():
        # 搜索所有版本目录
        for version_dir in sorted(wp_dir.iterdir(), reverse=True):
            if version_dir.is_dir() and version_dir.name.startswith("WP"):
                candidate = version_dir / f"wikipathways-{version_dir.name[2:]}-gmt-{latin_name}.gmt"
                if candidate.exists():
                    gmt_file = candidate
                    break

    if not gmt_file:
        print(f"|--- [跳过] WikiPathways 未找到 {latin_name} 的 GMT 文件，请先执行 download -d wikipathways")
        return

    print(f"|--- 构建 WikiPathways: {species} ({latin_name})")
    print(f"|--- GMT 文件: {gmt_file}")

    # 获取有效基因集合（从 gene_info.gz）
    valid_genes = self._load_gene_symbols(taxid)

    # 构建标准格式
    tab_path, disc_path = WikiPathwaysParser.build_database(
        gmt_path=gmt_file,
        output_dir=outdir,
        species=species,
        valid_genes=valid_genes,
    )

    if tab_path.exists() and disc_path.exists():
        # 生成 GMT 文件（供 GSEA 使用）
        from .gmt_generator import GMTGenerator
        gen = GMTGenerator(outdir)
        gen.generate_wikipathways_gmt(species)
        print(f"|--- WikiPathways 构建完成: {tab_path.name}, {disc_path.name}")
    else:
        print(f"|--- [警告] WikiPathways 构建失败: 输出文件为空")
```

注意: `build_wikipathways` 需要访问 `self.basic_dir`，确认 `DatabaseBuilder.__init__` 中已有此属性。如果没有，需要通过参数传入。

- [ ] **Step 2: 在 gmt_generator.py 中添加 generate_wikipathways_gmt() 方法**

在 `GMTGenerator` 类中添加:
```python
def generate_wikipathways_gmt(self, species: str) -> Optional[str]:
    """从 WikiPathways2gene.tab.gz 生成标准 GMT 文件

    Args:
        species: 物种代码

    Returns:
        生成的 GMT 文件路径，或 None
    """
    tab_file = self.output_dir / f"{species}.WikiPathways2gene.tab.gz"
    disc_file = self.output_dir / f"{species}.WikiPathways2disc.gz"

    if not tab_file.exists():
        logger.warning(f"WikiPathways tab 文件不存在: {tab_file}")
        return None

    return self._generate_gmt_from_tab(
        tab_file=tab_file,
        disc_file=disc_file,
        db_name="WikiPathways",
        species=species,
    )
```

- [ ] **Step 3: 验证构建命令**

Run: `python -m allenricher build -s hsa -t 9606 -d WikiPathways --database-dir ./database 2>&1 | head -20`
Expected: 显示 "构建 WikiPathways: hsa (Homo_sapiens)"，生成 `.WikiPathways2gene.tab.gz` 和 `.WikiPathways2disc.gz`

---

## Task 6: 扩展运行时数据库加载器

**Files:**
- Modify: `allenricher/database/manager.py`

- [ ] **Step 1: 在 load_database() 的 name_to_prefix 中添加映射**

在 `load_database()` 方法中:
```python
name_to_prefix = {
    'GO': 'GO',
    'KEGG': 'kegg',
    'REACTOME': 'Reactome',
    'DO': 'DO',
    'DISGENET': 'CUI',
    'WIKIPATHWAYS': 'WikiPathways',  # 新增
}
```

- [ ] **Step 2: 在 _load_term_names() 的 name_to_prefix 中添加映射**

在 `_load_term_names()` 方法中同样添加:
```python
name_to_prefix = {
    'GO': 'GO',
    'KEGG': 'kegg',
    'REACTOME': 'Reactome',
    'DO': 'DO',
    'DISGENET': 'CUI',
    'WIKIPATHWAYS': 'WikiPathways',  # 新增
}
```

- [ ] **Step 3: 在 _format_term_name() 中添加 WikiPathways 格式化（如需要）**

检查 `_format_term_name()` 是否需要特殊处理 WikiPathways 的通路名称。WikiPathways 通路名通常不需要层级格式化（不像 KEGG 的三层分类），可以直接使用原始名称。

如果存在 `_format_term_name()` 方法，添加:
```python
elif db_name.upper() == 'WIKIPATHWAYS':
    # WikiPathways 通路名直接使用，无需格式化
    return name
```

- [ ] **Step 4: 验证加载**

Run: `python -c "
from allenricher.database.manager import DatabaseManager
dm = DatabaseManager('./database', 'hsa')
# 检查文件查找逻辑
print('WikiPathways prefix mapping OK')
"`
Expected: 无报错

---

## Task 7: 扩展 CLI 命令

**Files:**
- Modify: `allenricher/cli.py`

- [ ] **Step 1: list-species 添加 --wikipathways 过滤**

在 `list_species_parser` 的参数定义中添加:
```python
list_species_parser.add_argument('--wikipathways', action='store_true', default=False, help='Filter by WikiPathways support')
```

在 `_cmd_list_species()` 函数中，将 `--wikipathways` 传入 `registry.filter_by_databases()`:
```python
registry.filter_by_databases(
    go=args.go, kegg=args.kegg, reactome=args.reactome, do=args.do,
    wikipathways=args.wikipathways,  # 新增
)
```

同时需要修改 `SpeciesRegistry.filter_by_databases()` 方法，添加 `wikipathways` 参数:
```python
def filter_by_databases(self, go=False, kegg=False, reactome=False, do=False, wikipathways=False) -> List[SpeciesEntry]:
    results = self._entries.values()
    if wikipathways:
        results = [e for e in results if e.has_wikipathways]
    # ... 其余过滤逻辑不变
```

- [ ] **Step 2: build 命令默认数据库列表添加 WikiPathways**

修改 `build_parser` 的 `--databases` 默认值:
```python
build_parser.add_argument('-d', '--databases', default='GO,KEGG,Reactome,WikiPathways', ...)
```

注意: 这是一个可选改动。也可以不修改默认值，让用户手动指定 `-d WikiPathways`。

- [ ] **Step 3: analyze 命令无需修改**

`analyze` 命令通过 `DatabaseManager.load_databases()` 加载数据库，只要 `manager.py` 中的映射正确，`-d GO,KEGG,WikiPathways` 即可自动工作。

- [ ] **Step 4: 验证 CLI**

Run: `python -m allenricher list-species --wikipathways --format table 2>&1 | head -20`
Expected: 列出支持 WikiPathways 的物种

---

## Task 8: 端到端测试

**Files:**
- Create: `test_e2e_2026/test_wikipathways.py`

- [ ] **Step 1: 编写端到端测试脚本**

```python
"""
WikiPathways 端到端测试

测试 WikiPathways 数据库的完整流程：
1. 下载 → 2. 构建 → 3. 富集分析 → 4. 输出验证
"""
import subprocess
import sys
from pathlib import Path

PROJECT_DIR = Path(__file__).parent.parent
OUTPUT_DIR = PROJECT_DIR / "test_wikipathways_output"
SPECIES = "hsa"
TAXID = 9606


def run_cmd(cmd: str, cwd=None):
    """执行命令并返回结果"""
    result = subprocess.run(
        cmd, shell=True, capture_output=True, text=True, cwd=cwd or str(PROJECT_DIR)
    )
    return result.returncode, result.stdout, result.stderr


def test_download():
    """测试 WikiPathways 数据下载"""
    print("\n=== 测试 1: 下载 WikiPathways 数据 ===")
    rc, stdout, stderr = run_cmd(
        f"python -m allenricher download -d wikipathways --force"
    )
    print(stdout[-500:] if len(stdout) > 500 else stdout)
    assert rc == 0, f"下载失败: rc={rc}"
    print("✅ 下载成功")


def test_build():
    """测试 WikiPathways 物种数据库构建"""
    print("\n=== 测试 2: 构建 WikiPathways 数据库 ===")
    rc, stdout, stderr = run_cmd(
        f"python -m allenricher build -s {SPECIES} -t {TAXID} -d WikiPathways"
    )
    print(stdout[-500:] if len(stdout) > 500 else stdout)
    assert rc == 0, f"构建失败: rc={rc}"

    # 检查输出文件
    organism_dir = PROJECT_DIR / "database" / "organism"
    # 找到最新版本目录
    version_dirs = sorted(organism_dir.glob("v*"), reverse=True)
    assert len(version_dirs) > 0, "未找到任何版本目录"

    species_dir = version_dirs[0] / SPECIES
    assert species_dir.exists(), f"物种目录不存在: {species_dir}"

    expected_files = [
        f"{SPECIES}.WikiPathways2gene.tab.gz",
        f"{SPECIES}.WikiPathways2disc.gz",
        f"{SPECIES}.WikiPathways.gmt.gz",
    ]
    for f in expected_files:
        assert (species_dir / f).exists(), f"文件不存在: {species_dir / f}"
    print(f"✅ 构建成功，文件位于: {species_dir}")


def test_analyze():
    """测试 WikiPathways 富集分析"""
    print("\n=== 测试 3: WikiPathways 富集分析 ===")

    # 使用已有的基因列表
    gene_list = PROJECT_DIR / "test_e2e_2026" / "00_input_data" / "hsa_gene_list_500.txt"
    if not gene_list.exists():
        print("⚠️ 基因列表文件不存在，跳过分析测试")
        return

    if OUTPUT_DIR.exists():
        import shutil
        shutil.rmtree(OUTPUT_DIR)

    rc, stdout, stderr = run_cmd(
        f"python -m allenricher analyze -i {gene_list} -m hypergeometric "
        f"-s {SPECIES} -d WikiPathways -o {OUTPUT_DIR} "
        f"--use-version {find_latest_version()} -p 1 -q 1"
    )
    print(stdout[-1000:] if len(stdout) > 1000 else stdout)

    # 检查输出文件
    tsv_file = OUTPUT_DIR / "WikiPathways_enrichment.tsv"
    assert tsv_file.exists(), f"TSV 文件不存在: {tsv_file}"

    # 检查 TSV 列名
    with open(tsv_file, 'r', encoding='utf-8') as f:
        header_line = ""
        for line in f:
            if line.startswith('#'):
                continue
            header_line = line.strip()
            break

    expected_cols = ['Term_ID', 'Term_Name', 'P_Value']
    for col in expected_cols:
        assert col in header_line, f"缺少列: {col}, 实际列: {header_line}"

    print(f"✅ 分析成功，输出: {tsv_file}")


def find_latest_version() -> str:
    """找到最新的数据库版本"""
    organism_dir = PROJECT_DIR / "database" / "organism"
    version_dirs = sorted(organism_dir.glob("v*"), reverse=True)
    if version_dirs:
        return version_dirs[0].name
    return "v20260515"


if __name__ == "__main__":
    try:
        test_download()
        test_build()
        test_analyze()
        print("\n" + "=" * 60)
        print("✅ 所有 WikiPathways 端到端测试通过！")
        print("=" * 60)
    except AssertionError as e:
        print(f"\n❌ 测试失败: {e}")
        sys.exit(1)
    except Exception as e:
        print(f"\n❌ 测试异常: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)
```

- [ ] **Step 2: 运行端到端测试**

Run: `cd AllEnricher-v2 && python test_e2e_2026/test_wikipathways.py 2>&1`
Expected: 所有测试通过

---

## Task 9: 扩展四种分析方法的图表和报告支持

WikiPathways 数据库接入后，需要确保四种分析方法（ORA、GSEA、ssGSEA、GSVA）的图表生成和 HTML 报告都能正确显示 WikiPathways 结果。

### 9.1 图表系统架构回顾

**图表类型映射** (`_METHOD_PLOT_TYPES` in `cli.py`):
```python
_METHOD_PLOT_TYPES = {
    'gsea': {'enrichment', 'nes_barplot', 'dotplot'},
    'ssgsea': {'heatmap', 'group_comparison', 'dotplot', 'correlation'},
    'gsva': {'heatmap', 'group_comparison', 'dotplot', 'correlation'},
}
_COMMON_PLOT_TYPES = {'network', 'upset', 'volcano'}
```

**图表生成流程**:
1. CLI 根据 `--method` 和 `--plot-types` 决定生成哪些图表
2. `_generate_plots()` 函数分派到各方法的专用绘图函数
3. 所有数据库（GO/KEGG/Reactome/WikiPathways）共享同一套绘图逻辑

### 9.2 ORA 分析图表

**Files:**
- Modify: `allenricher/visualization/plotter.py`

**验证点**: WikiPathways 结果使用 ORA 格式（含 Rich_Factor, Gene_Ratio 等列），plotter 的 `_prepare_barplot_data()` 和 `_prepare_bubble_data()` 已支持动态列名映射，无需修改。

- [ ] **Step 1: 验证 ORA 图表对 WikiPathways 的兼容性**

Run: `python -c "
from allenricher.visualization.plotter import Plotter
import pandas as pd
# 模拟 WikiPathways ORA 结果
df = pd.DataFrame({
    'Term_ID': ['WP1', 'WP2'],
    'Term_Name': ['Pathway A', 'Pathway B'],
    'Gene_Count': [10, 20],
    'Rich_Factor': [0.5, 0.3],
    'P_Value': [0.01, 0.02],
    'Adjusted_P_Value': [0.02, 0.04],
    'Database': ['WikiPathways', 'WikiPathways']
})
# 验证列名映射
from allenricher.visualization.plotter import Plotter
p = Plotter('./test')
col_map = p._get_column_mapping(df)
print('Column mapping:', col_map)
assert 'gene_count' in col_map or 'setSize' in str(col_map)
print('✅ ORA 图表兼容')
"`
Expected: 正常打印列名映射，无报错

### 9.3 GSEA 分析图表

**Files:**
- Modify: `allenricher/visualization/gsea_plots.py` (如有需要)

**验证点**: GSEA 结果通过 `is_gsea = 'NES' in df.columns` 检测，WikiPathways GSEA 结果包含 NES 列，自动触发 GSEA 专用图表。

- [ ] **Step 2: 验证 GSEA 图表对 WikiPathways 的兼容性**

Run: `python -c "
import pandas as pd
# 模拟 WikiPathways GSEA 结果
df = pd.DataFrame({
    'Term_ID': ['WP1', 'WP2'],
    'Term_Name': ['Pathway A', 'Pathway B'],
    'NES': [1.5, -1.2],
    'p_value': [0.01, 0.02],
    'FDR': [0.02, 0.04],
    'setSize': [10, 20],
    'matched_genes': ['Gene1;Gene2', 'Gene3;Gene4']
})
is_gsea = 'NES' in df.columns or 'nes' in df.columns
print(f'is_gsea: {is_gsea}')
assert is_gsea
print('✅ GSEA 图表兼容')
"`
Expected: `is_gsea: True`

### 9.4 ssGSEA/GSVA 分析图表

**Files:**
- Modify: `allenricher/visualization/gsva_plots.py` (如有需要)

**验证点**: ssGSEA/GSVA 结果包含 `pathway` 和 `NES`（或 `score`）列，WikiPathways 通路 ID（如 WP254）作为 pathway 名称即可。

- [ ] **Step 3: 验证 ssGSEA/GSVA 图表对 WikiPathways 的兼容性**

WikiPathways 通路 ID 格式为 `WP` + 数字（如 WP254），与现有通路 ID 格式兼容，无需特殊处理。

### 9.5 HTML 报告扩展

**Files:**
- Modify: `allenricher/report/generator.py`

**报告结构**:
```
HTML报告:
├── Summary (统计摘要 - 按数据库分组)
├── Plots (通用图表 - barplot/bubble for ORA)
├── GSEA Visualization (条件渲染 - 如果含 NES 列)
├── GSVA Visualization (条件渲染 - 如果含 score 列)
├── Tables (数据表格 - 动态列名检测)
│   ├── GO table
│   ├── KEGG table
│   ├── Reactome table
│   └── WikiPathways table  <-- 新增
└── AI Interpretation (可选)
```

- [ ] **Step 4: 验证报告表格动态列名检测**

`generator.py` 第320-362行通过以下逻辑动态检测:
```python
is_gsea = 'NES' in df.columns or 'nes' in df.columns
# GSEA: 显示 NES 列
# ORA: 显示 Rich Factor 列
```

WikiPathways 结果无论是 ORA 还是 GSEA 格式，都通过相同的列名检测逻辑，无需修改。

- [ ] **Step 5: 验证报告中的数据库表格**

`generator.py` 第311-378行的 `_generate_tables()` 方法遍历所有数据库结果生成表格:
```python
for db_name, df in results.items():
    # db_name 可以是 'GO', 'KEGG', 'Reactome', 'WikiPathways'
    # 每个数据库独立生成一个表格 section
```

WikiPathways 作为新数据库，会自动生成独立的表格 section，无需修改。

### 9.6 图表样式扩展（可选）

**Files:**
- Modify: `allenricher/visualization/plot_theme.py` (可选)

如需为 WikiPathways 添加专属配色，可在 `COLOR_PALETTES` 中添加:
```python
COLOR_PALETTES = {
    'GO': {...},
    'KEGG': {...},
    'Reactome': {...},
    'WikiPathways': {  # 新增
        'primary': '#2E7D32',    # 绿色系
        'secondary': '#81C784',
        'accent': '#1B5E20',
    },
}
```

---

## Task 10: API 功能扩展

### 10.1 Python API 使用示例

WikiPathways 接入后，Python API 使用方式与其他数据库完全一致:

```python
from allenricher import EnrichmentAnalyzer, Config
from allenricher.database.manager import DatabaseManager

# 1. 初始化配置
config = Config(
    species='hsa',
    databases=['GO', 'KEGG', 'WikiPathways'],  # 添加 WikiPathways
    method='hypergeometric',  # 或 'gsea', 'ssgsea', 'gsva'
)

# 2. 加载数据库
db_manager = DatabaseManager('./database', 'hsa')
db_manager.load_databases(['GO', 'KEGG', 'WikiPathways'])
database_data = db_manager.get_all_term_data()

# 3. 执行分析
analyzer = EnrichmentAnalyzer(config)
results = analyzer.run_analysis(
    gene_set={'Gene1', 'Gene2', ...},
    background_set=db_manager.get_background_genes(),
    database_data=database_data,
)

# 4. 结果包含 WikiPathways
print(results['WikiPathways'])  # DataFrame with Term_ID, Term_Name, P_Value, ...
```

### 10.2 CLI API 使用示例

```bash
# ORA 分析
allenricher analyze -i genes.txt -m hypergeometric -s hsa -d GO,KEGG,WikiPathways -o ./results

# GSEA 分析
allenricher analyze -i genes.txt -r ranked_genes.rnk -m gsea -s hsa -d WikiPathways -o ./results

# ssGSEA 分析
allenricher analyze -i expression.tsv -m ssgsea -s hsa -d WikiPathways -o ./results

# GSVA 分析
allenricher analyze -i expression.tsv -m gsva -s hsa -d WikiPathways -o ./results
```

### 10.3 报告 API 使用示例

```python
from allenricher.report.generator import ReportGenerator

# 生成包含 WikiPathways 的 HTML 报告
gen = ReportGenerator(output_dir='./results')
gen.generate_report(
    results={
        'GO': go_df,
        'KEGG': kegg_df,
        'WikiPathways': wp_df,  # WikiPathways 结果
    },
    method='hypergeometric',
    species='hsa',
)
```

---

## Task 11: 生成完整测试报告

**Files:**
- Create: `test_wikipathways_output/WIKIPATHWAYS_TEST_REPORT.md`

- [ ] **Step 1: 测试完成后生成报告**

报告内容应包含:
- 测试日期、环境信息
- 下载测试结果（物种数量、文件大小）
- 构建测试结果（生成的文件列表、通路数量、基因数量）
- **四种分析方法测试结果**:
  - ORA: 柱状图、气泡图、HTML表格
  - GSEA: 富集曲线、NES条形图、气泡图、HTML表格
  - ssGSEA: 热图、组间比较、气泡图、相关性热图
  - GSVA: 同上
- API 使用示例验证
- 已知限制和后续改进方向

---

## 自检清单

1. **Spec 覆盖**: 下载 ✓, 物种注册 ✓, 构建 ✓, 加载 ✓, 分析 ✓, CLI ✓, **图表 ✓**, **报告 ✓**, **API ✓**
2. **占位符扫描**: 无 TBD/TODO
3. **类型一致性**: `WikiPathways` 大写在 `DatabaseType` 枚举中已定义；`name_to_prefix` 使用 `'WIKIPATHWAYS': 'WikiPathways'`；文件前缀 `WikiPathways` 与枚举值一致
4. **数据流完整性**: GMT → parse_gmt() → build_database() → .tab.gz + .disc.gz → DatabaseManager._parse_tab_file() → {term_id: {name, genes}} → EnrichmentAnalyzer → **Plotter/ReportGenerator**
