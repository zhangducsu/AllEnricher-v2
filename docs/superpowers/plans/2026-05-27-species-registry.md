# 全物种支持名单与查询功能 Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** 在 download 阶段自动从 NCBI (gene2go/gene_info) 和 KEGG (list/organisms) 提取全物种支持名单，生成统一的物种支持表（含拉丁名、TaxID、物种缩写、缩写来源），并提供 `allenricher list` CLI 命令供用户按拉丁名或 TaxID 查询各数据库的物种支持情况。

**Architecture:** 新增 `SpeciesRegistry` 模块作为物种注册中心。在 `download` 阶段，解析已下载的 gene2go.gz/gene_info.gz 提取 GO/Reactome 支持物种，调用 KEGG API `list/organisms` 获取 KEGG 全物种列表，合并后生成 `supported_species.tsv.gz`。新增 `allenricher list` CLI 命令提供查询接口。对于不在 KEGG 中的物种，按 KEGG 惯例自动生成缩写（属名前3字母+种名首字母）。

**Tech Stack:** Python 3.10, requests (KEGG API), gzip/pandas (文件解析), argparse (CLI)

---

## 现状分析（已验证）

| 数据库 | 物种数据来源 | 当前物种发现方式 | 实际支持物种数 |
|--------|-------------|-----------------|--------------|
| GO (v1) | NCBI gene2go.gz (2019) | 无动态发现 | **45 物种** |
| GO (v2) | NCBI gene2go.gz (2026) | 无动态发现 | **2,332 物种** |
| Reactome (v1/v2) | NCBI2Reactome | 无动态发现 | **16 物种** |
| gene_info (v2) | NCBI gene_info.gz | - | **53,590 物种** |
| KEGG | REST API (build 阶段) | 硬编码 16 个物种 | ~8,000 物种 |
| DO | 仅人类 | 硬编码 | 1 |
| DisGeNET | 已失效 | - | 0 |

### 关键约束
1. KEGG 数据在 build 阶段才获取，但 物种名单可以在 download 阶段通过 `list/organisms` API 一次性获取
2. GO/Reactome 物种名单从已下载的 gene2go.gz 和 NCBI2Reactome 文件中提取
3. 用户查询时输入拉丁名或 TaxID，返回各数据库支持情况
4. 物种缩写来源需标注：KEGG 原生 vs 自动生成
5. **列名使用英文**: `species_code`（物种缩写）、`species_code_source`（物种缩写来源）

---

## 文件结构

| 文件 | 职责 |
|------|------|
| `allenricher/database/species_registry.py` | 物种注册中心：提取、合并、保存、查询物种名单 |
| `allenricher/database/kegg_fetcher.py` | 修改：新增 `fetch_organism_list()` 方法 |
| `allenricher/cli.py` | 修改：新增 `list` 子命令 |
| `tests/test_species_registry.py` | 物种注册中心单元测试 |
| `tests/test_cli_list.py` | CLI list 命令测试 |

---

## Task 1: 新增 KEGG 全物种列表获取

**Files:**
- Modify: `allenricher/database/kegg_fetcher.py`

- [ ] **Step 1: 在 KEGGFetcher 中新增 `fetch_organism_list()` 方法**

在 `KEGGFetcher` 类中添加以下方法：

```python
def fetch_organism_list(self) -> pd.DataFrame:
    """获取 KEGG 全物种列表

    调用 KEGG REST API list/organisms 获取所有支持的物种。
    返回 DataFrame 包含列: species_code, latin_name, taxonomy_id, full_name

    KEGG list/organisms 返回格式:
        T01001    hsa    Homo sapiens (human)    Eukaryota;Animals;Vertebrata;Mammals;...
        T01001    mmu    Mus musculus (mouse)    ...

    Returns:
        pd.DataFrame: 全物种列表
    """
    cache_file = self.cache_dir / "kegg_all_organisms.txt"

    if cache_file.exists() and not self.overwrite:
        lines = cache_file.read_text(encoding='utf-8').strip().split('\n')
    else:
        url = f"{BASE_URL}/list/organisms"
        resp = self._request(url)
        lines = resp.strip().split('\n')
        cache_file.write_text(resp, encoding='utf-8')

    records = []
    for line in lines:
        if not line.strip():
            continue
        parts = line.split('\t')
        if len(parts) < 3:
            continue
        taxonomy_id_str = parts[0].lstrip('T')  # 去掉 T 前缀
        species_code = parts[1]
        full_name = parts[2]  # 如 "Homo sapiens (human)"

        # 提取拉丁学名（去掉括号中的常用名）
        latin_name = full_name.split('(')[0].strip()

        try:
            taxonomy_id = int(taxonomy_id_str)
        except ValueError:
            continue

        records.append({
            'species_code': species_code,
            'latin_name': latin_name,
            'taxonomy_id': taxonomy_id,
            'full_name': full_name,
            'species_code_source': 'kegg',
        })

    return pd.DataFrame(records)
```

- [ ] **Step 2: 运行现有 KEGG 测试确认无回归**

Run: `cd f:\OneDrive\Documents\TraeSOLO\AllEnricher\AllEnricher-v2 && python -m pytest tests/test_download.py -v --tb=short -k kegg`
Expected: 全部通过

---

## Task 2: 创建 SpeciesRegistry 物种注册中心

**Files:**
- Create: `allenricher/database/species_registry.py`

- [ ] **Step 1: 编写 SpeciesRegistry 类**

```python
"""
物种注册中心模块

在 download 阶段从各数据源提取支持的物种名单，
合并为统一的物种支持表，并提供查询接口。
"""

import gzip
import os
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple

import pandas as pd


@dataclass
class SpeciesSupport:
    """单个物种的数据库支持情况"""
    taxonomy_id: int
    latin_name: str
    species_code: str = ""
    species_code_source: str = ""  # "kegg" 或 "auto"
    has_go: bool = False
    has_reactome: bool = False
    has_kegg: bool = False
    go_gene_count: int = 0
    reactome_pathway_count: int = 0
    kegg_pathway_count: int = 0


class SpeciesRegistry:
    """物种注册中心"""

    def __init__(self, root_dir: str = "./database"):
        self.root_dir = Path(root_dir)
        self.basic_dir = self.root_dir / "basic"
        self._species_data: Dict[int, SpeciesSupport] = {}

    # ========== 数据提取方法 ==========

    def extract_go_species(self, gene2go_path: str) -> Dict[int, Tuple[str, int]]:
        """从 gene2go.gz 提取 GO 支持的物种列表

        扫描 gene2go.gz 第一列（taxid），统计每个物种的注释基因数。
        同时从 gene_info.gz 获取拉丁学名。

        Args:
            gene2go_path: gene2go.gz 文件路径

        Returns:
            Dict[taxid, (latin_name, gene_count)]
        """
        taxid_genes: Dict[int, Set[str]] = {}
        opener = gzip.open if gene2go_path.endswith('.gz') else open

        with opener(gene2go_path, 'rt', encoding='utf-8') as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith('#'):
                    continue
                parts = line.split('\t')
                try:
                    taxid = int(parts[0])
                    gene_id = parts[1]
                except (ValueError, IndexError):
                    continue
                if taxid not in taxid_genes:
                    taxid_genes[taxid] = set()
                taxid_genes[taxid].add(gene_id)

        return {tid: (None, len(genes)) for tid, genes in taxid_genes.items()}

    def extract_reactome_species(self, ncbi2reactome_path: str) -> Dict[str, int]:
        """从 NCBI2Reactome 文件提取 Reactome 支持的物种列表

        通过 pathway_id (如 R-HSA-12345) 中的物种代码提取物种，
        统计每个物种的通路数量。

        Args:
            ncbi2reactome_path: NCBI2Reactome 文件路径

        Returns:
            Dict[species_code_upper, pathway_count]
        """
        species_pathways: Dict[str, Set[str]] = {}
        opener = gzip.open if ncbi2reactome_path.endswith('.gz') else open

        with opener(ncbi2reactome_path, 'rt', encoding='utf-8') as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith('#'):
                    continue
                parts = line.split('\t')
                if len(parts) < 2:
                    continue
                pathway_id = parts[1]
                pathway_parts = pathway_id.split('-')
                if len(pathway_parts) < 2:
                    continue
                species_code = pathway_parts[1]  # 如 "HSA"
                if species_code not in species_pathways:
                    species_pathways[species_code] = set()
                species_pathways[species_code].add(pathway_id)

        return {code: len(paths) for code, paths in species_pathways.items()}

    def extract_gene_info_names(self, gene_info_path: str) -> Dict[int, str]:
        """从 gene_info.gz 提取 taxid 到拉丁学名的映射

        Args:
            gene_info_path: gene_info.gz 文件路径

        Returns:
            Dict[taxid, latin_name]
        """
        taxid_to_name: Dict[int, str] = {}
        opener = gzip.open if gene_info_path.endswith('.gz') else open

        with opener(gene_info_path, 'rt', encoding='utf-8') as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith('#'):
                    continue
                parts = line.split('\t')
                try:
                    taxid = int(parts[0])
                    # latin_name 在 gene_info.gz 中的位置
                    # 格式: taxid\tGeneID\tSymbol\t...\tlatin_name
                    if len(parts) >= 13:
                        latin_name = parts[12]  # 通常在倒数几列
                    else:
                        latin_name = ""
                except (ValueError, IndexError):
                    continue
                if taxid not in taxid_to_name:
                    taxid_to_name[taxid] = latin_name

        return taxid_to_name

    # ========== 物种缩写自动生成 ==========

    @staticmethod
    def generate_kegg_abbreviation(latin_name: str) -> str:
        """按 KEGG 惯例自动生成物种缩写

        KEGG 惯例: 属名前3字母(小写) + 种名首字母(小写)
        例如: "Homo sapiens" -> "hsa", "Mus musculus" -> "mmu"

        Args:
            latin_name: 拉丁学名

        Returns:
            3-4 字母的小写缩写
        """
        parts = latin_name.strip().split()
        if len(parts) >= 2:
            genus = parts[0].lower()[:3]
            species = parts[1].lower()[0]
            return f"{genus}{species}"
        # 单个词的情况（少见）
        return latin_name.strip().lower()[:4]

    # ========== 合并与保存 ==========

    def build_registry(
        self,
        gene2go_path: Optional[str] = None,
        gene_info_path: Optional[str] = None,
        ncbi2reactome_path: Optional[str] = None,
        kegg_organisms_df: Optional[pd.DataFrame] = None,
    ) -> pd.DataFrame:
        """构建完整的物种支持表

        合并 GO、Reactome、KEGG 三个数据源的物种信息，
        生成统一的 DataFrame 并保存为 TSV 文件。

        Args:
            gene2go_path: gene2go.gz 路径
            gene_info_path: gene_info.gz 路径
            ncbi2reactome_path: NCBI2Reactome 路径
            kegg_organisms_df: KEGG 全物种 DataFrame (来自 KEGGFetcher.fetch_organism_list)

        Returns:
            pd.DataFrame: 统一物种支持表
        """
        # 1. 提取 GO 物种
        go_species = {}
        if gene2go_path and os.path.exists(gene2go_path):
            go_species = self.extract_go_species(gene2go_path)

        # 2. 提取 Reactome 物种
        reactome_species = {}
        if ncbi2reactome_path and os.path.exists(ncbi2reactome_path):
            reactome_species = self.extract_reactome_species(ncbi2reactome_path)

        # 3. 提取拉丁学名
        taxid_to_name = {}
        if gene_info_path and os.path.exists(gene_info_path):
            taxid_to_name = self.extract_gene_info_names(gene_info_path)

        # 4. 构建 KEGG 映射
        species_by_taxid: Dict[int, Tuple[str, str]] = {}  # taxid -> (species_code, source)
        species_by_code: Dict[str, Tuple[int, str]] = {}   # species_code -> (taxid, latin_name)
        if kegg_organisms_df is not None and len(kegg_organisms_df) > 0:
            for _, row in kegg_organisms_df.iterrows():
                species_by_taxid[row['taxonomy_id']] = (row['species_code'], 'kegg')
                species_by_code[row['species_code']] = (row['taxonomy_id'], row['latin_name'])

        # 5. Reactome 物种代码到 taxid 的映射
        # Reactome pathway_id 格式 R-HSA-12345，物种代码为 HSA
        # 需要通过 gene_info 或 KEGG 映射来关联 taxid
        reactome_code_to_taxid: Dict[str, int] = {}
        for code in reactome_species:
            # 尝试从 KEGG 映射中查找（代码转小写）
            species_code = code.lower()
            if species_code in species_by_code:
                reactome_code_to_taxid[code] = species_by_code[species_code][0]

        # 6. 合并所有物种
        all_taxids = set(go_species.keys()) | set(species_by_taxid.keys()) | set(reactome_code_to_taxid.values())

        records = []
        for taxid in sorted(all_taxids):
            latin_name = taxid_to_name.get(taxid, "")

            # KEGG 信息
            has_kegg = taxid in species_by_taxid
            species_code = ""
            species_code_source = ""
            if has_kegg:
                species_code, species_code_source = species_by_taxid[taxid]
            else:
                # 自动生成 物种缩写
                if latin_name:
                    species_code = self.generate_kegg_abbreviation(latin_name)
                    species_code_source = "auto"

            # GO 信息
            has_go = taxid in go_species
            go_gene_count = go_species.get(taxid, (None, 0))[1] if has_go else 0

            # Reactome 信息
            has_reactome = False
            reactome_pathway_count = 0
            for rcode, rtaxid in reactome_code_to_taxid.items():
                if rtaxid == taxid:
                    has_reactome = True
                    reactome_pathway_count = reactome_species.get(rcode, 0)
                    break

            records.append({
                'taxonomy_id': taxid,
                'latin_name': latin_name,
                'species_code': species_code,
                'species_code_source': species_code_source,
                'has_go': has_go,
                'has_reactome': has_reactome,
                'has_kegg': has_kegg,
                'go_gene_count': go_gene_count,
                'reactome_pathway_count': reactome_pathway_count,
                'supported_db_count': int(has_go) + int(has_reactome) + int(has_kegg),
            })

        df = pd.DataFrame(records)

        # 7. 保存
        output_dir = self.basic_dir / "species_registry"
        output_dir.mkdir(parents=True, exist_ok=True)
        output_path = output_dir / "supported_species.tsv.gz"
        df.to_csv(output_path, sep='\t', index=False, compression='gzip')

        self._species_data = {
            row['taxonomy_id']: SpeciesSupport(
                taxonomy_id=row['taxonomy_id'],
                latin_name=row['latin_name'],
                species_code=row['species_code'],
                species_code_source=row['species_code_source'],
                has_go=row['has_go'],
                has_reactome=row['has_reactome'],
                has_kegg=row['has_kegg'],
                go_gene_count=row['go_gene_count'],
                reactome_pathway_count=row['reactome_pathway_count'],
            )
            for _, row in df.iterrows()
        }

        return df

    # ========== 查询方法 ==========

    def load_registry(self) -> Optional[pd.DataFrame]:
        """从已保存的文件加载物种注册表

        Returns:
            pd.DataFrame 或 None（如果文件不存在）
        """
        registry_path = self.basic_dir / "species_registry" / "supported_species.tsv.gz"
        if not registry_path.exists():
            return None
        df = pd.read_csv(registry_path, sep='\t', compression='gzip')
        self._species_data = {
            row['taxonomy_id']: SpeciesSupport(
                taxonomy_id=row['taxonomy_id'],
                latin_name=row['latin_name'],
                species_code=row['species_code'],
                species_code_source=row['species_code_source'],
                has_go=bool(row['has_go']),
                has_reactome=bool(row['has_reactome']),
                has_kegg=bool(row['has_kegg']),
                go_gene_count=row['go_gene_count'],
                reactome_pathway_count=row['reactome_pathway_count'],
            )
            for _, row in df.iterrows()
        }
        return df

    def query_by_taxid(self, taxonomy_id: int) -> Optional[SpeciesSupport]:
        """通过 TaxID 查询物种支持情况"""
        if not self._species_data:
            self.load_registry()
        return self._species_data.get(taxonomy_id)

    def query_by_latin_name(self, latin_name: str) -> List[SpeciesSupport]:
        """通过拉丁学名查询（模糊匹配）"""
        if not self._species_data:
            self.load_registry()
        latin_lower = latin_name.lower()
        return [
            info for info in self._species_data.values()
            if latin_lower in info.latin_name.lower()
        ]

    def query_by_species_code(self, species_code: str) -> Optional[SpeciesSupport]:
        """通过 物种代码查询"""
        if not self._species_data:
            self.load_registry()
        code_lower = species_code.lower()
        for info in self._species_data.values():
            if info.species_code.lower() == code_lower:
                return info
        return None

    def get_summary(self) -> Dict[str, int]:
        """获取物种支持统计摘要

        Returns:
            Dict 包含: total_species, go_species, reactome_species, kegg_species
        """
        if not self._species_data:
            self.load_registry()
        return {
            'total_species': len(self._species_data),
            'go_species': sum(1 for s in self._species_data.values() if s.has_go),
            'reactome_species': sum(1 for s in self._species_data.values() if s.has_reactome),
            'kegg_species': sum(1 for s in self._species_data.values() if s.has_kegg),
        }
```

- [ ] **Step 2: 运行测试确认模块可导入**

Run: `cd f:\OneDrive\Documents\TraeSOLO\AllEnricher\AllEnricher-v2 && python -c "from allenricher.database.species_registry import SpeciesRegistry; print('OK')"`
Expected: OK

---

## Task 3: 集成到 download 流程

**Files:**
- Modify: `allenricher/database/downloader.py`

- [ ] **Step 1: 在 DataDownloader.download_all() 末尾调用 SpeciesRegistry**

在 `download_all()` 方法的返回前，添加物种注册表构建逻辑：

```python
# 在 download_all() 方法末尾、return results 之前添加：

# 构建物种支持名单
try:
    from .species_registry import SpeciesRegistry
    from .kegg_fetcher import KEGGFetcher

    registry = SpeciesRegistry(root_dir=self.root_dir)

    # 获取 KEGG 全物种列表
    species_fetcher = KEGGFetcher(
        cache_dir=str(self.basic_dir / "kegg"),
        overwrite=self.overwrite
    )
    species_df = species_fetcher.fetch_organism_list()

    # 定位已下载的文件路径
    gene2go_path = None
    gene_info_path = None
    ncbi2reactome_path = None

    if 'go' in downloaded:
        go_dir = Path(downloaded['go'])
        gene2go_path = str(go_dir / "gene2go.gz")
        gene_info_path = str(go_dir / "gene_info.gz")

    if 'reactome' in downloaded:
        reactome_dir = Path(downloaded['reactome'])
        ncbi2reactome_path = str(reactome_dir / "NCBI2Reactome_All_Levels.txt.gz")
        if not gene_info_path:
            gene_info_path = str(reactome_dir / "gene_info.gz")

    # 构建并保存
    species_df = registry.build_registry(
        gene2go_path=gene2go_path,
        gene_info_path=gene_info_path,
        ncbi2reactome_path=ncbi2reactome_path,
        kegg_organisms_df=species_df,
    )
    summary = registry.get_summary()
    print(f"\n物种支持名单已生成:")
    print(f"  总物种数: {summary['total_species']}")
    print(f"  GO 支持: {summary['go_species']} 物种")
    print(f"  Reactome 支持: {summary['reactome_species']} 物种")
    print(f"  KEGG 支持: {summary['kegg_species']} 物种")
    print(f"  保存路径: {self.basic_dir / 'species_registry' / 'supported_species.tsv.gz'}")
except Exception as e:
    print(f"警告: 物种支持名单生成失败: {e}")
```

- [ ] **Step 2: 运行 download 测试确认无回归**

Run: `cd f:\OneDrive\Documents\TraeSOLO\AllEnricher\AllEnricher-v2 && python -m pytest tests/test_download.py -v --tb=short`
Expected: 全部通过

---

## Task 4: 新增 `allenricher list` CLI 命令

**Files:**
- Modify: `allenricher/cli.py`

- [ ] **Step 1: 添加 list 子命令参数定义**

在 CLI 的 subparser 定义区域添加：

```python
# list 子命令
list_parser = subparsers.add_parser('list', help='查询支持的物种及数据库覆盖情况')
list_parser.add_argument('--taxid', type=int, help='按 Taxonomy ID 查询')
list_parser.add_argument('--name', type=str, help='按拉丁学名查询（支持模糊匹配）')
list_parser.add_argument('--code', type=str, help='按 物种代码查询')
list_parser.add_argument('--database-dir', default='./database', help='数据库目录')
list_parser.add_argument('--format', choices=['table', 'json', 'tsv'], default='table',
                        help='输出格式')
list_parser.add_argument('--all', action='store_true', help='列出所有支持的物种')
```

- [ ] **Step 2: 添加 cmd_list() 处理函数**

```python
def cmd_list(args):
    """处理 list 命令"""
    from allenricher.database.species_registry import SpeciesRegistry

    registry = SpeciesRegistry(root_dir=args.database_dir)
    df = registry.load_registry()

    if df is None or len(df) == 0:
        print("物种支持名单尚未生成。请先运行: allenricher download -d go,reactome")
        return 1

    # 查询模式
    if args.taxid is not None:
        result = registry.query_by_taxid(args.taxid)
        if result is None:
            print(f"未找到 TaxID={args.taxid} 的物种")
            return 1
        _print_species_detail(result, args.format)
        return 0

    if args.name is not None:
        results = registry.query_by_latin_name(args.name)
        if not results:
            print(f"未找到匹配 '{args.name}' 的物种")
            return 1
        for r in results:
            _print_species_detail(r, args.format)
        return 0

    if args.code is not None:
        result = registry.query_by_species_code(args.code)
        if result is None:
            print(f"未找到 物种代码 '{args.code}' 的物种")
            return 1
        _print_species_detail(result, args.format)
        return 0

    # 默认: 显示统计摘要
    summary = registry.get_summary()
    print(f"=== AllEnricher 物种支持统计 ===")
    print(f"总物种数: {summary['total_species']}")
    print(f"GO 支持: {summary['go_species']} 物种")
    print(f"Reactome 支持: {summary['reactome_species']} 物种")
    print(f"KEGG 支持: {summary['kegg_species']} 物种")
    print()
    print("使用 --taxid/--name/--code 查询具体物种")
    print("使用 --all 列出所有物种")

    if args.all:
        print(f"\n=== 全部物种列表 ({len(df)} 物种) ===")
        if args.format == 'json':
            print(df.to_json(orient='records', indent=2, force_ascii=False))
        elif args.format == 'tsv':
            print(df.to_csv(sep='\t', index=False))
        else:
            # 表格格式，只显示关键列
            display_cols = ['taxonomy_id', 'latin_name', 'species_code', 'species_code_source',
                           'has_go', 'has_reactome', 'has_kegg']
            print(df[display_cols].to_string(index=False))

    return 0


def _print_species_detail(species, fmt='table'):
    """打印单个物种的详细信息"""
    if fmt == 'json':
        import json
        print(json.dumps({
            'taxonomy_id': species.taxonomy_id,
            'latin_name': species.latin_name,
            'species_code': species.species_code,
            'species_code_source': species.species_code_source,
            'has_go': species.has_go,
            'has_reactome': species.has_reactome,
            'has_kegg': species.has_kegg,
            'go_gene_count': species.go_gene_count,
            'reactome_pathway_count': species.reactome_pathway_count,
        }, indent=2, ensure_ascii=False))
    else:
        print(f"=== 物种信息 ===")
        print(f"  拉丁学名: {species.latin_name}")
        print(f"  Taxonomy ID: {species.taxonomy_id}")
        print(f"  物种缩写: {species.species_code}")
        print(f"  缩写来源: {species.species_code_source}")
        print(f"  数据库支持:")
        print(f"    GO:       {'✓' if species.has_go else '✗'} ({species.go_gene_count} 基因)")
        print(f"    Reactome: {'✓' if species.has_reactome else '✗'} ({species.reactome_pathway_count} 通路)")
        print(f"    KEGG:     {'✓' if species.has_kegg else '✗'}")
```

- [ ] **Step 3: 在 main() 中注册 list 子命令**

在 `main()` 函数中添加：
```python
list_parser.set_defaults(func=cmd_list)
```

---

## Task 5: 编写测试

**Files:**
- Create: `tests/test_species_registry.py`

- [ ] **Step 1: 编写 SpeciesRegistry 单元测试**

```python
"""SpeciesRegistry 单元测试"""

import gzip
import os
import tempfile

import pandas as pd
import pytest

from allenricher.database.species_registry import SpeciesRegistry, SpeciesSupport


class TestExtractGoSpecies:
    """GO 物种提取测试"""

    def test_extract_from_gene2go(self, tmp_path):
        """从 gene2go.gz 提取物种列表"""
        # 构造测试数据
        gene2go_path = tmp_path / "gene2go.gz"
        with gzip.open(gene2go_path, 'wt') as f:
            f.write("9606\t1\tGO:0000001\tIDA\t\tbiological_process\t\t\tBP\n")
            f.write("9606\t2\tGO:0000002\tIDA\t\tbiological_process\t\t\tBP\n")
            f.write("9606\t1\tGO:0000003\tIDA\t\tmolecular_function\t\t\tMF\n")
            f.write("10090\t3\tGO:0000001\tIDA\t\tbiological_process\t\t\tBP\n")
            f.write("10090\t4\tGO:0000004\tIDA\t\tbiological_process\t\t\tBP\n")

        registry = SpeciesRegistry()
        result = registry.extract_go_species(str(gene2go_path))

        assert 9606 in result
        assert 10090 in result
        assert result[9606][1] == 2  # 2 个唯一基因
        assert result[10090][1] == 2  # 2 个唯一基因

    def test_empty_file(self, tmp_path):
        """空文件返回空字典"""
        gene2go_path = tmp_path / "gene2go.gz"
        with gzip.open(gene2go_path, 'wt') as f:
            f.write("")

        registry = SpeciesRegistry()
        result = registry.extract_go_species(str(gene2go_path))
        assert result == {}


class TestExtractReactomeSpecies:
    """Reactome 物种提取测试"""

    def test_extract_from_ncbi2reactome(self, tmp_path):
        """从 NCBI2Reactome 提取物种列表"""
        ncbi2reactome_path = tmp_path / "NCBI2Reactome.txt.gz"
        with gzip.open(ncbi2reactome_path, 'wt') as f:
            f.write("1\tR-HSA-12345\t\tCell Cycle\t\n")
            f.write("2\tR-HSA-12345\t\tCell Cycle\t\n")
            f.write("3\tR-HSA-67890\t\tApoptosis\t\n")
            f.write("4\tR-MMU-11111\t\tMetabolism\t\n")

        registry = SpeciesRegistry()
        result = registry.extract_reactome_species(str(ncbi2reactome_path))

        assert "HSA" in result
        assert "MMU" in result
        assert result["HSA"] == 2  # 2 个唯一通路
        assert result["MMU"] == 1


class TestKeggAbbreviation:
    """物种缩写自动生成测试"""

    def test_standard_binomial(self):
        """标准双名法"""
        assert SpeciesRegistry.generate_kegg_abbreviation("Homo sapiens") == "hsa"
        assert SpeciesRegistry.generate_kegg_abbreviation("Mus musculus") == "mmu"
        assert SpeciesRegistry.generate_kegg_abbreviation("Rattus norvegicus") == "rno"

    def test_three_letter_genus(self):
        """属名恰好3个字母"""
        assert SpeciesRegistry.generate_kegg_abbreviation("Sus scrofa") == "ssc"

    def test_long_genus(self):
        """属名超过3个字母"""
        assert SpeciesRegistry.generate_kegg_abbreviation("Danio rerio") == "dan"
        assert SpeciesRegistry.generate_kegg_abbreviation("Xenopus tropicalis") == "xet"


class TestBuildRegistry:
    """注册表构建测试"""

    def test_build_with_all_sources(self, tmp_path):
        """使用所有数据源构建注册表"""
        # 构造 gene2go
        gene2go_path = tmp_path / "gene2go.gz"
        with gzip.open(gene2go_path, 'wt') as f:
            f.write("9606\t1\tGO:0000001\tIDA\t\tBP\t\t\tBP\n")
            f.write("10090\t2\tGO:0000002\tIDA\t\tBP\t\t\tBP\n")

        # 构造 KEGG DataFrame
        species_df = pd.DataFrame([
            {'species_code': 'hsa', 'latin_name': 'Homo sapiens', 'taxonomy_id': 9606,
             'full_name': 'Homo sapiens (human)', 'species_code_source': 'kegg'},
            {'species_code': 'mmu', 'latin_name': 'Mus musculus', 'taxonomy_id': 10090,
             'full_name': 'Mus musculus (mouse)', 'species_code_source': 'kegg'},
        ])

        registry = SpeciesRegistry(root_dir=str(tmp_path / "db"))
        df = registry.build_registry(
            gene2go_path=str(gene2go_path),
            kegg_organisms_df=species_df,
        )

        assert len(df) == 2
        assert df[df['taxonomy_id'] == 9606].iloc[0]['has_go'] == True
        assert df[df['taxonomy_id'] == 9606].iloc[0]['has_kegg'] == True
        assert df[df['taxonomy_id'] == 9606].iloc[0]['species_code'] == 'hsa'
        assert df[df['taxonomy_id'] == 9606].iloc[0]['species_code_source'] == 'kegg'

    def test_auto_species_code_for_non_kegg_species(self, tmp_path):
        """非 物种自动生成缩写"""
        gene2go_path = tmp_path / "gene2go.gz"
        with gzip.open(gene2go_path, 'wt') as f:
            f.write("3702\t5\tGO:0000001\tIDA\t\tBP\t\t\tBP\n")  # 拟南芥

        # gene_info 提供拉丁名
        gene_info_path = tmp_path / "gene_info.gz"
        with gzip.open(gene_info_path, 'wt') as f:
            f.write("3702\t5\tAG1\t\t\t\t\t\t\t\t\tArabidopsis thaliana\t\n")

        registry = SpeciesRegistry(root_dir=str(tmp_path / "db"))
        df = registry.build_registry(
            gene2go_path=str(gene2go_path),
            gene_info_path=str(gene_info_path),
        )

        ara_row = df[df['taxonomy_id'] == 3702]
        assert len(ara_row) == 1
        assert ara_row.iloc[0]['has_go'] == True
        assert ara_row.iloc[0]['has_kegg'] == False
        assert ara_row.iloc[0]['species_code'] == 'art'  # 自动生成
        assert ara_row.iloc[0]['species_code_source'] == 'auto'


class TestQueryMethods:
    """查询方法测试"""

    def test_query_by_taxid(self, tmp_path):
        """通过 TaxID 查询"""
        registry_path = tmp_path / "db" / "basic" / "species_registry"
        registry_path.mkdir(parents=True)
        df = pd.DataFrame([
            {'taxonomy_id': 9606, 'latin_name': 'Homo sapiens', 'species_code': 'hsa',
             'species_code_source': 'kegg', 'has_go': True, 'has_reactome': True,
             'has_kegg': True, 'go_gene_count': 100, 'reactome_pathway_count': 50},
        ])
        df.to_csv(registry_path / "supported_species.tsv.gz", sep='\t', index=False, compression='gzip')

        registry = SpeciesRegistry(root_dir=str(tmp_path / "db"))
        result = registry.query_by_taxid(9606)
        assert result is not None
        assert result.latin_name == "Homo sapiens"
        assert result.has_go == True

    def test_query_by_latin_name(self, tmp_path):
        """通过拉丁学名模糊查询"""
        registry_path = tmp_path / "db" / "basic" / "species_registry"
        registry_path.mkdir(parents=True)
        df = pd.DataFrame([
            {'taxonomy_id': 9606, 'latin_name': 'Homo sapiens', 'species_code': 'hsa',
             'species_code_source': 'kegg', 'has_go': True, 'has_reactome': True,
             'has_kegg': True, 'go_gene_count': 100, 'reactome_pathway_count': 50},
            {'taxonomy_id': 10090, 'latin_name': 'Mus musculus', 'species_code': 'mmu',
             'species_code_source': 'kegg', 'has_go': True, 'has_reactome': False,
             'has_kegg': True, 'go_gene_count': 80, 'reactome_pathway_count': 0},
        ])
        df.to_csv(registry_path / "supported_species.tsv.gz", sep='\t', index=False, compression='gzip')

        registry = SpeciesRegistry(root_dir=str(tmp_path / "db"))
        results = registry.query_by_latin_name("homo")
        assert len(results) == 1
        assert results[0].taxonomy_id == 9606

    def test_query_not_found(self, tmp_path):
        """查询不存在的物种"""
        registry = SpeciesRegistry(root_dir=str(tmp_path / "db"))
        result = registry.query_by_taxid(999999)
        assert result is None
```

- [ ] **Step 2: 运行测试**

Run: `cd f:\OneDrive\Documents\TraeSOLO\AllEnricher\AllEnricher-v2 && python -m pytest tests/test_species_registry.py -v --tb=short`
Expected: 全部通过

---

## Task 6: 端对端测试（使用已下载的真实数据）

**Files:**
- Create: `tests/test_species_registry_e2e.py`

- [ ] **Step 1: 编写 E2E 测试**

```python
"""SpeciesRegistry 端对端测试

使用项目已下载的真实数据库文件测试物种注册表构建。
"""

import os

import pandas as pd
import pytest

from allenricher.database.species_registry import SpeciesRegistry

TEST_DATA_DIR = os.path.join(os.path.dirname(__file__), "..", "test_data")
DATABASE_DIR = os.path.join(os.path.dirname(__file__), "..", "database")


class TestSpeciesRegistryE2E:
    """使用真实数据的端对端测试"""

    @pytest.fixture(scope="class")
    def registry(self):
        """创建 SpeciesRegistry 实例"""
        reg = SpeciesRegistry(root_dir=DATABASE_DIR)
        df = reg.load_registry()
        # 如果注册表不存在，尝试从已下载的文件构建
        if df is None:
            # 查找已下载的文件
            basic_dir = os.path.join(DATABASE_DIR, "basic")
            gene2go_path = None
            gene_info_path = None
            ncbi2reactome_path = None

            # 查找 GO 文件
            for root, dirs, files in os.walk(basic_dir):
                for f in files:
                    if f == "gene2go.gz":
                        gene2go_path = os.path.join(root, f)
                    elif f == "gene_info.gz":
                        gene_info_path = os.path.join(root, f)
                    elif f == "NCBI2Reactome_All_Levels.txt.gz":
                        ncbi2reactome_path = os.path.join(root, f)

            if gene2go_path:
                df = reg.build_registry(
                    gene2go_path=gene2go_path,
                    gene_info_path=gene_info_path,
                    ncbi2reactome_path=ncbi2reactome_path,
                )
        yield reg, df

    def test_registry_has_data(self, registry):
        """注册表有数据"""
        reg, df = registry
        if df is None:
            pytest.skip("无已下载数据库文件，跳过 E2E 测试")
        assert len(df) > 0

    def test_go_species_count(self, registry):
        """GO 物种数量合理（预期 > 100）"""
        reg, df = registry
        if df is None:
            pytest.skip("无已下载数据库文件，跳过 E2E 测试")
        go_count = df['has_go'].sum()
        assert go_count > 100, f"GO 物种数量过少: {go_count}"

    def test_query_human(self, registry):
        """查询人类（TaxID=9606）"""
        reg, df = registry
        if df is None:
            pytest.skip("无已下载数据库文件，跳过 E2E 测试")
        result = reg.query_by_taxid(9606)
        if result is None:
            pytest.skip("人类不在注册表中")
        assert result.has_go == True
        assert "sapiens" in result.latin_name.lower()

    def test_query_mouse(self, registry):
        """查询小鼠（TaxID=10090）"""
        reg, df = registry
        if df is None:
            pytest.skip("无已下载数据库文件，跳过 E2E 测试")
        result = reg.query_by_taxid(10090)
        if result is None:
            pytest.skip("小鼠不在注册表中")
        assert result.has_go == True

    def test_species_code_source_labeled(self, registry):
        """物种缩写来源已标注"""
        reg, df = registry
        if df is None:
            pytest.skip("无已下载数据库文件，跳过 E2E 测试")
        assert 'species_code_source' in df.columns
        sources = df['species_code_source'].unique()
        # 至少应该有 'auto' 或 'kegg'
        assert len(sources) > 0

    def test_summary(self, registry):
        """统计摘要正确"""
        reg, df = registry
        if df is None:
            pytest.skip("无已下载数据库文件，跳过 E2E 测试")
        summary = reg.get_summary()
        assert summary['total_species'] > 0
        assert summary['go_species'] > 0
```

- [ ] **Step 2: 运行 E2E 测试**

Run: `cd f:\OneDrive\Documents\TraeSOLO\AllEnricher\AllEnricher-v2 && python -m pytest tests/test_species_registry_e2e.py -v --tb=short`
Expected: 通过（如果已有下载数据）或 skip

---

## Task 7: 更新 __init__.py 导出和 README

**Files:**
- Modify: `allenricher/database/__init__.py`
- Modify: `README.md`

- [ ] **Step 1: 更新 __init__.py 导出**

在 `allenricher/database/__init__.py` 中添加：
```python
from .species_registry import SpeciesRegistry, SpeciesSupport
```

并更新 `__all__` 列表。

- [ ] **Step 2: 更新 README.md**

在 README.md 中添加 `list` 命令的文档说明，包括：
- 命令用法示例
- 查询输出格式说明
- 物种支持表字段说明

---

## Self-Review 检查清单

1. **Spec 覆盖**: 
   - ✅ 下载全物种数据库时整理物种名单表 → Task 3
   - ✅ KEGG 在 download 阶段获取物种名单 → Task 1
   - ✅ 自动生成 物种缩写 → Task 2 (generate_kegg_abbreviation)
   - ✅ 缩写来源标注 → Task 2 (species_code_source 字段)
   - ✅ 按拉丁名或 TaxID 查询 → Task 2 (query_by_taxid/query_by_latin_name)
   - ✅ 查询返回各数据库支持数量、TaxID、拉丁名、物种缩写、缩写来源 → Task 2 (SpeciesSupport)

2. **Placeholder 扫描**: 无 TBD/TODO

3. **类型一致性**: SpeciesSupport 在 Task 2 定义，Task 5 测试和 Task 6 E2E 中使用一致
