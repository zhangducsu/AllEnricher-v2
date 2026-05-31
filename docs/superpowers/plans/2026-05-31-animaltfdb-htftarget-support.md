# AnimalTFDB + hTFtarget 多物种 TF 富集分析 实现计划

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** 为 AllEnricher v2 添加 AnimalTFDB 4.0（183个动物物种）和 hTFtarget（人类TF-target关系）数据库支持，通过同源映射策略将人类TF-target关系推断到非模式物种，实现多物种TF富集分析、可视化图表和HTML报告。

**Architecture:**
- hTFtarget 提供人类 TF→靶基因关系（~134万条，ChIP-seq实验验证），作为映射源
- AnimalTFDB 4.0 提供 183 个动物物种的 TF 列表和人类直系同源映射，作为桥梁
- 同源映射流程：人类TF-target → 通过直系同源映射 → 目标物种TF-target
- 与现有 TRRUST/ChEA3 架构一致：Fetcher → Parser → Builder → Manager → Analyzer → Visualizer → CLI
- 复用现有 TFEnrichmentAnalyzer 进行 ORA/GSEA/ssGSEA/GSVA 分析

**Tech Stack:** Python 3.x, requests (下载), pandas (解析), numpy, plotly (可视化), Jinja2 (HTML报告), scipy/statsmodels (统计检验)

---

## 前置调研：数据源特性

### hTFtarget（人类TF-target关系）

| 属性 | 详情 |
|------|------|
| **URL** | https://guolab.wchscu.cn/hTFtarget/ |
| **下载URL** | `https://guolab.wchscu.cn/static/hTFtarget/file_download/tf-target-infomation.txt` |
| **文件大小** | ~56 MB |
| **数据量** | ~1,342,129 条 TF-target 关系 |
| **格式** | TSV (TF, target, tissue) |
| **数据来源** | ENCODE/SRA ChIP-Seq，30%阈值过滤（至少3个数据集） |
| **物种** | 仅人类 |

### AnimalTFDB 4.0（动物TF分类注释）

| 属性 | 详情 |
|------|------|
| **URL** | https://guolab.wchscu.cn/AnimalTFDB4/ |
| **物种覆盖** | 183 个动物物种 |
| **TF数量** | 274,633 个TF，73个家族 |
| **下载文件** | TF列表(TSV)、直系同源映射(TSV)、蛋白序列(FASTA) |

**关键下载URL模式：**
```
# TF列表
https://guolab.wchscu.cn/AnimalTFDB4_static/download/TF_list_final/{Species}_TF

# 直系同源映射（到人类）
https://guolab.wchscu.cn/AnimalTFDB4_static/download/ortholog_to_human_download/{Species}_ortholog_to_human
```

**TF列表列结构：** `Species, Symbol, Ensembl, Family, Protein, Entrez_ID`
**直系同源列结构：** 需下载确认，预期为 `{Species_Gene} \t {Human_Gene}`

### 同源映射策略

```
人类 TF → 靶基因 (hTFtarget)
    ↓ 直系同源映射 (AnimalTFDB ortholog_to_human)
目标物种 TF → 推断的靶基因
```

**映射规则：**
1. 从 AnimalTFDB 下载目标物种的 `ortholog_to_human` 文件
2. 构建 `{目标物种基因} → {人类基因}` 映射表
3. 从 hTFtarget 获取 `{人类TF} → {人类靶基因}` 关系
4. 反向映射：`{人类靶基因} → {目标物种基因}`
5. 合并得到 `{目标物种TF} → {目标物种靶基因}` 关系

---

## 与现有数据库对比

| 数据库 | 数据类型 | 物种数 | 核心内容 | 准确性 |
|--------|----------|--------|----------|--------|
| TRRUST | TF-target（实验） | 2 (人/鼠) | 调控网络(激活/抑制) | 高 |
| ChEA3 | TF-target（ChIP-seq整合） | 1 (人) | 6库整合 | 高 |
| **hTFtarget** | TF-target（ChIP-seq） | 1 (人) | ~134万条关系 | 高 |
| **AnimalTFDB** | TF注释+同源映射 | 183 (动物) | TF分类+直系同源 | 中-高 |

**AnimalTFDB独特价值：**
- 通过同源映射，将人类实验验证的TF-target关系扩展到183个动物物种
- 覆盖家畜（牛/猪/羊/鸡/狗/马等）、模式生物（斑马鱼/果蝇/线虫等）
- 与TRRUST/ChEA3互补：TRRUST/ChEA3直接用于人/鼠，AnimalTFDB扩展到其他物种

---

## 文件变更清单

| 操作 | 文件路径 | 职责 |
|------|---------|------|
| 新建 | `allenricher/database/htftarget_fetcher.py` | hTFtarget 数据下载器 |
| 新建 | `allenricher/database/animaltfdb_fetcher.py` | AnimalTFDB 数据下载器 |
| 新建 | `allenricher/database/parsers/htftarget.py` | hTFtarget TSV 解析器 |
| 新建 | `allenricher/database/parsers/animaltfdb.py` | AnimalTFDB TF列表+同源映射解析器 |
| 新建 | `allenricher/database/ortholog_mapper.py` | 同源映射引擎（核心） |
| 修改 | `allenricher/database/builder.py` | 添加 build_animaltfdb 方法 |
| 修改 | `allenricher/database/manager.py` | 添加 load_animaltfdb 方法 |
| 修改 | `allenricher/database/species_registry.py` | 添加 AnimalTFDB 字段 |
| 修改 | `allenricher/database/parsers/__init__.py` | 导出新 Parser |
| 修改 | `allenricher/cli.py` | 添加 AnimalTFDB 相关命令 |
| 修改 | `allenricher/report/templates/tf_report.html` | 添加映射来源标注 |

---

## Task 1: hTFtarget 数据下载器

**Files:**
- Create: `allenricher/database/htftarget_fetcher.py`

- [ ] **Step 1: 创建 hTFtarget 下载器类**

```python
"""
hTFtarget 数据下载器

从 guolab.wchscu.cn 下载人类转录因子-靶基因关系数据。

hTFtarget 特性：
- 基于 ENCODE/SRA 的 ChIP-Seq 数据
- 659 个人类 TF 的 ~134 万条 TF-target 关系
- 30% 阈值过滤（至少 3 个数据集）
- 包含组织来源信息

数据源：
- hTFtarget: https://guolab.wchscu.cn/hTFtarget/
"""

from pathlib import Path
from typing import Dict, List, Optional
import requests
import logging

logger = logging.getLogger(__name__)

# hTFtarget 下载链接
HTFTARGET_DOWNLOAD_URL = (
    "https://guolab.wchscu.cn/static/hTFtarget/file_download/tf-target-infomation.txt"
)


class HTFtargetFetcher:
    """hTFtarget 数据下载器

    下载人类转录因子-靶基因关系数据。

    Usage::

        fetcher = HTFtargetFetcher(basic_dir='./database/basic')
        fetcher.download()
    """

    REQUEST_TIMEOUT = 300  # 56MB 文件需要较长时间

    def __init__(self, basic_dir: str):
        """
        Args:
            basic_dir: 基础缓存目录
        """
        self.basic_dir = Path(basic_dir)
        self.basic_dir.mkdir(parents=True, exist_ok=True)

    def _get_cache_dir(self) -> Path:
        """获取缓存目录"""
        cache_dir = self.basic_dir / "htftarget"
        cache_dir.mkdir(parents=True, exist_ok=True)
        return cache_dir

    def download(self, overwrite: bool = False) -> Path:
        """下载 hTFtarget TF-target 关系文件

        Args:
            overwrite: 是否覆盖已存在的文件

        Returns:
            下载的 TSV 文件路径
        """
        cache_dir = self._get_cache_dir()
        local_path = cache_dir / "tf-target-information.txt"

        if local_path.exists() and not overwrite:
            logger.info(f"hTFtarget 已缓存，跳过: {local_path}")
            return local_path

        logger.info(f"下载 hTFtarget: {HTFTARGET_DOWNLOAD_URL}")

        try:
            resp = requests.get(
                HTFTARGET_DOWNLOAD_URL,
                timeout=self.REQUEST_TIMEOUT,
                stream=True,
                headers={"User-Agent": "Mozilla/5.0"},
            )
            resp.raise_for_status()
        except requests.RequestException as e:
            raise RuntimeError(f"hTFtarget 下载失败: {e}") from e

        # 流式写入
        with open(local_path, "wb") as f:
            for chunk in resp.iter_content(chunk_size=8192):
                if chunk:
                    f.write(chunk)

        size_mb = local_path.stat().st_size / 1024 / 1024
        logger.info(f"hTFtarget 已保存: {local_path} ({size_mb:.1f} MB)")
        return local_path

    @staticmethod
    def get_info() -> Dict[str, str]:
        """获取数据库基本信息"""
        return {
            "name": "hTFtarget",
            "version": "2020",
            "url": "https://guolab.wchscu.cn/hTFtarget/",
            "species": "Homo sapiens",
            "description": "人类TF-target关系，基于ChIP-Seq",
        }
```

- [ ] **Step 2: 验证下载器**

Run: `python -c "from allenricher.database.htftarget_fetcher import HTFtargetFetcher; print('HTFtargetFetcher OK')"`
Expected: `HTFtargetFetcher OK`

---

## Task 2: AnimalTFDB 数据下载器

**Files:**
- Create: `allenricher/database/animaltfdb_fetcher.py`

- [ ] **Step 1: 创建 AnimalTFDB 下载器类**

```python
"""
AnimalTFDB 4.0 数据下载器

从 guolab.wchscu.cn/AnimalTFDB4/ 下载动物转录因子注释数据。

AnimalTFDB 4.0 特性：
- 183 个动物物种
- 274,633 个 TF，73 个家族
- TF 列表、直系同源映射、蛋白序列

数据源：
- AnimalTFDB 4.0: https://guolab.wchscu.cn/AnimalTFDB4/
"""

from pathlib import Path
from typing import Dict, List, Optional, Set
import requests
import logging

logger = logging.getLogger(__name__)

# AnimalTFDB 下载基础 URL
ANIMALTFDB_BASE_URL = (
    "https://guolab.wchscu.cn/AnimalTFDB4_static/download"
)

# 主要模式生物和经济动物（优先支持）
ANIMALTFDB_PRIORITY_SPECIES = [
    "Homo_sapiens", "Mus_musculus", "Rattus_norvegicus",
    "Danio_rerio", "Drosophila_melanogaster", "Caenorhabditis_elegans",
    "Bos_taurus", "Sus_scrofa", "Ovis_aries", "Capra_hircus",
    "Gallus_gallus", "Canis_lupus_familiaris", "Equus_caballus",
    "Felis_catus", "Macaca_mulatta", "Gorilla_gorilla_gorilla",
    "Pan_troglodytes", "Oryctolagus_cuniculus", "Xenopus_tropicalis",
]


class AnimalTFDBFetcher:
    """AnimalTFDB 4.0 数据下载器

    下载动物转录因子注释和直系同源映射数据。

    Usage::

        fetcher = AnimalTFDBFetcher(basic_dir='./database/basic')
        # 下载物种注册信息（TF数量统计）
        fetcher.download_species_registry()
        # 下载指定物种的TF列表和同源映射
        fetcher.download_species_data("Bos_taurus")
    """

    REQUEST_TIMEOUT = 120

    def __init__(self, basic_dir: str):
        """
        Args:
            basic_dir: 基础缓存目录
        """
        self.basic_dir = Path(basic_dir)
        self.basic_dir.mkdir(parents=True, exist_ok=True)

    def _get_cache_dir(self) -> Path:
        """获取缓存目录"""
        cache_dir = self.basic_dir / "animaltfdb" / "AnimalTFDBv4.0"
        cache_dir.mkdir(parents=True, exist_ok=True)
        return cache_dir

    def _download_url(self, url: str, local_path: Path, overwrite: bool = False) -> Path:
        """通用下载方法"""
        if local_path.exists() and not overwrite:
            logger.info(f"已缓存，跳过: {local_path.name}")
            return local_path

        logger.info(f"下载: {url}")

        try:
            resp = requests.get(
                url, timeout=self.REQUEST_TIMEOUT,
                stream=True, headers={"User-Agent": "Mozilla/5.0"},
            )
            resp.raise_for_status()
        except requests.RequestException as e:
            raise RuntimeError(f"下载失败: {e}") from e

        with open(local_path, "wb") as f:
            for chunk in resp.iter_content(chunk_size=8192):
                if chunk:
                    f.write(chunk)

        logger.info(f"已保存: {local_path} ({local_path.stat().st_size / 1024:.1f} KB)")
        return local_path

    def download_tf_list(self, species: str, overwrite: bool = False) -> Path:
        """下载指定物种的 TF 列表

        Args:
            species: 物种拉丁名（下划线格式，如 Bos_taurus）
            overwrite: 是否覆盖

        Returns:
            TF 列表文件路径
        """
        cache_dir = self._get_cache_dir()
        url = f"{ANIMALTFDB_BASE_URL}/TF_list_final/{species}_TF"
        local_path = cache_dir / f"{species}_TF"
        return self._download_url(url, local_path, overwrite)

    def download_ortholog_to_human(self, species: str, overwrite: bool = False) -> Path:
        """下载指定物种到人类的直系同源映射

        Args:
            species: 物种拉丁名（下划线格式）
            overwrite: 是否覆盖

        Returns:
            同源映射文件路径
        """
        cache_dir = self._get_cache_dir()
        url = f"{ANIMALTFDB_BASE_URL}/ortholog_to_human_download/{species}_ortholog_to_human"
        local_path = cache_dir / f"{species}_ortholog_to_human"
        return self._download_url(url, local_path, overwrite)

    def download_species_data(self, species: str, overwrite: bool = False) -> Dict[str, Path]:
        """下载指定物种的全部数据（TF列表 + 同源映射）

        Args:
            species: 物种拉丁名（下划线格式）
            overwrite: 是否覆盖

        Returns:
            {'tf_list': Path, 'ortholog': Path}
        """
        results = {}
        try:
            results['tf_list'] = self.download_tf_list(species, overwrite)
        except Exception as e:
            logger.error(f"下载 {species} TF列表失败: {e}")

        try:
            results['ortholog'] = self.download_ortholog_to_human(species, overwrite)
        except Exception as e:
            logger.error(f"下载 {species} 同源映射失败: {e}")

        return results

    def download_htftarget(self, overwrite: bool = False) -> Path:
        """下载 hTFtarget 人类TF-target关系数据

        Args:
            overwrite: 是否覆盖

        Returns:
            hTFtarget 文件路径
        """
        from .htftarget_fetcher import HTFtargetFetcher
        fetcher = HTFtargetFetcher(str(self.basic_dir))
        return fetcher.download(overwrite)

    @staticmethod
    def get_priority_species() -> List[str]:
        """获取优先支持的物种列表"""
        return ANIMALTFDB_PRIORITY_SPECIES.copy()

    @staticmethod
    def get_info() -> Dict[str, str]:
        """获取数据库基本信息"""
        return {
            "name": "AnimalTFDB",
            "version": "4.0",
            "url": "https://guolab.wchscu.cn/AnimalTFDB4/",
            "species_count": 183,
            "tf_count": 274633,
            "description": "183个动物物种的TF分类注释和直系同源映射",
        }
```

- [ ] **Step 2: 验证下载器**

Run: `python -c "from allenricher.database.animaltfdb_fetcher import AnimalTFDBFetcher; print(f'Priority species: {len(AnimalTFDBFetcher.get_priority_species())}'); print('AnimalTFDBFetcher OK')"`
Expected: `Priority species: 19 AnimalTFDBFetcher OK`

---

## Task 3: hTFtarget 解析器

**Files:**
- Create: `allenricher/database/parsers/htftarget.py`

- [ ] **Step 1: 创建 hTFtarget 解析器类**

```python
"""
hTFtarget 数据库解析器

解析 hTFtarget 的 TF-target 关系 TSV 文件。

输入文件：
- tf-target-information.txt: TF, target, tissue

输出格式：
- {species}.hTF_2tf.tab.gz: TF 信息表（TF, target_count, tissues）
- {species}.hTF_2gene.tab.gz: Gene x TF 0/1 矩阵
- {species}.hTF_2disc.gz: TF 描述信息
"""

import gzip
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple
from collections import defaultdict
import logging

import pandas as pd

logger = logging.getLogger(__name__)


class HTFtargetParser:
    """hTFtarget 数据库解析器

    解析 hTFtarget 的 TF-target 关系文件。
    """

    @staticmethod
    def parse_tsv(tsv_path: str) -> Tuple[Dict[str, Set[str]], Dict[str, Set[str]], Dict[str, Set[str]]]:
        """解析 hTFtarget TSV 文件

        Args:
            tsv_path: tf-target-information.txt 文件路径

        Returns:
            (tf_to_targets, gene_to_tfs, tf_to_tissues)
            - tf_to_targets: {TF: {target1, target2, ...}}
            - gene_to_tfs: {gene: {TF1, TF2, ...}}
            - tf_to_tissues: {TF: {tissue1, tissue2, ...}}
        """
        tf_to_targets: Dict[str, Set[str]] = defaultdict(set)
        gene_to_tfs: Dict[str, Set[str]] = defaultdict(set)
        tf_to_tissues: Dict[str, Set[str]] = defaultdict(set)

        count = 0
        with open(tsv_path, 'r') as f:
            header = f.readline()  # 跳过表头: TF\ttarget\ttissue
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) < 3:
                    continue

                tf = parts[0].strip()
                target = parts[1].strip()
                tissues = parts[2].strip()

                if not tf or not target:
                    continue

                tf_to_targets[tf].add(target)
                gene_to_tfs[target].add(tf)

                # tissue 可能是逗号分隔的多个值
                for tissue in tissues.split(','):
                    tissue = tissue.strip()
                    if tissue:
                        tf_to_tissues[tf].add(tissue)

                count += 1

        logger.info(f"hTFtarget: 解析 {count} 条 TF-target 关系, "
                     f"{len(tf_to_targets)} 个 TF, {len(gene_to_tfs)} 个靶基因")
        return dict(tf_to_targets), dict(gene_to_tfs), dict(tf_to_tissues)

    @staticmethod
    def build_database(tsv_path: str, output_dir: str, species: str,
                       valid_genes: Optional[Set[str]] = None) -> None:
        """构建 hTFtarget 数据库

        Args:
            tsv_path: tf-target-information.txt 路径
            output_dir: 输出目录
            species: 物种代码（如 hsa）
            valid_genes: 有效基因集合（可选）
        """
        outdir = Path(output_dir)
        outdir.mkdir(parents=True, exist_ok=True)

        logger.info(f"HTFtargetParser: 开始构建数据库 (species={species})")

        tf_to_targets, gene_to_tfs, tf_to_tissues = HTFtargetParser.parse_tsv(tsv_path)

        # 过滤有效基因
        if valid_genes:
            tf_to_targets = {
                tf: {g for g in targets if g in valid_genes}
                for tf, targets in tf_to_targets.items()
            }
            tf_to_targets = {tf: targets for tf, targets in tf_to_targets.items() if targets}
            logger.info(f"过滤后: {len(tf_to_targets)} 个 TF")

        # 获取所有基因
        all_genes = set()
        for targets in tf_to_targets.values():
            all_genes.update(targets)
        all_tfs = set(tf_to_targets.keys())

        # 构建 Gene x TF 矩阵
        gene_list = sorted(all_genes)
        tf_list = sorted(all_tfs)

        gene2tf_file = outdir / f"{species}.hTF_2gene.tab.gz"
        logger.info(f"写入文件: {gene2tf_file}")

        with gzip.open(gene2tf_file, 'wt') as f:
            f.write('Gene\t' + '\t'.join(tf_list) + '\n')
            for gene in gene_list:
                regulating_tfs = gene_to_tfs.get(gene, set())
                values = ['1' if tf in regulating_tfs else '0' for tf in tf_list]
                f.write(gene + '\t' + '\t'.join(values) + '\n')

        # 构建 TF 描述文件
        disc_file = outdir / f"{species}.hTF_2disc.gz"
        logger.info(f"写入文件: {disc_file}")

        with gzip.open(disc_file, 'wt') as f:
            f.write("TF\ttarget_count\ttissues\tsource\n")
            for tf in tf_list:
                target_count = len(tf_to_targets[tf])
                tissues = ','.join(sorted(tf_to_tissues.get(tf, set())))
                f.write(f"{tf}\t{target_count}\t{tissues}\thTFtarget\n")

        logger.info(f"HTFtargetParser: 数据库构建完成")
```

- [ ] **Step 2: 导出解析器**

在 `allenricher/database/parsers/__init__.py` 中添加:
```python
from .htftarget import HTFtargetParser

__all__ = [..., 'HTFtargetParser']
```

- [ ] **Step 3: 验证解析器**

Run: `python -c "from allenricher.database.parsers import HTFtargetParser; print('HTFtargetParser OK')"`
Expected: `HTFtargetParser OK`

---

## Task 4: AnimalTFDB 解析器

**Files:**
- Create: `allenricher/database/parsers/animaltfdb.py`

- [ ] **Step 1: 创建 AnimalTFDB 解析器类**

```python
"""
AnimalTFDB 4.0 数据库解析器

解析 AnimalTFDB 的 TF 列表和直系同源映射文件。

输入文件：
- {Species}_TF: TF 列表 (Species, Symbol, Ensembl, Family, Protein, Entrez_ID)
- {Species}_ortholog_to_human: 直系同源映射

输出格式：
- {species}.AnimalTFDB_2tf.tab.gz: TF 信息表
- {species}.AnimalTFDB_2disc.gz: TF 描述信息
- {species}.AnimalTFDB_ortholog.gz: 同源映射表
"""

import gzip
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple
from collections import defaultdict
import logging

import pandas as pd

logger = logging.getLogger(__name__)


class AnimalTFDBParser:
    """AnimalTFDB 4.0 数据库解析器

    解析 TF 列表和直系同源映射。
    """

    @staticmethod
    def parse_tf_list(tf_list_path: str) -> pd.DataFrame:
        """解析 TF 列表文件

        Args:
            tf_list_path: {Species}_TF 文件路径

        Returns:
            TF 信息 DataFrame，列: Species, Symbol, Ensembl, Family, Protein, Entrez_ID
        """
        df = pd.read_csv(tf_list_path, sep='\t', low_memory=False)

        # 标准化列名（去除可能的空格）
        df.columns = df.columns.str.strip()

        logger.info(f"AnimalTFDB TF列表: {len(df)} 个 TF")
        return df

    @staticmethod
    def parse_ortholog_to_human(ortholog_path: str) -> Dict[str, str]:
        """解析直系同源映射文件

        Args:
            ortholog_path: {Species}_ortholog_to_human 文件路径

        Returns:
            {物种基因Symbol: 人类基因Symbol} 映射字典
        """
        ortholog_map: Dict[str, str] = {}

        with open(ortholog_path, 'r') as f:
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) >= 2:
                    species_gene = parts[0].strip()
                    human_gene = parts[1].strip()
                    if species_gene and human_gene:
                        ortholog_map[species_gene] = human_gene

        logger.info(f"AnimalTFDB 同源映射: {len(ortholog_map)} 对")
        return ortholog_map

    @staticmethod
    def build_database(
        tf_list_path: str,
        ortholog_path: str,
        output_dir: str,
        species: str,
        valid_genes: Optional[Set[str]] = None
    ) -> Tuple[pd.DataFrame, Dict[str, str]]:
        """构建 AnimalTFDB 数据库

        Args:
            tf_list_path: TF 列表文件路径
            ortholog_path: 同源映射文件路径
            output_dir: 输出目录
            species: 物种代码（如 bta 代表牛）
            valid_genes: 有效基因集合（可选）

        Returns:
            (tf_df, ortholog_map) 用于后续同源映射
        """
        outdir = Path(output_dir)
        outdir.mkdir(parents=True, exist_ok=True)

        logger.info(f"AnimalTFDBParser: 开始构建数据库 (species={species})")

        # 解析 TF 列表
        tf_df = AnimalTFDBParser.parse_tf_list(tf_list_path)

        # 过滤有效基因
        if valid_genes:
            tf_df = tf_df[tf_df['Symbol'].isin(valid_genes)]
            logger.info(f"过滤后: {len(tf_df)} 个 TF")

        # 解析同源映射
        ortholog_map = AnimalTFDBParser.parse_ortholog_to_human(ortholog_path)

        # 保存 TF 信息表
        tf_file = outdir / f"{species}.AnimalTFDB_2tf.tab.gz"
        logger.info(f"写入文件: {tf_file}")

        with gzip.open(tf_file, 'wt') as f:
            f.write('\t'.join(tf_df.columns) + '\n')
            for _, row in tf_df.iterrows():
                f.write('\t'.join(str(v) for v in row.values) + '\n')

        # 保存 TF 描述文件
        disc_file = outdir / f"{species}.AnimalTFDB_2disc.gz"
        logger.info(f"写入文件: {disc_file}")

        with gzip.open(disc_file, 'wt') as f:
            f.write("TF\tFamily\tEntrez_ID\tEnsembl\tsource\n")
            for _, row in tf_df.iterrows():
                symbol = row.get('Symbol', '')
                family = row.get('Family', 'Unknown')
                entrez = row.get('Entrez_ID', 'NA')
                ensembl = row.get('Ensembl', 'NA')
                f.write(f"{symbol}\t{family}\t{entrez}\t{ensembl}\tAnimalTFDB\n")

        # 保存同源映射
        ortholog_file = outdir / f"{species}.AnimalTFDB_ortholog.gz"
        logger.info(f"写入文件: {ortholog_file}")

        with gzip.open(ortholog_file, 'wt') as f:
            f.write("Species_Gene\tHuman_Gene\n")
            for sp_gene, hu_gene in ortholog_map.items():
                f.write(f"{sp_gene}\t{hu_gene}\n")

        logger.info(f"AnimalTFDBParser: 数据库构建完成")
        return tf_df, ortholog_map
```

- [ ] **Step 2: 导出解析器**

在 `allenricher/database/parsers/__init__.py` 中添加:
```python
from .animaltfdb import AnimalTFDBParser

__all__ = [..., 'AnimalTFDBParser']
```

- [ ] **Step 3: 验证解析器**

Run: `python -c "from allenricher.database.parsers import AnimalTFDBParser; print('AnimalTFDBParser OK')"`
Expected: `AnimalTFDBParser OK`

---

## Task 5: 同源映射引擎（核心）

**Files:**
- Create: `allenricher/database/ortholog_mapper.py`

- [ ] **Step 1: 创建同源映射引擎**

```python
"""
同源映射引擎

通过 AnimalTFDB 的直系同源映射，将人类 TF-target 关系推断到目标物种。

映射流程：
1. 人类 TF → 靶基因 (hTFtarget)
2. 目标物种基因 → 人类基因 (AnimalTFDB ortholog_to_human, 反向)
3. 目标物种 TF → 人类 TF (AnimalTFDB ortholog_to_human)
4. 合并得到：目标物种 TF → 目标物种靶基因
"""

import gzip
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple
from collections import defaultdict
import logging

import pandas as pd

logger = logging.getLogger(__name__)


class OrthologMapper:
    """同源映射引擎

    将人类 TF-target 关系通过直系同源映射推断到目标物种。
    """

    def __init__(
        self,
        human_tf_to_targets: Dict[str, Set[str]],
        species_to_human: Dict[str, str],
        species_tf_set: Optional[Set[str]] = None,
    ):
        """初始化同源映射引擎

        Args:
            human_tf_to_targets: 人类 TF→靶基因关系 {TF: {target1, target2, ...}}
            species_to_human: 目标物种基因→人类基因映射 {sp_gene: human_gene}
            species_tf_set: 目标物种的 TF 集合（可选，用于过滤）
        """
        self.human_tf_to_targets = human_tf_to_targets
        self.species_to_human = species_to_human
        self.species_tf_set = species_tf_set

        # 构建反向映射：人类基因 → 目标物种基因
        self.human_to_species: Dict[str, Set[str]] = defaultdict(set)
        for sp_gene, hu_gene in species_to_human.items():
            self.human_to_species[hu_gene].add(sp_gene)

        # 构建目标物种 TF → 人类 TF 映射
        self.species_tf_to_human_tf: Dict[str, str] = {}
        if species_tf_set:
            for sp_tf in species_tf_set:
                if sp_tf in species_to_human:
                    self.species_tf_to_human_tf[sp_tf] = species_to_human[sp_tf]

    def map_tf_targets(self) -> Tuple[Dict[str, Set[str]], Dict[str, Set[str]]]:
        """执行同源映射，推断目标物种的 TF-target 关系

        Returns:
            (tf_to_targets, gene_to_tfs)
            - tf_to_targets: {目标物种TF: {目标物种靶基因1, ...}}
            - gene_to_tfs: {目标物种基因: {目标物种TF1, ...}}
        """
        tf_to_targets: Dict[str, Set[str]] = defaultdict(set)
        gene_to_tfs: Dict[str, Set[str]] = defaultdict(set)

        mapped_tf_count = 0
        unmapped_tf_count = 0

        for sp_tf, human_tf in self.species_tf_to_human_tf.items():
            # 获取人类 TF 的靶基因
            human_targets = self.human_tf_to_targets.get(human_tf, set())
            if not human_targets:
                unmapped_tf_count += 1
                continue

            # 将人类靶基因映射回目标物种
            for human_target in human_targets:
                sp_targets = self.human_to_species.get(human_target, set())
                for sp_target in sp_targets:
                    # 排除自映射（TF不能调控自己）
                    if sp_target != sp_tf:
                        tf_to_targets[sp_tf].add(sp_target)
                        gene_to_tfs[sp_target].add(sp_tf)

            mapped_tf_count += 1

        logger.info(
            f"同源映射完成: {mapped_tf_count} 个TF成功映射, "
            f"{unmapped_tf_count} 个TF无人类对应关系, "
            f"共 {len(tf_to_targets)} 个TF有靶基因, "
            f"共 {len(gene_to_tfs)} 个靶基因"
        )

        return dict(tf_to_targets), dict(gene_to_tfs)

    @staticmethod
    def build_mapped_database(
        tf_to_targets: Dict[str, Set[str]],
        gene_to_tfs: Dict[str, Set[str]],
        species_tf_df: pd.DataFrame,
        output_dir: str,
        species: str,
    ) -> None:
        """构建同源映射后的数据库文件

        Args:
            tf_to_targets: 映射后的 TF→靶基因关系
            gene_to_tfs: 映射后的 基因→TF 关系
            species_tf_df: 目标物种 TF 信息 DataFrame
            output_dir: 输出目录
            species: 物种代码
        """
        outdir = Path(output_dir)
        outdir.mkdir(parents=True, exist_ok=True)

        all_genes = set(gene_to_tfs.keys())
        all_tfs = set(tf_to_targets.keys())
        gene_list = sorted(all_genes)
        tf_list = sorted(all_tfs)

        # 构建 Gene x TF 矩阵
        gene2tf_file = outdir / f"{species}.AnimalTFDB_2gene.tab.gz"
        logger.info(f"写入文件: {gene2tf_file}")

        with gzip.open(gene2tf_file, 'wt') as f:
            f.write('Gene\t' + '\t'.join(tf_list) + '\n')
            for gene in gene_list:
                regulating_tfs = gene_to_tfs.get(gene, set())
                values = ['1' if tf in regulating_tfs else '0' for tf in tf_list]
                f.write(gene + '\t' + '\t'.join(values) + '\n')

        # 构建 TF 描述文件
        disc_file = outdir / f"{species}.AnimalTFDB_mapped_2disc.gz"
        logger.info(f"写入文件: {disc_file}")

        # 构建 TF family 查找表
        tf_family_map = {}
        if species_tf_df is not None and 'Symbol' in species_tf_df.columns:
            for _, row in species_tf_df.iterrows():
                tf_family_map[row['Symbol']] = row.get('Family', 'Unknown')

        with gzip.open(disc_file, 'wt') as f:
            f.write("TF\ttarget_count\tFamily\tsource\n")
            for tf in tf_list:
                target_count = len(tf_to_targets[tf])
                family = tf_family_map.get(tf, 'Unknown')
                f.write(f"{tf}\t{target_count}\t{family}\tAnimalTFDB_ortholog\n")

        logger.info(f"OrthologMapper: 数据库构建完成")
```

- [ ] **Step 2: 验证映射引擎**

Run: `python -c "from allenricher.database.ortholog_mapper import OrthologMapper; print('OrthologMapper OK')"`
Expected: `OrthologMapper OK`

---

## Task 6: 更新 DatabaseBuilder

**Files:**
- Modify: `allenricher/database/builder.py`

- [ ] **Step 1: 添加 build_htftarget 方法**

在 `DatabaseBuilder` 类中添加:

```python
def build_htftarget(self, species: str, taxid: int) -> str:
    """构建 hTFtarget 数据库（人类专用）

    Args:
        species: 物种代码（如 hsa）
        taxid: NCBI TaxID

    Returns:
        输出目录路径
    """
    from .htftarget_fetcher import HTFtargetFetcher
    from .parsers.htftarget import HTFtargetParser

    # 查找已下载的数据
    htftarget_file = self.basic_dir / "htftarget" / "tf-target-information.txt"

    if not htftarget_file.exists():
        print(f"|--- [跳过] 未找到 hTFtarget 文件: {htftarget_file}")
        print("|--- 请先运行: allenricher download --animaltfdb")
        return ""

    # 输出目录
    date_str = datetime.now().strftime("%Y%m%d")
    outdir = self.organism_dir / f"v{date_str}" / species
    outdir.mkdir(parents=True, exist_ok=True)

    print(f"\n{'='*60}")
    print(f"构建 hTFtarget 数据库: {species} (taxid={taxid})")
    print(f"数据源: {htftarget_file}")
    print(f"输出目录: {outdir}")
    print(f"{'='*60}")

    # 获取有效基因集合
    gene_info_path = self._get_gene_info_path()
    valid_genes = None
    if gene_info_path:
        valid_genes = self._load_valid_genes(str(gene_info_path), taxid)

    # 构建数据库
    HTFtargetParser.build_database(
        tsv_path=str(htftarget_file),
        output_dir=str(outdir),
        species=species,
        valid_genes=valid_genes,
    )

    # 验证输出
    expected_files = [
        f"{species}.hTF_2gene.tab.gz",
        f"{species}.hTF_2disc.gz",
    ]
    for fname in expected_files:
        fpath = outdir / fname
        if fpath.exists():
            print(f"    ✅ {fname}")
        else:
            print(f"    ❌ {fname} - 未生成")

    print(f"\nhTFtarget 数据库构建完成 → {outdir}")
    return str(outdir)
```

- [ ] **Step 2: 添加 build_animaltfdb 方法**

```python
def build_animaltfdb(self, species: str, taxid: int,
                     species_latin: str = "") -> str:
    """构建 AnimalTFDB 数据库（通过同源映射）

    对于人类(hsa)：直接构建 hTFtarget 数据库
    对于其他物种：下载TF列表+同源映射，通过同源映射推断TF-target关系

    Args:
        species: 物种代码（如 bta 代表牛）
        taxid: NCBI TaxID
        species_latin: 物种拉丁名（下划线格式，如 Bos_taurus）

    Returns:
        输出目录路径
    """
    from .animaltfdb_fetcher import AnimalTFDBFetcher
    from .parsers.animaltfdb import AnimalTFDBParser
    from .parsers.htftarget import HTFtargetParser
    from .ortholog_mapper import OrthologMapper

    # 人类直接使用 hTFtarget
    if species == 'hsa':
        return self.build_htftarget(species, taxid)

    # 确定物种拉丁名
    if not species_latin:
        # 尝试从物种注册表查找
        from .species_registry import SpeciesRegistry
        registry = SpeciesRegistry.load_default(str(self.root_dir))
        entry = registry.query_by_kegg_code(species)
        if entry and entry.latin_name:
            species_latin = entry.latin_name.replace(' ', '_')
        else:
            print(f"|--- [错误] 无法确定物种 {species} 的拉丁名，请使用 --latin-name 参数")
            return ""

    cache_dir = self.basic_dir / "animaltfdb" / "AnimalTFDBv4.0"

    # 查找已下载的数据
    tf_list_file = cache_dir / f"{species_latin}_TF"
    ortholog_file = cache_dir / f"{species_latin}_ortholog_to_human"
    htftarget_file = self.basic_dir / "htftarget" / "tf-target-information.txt"

    missing = []
    if not tf_list_file.exists():
        missing.append(f"TF列表: {tf_list_file}")
    if not ortholog_file.exists():
        missing.append(f"同源映射: {ortholog_file}")
    if not htftarget_file.exists():
        missing.append(f"hTFtarget: {htftarget_file}")

    if missing:
        print(f"|--- [跳过] 缺少必要文件:")
        for m in missing:
            print(f"|---   - {m}")
        print(f"|--- 请先运行: allenricher download --animaltfdb --species {species_latin}")
        return ""

    # 输出目录
    date_str = datetime.now().strftime("%Y%m%d")
    outdir = self.organism_dir / f"v{date_str}" / species
    outdir.mkdir(parents=True, exist_ok=True)

    print(f"\n{'='*60}")
    print(f"构建 AnimalTFDB 数据库: {species} ({species_latin}, taxid={taxid})")
    print(f"策略: 同源映射 (hTFtarget → {species_latin})")
    print(f"输出目录: {outdir}")
    print(f"{'='*60}")

    # 获取有效基因集合
    gene_info_path = self._get_gene_info_path()
    valid_genes = None
    if gene_info_path:
        valid_genes = self._load_valid_genes(str(gene_info_path), taxid)

    # Step 1: 解析 AnimalTFDB TF列表和同源映射
    print("\n[1/3] 解析 AnimalTFDB 数据...")
    tf_df, ortholog_map = AnimalTFDBParser.build_database(
        tf_list_path=str(tf_list_file),
        ortholog_path=str(ortholog_file),
        output_dir=str(outdir),
        species=species,
        valid_genes=valid_genes,
    )

    # Step 2: 解析 hTFtarget
    print("[2/3] 解析 hTFtarget 人类TF-target关系...")
    human_tf_to_targets, _, _ = HTFtargetParser.parse_tsv(str(htftarget_file))

    # Step 3: 同源映射
    print("[3/3] 执行同源映射...")
    species_tf_set = set(tf_df['Symbol'].values) if 'Symbol' in tf_df.columns else None

    mapper = OrthologMapper(
        human_tf_to_targets=human_tf_to_targets,
        species_to_human=ortholog_map,
        species_tf_set=species_tf_set,
    )

    tf_to_targets, gene_to_tfs = mapper.map_tf_targets()

    if not tf_to_targets:
        print("|--- [警告] 同源映射结果为空，无法构建TF-target关系")
        return str(outdir)

    # 构建映射后的数据库文件
    OrthologMapper.build_mapped_database(
        tf_to_targets=tf_to_targets,
        gene_to_tfs=gene_to_tfs,
        species_tf_df=tf_df,
        output_dir=str(outdir),
        species=species,
    )

    # 验证输出
    expected_files = [
        f"{species}.AnimalTFDB_2tf.tab.gz",
        f"{species}.AnimalTFDB_2disc.gz",
        f"{species}.AnimalTFDB_ortholog.gz",
        f"{species}.AnimalTFDB_2gene.tab.gz",
        f"{species}.AnimalTFDB_mapped_2disc.gz",
    ]
    for fname in expected_files:
        fpath = outdir / fname
        if fpath.exists():
            print(f"    ✅ {fname}")
        else:
            print(f"    ❌ {fname} - 未生成")

    print(f"\nAnimalTFDB 数据库构建完成 → {outdir}")
    return str(outdir)
```

- [ ] **Step 3: 更新 build_species_db 方法**

在 `build_species_db` 的数据库列表处理中添加:

```python
elif db_upper == 'ANIMALTFDB':
    latin_name = kwargs.get('latin_name', '')
    self.build_animaltfdb(species, taxid, species_latin=latin_name)
elif db_upper == 'HTFTARGET':
    self.build_htftarget(species, taxid)
```

---

## Task 7: 更新 DatabaseManager

**Files:**
- Modify: `allenricher/database/manager.py`

- [ ] **Step 1: 添加 load_animaltfdb 和 load_htftarget 方法**

```python
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

    # 优先加载同源映射结果，回退到 hTFtarget（人类）
    gene2tf_file = db_dir / f"{sp}.AnimalTFDB_2gene.tab.gz"
    disc_file = db_dir / f"{sp}.AnimalTFDB_mapped_2disc.gz"

    if not gene2tf_file.exists():
        # 回退到 hTFtarget
        return self.load_htftarget(sp)

    result = {}
    result['gene2tf'] = pd.read_csv(gene2tf_file, sep='\t', compression='gzip', low_memory=False)

    if disc_file.exists():
        result['tf_info'] = pd.read_csv(disc_file, sep='\t', compression='gzip')

    return result
```

- [ ] **Step 2: 更新 load_database 方法**

在 `name_to_prefix` 映射中添加:

```python
'ANIMALTFDB': 'AnimalTFDB_2gene',
'HTFTARGET': 'hTF_2gene',
```

---

## Task 8: 更新物种注册表

**Files:**
- Modify: `allenricher/database/species_registry.py`

- [ ] **Step 1: 添加 AnimalTFDB 字段**

在 `_FIELD_NAMES` 列表中添加:
```python
"has_animaltfdb", "animaltfdb_tf_count", "animaltfdb_mapped_target_count",
```

在 `SpeciesEntry` dataclass 中添加:
```python
# AnimalTFDB 相关字段
has_animaltfdb: bool = False
animaltfdb_tf_count: Optional[int] = None
animaltfdb_mapped_target_count: Optional[int] = None
```

在 `filter_by_databases` 方法中添加 `animaltfdb` 参数。

---

## Task 9: 更新 CLI

**Files:**
- Modify: `allenricher/cli.py`

- [ ] **Step 1: 添加 AnimalTFDB 下载命令**

```python
def _cmd_download_animaltfdb(args) -> int:
    """下载 AnimalTFDB 和 hTFtarget 数据"""
    from allenricher.database.animaltfdb_fetcher import AnimalTFDBFetcher

    fetcher = AnimalTFDBFetcher(basic_dir=args.database_dir + "/basic")

    # 始终下载 hTFtarget（映射源）
    print("下载 hTFtarget 人类TF-target关系...")
    fetcher.download_htftarget(overwrite=args.force)

    # 下载指定物种的TF列表和同源映射
    species_list = args.species.split(',') if args.species else []

    if not species_list:
        print("未指定物种，仅下载 hTFtarget 映射源。")
        print("使用 --species Bos_taurus,Sus_scrofa 指定要下载的物种。")
        return 0

    print(f"下载 {len(species_list)} 个物种的 AnimalTFDB 数据...")
    for sp in species_list:
        print(f"\n--- 下载 {sp} ---")
        fetcher.download_species_data(sp, overwrite=args.force)

    print("\n下载完成")
    return 0
```

- [ ] **Step 2: 添加 AnimalTFDB 构建支持**

在 `_cmd_build` 处理函数中添加 `--latin-name` 参数支持，并添加 `ANIMALTFDB` 和 `HTFTARGET` 数据库选项。

- [ ] **Step 3: 更新 tf-enrich 命令**

在 `tf-enrich` 命令的 `--database` choices 中添加 `'animaltfdb'` 和 `'htftarget'`。

- [ ] **Step 4: 注册命令**

在 CLI 主函数中添加:

```python
# download 子命令
download_parser.add_argument('--animaltfdb', action='store_true',
                            help='下载 AnimalTFDB 和 hTFtarget 数据')
download_parser.add_argument('--species', type=str, default='',
                            help='要下载的物种列表（逗号分隔，如 Bos_taurus,Sus_scrofa）')

# build 子命令
build_parser.add_argument('--latin-name', type=str, default='',
                          help='物种拉丁名（下划线格式，如 Bos_taurus）')

# tf-enrich 子命令
tf_enrich_parser.add_argument('-d', '--database', default='trrust',
                              choices=['trrust', 'chea3', 'animaltfdb', 'htftarget'],
                              help='TF数据库')
```

---

## Task 10: 更新 HTML 报告模板

**Files:**
- Modify: `allenricher/report/templates/tf_report.html`

- [ ] **Step 1: 添加映射来源标注**

在报告模板中添加数据来源说明区域:

```html
{% if source == 'AnimalTFDB_ortholog' %}
<div class="alert alert-info">
    <strong>数据来源说明：</strong>
    本分析基于 AnimalTFDB 同源映射策略。TF-target 关系通过人类 ChIP-Seq 实验数据
    (hTFtarget) 经直系同源映射推断得到。映射准确性取决于物种间的直系同源关系保守程度。
    人类 TF-target 关系来源: hTFtarget (基于 ENCODE/SRA ChIP-Seq, 30% 阈值过滤)。
</div>
{% elif source == 'hTFtarget' %}
<div class="alert alert-info">
    <strong>数据来源：</strong>
    hTFtarget - 人类转录因子靶基因数据库，基于 ENCODE/SRA ChIP-Seq 数据。
</div>
{% endif %}
```

- [ ] **Step 2: 更新报告生成逻辑**

在 `_generate_tf_enrichment_report()` 函数中，根据数据库类型传递 `source` 变量:

```python
source_map = {
    'trrust': 'TRRUST',
    'chea3': 'ChEA3',
    'animaltfdb': 'AnimalTFDB_ortholog',
    'htftarget': 'hTFtarget',
}
```

---

## Task 11: 端到端测试

**Files:**
- Create: `test_e2e_2026/test_animaltfdb.py`

- [ ] **Step 1: 创建测试文件**

```python
"""
AnimalTFDB + hTFtarget 端到端测试
"""

import pytest
import gzip
import tempfile
from pathlib import Path
from unittest.mock import patch, MagicMock
import numpy as np
import pandas as pd


class TestHTFtargetFetcher:
    """hTFtarget 下载器测试"""

    def test_import(self):
        from allenricher.database.htftarget_fetcher import HTFtargetFetcher
        assert HTFtargetFetcher is not None

    def test_get_info(self):
        from allenricher.database.htftarget_fetcher import HTFtargetFetcher
        info = HTFtargetFetcher.get_info()
        assert info['name'] == 'hTFtarget'
        assert info['species'] == 'Homo sapiens'


class TestAnimalTFDBFetcher:
    """AnimalTFDB 下载器测试"""

    def test_import(self):
        from allenricher.database.animaltfdb_fetcher import AnimalTFDBFetcher
        assert AnimalTFDBFetcher is not None

    def test_priority_species(self):
        from allenricher.database.animaltfdb_fetcher import AnimalTFDBFetcher
        species = AnimalTFDBFetcher.get_priority_species()
        assert len(species) >= 19
        assert "Bos_taurus" in species
        assert "Sus_scrofa" in species
        assert "Homo_sapiens" in species

    def test_get_info(self):
        from allenricher.database.animaltfdb_fetcher import AnimalTFDBFetcher
        info = AnimalTFDBFetcher.get_info()
        assert info['name'] == 'AnimalTFDB'
        assert info['species_count'] == 183


class TestHTFtargetParser:
    """hTFtarget 解析器测试"""

    def test_import(self):
        from allenricher.database.parsers import HTFtargetParser
        assert HTFtargetParser is not None

    def test_parse_tsv(self, tmp_path):
        """测试解析 hTFtarget TSV 文件"""
        from allenricher.database.parsers.htftarget import HTFtargetParser

        # 创建模拟 TSV 文件
        tsv_file = tmp_path / "test_htftarget.txt"
        tsv_file.write_text(
            "TF\ttarget\ttissue\n"
            "AEBP2\tTMEM53\tcolon\n"
            "AEBP2\tC1orf228\tcolon\n"
            "AEBP2\tFBXO31\tcolon,blood\n"
            "AFF4\tMIR31HG\tcolon\n"
            "AFF4\tMIR29A\tunclear\n"
        )

        tf_to_targets, gene_to_tfs, tf_to_tissues = HTFtargetParser.parse_tsv(str(tsv_file))

        assert len(tf_to_targets) == 2
        assert len(tf_to_targets['AEBP2']) == 3
        assert 'TMEM53' in tf_to_targets['AEBP2']
        assert 'MIR31HG' in tf_to_targets['AFF4']

        assert len(gene_to_tfs) == 5
        assert 'AEBP2' in gene_to_tfs['TMEM53']

        assert 'colon' in tf_to_tissues['AEBP2']
        assert 'blood' in tf_to_tissues['AEBP2']


class TestAnimalTFDBParser:
    """AnimalTFDB 解析器测试"""

    def test_import(self):
        from allenricher.database.parsers import AnimalTFDBParser
        assert AnimalTFDBParser is not None

    def test_parse_tf_list(self, tmp_path):
        """测试解析 TF 列表"""
        from allenricher.database.parsers.animaltfdb import AnimalTFDBParser

        tf_file = tmp_path / "test_TF"
        tf_file.write_text(
            "Species\tSymbol\tEnsembl\tFamily\tProtein\tEntrez_ID\n"
            "Bos_taurus\tZFP69\tENSBTAG00000039287\tzf-C2H2\tENSBTAP00000061159.1\t521549\n"
            "Bos_taurus\tEN2\tENSBTAG00000012345\tHomeobox\tENSBTAP00000012345.1\t613654\n"
        )

        df = AnimalTFDBParser.parse_tf_list(str(tf_file))
        assert len(df) == 2
        assert 'Symbol' in df.columns
        assert 'Family' in df.columns

    def test_parse_ortholog(self, tmp_path):
        """测试解析同源映射"""
        from allenricher.database.parsers.animaltfdb import AnimalTFDBParser

        ortholog_file = tmp_path / "test_ortholog"
        ortholog_file.write_text(
            "BTA12345\tTP53\n"
            "BTA67890\tMYC\n"
            "BTA11111\tBRCA1\n"
        )

        ortholog_map = AnimalTFDBParser.parse_ortholog_to_human(str(ortholog_file))
        assert len(ortholog_map) == 3
        assert ortholog_map['BTA12345'] == 'TP53'
        assert ortholog_map['BTA67890'] == 'MYC'


class TestOrthologMapper:
    """同源映射引擎测试"""

    def test_import(self):
        from allenricher.database.ortholog_mapper import OrthologMapper
        assert OrthologMapper is not None

    def test_map_tf_targets(self):
        """测试同源映射"""
        from allenricher.database.ortholog_mapper import OrthologMapper

        # 人类 TF-target 关系
        human_tf_to_targets = {
            'TP53': {'CDKN1A', 'MDM2', 'BAX'},
            'MYC': {'CDK4', 'CCND1', 'BAX'},
        }

        # 目标物种基因→人类基因映射
        species_to_human = {
            'BTA12345': 'TP53',  # 牛TP53同源
            'BTA67890': 'MYC',   # 牛MYC同源
            'BTA11111': 'CDKN1A', # 牛CDKN1A同源
            'BTA22222': 'MDM2',   # 牛MDM2同源
            'BTA33333': 'BAX',   # 牛BAX同源
            'BTA44444': 'CDK4',   # 牛CDK4同源
            'BTA55555': 'CCND1',  # 牛CCND1同源
        }

        # 目标物种 TF 集合
        species_tf_set = {'BTA12345', 'BTA67890'}

        mapper = OrthologMapper(
            human_tf_to_targets=human_tf_to_targets,
            species_to_human=species_to_human,
            species_tf_set=species_tf_set,
        )

        tf_to_targets, gene_to_tfs = mapper.map_tf_targets()

        # BTA12345(牛TP53) 应映射到 CDKN1A, MDM2, BAX 的牛同源
        assert 'BTA12345' in tf_to_targets
        assert 'BTA11111' in tf_to_targets['BTA12345']  # CDKN1A
        assert 'BTA22222' in tf_to_targets['BTA12345']  # MDM2
        assert 'BTA33333' in tf_to_targets['BTA12345']  # BAX

        # BTA67890(牛MYC) 应映射到 CDK4, CCND1, BAX 的牛同源
        assert 'BTA67890' in tf_to_targets
        assert 'BTA44444' in tf_to_targets['BTA67890']  # CDK4
        assert 'BTA55555' in tf_to_targets['BTA67890']  # CCND1

        # BAX 被两个 TF 共同调控
        assert 'BTA12345' in gene_to_tfs.get('BTA33333', set())
        assert 'BTA67890' in gene_to_tfs.get('BTA33333', set())

    def test_no_self_regulation(self):
        """测试TF不能调控自己"""
        from allenricher.database.ortholog_mapper import OrthologMapper

        human_tf_to_targets = {'TP53': {'TP53_TARGET'}}
        species_to_human = {'GENE_A': 'TP53', 'GENE_A_TARGET': 'TP53_TARGET'}
        species_tf_set = {'GENE_A'}

        mapper = OrthologMapper(
            human_tf_to_targets=human_tf_to_targets,
            species_to_human=species_to_human,
            species_tf_set=species_tf_set,
        )

        tf_to_targets, _ = mapper.map_tf_targets()
        # GENE_A 不应出现在自己的靶基因列表中
        assert 'GENE_A' not in tf_to_targets.get('GENE_A', set())


class TestIntegration:
    """集成测试：完整流程"""

    def test_full_pipeline(self, tmp_path):
        """测试完整解析+映射流程"""
        from allenricher.database.parsers.htftarget import HTFtargetParser
        from allenricher.database.parsers.animaltfdb import AnimalTFDBParser
        from allenricher.database.ortholog_mapper import OrthologMapper

        # 1. 创建模拟 hTFtarget 文件
        htftarget_file = tmp_path / "tf-target-information.txt"
        htftarget_file.write_text(
            "TF\ttarget\ttissue\n"
            "TP53\tCDKN1A\tblood\n"
            "TP53\tMDM2\tcolon\n"
            "TP53\tBAX\tblood,colon\n"
            "MYC\tCDK4\tblood\n"
            "MYC\tCCND1\tcolon\n"
        )

        # 2. 创建模拟 AnimalTFDB 文件
        tf_list_file = tmp_path / "Test_species_TF"
        tf_list_file.write_text(
            "Species\tSymbol\tEnsembl\tFamily\tProtein\tEntrez_ID\n"
            "Test_species\tGENE_A\tENSTEST001\tzf-C2H2\tENSTESTP001\t100\n"
            "Test_species\tGENE_B\tENSTEST002\tHomeobox\tENSTESTP002\t200\n"
        )

        ortholog_file = tmp_path / "Test_species_ortholog_to_human"
        ortholog_file.write_text(
            "GENE_A\tTP53\n"
            "GENE_B\tMYC\n"
            "TARGET_1\tCDKN1A\n"
            "TARGET_2\tMDM2\n"
            "TARGET_3\tBAX\n"
            "TARGET_4\tCDK4\n"
            "TARGET_5\tCCND1\n"
        )

        # 3. 解析
        human_tf_to_targets, _, _ = HTFtargetParser.parse_tsv(str(htftarget_file))
        _, ortholog_map = AnimalTFDBParser.build_database(
            tf_list_path=str(tf_list_file),
            ortholog_path=str(ortholog_file),
            output_dir=str(tmp_path / "db"),
            species="tst",
        )

        # 4. 映射
        mapper = OrthologMapper(
            human_tf_to_targets=human_tf_to_targets,
            species_to_human=ortholog_map,
            species_tf_set={'GENE_A', 'GENE_B'},
        )

        tf_to_targets, gene_to_tfs = mapper.map_tf_targets()

        assert 'GENE_A' in tf_to_targets  # TP53 同源
        assert 'GENE_B' in tf_to_targets  # MYC 同源
        assert len(tf_to_targets['GENE_A']) == 3  # TP53 有 3 个靶基因
        assert len(tf_to_targets['GENE_B']) == 2  # MYC 有 2 个靶基因
```

- [ ] **Step 2: 运行测试**

Run: `cd f:\OneDrive\Documents\TraeSOLO\AllEnricher\AllEnricher-v2 && python -m pytest test_e2e_2026/test_animaltfdb.py -v`

---

## 自检清单

1. **Spec 覆盖**:
   - hTFtarget 数据下载 ✓
   - AnimalTFDB 数据下载 ✓
   - hTFtarget 解析 ✓
   - AnimalTFDB 解析 ✓
   - 同源映射引擎 ✓
   - 数据库构建 ✓
   - 数据库加载 ✓
   - 物种注册表更新 ✓
   - CLI 命令 ✓
   - HTML 报告更新 ✓
   - 端到端测试 ✓

2. **类型一致性**:
   - HTFtargetFetcher/AnimalTFDBFetcher 接口与 TRRUSTFetcher/ChEA3Fetcher 一致 ✓
   - HTFtargetParser/AnimalTFDBParser 接口与 TRRUSTParser/ChEA3Parser 一致 ✓
   - DatabaseManager.load_animaltfdb() 返回格式与 load_trrust()/load_chea3() 一致 ✓
   - 复用 TFEnrichmentAnalyzer 进行富集分析 ✓

3. **数据流完整性**:
   - Download(hTFtarget+AnimalTFDB) → Parse → OrthologMap → Build → Load → Analyze → Visualize → CLI ✓

4. **同源映射正确性**:
   - 排除自映射（TF不能调控自己） ✓
   - 人类TF→靶基因 → 目标物种TF→靶基因 ✓
   - 映射来源标注（HTML报告） ✓

---

## 执行选项

计划保存到 `docs/superpowers/plans/2026-05-31-animaltfdb-htftarget-support.md`。

**两种执行选项：**

**1. Subagent-Driven (推荐)** - 每个 Task 启动独立 subagent，任务间有检查点，快速迭代

**2. Inline Execution** - 在当前会话中批量执行，有检查点供审查

**选择哪种方式？**
