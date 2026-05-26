# Spec: KEGG 在线数据获取优化

> **状态**: 待审批
> **优先级**: 高
> **影响范围**: `allenricher/database/` 新增模块, `builder.py`, `cli.py`, `downloader.py`

---

## 1. 背景与问题

### 1.1 当前状态

KEGG 在 v2 中是**唯一缺少自动数据获取**的数据库：

| 数据库 | 自动下载 | 自动构建 | 一键可用 |
|--------|----------|----------|----------|
| GO | ✅ | ✅ | ✅ |
| Reactome | ✅ | ✅ | ✅ |
| DO | ✅ | ✅ | ✅ |
| KEGG | ❌ | ⚠️ 需手动传参 | ❌ |

当前 `build_kegg()` 在缺少 `gene2pathway_path` 时**静默跳过**，不产出任何数据。

### 1.2 v1 的实现方式

v1 通过 `keggMapGrab.R` + `pathway2tab.pl` + `makeDB.kegg.v1.1.sh` 三层脚本协作：
1. **keggMapGrab.R**: 用 RCurl 爬取 KEGG 网页 HTML
2. **makeDB.kegg.v1.1.sh**: Perl 正则解析 HTML → 提取通路列表 → 为每个通路动态生成 R 代码 → 逐个下载通路详情 HTML → 正则提取基因 ID
3. **pathway2tab.pl**: 合并所有通路的基因列表为 0/1 矩阵

**v1 的瓶颈**：
- 串行逐通路爬取（人类 ~340 个通路，每个需一次 HTTP 请求）
- 网页 HTML 解析脆弱（依赖正则匹配 HTML 结构）
- R + Perl 多语言调用链复杂

### 1.3 KEGG REST API

KEG 官方提供 REST API（`https://rest.kegg.jp/`），比网页爬取更高效稳定：

| API 端点 | 用途 | 响应格式 |
|----------|------|----------|
| `list/pathway/{org}` | 获取物种所有通路 ID + 名称 | `pathway_id\tpathway_name` |
| `get/{org}{pathway_id}` | 获取通路详情（含基因列表） | 多行文本 |
| `link/{org}/pathway/{gene_id}` | 获取基因关联的通路 | `gene_id\tpathway_id` |
| `conv/{org}/ncbi-geneid` | KEGG 基因 ID ↔ NCBI Gene ID 映射 | `kegg_id\tncbi_id` |

**API 限制**：每秒不超过 10 次请求（约 0.1s/次），否则可能被临时封禁。

---

## 2. 优化方案

### 2.1 核心策略：REST API 替代网页爬取

用 KEGG REST API 替代 v1 的 HTML 爬取，优势：
- **速度提升 3-5x**：REST API 响应体积远小于 HTML 页面
- **解析简单可靠**：纯文本 TSV 格式，无需正则解析 HTML
- **代码精简**：Python 原生实现，消除 R + Perl 依赖

### 2.2 数据获取流程（新设计）

```
KEGG REST API
  │
  ├─ list/pathway/hsa          → 通路列表 (ID + 名称)
  │
  ├─ get/hsa{pathway_id} × N   → 每个通路的基因列表
  │   (串行，每请求间隔 0.12s)
  │
  └─ conv/hsa/ncbi-geneid      → KEGG ID ↔ NCBI Gene ID 映射
      (用于将 KEGG 基因 ID 转为 Gene Symbol)
      ↓
  gene2pathway.txt             → KEGGParser.build_database()
      ↓
  hsa.kegg2gene.tab.gz + hsa.kegg2disc.gz
```

### 2.3 速度优化措施

| 措施 | 说明 | 预期效果 |
|------|------|----------|
| **REST API 替代 HTML** | 响应体积减少 ~80% | 单请求时间缩短 60% |
| **批量 link 接口** | `link/hsa/pathway` 一次获取所有基因-通路关联 | 340 请求 → 2 请求 |
| **请求间隔控制** | 0.12s/次（略低于 10次/秒限制） | 避免被封禁 |
| **断点续传** | 已下载的通路跳过 | 中断后快速恢复 |
| **本地缓存** | 已获取数据存为 .tsv.gz | 重复构建零网络开销 |

### 2.4 批量 API 方案（最优）

**核心发现**：KEGG REST API 提供批量接口，可以大幅减少请求数：

```
# 方案 A：逐通路获取（v1 思路，~340 请求）
get/hsa00010  → 基因列表
get/hsa00020  → 基因列表
...（重复 340 次）

# 方案 B：批量 link 接口（新方案，仅 3 请求）
list/pathway/hsa                    → 所有通路 ID + 名称（1 请求）
conv/hsa/ncbi-geneid                → KEGG ID ↔ NCBI ID 映射（1 请求）
link/hsa/pathway                    → 所有基因-通路关联（1 请求）
```

**方案 B 仅需 3 次 API 请求**，总耗时 < 1 秒（vs v1 的 5-10 分钟）。

---

## 3. 技术设计

### 3.1 新增文件

```
allenricher/database/
├── kegg_fetcher.py      # 新增 - KEGG REST API 数据获取器
```

### 3.2 kegg_fetcher.py 接口设计

```python
class KEGGFetcher:
    """KEGG REST API 数据获取器"""

    BASE_URL = "https://rest.kegg.jp"
    REQUEST_INTERVAL = 0.12  # 秒，略低于 10次/秒限制

    def __init__(self, cache_dir: str, overwrite: bool = False):
        ...

    def fetch_species_data(
        self,
        species: str,          # KEGG 物种代码，如 "hsa"
        gene_info_path: str,   # gene_info.gz 路径（用于 ID 映射）
    ) -> Tuple[str, str]:
        """获取物种 KEGG 数据

        Returns:
            (gene2pathway_path, pathway_summary_path)
            两个文件的路径，可直接传给 KEGGParser.build_database()
        """
        ...

    def _list_pathways(self, species: str) -> List[Tuple[str, str]]:
        """获取物种所有通路列表

        API: list/pathway/{species}
        Returns: [(pathway_id, pathway_name), ...]
        """
        ...

    def _get_gene_pathway_links(self, species: str) -> Dict[str, List[str]]:
        """获取所有基因-通路关联（批量接口）

        API: link/{species}/pathway
        Returns: {gene_id: [pathway_id, ...], ...}
        """
        ...

    def _get_kegg_ncbi_mapping(self, species: str) -> Dict[str, str]:
        """获取 KEGG ID → NCBI Gene ID 映射

        API: conv/{species}/ncbi-geneid
        Returns: {kegg_id: ncbi_gene_id, ...}
        """
        ...

    def _ncbi_id_to_symbol(self, gene_info_path: str) -> Dict[str, str]:
        """从 gene_info.gz 构建 NCBI Gene ID → Symbol 映射"""
        ...
```

### 3.3 修改文件

#### builder.py - build_kegg() 方法

```python
def build_kegg(self, species: str, taxid: int, ...):
    # 新增：如果 gene2pathway_path 为 None，自动调用 KEGGFetcher
    if gene2pathway_path is None:
        from .kegg_fetcher import KEGGFetcher
        fetcher = KEGGFetcher(cache_dir=str(self.basic_dir / "kegg"))
        gene2pathway_path, pathway_summary_path = fetcher.fetch_species_data(
            species=species,
            gene_info_path=str(gene_info_path),
        )
    # 后续逻辑不变...
```

#### downloader.py - 新增 download_kegg()

```python
def download_kegg(self, species: str) -> str:
    """下载 KEGG 数据（通过 REST API）"""
    ...
```

#### cli.py - download 子命令支持 kegg

```python
# download -d kegg 自动触发 KEGG REST API 获取
```

### 3.4 输出文件格式

#### gene2pathway.txt（与现有 KEGGParser 兼容）
```
#gene_symbol\tentrez_id\tpathway_id\tpathway_name
TP53\t7157\thsa04110\tCell cycle
BRCA1\t672\thsa05212\tBreast cancer
...
```

#### pathway_summary.txt（可选，与现有 KEGGParser 兼容）
```
#Category\tSubcategory\tpathway_id\tpathway_name\turl
Metabolism\tCarbohydrate_metabolism\t00010\tGlycolysis_/_Gluconeogenesis\thttps://www.kegg.jp/entry/hsa00010
...
```

> **注意**：KEGG REST API 的 `list/pathway/{org}` 返回的通路名称中包含分类信息（以缩进表示层级），需要解析为 Category + Subcategory。

---

## 4. 性能预期

| 指标 | v1 (HTML爬取) | v2 (REST API) | 提升 |
|------|---------------|---------------|------|
| API 请求数 | ~340+ | 3 | **113x** |
| 总耗时（人类） | 5-10 分钟 | < 5 秒 | **60-120x** |
| 网络流量 | ~50MB HTML | ~2MB 文本 | **25x** |
| 代码依赖 | R + Perl | 纯 Python | 消除外部依赖 |
| 解析可靠性 | HTML 正则（脆弱） | TSV 文本（稳定） | 质的飞跃 |

---

## 5. 风险与缓解

| 风险 | 概率 | 影响 | 缓解措施 |
|------|------|------|----------|
| KEGG REST API 不可用 | 低 | 高 | 保留手动传入 gene2pathway_path 的回退机制 |
| API 频率限制触发 | 中 | 中 | 0.12s 间隔 + 自动重试 + 指数退避 |
| 批量 link 接口数据不全 | 低 | 中 | 对比 list/pathway 结果，缺失通路逐个补抓 |
| KEGG ID → Symbol 映射缺失 | 低 | 低 | 使用 gene_info.gz 中的 Entrez ID → Symbol 映射 |

---

## 6. 实施任务清单

1. **创建 `kegg_fetcher.py`** - KEGG REST API 数据获取器
   - `_list_pathways()`: 获取通路列表
   - `_get_gene_pathway_links()`: 批量获取基因-通路关联
   - `_get_kegg_ncbi_mapping()`: KEGG ID ↔ NCBI ID 映射
   - `fetch_species_data()`: 编排上述步骤，输出标准格式文件

2. **修改 `builder.py`** - `build_kegg()` 自动调用 KEGGFetcher

3. **修改 `downloader.py`** - 新增 `download_kegg()` 方法

4. **修改 `cli.py`** - download 子命令支持 `kegg`

5. **编写测试** - 验证 API 调用、数据解析、格式兼容性

6. **集成验证** - `allenricher download -d kegg && allenricher build -s hsa -t 9606 -d KEGG`
