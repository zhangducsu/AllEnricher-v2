# AllEnricher v2 统一物种注册表架构 — 实施总结报告

**日期**: 2026-05-27  
**测试结果**: 628/628 通过 ✅  
**新增测试**: 16 个（species_registry: 11, goa_fetcher: 5）

---

## 一、实施概要

本次实施完成了统一物种注册表架构的全部核心模块，实现了以下目标：

1. **统一物种注册表** (`SpeciesRegistry`): 整合 GO、KEGG、Reactome、DO 四个数据库的物种支持信息
2. **UniProt GOA proteomes 集成** (`GOAFetcher`): 作为 GO 注释的扩展数据源，按需下载
3. **分层构建策略** (`DatabaseBuilder`): gene2go 优先，GOA 自动回退
4. **CLI 查询命令**: `list-species` 和 `query-species` 供用户查询物种支持情况
5. **背景基因集策略更新** (`EnrichmentAnalyzer`): 默认使用有注释的基因（clusterProfiler 方案）

---

## 二、文件变更清单

### 新增文件（2 个）

| 文件 | 说明 | 代码行数 |
|------|------|----------|
| `allenricher/database/species_registry.py` | 统一物种注册表管理模块 | ~350 行 |
| `allenricher/database/goa_fetcher.py` | UniProt GOA 按需获取模块 | ~250 行 |
| `tests/test_species_registry.py` | SpeciesRegistry 单元测试 | ~180 行 |
| `tests/test_goa_fetcher.py` | GOAFetcher 单元测试 | ~130 行 |

### 修改文件（5 个）

| 文件 | 修改内容 | 新增方法数 |
|------|----------|-----------|
| `allenricher/database/downloader.py` | 新增 8 个注册表构建方法 + 1 个辅助方法 | 9 |
| `allenricher/database/builder.py` | 新增 GOA 构建路径 + 分层回退逻辑 | 8 |
| `allenricher/database/kegg_fetcher.py` | 新增 `fetch_organism_list()` 方法 | 1 |
| `allenricher/cli.py` | 新增 `list-species` 和 `query-species` 命令 | 2 函数 + 2 子命令 |
| `allenricher/core/enrichment.py` | 新增 `get_annotated_genes()` 和 `resolve_background()` | 2 |

---

## 三、核心模块说明

### 3.1 SpeciesRegistry（统一物种注册表）

**数据模型**: `SpeciesEntry` 包含 21 个字段，覆盖四个数据库的支持状态和统计信息。

**核心能力**:
- 按 TaxID / 拉丁名 / KEGG 代码查询
- 按数据库支持状态筛选
- 自动生成 KEGG 缩写（属名前3字母 + 种名首字母）
- 统计汇总（各数据库物种数、基因数、条目数）

**文件格式**: `supported_species.tsv`（制表符分隔，21 列）

### 3.2 GOAFetcher（GOA 按需获取）

**核心能力**:
- 从 EBI FTP 按需下载指定物种的 GOA 文件
- 解析 GAF 2.2 格式（17 列），提取 Gene Symbol + GO ID
- 生成与现有 `GOParser.parse_gene2go` 格式一致的输出文件
- 缓存机制避免重复下载

### 3.3 DataDownloader（注册表构建）

**新增方法**:

| 方法 | 功能 |
|------|------|
| `_build_go_registry()` | 从 gene2go.gz 提取 GO 物种注册表 |
| `_download_goa_index()` | 从 EBI FTP 获取 UniProt GOA 物种索引 |
| `_merge_go_registries()` | 合并 gene2go 和 GOA 注册表 |
| `_build_kegg_registry()` | 从 KEGG API 获取 KEGG 物种注册表 |
| `_build_reactome_registry()` | 从 NCBI2Reactome 提取 Reactome 物种注册表 |
| `_build_do_registry()` | 生成 DO 物种注册表（仅人类） |
| `_merge_all_registries()` | 合并为统一 `supported_species.tsv` |
| `_report_download_summary()` | 打印统计摘要 |

### 3.4 DatabaseBuilder（分层构建）

**新增方法**:

| 方法 | 功能 |
|------|------|
| `build_go_from_goa()` | 从 GOA 文件构建 GO 数据库 |
| `build_go_with_fallback()` | 先尝试 gene2go，失败自动切换 GOA |
| `_check_gene2go_has_taxid()` | 快速检查 gene2go 是否包含指定物种 |
| `_get_species_dir()` | 生成 `taxid.Latin_Name` 格式目录名 |
| `_get_species_prefix()` | 生成文件名前缀 |

### 3.5 CLI 新增命令

**`allenricher list-species`**:
```
allenricher list-species                    # 列出所有物种
allenricher list-species --go               # 只列出支持 GO 的物种
allenricher list-species --summary          # 统计摘要
allenricher list-species --format tsv       # TSV 格式输出
```

**`allenricher query-species`**:
```
allenricher query-species --taxid 9606      # 按 TaxID 查询
allenricher query-species --name "Mus"      # 按拉丁名模糊查询
allenricher query-species --kegg hsa        # 按 KEGG 代码查询
```

### 3.6 背景基因集策略

| 模式 | 说明 |
|------|------|
| `annotated`（默认） | 使用数据库中有注释的基因（clusterProfiler 方案） |
| `genome` | 使用全基因组基因（v1/v2 原方案） |
| `custom` | 使用用户自定义背景基因列表 |

---

## 四、目录结构（目标）

```
database/
├── supported_species.tsv             # 统一物种注册表（主文件）
├── basic/
│   ├── go/GO{date}/
│   │   ├── go_species_registry.tsv   # GO 专用注册表
│   │   ├── gene2go.gz / gene_info.gz / go-basic.obo
│   │   └── goa/                      # GOA 缓存（build 阶段填充）
│   ├── kegg/
│   │   ├── kegg_species_registry.tsv # KEGG 专用注册表
│   │   └── cache/                    # KEGG 缓存
│   ├── reactome/Reactome{date}/
│   │   └── reactome_species_registry.tsv
│   └── do/
│       └── do_species_registry.tsv
└── organism/v{date}/
    ├── 9606.Homo_sapiens/            # taxid.拉丁名 命名
    │   ├── 9606.Homo_sapiens.GO2gene.tab.gz
    │   ├── 9606.Homo_sapiens.KEGG2gene.tab.gz
    │   └── ...
    └── ...
```

---

## 五、测试结果

### 新增测试（16 个）

| 测试文件 | 测试数 | 状态 |
|----------|--------|------|
| `test_species_registry.py` | 11 | ✅ 全部通过 |
| `test_goa_fetcher.py` | 5 | ✅ 全部通过 |

### 回归测试

| 指标 | 结果 |
|------|------|
| 总测试数 | 628 |
| 通过 | 628 |
| 失败 | 0 |
| 跳过 | 0 |
| 耗时 | 368 秒 |

---

## 六、待完成事项

以下功能已实现代码框架，但尚未集成到主流程中：

1. **download 命令集成**: `download_all()` 方法尚未调用新增的注册表构建方法
2. **build 命令集成**: `build_species_db()` 尚未调用 `build_go_with_fallback()`
3. **CLI --background-mode 参数**: `analyze` 命令尚未暴露 `background_mode` 参数
4. **species_lookup.py 废弃**: 现有 `SpeciesLookup` 类的功能应逐步迁移到 `SpeciesRegistry`
5. **config.py SPECIES_CONFIGS**: 硬编码的 16 物种配置应改为从注册表动态加载
6. **E2E 集成测试**: 需要编写从 download → registry → build → analyze 的完整 E2E 测试

这些待完成事项属于集成层面的工作，不影响现有功能的正常运行。
