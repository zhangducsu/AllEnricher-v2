# AllEnricher v2 统一物种注册表 — 集成完善实施总结报告

**日期**: 2026-05-27  
**测试结果**: 628/628 通过 ✅  
**实施耗时**: ~30 分钟

---

## 一、实施概要

本次实施完成了上一轮遗留的 8 个集成缺失点（I-1 到 I-8）的全部修复工作，使统一物种注册表架构真正可用。

### 核心成果

| 编号 | 问题 | 状态 | 关键修改 |
|------|------|------|----------|
| I-1 | `download_all()` 未编排注册表构建流水线 | ✅ 已修复 | downloader.py: 下载后自动构建所有注册表 |
| I-2 | `build_species_db()` 未调用 `build_go_with_fallback()` | ✅ 已修复 | builder.py: GO 构建自动回退到 GOA 数据源 |
| I-3 | `analyze` 命令未暴露 `--background-mode` 参数 | ✅ 已修复 | cli.py: 新增 `--background-mode` 参数 |
| I-4 | `cmd_analyze()` 未使用 `resolve_background()` | ✅ 已修复 | cli.py: 背景基因集逻辑使用新机制 |
| I-5 | `Config.validate()` 硬编码依赖 SPECIES_CONFIGS | ✅ 已修复 | config.py: 验证时同时查询注册表 |
| I-6 | `cmd_list()` 硬编码使用 SPECIES_CONFIGS | ✅ 已修复 | cli.py: 优先从注册表读取物种列表 |
| I-7 | `species_lookup.py` 是死代码 | ✅ 已修复 | species_lookup.py: 添加废弃警告 |
| I-8 | `cmd_build()` 未使用 SpeciesRegistry | ✅ 已修复 | cli.py: 构建前查询注册表 |

---

## 二、详细变更

### 2.1 downloader.py — 注册表构建流水线集成

**修改位置**: `download_all()` 方法

**新增逻辑**:
```python
# 下载完成后，编排注册表构建
if go_dir and gene2go_path.exists():
    go_registry = self._build_go_registry(gene2go_path, go_dir)
    goa_index = self._download_goa_index(go_dir)
    if go_registry and goa_index:
        self._merge_go_registries(go_registry, goa_index, go_dir)

if kegg_dir:
    self._build_kegg_registry(kegg_dir)

if reactome_dir and ncbi2reactome_path.exists():
    self._build_reactome_registry(ncbi2reactome_path, reactome_dir)

if do_dir:
    self._build_do_registry(do_dir)

# 合并为统一注册表
self._merge_all_registries(...)
self._report_download_summary(...)
```

**特性**:
- 每步添加 try/except，单步失败不阻塞
- 仅在对应数据库下载成功时才构建其注册表
- 最后打印统计摘要到控制台

### 2.2 builder.py — GOA 回退逻辑集成

**修改位置**: `build_species_db()` 方法

**变更内容**:
- 将 `self.build_go(species, taxid, go_version)` 改为 `self.build_go_with_fallback(taxid, latin_name, go_version)`
- 自动获取 latin_name：优先从 SpeciesRegistry，回退到 SPECIES_CONFIGS

**回退逻辑**:
1. 检查 supported_species.tsv 中的 go_source
2. 如果为 "ncbi_gene2go" 或 "both"，使用 gene2go 路径
3. 如果为 "uniprot_goa"，使用 GOA 路径
4. 如果注册表不可用，快速扫描 gene2go.gz 前 10000 行判断

### 2.3 config.py — 物种验证放宽

**修改位置**: `validate()` 方法

**变更内容**:
```python
# 原有检查
if self.species in SPECIES_CONFIGS:
    return

# 新增：尝试从注册表查询
try:
    from ..database.species_registry import SpeciesRegistry
    registry = SpeciesRegistry.load_default()
    if registry.query_by_kegg_code(self.species) or \
       (self.species.isdigit() and registry.query_by_taxid(int(self.species))):
        return
except:
    pass

# 都不存在才报错
raise ValueError(f"Unknown species: {self.species}")
```

### 2.4 cli.py — 多项 CLI 集成

#### 2.4.1 analyze 命令 --background-mode 参数

**新增参数**:
```bash
allenricher analyze --background-mode {annotated,genome,custom}
```

**逻辑**:
- `annotated`（默认）：使用有注释的基因作为背景
- `genome`：使用全基因组基因作为背景
- `custom`：要求用户必须提供 `--background` 参数

#### 2.4.2 list 命令使用 SpeciesRegistry

**变更**: `cmd_list()` 优先从 SpeciesRegistry 加载物种列表，失败时回退到 SPECIES_CONFIGS

#### 2.4.3 build 命令使用 SpeciesRegistry

**变更**: `cmd_build()` 在构建前查询 SpeciesRegistry，打印物种在各数据库的支持状态

### 2.5 species_lookup.py — 标记废弃

**修改**: 模块顶部添加废弃警告
```python
import warnings
warnings.warn(
    "species_lookup is deprecated, use species_registry instead",
    DeprecationWarning,
    stacklevel=2
)
```

### 2.6 test_enrichment_extended.py — 测试修复

**修改**: 更新 `TestSpeciesLookup` 类，使用临时注册表而非依赖 `load_default()`

---

## 三、测试结果

### 全量回归测试

```
================ 628 passed, 66 warnings in 322.86s (0:05:22) =================
```

| 指标 | 结果 |
|------|------|
| 总测试数 | 628 |
| 通过 | 628 ✅ |
| 失败 | 0 |
| 跳过 | 0 |
| 耗时 | 5 分 22 秒 |

### 关键测试覆盖

- `test_species_registry.py` — 11 个测试 ✅
- `test_goa_fetcher.py` — 5 个测试 ✅
- `test_enrichment_extended.py::TestSpeciesLookup` — 6 个测试 ✅

---

## 四、使用示例

### 下载并构建注册表

```bash
# 下载数据（会自动构建注册表）
allenricher download -d go,kegg,reactome,do

# 查看统计摘要
allenricher list-species --summary
```

### 查询物种支持情况

```bash
# 列出所有支持 GO 的物种
allenricher list-species --go

# 查询特定物种
allenricher query-species --taxid 9606
allenricher query-species --kegg hsa
```

### 构建物种数据库（自动回退）

```bash
# 对于 gene2go 支持的物种，使用 gene2go
allenricher build -s hsa -t 9606 -d GO

# 对于仅 GOA 支持的物种，自动回退到 GOA
allenricher build -s ath -t 3702 -d GO
```

### 使用新的背景基因集策略

```bash
# 默认：使用有注释的基因（clusterProfiler 方案）
allenricher analyze -i degs.txt -s hsa -d GO

# 显式指定 annotated 模式
allenricher analyze -i degs.txt -s hsa -d GO --background-mode annotated

# 使用全基因组背景
allenricher analyze -i degs.txt -s hsa -d GO --background-mode genome

# 使用自定义背景
allenricher analyze -i degs.txt -s hsa -d GO --background-mode custom --background my_genes.txt
```

---

## 五、架构完整性

### 数据流

```
下载阶段 (allenricher download)
├── 下载基础数据 (gene2go, gene_info, obo, ...)
└── 构建注册表
    ├── go_species_registry.tsv
    ├── kegg_species_registry.tsv
    ├── reactome_species_registry.tsv
    ├── do_species_registry.tsv
    └── supported_species.tsv (统一注册表)

查询阶段 (allenricher list-species / query-species)
└── 读取 supported_species.tsv

构建阶段 (allenricher build)
├── 查询 supported_species.tsv
├── 确定数据源 (gene2go vs GOA)
└── 构建物种数据库

分析阶段 (allenricher analyze)
├── 验证物种 (SPECIES_CONFIGS + supported_species.tsv)
├── 确定背景基因集 (annotated/genome/custom)
└── 执行富集分析
```

### 向后兼容

| 场景 | 行为 |
|------|------|
| 无 supported_species.tsv | 回退到 SPECIES_CONFIGS |
| SpeciesRegistry 加载失败 | 静默回退，不影响功能 |
| 原有 16 种物种 | 完全兼容，验证通过 |
| 旧版 CLI 用法 | 完全兼容，参数不变 |

---

## 六、文件变更清单

| 文件 | 操作 | 说明 |
|------|------|------|
| `allenricher/database/downloader.py` | 修改 | 集成注册表构建流水线 |
| `allenricher/database/builder.py` | 修改 | 集成 GOA 回退逻辑 |
| `allenricher/core/config.py` | 修改 | 放宽物种验证 |
| `allenricher/cli.py` | 修改 | 集成 CLI 命令 |
| `allenricher/database/species_lookup.py` | 修改 | 添加废弃警告 |
| `tests/test_enrichment_extended.py` | 修改 | 修复测试 |

---

## 七、总结

本次实施完成了统一物种注册表架构的全部集成工作，实现了：

1. **全自动注册表构建**: `download` 命令下载完成后自动生成所有注册表
2. **智能数据源回退**: `build` 命令自动选择 gene2go 或 GOA 数据源
3. **灵活的背景基因集**: `analyze` 命令支持 annotated/genome/custom 三种模式
4. **完整的 CLI 查询**: `list-species` 和 `query-species` 命令可用
5. **向后兼容**: 所有修改保持向后兼容，628 个测试全部通过

统一物种注册表架构现已完全可用。
