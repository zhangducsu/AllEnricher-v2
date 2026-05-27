# 数据库版本管理系统 — 实施完成报告

> **项目**: AllEnricher-v2
> **日期**: 2026-05-28
> **计划文档**: `docs/superpowers/plans/2026-05-27-database-version-management.md`
> **执行方式**: Subagent-Driven (executing-plans skill)
> **测试结果**: 637 passed, 2 deselected (slow), 0 failed

---

## 一、实施概览

本次实施完成了 14 项任务（合并为 8 个执行单元），为 AllEnricher-v2 建立了完整的数据库版本管理体系，覆盖远程更新检测、本地版本追溯、构建血缘追踪、版本锁定/切换、旧版本清理、分析结果版本记录等全链路功能。

### 完成状态总览

| 任务 | 描述 | 状态 |
|------|------|------|
| Task 1 | 创建版本管理核心模块（DatabaseVersion/VersionManifest） | ✅ 完成 |
| Task 2 | 实现远程版本检测器（RemoteVersionChecker） | ✅ 完成 |
| Task 3 | 集成版本记录到下载流程（gene2go/go_obo/gene_info/reactome/taxonomy） | ✅ 完成 |
| Task 4 | 构建流程写入 build_manifest.json | ✅ 完成 |
| Task 5 | 支持版本锁定（config.py + manager.py） | ✅ 完成 |
| Task 6 | 新增 check-update CLI 命令 | ✅ 完成 |
| Task 7 | 新增 cleanup CLI 命令 | ✅ 完成 |
| Task 8 | download 命令集成更新检测 + --force 参数 | ✅ 完成 |
| Task 9 | 端到端验证 | ✅ 通过 |
| Task 10 | 新增 list-versions CLI 命令（--lineage/--json） | ✅ 完成 |
| Task 11 | 构建血缘追踪增强（build_manifest source_versions 填充） | ✅ 完成 |
| Task 12 | 分析结果版本记录（TSV 注释头 + HTML 报告嵌入） | ✅ 完成 |
| Task 13 | 版本切换与回退（--use-version 参数） | ✅ 完成 |
| Task 14 | 端到端验证（完整版） | ✅ 通过 |

---

## 二、新增/修改文件清单

### 新增文件

| 文件 | 说明 |
|------|------|
| `allenricher/database/version.py` | 版本管理核心模块（~400行），包含 DatabaseVersion、VersionManifest、DatabaseVersionManager、RemoteVersionChecker |
| `tests/test_version.py` | 版本管理单元测试（11个测试），覆盖序列化、版本比较、清单读写、远程检测 |

### 修改文件

| 文件 | 修改内容 |
|------|---------|
| `allenricher/database/downloader.py` | 下载完成后自动记录版本元数据到 versions.json（gene2go/go_obo/gene_info/reactome/taxonomy） |
| `allenricher/database/builder.py` | 构建完成后写入 build_manifest.json（含完整依赖链和源数据版本） |
| `allenricher/database/manager.py` | 支持版本锁定（version 参数）、active_version 属性、get_build_metadata() 方法 |
| `allenricher/core/config.py` | 新增 `use_version: Optional[str]` 字段 |
| `allenricher/cli.py` | 新增 3 个子命令（check-update/cleanup/list-versions）、download --force 参数、analyze --use-version 参数、分析结果版本元数据传递 |
| `allenricher/core/enrichment.py` | save_results() 新增 metadata 参数，TSV 文件头部写入 `#` 注释版本信息 |
| `allenricher/report/generator.py` | HTML 报告动态版本号（allenricher.__version__）、可选嵌入数据库版本信息 |
| `database/organism/v20260515/hsa/build_manifest.json` | 回填 source_versions 字段（go_obo/gene2go/reactome 版本号） |

---

## 三、功能详解

### 3.1 远程版本检测

`RemoteVersionChecker` 类支持 7 个数据源的远程版本检测：

| 数据源 | 检测方式 | 版本格式示例 |
|--------|---------|-------------|
| NCBI gene2go.gz | HTTP HEAD Last-Modified | `Thu, 02 Apr 2026 02:53:36 GMT` |
| NCBI gene_info.gz | HTTP HEAD Last-Modified | 同上 |
| GO go-basic.obo | 文件内 data-version 字段 | `releases/2026-03-25` |
| EBI GOA proteomes | HTTP HEAD 目录页 | Last-Modified |
| KEGG | REST API info/kegg | `Release 118.0+/05-28, May 26` |
| Reactome | 下载页面 URL 解析 | `v96` |
| NCBI Taxonomy | HTTP HEAD Last-Modified | Last-Modified |

### 3.2 本地版本管理

`DatabaseVersionManager` 类提供：

- **versions.json 读写**: 记录每次下载的远程版本号、Last-Modified、下载时间
- **已安装版本扫描**: 自动扫描 `database/basic/{source}/` 和 `database/organism/v{date}/`
- **旧版本清理**: `find_stale_versions()` / `remove_stale_versions()` 支持 dry-run
- **构建血缘查询**: `get_build_lineage()` 读取 build_manifest.json
- **版本清单格式化**: `get_summary_table()` / `get_summary_json()` / `get_full_lineage_report()`

### 3.3 构建血缘追踪

每次 `build` 操作完成后，在物种数据库目录下生成 `build_manifest.json`：

```json
{
  "schema_version": "1.0",
  "built_at": "2026-05-15T00:00:00+00:00",
  "allenricher_version": "0.2.0",
  "species": "hsa",
  "taxid": "9606",
  "databases": ["DO", "GO", "KEGG", "Reactome"],
  "dependencies": {
    "GO": { "basic_dir": "basic/go/GO20260515", "files": ["gene2go.gz", ...] },
    "Reactome": { "basic_dir": "basic/reactome/Reactome20260515", ... },
    "KEGG": { "source": "REST API (real-time)", ... },
    "DO": { "basic_dir": "basic/do", ... }
  },
  "source_versions": {
    "go_obo": "releases/2026-03-25",
    "gene2go": "Fri, 15 May 2026 03:30:01 GMT",
    "reactome": "v96"
  }
}
```

### 3.4 新增 CLI 命令

```
# 检查远程数据源更新
allenricher check-update [--database-dir DIR] [--json]

# 查看本地已安装版本
allenricher list-versions [--database-dir DIR] [--json] [--lineage]

# 清理旧版本
allenricher cleanup [--keep N] [--dry-run] [--database-dir DIR]

# 强制重新下载（跳过更新检查）
allenricher download -d GO,KEGG --force

# 使用指定版本分析
allenricher analyze -i genes.txt -s hsa -d GO --use-version v20260515 -o results/
```

### 3.5 分析结果版本记录

TSV 输出文件头部自动嵌入版本注释：

```
# AllEnricher version: 2.0.3
# Analysis date: 2026-05-28T06:00:00+00:00
# Database version: v20260515
# Species: hsa
# Source data versions:
#   go_obo: releases/2026-03-25
#   gene2go: Fri, 15 May 2026 03:30:01 GMT
#   reactome: v96
#
Term_ID    Term_Name    ...
```

HTML 报告头部动态显示版本号，包含数据库版本信息。

---

## 四、已知限制与待改进项

| 项目 | 说明 | 优先级 |
|------|------|--------|
| DO/DisGeNET 下载未记录版本 | downloader.py 中 DO 数据下载后未调用 record_download() | 中 |
| check-update 无本地版本时全部显示"有更新" | 因为 versions.json 不在项目数据库目录中（运行时在 CWD 生成），首次使用需先 download | 低 |
| 单元测试覆盖不足 | find_stale_versions/remove_stale_versions/get_build_lineage 等方法缺少专门测试 | 中 |
| Reactome 版本记录代码重复 | downloader.py 中 Reactome 版本记录逻辑在两个分支中各写了一遍 | 低 |
| versions.json 路径依赖 CWD | DatabaseVersionManager 默认使用 `./database`，非项目安装目录时需 --database-dir | 低 |

---

## 五、测试验证结果

### 单元测试
```
tests/test_version.py: 11 passed (含 2 个 slow 网络测试)
```

### 全量测试
```
637 passed, 2 deselected (slow), 66 warnings, 482.29s
```

### CLI E2E 验证

| 命令 | 结果 |
|------|------|
| `allenricher --help` | ✅ 显示全部 11 个子命令（含新增的 check-update/cleanup/list-versions） |
| `allenricher list-versions` | ✅ 显示基础数据版本（5个数据源）+ 物种数据库版本（2个版本） |
| `allenricher list-versions --lineage` | ✅ 显示 v20260515/hsa 的完整构建血缘（依赖链 + 源数据版本） |
| `allenricher cleanup --dry-run --keep 1` | ✅ 预览将删除 GO20260515 和 v20260515 |
| `allenricher check-update` | ✅ 检测到 5 个数据源的远程版本信息 |
