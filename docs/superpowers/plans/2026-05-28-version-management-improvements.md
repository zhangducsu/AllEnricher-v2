# 版本管理改进计划

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** 修复数据库版本管理系统中已知的 5 项限制，补充单元测试覆盖，完成全量端到端验证

**Architecture:** 针对每个限制进行精准修复：DO 下载补版本记录、Reactome 重复代码提取、补充 12 个未测试方法的单元测试、统一 CLI 的 db_dir 解析逻辑

**Tech Stack:** Python, pytest, requests (mock)

---

## 现状问题清单

| # | 问题 | 文件 | 行号 | 优先级 |
|---|------|------|------|--------|
| 1 | DO 下载未记录版本 | downloader.py | 294-340 | 中 |
| 2 | Reactome 版本记录代码重复 | downloader.py | 241-255 vs 273-287 | 低 |
| 3 | 单元测试覆盖不足（12 个方法无测试） | test_version.py | - | 中 |
| 4 | versions.json 路径依赖 CWD | cli.py | 1214, 1267, 1314 | 低 |
| 5 | check-update 无本地版本时全部显示"有更新" | cli.py + version.py | - | 低 |

---

## 文件结构

```
allenricher/database/
├── downloader.py           ← 修改：DO 版本记录 + Reactome 去重
├── version.py              ← 不变
allenricher/cli.py           ← 修改：统一 db_dir 解析
tests/test_version.py        ← 修改：补充 ~15 个测试
```

---

## Task 1: DO 下载版本记录 + Reactome 去重

**Files:**
- Modify: `allenricher/database/downloader.py:241-287, 294-340`

- [ ] **Step 1: 提取 Reactome 版本记录为私有方法**

在 `DataDownloader` 类中新增方法（建议放在 `download_reactome_basic()` 之后）：

```python
def _record_reactome_version(self, version: str) -> None:
    """记录 Reactome 版本元数据到 versions.json"""
    try:
        from allenricher.database.version import DatabaseVersionManager, RemoteVersionChecker
        _vm = DatabaseVersionManager(database_dir=str(self.root_dir))
        _checker = RemoteVersionChecker()
        _re_info = _checker.check_reactome_version()
        if _re_info:
            _vm.record_download(
                source="reactome", local_version=version,
                local_path=f"basic/reactome/{version}",
                remote_version=_re_info.get("remote_version"),
                remote_last_modified=_re_info.get("last_modified"),
            )
    except Exception as _e:
        logger.warning("记录 Reactome 版本元数据失败: %s", _e)
```

- [ ] **Step 2: 替换 download_reactome_basic() 中的两处重复代码**

将第 241-255 行的 try 块替换为：
```python
self._record_reactome_version(version)
```

将第 273-287 行的 try 块替换为：
```python
self._record_reactome_version(version)
```

- [ ] **Step 3: 在 download_do_files() 末尾添加版本记录**

在 `download_do_files()` 的 `return files` 之前添加：

```python
# 记录 DO 版本元数据到 versions.json
try:
    from allenricher.database.version import DatabaseVersionManager, RemoteVersionChecker
    _vm = DatabaseVersionManager(database_dir=str(self.root_dir))
    _vm.record_download(
        source="do",
        local_version="cached",
        local_path="basic/do",
        downloaded_at=datetime.now(timezone.utc).isoformat(),
    )
except Exception as _e:
    logger.warning("记录 DO 版本元数据失败: %s", _e)
```

注意：需要在文件顶部确认 `datetime` 和 `timezone` 已导入（检查现有 import）。

- [ ] **Step 4: 验证语法正确**

Run: `cd f:\OneDrive\Documents\TraeSOLO\AllEnricher\AllEnricher-v2 && python -c "from allenricher.database.downloader import DataDownloader; print('OK')"`
Expected: `OK`

- [ ] **Step 5: Commit**

```bash
git add allenricher/database/downloader.py
git commit -m "fix(download): DO 版本记录 + Reactome 版本记录去重"
```

---

## Task 2: 统一 CLI db_dir 解析

**Files:**
- Modify: `allenricher/cli.py:1214, 1267, 1314`

- [ ] **Step 1: 添加统一的 db_dir 解析辅助函数**

在 `cli.py` 中（函数定义区域，`cmd_analyze` 之前）添加：

```python
def _resolve_db_dir(args) -> str:
    """统一解析数据库目录路径

    优先级：CLI --database-dir > config.database_dir > 默认 ./database
    """
    if hasattr(args, 'database_dir') and args.database_dir:
        return args.database_dir
    # 对于非 analyze 命令，没有 config 对象，直接用默认值
    return "./database"
```

- [ ] **Step 2: 替换 cmd_check_update 中的 db_dir 解析**

第 1214 行：
```python
# 原: db_dir = args.database_dir or "./database"
db_dir = _resolve_db_dir(args)
```

- [ ] **Step 3: 替换 cmd_cleanup 中的 db_dir 解析**

第 1267 行：
```python
# 原: db_dir = args.database_dir or "./database"
db_dir = _resolve_db_dir(args)
```

- [ ] **Step 4: 替换 cmd_list_versions 中的 db_dir 解析**

第 1314 行：
```python
# 原: db_dir = args.database_dir or "./database"
db_dir = _resolve_db_dir(args)
```

- [ ] **Step 5: 验证**

Run: `cd f:\OneDrive\Documents\TraeSOLO\AllEnricher\AllEnricher-v2 && python -m allenricher check-update --help`
Expected: 正常显示帮助信息

- [ ] **Step 6: Commit**

```bash
git add allenricher/cli.py
git commit -m "refactor(cli): 统一 db_dir 解析逻辑"
```

---

## Task 3: 补充单元测试

**Files:**
- Modify: `tests/test_version.py`

- [ ] **Step 1: 补充 DatabaseVersionManager 测试**

在 `TestDatabaseVersionManager` 类中追加以下测试：

```python
def test_list_local_versions(self, tmp_path):
    """测试 list_local_versions 返回所有已记录版本"""
    mgr = DatabaseVersionManager(str(tmp_path))
    mgr.record_download("go", "GO20260527", "basic/go/GO20260527")
    mgr.record_download("kegg", "cached", "basic/kegg")
    versions = mgr.list_local_versions()
    assert "go" in versions
    assert "kegg" in versions
    assert len(versions) == 2

def test_get_organism_build_info(self, tmp_path):
    """测试 get_organism_build_info 返回版本下的物种列表"""
    org_dir = tmp_path / "organism" / "v20260515"
    org_dir.mkdir(parents=True)
    (org_dir / "hsa").mkdir()
    (org_dir / "mmu").mkdir()
    mgr = DatabaseVersionManager(str(tmp_path))
    info = mgr.get_organism_build_info("v20260515")
    assert "v20260515" in info
    assert set(info["v20260515"]) == {"hsa", "mmu"}

def test_find_stale_versions(self, tmp_path):
    """测试 find_stale_versions 正确识别旧版本"""
    # 创建 basic/go/ 下 3 个版本目录
    for ver in ["GO20260501", "GO20260515", "GO20260527"]:
        (tmp_path / "basic" / "go" / ver).mkdir(parents=True)
    # 创建 organism/ 下 3 个版本目录
    for ver in ["v20260501", "v20260515", "v20260527"]:
        org_dir = tmp_path / "organism" / ver
        org_dir.mkdir(parents=True)
        (org_dir / "hsa").mkdir()

    mgr = DatabaseVersionManager(str(tmp_path))
    stale = mgr.find_stale_versions(keep_count=2)
    assert "go" in stale
    assert "GO20260501" in stale["go"]
    assert "GO20260527" not in stale["go"]
    assert "organism" in stale
    assert "v20260501" in stale["organism"]

def test_remove_stale_versions_dry_run(self, tmp_path):
    """测试 remove_stale_versions dry-run 不删除文件"""
    for ver in ["GO20260501", "GO20260515"]:
        (tmp_path / "basic" / "go" / ver).mkdir(parents=True)
    mgr = DatabaseVersionManager(str(tmp_path))
    removed = mgr.remove_stale_versions(keep_count=1, dry_run=True)
    assert "GO20260501" in removed["go"]
    # 文件仍存在
    assert (tmp_path / "basic" / "go" / "GO20260501").exists()

def test_remove_stale_versions_actual(self, tmp_path):
    """测试 remove_stale_versions 实际删除"""
    for ver in ["GO20260501", "GO20260515"]:
        (tmp_path / "basic" / "go" / ver).mkdir(parents=True)
    mgr = DatabaseVersionManager(str(tmp_path))
    removed = mgr.remove_stale_versions(keep_count=1, dry_run=False)
    assert "GO20260501" in removed["go"]
    # 文件已删除
    assert not (tmp_path / "basic" / "go" / "GO20260501").exists()
    # 最新版保留
    assert (tmp_path / "basic" / "go" / "GO20260515").exists()

def test_get_build_lineage(self, tmp_path):
    """测试 get_build_lineage 读取 build_manifest.json"""
    manifest_dir = tmp_path / "organism" / "v20260515" / "hsa"
    manifest_dir.mkdir(parents=True)
    manifest_data = {
        "built_at": "2026-05-15T00:00:00+00:00",
        "species": "hsa",
        "source_versions": {"go_obo": "releases/2026-03-25"},
    }
    import json
    with open(manifest_dir / "build_manifest.json", "w") as f:
        json.dump(manifest_data, f)

    mgr = DatabaseVersionManager(str(tmp_path))
    lineage = mgr.get_build_lineage("v20260515", "hsa")
    assert lineage is not None
    assert lineage["species"] == "hsa"
    assert lineage["source_versions"]["go_obo"] == "releases/2026-03-25"

def test_get_build_lineage_not_found(self, tmp_path):
    """测试 get_build_lineage 不存在时返回 None"""
    mgr = DatabaseVersionManager(str(tmp_path))
    lineage = mgr.get_build_lineage("v99999999", "hsa")
    assert lineage is None

def test_get_summary_json(self, tmp_path):
    """测试 get_summary_json 返回结构化数据"""
    mgr = DatabaseVersionManager(str(tmp_path))
    mgr.record_download("go", "GO20260527", "basic/go/GO20260527")
    summary = mgr.get_summary_json()
    assert "basic_versions" in summary
    assert "organism_versions" in summary
    assert "version_records" in summary
    assert "go" in summary["version_records"]

def test_get_summary_table(self, tmp_path):
    """测试 get_summary_table 返回非空字符串"""
    for ver in ["GO20260527"]:
        (tmp_path / "basic" / "go" / ver).mkdir(parents=True)
    mgr = DatabaseVersionManager(str(tmp_path))
    table = mgr.get_summary_table()
    assert "本地数据库版本清单" in table
    assert "go" in table

def test_get_full_lineage_report(self, tmp_path):
    """测试 get_full_lineage_report 生成血缘报告"""
    manifest_dir = tmp_path / "organism" / "v20260515" / "hsa"
    manifest_dir.mkdir(parents=True)
    import json
    with open(manifest_dir / "build_manifest.json", "w") as f:
        json.dump({"built_at": "2026-05-15", "species": "hsa", "databases": ["GO"]}, f)

    mgr = DatabaseVersionManager(str(tmp_path))
    report = mgr.get_full_lineage_report()
    assert "构建血缘追踪报告" in report
    assert "v20260515/hsa" in report
```

- [ ] **Step 2: 补充 RemoteVersionChecker 测试**

在 `TestRemoteVersionChecker` 类中追加以下测试（使用 mock 避免网络依赖）：

```python
from unittest.mock import patch, MagicMock

def test_check_head_success(self):
    """测试 check_head 成功解析 Last-Modified"""
    checker = RemoteVersionChecker()
    mock_resp = MagicMock()
    mock_resp.headers = {"Last-Modified": "Thu, 02 May 2026 02:53:36 GMT", "ETag": '"abc123"'}
    mock_resp.raise_for_status = MagicMock()
    with patch("requests.head", return_value=mock_resp):
        result = checker.check_head("https://example.com/file.gz")
    assert result is not None
    assert result["last_modified"] == "Thu, 02 May 2026 02:53:36 GMT"
    assert result["etag"] == '"abc123"'

def test_check_head_failure(self):
    """测试 check_head 失败返回 None"""
    checker = RemoteVersionChecker()
    with patch("requests.head", side_effect=Exception("Connection error")):
        result = checker.check_head("https://example.com/file.gz")
    assert result is None

def test_check_updates_no_local(self, tmp_path):
    """测试 check_updates 本地无记录时全部显示有更新"""
    checker = RemoteVersionChecker()
    mgr = DatabaseVersionManager(str(tmp_path))
    # mock check_all_sources 返回固定结果
    mock_results = {
        "gene2go": {"last_modified": "Thu, 02 May 2026 02:53:36 GMT"},
        "go_obo": {"remote_version": "releases/2026-05-01"},
    }
    with patch.object(checker, "check_all_sources", return_value=mock_results):
        status = checker.check_updates(mgr)
    for source, info in status.items():
        assert info["has_update"] is True

def test_check_updates_up_to_date(self, tmp_path):
    """测试 check_updates 本地已是最新时显示无更新"""
    mgr = DatabaseVersionManager(str(tmp_path))
    mgr.record_download(
        source="gene2go", local_version="GO20260527",
        local_path="basic/go/GO20260527",
        remote_last_modified="Thu, 02 May 2026 02:53:36 GMT",
    )
    checker = RemoteVersionChecker()
    mock_results = {
        "gene2go": {"last_modified": "Thu, 02 May 2026 02:53:36 GMT"},
    }
    with patch.object(checker, "check_all_sources", return_value=mock_results):
        status = checker.check_updates(mgr)
    assert status["gene2go"]["has_update"] is False
```

注意：需要在文件顶部添加 `from unittest.mock import patch, MagicMock`。

- [ ] **Step 3: 运行测试验证**

Run: `cd f:\OneDrive\Documents\TraeSOLO\AllEnricher\AllEnricher-v2 && python -m pytest tests/test_version.py -v`
Expected: 全部通过（原 11 个 + 新增 ~15 个）

- [ ] **Step 4: Commit**

```bash
git add tests/test_version.py
git commit -m "test(version): 补充版本管理单元测试覆盖"
```

---

## Task 4: 全量端到端测试

**Files:**
- Test: 手动 E2E 测试

- [ ] **Step 1: 运行全量 pytest**

Run: `cd f:\OneDrive\Documents\TraeSOLO\AllEnricher\AllEnricher-v2 && python -m pytest tests/ -v --tb=short -m "not slow" 2>&1 | tail -30`
Expected: 全部通过（0 failed）

- [ ] **Step 2: E2E 验证 list-versions**

Run: `cd f:\OneDrive\Documents\TraeSOLO\AllEnricher\AllEnricher-v2 && python -m allenricher list-versions`
Expected: 显示基础数据版本 + 物种数据库版本

- [ ] **Step 3: E2E 验证 list-versions --lineage**

Run: `cd f:\OneDrive\Documents\TraeSOLO\AllEnricher\AllEnricher-v2 && python -m allenricher list-versions --lineage`
Expected: 显示构建血缘追踪报告

- [ ] **Step 4: E2E 验证 cleanup --dry-run**

Run: `cd f:\OneDrive\Documents\TraeSOLO\AllEnricher\AllEnricher-v2 && python -m allenricher cleanup --dry-run --keep 1`
Expected: 预览将删除的旧版本

- [ ] **Step 5: E2E 验证 check-update**

Run: `cd f:\OneDrive\Documents\TraeSOLO\AllEnricher\AllEnricher-v2 && python -m allenricher check-update`
Expected: 显示各数据源更新状态

- [ ] **Step 6: E2E 验证 analyze --use-version**

Run: `cd f:\OneDrive\Documents\TraeSOLO\AllEnricher\AllEnricher-v2 && python -m allenricher analyze -i "../AllEnricher-v1/example/example.glist" -s hsa -d GO --use-version v20260515 -o test_output/e2e_improvement`
Expected: 分析成功，TSV 文件头部包含版本注释

- [ ] **Step 7: 验证 TSV 版本注释**

Run: `head -10 test_output/e2e_improvement/GO_enrichment.tsv`
Expected: 前几行包含 `# AllEnricher version:`, `# Database version: v20260515`

- [ ] **Step 8: 运行版本管理专项测试**

Run: `cd f:\OneDrive\Documents\TraeSOLO\AllEnricher\AllEnricher-v2 && python -m pytest tests/test_version.py -v`
Expected: 全部通过（含新增测试）

---

## Task 5: 生成改进完成报告

**Files:**
- Create: `docs/version-management-improvement-report.md`

- [ ] **Step 1: 汇总所有改进项和测试结果，生成报告**

报告包含：
1. 改进项清单及完成状态
2. 代码变更摘要（文件、行号、变更内容）
3. 测试覆盖变化（改进前 vs 改进后）
4. E2E 测试结果
5. 已知遗留问题（如有）
