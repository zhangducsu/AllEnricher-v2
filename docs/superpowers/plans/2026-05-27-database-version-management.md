# 数据库版本管理方案 Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** 建立完整的数据库版本管理体系，支持远程更新检测、本地版本追溯、版本锁定和冗余清理

**Architecture:** 新增 `DatabaseVersionManager` 模块统一管理版本元数据；为每个数据源实现远程版本探测器（RemoteVersionChecker）；在 `download`/`build`/`analyze` 三个阶段分别集成版本感知逻辑；新增 `check-update` CLI 命令供用户主动检查更新。

**Tech Stack:** Python, requests (HTTP HEAD), json, dataclasses, pytest

---

## 现状问题总结

| 问题 | 现状 | 影响 |
|------|------|------|
| 版本号仅基于日期 | `GO20260527` 由下载当天日期决定，不反映数据源实际版本 | 同一天重下覆盖旧版，无法区分数据源是否真正更新 |
| 无版本元数据 | 构建产物不记录源数据版本号 | 无法追溯某次分析使用了哪个版本的源数据 |
| 无更新检测 | download 命令不检查远程是否有新版本 | 用户无法判断是否需要重新下载 |
| 无版本锁定 | analyze 始终使用最新版本 | 无法复现历史分析结果 |
| 无冗余清理 | 旧版本目录永久保留 | 磁盘空间持续增长 |
| KEGG/DO 无版本管理 | 直接存放文件，无版本目录 | 无法管理这两个数据源的版本 |
| auto_update 未实现 | config.py 中有此字段但无代码 | 配置项是空壳 |

## 远程数据源版本检测方案

| 数据源 | 版本获取方式 | 官方版本号 | 更新频率 |
|--------|-------------|-----------|---------|
| NCBI gene2go.gz | HTTP HEAD `Last-Modified` | 无 | 每日 |
| NCBI gene_info.gz | HTTP HEAD `Last-Modified` | 无 | 每日 |
| GO go-basic.obo | 文件内 `data-version:` 字段 | `releases/YYYY-MM-DD` | 每月 |
| EBI GOA proteomes | HTTP HEAD 目录页或具体文件 `Last-Modified` | 无 | 每月 |
| KEGG | REST API `GET info/kegg` | `Release N.M+/MM-DD` | 每日(工作日) |
| Reactome | 下载页面解析 `download.reactome.org/{N}/` | `v{N}` | 每季度 |
| NCBI Taxonomy | HTTP HEAD `Last-Modified` | 无 | 每日 |

---

## 文件结构

```
allenricher/database/
├── version.py              ← 新增：版本管理核心模块
│   ├── DatabaseVersion      数据类：单个数据源的版本信息
│   ├── VersionManifest      数据类：全局版本清单
│   ├── DatabaseVersionManager 类：版本清单读写、比较、清理
│   └── RemoteVersionChecker 类：远程版本检测
├── downloader.py           ← 修改：集成版本检测和记录
├── builder.py              ← 修改：构建时写入版本元数据
├── manager.py              ← 修改：支持版本锁定
├── goa_fetcher.py          ← 修改：无重大改动
└── kegg_fetcher.py         ← 修改：添加版本查询方法
allenricher/cli.py           ← 修改：新增 check-update 命令、download/build 版本参数
allenricher/core/config.py   ← 修改：新增版本锁定配置
tests/test_version.py        ← 新增：版本管理测试
```

---

## Task 1: 创建版本管理核心模块

**Files:**
- Create: `allenricher/database/version.py`
- Test: `tests/test_version.py`

- [ ] **Step 1: 定义版本数据类**

```python
# allenricher/database/version.py
"""数据库版本管理模块

管理 AllEnricher 各数据源的版本信息，支持：
- 远程版本检测（检查数据源是否有更新）
- 本地版本清单读写（versions.json）
- 版本比较和更新判断
- 旧版本清理
"""
from __future__ import annotations

import json
import logging
from dataclasses import dataclass, field, asdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List, Optional

logger = logging.getLogger(__name__)


@dataclass
class DatabaseVersion:
    """单个数据源的版本信息"""

    source: str                          # 数据源标识，如 "go", "kegg", "reactome", "taxonomy"
    remote_version: Optional[str] = None # 远程数据源版本号（如 "releases/2026-03-25", "Release 118.0+", "v96"）
    remote_last_modified: Optional[str] = None  # 远程文件最后修改时间（HTTP Last-Modified）
    local_version: Optional[str] = None  # 本地版本目录名（如 "GO20260527"）
    local_path: Optional[str] = None     # 本地文件/目录路径（相对于 database/）
    downloaded_at: Optional[str] = None  # 下载时间（ISO 8601）
    file_hash: Optional[str] = None      # 文件 MD5/SHA256（可选，用于完整性校验）

    def is_newer_than(self, other: "DatabaseVersion") -> bool:
        """判断远程版本是否比本地版本更新"""
        if not self.remote_last_modified or not other.remote_last_modified:
            return False
        try:
            remote_dt = datetime.strptime(
                self.remote_last_modified, "%a, %d %b %Y %H:%M:%S GMT"
            )
            local_dt = datetime.strptime(
                other.remote_last_modified, "%a, %d %b %Y %H:%M:%S GMT"
            )
            return remote_dt > local_dt
        except ValueError:
            return False

    def to_dict(self) -> dict:
        return asdict(self)

    @classmethod
    def from_dict(cls, data: dict) -> "DatabaseVersion":
        return cls(**{k: v for k, v in data.items() if k in cls.__dataclass_fields__})


@dataclass
class VersionManifest:
    """全局版本清单"""

    created_at: str = field(default_factory=lambda: datetime.now(timezone.utc).isoformat())
    updated_at: str = field(default_factory=lambda: datetime.now(timezone.utc).isoformat())
    versions: Dict[str, DatabaseVersion] = field(default_factory=dict)

    def get(self, source: str) -> Optional[DatabaseVersion]:
        return self.versions.get(source)

    def set(self, source: str, version: DatabaseVersion) -> None:
        self.versions[source] = version
        self.updated_at = datetime.now(timezone.utc).isoformat()

    def to_dict(self) -> dict:
        return {
            "created_at": self.created_at,
            "updated_at": self.updated_at,
            "versions": {k: v.to_dict() for k, v in self.versions.items()},
        }

    @classmethod
    def from_dict(cls, data: dict) -> "VersionManifest":
        manifest = cls(
            created_at=data.get("created_at", ""),
            updated_at=data.get("updated_at", ""),
        )
        for source, ver_data in data.get("versions", {}).items():
            manifest.versions[source] = DatabaseVersion.from_dict(ver_data)
        return manifest

    def save(self, path: Path) -> None:
        path.parent.mkdir(parents=True, exist_ok=True)
        with open(path, "w", encoding="utf-8") as f:
            json.dump(self.to_dict(), f, indent=2, ensure_ascii=False)
        logger.info("版本清单已保存: %s", path)

    @classmethod
    def load(cls, path: Path) -> "VersionManifest":
        if not path.exists():
            logger.info("版本清单不存在，创建新清单: %s", path)
            return cls()
        with open(path, "r", encoding="utf-8") as f:
            return cls.from_dict(json.load(f))
```

- [ ] **Step 2: 定义 DatabaseVersionManager**

```python
class DatabaseVersionManager:
    """数据库版本管理器

    负责：
    - 读写 versions.json 版本清单
    - 查询本地已安装的版本
    - 比较本地与远程版本差异
    """

    MANIFEST_FILENAME = "versions.json"

    def __init__(self, database_dir: str = "./database"):
        self.database_dir = Path(database_dir)
        self.manifest_path = self.database_dir / self.MANIFEST_FILENAME
        self._manifest: Optional[VersionManifest] = None

    @property
    def manifest(self) -> VersionManifest:
        if self._manifest is None:
            self._manifest = VersionManifest.load(self.manifest_path)
        return self._manifest

    def save_manifest(self) -> None:
        self.manifest.save(self.manifest_path)

    def record_download(
        self,
        source: str,
        local_version: str,
        local_path: str,
        remote_version: Optional[str] = None,
        remote_last_modified: Optional[str] = None,
    ) -> None:
        """记录一次下载操作"""
        ver = DatabaseVersion(
            source=source,
            remote_version=remote_version,
            remote_last_modified=remote_last_modified,
            local_version=local_version,
            local_path=local_path,
            downloaded_at=datetime.now(timezone.utc).isoformat(),
        )
        self.manifest.set(source, ver)
        self.save_manifest()

    def get_local_version(self, source: str) -> Optional[DatabaseVersion]:
        """获取某个数据源的本地版本信息"""
        return self.manifest.get(source)

    def list_local_versions(self) -> Dict[str, DatabaseVersion]:
        """列出所有本地已安装的数据源版本"""
        return dict(self.manifest.versions)

    def list_installed_basic_versions(self) -> Dict[str, List[str]]:
        """扫描 database/basic/ 目录，列出已安装的基础数据版本

        Returns:
            {source_name: [version_dir_names]}
            例如 {"go": ["GO20260515", "GO20260527"], "reactome": ["Reactome20260515"]}
        """
        result: Dict[str, List[str]] = {}
        basic_dir = self.database_dir / "basic"

        # GO
        go_dir = basic_dir / "go"
        if go_dir.exists():
            result["go"] = sorted(
                d.name for d in go_dir.iterdir() if d.is_dir() and d.name.startswith("GO")
            )

        # Reactome
        reactome_dir = basic_dir / "reactome"
        if reactome_dir.exists():
            result["reactome"] = sorted(
                d.name for d in reactome_dir.iterdir()
                if d.is_dir() and d.name.startswith("Reactome")
            )

        # KEGG（无版本目录，检查是否有文件）
        kegg_dir = basic_dir / "kegg"
        if kegg_dir.exists() and any(kegg_dir.iterdir()):
            result["kegg"] = ["cached"]

        # DO
        do_dir = basic_dir / "do"
        if do_dir.exists() and any(do_dir.iterdir()):
            result["do"] = ["cached"]

        # Taxonomy
        tax_dir = basic_dir / "taxonomy"
        if tax_dir.exists() and (tax_dir / "names.dmp").exists():
            result["taxonomy"] = ["cached"]

        return result

    def list_installed_organism_versions(self) -> List[str]:
        """扫描 database/organism/ 目录，列出已构建的物种数据库版本

        Returns:
            版本目录名列表，如 ["v20260515"]
        """
        organism_dir = self.database_dir / "organism"
        if not organism_dir.exists():
            return []
        return sorted(
            d.name for d in organism_dir.iterdir()
            if d.is_dir() and d.name.startswith("v")
        )

    def get_organism_build_info(self, version: str) -> Dict[str, List[str]]:
        """获取某个构建版本下包含哪些物种

        Returns:
            {version: [species_codes]}
        """
        version_dir = self.database_dir / "organism" / version
        if not version_dir.exists():
            return {}
        return {
            version: sorted(
                d.name for d in version_dir.iterdir() if d.is_dir()
            )
        }

    def find_stale_versions(self, keep_count: int = 2) -> Dict[str, List[str]]:
        """找出可以清理的旧版本

        Args:
            keep_count: 每个数据源保留的最新版本数量

        Returns:
            {source: [stale_version_names]}
        """
        stale: Dict[str, List[str]] = {}

        # 基础数据
        basic_versions = self.list_installed_basic_versions()
        for source, versions in basic_versions.items():
            if len(versions) > keep_count:
                stale[source] = versions[:-keep_count]

        # 物种数据库
        organism_versions = self.list_installed_organism_versions()
        if len(organism_versions) > keep_count:
            stale["organism"] = organism_versions[:-keep_count]

        return stale

    def remove_stale_versions(self, keep_count: int = 2, dry_run: bool = True) -> Dict[str, List[str]]:
        """清理旧版本

        Args:
            keep_count: 保留的最新版本数量
            dry_run: 仅预览不实际删除

        Returns:
            {source: [removed_version_names]}
        """
        stale = self.find_stale_versions(keep_count)
        removed: Dict[str, List[str]] = {}

        for source, versions in stale.items():
            removed[source] = []
            for ver in versions:
                if source == "organism":
                    dir_path = self.database_dir / "organism" / ver
                else:
                    dir_path = self.database_dir / "basic" / source / ver

                if dir_path.exists():
                    if not dry_run:
                        import shutil
                        shutil.rmtree(dir_path)
                        logger.info("已删除: %s", dir_path)
                    else:
                        logger.info("[dry-run] 将删除: %s", dir_path)
                    removed[source].append(ver)

        return removed
```

- [ ] **Step 3: 编写测试**

```python
# tests/test_version.py
import json
import pytest
from pathlib import Path
from allenricher.database.version import DatabaseVersion, VersionManifest, DatabaseVersionManager


class TestDatabaseVersion:
    def test_from_dict_roundtrip(self):
        data = {
            "source": "go",
            "remote_version": "releases/2026-03-25",
            "remote_last_modified": "Thu, 02 Apr 2026 02:53:36 GMT",
            "local_version": "GO20260527",
            "local_path": "basic/go/GO20260527",
            "downloaded_at": "2026-05-27T00:00:00+00:00",
        }
        ver = DatabaseVersion.from_dict(data)
        assert ver.source == "go"
        assert ver.remote_version == "releases/2026-03-25"
        assert ver.to_dict() == data

    def test_is_newer_than(self):
        old = DatabaseVersion(
            source="go",
            remote_last_modified="Thu, 02 Apr 2026 02:53:36 GMT",
        )
        new = DatabaseVersion(
            source="go",
            remote_last_modified="Thu, 02 May 2026 02:53:36 GMT",
        )
        assert new.is_newer_than(old) is True
        assert old.is_newer_than(new) is False

    def test_is_newer_than_missing_dates(self):
        a = DatabaseVersion(source="go")
        b = DatabaseVersion(source="go", remote_last_modified="Thu, 02 May 2026 02:53:36 GMT")
        assert a.is_newer_than(b) is False
        assert b.is_newer_than(a) is False


class TestVersionManifest:
    def test_save_and_load(self, tmp_path):
        manifest = VersionManifest()
        manifest.set("go", DatabaseVersion(source="go", local_version="GO20260527"))
        manifest.set("kegg", DatabaseVersion(source="kegg", local_version="cached"))

        path = tmp_path / "versions.json"
        manifest.save(path)

        loaded = VersionManifest.load(path)
        assert loaded.get("go").local_version == "GO20260527"
        assert loaded.get("kegg").local_version == "cached"

    def test_load_nonexistent(self, tmp_path):
        manifest = VersionManifest.load(tmp_path / "missing.json")
        assert len(manifest.versions) == 0


class TestDatabaseVersionManager:
    def test_record_and_retrieve(self, tmp_path):
        mgr = DatabaseVersionManager(str(tmp_path))
        mgr.record_download(
            source="go",
            local_version="GO20260527",
            local_path="basic/go/GO20260527",
            remote_version="releases/2026-03-25",
            remote_last_modified="Thu, 02 Apr 2026 02:53:36 GMT",
        )
        ver = mgr.get_local_version("go")
        assert ver is not None
        assert ver.local_version == "GO20260527"
        assert ver.remote_version == "releases/2026-03-25"

    def test_manifest_persists(self, tmp_path):
        mgr = DatabaseVersionManager(str(tmp_path))
        mgr.record_download("go", "GO20260527", "basic/go/GO20260527")
        mgr2 = DatabaseVersionManager(str(tmp_path))
        assert mgr2.get_local_version("go").local_version == "GO20260527"
```

- [ ] **Step 4: 运行测试验证**

Run: `cd f:\OneDrive\Documents\TraeSOLO\AllEnricher\AllEnricher-v2 && python -m pytest tests/test_version.py -v`
Expected: 7 passed

- [ ] **Step 5: Commit**

```bash
git add allenricher/database/version.py tests/test_version.py
git commit -m "feat(version): 新增数据库版本管理核心模块"
```

---

## Task 2: 实现远程版本检测器

**Files:**
- Modify: `allenricher/database/version.py` (追加 RemoteVersionChecker)
- Modify: `tests/test_version.py` (追加测试)

- [ ] **Step 1: 实现 RemoteVersionChecker**

在 `version.py` 末尾追加：

```python
class RemoteVersionChecker:
    """远程数据源版本检测器

    通过 HTTP HEAD / API 查询检测远程数据源是否有更新。
    """

    TIMEOUT = 30  # 秒

    # 各数据源的关键文件 URL
    SOURCE_URLS = {
        "gene2go": "https://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2go.gz",
        "gene_info": "https://ftp.ncbi.nlm.nih.gov/gene/DATA/gene_info.gz",
        "go_obo": "http://purl.obolibrary.org/obo/go/go-basic.obo",
        "goa_proteomes": "https://ftp.ebi.ac.uk/pub/databases/GO/goa/proteomes/",
        "kegg": "https://rest.kegg.jp/info/kegg",
        "reactome": "https://reactome.org/download/",
        "taxonomy": "https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz",
    }

    def __init__(self, timeout: int = TIMEOUT):
        self.timeout = timeout

    def check_head(self, url: str) -> Optional[Dict[str, str]]:
        """发送 HTTP HEAD 请求获取 Last-Modified / ETag

        Returns:
            {"last_modified": "...", "etag": "...", "content_length": "..."} 或 None
        """
        import requests
        try:
            resp = requests.head(url, timeout=self.timeout, allow_redirects=True)
            resp.raise_for_status()
            result = {}
            if "Last-Modified" in resp.headers:
                result["last_modified"] = resp.headers["Last-Modified"]
            if "ETag" in resp.headers:
                result["etag"] = resp.headers["ETag"]
            if "Content-Length" in resp.headers:
                result["content_length"] = resp.headers["Content-Length"]
            return result if result else None
        except Exception as e:
            logger.warning("HEAD 请求失败 %s: %s", url, e)
            return None

    def check_go_obo_version(self) -> Optional[Dict[str, str]]:
        """检测 GO Ontology 版本

        下载 go-basic.obo 前 10 行，解析 data-version 字段。
        """
        import requests
        url = self.SOURCE_URLS["go_obo"]
        try:
            # 使用 Range 请求只获取前 1024 字节
            resp = requests.get(
                url, timeout=self.timeout, headers={"Range": "bytes=0-1023"}, stream=True
            )
            resp.raise_for_status()
            first_chunk = resp.content.decode("utf-8", errors="ignore")
            for line in first_chunk.split("\n"):
                if line.startswith("data-version:"):
                    version = line.split(":", 1)[1].strip()
                    return {
                        "remote_version": version,
                        "last_modified": resp.headers.get("Last-Modified", ""),
                    }
        except Exception as e:
            logger.warning("检测 GO 版本失败: %s", e)
        return None

    def check_kegg_version(self) -> Optional[Dict[str, str]]:
        """检测 KEGG 版本

        调用 REST API info/kegg 获取版本号。
        """
        import requests
        url = self.SOURCE_URLS["kegg"]
        try:
            resp = requests.get(url, timeout=self.timeout)
            resp.raise_for_status()
            for line in resp.text.strip().split("\n"):
                if line.startswith("Release"):
                    return {
                        "remote_version": line.strip(),
                        "last_modified": "",
                    }
        except Exception as e:
            logger.warning("检测 KEGG 版本失败: %s", e)
        return None

    def check_reactome_version(self) -> Optional[Dict[str, str]]:
        """检测 Reactome 版本

        解析下载页面中的版本化 URL。
        """
        import requests
        import re
        url = self.SOURCE_URLS["reactome"]
        try:
            resp = requests.get(url, timeout=self.timeout)
            resp.raise_for_status()
            # 匹配 download.reactome.org/{N}/ 或 /{N}/
            match = re.search(r"download\.reactome\.org/(\d+)/", resp.text)
            if match:
                version_num = match.group(1)
                return {
                    "remote_version": f"v{version_num}",
                    "last_modified": resp.headers.get("Last-Modified", ""),
                }
        except Exception as e:
            logger.warning("检测 Reactome 版本失败: %s", e)
        return None

    def check_all_sources(self) -> Dict[str, Dict[str, str]]:
        """检测所有数据源的远程版本

        Returns:
            {
                "gene2go": {"last_modified": "...", ...},
                "go_obo": {"remote_version": "releases/2026-03-25", ...},
                "kegg": {"remote_version": "Release 118.0+...", ...},
                ...
            }
        """
        results: Dict[str, Dict[str, str]] = {}

        # HTTP HEAD 类（gene2go, gene_info, taxonomy）
        for source in ["gene2go", "gene_info", "taxonomy"]:
            url_key = source
            info = self.check_head(self.SOURCE_URLS[url_key])
            if info:
                results[source] = info

        # GO OBO（文件内容解析）
        go_info = self.check_go_obo_version()
        if go_info:
            results["go_obo"] = go_info

        # GOA Proteomes（目录 HEAD）
        goa_info = self.check_head(self.SOURCE_URLS["goa_proteomes"])
        if goa_info:
            results["goa_proteomes"] = goa_info

        # KEGG（API 查询）
        kegg_info = self.check_kegg_version()
        if kegg_info:
            results["kegg"] = kegg_info

        # Reactome（页面解析）
        reactome_info = self.check_reactome_version()
        if reactome_info:
            results["reactome"] = reactome_info

        return results

    def check_updates(self, local_manager: DatabaseVersionManager) -> Dict[str, Dict]:
        """检查所有数据源是否有更新

        比较远程版本与本地记录，返回更新状态。

        Returns:
            {
                "gene2go": {
                    "has_update": True,
                    "local": {"last_modified": "...", ...},
                    "remote": {"last_modified": "...", ...},
                },
                ...
            }
        """
        remote_versions = self.check_all_sources()
        update_status: Dict[str, Dict] = {}

        for source, remote_info in remote_versions.items():
            local_ver = local_manager.get_local_version(source)
            has_update = False

            if local_ver is None:
                has_update = True  # 本地从未下载过
            elif "last_modified" in remote_info and local_ver.remote_last_modified:
                # 比较 Last-Modified
                try:
                    from email.utils import parsedate_to_datetime
                    remote_dt = parsedate_to_datetime(remote_info["last_modified"])
                    local_dt = parsedate_to_datetime(local_ver.remote_last_modified)
                    has_update = remote_dt > local_dt
                except Exception:
                    has_update = False
            elif "remote_version" in remote_info and local_ver.remote_version:
                has_update = remote_info["remote_version"] != local_ver.remote_version

            update_status[source] = {
                "has_update": has_update,
                "local": {
                    "version": local_ver.local_version if local_ver else None,
                    "remote_version": local_ver.remote_version if local_ver else None,
                    "last_modified": local_ver.remote_last_modified if local_ver else None,
                    "downloaded_at": local_ver.downloaded_at if local_ver else None,
                },
                "remote": remote_info,
            }

        return update_status
```

- [ ] **Step 2: 追加测试**

```python
# 追加到 tests/test_version.py
from allenricher.database.version import RemoteVersionChecker


class TestRemoteVersionChecker:
    def test_check_kegg_version(self):
        checker = RemoteVersionChecker()
        result = checker.check_kegg_version()
        assert result is not None
        assert "remote_version" in result
        assert "Release" in result["remote_version"]

    def test_check_go_obo_version(self):
        checker = RemoteVersionChecker()
        result = checker.check_go_obo_version()
        assert result is not None
        assert "remote_version" in result
        assert "releases/" in result["remote_version"]

    def test_check_all_sources(self):
        checker = RemoteVersionChecker()
        results = checker.check_all_sources()
        # 至少应检测到 gene2go, go_obo, kegg
        assert len(results) >= 3

    def test_check_updates_no_local(self, tmp_path):
        checker = RemoteVersionChecker()
        mgr = DatabaseVersionManager(str(tmp_path))
        status = checker.check_updates(mgr)
        # 本地无任何记录，所有数据源都应显示有更新
        for source, info in status.items():
            assert info["has_update"] is True
```

- [ ] **Step 3: 运行测试**

Run: `cd f:\OneDrive\Documents\TraeSOLO\AllEnricher\AllEnricher-v2 && python -m pytest tests/test_version.py -v`
Expected: 全部通过（注意网络测试可能超时，可标记为 slow）

- [ ] **Step 4: Commit**

```bash
git add allenricher/database/version.py tests/test_version.py
git commit -m "feat(version): 实现远程版本检测器"
```

---

## Task 3: 集成版本记录到下载流程

**Files:**
- Modify: `allenricher/database/downloader.py`
- Modify: `allenricher/database/version.py`

- [ ] **Step 1: 在 download_go_basic() 中记录版本**

在 `download_go_basic()` 方法末尾（所有文件下载完成后），添加版本记录：

```python
# 在 download_go_basic() 方法末尾追加
from allenricher.database.version import DatabaseVersionManager

# 记录 gene2go.gz 的版本
vm = DatabaseVersionManager(root_dir=str(self.root_dir))
gene2go_info = RemoteVersionChecker().check_head(
    "https://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2go.gz"
)
if gene2go_info:
    vm.record_download(
        source="gene2go",
        local_version=version,
        local_path=f"basic/go/{version}",
        remote_last_modified=gene2go_info.get("last_modified"),
    )

# 记录 go-basic.obo 的版本
go_info = RemoteVersionChecker().check_go_obo_version()
if go_info:
    vm.record_download(
        source="go_obo",
        local_version=version,
        local_path=f"basic/go/{version}/go-basic.obo",
        remote_version=go_info.get("remote_version"),
        remote_last_modified=go_info.get("last_modified"),
    )
```

- [ ] **Step 2: 在 download_reactome_basic() 中记录版本**

类似地，在 Reactome 下载完成后记录版本。

- [ ] **Step 3: 在 _download_taxonomy_names() 中记录版本**

类似地，在 Taxonomy 下载完成后记录版本。

- [ ] **Step 4: 验证下载后生成 versions.json**

Run: `cd f:\OneDrive\Documents\TraeSOLO\AllEnricher\AllEnricher-v2 && python -c "from allenricher.database.version import DatabaseVersionManager; vm = DatabaseVersionManager(); print(vm.manifest.to_dict())" | python -m json.tool`
Expected: 输出包含 gene2go, go_obo, taxonomy 等数据源的版本信息

- [ ] **Step 5: Commit**

```bash
git add allenricher/database/downloader.py
git commit -m "feat(download): 下载时自动记录版本元数据"
```

---

## Task 4: 集成版本记录到构建流程

**Files:**
- Modify: `allenricher/database/builder.py`

- [ ] **Step 1: 在构建输出目录中写入 build_manifest.json**

在 `build_species_db()` 方法末尾，写入构建清单：

```python
# 在 build_species_db() 方法末尾追加
import json
from datetime import datetime, timezone

build_manifest = {
    "built_at": datetime.now(timezone.utc).isoformat(),
    "species": species,
    "taxid": taxid,
    "databases": databases,
    "source_versions": {},
}

# 记录各数据源版本
downloader = DataDownloader(root_dir=str(self.root_dir))
vm = DatabaseVersionManager(root_dir=str(self.root_dir))

if "GO" in databases:
    go_ver = downloader.get_latest_go_version()
    build_manifest["source_versions"]["go"] = go_ver
    go_info = vm.get_local_version("go_obo")
    if go_info:
        build_manifest["source_versions"]["go_remote"] = go_info.remote_version

if "KEGG" in databases:
    kegg_info = vm.get_local_version("kegg")
    if kegg_info:
        build_manifest["source_versions"]["kegg"] = kegg_info.remote_version

# 写入构建清单
manifest_path = outdir / "build_manifest.json"
with open(manifest_path, "w", encoding="utf-8") as f:
    json.dump(build_manifest, f, indent=2, ensure_ascii=False)
logger.info("构建清单已保存: %s", manifest_path)
```

- [ ] **Step 2: 验证构建后生成 build_manifest.json**

Run: 检查 `database/organism/v{date}/hsa/build_manifest.json` 是否存在

- [ ] **Step 3: Commit**

```bash
git add allenricher/database/builder.py
git commit -m "feat(build): 构建时写入版本清单 build_manifest.json"
```

---

## Task 5: 支持版本锁定

**Files:**
- Modify: `allenricher/database/manager.py`
- Modify: `allenricher/core/config.py`

- [ ] **Step 1: 在 config.py 中添加版本锁定配置**

```python
# 在 EnrichmentConfig 数据类中添加
database_version: Optional[str] = None  # 锁定使用的数据库版本（如 "v20260515"），None 表示使用最新
```

- [ ] **Step 2: 修改 DatabaseManager 支持版本锁定**

修改 `_find_species_dir()` 方法，添加 `version` 参数：

```python
def _find_species_dir(
    self, database_dir: Path, species: str, version: Optional[str] = None
) -> Path:
    """自动查找物种数据库目录

    Args:
        database_dir: 数据库根目录
        species: 物种代码
        version: 锁定版本（如 "v20260515"），None 表示使用最新

    Returns:
        物种数据库目录路径
    """
    organism_dir = database_dir / "organism"

    # 如果指定了版本，直接使用
    if version and organism_dir.exists():
        species_dir = organism_dir / version / species
        if species_dir.exists():
            return species_dir
        logger.warning("指定版本 %s 的物种 %s 不存在，回退到最新版本", version, species)

    # 自动查找最新版本
    if organism_dir.exists():
        for version_dir in sorted(organism_dir.iterdir(), reverse=True):
            if version_dir.is_dir():
                species_dir = version_dir / species
                if species_dir.exists():
                    return species_dir

    # v1 兼容
    if (database_dir / f"{species}.GO2gene.tab.gz").exists():
        return database_dir

    return database_dir
```

- [ ] **Step 3: 修改 CLI analyze 命令传递版本参数**

在 `cmd_analyze()` 中：

```python
# 在创建 DatabaseManager 时传递版本
db_manager = DatabaseManager(db_dir, config.species)
# load_databases 时传递版本
db_manager.load_databases(config.databases, version=config.database_version)
```

- [ ] **Step 4: 在 analyze 输出中显示使用的数据库版本**

在分析报告的元数据中添加 `database_version` 字段。

- [ ] **Step 5: Commit**

```bash
git add allenricher/database/manager.py allenricher/core/config.py allenricher/cli.py
git commit -m "feat(version): 支持版本锁定，analyze 可指定数据库版本"
```

---

## Task 6: 新增 check-update CLI 命令

**Files:**
- Modify: `allenricher/cli.py`

- [ ] **Step 1: 添加 check-update 子命令**

在 CLI 参数解析中添加：

```python
# 在 cli.py 的 subparsers 中添加
sub_check = subparsers.add_parser(
    "check-update",
    help="检查远程数据源是否有更新",
)
sub_check.add_argument(
    "--database-dir",
    type=str,
    default=None,
    help="数据库目录（默认 ./database）",
)
sub_check.add_argument(
    "--json",
    action="store_true",
    help="以 JSON 格式输出",
)
```

- [ ] **Step 2: 实现 cmd_check_update() 处理函数**

```python
def cmd_check_update(args):
    """检查远程数据源是否有更新"""
    from allenricher.database.version import RemoteVersionChecker, DatabaseVersionManager

    db_dir = args.database_dir or "./database"
    checker = RemoteVersionChecker()
    vm = DatabaseVersionManager(db_dir)

    print("正在检查远程数据源更新...")
    print("=" * 60)

    status = checker.check_updates(vm)

    has_any_update = False
    for source, info in sorted(status.items()):
        local_ver = info["local"].get("version", "未安装")
        local_remote = info["local"].get("remote_version", "-")
        local_date = info["local"].get("last_modified", "-")
        remote_ver = info["remote"].get("remote_version", "-")
        remote_date = info["remote"].get("last_modified", "-")

        if info["has_update"]:
            status_icon = "🔄"
            has_any_update = True
        else:
            status_icon = "✅"

        print(f"{status_icon} {source:<20} 本地: {local_ver:<20} 远程: {remote_ver or remote_date}")

    print("=" * 60)
    if has_any_update:
        print("有可用的更新。运行 `allenricher download -d go,kegg,reactome` 下载最新数据。")
    else:
        print("所有数据源均为最新版本。")

    if args.json:
        import json
        print(json.dumps(status, indent=2, ensure_ascii=False))
```

- [ ] **Step 3: 验证 check-update 命令**

Run: `cd f:\OneDrive\Documents\TraeSOLO\AllEnricher\AllEnricher-v2 && python -m allenricher check-update`
Expected: 显示各数据源的更新状态

- [ ] **Step 4: Commit**

```bash
git add allenricher/cli.py
git commit -m "feat(cli): 新增 check-update 命令检查远程数据源更新"
```

---

## Task 7: 新增 cleanup CLI 命令

**Files:**
- Modify: `allenricher/cli.py`

- [ ] **Step 1: 添加 cleanup 子命令**

```python
sub_cleanup = subparsers.add_parser(
    "cleanup",
    help="清理旧版本的数据库文件",
)
sub_cleanup.add_argument(
    "--keep",
    type=int,
    default=2,
    help="保留的最新版本数量（默认 2）",
)
sub_cleanup.add_argument(
    "--dry-run",
    action="store_true",
    help="仅预览不实际删除",
)
sub_cleanup.add_argument(
    "--database-dir",
    type=str,
    default=None,
    help="数据库目录（默认 ./database）",
)
```

- [ ] **Step 2: 实现 cmd_cleanup() 处理函数**

```python
def cmd_cleanup(args):
    """清理旧版本的数据库文件"""
    from allenricher.database.version import DatabaseVersionManager

    db_dir = args.database_dir or "./database"
    vm = DatabaseVersionManager(db_dir)

    if args.dry_run:
        print("[预览模式] 以下旧版本将被删除：")
    else:
        print("正在清理旧版本...")

    removed = vm.remove_stale_versions(keep_count=args.keep, dry_run=args.dry_run)

    total = sum(len(v) for v in removed.values())
    if total == 0:
        print("没有需要清理的旧版本。")
    else:
        for source, versions in removed.items():
            for ver in versions:
                print(f"  删除: {source}/{ver}")
        print(f"共清理 {total} 个旧版本目录。")
```

- [ ] **Step 3: Commit**

```bash
git add allenricher/cli.py
git commit -m "feat(cli): 新增 cleanup 命令清理旧版本数据库"
```

---

## Task 8: 修改 download 命令集成更新检测

**Files:**
- Modify: `allenricher/cli.py`

- [ ] **Step 1: download 命令默认先检查更新**

修改 `cmd_download()` 函数，在下载前自动检查是否有更新：

```python
def cmd_download(args):
    # ... 现有参数解析 ...

    # 自动检查更新（除非用户指定 --force）
    if not args.force:
        from allenricher.database.version import RemoteVersionChecker, DatabaseVersionManager
        checker = RemoteVersionChecker()
        vm = DatabaseVersionManager(db_dir)
        status = checker.check_updates(vm)

        requested = set(d.strip().lower() for d in args.databases.split(","))
        relevant_updates = {
            s: info for s, info in status.items()
            if info["has_update"] and _source_matches_database(s, requested)
        }

        if relevant_updates:
            print("检测到以下数据源有更新：")
            for source, info in relevant_updates.items():
                print(f"  🔄 {source}")
            print("开始下载...")
        else:
            print("所有请求的数据源均为最新版本。使用 --force 强制重新下载。")
            return
```

- [ ] **Step 2: 添加 --force 参数到 download 命令**

```python
sub_download.add_argument(
    "--force",
    action="store_true",
    help="强制重新下载，即使本地已是最新版本",
)
```

- [ ] **Step 3: Commit**

```bash
git add allenricher/cli.py
git commit -m "feat(download): 下载前自动检查更新，新增 --force 参数"
```

---

## Task 9: 端到端验证

**Files:**
- Test: 手动 E2E 测试

- [ ] **Step 1: 验证 check-update 命令**

```bash
cd f:\OneDrive\Documents\TraeSOLO\AllEnricher\AllEnricher-v2
python -m allenricher check-update
```

Expected: 显示各数据源的更新状态

- [ ] **Step 2: 验证 download 后生成 versions.json**

```bash
python -c "from pathlib import Path; import json; p = Path('./database/versions.json'); print(json.dumps(json.load(p.open()), indent=2))" | head -30
```

Expected: 包含各数据源的版本信息

- [ ] **Step 3: 验证 build 后生成 build_manifest.json**

```bash
python -c "from pathlib import Path; import json; p = list(Path('./database/organism').glob('*/hsa/build_manifest.json')); print(json.dumps(json.load(p[0].open()), indent=2) if p else 'Not found')"
```

Expected: 包含构建时使用的源数据版本

- [ ] **Step 4: 验证 cleanup 命令**

```bash
python -m allenricher cleanup --dry-run --keep 1
```

Expected: 预览将删除的旧版本目录

- [ ] **Step 5: 运行完整测试套件**

```bash
python -m pytest tests/test_version.py -v
python -m pytest tests/ -v --tb=short 2>&1 | tail -20
```

Expected: 全部通过

---

## Task 10: 新增 list-versions CLI 命令（本地版本清单查看）

**Files:**
- Modify: `allenricher/cli.py`
- Modify: `allenricher/database/version.py`

- [ ] **Step 1: 在 version.py 中添加版本清单格式化方法**

在 `DatabaseVersionManager` 类中追加：

```python
def get_summary_table(self) -> str:
    """生成人类可读的本地版本清单表格

    Returns:
        格式化的表格字符串
    """
    lines = []
    lines.append("本地数据库版本清单")
    lines.append("=" * 80)

    # 基础数据版本
    basic_versions = self.list_installed_basic_versions()
    if basic_versions:
        lines.append("\n[基础数据 (basic/)]")
        lines.append(f"  {'数据源':<15} {'已安装版本':<40}")
        lines.append(f"  {'-'*15} {'-'*40}")
        for source, versions in basic_versions.items():
            ver_str = ", ".join(versions) if versions else "无"
            # 标记最新版本
            latest = f"← 最新" if versions else ""
            lines.append(f"  {source:<15} {ver_str:<40} {latest}")

    # 物种数据库版本
    organism_versions = self.list_installed_organism_versions()
    if organism_versions:
        lines.append("\n[物种数据库 (organism/)]")
        lines.append(f"  {'版本':<15} {'包含物种'}")
        lines.append(f"  {'-'*15} {'-'*40}")
        for ver in organism_versions:
            species_list = self.get_organism_build_info(ver).get(ver, [])
            species_str = ", ".join(species_list) if species_list else "空"
            latest = "← 最新" if ver == organism_versions[-1] else ""
            lines.append(f"  {ver:<15} {species_str:<40} {latest}")

    # versions.json 中记录的远程版本
    local_records = self.list_local_versions()
    if local_records:
        lines.append("\n[版本元数据 (versions.json)]")
        lines.append(f"  {'数据源':<20} {'本地版本':<20} {'远程版本':<25} {'下载时间'}")
        lines.append(f"  {'-'*20} {'-'*20} {'-'*25} {'-'*20}")
        for source, ver in sorted(local_records.items()):
            remote_ver = ver.remote_version or "-"
            downloaded = (ver.downloaded_at[:10] if ver.downloaded_at else "-")
            lines.append(f"  {source:<20} {ver.local_version or '-':<20} {remote_ver:<25} {downloaded}")

    lines.append("=" * 80)
    return "\n".join(lines)

def get_summary_json(self) -> dict:
    """生成 JSON 格式的版本清单"""
    return {
        "basic_versions": self.list_installed_basic_versions(),
        "organism_versions": self.list_installed_organism_versions(),
        "version_records": {
            source: ver.to_dict()
            for source, ver in self.list_local_versions().items()
        },
    }
```

- [ ] **Step 2: 在 cli.py 中添加 list-versions 子命令**

```python
# 在 cli.py 的 subparsers 中添加
sub_list_ver = subparsers.add_parser(
    "list-versions",
    help="查看本地已安装的数据库版本",
)
sub_list_ver.add_argument(
    "--database-dir",
    type=str,
    default=None,
    help="数据库目录（默认 ./database）",
)
sub_list_ver.add_argument(
    "--json",
    action="store_true",
    help="以 JSON 格式输出",
)
```

- [ ] **Step 3: 实现 cmd_list_versions() 处理函数**

```python
def cmd_list_versions(args):
    """查看本地已安装的数据库版本"""
    from allenricher.database.version import DatabaseVersionManager

    db_dir = args.database_dir or "./database"
    vm = DatabaseVersionManager(db_dir)

    if args.json:
        import json
        print(json.dumps(vm.get_summary_json(), indent=2, ensure_ascii=False))
    else:
        print(vm.get_summary_table())
```

- [ ] **Step 4: 验证**

Run: `cd f:\OneDrive\Documents\TraeSOLO\AllEnricher\AllEnricher-v2 && python -m allenricher list-versions`
Expected: 显示基础数据版本、物种数据库版本、版本元数据三段表格

- [ ] **Step 5: Commit**

```bash
git add allenricher/database/version.py allenricher/cli.py
git commit -m "feat(cli): 新增 list-versions 命令查看本地数据库版本"
```

---

## Task 11: 构建血缘追踪（build_manifest 增强）

**Files:**
- Modify: `allenricher/database/builder.py`
- Modify: `allenricher/database/version.py`

**目标**: 构建产物中记录完整的源数据依赖链，支持从 organism 版本反查依赖的 basic 版本。

- [ ] **Step 1: 增强 build_manifest.json 内容**

修改 Task 4 中定义的 `build_manifest.json` 结构，增加完整的依赖追踪：

```python
# 在 build_species_db() 方法末尾，写入增强版构建清单
import json
from datetime import datetime, timezone
from allenricher.database.version import DatabaseVersionManager

vm = DatabaseVersionManager(root_dir=str(self.root_dir))

build_manifest = {
    "schema_version": "1.0",
    "built_at": datetime.now(timezone.utc).isoformat(),
    "allenricher_version": __import__("allenricher").__version__,
    "species": species,
    "taxid": taxid,
    "databases": sorted(databases),
    "dependencies": {},   # 关键：记录每个数据库依赖的 basic 源版本
    "source_versions": {}, # 记录远程数据源版本号
}

# ---- 记录依赖链 ----

# GO 依赖
if "GO" in databases:
    go_basic_ver = go_version or downloader.get_latest_go_version()
    build_manifest["dependencies"]["GO"] = {
        "basic_dir": f"basic/go/{go_basic_ver}",
        "files": ["gene2go.gz", "gene_info.gz", "go-basic.obo"],
        "goa_fallback": goa_used if 'goa_used' in dir() else False,
    }
    # 从 versions.json 读取远程版本号
    go_obo_ver = vm.get_local_version("go_obo")
    if go_obo_ver:
        build_manifest["source_versions"]["go_obo"] = go_obo_ver.remote_version
    gene2go_ver = vm.get_local_version("gene2go")
    if gene2go_ver:
        build_manifest["source_versions"]["gene2go"] = gene2go_ver.remote_last_modified

# Reactome 依赖
if "Reactome" in databases:
    reactome_basic_ver = reactome_version or downloader.get_latest_reactome_version()
    build_manifest["dependencies"]["Reactome"] = {
        "basic_dir": f"basic/reactome/{reactome_basic_ver}",
        "files": ["NCBI2Reactome_All_Levels.txt.gz", "gene_info.gz"],
    }
    reactome_ver = vm.get_local_version("reactome")
    if reactome_ver:
        build_manifest["source_versions"]["reactome"] = reactome_ver.remote_version

# KEGG 依赖（REST API 实时获取，记录获取时间）
if "KEGG" in databases:
    kegg_ver = vm.get_local_version("kegg")
    build_manifest["dependencies"]["KEGG"] = {
        "source": "REST API (real-time)",
        "gene_info_from": f"basic/go/{go_version or downloader.get_latest_go_version()}/gene_info.gz",
    }
    if kegg_ver:
        build_manifest["source_versions"]["kegg"] = kegg_ver.remote_version

# DO 依赖
if "DO" in databases:
    build_manifest["dependencies"]["DO"] = {
        "basic_dir": "basic/do",
        "files": ["human_disease_knowledge_filtered.tsv.gz",
                   "human_disease_experiments_filtered.tsv.gz",
                   "human_disease_textmining_filtered.tsv.gz"],
        "gene_info_from": f"basic/go/{go_version or downloader.get_latest_go_version()}/gene_info.gz",
    }

# 写入
manifest_path = outdir / "build_manifest.json"
with open(manifest_path, "w", encoding="utf-8") as f:
    json.dump(build_manifest, f, indent=2, ensure_ascii=False)
logger.info("构建清单已保存: %s", manifest_path)
```

- [ ] **Step 2: 在 version.py 中添加血缘查询方法**

```python
class DatabaseVersionManager:
    # ... 现有代码 ...

    def get_build_lineage(self, organism_version: str, species: str) -> Optional[dict]:
        """查询某个 organism 版本的构建血缘

        Args:
            organism_version: 如 "v20260515"
            species: 如 "hsa"

        Returns:
            build_manifest.json 的内容，或 None
        """
        manifest_path = (
            self.database_dir / "organism" / organism_version / species / "build_manifest.json"
        )
        if not manifest_path.exists():
            return None
        with open(manifest_path, "r", encoding="utf-8") as f:
            return json.load(f)

    def get_full_lineage_report(self) -> str:
        """生成所有 organism 版本的血缘报告

        Returns:
            格式化的血缘报告字符串
        """
        lines = []
        lines.append("构建血缘追踪报告")
        lines.append("=" * 80)

        for org_ver in self.list_installed_organism_versions():
            for species in self.get_organism_build_info(org_ver).get(org_ver, []):
                lineage = self.get_build_lineage(org_ver, species)
                if not lineage:
                    continue

                lines.append(f"\n[{org_ver}/{species}]")
                lines.append(f"  构建时间: {lineage.get('built_at', '-')}")
                lines.append(f"  软件版本: {lineage.get('allenricher_version', '-')}")
                lines.append(f"  数据库: {', '.join(lineage.get('databases', []))}")

                deps = lineage.get("dependencies", {})
                if deps:
                    lines.append(f"  依赖链:")
                    for db_name, dep_info in deps.items():
                        basic_dir = dep_info.get("basic_dir", "-")
                        lines.append(f"    {db_name:<12} ← {basic_dir}")

                src_vers = lineage.get("source_versions", {})
                if src_vers:
                    lines.append(f"  源数据版本:")
                    for src_name, src_ver in src_vers.items():
                        lines.append(f"    {src_name:<12} = {src_ver}")

        lines.append("=" * 80)
        return "\n".join(lines)
```

- [ ] **Step 3: 在 list-versions 命令中支持 --lineage 参数**

```python
sub_list_ver.add_argument(
    "--lineage",
    action="store_true",
    help="显示构建血缘追踪信息",
)
```

```python
def cmd_list_versions(args):
    from allenricher.database.version import DatabaseVersionManager
    db_dir = args.database_dir or "./database"
    vm = DatabaseVersionManager(db_dir)

    if args.json:
        import json
        print(json.dumps(vm.get_summary_json(), indent=2, ensure_ascii=False))
    elif args.lineage:
        print(vm.get_full_lineage_report())
    else:
        print(vm.get_summary_table())
```

- [ ] **Step 4: 验证**

Run: `cd f:\OneDrive\Documents\TraeSOLO\AllEnricher\AllEnricher-v2 && python -m allenricher list-versions --lineage`
Expected: 显示每个 organism 版本依赖的 basic 数据源版本

- [ ] **Step 5: Commit**

```bash
git add allenricher/database/version.py allenricher/database/builder.py allenricher/cli.py
git commit -m "feat(version): 构建血缘追踪，build_manifest 记录完整依赖链"
```

---

## Task 12: 分析结果版本记录

**Files:**
- Modify: `allenricher/core/enrichment.py:155-188` (to_dict 方法)
- Modify: `allenricher/core/enrichment.py:1822-1842` (save_results 方法)
- Modify: `allenricher/report/generator.py:211,1167` (HTML 报告 header)

**目标**: 每次 analyze 的输出（TSV/HTML）中自动嵌入使用的数据库版本信息，确保结果可复现。

- [ ] **Step 1: 修改 EnrichmentResult.to_dict() 添加版本字段**

在 `allenricher/core/enrichment.py` 的 `to_dict()` 方法中，添加版本元数据：

```python
def to_dict(self) -> dict:
    """转换为字典（TSV 行格式）"""
    d = {
        "Term_ID": self.term_id,
        "Term_Name": self.term_name,
        "Database": self.database,
        # ... 现有字段 ...
    }
    return d
```

注意：TSV 的每一行不需要重复版本信息。版本信息应写入 TSV 文件的**注释行**（以 `#` 开头）。

- [ ] **Step 2: 修改 save_results() 在 TSV 头部写入版本注释**

在 `save_results()` 方法中，在 DataFrame 写入之前，先写入版本注释行：

```python
def save_results(self, filepath: str, metadata: dict = None) -> None:
    """保存富集结果到 TSV 文件

    Args:
        filepath: 输出文件路径
        metadata: 版本元数据字典（包含 database_version, source_versions 等）
    """
    df = pd.DataFrame([r.to_dict() for r in self.results])

    # 写入版本元数据注释行
    if metadata:
        with open(filepath, "w", encoding="utf-8") as f:
            f.write(f"# AllEnricher version: {metadata.get('allenricher_version', 'unknown')}\n")
            f.write(f"# Analysis date: {metadata.get('analysis_date', 'unknown')}\n")
            f.write(f"# Database version: {metadata.get('database_version', 'unknown')}\n")
            f.write(f"# Species: {metadata.get('species', 'unknown')}\n")

            source_versions = metadata.get("source_versions", {})
            if source_versions:
                f.write("# Source data versions:\n")
                for src_name, src_ver in source_versions.items():
                    f.write(f"#   {src_name}: {src_ver}\n")

            f.write("#\n")
            df.to_csv(f, sep="\t", index=False)
    else:
        df.to_csv(filepath, sep="\t", index=False)
```

- [ ] **Step 3: 修改 DatabaseManager 保留版本信息**

修改 `DatabaseManager._find_species_dir()` 使其同时返回版本号：

```python
def _find_species_dir(self, database_dir, species, version=None):
    """查找物种目录，同时记录匹配的版本号"""
    organism_dir = database_dir / "organism"

    if version and organism_dir.exists():
        species_dir = organism_dir / version / species
        if species_dir.exists():
            self._active_version = version  # 记录当前使用的版本
            return species_dir

    if organism_dir.exists():
        for version_dir in sorted(organism_dir.iterdir(), reverse=True):
            if version_dir.is_dir():
                species_dir = version_dir / species
                if species_dir.exists():
                    self._active_version = version_dir.name  # 记录版本号
                    return species_dir

    # v1 兼容
    if (database_dir / f"{species}.GO2gene.tab.gz").exists():
        self._active_version = "v1-legacy"
        return database_dir

    self._active_version = None
    return database_dir

@property
def active_version(self) -> Optional[str]:
    """当前加载的数据库版本号"""
    return getattr(self, "_active_version", None)

def get_build_metadata(self) -> Optional[dict]:
    """读取当前版本的 build_manifest.json"""
    if not self._active_version:
        return None
    manifest_path = (
        Path(self.database_dir) / "organism" / self._active_version /
        self.species / "build_manifest.json"
    )
    if not manifest_path.exists():
        return None
    with open(manifest_path, "r", encoding="utf-8") as f:
        return json.load(f)
```

- [ ] **Step 4: 修改 cmd_analyze() 传递版本元数据**

在 `cmd_analyze()` 中，构建 metadata 并传递给 save_results：

```python
# 在分析完成后、保存结果之前
import json
from datetime import datetime, timezone
import allenricher

metadata = {
    "allenricher_version": allenricher.__version__,
    "analysis_date": datetime.now(timezone.utc).isoformat(),
    "database_version": db_manager.active_version,
    "species": config.species,
    "databases": config.databases,
}

# 从 build_manifest 读取源数据版本
build_meta = db_manager.get_build_metadata()
if build_meta:
    metadata["source_versions"] = build_meta.get("source_versions", {})
    metadata["built_at"] = build_meta.get("built_at", "")

# 传递给 save_results
result.save_results(output_tsv, metadata=metadata)
```

- [ ] **Step 5: 修改 HTML 报告嵌入版本信息**

在 `allenricher/report/generator.py` 中，将硬编码的版本替换为动态获取：

```python
# 替换第 211 行和第 1167 行
# 原: Version 2.0 | {datetime.now().strftime("%Y-%m-%d")}
# 改:
header_text = f"Version {allenricher.__version__} | {datetime.now().strftime('%Y-%m-%d')}"
if metadata and metadata.get("database_version"):
    header_text += f" | DB: {metadata['database_version']}"
```

- [ ] **Step 6: 验证**

```bash
# 运行分析
cd f:\OneDrive\Documents\TraeSOLO\AllEnricher\AllEnricher-v2
python -m allenricher analyze -i '../AllEnricher-v1/example/example.glist' -s hsa -d GO -o test_output/version_test

# 检查 TSV 文件头部
head -10 test_output/version_test/GO_enrichment.tsv
```

Expected: TSV 文件前几行包含 `# AllEnricher version:`, `# Database version:`, `# Source data versions:` 等注释

- [ ] **Step 7: Commit**

```bash
git add allenricher/core/enrichment.py allenricher/database/manager.py allenricher/report/generator.py allenricher/cli.py
git commit -m "feat(analyze): 分析结果中记录数据库版本信息，确保可复现"
```

---

## Task 13: 版本切换与回退

**Files:**
- Modify: `allenricher/cli.py`
- Modify: `allenricher/core/config.py`

**目标**: 支持用户通过 `--use-version` 参数指定任意已安装的本地版本进行分析。

- [ ] **Step 1: 在 analyze 命令添加 --use-version 参数**

```python
# 在 analyze 子命令的参数定义中添加
sub_analyze.add_argument(
    "--use-version",
    type=str,
    default=None,
    help="指定使用的数据库版本（如 v20260515），默认使用最新版本",
)
```

- [ ] **Step 2: 在 config.py 中添加 use_version 字段**

```python
# 在 EnrichmentConfig 数据类中添加
use_version: Optional[str] = None  # 指定使用的数据库版本
```

- [ ] **Step 3: 修改 cmd_analyze() 传递版本参数**

```python
def cmd_analyze(args):
    # ... 现有代码 ...

    # 版本锁定
    if args.use_version:
        config.use_version = args.use_version

    # 创建 DatabaseManager 时传递版本
    db_manager = DatabaseManager(db_dir, config.species)
    db_manager.load_databases(
        config.databases,
        version=config.use_version  # 传递版本锁定参数
    )
```

- [ ] **Step 4: 添加版本不存在时的友好提示**

```python
# 在 _find_species_dir() 中，当指定版本不存在时
if version and organism_dir.exists():
    species_dir = organism_dir / version / species
    if not species_dir.exists():
        # 列出可用版本
        available = sorted(
            d.name for d in organism_dir.iterdir()
            if d.is_dir() and (d / species).exists()
        )
        if available:
            logger.error(
                "版本 '%s' 的物种 '%s' 不存在。可用版本: %s",
                version, species, ", ".join(available)
            )
        else:
            logger.error("物种 '%s' 没有任何已构建的版本。请先运行 allenricher build。", species)
        raise FileNotFoundError(f"Database version {version}/{species} not found")
```

- [ ] **Step 5: 验证版本切换**

```bash
# 使用指定版本分析
python -m allenricher analyze -i example.glist -s hsa -d GO --use-version v20260515 -o test_output/v20260515_test

# 检查输出中的版本信息
head -5 test_output/v20260515_test/GO_enrichment.tsv
```

Expected: TSV 注释中显示 `# Database version: v20260515`

- [ ] **Step 6: 验证版本不存在时的错误提示**

```bash
python -m allenricher analyze -i example.glist -s hsa -d GO --use-version v99999999 -o test_output/error_test
```

Expected: 显示 "版本 'v99999999' 的物种 'hsa' 不存在。可用版本: v20260515"

- [ ] **Step 7: Commit**

```bash
git add allenricher/cli.py allenricher/core/config.py allenricher/database/manager.py
git commit -m "feat(version): 支持 --use-version 版本切换与回退"
```

---

## Task 14: 端到端验证（完整版）

**Files:**
- Test: 手动 E2E 测试

- [ ] **Step 1: 验证 list-versions 命令**

```bash
python -m allenricher list-versions
python -m allenricher list-versions --lineage
python -m allenricher list-versions --json | python -m json.tool | head -30
```

- [ ] **Step 2: 验证 check-update 命令**

```bash
python -m allenricher check-update
```

- [ ] **Step 3: 验证版本切换**

```bash
python -m allenricher analyze -i example.glist -s hsa -d GO --use-version v20260515 -o test_output/e2e_versioned
head -8 test_output/e2e_versioned/GO_enrichment.tsv
```

Expected: TSV 注释行包含 `# Database version: v20260515` 和源数据版本信息

- [ ] **Step 4: 验证 cleanup 命令**

```bash
python -m allenricher cleanup --dry-run --keep 1
```

- [ ] **Step 5: 运行完整测试套件**

```bash
python -m pytest tests/test_version.py -v
python -m pytest tests/ -v --tb=short 2>&1 | tail -20
```

Expected: 全部通过

---

## 执行选项

**Plan complete and saved to `docs/superpowers/plans/2026-05-27-database-version-management.md`. Two execution options:**

**1. Subagent-Driven (recommended)** - I dispatch a fresh subagent per task, review between tasks, fast iteration

**2. Inline Execution** - Execute tasks in this session using executing-plans, batch execution with checkpoints

**Which approach?**
