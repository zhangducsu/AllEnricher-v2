"""
数据库版本管理模块单元测试

测试 DatabaseVersion、VersionManifest、DatabaseVersionManager、RemoteVersionChecker
"""

import json
import pytest
from pathlib import Path
from unittest.mock import patch, MagicMock

from allenricher.database.version import (
    DatabaseVersion,
    VersionManifest,
    DatabaseVersionManager,
    RemoteVersionChecker,
)


# ============================================================
# TestDatabaseVersion
# ============================================================

class TestDatabaseVersion:
    """DatabaseVersion 数据类测试"""

    def test_from_dict_roundtrip(self):
        """测试字典序列化/反序列化往返一致性"""
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
        # to_dict() 包含所有字段（含 file_hash=None），验证关键字段一致
        result = ver.to_dict()
        assert result["source"] == data["source"]
        assert result["remote_version"] == data["remote_version"]
        assert result["remote_last_modified"] == data["remote_last_modified"]
        assert result["local_version"] == data["local_version"]
        assert result["downloaded_at"] == data["downloaded_at"]
        assert result["file_hash"] is None

    def test_is_newer_than(self):
        """测试版本新旧比较"""
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
        """测试缺少日期时比较返回 False"""
        a = DatabaseVersion(source="go")
        b = DatabaseVersion(source="go", remote_last_modified="Thu, 02 May 2026 02:53:36 GMT")
        assert a.is_newer_than(b) is False
        assert b.is_newer_than(a) is False


# ============================================================
# TestVersionManifest
# ============================================================

class TestVersionManifest:
    """VersionManifest 数据类测试"""

    def test_save_and_load(self, tmp_path):
        """测试保存和加载版本清单"""
        manifest = VersionManifest()
        manifest.set("go", DatabaseVersion(source="go", local_version="GO20260527"))
        manifest.set("kegg", DatabaseVersion(source="kegg", local_version="cached"))

        path = tmp_path / "versions.json"
        manifest.save(path)

        loaded = VersionManifest.load(path)
        assert loaded.get("go").local_version == "GO20260527"
        assert loaded.get("kegg").local_version == "cached"

    def test_load_nonexistent(self, tmp_path):
        """测试加载不存在的文件返回空清单"""
        manifest = VersionManifest.load(tmp_path / "missing.json")
        assert len(manifest.versions) == 0


# ============================================================
# TestDatabaseVersionManager
# ============================================================

class TestDatabaseVersionManager:
    """DatabaseVersionManager 测试"""

    def test_record_and_retrieve(self, tmp_path):
        """测试记录下载并检索版本信息"""
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
        """测试版本清单持久化（跨实例读取）"""
        mgr = DatabaseVersionManager(str(tmp_path))
        mgr.record_download("go", "GO20260527", "basic/go/GO20260527")
        mgr2 = DatabaseVersionManager(str(tmp_path))
        assert mgr2.get_local_version("go").local_version == "GO20260527"

    def test_list_installed_basic_versions(self, tmp_path):
        """测试扫描 basic/ 目录列出已安装版本"""
        # 创建模拟目录结构
        (tmp_path / "basic" / "go" / "GO20260515").mkdir(parents=True)
        (tmp_path / "basic" / "go" / "GO20260527").mkdir(parents=True)
        (tmp_path / "basic" / "reactome" / "Reactome20260515").mkdir(parents=True)
        (tmp_path / "basic" / "kegg").mkdir(parents=True)
        (tmp_path / "basic" / "kegg" / "hsa_gene2pathway.txt").write_text("test")

        mgr = DatabaseVersionManager(str(tmp_path))
        versions = mgr.list_installed_basic_versions()

        assert "go" in versions
        assert versions["go"] == ["GO20260515", "GO20260527"]
        assert "reactome" in versions
        assert versions["reactome"] == ["Reactome20260515"]
        assert "kegg" in versions
        assert versions["kegg"] == ["cached"]

    def test_list_installed_organism_versions(self, tmp_path):
        """测试扫描 organism/ 目录列出已构建版本"""
        (tmp_path / "organism" / "v20260515" / "hsa").mkdir(parents=True)
        (tmp_path / "organism" / "v20260527" / "hsa").mkdir(parents=True)

        mgr = DatabaseVersionManager(str(tmp_path))
        versions = mgr.list_installed_organism_versions()

        assert versions == ["v20260515", "v20260527"]

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
        for ver in ["GO20260501", "GO20260515", "GO20260527"]:
            (tmp_path / "basic" / "go" / ver).mkdir(parents=True)
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
        assert (tmp_path / "basic" / "go" / "GO20260501").exists()

    def test_remove_stale_versions_actual(self, tmp_path):
        """测试 remove_stale_versions 实际删除"""
        for ver in ["GO20260501", "GO20260515"]:
            (tmp_path / "basic" / "go" / ver).mkdir(parents=True)
        mgr = DatabaseVersionManager(str(tmp_path))
        removed = mgr.remove_stale_versions(keep_count=1, dry_run=False)
        assert "GO20260501" in removed["go"]
        assert not (tmp_path / "basic" / "go" / "GO20260501").exists()
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
        with open(manifest_dir / "build_manifest.json", "w") as f:
            json.dump({"built_at": "2026-05-15", "species": "hsa", "databases": ["GO"]}, f)
        mgr = DatabaseVersionManager(str(tmp_path))
        report = mgr.get_full_lineage_report()
        assert "构建血缘追踪报告" in report
        assert "v20260515/hsa" in report


# ============================================================
# TestRemoteVersionChecker
# ============================================================

class TestRemoteVersionChecker:
    """RemoteVersionChecker 测试（涉及网络请求）"""

    @pytest.mark.slow
    def test_check_kegg_version(self):
        """测试 KEGG REST API 版本检测"""
        checker = RemoteVersionChecker()
        result = checker.check_kegg_version()
        assert result is not None
        assert "remote_version" in result
        assert "Release" in result["remote_version"]

    @pytest.mark.slow
    def test_check_go_obo_version(self):
        """测试 GO OBO 文件版本检测（标记为 slow）"""
        checker = RemoteVersionChecker()
        result = checker.check_go_obo_version()
        assert result is not None
        assert "remote_version" in result
        assert "releases/" in result["remote_version"]

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
