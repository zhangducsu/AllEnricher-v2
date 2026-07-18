"""
Test of module for database version management

Test DatabaseVersion, VersionManifest, DatabaseVersionManager, RemoteVersionChecker
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
    """Data class testing for DatabaseVersion"""

    def test_from_dict_roundtrip(self):
        """Test dictionaries/reverse sequenced returns"""
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
        # To_dict() contains all fields (including file_hash=Noe), verifying that key fields are consistent
        result = ver.to_dict()
        assert result["source"] == data["source"]
        assert result["remote_version"] == data["remote_version"]
        assert result["remote_last_modified"] == data["remote_last_modified"]
        assert result["local_version"] == data["local_version"]
        assert result["downloaded_at"] == data["downloaded_at"]
        assert result["file_hash"] is None

    def test_is_newer_than(self):
        """Test version of old comparison"""
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
        """Returns False when testing missing date"""
        a = DatabaseVersion(source="go")
        b = DatabaseVersion(source="go", remote_last_modified="Thu, 02 May 2026 02:53:36 GMT")
        assert a.is_newer_than(b) is False
        assert b.is_newer_than(a) is False


# ============================================================
# TestVersionManifest
# ============================================================

class TestVersionManifest:
    """VersionManifest Data Class Test"""

    def test_save_and_load(self, tmp_path):
        """Test to save and load the list of versions"""
        manifest = VersionManifest()
        manifest.set("go", DatabaseVersion(source="go", local_version="GO20260527"))
        manifest.set("kegg", DatabaseVersion(source="kegg", local_version="cached"))

        path = tmp_path / "versions.json"
        manifest.save(path)

        loaded = VersionManifest.load(path)
        assert loaded.get("go").local_version == "GO20260527"
        assert loaded.get("kegg").local_version == "cached"

    def test_load_nonexistent(self, tmp_path):
        """Test load non-existent file returns empty list"""
        manifest = VersionManifest.load(tmp_path / "missing.json")
        assert len(manifest.versions) == 0


# ============================================================
# TestDatabaseVersionManager
# ============================================================

class TestDatabaseVersionManager:
    """Database VersionManager Test"""

    def test_record_and_retrieve(self, tmp_path):
        """Test log download and retrieve version of information"""
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
        """Test version list duration (read across instance)"""
        mgr = DatabaseVersionManager(str(tmp_path))
        mgr.record_download("go", "GO20260527", "basic/go/GO20260527")
        mgr2 = DatabaseVersionManager(str(tmp_path))
        assert mgr2.get_local_version("go").local_version == "GO20260527"

    def test_list_installed_basic_versions(self, tmp_path):
        """Test scan Basic/ Directory Lists installed version"""
        # Create Simulate Directory Structure
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
        """Test scan organism/ directory list of constructed versions"""
        (tmp_path / "organism" / "v20260515" / "hsa").mkdir(parents=True)
        (tmp_path / "organism" / "v20260527" / "hsa").mkdir(parents=True)

        mgr = DatabaseVersionManager(str(tmp_path))
        versions = mgr.list_installed_organism_versions()

        assert versions == ["v20260515", "v20260527"]

    def test_list_local_versions(self, tmp_path):
        """Test list_local_versions to return all recorded versions"""
        mgr = DatabaseVersionManager(str(tmp_path))
        mgr.record_download("go", "GO20260527", "basic/go/GO20260527")
        mgr.record_download("kegg", "cached", "basic/kegg")
        versions = mgr.list_local_versions()
        assert "go" in versions
        assert "kegg" in versions
        assert len(versions) == 2

    def test_get_organism_build_info(self, tmp_path):
        """Test get_organism_build_info to return the list of species under the version"""
        org_dir = tmp_path / "organism" / "v20260515"
        org_dir.mkdir(parents=True)
        (org_dir / "hsa").mkdir()
        (org_dir / "mmu").mkdir()
        mgr = DatabaseVersionManager(str(tmp_path))
        info = mgr.get_organism_build_info("v20260515")
        assert "v20260515" in info
        assert set(info["v20260515"]) == {"hsa", "mmu"}

    def test_find_stale_versions(self, tmp_path):
        """Test find_stale_versions Correct recognition of old versions"""
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
        """Test not remove files"""
        for ver in ["GO20260501", "GO20260515"]:
            (tmp_path / "basic" / "go" / ver).mkdir(parents=True)
        mgr = DatabaseVersionManager(str(tmp_path))
        removed = mgr.remove_stale_versions(keep_count=1, dry_run=True)
        assert "GO20260501" in removed["go"]
        assert (tmp_path / "basic" / "go" / "GO20260501").exists()

    def test_remove_stale_versions_actual(self, tmp_path):
        """Test remove_stale_versions"""
        for ver in ["GO20260501", "GO20260515"]:
            (tmp_path / "basic" / "go" / ver).mkdir(parents=True)
        mgr = DatabaseVersionManager(str(tmp_path))
        removed = mgr.remove_stale_versions(keep_count=1, dry_run=False)
        assert "GO20260501" in removed["go"]
        assert not (tmp_path / "basic" / "go" / "GO20260501").exists()
        assert (tmp_path / "basic" / "go" / "GO20260515").exists()

    def test_get_build_lineage(self, tmp_path):
        """Test get built_lineage to read"""
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
        """Test get_Build_lineage returns None"""
        mgr = DatabaseVersionManager(str(tmp_path))
        lineage = mgr.get_build_lineage("v99999999", "hsa")
        assert lineage is None

    def test_get_summary_json(self, tmp_path):
        """Test get_summary_json returns structured data"""
        mgr = DatabaseVersionManager(str(tmp_path))
        mgr.record_download("go", "GO20260527", "basic/go/GO20260527")
        summary = mgr.get_summary_json()
        assert "basic_versions" in summary
        assert "organism_versions" in summary
        assert "version_records" in summary
        assert "go" in summary["version_records"]

    def test_get_summary_table(self, tmp_path):
        """Test get_sumary_table returns non-empty string"""
        for ver in ["GO20260527"]:
            (tmp_path / "basic" / "go" / ver).mkdir(parents=True)
        mgr = DatabaseVersionManager(str(tmp_path))
        table = mgr.get_summary_table()
        assert "List of local database versions" in table
        assert "go" in table

    def test_get_full_lineage_report(self, tmp_path):
        """Test get_full_line_report to generate blood report"""
        manifest_dir = tmp_path / "organism" / "v20260515" / "hsa"
        manifest_dir.mkdir(parents=True)
        with open(manifest_dir / "build_manifest.json", "w") as f:
            json.dump({"built_at": "2026-05-15", "species": "hsa", "databases": ["GO"]}, f)
        mgr = DatabaseVersionManager(str(tmp_path))
        report = mgr.get_full_lineage_report()
        assert "Database build provenance" in report
        assert "v20260515/hsa" in report


# ============================================================
# TestRemoteVersionChecker
# ============================================================

class TestRemoteVersionChecker:
    """RemoteVersionChecker test (in connection with network requests)"""

    @pytest.mark.slow
    def test_check_kegg_version(self):
        """Test KEGREST API version for testing"""
        checker = RemoteVersionChecker()
        result = checker.check_kegg_version()
        assert result is not None
        assert "remote_version" in result
        assert result["remote_version"]

    def test_check_kegg_version_without_release_line(self):
        checker = RemoteVersionChecker()
        response = MagicMock()
        response.text = "kegg\tKEGG\n\tpathway 587 2026/07/10\n\tgenes 1 2026/07/12\n"
        response.raise_for_status = MagicMock()
        with patch("requests.get", return_value=response):
            result = checker.check_kegg_version()
        assert result == {
            "remote_version": "KEGG data 2026/07/12",
            "last_modified": "2026/07/12",
        }

    @pytest.mark.slow
    def test_check_go_obo_version(self):
        """Test the GO OBO file version for detection (marked as slow)"""
        checker = RemoteVersionChecker()
        result = checker.check_go_obo_version()
        assert result is not None
        assert "remote_version" in result
        assert "releases/" in result["remote_version"]

    def test_check_head_success(self):
        """Test check_head successfully parsed Last-Modified"""
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
        """Test check_head failed to return None"""
        checker = RemoteVersionChecker()
        with patch("requests.head", side_effect=Exception("Connection error")):
            result = checker.check_head("https://example.com/file.gz")
        assert result is None

    def test_check_updates_no_local(self, tmp_path):
        """Test check_updates to display all updates when local records are not available"""
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
        """Test check_updates Local is updated and no updates are shown"""
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

    def test_check_updates_without_local_remote_metadata(self, tmp_path):
        mgr = DatabaseVersionManager(str(tmp_path))
        mgr.record_download(
            source="gene2go", local_version="GO20260527",
            local_path="basic/go/GO20260527",
        )
        checker = RemoteVersionChecker()
        remote = {"gene2go": {"last_modified": "Sun, 12 Jul 2026 05:29:38 GMT"}}
        with patch.object(checker, "check_all_sources", return_value=remote):
            status = checker.check_updates(mgr)
        assert status["gene2go"]["has_update"] is True
