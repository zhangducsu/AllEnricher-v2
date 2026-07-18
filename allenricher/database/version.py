"""Record local database versions, lineage, and remote update metadata."""

from __future__ import annotations

import json
import logging
import re
import shutil
from dataclasses import dataclass, field, asdict
from datetime import datetime, timezone
from email.utils import parsedate_to_datetime
from pathlib import Path
from typing import Dict, List, Optional

import requests

logger = logging.getLogger(__name__)


# ============================================================
# Data class
# ============================================================

@dataclass
class DatabaseVersion:
    """Describe one downloaded upstream data snapshot."""

    source: str
    remote_version: Optional[str] = None
    remote_last_modified: Optional[str] = None
    local_version: Optional[str] = None
    local_path: Optional[str] = None
    downloaded_at: Optional[str] = None
    file_hash: Optional[str] = None

    def is_newer_than(self, other: DatabaseVersion) -> bool:
        """Return whether this version is newer than another snapshot."""
        if not self.remote_last_modified or not other.remote_last_modified:
            return False
        try:
            remote_dt = parsedate_to_datetime(self.remote_last_modified)
            local_dt = parsedate_to_datetime(other.remote_last_modified)
            return remote_dt > local_dt
        except (ValueError, TypeError):
            return False

    def to_dict(self) -> dict:
        """Return a JSON-serializable representation."""
        return asdict(self)

    @classmethod
    def from_dict(cls, data: dict) -> DatabaseVersion:
        """Create an instance from a serialized mapping."""
        return cls(**{k: v for k, v in data.items() if k in cls.__dataclass_fields__})


@dataclass
class VersionManifest:
    """Store version metadata for all downloaded data sources."""

    created_at: str = field(default_factory=lambda: datetime.now(timezone.utc).isoformat())
    updated_at: str = field(default_factory=lambda: datetime.now(timezone.utc).isoformat())
    versions: Dict[str, DatabaseVersion] = field(default_factory=dict)

    def get(self, source: str) -> Optional[DatabaseVersion]:
        """Return version metadata for one data source."""
        return self.versions.get(source)

    def set(self, source: str, version: DatabaseVersion) -> None:
        """Record version metadata for one data source."""
        self.versions[source] = version
        self.updated_at = datetime.now(timezone.utc).isoformat()

    def to_dict(self) -> dict:
        """Return a JSON-serializable representation."""
        return {
            "created_at": self.created_at,
            "updated_at": self.updated_at,
            "versions": {k: v.to_dict() for k, v in self.versions.items()},
        }

    @classmethod
    def from_dict(cls, data: dict) -> VersionManifest:
        """Create an instance from a serialized mapping."""
        manifest = cls(
            created_at=data.get("created_at", ""),
            updated_at=data.get("updated_at", ""),
        )
        for source, ver_data in data.get("versions", {}).items():
            manifest.versions[source] = DatabaseVersion.from_dict(ver_data)
        return manifest

    def save(self, path: Path) -> None:
        """Write the species registry as deterministic TSV."""
        path.parent.mkdir(parents=True, exist_ok=True)
        with open(path, "w", encoding="utf-8") as f:
            json.dump(self.to_dict(), f, indent=2, ensure_ascii=False)
        logger.info("Version List saved: %s", path)

    @classmethod
    def load(cls, path: Path) -> VersionManifest:
        """Load a species registry from TSV."""
        if not path.exists():
            logger.info("Version list does not exist, creating a new list: %s", path)
            return cls()
        with open(path, "r", encoding="utf-8") as f:
            return cls.from_dict(json.load(f))


# ============================================================
# Version Manager
# ============================================================

class DatabaseVersionManager:
    """Record versions, builds, lineage, and cleanup candidates."""

    MANIFEST_FILENAME = "versions.json"

    def __init__(self, database_dir: str = "./database"):
        """
        Args:
database_dir: Path to the database root directory
        """
        self.database_dir = Path(database_dir)
        self.manifest_path = self.database_dir / self.MANIFEST_FILENAME
        self._manifest: Optional[VersionManifest] = None

    @property
    def manifest(self) -> VersionManifest:
        """Load the version manifest on first access."""
        if self._manifest is None:
            self._manifest = VersionManifest.load(self.manifest_path)
        return self._manifest

    def save_manifest(self) -> None:
        """Write the current version manifest."""
        self.manifest.save(self.manifest_path)

    def record_download(
        self,
        source: str,
        local_version: str,
        local_path: str,
        remote_version: Optional[str] = None,
        remote_last_modified: Optional[str] = None,
    ) -> None:
        """Record one completed source download."""
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
        """Return local version metadata for one data source."""
        ver = self.manifest.get(source)
        if ver is not None:
            return ver

        # Backs: Scan disk directory extrapolated version
        inferred = self._infer_version_from_disk(source)
        if inferred is not None:
            # Write back versions.json to follow up on reading directly
            self.manifest.set(source, inferred)
            self.save_manifest()
        return inferred

    def _infer_version_from_disk(self, source: str) -> Optional[DatabaseVersion]:
        """Infer legacy version metadata from downloaded directories."""
        basic_dir = self.database_dir / "basic"

        if source == "go":
            go_dir = basic_dir / "go"
            if not go_dir.exists():
                return None
            versions = sorted(
                d.name for d in go_dir.iterdir()
                if d.is_dir() and d.name.startswith("GO")
            )
            if not versions:
                return None
            latest = versions[-1]
            return DatabaseVersion(
                source="go",
                local_version=latest,
                local_path=f"basic/go/{latest}",
                downloaded_at="",
            )

        elif source == "reactome":
            re_dir = basic_dir / "reactome"
            if not re_dir.exists():
                return None
            versions = sorted(
                d.name for d in re_dir.iterdir()
                if d.is_dir() and d.name.startswith("Reactome")
            )
            if not versions:
                return None
            latest = versions[-1]
            return DatabaseVersion(
                source="reactome",
                local_version=latest,
                local_path=f"basic/reactome/{latest}",
                downloaded_at="",
            )

        elif source == "do":
            do_dir = basic_dir / "do"
            if do_dir.exists() and any(do_dir.iterdir()):
                return DatabaseVersion(
                    source="do",
                    local_version="cached",
                    local_path="basic/do",
                    downloaded_at="",
                )
            return None

        elif source == "kegg":
            kegg_dir = basic_dir / "kegg"
            if kegg_dir.exists() and any(kegg_dir.iterdir()):
                return DatabaseVersion(
                    source="kegg",
                    local_version="cached",
                    local_path="basic/kegg",
                    downloaded_at="",
                )
            return None

        elif source == "taxonomy":
            tax_dir = basic_dir / "taxonomy"
            if tax_dir.exists() and (tax_dir / "names.dmp").exists():
                return DatabaseVersion(
                    source="taxonomy",
                    local_version="cached",
                    local_path="basic/taxonomy",
                    downloaded_at="",
                )
            return None

        return None

    def list_local_versions(self) -> Dict[str, DatabaseVersion]:
        """Return locally installed source versions."""
        return dict(self.manifest.versions)

    def list_installed_basic_versions(self) -> Dict[str, List[str]]:
        """List installed shared-source snapshots."""
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

        # KEGG (no version of directory, check if file is available)
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
        """List installed species database builds."""
        organism_dir = self.database_dir / "organism"
        if not organism_dir.exists():
            return []
        return sorted(
            d.name for d in organism_dir.iterdir()
            if d.is_dir() and d.name.startswith("v")
        )

    def get_organism_build_info(self, version: str) -> Dict[str, List[str]]:
        """Return species included in one database build."""
        version_dir = self.database_dir / "organism" / version
        if not version_dir.exists():
            return {}
        return {
            version: sorted(
                d.name for d in version_dir.iterdir() if d.is_dir()
            )
        }

    def find_stale_versions(self, keep_count: int = 2) -> Dict[str, List[str]]:
        """Return old snapshots eligible for cleanup."""
        stale: Dict[str, List[str]] = {}

        # Basic data
        basic_versions = self.list_installed_basic_versions()
        for source, versions in basic_versions.items():
            if len(versions) > keep_count:
                stale[source] = versions[:-keep_count]

        # Species database
        organism_versions = self.list_installed_organism_versions()
        if len(organism_versions) > keep_count:
            stale["organism"] = organism_versions[:-keep_count]

        return stale

    def remove_stale_versions(self, keep_count: int = 2, dry_run: bool = True) -> Dict[str, List[str]]:
        """Remove selected stale snapshots after caller confirmation."""
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
                        shutil.rmtree(dir_path)
                        logger.info("Deleted: %s", dir_path)
                    else:
                        logger.info("[dry-run] Delete: %s", dir_path)
                    removed[source].append(ver)

        return removed

    def get_build_lineage(self, organism_version: str, species: str) -> Optional[dict]:
        """Return the source lineage of one species database build."""
        manifest_path = (
            self.database_dir / "organism" / organism_version / species / "build_manifest.json"
        )
        if not manifest_path.exists():
            return None
        with open(manifest_path, "r", encoding="utf-8") as f:
            return json.load(f)

    def get_full_lineage_report(self) -> str:
        """Return lineage metadata for all installed species builds."""
        lines = []
        lines.append("Database build provenance")
        lines.append("=" * 80)

        for org_ver in self.list_installed_organism_versions():
            for species in self.get_organism_build_info(org_ver).get(org_ver, []):
                lineage = self.get_build_lineage(org_ver, species)
                if not lineage:
                    continue

                lines.append(f"\n[{org_ver}/{species}]")
                lines.append(f"Build time: {lineage.get('built_at', '-')}")
                lines.append(f"Software version: {lineage.get('allenricher_version', '-')}")
                lines.append(f"Database: {', '.join(lineage.get('databases', []))}")

                deps = lineage.get("dependencies", {})
                if deps:
                    lines.append("Dependencies:")
                    for db_name, dep_info in deps.items():
                        basic_dir = dep_info.get("basic_dir", "-")
                        lines.append(f"    {db_name:<12} <- {basic_dir}")

                src_vers = lineage.get("source_versions", {})
                if src_vers:
                    lines.append("Source versions:")
                    for src_name, src_ver in src_vers.items():
                        lines.append(f"    {src_name:<12} = {src_ver}")

        lines.append("=" * 80)
        return "\n".join(lines)

    def get_summary_table(self) -> str:
        """Render installed version metadata as a readable table."""
        lines = []
        lines.append("List of local database versions")
        lines.append("=" * 80)

        # Basic Data Version
        basic_versions = self.list_installed_basic_versions()
        if basic_versions:
            lines.append("\n [Base Data (basic/)]")
            lines.append(f"  {'Data Sources':<15} {'Installed version':<40}")
            lines.append(f"  {'-'*15} {'-'*40}")
            for source, versions in basic_versions.items():
                ver_str = ", ".join(versions) if versions else "None"
                latest = "<- Latest" if versions else ""
                lines.append(f"  {source:<15} {ver_str:<40} {latest}")

        # Species database version
        organism_versions = self.list_installed_organism_versions()
        if organism_versions:
            lines.append("\n[Specific database (organism/)]")
            lines.append(f"  {'Version':<15} {'Include species'}")
            lines.append(f"  {'-'*15} {'-'*40}")
            for ver in organism_versions:
                species_list = self.get_organism_build_info(ver).get(ver, [])
                species_str = ", ".join(species_list) if species_list else "Empty"
                latest = "<- Latest" if ver == organism_versions[-1] else ""
                lines.append(f"  {ver:<15} {species_str:<40} {latest}")

        # Remote version recorded in version.json
        local_records = self.list_local_versions()
        if local_records:
            lines.append("\n[version metadata (versions.json)]")
            lines.append(f"  {'Data Sources':<20} {'Local version':<20} {'Remote version':<25} {'Download Time'}")
            lines.append(f"  {'-'*20} {'-'*20} {'-'*25} {'-'*20}")
            for source, ver in sorted(local_records.items()):
                remote_ver = ver.remote_version or "-"
                downloaded = (ver.downloaded_at[:10] if ver.downloaded_at else "-")
                lines.append(f"  {source:<20} {ver.local_version or '-':<20} {remote_ver:<25} {downloaded}")

        lines.append("=" * 80)
        return "\n".join(lines)

    def get_summary_json(self) -> dict:
        """Render installed version metadata as JSON."""
        return {
            "basic_versions": self.list_installed_basic_versions(),
            "organism_versions": self.list_installed_organism_versions(),
            "version_records": {
                source: ver.to_dict()
                for source, ver in self.list_local_versions().items()
            },
        }


# ============================================================
# Remote version detector
# ============================================================

class RemoteVersionChecker:
    """Inspect upstream metadata for newer data snapshots."""

    TIMEOUT = 30  # sec

    # Key files for data sources URL
    SOURCE_URLS: Dict[str, str] = {
        "gene2go": "https://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2go.gz",
        "gene_info": "https://ftp.ncbi.nlm.nih.gov/gene/DATA/gene_info.gz",
        "go_obo": "http://purl.obolibrary.org/obo/go/go-basic.obo",
        "goa_proteomes": "https://ftp.ebi.ac.uk/pub/databases/GO/goa/proteomes/",
        "kegg": "https://rest.kegg.jp/info/kegg",
        "reactome": "https://reactome.org/download/",
        "taxonomy": "https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz",
        "wikipathways": "https://data.wikipathways.org/",
    }

    # Remote source key... and possible key list in local versions. json
    _SOURCE_KEY_ALIASES: Dict[str, List[str]] = {
        "gene2go": ["gene2go", "go"],
        "gene_info": ["gene_info", "go"],
        "go_obo": ["go_obo", "go"],
        "goa_proteomes": ["goa_proteomes"],
        "kegg": ["kegg"],
        "reactome": ["reactome"],
        "taxonomy": ["taxonomy"],
    }

    def __init__(self, timeout: int = TIMEOUT):
        """
        Args:
timeout: HTTP Request timeout in seconds
        """
        self.timeout = timeout

    def check_head(self, url: str) -> Optional[Dict[str, str]]:
        """Read Last-Modified and ETag metadata from an upstream URL."""
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
            logger.warning("HEAD request failed%s: %s", url, e)
            return None

    def check_go_obo_version(self) -> Optional[Dict[str, str]]:
        """Inspect the current upstream Gene Ontology release."""
        url = self.SOURCE_URLS["go_obo"]
        try:
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
            logger.warning("Could not close temporary folder: %s%s", e)
        return None

    def check_kegg_version(self) -> Optional[Dict[str, str]]:
        """Inspect the current upstream KEGG release metadata."""
        url = self.SOURCE_URLS["kegg"]
        try:
            resp = requests.get(url, timeout=self.timeout)
            resp.raise_for_status()
            for line in resp.text.strip().split("\n"):
                if "Release" in line:
                    # Remove Release and subsequent contents (remove kegg prefix from the header)
                    idx = line.find("Release")
                    version_str = line[idx:].strip()
                    return {
                        "remote_version": version_str,
                        "last_modified": "",
                    }
            dates = re.findall(r"\b\d{4}/\d{2}/\d{2}\b", resp.text)
            if dates:
                latest = max(dates, key=lambda value: datetime.strptime(value, "%Y/%m/%d"))
                return {
                    "remote_version": f"KEGG data {latest}",
                    "last_modified": latest,
                }
        except Exception as e:
            logger.warning("Could not close temporary folder: %s%s", e)
        return None

    def check_reactome_version(self) -> Optional[Dict[str, str]]:
        """Inspect the current upstream Reactome release."""
        url = "https://reactome.org/ContentService/data/database/version"
        try:
            resp = requests.get(url, timeout=self.timeout, headers={
                "Accept": "text/plain",
                "User-Agent": "AllEnricher/2.0 (database version checker)"
            })
            resp.raise_for_status()
            version_num = resp.text.strip()
            if version_num and version_num.isdigit():
                return {
                    "remote_version": f"v{version_num}",
                    "last_modified": "",  # API does not return last-modified
                }
        except Exception as e:
            logger.warning("Could not close temporary folder: %s%s", e)
        return None

    def check_all_sources(self) -> Dict[str, Dict[str, str]]:
        """Check all supported upstream sources for newer versions."""
        results: Dict[str, Dict[str, str]] = {}

        # HTTP HEAD class (gene2go, gene_info, taxonomy)
        for source in ["gene2go", "gene_info", "taxonomy"]:
            info = self.check_head(self.SOURCE_URLS[source])
            if info:
                results[source] = info

        # GO OBO (file content interpretation)
        go_info = self.check_go_obo_version()
        if go_info:
            results["go_obo"] = go_info

        # GOA Proteomes (Dialer HEAD)
        goa_info = self.check_head(self.SOURCE_URLS["goa_proteomes"])
        if goa_info:
            results["goa_proteomes"] = goa_info

        # KEG (API Query)
        kegg_info = self.check_kegg_version()
        if kegg_info:
            results["kegg"] = kegg_info

        # Reactome (page parsing)
        reactome_info = self.check_reactome_version()
        if reactome_info:
            results["reactome"] = reactome_info

        return results

    def check_updates(self, local_manager: DatabaseVersionManager) -> Dict[str, Dict]:
        """Compare remote source metadata with the local manifest."""
        remote_versions = self.check_all_sources()
        update_status: Dict[str, Dict] = {}

        for source, remote_info in remote_versions.items():
            # Try multiple key aliases for local versions
            local_ver = None
            for alias in self._SOURCE_KEY_ALIASES.get(source, [source]):
                local_ver = local_manager.get_local_version(alias)
                if local_ver is not None:
                    break
            has_update = False

            if local_ver is None:
                has_update = True  # Locally never downloaded
            elif "last_modified" in remote_info and local_ver.remote_last_modified:
                try:
                    remote_dt = parsedate_to_datetime(remote_info["last_modified"])
                    local_dt = parsedate_to_datetime(local_ver.remote_last_modified)
                    has_update = remote_dt > local_dt
                except Exception:
                    has_update = False
            elif "remote_version" in remote_info and local_ver.remote_version:
                has_update = remote_info["remote_version"] != local_ver.remote_version
            elif remote_info.get("last_modified") or remote_info.get("remote_version"):
                has_update = True  # The old version of the record lacks a remote benchmark and cannot prove to be the latest

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
