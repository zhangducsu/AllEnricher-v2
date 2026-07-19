"""Download shared source data and maintain database-specific species registries."""

from __future__ import annotations

import csv
import gzip
import logging
import re
import shutil
import sys
import time
from pathlib import Path
from datetime import datetime
from typing import Dict, Iterable, List, Optional, Tuple

import requests

from .download_manager import DownloadManager
from .mirrors import get_mirrors, JENSEN_SOURCES
from .species_registry import SpeciesRegistry, SpeciesEntry

logger = logging.getLogger(__name__)


class DataDownloader:
    """Coordinate source downloads and species coverage registry generation."""

    def __init__(
        self,
        root_dir: str = "./database",
        overwrite: bool = False,
        max_workers: int = 4,
        use_multi_thread: bool = True,
        verify_integrity: bool = True,
    ):
        """Initialise Downloader

        Args:
root_dir: database root
overwrite: Whether to overwrite existing files
Max_workers: number of threads downloaded over multiple threads
use_multi_thread: Whether to enable multiple threads to download large files
vereffy_integrity: Verify the integrity of the download file
        """
        self.root_dir = Path(root_dir)
        self.basic_dir = self.root_dir / "basic"
        self.basic_dir.mkdir(parents=True, exist_ok=True)

        self.manager = DownloadManager(
            root_dir=root_dir,
            overwrite=overwrite,
            max_workers=max_workers,
            use_multi_thread=use_multi_thread,
            verify_integrity=verify_integrity,
            show_progress=True,
        )

    _DATABASE_REGISTRY_PATHS = {
        "disgenet": "disgenet/disgenet_species_registry.tsv",
        "trrust": "trrust/trrust_species_registry.tsv",
        "chea3": "chea3/chea3_species_registry.tsv",
        "animaltfdb": "animaltfdb/animaltfdb_species_registry.tsv",
        "htftarget": "htftarget/htftarget_species_registry.tsv",
    }

    def record_database_species(
        self,
        database: str,
        species: Iterable[Tuple[int, str]],
    ) -> Path:
        """Write TaxID-keyed coverage for one downloaded database."""
        key = database.strip().lower()
        relative_path = self._DATABASE_REGISTRY_PATHS.get(key)
        if relative_path is None:
            raise ValueError(f"The establishment of a database of species registers is not supported: {database}")

        output_path = self.basic_dir / relative_path
        output_path.parent.mkdir(parents=True, exist_ok=True)
        rows = sorted({int(taxid): str(name).strip() for taxid, name in species}.items())
        if not rows or any(not name for _, name in rows):
            raise ValueError(f"{database} species register is empty or contains empty species names")

        with output_path.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t")
            writer.writerow(["taxid", "latin_name"])
            writer.writerows(rows)
        logger.info("%sSpecies register has been generated: %s (%d(species)", database, output_path, len(rows))
        return output_path

    def _latest_registry(self, pattern: str) -> Optional[Path]:
        matches = sorted(self.basic_dir.glob(pattern))
        return matches[-1] if matches else None

    def refresh_supported_species_registry(self) -> Path:
        """Merge database coverage files into the unified species registry."""
        missing = Path("/__allenricher_missing_registry__")
        output_path = self.basic_dir / "supported_species.tsv"
        return self._merge_all_registries(
            go_registry=self._latest_registry("go/GO*/go_species_registry.tsv") or missing,
            kegg_registry=self._latest_registry("kegg/kegg_species_registry.tsv") or missing,
            reactome_registry=self._latest_registry(
                "reactome/Reactome*/reactome_species_registry.tsv"
            ) or missing,
            do_registry=self._latest_registry("do/do_species_registry.tsv") or missing,
            wikipathways_registry=self._latest_registry(
                "wikipathways_species_registry.tsv"
            ) or missing,
            output_path=output_path,
            disgenet_registry=self._latest_registry("disgenet/disgenet_species_registry.tsv"),
            trrust_registry=self._latest_registry("trrust/trrust_species_registry.tsv"),
            chea3_registry=self._latest_registry("chea3/chea3_species_registry.tsv"),
            animaltfdb_registry=self._latest_registry(
                "animaltfdb/animaltfdb_species_registry.tsv"
            ),
            htftarget_registry=self._latest_registry("htftarget/htftarget_species_registry.tsv"),
        )

    # ============================
    # GO Basic Data Download (All Species)
    # ============================
    def download_go_basic(self, version: Optional[str] = None) -> str:
        """Download shared Gene Ontology and NCBI annotation inputs."""
        if version is None:
            version = f"GO{datetime.now().strftime('%Y%m%d')}"

        go_dir = self.basic_dir / "go" / version
        go_dir.mkdir(parents=True, exist_ok=True)

        print(f"\n{'='*60}")
        print(f"Download GO Basic data -> {go_dir}")
        print(f"{'='*60}")

        ncbi_mirrors = get_mirrors('ncbi')

        self.manager.download_with_mirror_fallback(
            ncbi_mirrors, "gene2go.gz",
            go_dir / "gene2go.gz", desc="gene2go.gz"
        )
        self.manager.download_with_mirror_fallback(
            ncbi_mirrors, "gene_info.gz",
            go_dir / "gene_info.gz", desc="gene_info.gz"
        )

        go_mirrors = get_mirrors('go')
        self.manager.download_with_mirror_fallback(
            go_mirrors, "go-basic.obo",
            go_dir / "go-basic.obo", desc="go-basic.obo"
        )

        print(f"GO Base data download complete -> {go_dir}")

        # Recording version metadata to versions. json
        try:
            from allenricher.database.version import DatabaseVersionManager, RemoteVersionChecker
            _vm = DatabaseVersionManager(database_dir=str(self.root_dir))
            _checker = RemoteVersionChecker()

            # Record gene2go
            _g2g_info = _checker.check_head("https://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2go.gz")
            if _g2g_info:
                _vm.record_download(
                    source="gene2go", local_version=version,
                    local_path=f"basic/go/{version}",
                    remote_last_modified=_g2g_info.get("last_modified"),
                )

            # Record
            _obo_info = _checker.check_go_obo_version()
            if _obo_info:
                _vm.record_download(
                    source="go_obo", local_version=version,
                    local_path=f"basic/go/{version}/go-basic.obo",
                    remote_version=_obo_info.get("remote_version"),
                    remote_last_modified=_obo_info.get("last_modified"),
                )

            # Record Gene_info
            _gi_info = _checker.check_head("https://ftp.ncbi.nlm.nih.gov/gene/DATA/gene_info.gz")
            if _gi_info:
                _vm.record_download(
                    source="gene_info", local_version=version,
                    local_path=f"basic/go/{version}",
                    remote_last_modified=_gi_info.get("last_modified"),
                )

            _vm.record_download(
                source="go", local_version=version,
                local_path=f"basic/go/{version}",
                remote_version=(_obo_info or {}).get("remote_version"),
                remote_last_modified=(_obo_info or {}).get("last_modified"),
            )
        except Exception as _e:
            logger.warning("Recording of metadata from the GO version failed: %s", _e)

        return str(go_dir)

    # ============================
    # Reactome Basic Data Download (All Species)
    # ============================
    def download_reactome_basic(self, version: Optional[str] = None,
                                  go_version: Optional[str] = None) -> str:
        """Download shared Reactome mapping and hierarchy inputs."""
        if version is None:
            version = f"Reactome{datetime.now().strftime('%Y%m%d')}"

        re_dir = self.basic_dir / "reactome" / version
        re_dir.mkdir(parents=True, exist_ok=True)

        print(f"\n{'='*60}")
        print(f"Download Reactome Basic data -> {re_dir}")
        print(f"{'='*60}")

        # Gene_info.gz: priority for re-do downloaded files
        re_gene_info = re_dir / "gene_info.gz"
        if re_gene_info.exists() and not self.manager.overwrite:
            if self.manager.verify_integrity:
                from .download_utils import verify_gzip_integrity
                valid, _ = verify_gzip_integrity(re_gene_info)
                if valid:
                    print(f"--- exists and is valid and skips: gene_info.gz")
        else:
            # Try to copy from the Go directory
            if go_version is None:
                go_version = self.get_latest_go_version()
            if go_version:
                go_gene_info = self.basic_dir / "go" / go_version / "gene_info.gz"
                if go_gene_info.exists():
                    print(f"|---Reuse GO data: {go_gene_info}")
                    shutil.copy2(go_gene_info, re_gene_info)
                    print(f"--- has been copied: gene_info.gz ({re_gene_info.stat().st_size / 1024 / 1024: .1f} MB)")
                else:
                    # No GO directory, download from mirror.
                    ncbi_mirrors = get_mirrors('ncbi')
                    self.manager.download_with_mirror_fallback(
                        ncbi_mirrors, "gene_info.gz",
                        re_gene_info, desc="gene_info.gz"
                    )
            else:
                # No Go version, download from mirror
                ncbi_mirrors = get_mirrors('ncbi')
                self.manager.download_with_mirror_fallback(
                    ncbi_mirrors, "gene_info.gz",
                    re_gene_info, desc="gene_info.gz"
                )

        # NCBI2Reactome (gzip compression after download)
        reactome_mirrors = get_mirrors('reactome')
        raw_file = re_dir / "NCBI2Reactome_All_Levels.txt"
        gz_file = re_dir / "NCBI2Reactome_All_Levels.txt.gz"

        ncbi_ready = False
        if gz_file.exists() and not self.manager.overwrite:
            if self.manager.verify_integrity:
                from .download_utils import verify_gzip_integrity
                valid, _ = verify_gzip_integrity(gz_file)
                if valid:
                    print(f"|---Existence and validity; skipping: {gz_file.name}")
                    ncbi_ready = True
            else:
                ncbi_ready = True

        if not ncbi_ready:
            self.manager.download_with_mirror_fallback(
                reactome_mirrors, "NCBI2Reactome_All_Levels.txt",
                raw_file, desc="NCBI2Reactome"
            )

            if raw_file.exists():
                print(f"|--- Compression: {raw_file.name} -> {gz_file.name}")
                with open(raw_file, 'rb') as f_in:
                    with gzip.open(gz_file, 'wb') as f_out:
                        shutil.copyfileobj(f_in, f_out)
                raw_file.unlink()

        for filename in ("ReactomePathways.txt", "ReactomePathwaysRelation.txt"):
            self._download_small_reactome_metadata(
                reactome_mirrors, filename, re_dir / filename
            )

        print(f"Reactome Base data download complete -> {re_dir}")

        self._record_reactome_version(version)

        return str(re_dir)

    def _download_small_reactome_metadata(
        self,
        mirrors,
        filename: str,
        destination: Path,
    ) -> None:
        """Download small Reactome metadata files without range requests."""
        last_error = None
        for mirror in mirrors:
            url = f"{mirror.base_url.rstrip('/')}/{filename}"
            try:
                response = requests.head(url, allow_redirects=True, timeout=30)
                response.raise_for_status()
                expected = int(response.headers.get("Content-Length") or 0)
                if expected == 0:
                    etag_size = response.headers.get("ETag", "").strip('"').split("-", 1)[0]
                    if re.fullmatch(r"[0-9a-fA-F]+", etag_size):
                        expected = int(etag_size, 16)
                if (
                    destination.is_file()
                    and not self.manager.overwrite
                    and expected > 0
                    and destination.stat().st_size == expected
                ):
                    print(f"|---Existence and size correct, skip: {filename}")
                    return
                previous_multi_thread = self.manager.use_multi_thread
                previous_overwrite = self.manager.overwrite
                self.manager.use_multi_thread = False
                self.manager.overwrite = True
                try:
                    self.manager.download_file(
                        url,
                        destination,
                        expected_size=expected or None,
                        desc=filename,
                    )
                finally:
                    self.manager.use_multi_thread = previous_multi_thread
                    self.manager.overwrite = previous_overwrite
                if expected > 0 and destination.stat().st_size != expected:
                    raise RuntimeError(
                        f"Downloaded size mismatch for {filename}: "
                        f"{destination.stat().st_size} != {expected} bytes"
                    )
                return
            except Exception as exc:
                last_error = exc
                print(f"    [FAILED] {mirror.name}: {exc}")
        raise RuntimeError(f"{filename}All Reactome mirrors failed: {last_error}")

    def _record_reactome_version(self, version: str) -> None:
        """Record Reactome source version and retrieval metadata."""
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
            logger.warning("Failed to record the Reactome version metadata: %s", _e)

    # ============================
    # DO / DisGeNET (humans only)
    # ============================
    def download_do_files(self) -> Dict[str, str]:
        """Download the human Disease Ontology source files."""
        do_dir = self.basic_dir / "do"
        do_dir.mkdir(parents=True, exist_ok=True)

        print(f"\n{'='*60}")
        print(f"Download DO Data -> {do_dir}")
        print(f"{'='*60}")

        files = {}
        for url in JENSEN_SOURCES:
            fname = url.split('/')[-1]
            raw_dest = do_dir / fname
            gz_dest = do_dir / f"{fname}.gz"

            # Existing and valid compression file  Skip
            if gz_dest.exists() and not self.manager.overwrite:
                if self.manager.verify_integrity:
                    from .download_utils import verify_gzip_integrity
                    valid, _ = verify_gzip_integrity(gz_dest)
                    if valid:
                        print(f"|---Existence and validity; skipping: {gz_dest.name}")
                        files[fname] = str(gz_dest)
                        continue
                else:
                    print(f"|---Existence, Skipping: {gz_dest.name}")
                    files[fname] = str(gz_dest)
                    continue

            # Download original TSV
            self.manager.download_file(url, raw_dest, desc=fname)

            # gzip compression
            if raw_dest.exists():
                print(f"|--- Compression: {raw_dest.name} -> {gz_dest.name}")
                with open(raw_dest, 'rb') as f_in:
                    with gzip.open(gz_dest, 'wb') as f_out:
                        shutil.copyfileobj(f_in, f_out)
                raw_dest.unlink()
                files[fname] = str(gz_dest)
            else:
                files[fname] = str(raw_dest)

        ontology_path = do_dir / "doid.obo"
        if ontology_path.exists() and not self.manager.overwrite:
            print("|---Existence, Skip: doid.obo")
        else:
            self.manager.download_file(
                "https://purl.obolibrary.org/obo/doid.obo",
                ontology_path,
                desc="Disease Ontology doid.obo",
            )
        files[ontology_path.name] = str(ontology_path)

        # Recording the version of the Do metadata to versions. json
        try:
            from allenricher.database.version import DatabaseVersionManager
            _vm = DatabaseVersionManager(database_dir=str(self.root_dir))
            _vm.record_download(
                source="do",
                local_version="cached",
                local_path="basic/do",
            )
        except Exception as _e:
            logger.warning("Failed to record metadata for version DO: %s", _e)

        return files

    def download_disgenet(self) -> str:
        """Install the retained human DisGeNET snapshot."""
        dg_dir = self.basic_dir / "disgenet"
        dg_dir.mkdir(parents=True, exist_ok=True)

        print(f"\n{'='*60}")
        print(f"Download DisGeNET Data -> {dg_dir}")
        print(f"{'='*60}")

        url = (
            "http://www.disgenet.org/static/disgenet_ap1/files/"
            "downloads/all_gene_disease_associations.tsv.gz"
        )
        dest = dg_dir / "all_gene_disease_associations.tsv.gz"

        legacy_files = sorted(self.root_dir.glob("organism/v*/hsa/hsa.DisGeNET.gmt.gz"))
        if legacy_files:
            print("--- reminder: DisGeNET is still being updated, but the new edition is no longer available for public download free of charge")
            print("|---AllEnricher only repeats the constructed DisGeNET database with v1.")
            print(f"|---Current reuse: {legacy_files[-1]}")
            return str(legacy_files[-1])
        raise RuntimeError(
            "Current DisGeNET releases require authorization and are not available through the free downloader. "
            "No reusable AllEnricher v1 snapshot was found. Copy hsa.CUI2gene.tab.gz and "
            "hsa.CUI2disc.gz from the v1 species database into the v2 species database."
        )

    # ============================
    # WikiPathways Basic Data Download
    # ============================
    def download_wikipathways_basic(self, version: Optional[str] = None) -> str:
        """Download species-specific WikiPathways GMT files."""
        from .wikipathways_fetcher import WikiPathwaysFetcher

        fetcher = WikiPathwaysFetcher(str(self.basic_dir))

        print(f"\n{'='*60}")
        print(f"Download WikiPathways Basic Data")
        print(f"{'='*60}")

        if version is None:
            version = fetcher._detect_latest_version()

        results = fetcher.download_all_gmt(version=version, overwrite=self.manager.overwrite)

        version_str = version or "unknown"
        version_dir = self.basic_dir / "wikipathways" / f"WP{version_str}"

        # Recording version metadata to versions. json
        try:
            from allenricher.database.version import DatabaseVersionManager, RemoteVersionChecker
            _vm = DatabaseVersionManager(database_dir=str(self.root_dir))
            _checker = RemoteVersionChecker()
            _wp_info = _checker.check_head("https://data.wikipathways.org/")
            if _wp_info:
                _vm.record_download(
                    source="wikipathways",
                    local_version=f"WP{version_str}",
                    local_path=f"basic/wikipathways/WP{version_str}",
                    remote_last_modified=_wp_info.get("last_modified"),
                )
        except Exception as _e:
            logger.warning("Failed to record metadata from WikiPathways version: %s", _e)

        print(f"WikiPathways Base data download complete -> {version_dir}")
        print(f"|---Downloaded.{len(results)}GMT file for individual species")

        return str(version_dir)

    def _build_wikipathways_registry(self) -> Optional[Path]:
        """Build WikiPathways species coverage from downloaded metadata."""
        from .wikipathways_fetcher import WikiPathwaysFetcher, SPECIES_NAME_MAP

        logger.info("Build WikiPathways species register...")

        fetcher = WikiPathwaysFetcher(str(self.basic_dir))
        registry_path = self.basic_dir / "wikipathways_species_registry.tsv"

        # Get the list of species in the latest version
        try:
            species_list = fetcher.get_available_species()
        except Exception as e:
            logger.warning(f"Failed to obtain list of WikiPathways species: {e}")
            return None

        if not species_list:
            logger.warning("WikiPathways species list is empty")
            return None

        with open(registry_path, 'w', encoding='utf-8', newline='') as f:
            writer = csv.writer(f, delimiter='\t')
            writer.writerow(["species_latin_name", "species_code", "gene_count", "pathway_count"])
            for latin_name in sorted(species_list):
                code = fetcher.get_species_code(latin_name)
                code_str = code or "-"
                # gene_count and path_account fill in during build phase and leave empty during download phase
                writer.writerow([latin_name, code_str, "", ""])

        logger.info(
            "WikiPathways Species Register Generated: %s (%d Specimon)",
            registry_path, len(species_list)
        )
        return registry_path

    # ============================
    # Batch Downloads
    # ============================
    def download_all(self, db_types: List[str] = None) -> Dict[str, str]:
        """Download TRRUST data for all supported species."""
        if db_types is None:
            db_types = ['go', 'reactome']

        result = {}
        go_version = None  # Record GO version for Reactome reuse

        errors = []
        for db_type in db_types:
            db_type = db_type.lower().strip()
            try:
                if db_type in ('go',):
                    result['go'] = self.download_go_basic()
                    go_version = self.get_latest_go_version()  # Record version
                elif db_type in ('reactome',):
                    result['reactome'] = self.download_reactome_basic(go_version=go_version)
                elif db_type in ('do',):
                    result['do'] = str(self.basic_dir / "do")
                    self.download_do_files()
                elif db_type in ('kegg',):
                    # KEGG fetched via REST API, needing gene_info.gz
                    result['kegg'] = str(self.basic_dir / "kegg")
                    print("|---KEG data will be automatically retrieved through REST API during the build")
                elif db_type in ('taxonomy',):
                    taxonomy_dir = self.basic_dir / "taxonomy"
                    taxonomy_path = self._download_taxonomy_names(taxonomy_dir)
                    if taxonomy_path is None:
                        raise RuntimeError("NCBI Taxonomy Download or Parsing Failed")
                    result['taxonomy'] = str(taxonomy_dir)
                elif db_type in ('disgenet',):
                    result['disgenet'] = self.download_disgenet()
                elif db_type in ('wikipathways',):
                    result['wikipathways'] = self.download_wikipathways_basic()
                else:
                    print(f"|---[Warning] Unknown database type: {db_type}")
            except Exception as e:
                print(f"|---[Error] Failed to download {db_type}: {e}")
                errors.append(f"{db_type}: {e}")

        if errors:
            raise RuntimeError("; ".join(errors))

        registry_sources = {'go', 'kegg', 'reactome', 'do', 'disgenet', 'wikipathways'}
        if not registry_sources.intersection(result):
            return result

        # ============================
        # Form to build a current
        # ============================
        print(f"\n{'='*60}")
        print("Build species register")
        print(f"{'='*60}")

        go_registry: Optional[Path] = None
        goa_index: Optional[Path] = None
        kegg_registry: Optional[Path] = None
        reactome_registry: Optional[Path] = None
        do_registry: Optional[Path] = None

        # Build the GO species registry when GO data were downloaded.
        if 'go' in result:
            go_dir = Path(result['go'])
            gene2go_path = go_dir / "gene2go.gz"

            if gene2go_path.exists():
                try:
                    go_registry = self._build_go_registry(gene2go_path, go_dir)
                    logger.info(f"GO registry built: {go_registry}")
                except Exception as e:
                    logger.warning(f"Failed to build GO registry: {e}")

                try:
                    goa_index = self._download_goa_index(go_dir)
                    logger.info(f"GOA index downloaded: {goa_index}")
                except Exception as e:
                    logger.warning(f"Failed to download GOA index: {e}")

                if go_registry and goa_index:
                    try:
                        go_registry = self._merge_go_registries(go_registry, goa_index, go_dir)
                        logger.info(f"GO registries merged: {go_registry}")
                    except Exception as e:
                        logger.warning(f"Failed to merge GO registries: {e}")

        # Build the KEGG species registry when KEGG data were downloaded.
        if 'kegg' in result:
            kegg_dir = Path(result['kegg'])
            try:
                kegg_registry = self._build_kegg_registry(kegg_dir)
                logger.info(f"KEGG registry built: {kegg_registry}")
                from allenricher.database.version import DatabaseVersionManager, RemoteVersionChecker
                kegg_info = RemoteVersionChecker().check_kegg_version() or {}
                DatabaseVersionManager(database_dir=str(self.root_dir)).record_download(
                    source="kegg",
                    local_version="cached",
                    local_path="basic/kegg",
                    remote_version=kegg_info.get("remote_version"),
                    remote_last_modified=kegg_info.get("last_modified"),
                )
            except Exception as e:
                raise RuntimeError(f"Failed to build KEGG registry: {e}") from e

        # Record Reactome species coverage when Reactome data were downloaded.
        if 'reactome' in result:
            reactome_dir = Path(result['reactome'])
            ncbi2reactome_path = reactome_dir / "NCBI2Reactome_All_Levels.txt.gz"

            if ncbi2reactome_path.exists():
                try:
                    reactome_registry = self._build_reactome_registry(ncbi2reactome_path, reactome_dir)
                    logger.info(f"Reactome registry built: {reactome_registry}")
                except Exception as e:
                    logger.warning(f"Failed to build Reactome registry: {e}")

        # Add DO support to the species registry when DO was downloaded.
        if 'do' in result:
            do_dir = Path(result['do'])
            try:
                do_registry = self._build_do_registry(do_dir)
                logger.info(f"DO registry built: {do_registry}")
            except Exception as e:
                logger.warning(f"Failed to build DO registry: {e}")

        # Record WikiPathways species coverage when data were downloaded.
        wikipathways_registry: Optional[Path] = None
        if 'wikipathways' in result:
            try:
                wikipathways_registry = self._build_wikipathways_registry()
                logger.info(f"WikiPathways registry built: {wikipathways_registry}")
            except Exception as e:
                logger.warning(f"Failed to build WikiPathways registry: {e}")

        if 'disgenet' in result:
            self.record_database_species("DisGeNET", [(9606, "Homo sapiens")])

        # Merge database-specific species registries while preserving entries
        # for databases that were not part of this download operation.
        try:
            supported_species_path = self.refresh_supported_species_registry()
            logger.info(f"All registries merged: {supported_species_path}")
            result['supported_species'] = str(supported_species_path)
        except Exception as e:
            raise RuntimeError(f"Failed to merge all registries: {e}") from e

        # Print download statistical summary
        try:
            self._report_download_summary(supported_species_path)
        except Exception as e:
            logger.warning(f"Failed to report download summary: {e}")

        return result

    # ============================
    # Version Management
    # ============================
    def list_go_versions(self) -> list:
        """Return locally downloaded Gene Ontology versions."""
        go_basic = self.basic_dir / "go"
        if not go_basic.exists():
            return []
        return sorted([d.name for d in go_basic.iterdir() if d.is_dir()])

    def list_reactome_versions(self) -> list:
        """Return locally downloaded Reactome versions."""
        re_basic = self.basic_dir / "reactome"
        if not re_basic.exists():
            return []
        return sorted([d.name for d in re_basic.iterdir() if d.is_dir()])

    def get_latest_go_version(self) -> Optional[str]:
        """Return the newest local Gene Ontology version."""
        versions = self.list_go_versions()
        return versions[-1] if versions else None

    def get_latest_reactome_version(self) -> Optional[str]:
        """Return the newest local Reactome version."""
        versions = self.list_reactome_versions()
        return versions[-1] if versions else None

    def _download_taxonomy_names(self, output_dir: Path) -> Optional[Path]:
        """Download NCBI taxonomy names used by the species registry."""
        logger.info("Download NCBI Taxony species names...")
        taxonomy_dir = self.basic_dir / "taxonomy"
        taxonomy_dir.mkdir(parents=True, exist_ok=True)

        # Check if names.dmp
        names_dmp = taxonomy_dir / "names.dmp"
        if names_dmp.exists() and not self.manager.overwrite:
            logger.info("NCBI Taxonomy exists and skips download")
        else:
            # Downloads
            taxdump_url = "https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz"
            taxdump_tgz = taxonomy_dir / "taxdump.tar.gz"

            logger.info("Download data from NCBI from taxonomy: %s", taxdump_url)
            try:
                response = requests.get(taxdump_url, stream=True, timeout=300)
                response.raise_for_status()

                with open(taxdump_tgz, 'wb') as f:
                    for chunk in response.iter_content(chunk_size=8192):
                        f.write(chunk)
                logger.info("Download completed: %s", taxdump_tgz)

                # Unlock name. dmp
                import tarfile
                with tarfile.open(taxdump_tgz, 'r:gz') as tar:
                    names_member = tar.getmember('names.dmp')
                    tar.extract(names_member, taxonomy_dir)
                    logger.info("Unlock names.dmp succeeded")
            except Exception as e:
                logger.warning("Failed to download NCBI Taxonomy: %s", e)
                return None

        if not names_dmp.exists():
            return None

        # Parsing name names.dmp to generate atxid & latin_name map
        # Names.dmp format: tax_id\t name \t unique name \t name class
        taxid_to_name: Dict[int, str] = {}
        seen_taxids: set = set()

        with open(names_dmp, 'r', encoding='utf-8', errors='ignore') as f:
            for line in f:
                line = line.rstrip('\t|\n')
                if not line:
                    continue
                parts = line.split('\t|\t')
                if len(parts) < 4:
                    continue
                try:
                    taxid = int(parts[0])
                    name_class = parts[3].rstrip('\t|')

                    # Only scientifier name
                    if name_class == 'scientific name':
                        taxid_to_name[taxid] = parts[1]
                except (ValueError, IndexError):
                    continue

        # Save as TSV
        output_path = output_dir / "taxid_to_name.tsv"
        with open(output_path, 'w', encoding='utf-8', newline='') as f:
            writer = csv.writer(f, delimiter='\t')
            writer.writerow(['taxid', 'latin_name'])
            for taxid in sorted(taxid_to_name.keys()):
                writer.writerow([taxid, taxid_to_name[taxid]])

        logger.info("NCBI Taxonomy species name has been created: %s (%d(species)", output_path, len(taxid_to_name))

        # Recording version metadata to versions. json
        try:
            from allenricher.database.version import DatabaseVersionManager, RemoteVersionChecker
            _vm = DatabaseVersionManager(database_dir=str(self.root_dir))
            _checker = RemoteVersionChecker()
            _tax_info = _checker.check_head("https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz")
            if _tax_info:
                _vm.record_download(
                    source="taxonomy", local_version="cached",
                    local_path="basic/taxonomy",
                    remote_last_modified=_tax_info.get("last_modified"),
                )
        except Exception as _e:
            logger.warning("Failed to record metadata for Taxonomy version: %s", _e)

        return output_path

    def _load_taxid_to_name_map(self, taxonomy_dir: Path) -> Dict[int, str]:
        """Load the local NCBI TaxID-to-scientific-name mapping."""
        taxid_map_path = taxonomy_dir / "taxid_to_name.tsv"
        if not taxid_map_path.exists():
            return {}

        result: Dict[int, str] = {}
        with open(taxid_map_path, 'r', encoding='utf-8') as f:
            reader = csv.DictReader(f, delimiter='\t')
            for row in reader:
                try:
                    result[int(row['taxid'])] = row['latin_name']
                except (ValueError, KeyError):
                    continue
        return result

    # ============================
    # Species Register Construction
    # ============================

    def _build_go_registry(
        self, gene2go_path: Path, output_dir: Path
    ) -> Path:
        """Derive Gene Ontology species coverage from NCBI gene2go."""
        logger.info("Build GO species register: gene2go_path=%s", gene2go_path)

        # ----1. Remove all the only taxid and statistics from gene2go.gz----
        taxid_stats: Dict[int, Dict[str, int]] = {}  # taxid -> {gene_count, term_count}
        with gzip.open(gene2go_path, "rt", encoding="utf-8") as fh:
            for line in fh:
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                parts = line.split("\t")
                if len(parts) < 6:
                    continue
                try:
                    taxid = int(parts[0])
                    gene_id = parts[1]
                    go_term = parts[2]
                except (ValueError, IndexError):
                    continue
                if taxid not in taxid_stats:
                    taxid_stats[taxid] = {"gene_count": set(), "term_count": set()}
                taxid_stats[taxid]["gene_count"].add(gene_id)
                taxid_stats[taxid]["term_count"].add(go_term)

        # ----2. Fetch taxid  Latin_name map from NCBI Taxonomy----
        # Priority for using the sciencefic name in the taxony database
        taxonomy_dir = self.basic_dir / "taxonomy"
        taxid_to_name = self._load_taxid_to_name_map(taxonomy_dir)

        # If the taxony map is empty, download first
        if not taxid_to_name:
            logger.info("NCBI Taxony map is empty, tries to download...")
            taxonomy_tsv = self._download_taxonomy_names(taxonomy_dir)
            if taxonomy_tsv:
                taxid_to_name = self._load_taxid_to_name_map(taxonomy_dir)

        # Back to gene_info if still empty (for cases where taxonomy is not available only)
        # Note: Gene_info does not contain Latin names of species only for the purpose of recording taxid presence
        if not taxid_to_name:
            logger.warning("No scientific names could be resolved from NCBI Taxonomy or gene_info")
            logger.warning("The GO species registry will contain blank scientific names; review the source data")

        # Step 3: Write go_species_registry.tsv.
        output_dir.mkdir(parents=True, exist_ok=True)
        output_path = output_dir / "go_species_registry.tsv"

        with open(output_path, "w", encoding="utf-8", newline="") as fh:
            writer = csv.writer(fh, delimiter="\t")
            writer.writerow(["taxid", "latin_name", "source", "gene_count", "term_count"])
            for taxid in sorted(taxid_stats):
                latin_name = taxid_to_name.get(taxid, "")
                stats = taxid_stats[taxid]
                writer.writerow([
                    taxid,
                    latin_name,
                    "ncbi_gene2go",
                    len(stats["gene_count"]),
                    len(stats["term_count"]),
                ])

        logger.info(
            "GO species register has been generated: %s (%d(species)", output_path, len(taxid_stats)
        )
        return output_path

    def _download_goa_index(self, output_dir: Path) -> Path:
        """Retrieve the UniProt GOA species file index."""
        logger.info("Download GOA species index...")

        url = "https://ftp.ebi.ac.uk/pub/databases/GO/goa/proteomes/"
        headers = {
            "User-Agent": (
                "AllEnricher/2.0 (https://github.com/allenricher; "
                "data download pipeline)"
            ),
        }

        resp = requests.get(url, headers=headers, timeout=30)
        resp.raise_for_status()

        # Parsing HTML, extracting.goa file links
        goa_entries: List[Dict[str, str]] = []
        taxid_list: List[int] = []  # Collect all taxids for follow-up queries
        # Match href="xx.goa" or href="xx.goa.gz"
        for match in re.finditer(r'href="([^"]+\.goa(?:\.gz)?)"', resp.text):
            filename = match.group(1)
            # Remove.gz suffix for uniform treatment
            base_name = filename.removesuffix(".gz")
            # Filename format: {taxid}.{species_name}.goa
            # For example: 9606. Homo_sapiens.goa
            if not base_name.endswith(".goa"):
                continue
            stem = base_name[:-4]  # Get rid of it. Goa.
            dot_pos = stem.find(".")
            if dot_pos < 0:
                continue
            taxid_str = stem[:dot_pos]
            species_part = stem[dot_pos + 1:]
            # Taxid must be a pure number.
            if not taxid_str.isdigit():
                continue
            taxid = int(taxid_str)
            taxid_list.append(taxid)
            # Use the name of the file name for the time being. Next time you will get it from NCBI Taxonomy
            latin_name = species_part.replace("_", " ")
            goa_entries.append({
                "taxid": taxid,
                "latin_name": latin_name,
                "filename": filename,
            })

        # Fetch Latin_name from NCBI Taxony
        taxonomy_dir = self.basic_dir / "taxonomy"
        taxid_to_name = self._load_taxid_to_name_map(taxonomy_dir)

        # Update Latin_name: use the name in NCBI Taxonomy as a priority
        for entry in goa_entries:
            taxid = entry["taxid"]
            if taxid in taxid_to_name:
                entry["latin_name"] = taxid_to_name[taxid]

        # Write goa_species_index.tsv
        output_dir.mkdir(parents=True, exist_ok=True)
        output_path = output_dir / "goa_species_index.tsv"

        with open(output_path, "w", encoding="utf-8", newline="") as fh:
            writer = csv.writer(fh, delimiter="\t")
            writer.writerow(["taxid", "latin_name", "filename"])
            for entry in sorted(goa_entries, key=lambda e: e["taxid"]):
                writer.writerow([
                    entry["taxid"],
                    entry["latin_name"],
                    entry["filename"],
                ])

        logger.info(
            "GOA species index has been generated: %s (%d(species)", output_path, len(goa_entries)
        )
        return output_path

    def _merge_go_registries(
        self,
        gene2go_registry: Path,
        goa_index: Path,
        output_dir: Path,
    ) -> Path:
        """Merge NCBI gene2go and UniProt GOA coverage by TaxID."""
        logger.info("Merging GO species registries...")

        # Read species coverage derived from NCBI gene2go.
        g2g_data: Dict[int, Dict[str, str]] = {}
        if gene2go_registry.exists():
            with open(gene2go_registry, "r", encoding="utf-8") as fh:
                reader = csv.DictReader(fh, delimiter="\t")
                for row in reader:
                    taxid = int(row["taxid"])
                    g2g_data[taxid] = {
                        "latin_name": row.get("latin_name", ""),
                        "source": row.get("source", "ncbi_gene2go"),
                        "gene_count": row.get("gene_count", ""),
                        "term_count": row.get("term_count", ""),
                    }

        # Read GOA Index
        goa_data: Dict[int, Dict[str, str]] = {}
        if goa_index.exists():
            with open(goa_index, "r", encoding="utf-8") as fh:
                reader = csv.DictReader(fh, delimiter="\t")
                for row in reader:
                    taxid = int(row["taxid"])
                    goa_data[taxid] = {
                        "latin_name": row.get("latin_name", ""),
                        "filename": row.get("filename", ""),
                    }

        # Merge - prioritize the use of more reliable species names
        all_taxids = sorted(set(g2g_data) | set(goa_data))
        merged: List[Dict[str, str]] = []
        for taxid in all_taxids:
            in_g2g = taxid in g2g_data
            in_goa = taxid in goa_data

            if in_g2g and in_goa:
                source = "both"
                # Use ling_name of gene2go preferred (from NCBI Taxonomy)
                # Only use GOA names when the name of the gene2go is empty
                g2g_name = g2g_data[taxid]["latin_name"]
                goa_name = goa_data[taxid]["latin_name"]
                if g2g_name:
                    latin_name = g2g_name
                else:
                    latin_name = goa_name
                gene_count = g2g_data[taxid]["gene_count"]
                term_count = g2g_data[taxid]["term_count"]
            elif in_g2g:
                source = "ncbi_gene2go"
                latin_name = g2g_data[taxid]["latin_name"]
                gene_count = g2g_data[taxid]["gene_count"]
                term_count = g2g_data[taxid]["term_count"]
            else:
                source = "uniprot_goa"
                latin_name = goa_data[taxid]["latin_name"]
                gene_count = ""
                term_count = ""

            merged.append({
                "taxid": str(taxid),
                "latin_name": latin_name,
                "source": source,
                "gene_count": gene_count,
                "term_count": term_count,
            })

        # Write the merged GO species registry.
        output_dir.mkdir(parents=True, exist_ok=True)
        output_path = output_dir / "go_species_registry.tsv"

        with open(output_path, "w", encoding="utf-8", newline="") as fh:
            writer = csv.DictWriter(
                fh, fieldnames=["taxid", "latin_name", "source", "gene_count", "term_count"],
                delimiter="\t",
            )
            writer.writeheader()
            writer.writerows(merged)

        logger.info(
            "GO species registry written: %s (%d species; gene2go=%d, GOA=%d, both=%d)",
            output_path,
            len(merged),
            sum(1 for m in merged if m["source"] in ("ncbi_gene2go", "both")),
            sum(1 for m in merged if m["source"] in ("uniprot_goa", "both")),
            sum(1 for m in merged if m["source"] == "both"),
        )
        return output_path

    def _build_kegg_registry(self, output_dir: Path) -> Path:
        """Build KEGG species coverage from the KEGG organism API."""
        logger.info("Build KEGG species register...")

        url = "https://rest.kegg.jp/list/genome"
        headers = {
            "User-Agent": (
                "AllEnricher/2.0 (https://github.com/allenricher; "
                "data download pipeline)"
            ),
        }
        output_dir.mkdir(parents=True, exist_ok=True)
        output_path = output_dir / "kegg_species_registry.tsv"
        try:
            text = self._api_get_with_retry(
                url, headers=headers, timeout=30, max_retries=3
            )
        except RuntimeError:
            if output_path.is_file() and output_path.stat().st_size > 0:
                logger.warning("The KEGG API is unavailable; using the existing species registry: %s", output_path)
                return output_path
            raise

        taxonomy = self._load_taxid_to_name_map(self.basic_dir / "taxonomy")
        if not taxonomy:
            raise RuntimeError("Lack of NCBI Taxonomy name map to match taxid for KEG species")
        name_to_taxid = {name.casefold(): taxid for taxid, name in taxonomy.items()}

        # Current interface format: gender_id\tkegg_code; definition
        # Compatibility with the old four-column formatting.
        entries: List[Dict[str, str]] = []
        for line in text.strip().split("\n"):
            if not line.strip():
                continue
            parts = line.split("\t")
            if len(parts) >= 4:
                kegg_code = parts[1].strip()
                definition = parts[2].strip()
            elif len(parts) >= 2 and ";" in parts[1]:
                kegg_code, definition = (value.strip() for value in parts[1].split(";", 1))
            else:
                continue
            latin_name = re.sub(r"\s+\([^()]*\)\s*$", "", definition).strip()
            taxid = name_to_taxid.get(latin_name.casefold())
            if taxid is None:
                continue

            entries.append({
                "taxid": str(taxid),
                "latin_name": latin_name,
                "kegg_code": kegg_code,
                "kegg_code_source": "kegg",
                "gene_count": "",
            })

        # Writing kegg_species_registry.tsv
        with open(output_path, "w", encoding="utf-8", newline="") as fh:
            writer = csv.DictWriter(
                fh,
                fieldnames=["taxid", "latin_name", "kegg_code", "kegg_code_source", "gene_count"],
                delimiter="\t",
            )
            writer.writeheader()
            writer.writerows(entries)

        logger.info(
            "KEGG species register has been generated: %s (%d species)", output_path, len(entries)
        )
        return output_path

    def _build_reactome_registry(
        self, ncbi2reactome_path: Path, output_dir: Path
    ) -> Path:
        """Build Reactome species coverage from its NCBI mapping file."""
        logger.info("Construct Reactome species register...")

        # Reactome species code * Taxid internal map
        REACTOME_CODE_TO_TAXID: Dict[str, int] = {
            "HSA": 9606, "MMU": 10090, "RNO": 10116, "CEL": 6239,
            "DME": 7227, "SCE": 4932, "ATH": 3702, "DDI": 44689,
            "GGA": 9031, "SSC": 9823, "BTA": 9913, "XTR": 8364,
            "CFA": 9615, "DRE": 7955, "PFA": 5833, "SPO": 4896,
            "MTU": 1772,
        }

        # Reverse Map:
        taxid_to_code: Dict[int, str] = {}

        # Remove species code from file path_id
        pathway_pattern = re.compile(r"R-([A-Z]{3})-\d+")
        with gzip.open(ncbi2reactome_path, "rt", encoding="utf-8") as fh:
            for line in fh:
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                parts = line.split("\t")
                if len(parts) < 2:
                    continue
                pathway_id = parts[1]
                m = pathway_pattern.match(pathway_id)
                if m:
                    code = m.group(1)
                    if code in REACTOME_CODE_TO_TAXID:
                        taxid = REACTOME_CODE_TO_TAXID[code]
                        taxid_to_code[taxid] = code

        # Write_registry.tsv
        output_dir.mkdir(parents=True, exist_ok=True)
        output_path = output_dir / "reactome_species_registry.tsv"

        with open(output_path, "w", encoding="utf-8", newline="") as fh:
            writer = csv.writer(fh, delimiter="\t")
            writer.writerow(["taxid", "reactome_code"])
            for taxid in sorted(taxid_to_code):
                writer.writerow([taxid, taxid_to_code[taxid]])

        logger.info(
            "The Reactome species register has been created: %s (%d(species)",
            output_path, len(taxid_to_code),
        )
        return output_path

    def _build_do_registry(self, output_dir: Path) -> Path:
        """Record the fixed human coverage of Disease Ontology."""
        logger.info("Build DO species register...")

        output_dir.mkdir(parents=True, exist_ok=True)
        output_path = output_dir / "do_species_registry.tsv"

        with open(output_path, "w", encoding="utf-8", newline="") as fh:
            writer = csv.writer(fh, delimiter="\t")
            writer.writerow(["taxid", "latin_name"])
            writer.writerow([9606, "Homo sapiens"])

        logger.info("DO species register has been generated: %s", output_path)
        return output_path

    def _merge_all_registries(
        self,
        go_registry: Path,
        kegg_registry: Path,
        reactome_registry: Path,
        do_registry: Path,
        wikipathways_registry: Path,
        output_path: Path,
        disgenet_registry: Optional[Path] = None,
        trrust_registry: Optional[Path] = None,
        chea3_registry: Optional[Path] = None,
        animaltfdb_registry: Optional[Path] = None,
        htftarget_registry: Optional[Path] = None,
    ) -> Path:
        """Merge all database coverage records into supported_species.tsv."""
        logger.info("Merge all species registers...")

        # Read GO species coverage.
        go_data: Dict[int, Dict[str, str]] = {}
        if go_registry.exists():
            with open(go_registry, "r", encoding="utf-8") as fh:
                reader = csv.DictReader(fh, delimiter="\t")
                for row in reader:
                    taxid = int(row["taxid"])
                    go_data[taxid] = {
                        "latin_name": row.get("latin_name", ""),
                        "source": row.get("source", ""),
                        "gene_count": row.get("gene_count", ""),
                        "term_count": row.get("term_count", ""),
                    }

        # Read KEGG species coverage.
        kegg_data: Dict[int, Dict[str, str]] = {}
        if kegg_registry and kegg_registry.exists():
            with open(kegg_registry, "r", encoding="utf-8") as fh:
                reader = csv.DictReader(fh, delimiter="\t")
                for row in reader:
                    taxid = int(row["taxid"])
                    kegg_data[taxid] = {
                        "latin_name": row.get("latin_name", ""),
                        "kegg_code": row.get("kegg_code", ""),
                        "kegg_code_source": row.get("kegg_code_source", "kegg"),
                        "gene_count": row.get("gene_count", ""),
                    }
        latin_name_to_taxid: Dict[str, int] = {}
        for source in (go_data, kegg_data):
            for taxid, data in source.items():
                latin_name = data.get("latin_name", "").strip().casefold()
                if latin_name:
                    latin_name_to_taxid.setdefault(latin_name, taxid)

        # Read Reactome species-registry records.
        reactome_data: Dict[int, Dict[str, str]] = {}
        if reactome_registry and reactome_registry.exists():
            with open(reactome_registry, "r", encoding="utf-8") as fh:
                reader = csv.DictReader(fh, delimiter="\t")
                for row in reader:
                    taxid = int(row["taxid"])
                    reactome_data[taxid] = {
                        "reactome_code": row.get("reactome_code", ""),
                    }

        # Read Disease Ontology species coverage.
        do_taxids: set = set()
        if do_registry and do_registry.exists():
            with open(do_registry, "r", encoding="utf-8") as fh:
                reader = csv.DictReader(fh, delimiter="\t")
                for row in reader:
                    try:
                        do_taxids.add(int(row["taxid"]))
                    except (ValueError, KeyError):
                        continue

        # Read WikiPathways species-registry records.
        wikipathways_data: Dict[int, Dict[str, str]] = {}
        if wikipathways_registry and wikipathways_registry.exists():
            with open(wikipathways_registry, "r", encoding="utf-8") as fh:
                reader = csv.DictReader(fh, delimiter="\t")
                rows = list(reader)
            unresolved = {
                row.get("species_latin_name", "").strip().casefold()
                for row in rows
                if row.get("species_latin_name", "").strip().casefold() not in latin_name_to_taxid
            }
            taxonomy_matches: Dict[str, set] = {name: set() for name in unresolved if name}
            names_path = self.basic_dir / "taxonomy" / "names.dmp"
            if taxonomy_matches and names_path.exists():
                with open(names_path, "r", encoding="utf-8", errors="ignore") as fh:
                    for line in fh:
                        parts = [part.strip() for part in line.split("|")]
                        if len(parts) >= 2 and parts[1].casefold() in taxonomy_matches:
                            taxonomy_matches[parts[1].casefold()].add(int(parts[0]))
            for row in rows:
                species_name = row.get("species_latin_name", "").strip().casefold()
                taxid = latin_name_to_taxid.get(species_name)
                if not taxid and len(taxonomy_matches.get(species_name, ())) == 1:
                    taxid = next(iter(taxonomy_matches[species_name]))
                if not taxid:
                    logger.warning("The WikiPathways species cannot be uniquely mapped at TaxID: %s", row.get("species_latin_name", ""))
                    continue
                wikipathways_data[taxid] = {
                    "species_code": row.get("species_code", ""),
                    "gene_count": row.get("gene_count", ""),
                    "pathway_count": row.get("pathway_count", ""),
                }

        def read_simple_registry(path: Optional[Path]) -> Dict[int, str]:
            data: Dict[int, str] = {}
            if path and path.exists():
                with path.open("r", encoding="utf-8") as handle:
                    for row in csv.DictReader(handle, delimiter="\t"):
                        try:
                            taxid = int(row["taxid"])
                        except (KeyError, TypeError, ValueError):
                            continue
                        latin_name = row.get("latin_name", "").strip()
                        if latin_name:
                            data[taxid] = latin_name
            return data

        disgenet_data = read_simple_registry(disgenet_registry)
        trrust_data = read_simple_registry(trrust_registry)
        chea3_data = read_simple_registry(chea3_registry)
        animaltfdb_data = read_simple_registry(animaltfdb_registry)
        htftarget_data = read_simple_registry(htftarget_registry)

        # ----Merge----
        all_taxids = sorted(
            set(go_data)
            | set(kegg_data)
            | set(reactome_data)
            | do_taxids
            | set(wikipathways_data)
            | set(disgenet_data)
            | set(trrust_data)
            | set(chea3_data)
            | set(animaltfdb_data)
            | set(htftarget_data)
        )

        entries: List[SpeciesEntry] = []
        for taxid in all_taxids:
            # Determined latin_name: GO > Kegg > empty
            latin_name = ""
            if taxid in go_data and go_data[taxid]["latin_name"]:
                latin_name = go_data[taxid]["latin_name"]
            elif taxid in kegg_data and kegg_data[taxid]["latin_name"]:
                latin_name = kegg_data[taxid]["latin_name"]
            else:
                for source in (
                    disgenet_data,
                    trrust_data,
                    chea3_data,
                    animaltfdb_data,
                    htftarget_data,
                ):
                    if taxid in source:
                        latin_name = source[taxid]
                        break

            entry = SpeciesEntry(taxid=taxid, latin_name=latin_name)

            # GO
            if taxid in go_data:
                entry.has_go = True
                entry.go_source = go_data[taxid].get("source")
                gc = go_data[taxid].get("gene_count", "")
                tc = go_data[taxid].get("term_count", "")
                entry.go_gene_count = int(gc) if gc and gc.isdigit() else None
                entry.go_term_count = int(tc) if tc and tc.isdigit() else None

            # KEGG
            if taxid in kegg_data:
                entry.has_kegg = True
                entry.kegg_code = kegg_data[taxid].get("kegg_code", "")
                entry.kegg_code_source = kegg_data[taxid].get(
                    "kegg_code_source", "kegg"
                )
                gc = kegg_data[taxid].get("gene_count", "")
                entry.kegg_gene_count = int(gc) if gc and gc.isdigit() else None
            elif latin_name:
                # Auto Generate KEGG abbreviations
                entry.kegg_code = SpeciesRegistry.generate_kegg_abbreviation(
                    latin_name
                )
                entry.kegg_code_source = "auto"

            # Reactome
            if taxid in reactome_data:
                entry.has_reactome = True
                entry.reactome_code = reactome_data[taxid].get("reactome_code", "")

            # DO
            if taxid in do_taxids:
                entry.has_do = True

            if taxid in disgenet_data:
                entry.has_disgenet = True

            # WikiPathways
            if taxid in wikipathways_data:
                entry.has_wikipathways = True
                gc = wikipathways_data[taxid].get("gene_count", "")
                pc = wikipathways_data[taxid].get("pathway_count", "")
                entry.wikipathways_gene_count = int(gc) if gc and gc.isdigit() else None
                entry.wikipathways_pathway_count = int(pc) if pc and pc.isdigit() else None

            entry.has_trrust = taxid in trrust_data
            entry.has_chea3 = taxid in chea3_data
            entry.has_animaltfdb = taxid in animaltfdb_data and taxid != 9606
            entry.has_htftarget = taxid in htftarget_data

            entries.append(entry)

        # ----Writing supported_species.tsv----
        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)

        registry = SpeciesRegistry(registry_path=output_path)
        for entry in entries:
            registry.add_entry(entry)
        registry.save()

        logger.info(
            "Uniform species register has been generated: %s (%d species)", output_path, len(entries)
        )
        return output_path

    def _report_download_summary(self, registry_path: Path) -> None:
        """Print a concise summary of completed and failed downloads."""
        registry = SpeciesRegistry(registry_path=registry_path)
        registry.load()

        if not registry.entries:
            print("The species registry is empty; no coverage statistics are available.", file=sys.stderr)
            return

        summary = registry.get_summary()
        total = summary["total_species"]

        go_info = summary["go"]
        kegg_info = summary["kegg"]
        reactome_info = summary["reactome"]
        do_info = summary["do"]

        sep = "=" * 60
        print(f"\n{sep}", file=sys.stderr)
        print("Statistical summary of data downloads", file=sys.stderr)
        print(sep, file=sys.stderr)
        print(f"Total species: {total}", file=sys.stderr)
        print("-" * 60, file=sys.stderr)
        print(f"GO supported species: {go_info['count']}", file=sys.stderr)
        print(f"- Species with gene data: {go_info['with_gene_count']}", file=sys.stderr)
        print(f"- Species with term data: {go_info['with_term_count']}", file=sys.stderr)
        print("-" * 60, file=sys.stderr)
        print(f"KEGG supported species: {kegg_info['count']}", file=sys.stderr)
        print(f"- Species with gene data: {kegg_info['with_gene_count']}", file=sys.stderr)
        print(f"- Species with pathway data: {kegg_info['with_pathway_count']}", file=sys.stderr)
        print("-" * 60, file=sys.stderr)
        print(f"Reactome supported species: {reactome_info['count']}", file=sys.stderr)
        print(f"- Species with gene data: {reactome_info['with_gene_count']}", file=sys.stderr)
        print(f"- Species with pathway data: {reactome_info['with_pathway_count']}", file=sys.stderr)
        print("-" * 60, file=sys.stderr)
        print(f"DO supported species: {do_info['count']}", file=sys.stderr)
        print(f"- Species with gene data: {do_info['with_gene_count']}", file=sys.stderr)
        print(f"- Species with term data: {do_info['with_term_count']}", file=sys.stderr)
        print(sep, file=sys.stderr)

    # ============================
    # Internal support methods
    # ============================

    def _api_get_with_retry(
        self,
        url: str,
        headers: Optional[Dict[str, str]] = None,
        timeout: int = 30,
        max_retries: int = 3,
    ) -> str:
        """Perform an HTTP GET request with bounded retries."""
        if headers is None:
            headers = {
                "User-Agent": (
                    "AllEnricher/2.0 (https://github.com/allenricher; "
                    "data download pipeline)"
                ),
            }

        last_error: Optional[Exception] = None
        for attempt in range(1, max_retries + 1):
            try:
                resp = requests.get(url, headers=headers, timeout=timeout)
                resp.raise_for_status()
                return resp.text
            except requests.RequestException as exc:
                last_error = exc
                logger.warning(
                    "Request failed (No. 1)%d/%d(b) Seconds: %s - %s",
                    attempt, max_retries, url, exc,
                )
                if attempt < max_retries:
                    wait = 2 ** attempt
                    logger.info("Wait%dRetry in seconds...", wait)
                    time.sleep(wait)

        raise RuntimeError(
            f"API request failed (retryed){max_retries}(b) Seconds: {url}"
        ) from last_error
