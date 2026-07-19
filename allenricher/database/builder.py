"""Build species-specific enrichment databases from downloaded source data."""

import os
import gzip
import logging
import shutil
from pathlib import Path
from datetime import datetime
from typing import Dict, List, Optional

from .parsers.go import GOParser
from .parsers.kegg import KEGGParser
from .parsers.reactome import ReactomeParser
from .parsers.do import DOParser
from .parsers.disgenet import DisGeNETParser
from .parsers.wikipathways import WikiPathwaysParser
from .parsers.trrust import TRRUSTParser
from .parsers.chea3 import ChEA3Parser
from .downloader import DataDownloader
from .gmt_generator import GMTGenerator
from .species_registry import SpeciesRegistry, SpeciesEntry
from .goa_fetcher import GOAFetcher
from .wikipathways_fetcher import WikiPathwaysFetcher
from .trrust_fetcher import TRRUSTFetcher
from .manager import validate_tf_database_species

logger = logging.getLogger(__name__)


class DatabaseBuilder:
    """Build analysis-ready database artifacts for one or more species."""

    def __init__(self, root_dir: str = "./database"):
        """Initialize Builder

        Args:
root_dir: Database Root Directory, Must contain basic/ Subdirectorate
        """
        self.root_dir = Path(root_dir)
        self.basic_dir = self.root_dir / "basic"
        self.organism_dir = self.root_dir / "organism"

    # ============================
    # GO database construction
    # ============================
    def build_go(self, species: str, taxid: int,
                 go_version: Optional[str] = None) -> str:
        """Build the Gene Ontology database for one species."""
        # Resolve the most recent downloaded GO snapshot unless explicitly selected.
        if go_version is None:
            downloader = DataDownloader(root_dir=str(self.root_dir))
            go_version = downloader.get_latest_go_version()
            if go_version is None:
                raise FileNotFoundError(
                    "No Gene Ontology source data were found. Download them first:\n"
                    "  allenricher download go\n"
                    "or\n"
                    "  python -m allenricher.cli download -d go"
                )

        # Validate the selected source snapshot.
        go_basic = self.basic_dir / "go" / go_version
        if not go_basic.exists():
            raise FileNotFoundError(
                f"Gene Ontology source directory does not exist: {go_basic}\n"
                f"Available GO versions: {DataDownloader(root_dir=str(self.root_dir)).list_go_versions()}"
            )

        # Keep all databases built in the same run under one versioned directory.
        date_str = datetime.now().strftime("%Y%m%d")
        outdir = self.organism_dir / f"v{date_str}" / species
        outdir.mkdir(parents=True, exist_ok=True)

        gene2go_path = go_basic / "gene2go.gz"
        gene_info_path = go_basic / "gene_info.gz"
        obo_path = go_basic / "go-basic.obo"

        for fpath in [gene2go_path, gene_info_path, obo_path]:
            if not fpath.exists():
                raise FileNotFoundError(f"Missing file: {fpath}")

        print(f"\n{'='*60}")
        print(f"Build GO database: {species} (taxid={taxid})")
        print(f"Data source: {go_basic}")
        print(f"Output directory: {outdir}")
        print(f"{'='*60}")

        print("|--- Step 1/2: extract species annotations and build the GO membership matrix...")
        GOParser.parse_gene2go(
            gene2go_path=str(gene2go_path),
            gene_info_path=str(gene_info_path),
            taxid=taxid,
            species=species,
            outdir=str(outdir)
        )

        print("|--- Step 2/2: write GO term names from go-basic.obo...")
        GOParser.parse_obo(
            obo_path=str(obo_path),
            outdir=str(outdir)
        )

        # Report every required output explicitly.
        expected_files = [
            f"{species}.GO2gene.tab.gz",
            f"{species}.gene2go.txt",
            "GO2disc.gz"
        ]
        for fname in expected_files:
            fpath = outdir / fname
            if fpath.exists():
                print(f"    [OK] {fname}")
            else:
                print(f"    [MISSING] {fname}")

        print(f"\nGO Database build complete -> {outdir}")
        # Generate the canonical GMT representation used by all analysis methods.
        self.generate_gmt_files(species, str(outdir))
        return str(outdir)

    # ============================
    # GO Database Build (GoA Source)
    # ============================
    def _get_species_dir(self, taxid: int, latin_name: str) -> str:
        """Return the canonical TaxID-and-name directory for a species."""
        return f"{taxid}.{latin_name.replace(' ', '_')}"

    def _get_species_prefix(self, taxid: int, latin_name: str) -> str:
        """Return the filename-safe scientific-name prefix for a species."""
        return f"{taxid}.{latin_name.replace(' ', '_')}"

    def build_go_from_goa(self, taxid: int, latin_name: str,
                          goa_filename: str, go_version: str = None) -> str:
        """Build a Gene Ontology database from UniProt GOA annotations."""
        # Resolve a GO source version for provenance and ontology term names.
        if go_version is None:
            downloader = DataDownloader(root_dir=str(self.root_dir))
            go_version = downloader.get_latest_go_version()
            if go_version is None:
                date_str = datetime.now().strftime("%Y%m%d")
                go_version = f"GO{date_str}"
                logger.warning("No GO basic data version found, using default version: %s", go_version)

        # Use TaxID in paths so species identity is unambiguous.
        species_dir = self._get_species_dir(taxid, latin_name)
        prefix = self._get_species_prefix(taxid, latin_name)

        # Create the versioned species output directory.
        date_str = datetime.now().strftime("%Y%m%d")
        outdir = self.organism_dir / f"v{date_str}" / species_dir
        outdir.mkdir(parents=True, exist_ok=True)

        # GOA provides memberships; go-basic.obo provides stable term names.
        go_basic = self.basic_dir / "go" / go_version
        obo_path = go_basic / "go-basic.obo"

        print(f"\n{'='*60}")
        print(f"Build GO database (GOA source): {species_dir} (taxid={taxid})")
        print(f"GOA file: {goa_filename}")
        print(f"GO Version: {go_version}")
        print(f"Output directory: {outdir}")
        print(f"{'='*60}")

        # Retrieve the species-specific GOA source file.
        goa_date = go_version.replace("GO", "") if go_version.startswith("GO") else date_str
        goa_cache_dir = self.basic_dir / "goa" / f"GOA{goa_date}"
        fetcher = GOAFetcher(cache_dir=str(goa_cache_dir), overwrite=False)

        latin_name_underscore = latin_name.replace(" ", "_")
        goa_file = fetcher.fetch_species_data(
            taxid=taxid,
            latin_name=latin_name_underscore,
            goa_filename=goa_filename,
        )
        logger.info("GOA file ready: %s", goa_file)

        print("|--- Step 1/4: parse the GOA file...")
        gene_to_go, all_genes = fetcher.parse_goa_file(goa_file, taxid)

        if not all_genes:
            raise ValueError(
                f"The GOA file contains no valid annotations for TaxID {taxid}."
            )

        # Collect all terms represented by the species annotations.
        all_go_terms: set = set()
        for go_set in gene_to_go.values():
            all_go_terms.update(go_set)

        logger.info("Parsed %d genes and %d GO terms", len(all_genes), len(all_go_terms))

        # Resolve descriptive term names from the ontology snapshot.
        go_names: Dict[str, str] = {}
        if obo_path.exists():
            print("|--- Step 2/4: read GO term names from go-basic.obo...")
            go_names = self._extract_go_names_from_obo(str(obo_path))
            logger.info("Loaded %d GO term names from the ontology", len(go_names))
        else:
            logger.warning("go-basic.obo is missing at %s; term names will be unavailable", obo_path)

        # Write membership and description artifacts expected by DatabaseManager.
        print("|--- Step 3/4: Generate GO2gene.tab.gz...")
        go2gene_path = outdir / f"{prefix}.GO2gene.tab.gz"
        GOAFetcher.build_go2gene_matrix(
            gene_to_go=gene_to_go,
            all_genes=all_genes,
            all_go_terms=all_go_terms,
            output_path=go2gene_path,
        )

        print("|--- Step 3/4: write gene2go.txt...")
        gene2go_path = outdir / f"{prefix}.gene2go.txt"
        GOAFetcher.build_gene2go_list(
            gene_to_go=gene_to_go,
            go_names=go_names,
            output_path=gene2go_path,
        )

        print("|--- Step 4/4: write GO2disc.gz...")
        if obo_path.exists():
            GOParser.parse_obo(obo_path=str(obo_path), outdir=str(outdir))
        else:
            logger.warning("Skipped GO2disc.gz because go-basic.obo is unavailable")

        # Report every required output explicitly.
        expected_files = [
            f"{prefix}.GO2gene.tab.gz",
            f"{prefix}.gene2go.txt",
            "GO2disc.gz",
        ]
        for fname in expected_files:
            fpath = outdir / fname
            if fpath.exists():
                print(f"    [OK] {fname}")
            else:
                print(f"    [MISSING] {fname}")

        print(f"\nGO Database build complete (GOA Source) -> {outdir}")

        # Also provide GMT for GSEA, ssGSEA, GSVA, and custom workflows.
        self._generate_gmt_for_species_dir(species_dir, str(outdir), prefix)

        return str(outdir)

    def _extract_go_names_from_obo(self, obo_path: str) -> Dict[str, str]:
        """Read Gene Ontology term identifiers and names from an OBO file."""
        import re

        go_names: Dict[str, str] = {}
        go_id_pattern = re.compile(r'^id:\s(GO:\d+)')
        name_pattern = re.compile(r'^name:\s(.*)')

        with open(obo_path, 'r', encoding='utf-8') as f:
            current_id = None
            for line in f:
                line = line.strip()
                m = go_id_pattern.match(line)
                if m:
                    current_id = m.group(1)
                    continue
                m = name_pattern.match(line)
                if m and current_id:
                    go_names[current_id] = m.group(1)
                    current_id = None

        return go_names

    def _generate_gmt_for_species_dir(self, species_dir: str, output_dir: str,
                                       prefix: str) -> Dict[str, str]:
        """Generate GMT files for all built databases in one species directory."""
        generator = GMTGenerator(organism_dir=output_dir)
        results: Dict[str, str] = {}

        # Generate GO GMT only when both membership and description files exist.
        tab_path = Path(output_dir) / f"{prefix}.GO2gene.tab.gz"
        disc_path = Path(output_dir) / "GO2disc.gz"
        if tab_path.exists() and disc_path.exists():
            try:
                terms, term_to_genes = generator._read_tab_matrix(str(tab_path))
                descriptions = generator._read_description(str(disc_path))
                gmt_path = str(Path(output_dir) / f"{prefix}.GO.gmt.gz")
                generator._write_gmt(term_to_genes, descriptions, gmt_path)
                results["GO"] = gmt_path
            except Exception as e:
                logger.warning("Failed to generate the GO GMT file: %s", e)

        return results

    # ============================
    # GO database build with an explicit NCBI-to-UniProt fallback.
    # ============================
    def build_go_with_fallback(self, taxid: int, latin_name: str,
                               go_version: str = None) -> str:
        """Build GO from NCBI gene2go, falling back to UniProt GOA when necessary."""
        logger.info("build_go_with_fallback: taxid=%d, latin_name=%s", taxid, latin_name)

        # Prefer NCBI gene2go when it contains the requested TaxID.
        use_gene2go = False
        species_abbr = None

        # The registry records the source selected during download preparation.
        registry_path = self.root_dir / "supported_species.tsv"
        if registry_path.exists():
            registry = SpeciesRegistry(registry_path=registry_path)
            registry.load()
            entry = registry.query_by_taxid(taxid)
            if entry and entry.has_go and entry.go_source:
                if entry.go_source.lower() == "gene2go":
                    use_gene2go = True
                    logger.info("TaxID %d is registered with NCBI gene2go coverage", taxid)

        # Legacy registries may lack source metadata, so inspect gene2go directly.
        if not use_gene2go:
            use_gene2go = self._check_gene2go_has_taxid(taxid, go_version)

        # Build from NCBI when a project species code can be resolved.
        if use_gene2go:
            # Try to get abbreviations
            species_abbr = self._get_species_abbr(taxid, latin_name)
            if species_abbr:
                logger.info("Building GO from NCBI gene2go for species=%s", species_abbr)
                return self.build_go(species=species_abbr, taxid=taxid,
                                     go_version=go_version)
            else:
                logger.warning("Could not resolve a project species code for TaxID %d", taxid)

        # Fall back to the UniProt GOA file recorded for the species.
        goa_filename = self._find_goa_filename(taxid, latin_name)
        if goa_filename is None:
            raise ValueError(
                f"No GO annotation source is available for TaxID {taxid} ({latin_name}).\n"
                "The species is absent from NCBI gene2go and has no UniProt GOA registry entry."
            )

        logger.info("Building GO from UniProt GOA file %s", goa_filename)
        return self.build_go_from_goa(taxid, latin_name, goa_filename, go_version)

    def _check_gene2go_has_taxid(self, taxid: int,
                                  go_version: str = None) -> bool:
        """Return whether gene2go contains annotations for the requested TaxID."""
        if go_version is None:
            downloader = DataDownloader(root_dir=str(self.root_dir))
            go_version = downloader.get_latest_go_version()
            if go_version is None:
                return False

        gene2go_path = self.basic_dir / "go" / go_version / "gene2go.gz"
        if not gene2go_path.exists():
            return False

        taxid_str = str(taxid)

        try:
            with gzip.open(gene2go_path, 'rt', encoding='utf-8') as f:
                for i, line in enumerate(f, 1):
                    line = line.strip()
                    if not line or line.startswith('#'):
                        continue
                    parts = line.split('\t')
                    if parts[0] == taxid_str:
                        logger.info("Found TaxID %s in gene2go.gz at data line %d", taxid_str, i)
                        return True
        except Exception as e:
            logger.warning("Checking gene2go.gz failed: %s", e)
            return False

        logger.info("TaxID %s is absent from gene2go.gz", taxid_str)
        return False

    def _get_species_abbr(self, taxid: int, latin_name: str) -> Optional[str]:
        """Resolve the project species code for a scientific name or TaxID."""
        # Prefer the downloaded TaxID-keyed registry.
        registry_path = self.root_dir / "supported_species.tsv"
        if registry_path.exists():
            registry = SpeciesRegistry(registry_path=registry_path)
            registry.load()
            entry = registry.query_by_taxid(taxid)
            if entry and entry.kegg_code:
                return entry.kegg_code

        # Derive a candidate code only when the registry lacks one.
        return SpeciesRegistry.generate_kegg_abbreviation(latin_name)

    def _find_goa_filename(self, taxid: int, latin_name: str) -> Optional[str]:
        """Find the UniProt GOA filename associated with a species."""
        # Prefer the downloaded TaxID-keyed registry.
        registry_path = self.root_dir / "supported_species.tsv"
        if registry_path.exists():
            registry = SpeciesRegistry(registry_path=registry_path)
            registry.load()
            entry = registry.query_by_taxid(taxid)
            if entry and entry.has_go and entry.go_filename:
                return entry.go_filename

        # UniProt GOA proteome files conventionally use the TaxID as the stem.
        return f"{taxid}.goa"

    # ============================
    # Reactome Database Construction
    # ============================
    def build_reactome(self, species: str, taxid: int,
                       reactome_version: Optional[str] = None) -> str:
        """Build the Reactome pathway database for one species."""
        if reactome_version is None:
            downloader = DataDownloader(root_dir=str(self.root_dir))
            reactome_version = downloader.get_latest_reactome_version()
            if reactome_version is None:
                raise FileNotFoundError(
                    "No basic data for Reactome found. Run download first: \n"
                    "  allenricher download reactome\n"
                    "or\n"
                    "  python -m allenricher.cli download -d reactome"
                )

        re_basic = self.basic_dir / "reactome" / reactome_version
        if not re_basic.exists():
            raise FileNotFoundError(
                f"Reactome Basic Data Directory does not exist: {re_basic}"
            )

        date_str = datetime.now().strftime("%Y%m%d")
        outdir = self.organism_dir / f"v{date_str}" / species
        outdir.mkdir(parents=True, exist_ok=True)

        ncbi2reactome_path = re_basic / "NCBI2Reactome_All_Levels.txt.gz"
        gene_info_path = re_basic / "gene_info.gz"
        pathways_path = re_basic / "ReactomePathways.txt"
        relations_path = re_basic / "ReactomePathwaysRelation.txt"

        for fpath in [ncbi2reactome_path, gene_info_path, pathways_path, relations_path]:
            if not fpath.exists():
                raise FileNotFoundError(f"Missing file: {fpath}")

        print(f"\n{'='*60}")
        print(f"Construct Reactome database: {species} (taxid={taxid})")
        print(f"Data sources: {re_basic}")
        print(f"Output directory: {outdir}")
        print(f"{'='*60}")

        print("|--- Extract species pathways from the NCBI-to-Reactome mapping...")
        ReactomeParser.parse_ncbi2reactome(
            ncbi2reactome_path=str(ncbi2reactome_path),
            gene_info_path=str(gene_info_path),
            taxid=taxid,
            species=species,
            outdir=str(outdir),
            pathways_path=str(pathways_path),
            relations_path=str(relations_path),
        )

        # Report every required output explicitly.
        expected_files = [
            f"{species}.Reactome2gene.tab.gz",
            f"{species}.Reactome2disc.gz",
            f"{species}.gene2pathway.txt"
        ]
        for fname in expected_files:
            fpath = outdir / fname
            if fpath.exists():
                print(f"    [OK] {fname}")
            else:
                print(f"    [MISSING] {fname}")

        print(f"\nReactome Database build complete -> {outdir}")
        # Generate the canonical GMT representation.
        self.generate_gmt_files(species, str(outdir))
        return str(outdir)

    # ============================
    # KEGG database construction.
    # ============================
    def build_kegg(self, species: str, taxid: int,
                   go_version: Optional[str] = None,
                   gene2pathway_path: Optional[str] = None,
                   pathway_summary_path: Optional[str] = None) -> str:
        """Build the KEGG pathway database for one species."""
        date_str = datetime.now().strftime("%Y%m%d")
        outdir = self.organism_dir / f"v{date_str}" / species
        outdir.mkdir(parents=True, exist_ok=True)

        # Get gene_info from GO basic (if not specifically provided)
        if go_version is None:
            downloader = DataDownloader(root_dir=str(self.root_dir))
            go_version = downloader.get_latest_go_version()
            if go_version is None:
                go_version = f"GO{date_str}"

        gene_info_path = self.basic_dir / "go" / go_version / "gene_info.gz"
        if not gene_info_path.exists():
            # Try to get from real base
            re_ver = DataDownloader(root_dir=str(self.root_dir)).get_latest_reactome_version()
            if re_ver:
                alt_path = self.basic_dir / "reactome" / re_ver / "gene_info.gz"
                if alt_path.exists():
                    gene_info_path = alt_path

        print(f"\n{'='*60}")
        print(f"Build KEGG database: {species} (taxid={taxid})")
        print(f"Output directory: {outdir}")
        print(f"{'='*60}")

        if gene2pathway_path is None:
            # Retrieve missing species data through the KEGG REST API.
            from .kegg_fetcher import KEGGFetcher
            fetcher = KEGGFetcher(
                cache_dir=str(self.basic_dir / "kegg"),
                overwrite=False,
            )
            try:
                gene2pathway_path, pathway_summary_path = fetcher.fetch_species_data(
                    species=species,
                    gene_info_path=str(gene_info_path),
                    taxid=taxid,
                )
            except Exception as e:
                print(f"|--- [ERROR] KEGG REST retrieval failed: {e}")
                print("Check the network connection and try again.")
                raise RuntimeError(f"KEGG REST retrieval failed: {e}") from e

        KEGGParser.build_database(
            species=species,
            gene_info_path=str(gene_info_path),
            gene2pathway_path=gene2pathway_path,
            outdir=str(outdir),
            pathway_summary_path=pathway_summary_path,
            taxid=taxid,
        )

        print(f"\nKEGG Database build complete -> {outdir}")
        # Generate the canonical GMT representation.
        self.generate_gmt_files(species, str(outdir))
        return str(outdir)

    # ============================
    # DO / DisGeNET (humans only)
    #
    # v1 make_speciesDB L127-135:
    #   if [ $organism == "hsa" ]; then
    #     sh $bin/src/makeDB.do.v1.0.sh ...
    #     sh $bin/src/makeDB.DisGeNET.v1.0.sh ...
    #   fi
    #
    # Both sources are human-specific; build_species_db enforces species=hsa.
    # ============================
    def build_do(self, taxid: int,
                 go_version: Optional[str] = None) -> str:
        """Build the human Disease Ontology database."""
        species = "hsa"
        date_str = datetime.now().strftime("%Y%m%d")
        outdir = self.organism_dir / f"v{date_str}" / species
        outdir.mkdir(parents=True, exist_ok=True)

        # Accept compressed or uncompressed Jensen Lab disease association files.
        do_dir = self.basic_dir / "do"
        disease_files = sorted(do_dir.glob("human_disease_*_filtered.tsv.gz"))
        if not disease_files:
            disease_files = sorted(do_dir.glob("human_disease_*_filtered.tsv"))
        if not disease_files:
            disease_files = sorted(do_dir.glob("human_disease_*.tsv.gz"))
        if not disease_files:
            disease_files = sorted(do_dir.glob("human_disease_*.tsv"))
        if not disease_files:
            raise FileNotFoundError(
                f"No Disease Ontology association files were found in {do_dir}. Download them first:\n"
                f"  allenricher download do"
            )

        # NCBI gene_info constrains associations to valid human gene symbols.
        if go_version is None:
            downloader = DataDownloader(root_dir=str(self.root_dir))
            go_version = downloader.get_latest_go_version()
        gene_info_path = self.basic_dir / "go" / go_version / "gene_info.gz"
        ontology_path = do_dir / "doid.obo"
        if not ontology_path.is_file():
            raise FileNotFoundError(
                f"Disease Ontology file is missing: {ontology_path}. Rerun `allenricher download do`."
            )

        print(f"\n{'='*60}")
        print(f"Build DO database (taxid={taxid})")
        print(f"Data source: {do_dir}")
        print(f"Output directory: {outdir}")
        print(f"{'='*60}")

        DOParser.parse_disease_files(
            disease_files=[str(f) for f in disease_files],
            gene_info_path=str(gene_info_path),
            taxid=taxid,
            outdir=str(outdir),
            ontology_path=str(ontology_path),
        )

        print(f"\nDO Database build complete -> {outdir}")
        # Generate the canonical GMT representation.
        self.generate_gmt_files("hsa", str(outdir))
        return str(outdir)

    def build_disgenet(self, taxid: int,
                       go_version: Optional[str] = None) -> str:
        """Build the human DisGeNET database from the licensed snapshot."""
        species = "hsa"
        date_str = datetime.now().strftime("%Y%m%d")
        outdir = self.organism_dir / f"v{date_str}" / species
        outdir.mkdir(parents=True, exist_ok=True)

        legacy_dirs = sorted(self.organism_dir.glob("v*/hsa"))
        for legacy_dir in reversed(legacy_dirs):
            gene_file = legacy_dir / "hsa.CUI2gene.tab.gz"
            disc_file = legacy_dir / "hsa.CUI2disc.gz"
            if legacy_dir != outdir and gene_file.is_file() and disc_file.is_file():
                print("NOTE: DisGeNET is still maintained, but current releases are no longer freely downloadable.")
                print(f"|--- Reusing the retained AllEnricher v1 snapshot from {legacy_dir}")
                for source in (gene_file, disc_file, legacy_dir / "hsa.DisGeNET.gmt.gz"):
                    if source.is_file():
                        shutil.copy2(source, outdir / source.name)
                return str(outdir)

        assoc_path = self.basic_dir / "disgenet" / "all_gene_disease_associations.tsv.gz"
        if not assoc_path.exists():
            raise FileNotFoundError(
                f"No licensed DisGeNET association file was found at {assoc_path}. "
                "Current DisGeNET downloads require authorization. AllEnricher can instead reuse "
                "the retained v1 hsa.CUI2gene.tab.gz and hsa.CUI2disc.gz snapshot when those files "
                "are present in a legacy human database directory."
            )

        if go_version is None:
            downloader = DataDownloader(root_dir=str(self.root_dir))
            go_version = downloader.get_latest_go_version()
        gene_info_path = self.basic_dir / "go" / go_version / "gene_info.gz"

        print(f"\n{'='*60}")
        print(f"Build DisGeNET database (taxid={taxid})")
        print(f"Data source: {assoc_path}")
        print(f"Output directory: {outdir}")
        print(f"{'='*60}")

        DisGeNETParser.parse_associations(
            assoc_path=str(assoc_path),
            gene_info_path=str(gene_info_path),
            taxid=taxid,
            outdir=str(outdir)
        )

        print(f"\nDisGeNET Database build complete -> {outdir}")
        # Generate the canonical GMT representation.
        self.generate_gmt_files("hsa", str(outdir))
        return str(outdir)

    # ============================
    # Shared GO-version and NCBI gene-info helpers.
    # ============================
    def _get_go_version(self) -> Optional[str]:
        """Return the version recorded for the local Gene Ontology snapshot."""
        go_root = self.basic_dir / "go"
        if go_root.is_dir():
            versions = sorted(
                path.name for path in go_root.iterdir()
                if path.is_dir() and path.name.startswith("GO")
            )
            if versions:
                return versions[-1]
        downloader = DataDownloader(root_dir=str(self.root_dir))
        return downloader.get_latest_go_version()

    def _get_gene_info_path(self) -> Optional[Path]:
        """Locate the NCBI gene_info file used for gene identifier validation."""
        go_version = self._get_go_version()
        if go_version is None:
            return None

        gene_info_path = self.basic_dir / "go" / go_version / "gene_info.gz"
        if gene_info_path.exists():
            return gene_info_path

        return None

    def _load_valid_genes(self, gene_info_path: str, taxid: int) -> Optional[set]:
        """Load valid gene symbols for one NCBI TaxID."""
        taxid_str = str(taxid)
        valid_genes = set()
        try:
            with gzip.open(gene_info_path, 'rt', encoding='utf-8') as f:
                for line in f:
                    line = line.strip()
                    if not line or line.startswith('#'):
                        continue
                    parts = line.split('\t')
                    if len(parts) < 3:
                        continue
                    if parts[0] == taxid_str:
                        gene_symbol = parts[2]
                        if gene_symbol and gene_symbol != '-':
                            valid_genes.add(gene_symbol)
            logger.info("Loaded %d valid gene symbols from gene_info.gz for TaxID %s",
                        len(valid_genes), taxid_str)
            return valid_genes
        except Exception as e:
            logger.warning("Failed to load gene_info.gz: %s", e)
            return None

    # ============================
    # WikiPathways database construction
    # ============================
    def build_wikipathways(self, species: str, taxid: int,
                           wikipathways_version: Optional[str] = None) -> str:
        """Build the WikiPathways database for one species."""
        # WikiPathways filenames are keyed by scientific name.
        latin_name = WikiPathwaysFetcher.get_latin_name(species)
        if latin_name is None:
            raise ValueError(
                f"No WikiPathways scientific name is registered for species '{species}'."
            )

        # Use the newest cached release unless a version was requested.
        if wikipathways_version is None:
            fetcher = WikiPathwaysFetcher(basic_dir=str(self.basic_dir))
            versions = fetcher.list_cached_versions()
            if not versions:
                raise FileNotFoundError(
                    "No WikiPathways source data were found. Download them first:\n"
                    "  allenricher download wikipathways\n"
                    "or\n"
                    "  python -m allenricher.cli download -d wikipathways"
                )
            wikipathways_version = versions[-1]  # Use the latest version

        # Validate the selected WikiPathways release directory.
        wp_basic = self.basic_dir / "wikipathways" / f"WP{wikipathways_version}"
        if not wp_basic.exists():
            raise FileNotFoundError(
                f"The WikiPathways base data directory does not exist: {wp_basic}\n"
                f"Available versions: "
                f"{WikiPathwaysFetcher(basic_dir=str(self.basic_dir)).list_cached_versions()}"
            )

        # Resolve the official species-specific GMT filename.
        species_filename = latin_name.replace(" ", "_")
        gmt_filename = f"wikipathways-{wikipathways_version}-gmt-{species_filename}.gmt"
        gmt_path = wp_basic / gmt_filename

        if not gmt_path.exists():
            raise FileNotFoundError(
                f"The WikiPathways GMT file does not exist: {gmt_path}\n"
                f"Download WikiPathways data for {species} ({latin_name}) before building the database."
            )

        # Output Directory
        date_str = datetime.now().strftime("%Y%m%d")
        outdir = self.organism_dir / f"v{date_str}" / species
        outdir.mkdir(parents=True, exist_ok=True)

        # Use NCBI gene_info for identifier conversion when available.
        gene_info_path = self._get_gene_info_path()

        print(f"\n{'='*60}")
        print(f"Build WikiPathways database: {species} (taxid={taxid})")
        print(f"Latin name: {latin_name}")
        print(f"Data source: {gmt_path}")
        print(f"Output directory: {outdir}")
        print(f"{'='*60}")

        print("|--- Parse WikiPathways GMT gene sets and write database artifacts...")
        WikiPathwaysParser.build_database(
            gmt_path=str(gmt_path),
            output_dir=str(outdir),
            species=species,
            taxid=taxid,
            gene_info_path=str(gene_info_path) if gene_info_path else None
        )

        # Report every required output explicitly.
        expected_files = [
            f"{species}.WikiPathways2gene.tab.gz",
            f"{species}.WikiPathways2disc.gz",
        ]
        for fname in expected_files:
            fpath = outdir / fname
            if fpath.exists():
                print(f"    [OK] {fname}")
            else:
                print(f"    [MISSING] {fname}")

        print(f"\nWikiPathways Database build complete -> {outdir}")
        # Generate the canonical GMT representation.
        self.generate_gmt_files(species, str(outdir))
        return str(outdir)

    # ============================
    # TRRUST database build
    # ============================
    def build_trrust(self, species: str, taxid: int) -> str:
        """Build the TRRUST regulatory database for human or mouse."""
        # TRRUST publishes human and mouse datasets only.
        latin_name = TRRUSTFetcher.get_latin_name(species)
        if latin_name is None:
            raise ValueError(
                f"TRRUST does not support species '{species}'. "
                f"Supported species: {TRRUSTFetcher.get_supported_species()}"
            )

        # Locate the downloaded TRRUST v2 source table.
        trrust_raw = self.basic_dir / "trrust" / "TRRUSTv2" / f"trrust_rawdata.{species}.tsv"
        if not trrust_raw.exists():
            raise FileNotFoundError(
                f"The TRRUST source file does not exist: {trrust_raw}\n"
                "Download TRRUST data first:\n"
                f"  allenricher download trrust"
            )

        # Output Directory
        date_str = datetime.now().strftime("%Y%m%d")
        outdir = self.organism_dir / f"v{date_str}" / species
        outdir.mkdir(parents=True, exist_ok=True)

        # Filter targets to valid species gene symbols when possible.
        gene_info_path = self._get_gene_info_path()
        valid_genes = None
        if gene_info_path:
            valid_genes = self._load_valid_genes(str(gene_info_path), taxid)

        print(f"\n{'='*60}")
        print(f"Build TRRUST database: {species} (taxid={taxid})")
        print(f"Latin name: {latin_name}")
        print(f"Data source: {trrust_raw}")
        print(f"Output directory: {outdir}")
        print(f"{'='*60}")

        print("|--- Parse TRRUST TF-target interactions...")
        TRRUSTParser.build_database(
            tsv_path=str(trrust_raw),
            output_dir=str(outdir),
            species=species,
            valid_genes=valid_genes
        )

        # Report every required output explicitly.
        expected_files = [
            f"{species}.TF2target.tab.gz",
            f"{species}.gene2TF.tab.gz",
            f"{species}.TF2disc.gz",
            f"{species}.TRRUST_edges.tsv.gz",
        ]
        for fname in expected_files:
            fpath = outdir / fname
            if fpath.exists():
                print(f"    [OK] {fname}")
            else:
                print(f"    [MISSING] {fname}")

        print(f"\nTRRUST Database build complete -> {outdir}")
        # Generate the canonical GMT representation.
        self.generate_gmt_files(species, str(outdir))
        return str(outdir)

    # ============================
    # ChEA3 database construction
    # ============================
    def build_chea3(self, species: str, taxid: int,
                    merge_method: str = "separate") -> str:
        """Build the human ChEA3 regulatory database."""
        validate_tf_database_species("ChEA3", species)

        # ChEA3 is a human-specific collection of curated TF libraries.
        chea3_dir = self.basic_dir / "chea3" / "ChEA3v2024"
        if not chea3_dir.exists():
            raise FileNotFoundError(
                f"The ChEA3 base data directory does not exist: {chea3_dir}\n"
                "Download ChEA3 data first:\n"
                f"  allenricher download chea3"
            )

        gmt_files = sorted(chea3_dir.glob("*_tf.gmt"))
        if not gmt_files:
            raise FileNotFoundError(
                f"No ChEA3 *_tf.gmt libraries were found in {chea3_dir}.\n"
                "Download ChEA3 data first:\n"
                f"  allenricher download chea3"
            )

        # Output Directory
        date_str = datetime.now().strftime("%Y%m%d")
        outdir = self.organism_dir / f"v{date_str}" / species
        outdir.mkdir(parents=True, exist_ok=True)

        # Restrict target genes to valid human symbols when gene_info is available.
        gene_info_path = self._get_gene_info_path()
        valid_genes = None
        if gene_info_path:
            valid_genes = self._load_valid_genes(str(gene_info_path), taxid)
            if valid_genes:
                print(f"|--- Valid gene-symbol filter: {len(valid_genes)} genes")

        print(f"\n{'='*60}")
        print(f"Build ChEA3 database: {species} (taxid={taxid})")
        print(f"Data source: {chea3_dir} ({len(gmt_files)} GMT libraries)")
        print(f"Library handling: {merge_method}")
        print(f"Output directory: {outdir}")
        print(f"{'='*60}")

        # Build a database using ChEA3Parser
        print("|---Extract TF-target relationships from ChEA3 GMT files...")
        ChEA3Parser.build_database(
            gmt_paths=[str(f) for f in gmt_files],
            output_dir=str(outdir),
            species=species,
            merge_method=merge_method,
            valid_genes=valid_genes,
        )

        # Report every required output explicitly.
        expected_files = [
            f"{species}.ChEA3_2gene.tab.gz",
            f"{species}.ChEA3_2disc.gz",
        ]
        for fname in expected_files:
            fpath = outdir / fname
            if fpath.exists():
                print(f"    [OK] {fname}")
            else:
                print(f"    [MISSING] {fname}")

        print(f"\nChEA3 Database build complete -> {outdir}")
        # Generate the canonical GMT representation.
        self.generate_gmt_files(species, str(outdir))
        return str(outdir)

    # ============================
    # GMT Gene Set File Generation
    # ============================
    def generate_gmt_files(self, species: str,
                           output_dir: str = None) -> Dict[str, str]:
        """Generate GMT files from all available database artifacts."""
        if output_dir is None:
            # Automatically find the most recent species database directory
            if not self.organism_dir.exists():
                print("|---[Warning] The species database directory does not exist, skips GMT generation")
                return {}
            species_dirs = sorted(
                self.organism_dir.glob(f"*/{species}"),
                key=lambda p: p.parent.name,
                reverse=True
            )
            if not species_dirs:
                print(f"|---[Warning] No species found{species}Database directory, skip GMT generation")
                return {}
            output_dir = str(species_dirs[0])

        generator = GMTGenerator(organism_dir=output_dir)
        return generator.generate_all_gmt(species)

    # ============================
    # hTFtarget database construction (human exclusive)
    # ============================
    def build_htftarget(self, species: str, taxid: int) -> str:
        """Build the human hTFtarget regulatory database."""
        from .parsers.htftarget import HTFtargetParser

        validate_tf_database_species("hTFtarget", species)

        htftarget_file = self.basic_dir / "htftarget" / "tf-target-information.txt"

        if not htftarget_file.exists():
            print(f"|---[ Skipped] No HTFtarget file found: {htftarget_file}")
            print("|---Run first: --animaltfdb")
            return ""

        date_str = datetime.now().strftime("%Y%m%d")
        outdir = self.organism_dir / f"v{date_str}" / species
        outdir.mkdir(parents=True, exist_ok=True)

        print(f"\n{'='*60}")
        print(f"Build hTFtarget database: {species} (taxid={taxid})")
        print(f"Data source: {htftarget_file}")
        print(f"Output directory: {outdir}")
        print(f"{'='*60}")

        gene_info_path = self._get_gene_info_path()
        valid_genes = None
        if gene_info_path:
            valid_genes = self._load_valid_genes(str(gene_info_path), taxid)

        HTFtargetParser.build_database(
            tsv_path=str(htftarget_file),
            output_dir=str(outdir),
            species=species,
            valid_genes=valid_genes,
        )

        expected_files = [
            f"{species}.hTF_2gene.tab.gz",
            f"{species}.hTF_2disc.gz",
        ]
        for fname in expected_files:
            fpath = outdir / fname
            if fpath.exists():
                print(f"    [OK] {fname}")
            else:
                print(f"    [MISSING] {fname}")

        # Record human hTFtarget coverage in the unified species registry.
        if species == 'hsa':
            try:
                # Fetch TF statistics from hTFtarget files
                human_tf_to_targets, _, _ = HTFtargetParser.parse_tsv(str(htftarget_file))
                tf_count = len(human_tf_to_targets) if human_tf_to_targets else 0
                
                registry = SpeciesRegistry.load_default(str(self.root_dir))
                registry.update_animaltfdb_stats(
                    species_code='hsa',
                    tf_count=tf_count,
                    mapped_target_count=0,  # hTFtarget direct data, non-map
                    has_data=True
                )
                registry.save()
                print(f"[OK] Species Register updated: hsa -> hTFtarget")
            except Exception as e:
                print(f"[WARNING] Failed to update species register: {e}")

        print(f"\nhTFtarget Database build complete -> {outdir}")
        return str(outdir)

    # ============================
    # AnimalTFDB database build through target-species orthology.
    # ============================
    def build_animaltfdb(self, species: str, taxid: int,
                         species_latin: str = "",
                         gene_info_path: Optional[str] = None) -> str:
        """Build a species-specific AnimalTFDB database through orthology mapping."""
        from .parsers.animaltfdb import AnimalTFDBParser
        from .parsers.htftarget import HTFtargetParser
        from .ortholog_mapper import OrthologMapper

        validate_tf_database_species("AnimalTFDB", species)

        if not species_latin:
            registry = SpeciesRegistry.load_default(str(self.root_dir))
            entry = registry.query_by_kegg_code(species)
            if entry and entry.latin_name:
                species_latin = entry.latin_name.replace(' ', '_')
            else:
                print(f"|--- [ERROR] No scientific name is registered for species '{species}'. Use --latin-name.")
                return ""

        cache_dir = self.basic_dir / "animaltfdb" / "AnimalTFDBv4.0"
        tf_list_file = cache_dir / f"{species_latin}_TF"
        ortholog_file = cache_dir / f"{species_latin}_ortholog_to_human"
        htftarget_file = self.basic_dir / "htftarget" / "tf-target-information.txt"

        missing = []
        if not tf_list_file.exists():
            missing.append(f"TF List: {tf_list_file}")
        if not ortholog_file.exists():
            missing.append(f"Orthology mapping: {ortholog_file}")
        if not htftarget_file.exists():
            missing.append(f"hTFtarget: {htftarget_file}")

        if missing:
            print("|--- [SKIPPED] Required AnimalTFDB inputs are missing:")
            for m in missing:
                print(f"|---   - {m}")
            print(f"|--- Download AnimalTFDB data first for species {species_latin}")
            return ""

        date_str = datetime.now().strftime("%Y%m%d")
        outdir = self.organism_dir / f"v{date_str}" / species
        outdir.mkdir(parents=True, exist_ok=True)

        print(f"\n{'='*60}")
        print(f"Build AnimalTFDB database: {species} ({species_latin}, taxid={taxid})")
        print(f"Mapping strategy: human hTFtarget targets -> {species_latin} orthologs")
        print(f"Output directory: {outdir}")
        print(f"{'='*60}")

        gene_info_path = Path(gene_info_path) if gene_info_path else self._get_gene_info_path()
        valid_genes = None
        if gene_info_path:
            valid_genes = self._load_valid_genes(str(gene_info_path), taxid)

        print("\n[1/3] Parse AnimalTFDB TF and orthology data...")
        tf_df, ortholog_map = AnimalTFDBParser.build_database(
            tf_list_path=str(tf_list_file),
            ortholog_path=str(ortholog_file),
            output_dir=str(outdir),
            species=species,
            valid_genes=valid_genes,
            gene_info_path=str(gene_info_path) if gene_info_path else None,
            species_taxid=taxid,
        )

        print("[2/3] Parse human hTFtarget regulatory contexts...")
        human_tf_to_targets, _, human_term_metadata = HTFtargetParser.parse_context_terms(
            str(htftarget_file)
        )

        print("[3/3] Map TF-target gene sets through orthology...")
        species_tf_set = set(tf_df['Symbol'].values) if 'Symbol' in tf_df.columns else None

        mapper = OrthologMapper(
            human_tf_to_targets=human_tf_to_targets,
            species_to_human=ortholog_map,
            species_tf_set=species_tf_set,
            human_term_metadata=human_term_metadata,
        )

        tf_to_targets, gene_to_tfs = mapper.map_tf_targets()

        # Summarize mapping coverage before writing database files.
        dedup_stats = mapper.get_duplicate_stats()
        mapped_tf_names = {item['TF'] for item in mapper.mapped_term_metadata.values()}
        mapped_targets_by_tf = {tf: set() for tf in mapped_tf_names}
        for term_id, targets in tf_to_targets.items():
            tf = mapper.mapped_term_metadata[term_id]['TF']
            mapped_targets_by_tf[tf].update(targets)
        mapping_stats = {
            'species': species,
            'species_latin': species_latin,
            'total_species_genes': len(ortholog_map),
            'total_human_tfs': len({item['TF'] for item in human_term_metadata.values()}),
            'species_tfs_found': len(species_tf_set) if species_tf_set else 0,
            'mapped_tfs': len(mapped_tf_names),
            'mapped_terms': len(tf_to_targets),
            'mapped_targets': len(gene_to_tfs),
            'avg_targets_per_tf': (
                sum(len(targets) for targets in mapped_targets_by_tf.values()) / len(mapped_targets_by_tf)
                if mapped_targets_by_tf else 0
            ),
            'multi_mapping_human_genes': dedup_stats.get('multi_mapping_count', 0),
            'coverage_ratio': len(mapped_tf_names) / (len(species_tf_set) if species_tf_set else 1) * 100,
        }

        # Print a compact quality summary for manual review.
        print("\n" + "="*60)
        print("Orthology Mapping Summary")
        print("="*60)
        print(f"Species: {species_latin} ({species})")
        print(f"Species genes with human orthologs: {mapping_stats['total_species_genes']}")
        print(f"Human TFs in the source library: {mapping_stats['total_human_tfs']}")
        print(f"TFs annotated for the target species: {mapping_stats['species_tfs_found']}")
        print(f"Mapped TFs: {mapping_stats['mapped_tfs']} ({mapping_stats['coverage_ratio']:.1f}% coverage)")
        print(f"Mapped target genes: {mapping_stats['mapped_targets']}")
        print(f"Mean mapped targets per TF: {mapping_stats['avg_targets_per_tf']:.1f}")
        print(f"Human genes with ambiguous mappings: {mapping_stats['multi_mapping_human_genes']}")
        print("="*60)

        # Low coverage makes the inferred network unsuitable as a complete TF catalogue.
        if mapping_stats['coverage_ratio'] < 50:
            print(f"\nWARNING: TF mapping coverage is only {mapping_stats['coverage_ratio']:.1f}%.")
            print("The available direct human orthologs do not cover most annotated TFs in this species.")
            print("Interpret the inferred TF-target database as incomplete.")

        # Persist the same mapping summary for reproducible review.
        import json
        stats_file = outdir / f"{species}.AnimalTFDB_mapping_stats.json"
        with open(stats_file, 'w', encoding='utf-8') as f:
            json.dump(mapping_stats, f, indent=2, ensure_ascii=False)
        print(f"Mapping summary saved: {stats_file}")

        if not tf_to_targets:
            print("|--- [WARNING] Orthology mapping produced no TF-target gene sets")
            return str(outdir)

        OrthologMapper.build_mapped_database(
            tf_to_targets=tf_to_targets,
            gene_to_tfs=gene_to_tfs,
            species_tf_df=tf_df,
            output_dir=str(outdir),
            species=species,
            term_metadata=mapper.mapped_term_metadata,
        )

        expected_files = [
            f"{species}.AnimalTFDB_2tf.tab.gz",
            f"{species}.AnimalTFDB_2disc.gz",
            f"{species}.AnimalTFDB_ortholog.gz",
            f"{species}.AnimalTFDB_2gene.tab.gz",
            f"{species}.AnimalTFDB_mapped_2disc.gz",
        ]
        for fname in expected_files:
            fpath = outdir / fname
            if fpath.exists():
                print(f"    [OK] {fname}")
            else:
                print(f"    [MISSING] {fname}")

        # Publish successful coverage to the unified TaxID-keyed registry.
        try:
            registry = SpeciesRegistry.load_default(str(self.root_dir))
            
            tf_count = len(mapped_tf_names)
            mapped_target_count = len(gene_to_tfs) if gene_to_tfs else 0
            
            registry.update_animaltfdb_stats(
                species_code=species,
                tf_count=tf_count,
                mapped_target_count=mapped_target_count,
                has_data=True
            )
            registry.save()
            print(f"[OK] Species registry updated: {species} -> AnimalTFDB")
        except Exception as e:
            print(f"[WARNING] Failed to update the species registry: {e}")

        print(f"\nAnimalTFDB Database build complete -> {outdir}")
        return str(outdir)

    # ============================
    # Build a complete species database set in one command.
    # ============================
    def build_species_db(self, species: str, taxid: int,
                         databases: List[str] = None,
                         go_version: Optional[str] = None,
                         reactome_version: Optional[str] = None,
                         **kwargs) -> str:
        """Build every requested database for one species."""
        if databases is None:
            databases = ['GO', 'Reactome']

        date_str = datetime.now().strftime("%Y%m%d")
        outdir = self.organism_dir / f"v{date_str}" / species
        outdir.mkdir(parents=True, exist_ok=True)

        print(f"\n{'#'*60}")
        print(f"# AllEnricher v2 species database (make_speciesDB)")
        print(f"# Species: {species} (taxid={taxid})")
        print(f"# Database: {', '.join(databases)}")
        print(f"# Output directory: {outdir}")
        print(f"{'#'*60}")

        built_databases = []
        skipped_databases = []
        failures = []
        for db_name in databases:
            db_upper = db_name.upper().strip()
            try:
                if db_upper == 'GO':
                    # Resolve the scientific name from project configuration or the registry.
                    latin_name = species
                    try:
                        from ..core.config import SPECIES_CONFIGS
                        if species in SPECIES_CONFIGS:
                            latin_name = SPECIES_CONFIGS[species].name
                    except:
                        pass

                    try:
                        from .species_registry import SpeciesRegistry
                        registry = SpeciesRegistry.load_default()
                        entry = registry.query_by_taxid(taxid)
                        if entry:
                            latin_name = entry.latin_name
                    except:
                        pass

                    # The GO builder selects NCBI gene2go or UniProt GOA deterministically.
                    self.build_go_with_fallback(taxid, latin_name, go_version)
                    built_databases.append(db_upper)
                elif db_upper == 'REACTOME':
                    self.build_reactome(species, taxid, reactome_version)
                    built_databases.append(db_upper)
                elif db_upper == 'KEGG':
                    self.build_kegg(species, taxid, go_version)
                    built_databases.append(db_upper)
                elif db_upper == 'DO':
                    if species.lower() == 'hsa':
                        self.build_do(taxid, go_version)
                        built_databases.append(db_upper)
                    else:
                        print("|--- [SKIPPED] Disease Ontology supports human only (hsa)")
                        skipped_databases.append(db_upper)
                elif db_upper == 'DISGENET':
                    if species.lower() == 'hsa':
                        self.build_disgenet(taxid, go_version)
                        built_databases.append(db_upper)
                    else:
                        print("|--- [SKIPPED] DisGeNET supports human only (hsa)")
                        skipped_databases.append(db_upper)
                elif db_upper == 'WIKIPATHWAYS':
                    self.build_wikipathways(species, taxid)
                    built_databases.append(db_upper)
                elif db_upper == 'TRRUST':
                    self.build_trrust(species, taxid)
                    built_databases.append(db_upper)
                elif db_upper == 'CHEA3':
                    self.build_chea3(species, taxid)
                    built_databases.append(db_upper)
                elif db_upper == 'ANIMALTFDB':
                    latin_name = kwargs.get('latin_name', '')
                    self.build_animaltfdb(species, taxid, species_latin=latin_name)
                    built_databases.append(db_upper)
                elif db_upper == 'HTFTARGET':
                    self.build_htftarget(species, taxid)
                    built_databases.append(db_upper)
                else:
                    print(f"|---[Warning] Unknown database: {db_name}")
                    skipped_databases.append(db_upper)
            except Exception as e:
                print(f"|--- [ERROR] Failed to build {db_name}: {e}")
                failures.append(f"{db_name}: {e}")

        print(f"\n{'#'*60}")
        print(f"# Species database constructed - > {outdir}")
        print(f"{'#'*60}")

        # Record build lineage and source versions for provenance reporting.
        import json as _json
        from datetime import timezone as _tz

        try:
            _downloader = DataDownloader(root_dir=str(self.root_dir))

            _build_manifest = {
                "schema_version": "1.0",
                "built_at": datetime.now(_tz.utc).isoformat(),
                "allenricher_version": __import__("allenricher").__version__,
                "species": species,
                "taxid": taxid,
                "databases": sorted(databases),
                "built_databases": sorted(built_databases),
                "skipped_databases": sorted(skipped_databases),
                "failed_databases": sorted(failures),
                "dependencies": {},
                "source_versions": {},
            }

            # Fill source versions from the download manifest when available.
            try:
                from allenricher.database.version import DatabaseVersionManager
                _vm = DatabaseVersionManager(database_dir=str(self.root_dir))

                if "GO" in built_databases:
                    _go_obo_ver = _vm.get_local_version("go_obo")
                    if _go_obo_ver and _go_obo_ver.remote_version:
                        _build_manifest["source_versions"]["go_obo"] = _go_obo_ver.remote_version
                    _gene2go_ver = _vm.get_local_version("gene2go")
                    if _gene2go_ver and _gene2go_ver.remote_last_modified:
                        _build_manifest["source_versions"]["gene2go"] = _gene2go_ver.remote_last_modified

                if "REACTOME" in built_databases:
                    _reactome_ver = _vm.get_local_version("reactome")
                    if _reactome_ver and _reactome_ver.remote_version:
                        _build_manifest["source_versions"]["reactome"] = _reactome_ver.remote_version

                if "KEGG" in built_databases:
                    _kegg_ver = _vm.get_local_version("kegg")
                    if _kegg_ver and _kegg_ver.remote_version:
                        _build_manifest["source_versions"]["kegg"] = _kegg_ver.remote_version
            except Exception as _sv_err:
                logger.warning("Failed to fill source_versions: %s", _sv_err)

            # GO Dependence
            if "GO" in built_databases:
                _go_ver = go_version if go_version else _downloader.get_latest_go_version()
                _build_manifest["dependencies"]["GO"] = {
                    "basic_dir": f"basic/go/{_go_ver}",
                    "files": ["gene2go.gz", "gene_info.gz", "go-basic.obo"],
                }

            # Reactome Dependencies
            if "REACTOME" in built_databases:
                _reactome_ver = reactome_version if reactome_version else _downloader.get_latest_reactome_version()
                _build_manifest["dependencies"]["Reactome"] = {
                    "basic_dir": f"basic/reactome/{_reactome_ver}",
                    "files": ["NCBI2Reactome_All_Levels.txt.gz", "gene_info.gz"],
                }

            # Kegg Dependence
            if "KEGG" in built_databases:
                _build_manifest["dependencies"]["KEGG"] = {
                    "source": "REST API (real-time)",
                    "gene_info_from": f"basic/go/{go_version or _downloader.get_latest_go_version()}/gene_info.gz",
                }

            # DO Dependence
            if "DO" in built_databases:
                _build_manifest["dependencies"]["DO"] = {
                    "basic_dir": "basic/do",
                    "files": ["human_disease_knowledge_filtered.tsv.gz",
                               "human_disease_experiments_filtered.tsv.gz",
                               "human_disease_textmining_filtered.tsv.gz"],
                }

            _manifest_path = outdir / "build_manifest.json"
            with open(_manifest_path, "w", encoding="utf-8") as _mf:
                _json.dump(_build_manifest, _mf, indent=2, ensure_ascii=False)
            logger.info("Build manifest saved: %s", _manifest_path)
        except Exception as _e:
            logger.warning("Failed to write the build manifest: %s", _e)

        if failures:
            raise RuntimeError("; ".join(failures))
        return str(outdir)
