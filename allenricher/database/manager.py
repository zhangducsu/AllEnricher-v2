"""Discover, validate, and load enrichment databases for one species."""

import gzip
import csv
import json
import logging
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Set

import pandas as pd

logger = logging.getLogger(__name__)

_TF_DATABASE_SPECIES = {
    "TRRUST": {"hsa", "mmu"},
    "CHEA3": {"hsa"},
    "HTFTARGET": {"hsa"},
}
_TF_DATABASES = {*_TF_DATABASE_SPECIES, "ANIMALTFDB"}


def validate_tf_database_species(database: str, species: str) -> None:
    """Reject TF databases that do not support the requested species."""
    database_name = database.upper()
    species_code = species.lower()
    supported = _TF_DATABASE_SPECIES.get(database_name)
    if supported is not None and species_code not in supported:
        raise ValueError(
            f"{database} does not support species '{species}'; supported species codes: "
            f"{', '.join(sorted(supported))}"
        )
    if database_name == "ANIMALTFDB" and species_code == "hsa":
        raise ValueError(
            "AnimalTFDB does not provide an independent human TF-target network; "
            "use hTFtarget for human analyses"
        )


# KEGG species code map for NCBI Gene Taxid
# NCBI gene_info and GOA sources are matched by NCBI TaxID.
KEGG_CODE_TO_TAXID: Dict[str, int] = {
    # Common pattern organisms
    'hsa': 9606,   # Homo sapiens (Human)
    'hpy': 835,    # Helicobacter pylori 26695
    'mta': 4530,   # Oryza sativa (Rice)
    'ath': 3702,   # Arabidopsis thaliana
    'bta': 9913,   # Bos taurus (Bovine)
    'cel': 6239,   # Caenorhabditis elegans
    'cfa': 9615,   # Canis familiaris (Dog)
    'dre': 7955,   # Danio rerio (Zebrafish)
    'dme': 7227,   # Drosophila melanogaster (Fruit fly)
    'gga': 9031,   # Gallus gallus (Chicken)
    'mcf': 594,    # Mycobacterium tuberculosis CDC1551
    'mmu': 10090,  # Mus musculus (Mouse)
    'rno': 10116,  # Rattus norvegicus (Rat)
    'sce': 4932,   # Saccharomyces cerevisiae (Yeast)
    'spo': 4896,   # Schizosaccharomyces pombe (Fission yeast)
    'xla': 8355,   # Xenopus laevis (African clawed frog)
    'xtr': 8364,   # Xenopus tropicalis
    # Other common species
    'eco': 562,    # Escherichia coli K-12
    'bsu': 224308, # Bacillus subtilis 168
    'pae': 208964, # Pseudomonas aeruginosa PAO1
    'syf': 1148,   # Synechococcus elongatus PCC 7942
    'syn': 1140,   # Synechocystis sp. PCC 6803
    'mtu': 83332,  # Mycobacterium tuberculosis H37Rv
}


class DatabaseManager:
    """Load validated gene sets and metadata for a selected species."""

    def __init__(self, database_dir: str, species: str):
        """Initialize Database Manager

        Args:
database_dir: Database directory path
Spheres: species code (e.g. hsa)
        """
        self.database_dir = Path(database_dir)
        self.species = species
        self.databases: Dict[str, Dict] = {}
        self.term_names: Dict[str, Dict[str, str]] = {}  # {db_name: {term_id: term_name}}
        self.term_hierarchies: Dict[str, Dict[str, str]] = {}
        self._active_version: Optional[str] = None
        self._active_versions: Dict[str, str] = {}

    def _find_species_dir(self, database_dir: Path, species: str, version: Optional[str] = None) -> Path:
        """Locate the canonical or legacy database directory for a species."""
        # Mode 1: database/organism/v{date}/{species}/(v2 Structure)
        organism_dir = database_dir / "organism"

        # If a version is specified, use directly
        if version and organism_dir.exists():
            species_dir = organism_dir / version / species
            if species_dir.exists():
                self._active_version = version
                return species_dir
            # List available versions when they do not exist
            available = sorted(
                d.name for d in organism_dir.iterdir()
                if d.is_dir() and (d / species).exists()
            )
            if available:
                logger.error("Version '%s\"Specific species \"%s' Does not exist. Available version: %s", version, species, ", ".join(available))
            else:
                logger.error("Species '%s' There is no version of the building.", species)

        # Auto-Find the latest version
        if organism_dir.exists():
            for version_dir in sorted(organism_dir.iterdir(), reverse=True):
                if version_dir.is_dir():
                    species_dir = version_dir / species
                    if species_dir.exists():
                        self._active_version = version_dir.name
                        return species_dir

        # v1 Compatibility
        if (database_dir / f"{species}.GO2gene.tab.gz").exists():
            self._active_version = "v1-legacy"
            return database_dir

        self._active_version = None
        return database_dir

    @property
    def active_version(self) -> Optional[str]:
        """Return the active build version for this species database."""
        versions = set(self._active_versions.values())
        if len(versions) == 1:
            return next(iter(versions))
        if len(versions) > 1:
            return "mixed"
        return self._active_version

    @property
    def database_versions(self) -> Dict[str, str]:
        """Return source versions for each loaded database."""
        return dict(self._active_versions)

    def _iter_species_dirs(self, species: str, version: Optional[str] = None) -> Iterable[Path]:
        """Yield canonical and legacy species database directories."""
        root = self.database_dir
        organism_dir = root / "organism"
        if version and organism_dir.exists():
            candidate = organism_dir / version / species
            if candidate.is_dir():
                yield candidate
            return
        if root.name == species and root.is_dir():
            yield root
            return
        if organism_dir.exists():
            for version_dir in sorted(organism_dir.iterdir(), reverse=True):
                candidate = version_dir / species
                if version_dir.is_dir() and candidate.is_dir():
                    yield candidate
        if root.is_dir():
            yield root

    @staticmethod
    def _version_for_dir(path: Path, species: str) -> str:
        if path.name == species and path.parent.parent.name == "organism":
            return path.parent.name
        return "v1-legacy"

    def _find_species_dir_with_files(
        self,
        species: str,
        filenames: Iterable[str],
        version: Optional[str] = None,
    ) -> Optional[Path]:
        filenames = list(filenames)
        for directory in self._iter_species_dirs(species, version=version):
            if any((directory / filename).is_file() for filename in filenames):
                return directory
        return None

    def _matrix_filenames(self, name: str) -> List[str]:
        prefixes = {
            'GO': ['GO'],
            'KEGG': ['kegg', 'KEGG'],
            'REACTOME': ['Reactome', 'reactome'],
            'DO': ['DO'],
            'DISGENET': ['CUI', 'DisGeNET'],
            'WIKIPATHWAYS': ['WikiPathways', 'wikipathways'],
            'TRRUST': ['gene2TF'],
            'CHEA3': ['ChEA3_2gene'],
            'ANIMALTFDB': ['AnimalTFDB_2gene'],
            'HTFTARGET': ['hTF_2gene'],
        }.get(name.upper(), [name, name.lower()])
        complete_prefixes = {'gene2TF', 'ChEA3_2gene', 'AnimalTFDB_2gene', 'hTF_2gene'}
        return list(dict.fromkeys(
            f"{self.species}.{prefix}.tab.gz"
            if prefix in complete_prefixes
            else f"{self.species}.{prefix}2gene.tab.gz"
            for prefix in prefixes
        ))

    def _gmt_filenames(self, name: str) -> List[str]:
        labels = {
            'GO': ['GO'],
            'KEGG': ['KEGG', 'kegg'],
            'REACTOME': ['Reactome', 'reactome'],
            'DO': ['DO'],
            'DISGENET': ['DisGeNET', 'DISGENET'],
            'WIKIPATHWAYS': ['WikiPathways', 'wikipathways'],
        }.get(name.upper(), [name, name.lower()])
        return list(dict.fromkeys(
            filename
            for label in labels
            for filename in (
                f"{self.species}.{label}.gmt.gz",
                f"{self.species}.{label}.gmt",
            )
        ))

    def get_build_metadata(self) -> Optional[Dict]:
        """Return build metadata for the active species database."""
        metadata: Dict = {"database_versions": self.database_versions, "source_versions": {}}
        found = False
        for version in sorted(set(self._active_versions.values())):
            if version == "v1-legacy":
                continue
            manifest_path = self.database_dir / "organism" / version / self.species / "build_manifest.json"
            if not manifest_path.exists():
                continue
            try:
                with open(manifest_path, "r", encoding="utf-8") as f:
                    manifest = json.load(f)
                metadata["source_versions"].update(manifest.get("source_versions", {}))
                if manifest.get("built_at"):
                    metadata["built_at"] = manifest["built_at"]
                found = True
            except Exception as e:
                logger.warning("Reading build_manifest. json failed: %s", e)
        return metadata if found else None

    def load_databases(self, database_names: List[str], version: Optional[str] = None) -> None:
        """Load each requested database."""
        for name in database_names:
            self.load_database(name, version=version)

    def has_database(self, name: str, version: Optional[str] = None) -> bool:
        """Return whether a valid database is available for this species."""
        try:
            validate_tf_database_species(name, self.species)
        except ValueError:
            return False
        filenames = self._gmt_filenames(name) + self._matrix_filenames(name)
        return self._find_species_dir_with_files(
            self.species, filenames, version=version
        ) is not None

    def load_database(self, name: str, version: Optional[str] = None) -> None:
        """Load one database into the normalized in-memory schema."""
        validate_tf_database_species(name, self.species)
        gmt_filenames = self._gmt_filenames(name)
        matrix_filenames = self._matrix_filenames(name)
        filenames = gmt_filenames + matrix_filenames
        database_dir = self._find_species_dir_with_files(
            self.species, filenames, version=version
        )
        if database_dir is None:
            raise FileNotFoundError(
                f"Database file does not exist: {', '.join(filenames)} not found"
            )
        active_version = self._version_for_dir(database_dir, self.species)
        self._active_versions[name] = active_version
        self._active_version = active_version
        self._validate_tf_snapshot_schema(name, database_dir)

        # Load Term Name Map (from. tab.id.gz or 2disc.gz)
        self._load_term_names(name, database_dir)

        gmt_path = next(
            (database_dir / filename for filename in gmt_filenames
             if (database_dir / filename).is_file()),
            None,
        )
        if gmt_path is not None:
            self.databases[name] = self._drop_terms_without_names(
                name, self._parse_gmt_file(gmt_path, name)
            )
            return

        filepath = next(
            (database_dir / filename for filename in matrix_filenames
             if (database_dir / filename).is_file()),
            None,
        )
        if filepath is None:
            raise FileNotFoundError(
                f"Database file does not exist: {database_dir}Not found below{', '.join(matrix_filenames)}"
            )

        # Parsing database files (the name will use term_names map at this time)
        term_data = self._parse_tab_file(filepath, name)
        self.databases[name] = self._drop_terms_without_names(name, term_data)

    @staticmethod
    def _drop_terms_without_names(db_name: str, term_data: Dict[str, Dict]) -> Dict[str, Dict]:
        """Remove non-TF terms that lack a descriptive name."""
        if db_name.upper() in _TF_DATABASES:
            return term_data
        valid = {
            term_id: info
            for term_id, info in term_data.items()
            if str(info.get("name") or "").strip().casefold() != str(term_id).strip().casefold()
        }
        removed = len(term_data) - len(valid)
        if removed:
            logger.info(
                "Database%sYes.%dThe entry is missing a specific name and has been excluded from the analysis",
                db_name, removed,
            )
        return valid

    def _validate_tf_snapshot_schema(self, name: str, database_dir: Path) -> None:
        """Reject obsolete TF snapshots that lose library or context provenance."""
        requirements = {
            "CHEA3": (f"{self.species}.ChEA3_2disc.gz", {"Term_ID", "Library"}),
            "HTFTARGET": (f"{self.species}.hTF_2disc.gz", {"Term_ID", "Context"}),
            "ANIMALTFDB": (
                f"{self.species}.AnimalTFDB_mapped_2disc.gz",
                {"Term_ID", "Context", "Inference_Type"},
            ),
        }
        requirement = requirements.get(name.upper())
        if requirement is None:
            return
        filename, required = requirement
        path = database_dir / filename
        if not path.is_file():
            raise ValueError(
                f"{name} is missing evidence metadata file {filename}; rebuild the database"
            )
        with gzip.open(path, 'rt', encoding='utf-8') as handle:
            columns = set(handle.readline().rstrip('\n').split('\t'))
        if not required.issubset(columns):
            raise ValueError(
                f"The {name} database uses a legacy metadata schema; rebuild it to preserve evidence context"
            )

    @staticmethod
    def _capitalize(text: str) -> str:
        """Capitalize display labels while preserving established abbreviations."""
        # Keep a full capitalisation map (e. g. DNA-> DNA)
        upper_map = {
            'DNA': 'DNA', 'RNA': 'RNA', 'mRNA': 'mRNA', 'tRNA': 'tRNA',
            'rRNA': 'rRNA', 'ATP': 'ATP', 'ADP': 'ADP', 'GTP': 'GTP',
            'NAD': 'NAD', 'NADH': 'NADH', 'FAD': 'FAD', 'CoA': 'CoA',
            'AMP': 'AMP', 'GMP': 'GMP', 'UMP': 'UMP', 'CMP': 'CMP',
            'MAPK': 'MAPK', 'PI3K': 'PI3K', 'AKT': 'AKT', 'EGF': 'EGF',
            'TNF': 'TNF', 'IL': 'IL', 'IFN': 'IFN', 'TGF': 'TGF',
            'VEGF': 'VEGF', 'PDGF': 'PDGF', 'FGF': 'FGF', 'IGF': 'IGF',
            'JAK': 'JAK', 'STAT': 'STAT', 'NF': 'NF', 'AP': 'AP',
            'HIF': 'HIF', 'PPAR': 'PPAR', 'RXR': 'RXR', 'LXR': 'LXR',
            'FXR': 'FXR', 'CAR': 'CAR', 'PXR': 'PXR', 'SHP': 'SHP',
            'SREBP': 'SREBP', 'PGC': 'PGC', 'cAMP': 'cAMP', 'cGMP': 'cGMP',
            'PKA': 'PKA', 'PKC': 'PKC', 'PKG': 'PKG', 'AMPk': 'AMPK',
        }
        words = text.split()
        result = []
        for word in words:
            upper_word = word.upper()
            if upper_word in upper_map:
                result.append(upper_map[upper_word])
            else:
                result.append(word.capitalize())
        return ' '.join(result)

    def _format_term_name(self, db_name: str, raw_name: str) -> str:
        """Normalize a term name and optional hierarchy for display."""
        if db_name.upper() in _TF_DATABASES:
            return raw_name

        if db_name.upper() not in {
            "GO", "KEGG", "REACTOME", "DO", "DISGENET", "WIKIPATHWAYS",
        }:
            return raw_name

        # Replace underlined with space (pathway_name_with_underserves for KEG)
        name = raw_name.replace('_', ' ')

        if db_name.upper() == 'GO':
            # GO format: "biological_process: mitochondion inheritance"
            # The following is a translation: "Biological Problem of Mitochondion Inheritance"
            if ':' in name:
                namespace, term = name.split(':', 1)
                return f"{namespace.title()}|{self._capitalize(term)}"
            return self._capitalize(name)

        elif db_name.upper() == 'KEGG':
            # KEGG format is probably:
            #   "Category SubCategory PathwayName" (Category, 3th floor)
            #   "Uncategorized Uncategorized PathwayName" (unclassified, third floor but useless)
            #   "PathwayName" (name only)
            # The output format is based on actual hierarchical numbers
            if '|' in name:
                parts = name.split('|')
                pathway_name = self._capitalize(parts[-1])

                # If the first level is classified as Uncategorized, only the pass name is shown (level 1)
                if parts[0].lower() == 'uncategorized':
                    return pathway_name

                # There is a valid classification showing Category SubCategory PathwayName (third floor)
                if len(parts) >= 3:
                    cat = self._capitalize(parts[0])
                    subcat = self._capitalize(parts[1])
                    return f"{cat}|{subcat}|{pathway_name}"
                elif len(parts) == 2:
                    return f"{parts[0].title()}|{pathway_name}"
            return self._capitalize(name)

        # Reactome, DO and all hierarchical structures, capital letters.
        return self._capitalize(name)

    @staticmethod
    def _normalize_tf_term_name(db_name: str, term_id: str, term_name: str) -> str:
        """Create a readable TF label without changing its stable identifier."""
        if db_name.upper() not in _TF_DATABASES:
            return term_name
        cleaned = str(term_name or "").strip()
        if not cleaned or cleaned.casefold() == str(term_id).strip().casefold():
            return f"{term_id} targets [{db_name}]"
        return cleaned

    def _load_term_names(self, db_name: str, database_dir: Optional[Path] = None) -> None:
        """Load term identifiers, names, and hierarchy descriptions."""
        self.term_names[db_name] = {}
        self.term_hierarchies[db_name] = {}
        database_dir = database_dir or self.database_dir

        # Map of database name to description prefix (same as true product of various builder)
        name_to_prefix = {
            'GO': 'GO',
            'KEGG': 'kegg',
            'REACTOME': 'Reactome',
            'DO': 'DO',
            'DISGENET': 'CUI',
            'WIKIPATHWAYS': 'WikiPathways',
            'TRRUST': 'TF',
            'CHEA3': 'ChEA3_',
            'HTFTARGET': 'hTF_',
        }
        prefix = name_to_prefix.get(db_name.upper(), db_name)

        # Method 1: Try Loading{species}.{db}.tab.id.gz file
        id_file = database_dir / f"{self.species}.{db_name}.tab.id.gz"
        if not id_file.exists():
            id_file = database_dir / f"{self.species}.{db_name.lower()}.tab.id.gz"

        if id_file.exists():
            try:
                with gzip.open(id_file, 'rt', encoding='utf-8') as f:
                    for line in f:
                        line = line.strip()
                        if not line:
                            continue
                        # Format: TermID\tTermName\tParentTerms
                        parts = line.split('\t')
                        if len(parts) >= 2:
                            term_id = parts[0]
                            term_name = (
                                term_id if db_name.upper() in _TF_DATABASES else parts[1]
                            )
                            self.term_names[db_name][term_id] = self._normalize_tf_term_name(
                                db_name, term_id, term_name
                            )
                            if len(parts) >= 3 and '|' in parts[2]:
                                self.term_hierarchies[db_name][term_id] = parts[2]
                return  # If loaded successfully, return directly
            except Exception as e:
                print(f"Warning: Failed to load {id_file}: {e}")

        # Option 2: Try loading from *.2disc.gz.
        special_disc = {
            'ANIMALTFDB': [
                f"{self.species}.AnimalTFDB_mapped_2disc.gz",
                f"{self.species}.AnimalTFDB_2disc.gz",
            ],
        }.get(db_name.upper(), [])
        disc_candidates = special_disc + [
            f"{self.species}.{prefix}2disc.gz",
            f"{self.species}.{prefix.lower()}2disc.gz",
            f"{self.species}.{db_name}2disc.gz",
            f"{self.species}.{db_name.lower()}2disc.gz",
            f"{prefix}2disc.gz",
            f"{prefix.lower()}2disc.gz",
        ]
        disc_file = next(
            (database_dir / filename for filename in dict.fromkeys(disc_candidates)
             if (database_dir / filename).is_file()),
            database_dir / disc_candidates[0],
        )

        if disc_file.exists():
            try:
                count = 0
                has_named_schema = False
                with gzip.open(disc_file, 'rt', encoding='utf-8') as f:
                    for line in f:
                        line = line.strip()
                        if not line:
                            continue
                        # Format: TermID\tTermName (Term Name with a stony partition level)
                        parts = line.split('\t')
                        if len(parts) >= 2:
                            term_id = parts[0]
                            raw_name = parts[1]
                            if term_id == "Term_ID":
                                has_named_schema = raw_name == "Term_Name"
                                continue
                            if term_id == "TF":
                                continue
                            if db_name.upper() in _TF_DATABASES and not has_named_schema:
                                # Legacy TF descriptions used the second column for counts/modes,
                                # not a display name. Keep the TF identifier as the label.
                                term_name = term_id
                            else:
                                term_name = self._format_term_name(db_name, raw_name)
                            term_name = self._normalize_tf_term_name(
                                db_name, term_id, term_name
                            )
                            self.term_names[db_name][term_id] = term_name
                            hierarchy = parts[2].strip() if len(parts) >= 3 else ""
                            if '|' not in hierarchy and '|' in term_name:
                                hierarchy = term_name
                            if '|' in hierarchy:
                                self.term_hierarchies[db_name][term_id] = hierarchy
                            count += 1
                if count > 0:
                    print(f"Loaded {count} disease terms from {disc_file.name}")
            except Exception as e:
                print(f"Warning: failed to load {disc_file}: {e}")

    def _parse_tab_file(self, filepath: Path, db_name: str) -> Dict[str, Dict]:
        """Parse a compressed membership matrix and its descriptions."""
        term_data: Dict[str, Dict] = {}

        # Map to get the name of the database
        name_map = self.term_names.get(db_name, {})
        hierarchy_map = self.term_hierarchies.get(db_name, {})

        with gzip.open(filepath, 'rt', encoding='utf-8') as f:
            reader = csv.reader(f, delimiter='\t')

            # Read Table Headers
            header = next(reader)
            # Header format: Gene\tTermID1\tTermID2\t...
            term_ids = header[1:]  # The first column is Gene, followed by Term ID.

            # Initialize term_data, using name in term_names or Term_ID as name
            for term_id in term_ids:
                # If a name is map, use the name of the map; otherwise, use Term ID
                term_name = name_map.get(term_id, term_id)
                term_data[term_id] = {"name": term_name, "genes": []}
                if term_id in hierarchy_map:
                    term_data[term_id]["hierarchy"] = hierarchy_map[term_id]

            # Read Data Line
            for row in reader:
                if len(row) < 2:
                    continue

                gene = row[0]
                # The column that follows indicates whether the gene belongs to the corresponding entry (1 = yes, 0 = not)
                for i, value in enumerate(row[1:]):
                    if i < len(term_ids) and value == '1':
                        term_id = term_ids[i]
                        term_data[term_id]["genes"].append(gene)

        return term_data

    def _parse_gmt_file(self, filepath: Path, db_name: str) -> Dict[str, Dict]:
        """Parse GMT records into normalized term metadata and genes."""
        term_data: Dict[str, Dict] = {}
        name_map = self.term_names.get(db_name, {})
        hierarchy_map = self.term_hierarchies.get(db_name, {})
        opener = gzip.open if filepath.suffix.lower() == '.gz' else open
        with opener(filepath, 'rt', encoding='utf-8') as handle:
            for line_number, line in enumerate(handle, 1):
                parts = line.rstrip('\r\n').split('\t')
                if len(parts) < 3 or not parts[0].strip():
                    logger.warning("Skip invalid GMT lines%s: %d", filepath, line_number)
                    continue
                term_id = parts[0].strip()
                description = parts[1].strip()
                genes = [gene.strip() for gene in parts[2:] if gene.strip()]
                if not genes:
                    continue
                if term_id not in term_data:
                    mapped_name = name_map.get(term_id)
                    if mapped_name:
                        term_name = mapped_name
                    else:
                        raw_name = term_id if db_name.upper() in _TF_DATABASES else description or term_id
                        term_name = self._format_term_name(db_name, raw_name)
                    term_data[term_id] = {
                        "name": term_name,
                        "genes": [],
                    }
                    if term_id in hierarchy_map:
                        term_data[term_id]["hierarchy"] = hierarchy_map[term_id]
                term_data[term_id]["genes"] = list(dict.fromkeys(
                    [*term_data[term_id]["genes"], *genes]
                ))
        return term_data

    def get_all_term_data(self) -> Dict[str, Dict]:
        """Return normalized term data for all loaded databases."""
        return self.databases

    def get_background_genes(self) -> Set[str]:
        """Return the union of genes annotated by loaded databases."""
        background = set()
        for db_data in self.databases.values():
            for term_info in db_data.values():
                background.update(term_info["genes"])
        return background

    def get_database_genes(self, db_name: str) -> Set[str]:
        """Return all genes annotated by one database."""
        if db_name not in self.databases:
            return set()

        genes = set()
        for term_info in self.databases[db_name].values():
            genes.update(term_info["genes"])
        return genes

    def get_genome_genes(self, taxid: Optional[int] = None, species_code: Optional[str] = None) -> Set[str]:
        """Return all gene symbols recorded for the selected TaxID."""
        # Resolve the NCBI TaxID from the species code when available.
        if taxid is None and species_code is not None:
            taxid = KEGG_CODE_TO_TAXID.get(species_code.lower())

        if taxid is None:
            # If you still don't have a taxid, go back to the air.
            return set()

        genome_file = self.database_dir / "gene_info.gz"
        if not genome_file.exists():
            db_root = next(
                (path for path in (self.database_dir, *self.database_dir.parents)
                 if path.name == "database"),
                self.database_dir,
            )
            basic_go_dir = db_root / "basic" / "go"
            candidates = (
                sorted(basic_go_dir.glob("*/gene_info.gz"), reverse=True)
                if basic_go_dir.exists() else []
            )
            if not candidates:
                return set()
            genome_file = candidates[0]

        genes = set()
        # Keep gene records; exclude biological regions, pseudogenes, and other non-gene entries.
        _VALID_GENE_TYPES = frozenset({
            "protein-coding", "ncRNA", "snoRNA", "rRNA", "tRNA", "snRNA",
            "scRNA", "other",
        })
        with gzip.open(genome_file, "rt", encoding="utf-8") as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                parts = line.split("\t")
                if len(parts) >= 10:
                    file_taxid = parts[0]
                    gene_symbol = parts[2]  # Column 3: Gene Symbol (e.g. A1BG)
                    gene_type = parts[9]    # Column 10: gene type
                    # Retain valid gene records for the requested NCBI TaxID.
                    try:
                        if int(file_taxid) == taxid and gene_type in _VALID_GENE_TYPES:
                            if gene_symbol and gene_symbol != "-" and not gene_symbol.startswith("NEW|"):
                                genes.add(gene_symbol)
                    except ValueError:
                        continue
        return genes

    @staticmethod
    def _read_tf_info_file(path: Path, legacy_columns: List[str]) -> pd.DataFrame:
        with gzip.open(path, 'rt', encoding='utf-8') as handle:
            first_cell = handle.readline().split('\t', 1)[0]
        if first_cell in {"Term_ID", "TF"}:
            return pd.read_csv(path, sep='\t', compression='gzip', dtype=str)
        return pd.read_csv(
            path, sep='\t', compression='gzip', header=None,
            names=legacy_columns, dtype=str,
        )

    def load_trrust(self, species: Optional[str] = None) -> Optional[Dict[str, 'pd.DataFrame']]:
        """Load TRRUST TF-target gene sets."""
        sp = species or self.species
        validate_tf_database_species("TRRUST", sp)
        filenames = [f"{sp}.TF2target.tab.gz", f"{sp}.gene2TF.tab.gz"]
        base_dir = self._find_species_dir_with_files(sp, filenames)
        if base_dir is None:
            logger.warning("The TRRUST database file does not exist")
            return None

        tf2target_file = base_dir / filenames[0]
        gene2tf_file = base_dir / filenames[1]
        tf2disc_file = base_dir / f"{sp}.TF2disc.gz"

        # Check whether the necessary documents exist
        if not tf2target_file.exists() and not gene2tf_file.exists():
            logger.warning("The TRRUST database file does not exist: %s", base_dir)
            return None

        result: Dict[str, 'pd.DataFrame'] = {}

        if tf2target_file.exists():
            result['tf2target'] = pd.read_csv(tf2target_file, sep='\t', compression='gzip')
        else:
            result['tf2target'] = pd.DataFrame()

        if gene2tf_file.exists():
            result['gene2tf'] = pd.read_csv(gene2tf_file, sep='\t', compression='gzip')
        else:
            result['gene2tf'] = pd.DataFrame()

        if tf2disc_file.exists():
            result['tf_info'] = self._read_tf_info_file(
                tf2disc_file, ["TF", "Mode", "Target_Set_Size"]
            )
        else:
            result['tf_info'] = pd.DataFrame()

        edge_file = base_dir / f"{sp}.TRRUST_edges.tsv.gz"
        result['edges'] = (
            pd.read_csv(edge_file, sep='\t', compression='gzip', dtype=str)
            if edge_file.exists() else pd.DataFrame()
        )

        return result

    def load_chea3(self, species: Optional[str] = None) -> Optional[Dict[str, 'pd.DataFrame']]:
        """Load ChEA3 TF-target gene sets."""
        sp = species or self.species
        validate_tf_database_species("ChEA3", sp)
        gene2tf_name = f"{sp}.ChEA3_2gene.tab.gz"
        base_dir = self._find_species_dir_with_files(sp, [gene2tf_name])
        if base_dir is None:
            logger.warning("ChEA3 database file does not exist")
            return None

        gene2tf_file = base_dir / gene2tf_name
        tf2disc_file = base_dir / f"{sp}.ChEA3_2disc.gz"

        # Check whether the necessary documents exist
        if not gene2tf_file.exists():
            logger.warning("The ChEA3 database file does not exist: %s", base_dir)
            return None

        result: Dict[str, 'pd.DataFrame'] = {}

        result['gene2tf'] = pd.read_csv(gene2tf_file, sep='\t', compression='gzip')

        if tf2disc_file.exists():
            result['tf_info'] = self._read_tf_info_file(
                tf2disc_file, ["TF", "lib_count", "Target_Set_Size"]
            )
        else:
            result['tf_info'] = pd.DataFrame()

        if "Term_ID" not in result['tf_info'].columns:
            raise ValueError("ChEA3 uses a legacy metadata schema; rebuild the database")

        return result

    def load_htftarget(self, species: Optional[str] = None) -> Optional[Dict[str, pd.DataFrame]]:
        """Load hTFtarget regulatory-context gene sets."""
        sp = species or self.species
        validate_tf_database_species("hTFtarget", sp)
        gene2tf_name = f"{sp}.hTF_2gene.tab.gz"
        db_dir = self._find_species_dir_with_files(sp, [gene2tf_name])
        if db_dir is None:
            return None

        gene2tf_file = db_dir / gene2tf_name
        disc_file = db_dir / f"{sp}.hTF_2disc.gz"

        if not gene2tf_file.exists():
            return None

        result = {}
        result['gene2tf'] = pd.read_csv(gene2tf_file, sep='\t', compression='gzip', low_memory=False)

        if disc_file.exists():
            result['tf_info'] = self._read_tf_info_file(
                disc_file, ["TF", "target_count", "tissues", "source"]
            )
        else:
            result['tf_info'] = pd.DataFrame()

        if "Term_ID" not in result['tf_info'].columns:
            raise ValueError("hTFtarget uses a legacy metadata schema; rebuild the database")

        return result

    def load_animaltfdb(self, species: Optional[str] = None) -> Optional[Dict[str, pd.DataFrame]]:
        """Load species-specific AnimalTFDB TF-target gene sets."""
        sp = species or self.species
        validate_tf_database_species("AnimalTFDB", sp)
        gene2tf_name = f"{sp}.AnimalTFDB_2gene.tab.gz"
        db_dir = self._find_species_dir_with_files(sp, [gene2tf_name])
        if db_dir is None:
            return None

        gene2tf_file = db_dir / gene2tf_name
        disc_file = db_dir / f"{sp}.AnimalTFDB_mapped_2disc.gz"

        if not gene2tf_file.exists():
            return None

        result = {}
        result['gene2tf'] = pd.read_csv(gene2tf_file, sep='\t', compression='gzip', low_memory=False)

        if disc_file.exists():
            result['tf_info'] = self._read_tf_info_file(
                disc_file, ["TF", "target_count", "Family", "source"]
            )
        else:
            result['tf_info'] = pd.DataFrame()

        if not {"Term_ID", "Inference_Type"}.issubset(result['tf_info'].columns):
            raise ValueError("AnimalTFDB uses a legacy metadata schema; rebuild the database")

        return result
