"""Read, query, merge, and write the unified species coverage registry."""

from __future__ import annotations

import csv
import difflib
import logging
from dataclasses import dataclass, asdict
from pathlib import Path
from typing import Optional, Dict, List, Any, Tuple


logger = logging.getLogger(__name__)


# TSV column definition, consistent with file format
_FIELD_NAMES: List[str] = [
    "taxid", "latin_name", "common_name",
    "has_go", "go_source", "go_filename", "go_file_size", "go_gene_count", "go_term_count",
    "has_kegg", "kegg_code", "kegg_code_source", "kegg_gene_count", "kegg_pathway_count",
    "has_reactome", "reactome_code", "reactome_gene_count", "reactome_pathway_count",
    "has_do", "do_gene_count", "do_term_count",
    "has_disgenet", "disgenet_gene_count", "disgenet_term_count",
    "has_wikipathways", "wikipathways_data_type", "wikipathways_gene_count", "wikipathways_pathway_count",
    "has_trrust", "trrust_tf_count", "trrust_target_count",
    "has_chea3", "chea3_tf_count", "chea3_target_count",
    "has_animaltfdb", "animaltfdb_tf_count", "animaltfdb_mapped_target_count",
    "has_htftarget", "htftarget_tf_count", "htftarget_target_count",
    "synonyms",
]

@dataclass
class SpeciesEntry:
    """Represent one TaxID-keyed row in the unified species registry."""
    taxid: int
    latin_name: str
    common_name: Optional[str] = None
    # GO related field
    has_go: bool = False
    go_source: Optional[str] = None
    go_filename: Optional[str] = None
    go_file_size: Optional[int] = None
    go_gene_count: Optional[int] = None
    go_term_count: Optional[int] = None
    # KEGG-related fields
    has_kegg: bool = False
    kegg_code: Optional[str] = None
    kegg_code_source: Optional[str] = None
    kegg_gene_count: Optional[int] = None
    kegg_pathway_count: Optional[int] = None
    # Reactome Related Fields
    has_reactome: bool = False
    reactome_code: Optional[str] = None
    reactome_gene_count: Optional[int] = None
    reactome_pathway_count: Optional[int] = None
    # DO related fields
    has_do: bool = False
    do_gene_count: Optional[int] = None
    do_term_count: Optional[int] = None
    # DisGeNET-related fields
    has_disgenet: bool = False
    disgenet_gene_count: Optional[int] = None
    disgenet_term_count: Optional[int] = None
    # WikiPathways related fields
    has_wikipathways: bool = False
    wikipathways_data_type: Optional[str] = None  # 'gmt', 'gpml', or None
    wikipathways_gene_count: Optional[int] = None
    wikipathways_pathway_count: Optional[int] = None
    # TRUST-related fields
    has_trrust: bool = False
    trrust_tf_count: Optional[int] = None
    trrust_target_count: Optional[int] = None
    # ChEA3 related fields
    has_chea3: bool = False
    chea3_tf_count: Optional[int] = None
    chea3_target_count: Optional[int] = None
    # AnimalTFDB-related field
    has_animaltfdb: bool = False
    animaltfdb_tf_count: Optional[int] = None
    animaltfdb_mapped_target_count: Optional[int] = None
    # hTFtarget related field
    has_htftarget: bool = False
    htftarget_tf_count: Optional[int] = None
    htftarget_target_count: Optional[int] = None
    synonyms: Optional[str] = None  # All retraceable aliases, semicolons separated


class SpeciesRegistry:
    """Query and maintain database coverage keyed by NCBI TaxID."""

    def __init__(self, registry_path: Path) -> None:
        """Initialize a species registry.

        Args:
            registry_path: Path to the species registry TSV file.
        """
        self.registry_path = Path(registry_path)
        self.entries: Dict[int, SpeciesEntry] = {}
        self._synonym_index: Dict[str, int] = {}  # Name (lower) *taxid

    def load(self) -> None:
        """Load a species registry from TSV."""
        self.entries.clear()

        if not self.registry_path.exists():
            return

        with open(self.registry_path, "r", encoding="utf-8") as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            for row in reader:
                entry = self._parse_row(row)
                if entry is not None:
                    self.entries[entry.taxid] = entry

        self._load_synonyms_from_names_dmp()

    def save(self) -> None:
        """Write the species registry as deterministic TSV."""
        self.registry_path.parent.mkdir(parents=True, exist_ok=True)

        field_names = _FIELD_NAMES

        with open(self.registry_path, "w", encoding="utf-8", newline="") as fh:
            writer = csv.DictWriter(fh, fieldnames=field_names, delimiter="\t")
            writer.writeheader()

            for taxid in sorted(self.entries):
                writer.writerow(self._format_row(self.entries[taxid]))

    def add_entry(self, entry: SpeciesEntry) -> None:
        """Add or replace one TaxID-keyed species record."""
        self.entries[entry.taxid] = entry

    def query_by_taxid(self, taxid: int) -> Optional[SpeciesEntry]:
        """Return the exact NCBI TaxID match."""
        return self.entries.get(taxid)

    def query_by_latin_name(self, name: str) -> List[SpeciesEntry]:
        """Find species by case-insensitive scientific name."""
        keyword = name.strip().lower()
        if not keyword:
            return []
        return [
            entry for entry in self.entries.values()
            if keyword in entry.latin_name.lower()
        ]

    def query_by_kegg_code(self, code: str) -> Optional[SpeciesEntry]:
        """Return the exact KEGG organism-code match."""
        code_normalized = code.strip().lower()
        matches: List[SpeciesEntry] = []
        for entry in self.entries.values():
            if entry.has_kegg and entry.kegg_code is not None and entry.kegg_code.lower() == code_normalized:
                matches.append(entry)

        if not matches:
            return None
        if len(matches) == 1:
            return matches[0]

        # Rank ambiguous matches by database coverage, then prefer a full
        # scientific name over an abbreviated genus.
        def _completeness_score(e: SpeciesEntry) -> int:
            return sum([int(e.has_go), int(e.has_kegg), int(e.has_reactome), int(e.has_do)])
        def _name_fullness(e: SpeciesEntry) -> int:
            parts = e.latin_name.split()
            # A full genus name outranks an abbreviated genus.
            genus_full = 2 if len(parts) > 0 and len(parts[0]) > 1 else 1
            return genus_full * 100 + len(parts)
        matches.sort(key=lambda e: (_completeness_score(e), _name_fullness(e)), reverse=True)
        logger.info(
            "KEGG code '%s' matched %d registry entries; selected TaxID %d (%s) "
            "with database coverage score %d",
            code,
            len(matches),
            matches[0].taxid,
            matches[0].latin_name,
            _completeness_score(matches[0]),
        )
        return matches[0]

    def fuzzy_search(self, query: str, cutoff: float = 0.6) -> List[Tuple[SpeciesEntry, float, str]]:
        """Search scientific names, common names, synonyms, and identifiers."""
        query_lower = query.strip().lower()
        if not query_lower:
            return []

        results: List[Tuple[SpeciesEntry, float, str]] = []
        seen_taxids = set()

        # 0 Synonym/old name map matching (highest priority, second only to accurate matching)
        synonym_matches = self._check_synonyms(query_lower)
        for entry, syn_name in synonym_matches:
            if entry.taxid not in seen_taxids:
                results.append((entry, 0.95, f'synonym:{syn_name}'))
                seen_taxids.add(entry.taxid)

        # 1. Precise matching
        for entry in self.entries.values():
            if entry.taxid in seen_taxids:
                continue
            if entry.latin_name.lower() == query_lower:
                results.append((entry, 1.0, 'exact'))
                seen_taxids.add(entry.taxid)

        # 2. Subsequent matching
        for entry in self.entries.values():
            if entry.taxid in seen_taxids:
                continue
            if query_lower in entry.latin_name.lower():
                score = len(query_lower) / len(entry.latin_name.lower())
                results.append((entry, score, 'substring'))
                seen_taxids.add(entry.taxid)

        # 3. Fuzzy Match (using difflib.SequienceMatcher)
        for entry in self.entries.values():
            if entry.taxid in seen_taxids:
                continue
            similarity = difflib.SequenceMatcher(None, query_lower, entry.latin_name.lower()).ratio()
            if similarity >= cutoff:
                results.append((entry, similarity, 'fuzzy'))
                seen_taxids.add(entry.taxid)

        # Sort by fraction
        results.sort(key=lambda x: x[1], reverse=True)
        return results

    _NAME_CLASSES_FOR_SEARCH = {
        "scientific name", "synonym", "common name",
        "genbank common name", "blast name", "equivalent name", "acronym",
    }

    def _load_synonyms_from_names_dmp(self) -> None:
        """Build TaxID-keyed synonyms from NCBI taxonomy names.dmp."""
        self._synonym_index.clear()

        # Determine names.dmp path
        candidates = [
            self.registry_path.parent / "basic" / "taxonomy" / "names.dmp",
            self.registry_path.parent.parent / "basic" / "taxonomy" / "names.dmp",
        ]
        names_dmp = None
        for p in candidates:
            if p.exists():
                names_dmp = p
                break

        if names_dmp is None:
            return

        target_taxids = set(self.entries.keys())
        all_synonyms: Dict[int, List[str]] = {}

        with open(names_dmp, "r", encoding="utf-8") as fh:
            for line in fh:
                parts = line.rstrip("\n").split("\t|\t")
                if len(parts) < 4:
                    continue
                taxid_str, name, _, name_class = parts[0], parts[1], parts[2], parts[3].strip().rstrip("\t|").strip()
                try:
                    taxid = int(taxid_str)
                except ValueError:
                    continue

                if taxid not in target_taxids:
                    continue

                if name_class not in self._NAME_CLASSES_FOR_SEARCH:
                    continue

                name = name.strip()
                if not name:
                    continue

                name_lower = name.lower()
                self._synonym_index[name_lower] = taxid

                all_synonyms.setdefault(taxid, []).append(name)

        for taxid, names in all_synonyms.items():
            entry = self.entries.get(taxid)
            if entry:
                unique_names = list(dict.fromkeys(names))
                entry.synonyms = ";".join(unique_names)

    def _check_synonyms(self, query_lower: str) -> List[Tuple[SpeciesEntry, str]]:
        """Return whether a query matches a recorded taxonomy synonym."""
        results = []
        taxid = self._synonym_index.get(query_lower)
        if taxid is not None:
            entry = self.entries.get(taxid)
            if entry:
                results.append((entry, query_lower))
        return results

    def filter_by_databases(
        self,
        go: Optional[bool] = None,
        kegg: Optional[bool] = None,
        reactome: Optional[bool] = None,
        do: Optional[bool] = None,
        disgenet: Optional[bool] = None,
        wikipathways: Optional[bool] = None,
        trrust: Optional[bool] = None,
        chea3: Optional[bool] = None,
        animaltfdb: Optional[bool] = None,
        htftarget: Optional[bool] = None,
    ) -> List[SpeciesEntry]:
        """Return species covered by every requested database."""
        results: List[SpeciesEntry] = []
        for entry in self.entries.values():
            if go is not None and entry.has_go != go:
                continue
            if kegg is not None and entry.has_kegg != kegg:
                continue
            if reactome is not None and entry.has_reactome != reactome:
                continue
            if do is not None and entry.has_do != do:
                continue
            if disgenet is not None and entry.has_disgenet != disgenet:
                continue
            if wikipathways is not None and entry.has_wikipathways != wikipathways:
                continue
            if trrust is not None and entry.has_trrust != trrust:
                continue
            if chea3 is not None and entry.has_chea3 != chea3:
                continue
            if animaltfdb is not None and entry.has_animaltfdb != animaltfdb:
                continue
            if htftarget is not None and entry.has_htftarget != htftarget:
                continue
            results.append(entry)
        return results

    def get_summary(self) -> Dict[str, Any]:
        """Summarize TaxID-keyed species coverage by database."""
        go_count = 0
        go_with_genes = 0
        go_with_terms = 0
        kegg_count = 0
        kegg_with_genes = 0
        kegg_with_pathways = 0
        reactome_count = 0
        reactome_with_genes = 0
        reactome_with_pathways = 0
        do_count = 0
        do_with_genes = 0
        do_with_terms = 0
        disgenet_count = 0
        disgenet_with_genes = 0
        disgenet_with_terms = 0
        wikipathways_count = 0
        wikipathways_with_genes = 0
        wikipathways_with_pathways = 0
        trrust_count = 0
        trrust_with_tfs = 0
        trrust_with_targets = 0
        chea3_count = 0
        chea3_with_tfs = 0
        chea3_with_targets = 0
        animaltfdb_count = 0
        animaltfdb_with_tfs = 0
        animaltfdb_with_targets = 0
        htftarget_count = 0
        htftarget_with_tfs = 0
        htftarget_with_targets = 0

        for entry in self.entries.values():
            if entry.has_go:
                go_count += 1
                if entry.go_gene_count is not None:
                    go_with_genes += 1
                if entry.go_term_count is not None:
                    go_with_terms += 1
            if entry.has_kegg:
                kegg_count += 1
                if entry.kegg_gene_count is not None:
                    kegg_with_genes += 1
                if entry.kegg_pathway_count is not None:
                    kegg_with_pathways += 1
            if entry.has_reactome:
                reactome_count += 1
                if entry.reactome_gene_count is not None:
                    reactome_with_genes += 1
                if entry.reactome_pathway_count is not None:
                    reactome_with_pathways += 1
            if entry.has_do:
                do_count += 1
                if entry.do_gene_count is not None:
                    do_with_genes += 1
                if entry.do_term_count is not None:
                    do_with_terms += 1
            if entry.has_disgenet:
                disgenet_count += 1
                if entry.disgenet_gene_count is not None:
                    disgenet_with_genes += 1
                if entry.disgenet_term_count is not None:
                    disgenet_with_terms += 1
            if entry.has_wikipathways:
                wikipathways_count += 1
                if entry.wikipathways_gene_count is not None:
                    wikipathways_with_genes += 1
                if entry.wikipathways_pathway_count is not None:
                    wikipathways_with_pathways += 1
            if entry.has_trrust:
                trrust_count += 1
                if entry.trrust_tf_count is not None:
                    trrust_with_tfs += 1
                if entry.trrust_target_count is not None:
                    trrust_with_targets += 1
            if entry.has_chea3:
                chea3_count += 1
                if entry.chea3_tf_count is not None:
                    chea3_with_tfs += 1
                if entry.chea3_target_count is not None:
                    chea3_with_targets += 1
            if entry.has_animaltfdb and entry.taxid != 9606:
                animaltfdb_count += 1
                if entry.animaltfdb_tf_count is not None:
                    animaltfdb_with_tfs += 1
                if entry.animaltfdb_mapped_target_count is not None:
                    animaltfdb_with_targets += 1
            if entry.has_htftarget:
                htftarget_count += 1
                if entry.htftarget_tf_count is not None:
                    htftarget_with_tfs += 1
                if entry.htftarget_target_count is not None:
                    htftarget_with_targets += 1

        return {
            "total_species": len(self.entries),
            "go": {
                "count": go_count,
                "with_gene_count": go_with_genes,
                "with_term_count": go_with_terms,
            },
            "kegg": {
                "count": kegg_count,
                "with_gene_count": kegg_with_genes,
                "with_pathway_count": kegg_with_pathways,
            },
            "reactome": {
                "count": reactome_count,
                "with_gene_count": reactome_with_genes,
                "with_pathway_count": reactome_with_pathways,
            },
            "do": {
                "count": do_count,
                "with_gene_count": do_with_genes,
                "with_term_count": do_with_terms,
            },
            "disgenet": {
                "count": disgenet_count,
                "with_gene_count": disgenet_with_genes,
                "with_term_count": disgenet_with_terms,
            },
            "wikipathways": {
                "count": wikipathways_count,
                "with_gene_count": wikipathways_with_genes,
                "with_pathway_count": wikipathways_with_pathways,
            },
            "trrust": {
                "count": trrust_count,
                "with_tf_count": trrust_with_tfs,
                "with_target_count": trrust_with_targets,
            },
            "chea3": {
                "count": chea3_count,
                "with_tf_count": chea3_with_tfs,
                "with_target_count": chea3_with_targets,
            },
            "animaltfdb": {
                "count": animaltfdb_count,
                "with_tf_count": animaltfdb_with_tfs,
                "with_target_count": animaltfdb_with_targets,
            },
            "htftarget": {
                "count": htftarget_count,
                "with_tf_count": htftarget_with_tfs,
                "with_target_count": htftarget_with_targets,
            },
        }

    def get_species_detail(self, taxid: int) -> Optional[Dict[str, Any]]:
        """Return one species record as a serializable dictionary."""
        entry = self.entries.get(taxid)
        if entry is None:
            return None
        return asdict(entry)

    @staticmethod
    def generate_kegg_abbreviation(latin_name: str) -> str:
        """Derive a candidate KEGG code from a scientific name."""
        parts = latin_name.strip().split()
        if len(parts) < 2:
            # Only by name, first 3 letters
            genus = parts[0].lower()[:3]
            return genus.ljust(3, "x")[:3]
        genus_initial = parts[0].lower()[0]
        species_prefix = parts[1].lower()[:2]
        return genus_initial + species_prefix

    @classmethod
    def load_default(cls, root_dir: str = "./database") -> "SpeciesRegistry":
        """Load the unified registry from the default data directory."""
        root = Path(root_dir)
        current_path = root / "basic" / "supported_species.tsv"
        registry = cls(
            registry_path=current_path if current_path.exists() else root / "supported_species.tsv"
        )
        registry.load()
        return registry

    def update_animaltfdb_stats(self, species_code: str, tf_count: int,
                                mapped_target_count: int, has_data: bool = True) -> None:
        """Refresh AnimalTFDB coverage statistics in the registry."""
        entry = self.query_by_kegg_code(species_code)
        if entry:
            entry.has_animaltfdb = has_data
            entry.animaltfdb_tf_count = tf_count
            entry.animaltfdb_mapped_target_count = mapped_target_count
            logger.info(f"Update of species register: {species_code} - AnimalTFDB TF={tf_count}, targets={mapped_target_count}")
        else:
            logger.warning("Species %s was not found in the registry", species_code)

    # ------------------------------------------------------------------
    # Internal support methods
    # ------------------------------------------------------------------

    @staticmethod
    def _parse_bool(value: str) -> bool:
        """Parse a serialized boolean value."""
        return value.strip() == "True"

    @staticmethod
    def _parse_optional_int(value: str) -> Optional[int]:
        """Parse an optional integer value."""
        stripped = value.strip()
        if stripped in ("", "-"):
            return None
        return int(stripped)

    @staticmethod
    def _parse_optional_str(value: str) -> Optional[str]:
        """Parse an optional non-empty string."""
        stripped = value.strip()
        if stripped in ("", "-"):
            return None
        return stripped

    def _parse_row(self, row: Dict[str, str]) -> Optional[SpeciesEntry]:
        """Convert one TSV row into a species entry."""
        try:
            taxid = int(row["taxid"].strip())
        except (ValueError, KeyError):
            return None

        return SpeciesEntry(
            taxid=taxid,
            latin_name=row["latin_name"].strip(),
            common_name=self._parse_optional_str(row.get("common_name", "")),
            has_go=self._parse_bool(row.get("has_go", "False")),
            go_source=self._parse_optional_str(row.get("go_source", "")),
            go_filename=self._parse_optional_str(row.get("go_filename", "")),
            go_file_size=self._parse_optional_int(row.get("go_file_size", "")),
            go_gene_count=self._parse_optional_int(row.get("go_gene_count", "")),
            go_term_count=self._parse_optional_int(row.get("go_term_count", "")),
            has_kegg=self._parse_bool(row.get("has_kegg", "False")),
            kegg_code=self._parse_optional_str(row.get("kegg_code", "")),
            kegg_code_source=self._parse_optional_str(row.get("kegg_code_source", "")),
            kegg_gene_count=self._parse_optional_int(row.get("kegg_gene_count", "")),
            kegg_pathway_count=self._parse_optional_int(row.get("kegg_pathway_count", "")),
            has_reactome=self._parse_bool(row.get("has_reactome", "False")),
            reactome_code=self._parse_optional_str(row.get("reactome_code", "")),
            reactome_gene_count=self._parse_optional_int(row.get("reactome_gene_count", "")),
            reactome_pathway_count=self._parse_optional_int(row.get("reactome_pathway_count", "")),
            has_do=self._parse_bool(row.get("has_do", "False")),
            do_gene_count=self._parse_optional_int(row.get("do_gene_count", "")),
            do_term_count=self._parse_optional_int(row.get("do_term_count", "")),
            has_disgenet=self._parse_bool(row.get("has_disgenet", "False")),
            disgenet_gene_count=self._parse_optional_int(row.get("disgenet_gene_count", "")),
            disgenet_term_count=self._parse_optional_int(row.get("disgenet_term_count", "")),
            has_wikipathways=self._parse_bool(row.get("has_wikipathways", "False")),
            wikipathways_data_type=self._parse_optional_str(row.get("wikipathways_data_type", "")),
            wikipathways_gene_count=self._parse_optional_int(row.get("wikipathways_gene_count", "")),
            wikipathways_pathway_count=self._parse_optional_int(row.get("wikipathways_pathway_count", "")),
            has_trrust=self._parse_bool(row.get("has_trrust", "False")),
            trrust_tf_count=self._parse_optional_int(row.get("trrust_tf_count", "")),
            trrust_target_count=self._parse_optional_int(row.get("trrust_target_count", "")),
            has_chea3=self._parse_bool(row.get("has_chea3", "False")),
            chea3_tf_count=self._parse_optional_int(row.get("chea3_tf_count", "")),
            chea3_target_count=self._parse_optional_int(row.get("chea3_target_count", "")),
            has_animaltfdb=self._parse_bool(row.get("has_animaltfdb", "False")),
            animaltfdb_tf_count=self._parse_optional_int(row.get("animaltfdb_tf_count", "")),
            animaltfdb_mapped_target_count=self._parse_optional_int(row.get("animaltfdb_mapped_target_count", "")),
            has_htftarget=self._parse_bool(row.get("has_htftarget", "False")),
            htftarget_tf_count=self._parse_optional_int(row.get("htftarget_tf_count", "")),
            htftarget_target_count=self._parse_optional_int(row.get("htftarget_target_count", "")),
            synonyms=self._parse_optional_str(row.get("synonyms", "")),
        )

    @staticmethod
    def _format_optional(value: Optional[Any]) -> str:
        """Serialize an optional registry value."""
        if value is None:
            return "-"
        if isinstance(value, bool):
            return "True" if value else "False"
        return str(value)

    def _format_row(self, entry: SpeciesEntry) -> Dict[str, str]:
        """Convert one species entry into a TSV row."""
        entry_dict = asdict(entry)
        return {
            key: self._format_optional(entry_dict[key])
            for key in _FIELD_NAMES
        }
