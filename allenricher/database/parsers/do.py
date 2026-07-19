"""Build Disease Ontology gene sets with names and hierarchy paths."""

import gzip
from functools import lru_cache
from pathlib import Path
from typing import Dict, List, Optional, Set


class DOParser:
    """Build Disease Ontology gene sets and deterministic hierarchy paths."""

    @staticmethod
    def _open_gz_or_text(filepath: str):
        """Open plain text or gzip-compressed input transparently."""
        if filepath.endswith('.gz'):
            return gzip.open(filepath, 'rt', encoding='utf-8')
        else:
            return open(filepath, 'r', encoding='utf-8')

    @staticmethod
    def parse_disease_files(
        disease_files: List[str],
        gene_info_path: str,
        taxid: int,
        outdir: str,
        ontology_path: Optional[str] = None,
    ) -> None:
        """Merge disease association files into Disease Ontology gene sets."""
        outdir_path = Path(outdir)
        outdir_path.mkdir(parents=True, exist_ok=True)

        print(f"--- Disease Ontology parser: reading source data for TaxID {taxid}")

        # Step 1: Read the gene_info.gz to create an effective gene pool
        # v1 Logical: gene_filter.pl Read gene_info, create symbol map
        valid_genes: Set[str] = set()
        print(f"|---Read file: {gene_info_path}")

        with DOParser._open_gz_or_text(gene_info_path) as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith('#'):
                    continue
                parts = line.split('\t')
                if len(parts) >= 3:
                    file_taxid = parts[0]
                    gene_id = parts[1]
                    symbol = parts[2]
                    if int(file_taxid) == taxid:
                        valid_genes.add(symbol.upper())

        print(f"|---Found it.{len(valid_genes)}A valid gene (taxid={taxid})")

        # Step 2: Read all files, extract the gene-DOID connections
        # v1 Logical:
        #   Column [1] = gene_symbol, [2] = DOID, [3] =disease_name
        #   Filter Doid: Line at the beginning
        #   Replace space and hyphenation with underlined
        #   Remove single quotes
        all_doids: Set[str] = set()
        do_data: Dict[str, Dict[str, str]] = {}  # {symbol: {doid: disease_name}}
        all_symbols: Set[str] = set()
        n = 0

        for disease_file in disease_files:
            print(f"|---Read file: {disease_file}")

            with DOParser._open_gz_or_text(disease_file) as f:
                for line in f:
                    line = line.strip()
                    if not line or line.startswith('#'):
                        continue

                    # Remove single quote sign (v1) Logical: sed "s/"//g")
                    line = line.replace("'", "")

                    parts = line.split('\t')
                    if len(parts) < 4:
                        continue

                    gene_symbol = parts[1]
                    doid = parts[2]
                    disease_name = parts[3]

                    # v1 Logic: Filter DoID: The line at the beginning
                    if not doid.startswith('DOID:'):
                        continue

                    # v1 Logical: deposition_name spaces and hyphenation to underline
                    disease_name = disease_name.replace(' ', '_').replace('-', '_')

                    # Standarded gene symbol (capsulation comparison)
                    symbol_upper = gene_symbol.upper()

                    # v1 Logic: filter valid genes with gene_info
                    if symbol_upper not in valid_genes:
                        continue

                    # Use original symbol (size)
                    all_doids.add(doid)
                    all_symbols.add(gene_symbol)

                    if gene_symbol not in do_data:
                        do_data[gene_symbol] = {}
                    do_data[gene_symbol][doid] = disease_name
                    n += 1

        if n == 0:
            raise ValueError(
                "No valid gene-Disease Ontology associations were found"
            )

        print(f"|--- Parsed {n} Disease Ontology term-gene associations.")

        ontology_names: Dict[str, str] = {}
        hierarchy_paths: Dict[str, str] = {}
        obsolete_terms: Set[str] = set()
        if ontology_path:
            ontology_names, hierarchy_paths, obsolete_terms = DOParser._load_ontology(
                ontology_path, all_doids
            )
            all_doids.difference_update(obsolete_terms)
            for symbol in list(do_data):
                do_data[symbol] = {
                    doid: name
                    for doid, name in do_data[symbol].items()
                    if doid not in obsolete_terms
                }

        # Step 3: Write Do 2gene.tab.gz
        # v1 Logic: Header Gene\\tDOID1\\tDOID2\\t...
        sorted_doids = sorted(all_doids)
        tab_file = outdir_path / "hsa.DO2gene.tab.gz"
        print(f"|---Writing file: {tab_file}")

        # Retain DOIDs that have at least one gene association.
        uniq_disc: Dict[str, str] = {}  # {doid: disease_name}

        with gzip.open(tab_file, 'wt', encoding='utf-8') as f:
            header = ["Gene"] + sorted_doids
            f.write('\t'.join(header) + '\n')

            for symbol in sorted(all_symbols):
                row = [symbol]
                for doid in sorted_doids:
                    if doid in do_data.get(symbol, {}):
                        row.append('1')
                        uniq_disc[doid] = ontology_names.get(doid, do_data[symbol][doid])
                    else:
                        row.append('0')
                f.write('\t'.join(row) + '\n')

        # Step 4: Write Do2disc.gz
        # Preserve the v1 format: DOID, disease name, limited to associated diseases.
        disc_file = outdir_path / "hsa.DO2disc.gz"
        print(f"|---Writing file: {disc_file}")

        with gzip.open(disc_file, 'wt', encoding='utf-8') as f:
            for doid in sorted(uniq_disc.keys()):
                hierarchy = hierarchy_paths.get(doid, "")
                f.write(f"{doid}\t{uniq_disc[doid]}\t{hierarchy}\n")

        print(f"|---Total{len(uniq_disc)}DO Term")
        print(f"|---DOOPSer: Do databases constructed and completed")

    @staticmethod
    def _load_ontology(
        ontology_path: str,
        selected_terms: Set[str],
    ) -> tuple[Dict[str, str], Dict[str, str], Set[str]]:
        """Read names and deterministic hierarchy paths from an OBO ontology."""
        names: Dict[str, str] = {}
        parents: Dict[str, Set[str]] = {}
        obsolete: Set[str] = set()
        current: Dict[str, object] = {}

        def flush() -> None:
            term_id = str(current.get("id") or "")
            if not term_id:
                return
            if current.get("obsolete"):
                obsolete.add(term_id)
                return
            name = str(current.get("name") or "")
            if name:
                names[term_id] = name
            term_parents = current.get("parents") or set()
            if term_parents:
                parents[term_id] = set(term_parents)

        with DOParser._open_gz_or_text(ontology_path) as handle:
            for raw_line in handle:
                line = raw_line.strip()
                if line == "[Term]":
                    flush()
                    current = {"parents": set()}
                elif not current or not line or line.startswith("!"):
                    continue
                elif line.startswith("id: "):
                    current["id"] = line[4:].strip()
                elif line.startswith("name: "):
                    current["name"] = line[6:].strip()
                elif line.startswith("is_a: "):
                    parent_id = line[6:].split(" ! ", 1)[0].strip()
                    current.setdefault("parents", set()).add(parent_id)
                elif line == "is_obsolete: true":
                    current["obsolete"] = True
        flush()

        @lru_cache(maxsize=None)
        def best_path(term_id: str) -> tuple[str, ...]:
            candidates = []
            for parent_id in sorted(parents.get(term_id, set())):
                if parent_id == term_id:
                    continue
                parent_path = best_path(parent_id)
                if term_id not in parent_path:
                    candidates.append((*parent_path, term_id))
            if not candidates:
                return (term_id,)
            return sorted(candidates, key=lambda path: (-len(path), path))[0]

        hierarchies = {}
        for term_id in selected_terms - obsolete:
            path = best_path(term_id)
            labels = [names.get(path_id, path_id) for path_id in path]
            if len(labels) > 1:
                hierarchies[term_id] = "|".join(labels)
        return names, hierarchies, obsolete & selected_terms
