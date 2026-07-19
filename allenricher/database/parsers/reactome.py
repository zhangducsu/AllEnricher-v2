"""Build Reactome pathway gene sets with official names and hierarchy paths."""

import gzip
from functools import lru_cache
from pathlib import Path
from typing import Dict, Optional, Set


class ReactomeParser:
    """Build Reactome pathway files and hierarchy metadata."""

    @staticmethod
    def _open_gz_or_text(filepath: str):
        """Open plain text or gzip-compressed input transparently."""
        if filepath.endswith('.gz'):
            return gzip.open(filepath, 'rt', encoding='utf-8')
        else:
            return open(filepath, 'r', encoding='utf-8')

    @staticmethod
    def parse_ncbi2reactome(
        ncbi2reactome_path: str,
        gene_info_path: str,
        taxid: int,
        species: str,
        outdir: str,
        pathways_path: Optional[str] = None,
        relations_path: Optional[str] = None,
    ) -> None:
        """Build Reactome gene sets from the NCBI-to-Reactome mapping."""
        outdir_path = Path(outdir)
        outdir_path.mkdir(parents=True, exist_ok=True)

        print(f"|---ReactomeParser: Start deciphering NCBI2Reactome"
              f"(taxid={taxid}, species={species})")

        # Step 1: Read gene_info.gz, create geneid-> symbol map
        # Preserve the v1 mapping rule: filter gene_info by TaxID and extract GeneID-symbol pairs.
        gene_id_to_symbol: Dict[str, str] = {}
        print(f"|---Read file: {gene_info_path}")

        with ReactomeParser._open_gz_or_text(gene_info_path) as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith('#'):
                    continue
                parts = line.split('\t')
                file_taxid = parts[0]
                gene_id = parts[1]
                symbol = parts[2]
                if int(file_taxid) == taxid:
                    gene_id_to_symbol[gene_id] = symbol

        print(f"|---Found it.{len(gene_id_to_symbol)}Genome (taxid={taxid})")

        # Step 2: Read NCBI2Reactome files, filter the species specified
        # v1 Logical:
        #   Column [0] =geneid, [1] =pathway_id, [3] =pathway_name
        #   Pathway_id format "R-HSA-12345", with bb [1] (upsize) compared to uc (species)
        species_upper = species.upper()
        gene2pathway_file = outdir_path / f"{species}.gene2pathway.txt"
        tab_file = outdir_path / f"{species}.Reactome2gene.tab.gz"

        all_pathways: Set[str] = set()
        all_symbols: Set[str] = set()
        tab: Dict[str, Dict[str, int]] = {}  # {symbol: {pathway_id: 1}}
        pathway_names: Dict[str, str] = {}   # {pathway_id: pathway_name}

        print(f"|---Read file: {ncbi2reactome_path}")
        n = 0

        with ReactomeParser._open_gz_or_text(ncbi2reactome_path) as f_in, \
             open(gene2pathway_file, 'w', encoding='utf-8') as f_txt:
            for line in f_in:
                line = line.strip()
                if not line or line.startswith('#'):
                    continue
                parts = line.split('\t')
                gene_id = parts[0]
                pathway_id = parts[1]
                pathway_name = parts[3]

                # v1 Logical: Pathway_id format "R-HSA-12345", compared with bb[1] and uc(name)
                pathway_parts = pathway_id.split('-')
                if len(pathway_parts) < 2:
                    continue
                pathway_species_code = pathway_parts[1]
                if pathway_species_code != species_upper:
                    continue

                # Get agene symbol or use agene_id if not
                symbol = gene_id_to_symbol.get(gene_id, gene_id)

                all_symbols.add(symbol)
                all_pathways.add(pathway_id)
                pathway_names[pathway_id] = pathway_name

                # Write Gene2pathway.txt
                f_txt.write(f"{symbol}\t{gene_id}\t{pathway_id}\t{pathway_name}\n")

                # Build tab matrix data
                if symbol not in tab:
                    tab[symbol] = {}
                tab[symbol][pathway_id] = 1
                n += 1

        if n == 0:
            raise ValueError(
                f"[Error] No species= found in NCBI2Reactome file{species}The pass notes!"
            )

        print(f"|--- Parsed {n} Reactome pathway-gene associations.")

        # Step 3: Write to Reactome2gene.tab.gz
        # v1 Logic: Header Gene\\tpathway_id1\\tpathway_id2\\t...
        sorted_pathways = sorted(all_pathways)
        print(f"|---Writing file: {tab_file}")

        with gzip.open(tab_file, 'wt', encoding='utf-8') as f:
            header = ["Gene"] + sorted_pathways
            f.write('\t'.join(header) + '\n')

            for symbol in sorted(all_symbols):
                row = [symbol]
                for pid in sorted_pathways:
                    val = tab.get(symbol, {}).get(pid, 0)
                    row.append(str(val))
                f.write('\t'.join(row) + '\n')

        hierarchy_paths: Dict[str, str] = {}
        if pathways_path and relations_path:
            pathway_names, hierarchy_paths = ReactomeParser._load_hierarchies(
                pathways_path,
                relations_path,
                all_pathways,
                pathway_names,
            )

        # Step 4: Write Reactome2disc.gz
        disc_file = outdir_path / f"{species}.Reactome2disc.gz"
        print(f"|---Writing file: {disc_file}")

        with gzip.open(disc_file, 'wt', encoding='utf-8') as f:
            for pid in sorted_pathways:
                pname = pathway_names.get(pid, pid)
                hierarchy = hierarchy_paths.get(pid, "")
                f.write(f"{pid}\t{pname}\t{hierarchy}\n")

        print(f"|---ReactomeParser: Reactome database built")

    @staticmethod
    def _load_hierarchies(
        pathways_path: str,
        relations_path: str,
        selected_pathways: Set[str],
        fallback_names: Dict[str, str],
    ) -> tuple[Dict[str, str], Dict[str, str]]:
        """Return official pathway names and deterministic root-to-term paths."""
        names = dict(fallback_names)
        with ReactomeParser._open_gz_or_text(pathways_path) as handle:
            for line in handle:
                parts = line.rstrip("\r\n").split("\t")
                if len(parts) >= 2:
                    names[parts[0]] = parts[1]

        parents: Dict[str, Set[str]] = {}
        with ReactomeParser._open_gz_or_text(relations_path) as handle:
            for line in handle:
                parts = line.rstrip("\r\n").split("\t")
                if len(parts) >= 2:
                    parents.setdefault(parts[1], set()).add(parts[0])

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
        for term_id in selected_pathways:
            path = best_path(term_id)
            labels = [names.get(path_id, path_id) for path_id in path]
            if len(labels) > 1:
                hierarchies[term_id] = "|".join(labels)
        return names, hierarchies
