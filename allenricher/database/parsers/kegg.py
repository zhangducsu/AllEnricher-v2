"""Convert downloaded KEGG pathway annotations into database artifacts."""

import gzip
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple


class KEGGParser:
    """Write analysis-ready KEGG pathway files."""

    @staticmethod
    def build_database(species: str, gene_info_path: str,
                       gene2pathway_path: str, outdir: str,
                       pathway_summary_path: Optional[str] = None,
                       taxid: Optional[int] = None) -> None:
        """Build normalized database artifacts from parsed source data."""
        outdir_path = Path(outdir)
        outdir_path.mkdir(parents=True, exist_ok=True)

        print(f"|---KEGG parser: building the database for species {species}")

        # Step 1: read gene_info.gz and retain the requested TaxID.
        valid_genes: Set[str] = set()
        print(f"|---Read file: {gene_info_path}")

        if gene_info_path.endswith('.gz'):
            f_open = gzip.open(gene_info_path, 'rt', encoding='utf-8')
        else:
            f_open = open(gene_info_path, 'r', encoding='utf-8')

        with f_open:
            for line in f_open:
                line = line.strip()
                if not line or line.startswith('#'):
                    continue
                parts = line.split('\t')
                if len(parts) >= 3 and (taxid is None or parts[0] == str(taxid)):
                    valid_genes.add(parts[2])  # symbol column

        print(f"|---Valid genes for the requested TaxID: {len(valid_genes)}")

        # Step 2: read optional pathway hierarchy metadata.
        # Format: Category\\tSubcategory\\tpathway_id\\tpathway_name\\turl
        pathway_categories: Dict[str, str] = {}  # {pathway_id: "Category|Subcategory|PathwayName"}

        if pathway_summary_path and Path(pathway_summary_path).exists():
            print(f"|---Read file: {pathway_summary_path}")
            with open(pathway_summary_path, 'r', encoding='utf-8') as f:
                for line in f:
                    line = line.strip()
                    if not line or line.startswith('#'):
                        continue
                    parts = line.split('\t')
                    if len(parts) >= 4:
                        category = parts[0].replace(' ', '_')
                        subcategory = parts[1].replace(' ', '_')
                        pathway_id = parts[2]
                        pathway_name = parts[3].replace(' ', '_')
                        # Add species prefix
                        if not pathway_id.startswith(species):
                            pathway_id = f"{species}{pathway_id}"
                        pathway_categories[pathway_id] = (
                            f"{category}|{subcategory}|{pathway_name}"
                        )

        # Step 3: build the gene-by-pathway membership matrix.
        # Input format: gene_symbol, Entrez ID, pathway ID, pathway name.
        all_pathways: Set[str] = set()
        all_genes: Set[str] = set()
        tab: Dict[str, Dict[str, int]] = {}  # {gene: {pathway_id: 1}}
        pathway_names: Dict[str, str] = {}   # {pathway_id: pathway_name}

        print(f"|---Read file: {gene2pathway_path}")
        n = 0

        with open(gene2pathway_path, 'r', encoding='utf-8') as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith('#'):
                    continue
                parts = line.split('\t')
                if len(parts) < 4:
                    continue

                gene_symbol = parts[0]
                pathway_id = parts[2]
                pathway_name = parts[3]

                # Add species prefix (e.g. 04110-> hsa04110), in line with v1 database format
                if not pathway_id.startswith(species):
                    pathway_id = f"{species}{pathway_id}"

                # Preserve source symbols for backward compatibility. gene_info
                # is used to validate the requested species above.
                symbol = gene_symbol

                all_genes.add(symbol)
                all_pathways.add(pathway_id)
                pathway_names[pathway_id] = pathway_name

                if symbol not in tab:
                    tab[symbol] = {}
                tab[symbol][pathway_id] = 1
                n += 1

        if n == 0:
            raise ValueError(
                "No valid pathway-gene associations were found in gene2pathway.txt"
            )

        print(f"|--- Parsed {n} KEGG pathway-gene associations.")

        # Step 4: write to kegg2gene.tab.gz
        # v1 Logic: Header Gene\\tpathway_id1\\tpathway_id2\\t...
        sorted_pathways = sorted(all_pathways)
        tab_file = outdir_path / f"{species}.kegg2gene.tab.gz"
        print(f"|---Writing file: {tab_file}")

        with gzip.open(tab_file, 'wt', encoding='utf-8') as f:
            header = ["Gene"] + sorted_pathways
            f.write('\t'.join(header) + '\n')

            for gene in sorted(all_genes):
                row = [gene]
                for pid in sorted_pathways:
                    val = tab.get(gene, {}).get(pid, 0)
                    row.append(str(val))
                f.write('\t'.join(row) + '\n')

        # Step 5: write pathway descriptions and optional hierarchy labels.
        disc_file = outdir_path / f"{species}.kegg2disc.gz"
        print(f"|---Writing file: {disc_file}")

        with gzip.open(disc_file, 'wt', encoding='utf-8') as f:
            for pid in sorted_pathways:
                if pid in pathway_categories:
                    disc_name = pathway_categories[pid]
                else:
                    # Fall back to the pathway name when hierarchy data are unavailable.
                    pname = pathway_names.get(pid, pid)
                    disc_name = pname.replace(' ', '_')
                f.write(f"{pid}\t{disc_name}\n")

        print("|---KEGG parser: database build complete")
