"""Parse WikiPathways GMT records and normalize NCBI gene identifiers."""

import gzip
from pathlib import Path
from typing import Dict, Optional, Set, Tuple


class WikiPathwaysParser:
    """Build WikiPathways gene sets from GMT records."""

    @staticmethod
    def parse_gmt(gmt_path: str) -> Tuple[Dict[str, Set[str]], Dict[str, str]]:
        """Parse gene sets from GMT format."""
        gmt_file = Path(gmt_path)
        if not gmt_file.exists():
            raise FileNotFoundError(f"The GMT file does not exist: {gmt_path}")

        gene_sets: Dict[str, Set[str]] = {}
        descriptions: Dict[str, str] = {}

        open_func = gzip.open if gmt_path.endswith('.gz') else open
        mode = 'rt' if gmt_path.endswith('.gz') else 'r'

        with open_func(gmt_path, mode, encoding='utf-8') as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue

                parts = line.split('\t')
                if len(parts) < 3:
                    continue

                metadata = parts[0].split('%')
                wp_ids = [value for value in metadata if value.startswith('WP') and value[2:].isdigit()]
                if len(metadata) >= 3 and wp_ids:
                    pathway_id = wp_ids[0]
                    pathway_name = metadata[0]
                else:
                    pathway_id = parts[0]
                    pathway_name = parts[1]

                genes = {
                    gene
                    for field in parts[2:]
                    for gene in field.split('/')
                }

                # Filtering empty genes
                genes = {g.strip() for g in genes if g.strip()}

                gene_sets[pathway_id] = genes
                descriptions[pathway_id] = pathway_name

        return gene_sets, descriptions

    @staticmethod
    def load_gene_id_mapping(gene_info_path: str, taxid: int) -> Dict[str, str]:
        """Load NCBI Gene identifiers and symbols for one species."""
        gene_info_file = Path(gene_info_path)
        if not gene_info_file.exists():
            raise FileNotFoundError(f"The gene_info file does not exist: {gene_info_path}")

        id_mapping: Dict[str, str] = {}
        taxid_str = str(taxid)

        open_func = gzip.open if gene_info_path.endswith('.gz') else open
        mode = 'rt' if gene_info_path.endswith('.gz') else 'r'

        with open_func(gene_info_path, mode, encoding='utf-8') as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith('#'):
                    continue

                parts = line.split('\t')
                if len(parts) < 3:
                    continue

                # Check if the tabs_id match
                if parts[0] != taxid_str:
                    continue

                gene_id = parts[1]  # GeneID
                symbol = parts[2]   # Symbol

                if gene_id and symbol:
                    id_mapping[gene_id] = symbol

        return id_mapping

    @staticmethod
    def convert_ncbi_to_symbol(
        gene_sets: Dict[str, Set[str]],
        id_mapping: Dict[str, str]
    ) -> Dict[str, Set[str]]:
        """Convert namespaced NCBI Gene identifiers to symbols."""
        converted_sets: Dict[str, Set[str]] = {}

        for pathway_id, genes in gene_sets.items():
            converted_genes: Set[str] = set()

            for gene in genes:
                # Extracting the id from ncbigene: xx
                if gene.startswith('ncbigene:'):
                    ncbi_id = gene.split(':', 1)[1]
                elif gene.isdigit():
                    ncbi_id = gene
                else:
                    # Keep as if it were not for ncbigene: format
                    converted_genes.add(gene)
                    continue

                # Find the corresponding Symbol
                symbol = id_mapping.get(ncbi_id)
                if symbol:
                    converted_genes.add(symbol)
                else:
                    # Keep original ids for debug if map not found
                    converted_genes.add(gene)

            converted_sets[pathway_id] = converted_genes

        return converted_sets

    @staticmethod
    def build_database(
        gmt_path: str,
        output_dir: str,
        species: str,
        taxid: Optional[int] = None,
        gene_info_path: Optional[str] = None,
        valid_genes: Optional[Set[str]] = None
    ) -> None:
        """Build normalized database artifacts from parsed source data."""
        outdir_path = Path(output_dir)
        outdir_path.mkdir(parents=True, exist_ok=True)

        print(f"--- WikiPathwaysParser: Start building WikiPathways databases (SPECies={species})")

        # Step 1: Parsing GMT files
        print(f"|---Reading GMT files: {gmt_path}")
        gene_sets, descriptions = WikiPathwaysParser.parse_gmt(gmt_path)

        if not gene_sets:
            raise ValueError("No valid pathways were found in the WikiPathways GMT file")

        print(f"|--- Parsed {len(gene_sets)} WikiPathways gene sets.")

        # Step 2: Map NCBI Gene IDs to symbols when gene_info and TaxID are available.
        if gene_info_path and Path(gene_info_path).exists() and taxid is not None:
            print(f"|---Load NCBI Gene ID map (taxid={taxid})...")
            id_mapping = WikiPathwaysParser.load_gene_id_mapping(gene_info_path, taxid)
            print(f"|---Found it.{len(id_mapping)}Gene ID Map")

            print(f"|---Convert Gene ID to Symbol...")
            gene_sets = WikiPathwaysParser.convert_ncbi_to_symbol(gene_sets, id_mapping)
            print(f"|---ID conversion complete")

        # Step 3: Obtain effective gene pools (for filtering)
        valid_gene_set: Set[str] = set()
        if valid_genes:
            valid_gene_set = valid_genes
            print(f"|---Use of the available effective gene pool: {len(valid_gene_set)}Genome.")
        elif gene_info_path and Path(gene_info_path).exists():
            print(f"|---Load a list of active genes from the gene_info: {gene_info_path}")
            open_func = gzip.open if gene_info_path.endswith('.gz') else open
            mode = 'rt' if gene_info_path.endswith('.gz') else 'r'

            with open_func(gene_info_path, mode, encoding='utf-8') as f:
                for line in f:
                    line = line.strip()
                    if not line or line.startswith('#'):
                        continue
                    parts = line.split('\t')
                    if len(parts) >= 3 and (taxid is None or parts[0] == str(taxid)):
                        valid_gene_set.add(parts[2])  # symbol column

            print(f"|---Found it.{len(valid_gene_set)}A valid gene.")

        # Step 4: build the gene-by-pathway membership matrix.
        all_genes: Set[str] = set()
        tab: Dict[str, Dict[str, int]] = {}  # {gene: {pathway_id: 1}}

        for pathway_id, genes in gene_sets.items():
            for gene in genes:
                # Filter if effective gene pools are provided
                if valid_gene_set and gene not in valid_gene_set:
                    continue

                all_genes.add(gene)
                if gene not in tab:
                    tab[gene] = {}
                tab[gene][pathway_id] = 1

        if not all_genes:
            raise ValueError("No valid gene-pathway associations were found")

        all_pathways = {
            pathway_id
            for pathway_memberships in tab.values()
            for pathway_id in pathway_memberships
        }

        print(f"|---We found them together.{len(all_genes)}A gene.{sum(len(v) for v in tab.values())}Rhodes-access connections")

        # Step 5: Write WikiPathways2gene.tab.gz
        sorted_pathways = sorted(all_pathways)
        tab_file = outdir_path / f"{species}.WikiPathways2gene.tab.gz"
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

        # Step 6: Write WikiPathways2disc.gz
        disc_file = outdir_path / f"{species}.WikiPathways2disc.gz"
        print(f"|---Writing file: {disc_file}")

        with gzip.open(disc_file, 'wt', encoding='utf-8') as f:
            for pid in sorted_pathways:
                pname = descriptions.get(pid, pid)
                # Replace spaces with underlined (synchronous with Reactome, etc.)
                pname_underscore = pname.replace(' ', '_')
                f.write(f"{pid}\t{pname_underscore}\n")

        print(f"|---WikiPathwaysParser: WikiPathways database built")
