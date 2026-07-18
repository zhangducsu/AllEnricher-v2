"""Extract pathway genes from WikiPathways GPML archives."""

import gzip
import xml.etree.ElementTree as ET
import zipfile
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple


class WikiPathwaysGPMLParser:
    """Build WikiPathways gene sets from GPML XML archives."""

    GPML_NS = "http://pathvisio.org/GPML/2013a"

    def parse_gpml_zip(
        self,
        gpml_zip_path: str,
        gene_info_path: Optional[str] = None,
        taxid: Optional[str] = None
    ) -> Tuple[Dict[str, Set[str]], Dict[str, str]]:
        """Parse all GPML pathway files in one archive."""
        zip_path = Path(gpml_zip_path)
        if not zip_path.exists():
            raise FileNotFoundError(f"GPML ZIP file does not exist: {gpml_zip_path}")

        # Load the gene-ID-to-symbol map when gene_info is available.
        id_mapping: Dict[str, str] = {}
        if gene_info_path and taxid:
            id_mapping = self._load_gene_id_mapping(gene_info_path, taxid)

        gene_sets: Dict[str, Set[str]] = {}
        descriptions: Dict[str, str] = {}

        with zipfile.ZipFile(gpml_zip_path, 'r') as zf:
            for filename in zf.namelist():
                if not filename.endswith('.gpml'):
                    continue

                # Draw from file name path_id (e. g. WP1234.gpml->WP1234)
                pathway_id = Path(filename).stem

                try:
                    with zf.open(filename) as f:
                        tree = ET.parse(f)
                        root = tree.getroot()

                        # Get Access Name
                        pathway_name = root.get('Name', pathway_id)
                        descriptions[pathway_id] = pathway_name

                        # Gene extraction
                        genes = self._extract_genes_from_gpml(root, id_mapping)
                        gene_sets[pathway_id] = genes

                except ET.ParseError as e:
                    print(f"Warning: failed to parse {filename}: {e}")
                    continue

        return gene_sets, descriptions

    def _extract_genes_from_gpml(
        self,
        root: ET.Element,
        id_mapping: Dict[str, str]
    ) -> Set[str]:
        """Extract gene identifiers from one GPML XML document."""
        genes: Set[str] = set()
        ns = {'gpml': self.GPML_NS}

        # Find all DataNode elements
        for datanode in root.findall('.//gpml:DataNode', ns):
            # Get TextLabel
            text_label = datanode.get('TextLabel', '').strip()

            # Find Xref sub-elements
            xref = datanode.find('gpml:Xref', ns)
            if xref is not None:
                database = xref.get('Database', '')
                gene_id = xref.get('ID', '')

                # If it is Entrez Gene and provides a map
                if database == 'Entrez Gene' and gene_id and gene_id in id_mapping:
                    genes.add(id_mapping[gene_id])
                elif text_label:
                    # Use TextLabel
                    genes.add(text_label)
            elif text_label:
                # TextLabel when Xref is not available
                genes.add(text_label)

        return genes

    def _load_gene_id_mapping(
        self,
        gene_info_path: str,
        taxid: str
    ) -> Dict[str, str]:
        """Load the NCBI Gene ID-to-symbol mapping."""
        gene_info_file = Path(gene_info_path)
        if not gene_info_file.exists():
            raise FileNotFoundError(f"The gene_info file does not exist: {gene_info_path}")

        mapping: Dict[str, str] = {}

        open_func = gzip.open if gene_info_path.endswith('.gz') else open
        mode = 'rt' if gene_info_path.endswith('.gz') else 'r'

        with open_func(gene_info_path, mode, encoding='utf-8') as f:
            # Skip title line
            header = f.readline()

            for line in f:
                line = line.strip()
                if not line:
                    continue

                parts = line.split('\t')
                if len(parts) < 3:
                    continue

                # Gene_info format: tax_id Geneid Symbol...
                file_taxid = parts[0]
                gene_id = parts[1]
                symbol = parts[2]

                # Load only specified species
                if file_taxid == taxid:
                    mapping[gene_id] = symbol

        return mapping

    def build_database_from_gpml(
        self,
        gpml_zip_path: str,
        output_dir: str,
        species: str,
        taxid: Optional[str] = None,
        gene_info_path: Optional[str] = None
    ) -> Tuple[str, str]:
        """Build WikiPathways artifacts from a GPML archive."""
        outdir_path = Path(output_dir)
        outdir_path.mkdir(parents=True, exist_ok=True)

        print(f"--- WikiPathways GPMMLparser: Start building WikiPathways database (SPpecies={species})")

        # Step 1: parsing GPML ZIP files
        print(f"|---Read GPML ZIP file: {gpml_zip_path}")
        gene_sets, descriptions = self.parse_gpml_zip(
            gpml_zip_path,
            gene_info_path=gene_info_path,
            taxid=taxid
        )

        if not gene_sets:
            raise ValueError("No valid pathways were found in the WikiPathways GPML archive")

        print(f"|--- Parsed {len(gene_sets)} WikiPathways gene sets.")

        # Step 2: Build the pathway-by-gene membership matrix.
        all_pathways: Set[str] = set(gene_sets.keys())
        all_genes: Set[str] = set()
        tab: Dict[str, Dict[str, int]] = {}  # {gene: {pathway_id: 1}}

        for pathway_id, genes in gene_sets.items():
            for gene in genes:
                all_genes.add(gene)
                if gene not in tab:
                    tab[gene] = {}
                tab[gene][pathway_id] = 1

        if not all_genes:
            raise ValueError("No valid gene-pathway associations were found")

        print(f"|---We found them together.{len(all_genes)}A gene.{sum(len(v) for v in tab.values())}Rhodes-access connections")

        # Step 3: write WikiPathways2gene.tab.gz
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

        # Step 4: Write WikiPathways2disc.gz
        disc_file = outdir_path / f"{species}.WikiPathways2disc.gz"
        print(f"|---Writing file: {disc_file}")

        with gzip.open(disc_file, 'wt', encoding='utf-8') as f:
            for pid in sorted_pathways:
                pname = descriptions.get(pid, pid)
                # Replace spaces with underlined spaces (same as with other parser)
                pname_underscore = pname.replace(' ', '_')
                f.write(f"{pid}\t{pname_underscore}\n")

        print(f"--- WikiPathways GPPMLparser: WikiPathways database constructed")

        return str(tab_file), str(disc_file)
