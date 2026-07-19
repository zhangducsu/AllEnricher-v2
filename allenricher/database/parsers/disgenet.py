"""Convert DisGeNET gene-disease associations into enrichment gene sets."""

import gzip
import re
from pathlib import Path
from typing import Dict, List, Optional, Set


class DisGeNETParser:
    """Build disease gene sets from DisGeNET association records."""

    @staticmethod
    def _open_gz_or_text(filepath: str):
        """Open plain text or gzip-compressed input transparently."""
        if filepath.endswith('.gz'):
            return gzip.open(filepath, 'rt', encoding='utf-8')
        else:
            return open(filepath, 'r', encoding='utf-8')

    @staticmethod
    def parse_associations(assoc_path: str, gene_info_path: str,
                           taxid: int, outdir: str) -> None:
        """Parse DisGeNET associations into disease gene sets."""
        outdir_path = Path(outdir)
        outdir_path.mkdir(parents=True, exist_ok=True)

        print(f"|---DisGeNETParser: Start parsing DisGeNET files (taxid={taxid})")

        # Step 1: Read the gene_info.gz to create an effective gene pool
        # v1 Logic: DisGeNET_gene_filter.pl Read gene_info, create symbol map
        valid_genes: Set[str] = set()
        print(f"|---Read file: {gene_info_path}")

        with DisGeNETParser._open_gz_or_text(gene_info_path) as f:
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

        # Step 2: Read DisGeNET association file, extract the gene-cui association
        # v1 Logical:
        # Columns: gene symbol, gene ID, disease CUI, and disease name.
        #   Filter CUI: The beginning line (CUI format: Cxxxxx)
        #   Replace space and hyphenation with underlined
        #   Remove single quotes
        all_cuis: Set[str] = set()
        cui_data: Dict[str, Dict[str, str]] = {}  # {symbol: {cui: disease_name}}
        all_symbols: Set[str] = set()
        n = 0

        print(f"|---Read file: {assoc_path}")

        with DisGeNETParser._open_gz_or_text(assoc_path) as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith('#'):
                    continue

                # Remove single quote sign (v1) Logical: sed "s/"//g")
                line = line.replace("'", "")

                parts = line.split('\t')
                if len(parts) < 6:
                    continue

                gene_symbol = parts[0]
                gene_id = parts[1]
                disease_id = parts[4]
                disease_name = parts[5]

                # v1 Logic: Filter CUI: The beginning line (CUI format: Cxxxxx)
                if not disease_id.startswith('CUI:'):
                    continue

                # Standarded gene symbol (capsulation comparison)
                symbol_upper = gene_symbol.upper()

                # v1 Logic: filter valid genes with gene_info
                if symbol_upper not in valid_genes:
                    continue

                # v1 Logical: deposition_name spaces and hyphenation to underline
                disease_name = disease_name.replace(' ', '_').replace('-', '_')

                all_cuis.add(disease_id)
                all_symbols.add(gene_symbol)

                if gene_symbol not in cui_data:
                    cui_data[gene_symbol] = {}
                cui_data[gene_symbol][disease_id] = disease_name
                n += 1

        if n == 0:
            raise ValueError(
                "No valid disease-gene associations were found in the DisGeNET file"
            )

        print(f"|--- Parsed {n} DisGeNET disease-gene associations.")

        # Step 3: written into CUI2gene.tab.gz
        # v1 Logic: Header Gene\\tCUI1\\tCUI2\\t...
        sorted_cuis = sorted(all_cuis)
        tab_file = outdir_path / "hsa.CUI2gene.tab.gz"
        print(f"|---Writing file: {tab_file}")

        # Retain disease CUIs that have at least one gene association.
        uniq_disc: Dict[str, str] = {}  # {cui: disease_name}

        with gzip.open(tab_file, 'wt', encoding='utf-8') as f:
            header = ["Gene"] + sorted_cuis
            f.write('\t'.join(header) + '\n')

            for symbol in sorted(all_symbols):
                row = [symbol]
                for cui in sorted_cuis:
                    if cui in cui_data.get(symbol, {}):
                        row.append('1')
                        uniq_disc[cui] = cui_data[symbol][cui]
                    else:
                        row.append('0')
                f.write('\t'.join(row) + '\n')

        # Step 4: write CUI2disc.gz
        # Preserve the v1 format: CUI, disease name, limited to associated diseases.
        disc_file = outdir_path / "hsa.CUI2disc.gz"
        print(f"|---Writing file: {disc_file}")

        with gzip.open(disc_file, 'wt', encoding='utf-8') as f:
            for cui in sorted(uniq_disc.keys()):
                f.write(f"{cui}\t{uniq_disc[cui]}\n")

        print(f"|---Total{len(uniq_disc)}CUI Term")
        print(f"|---DisGeNETParser: DisGeNET database built")
