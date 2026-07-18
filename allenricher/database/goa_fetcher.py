"""Retrieve and convert species-specific UniProt Gene Ontology annotations."""

from __future__ import annotations

import gzip
import logging
from pathlib import Path
from typing import Dict, Optional, Set, Tuple

import requests

logger = logging.getLogger(__name__)


class GOAFetcher:
    """Retrieve and transform one species annotation file from UniProt GOA."""

    BASE_URL = "https://ftp.ebi.ac.uk/pub/databases/GO/goa/proteomes"
    UA = "Mozilla/5.0 (Windows NT 10.0; Win64; x64)"
    TIMEOUT = 30

    def __init__(self, cache_dir: str, overwrite: bool = False):
        """
        Args:
cache_dir: GOA File Cache Directory
overwrite: Whether to overwrite an existing cache
        """
        self.cache_dir = Path(cache_dir)
        self.cache_dir.mkdir(parents=True, exist_ok=True)
        self.overwrite = overwrite

    def fetch_species_data(
        self,
        taxid: int,
        latin_name: str,
        goa_filename: Optional[str] = None,
    ) -> Path:
        """Retrieve and convert source annotations for one species."""
        # EBI GOA filename format{taxid}.goa (like 9606.goa)
        if goa_filename is None:
            goa_filename = f"{taxid}.goa"

        # Local filenames match URL filenames
        local_file = self.cache_dir / f"{goa_filename}.gz"

        if local_file.exists() and not self.overwrite:
            logger.info("GOA file cached, skipping download: %s", local_file)
            return local_file

        url = f"{self.BASE_URL}/{goa_filename}"
        logger.info("Download GOA Documentation: %s -> %s", url, local_file)

        resp = requests.get(
            url,
            headers={"User-Agent": self.UA},
            timeout=self.TIMEOUT,
        )
        resp.raise_for_status()

        local_file.write_bytes(resp.content)
        logger.info("GoA file download completed: %s (%d bytes)", local_file, len(resp.content))

        return local_file

    def parse_goa_file(
        self,
        goa_file: Path,
        taxid: int,
    ) -> Tuple[Dict[str, Set[str]], Set[str]]:
        """Parse a Gene Association Format file into normalized records."""
        if not goa_file.exists():
            raise FileNotFoundError(f"GOA file does not exist: {goa_file}")

        opener: callable = gzip.open if str(goa_file).endswith(".gz") else open

        gene_to_go: Dict[str, Set[str]] = {}
        all_genes: Set[str] = set()
        skipped_not = 0
        skipped_wrong_taxid = 0
        total_lines = 0

        with opener(goa_file, "rt", encoding="utf-8") as f:
            for line in f:
                line = line.rstrip("\n\r")
                # Skip the comment line
                if line.startswith("!"):
                    continue

                parts = line.split("\t")
                if len(parts) < 17:
                    continue

                total_lines += 1

                # GAF 2.2 Column 13 (subscript 12 under zero base) is Taxon.
                taxon_field = parts[12]
                expected_taxon = f"taxon:{taxid}"
                if taxon_field.split('|')[0] != expected_taxon:
                    skipped_wrong_taxid += 1
                    continue

                relation = parts[3]
                if 'NOT' in relation.split('|'):
                    skipped_not += 1
                    continue

                # Column 2: DB_Object_ID, Column 3: DB_Object_Symbol
                db_object_id = parts[1]
                symbol = parts[2]

                symbol = symbol or db_object_id

                # Column 5: GO_ID
                go_id = parts[4]
                if not go_id.startswith("GO:"):
                    continue

                gene_to_go.setdefault(symbol, set()).add(go_id)
                all_genes.add(symbol)

        logger.info(
            "GOA parsed: TaxID=%d, total rows=%d, mapped genes=%d,"
            "Skip (Error taxid) =%d, Skip (NOT) =%d",
            taxid, total_lines, len(all_genes),
            skipped_wrong_taxid, skipped_not,
        )

        return gene_to_go, all_genes

    @staticmethod
    def build_go2gene_matrix(
        gene_to_go: Dict[str, Set[str]],
        all_genes: Set[str],
        all_go_terms: Set[str],
        output_path: Path,
    ) -> None:
        """Write a binary GO-term-by-gene membership matrix."""
        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)

        sorted_go = sorted(all_go_terms)
        sorted_genes = sorted(all_genes)

        with gzip.open(output_path, "wt", encoding="utf-8") as f:
            # Write
            header = ["Gene"] + sorted_go
            f.write("\t".join(header) + "\n")

            # Writing Data Line
            for symbol in sorted_genes:
                go_set = gene_to_go.get(symbol, set())
                row = [symbol] + ["1" if go_id in go_set else "0" for go_id in sorted_go]
                f.write("\t".join(row) + "\n")

        logger.info(
            "GO2gene.tab.gz Generate: %s (%d genes x %d GO terms)",
            output_path, len(sorted_genes), len(sorted_go),
        )

    @staticmethod
    def build_gene2go_list(
        gene_to_go: Dict[str, Set[str]],
        go_names: Dict[str, str],
        output_path: Path,
    ) -> None:
        """Write the normalized gene-to-GO association table."""
        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)

        with open(output_path, "w", encoding="utf-8") as f:
            for symbol in sorted(gene_to_go):
                for go_id in sorted(gene_to_go[symbol]):
                    go_name = go_names.get(go_id, "")
                    f.write(f"{symbol}\t{go_id}\t\t{go_name}\n")

        total_annotations = sum(len(v) for v in gene_to_go.values())
        logger.info(
            "Gene2go.txt Generated: %s (%d genes, %d annotations)",
            output_path, len(gene_to_go), total_annotations,
        )
