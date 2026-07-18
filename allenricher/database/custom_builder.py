"""Build a local enrichment database from user-provided annotations or GMT files."""

import gzip
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional

import pandas as pd

from .parsers.annotation_parser import AnnotationParser


class CustomDatabaseBuilder:
    """Build a custom database while preserving term names and hierarchy metadata."""

    def __init__(self, root_dir: str = "./database"):
        """Initialize a custom database builder.

        Args:
            root_dir: Root directory for versioned database builds.
        """
        self.root_dir = Path(root_dir)

    def build_from_annotation(
        self,
        annotation_file: str,
        species: str,
        taxid: int,
        db_name: str,
        format_type: Optional[str] = None,
        hierarchy_separator: str = '|'
    ) -> str:
        """Build a custom database from a gene-to-term annotation table."""
        annotation_path = Path(annotation_file)
        if not annotation_path.exists():
            raise FileNotFoundError(f"Annotation file does not exist: {annotation_file}")

        if annotation_path.name.lower().endswith((".gmt", ".gmt.gz")):
            return self.build_from_gmt(
                gmt_file=str(annotation_path),
                species=species,
                taxid=taxid,
                db_name=db_name,
            )

        # Parse and normalize the user annotation table.
        try:
            parser = AnnotationParser(
                filepath=str(annotation_path),
                format_type=format_type,
                hierarchy_separator=hierarchy_separator
            )
            parser.parse()
            term_genes = parser.get_term_genes()
            term_names = parser.get_term_names()
            term_hierarchies = parser.get_term_hierarchies()
        except (ValueError, FileNotFoundError):
            raise ValueError(
                f"No valid gene-to-term associations were found in: {annotation_file}"
            )

        if not term_genes:
            raise ValueError(f"No valid gene-to-term associations were found in: {annotation_file}")

        # Create a versioned species output directory.
        date_str = datetime.now().strftime("%Y%m%d")
        outdir = self.root_dir / "organism" / f"v{date_str}" / species
        outdir.mkdir(parents=True, exist_ok=True)

        print(f"\n{'='*60}")
        print(f"Building custom database: {species}.{db_name}")
        print(f"Annotation file: {annotation_file}")
        print(f"Output directory: {outdir}")
        print(f"Terms: {len(term_genes)}")
        print(f"{'='*60}")

        # Generate the gene-by-term membership matrix.
        print("|--- Step 1/3: writing gene-by-term membership matrix...")
        self._create_gene_matrix(term_genes, species, db_name, outdir)

        # Generate term descriptions and hierarchy metadata.
        print("|--- Step 2/3: writing term descriptions...")
        self._create_description_file(
            term_names, term_hierarchies, db_name, outdir
        )

        # Generate a portable GMT representation.
        print("|--- Step 3/3: writing GMT file...")
        self._create_gmt_file(term_genes, term_names, species, db_name, outdir)

        # Verify every required database artifact.
        expected_files = [
            f"{species}.{db_name}2gene.tab.gz",
            f"{db_name}2disc.gz",
            f"{species}.{db_name}.gmt.gz",
        ]
        for fname in expected_files:
            fpath = outdir / fname
            if fpath.exists():
                print(f"    [OK] {fname}")
            else:
                print(f"    [FAIL] {fname} was not generated")

        print(f"\nCustom database build completed -> {outdir}")
        return str(outdir)

    def build_from_gmt(
        self,
        gmt_file: str,
        species: str,
        taxid: int,
        db_name: str,
    ) -> str:
        """Build a custom database directly from a GMT file without a dense matrix."""
        gmt_path = Path(gmt_file)
        if not gmt_path.exists():
            raise FileNotFoundError(f"GMT file does not exist: {gmt_file}")

        opener = gzip.open if gmt_path.name.lower().endswith(".gz") else open
        term_genes: Dict[str, set[str]] = {}
        term_names: Dict[str, str] = {}
        with opener(gmt_path, "rt", encoding="utf-8") as handle:
            for line_number, line in enumerate(handle, start=1):
                parts = [part.strip() for part in line.rstrip("\r\n").split("\t")]
                if not parts or not parts[0] or parts[0].startswith("#"):
                    continue
                if len(parts) < 3:
                    raise ValueError(
                        f"GMT line {line_number} has fewer than three columns: {gmt_file}"
                    )
                term_id, term_name = parts[:2]
                genes = {gene for gene in parts[2:] if gene}
                if not genes:
                    continue
                term_genes.setdefault(term_id, set()).update(genes)
                term_names.setdefault(term_id, term_name or term_id)

        if not term_genes:
            raise ValueError(f"No valid gene sets were found in the GMT file: {gmt_file}")

        date_str = datetime.now().strftime("%Y%m%d")
        outdir = self.root_dir / "organism" / f"v{date_str}" / species
        outdir.mkdir(parents=True, exist_ok=True)
        self._create_description_file(term_names, {}, db_name, outdir)
        self._create_gmt_file(term_genes, term_names, species, db_name, outdir)
        return str(outdir)

    def _create_gene_matrix(
        self,
        term_genes: Dict[str, List[str]],
        species: str,
        db_name: str,
        outdir: Path
    ) -> str:
        """Write a binary gene-by-term membership matrix."""
        # Stable ordering makes repeated builds byte-for-byte reproducible.
        all_genes = sorted({gene for genes in term_genes.values() for gene in genes})

        # Collect all entries (sort)
        terms = sorted(term_genes.keys())

        # Write a binary membership matrix.
        data = {"Gene": all_genes}
        for term in terms:
            gene_set = set(term_genes[term])
            data[term] = [1 if g in gene_set else 0 for g in all_genes]

        df = pd.DataFrame(data)

        # Save as gzip Compressed TSV
        output_path = outdir / f"{species}.{db_name}2gene.tab.gz"
        df.to_csv(output_path, sep='\t', index=False, compression='gzip')

        print(f"    Membership matrix: {len(all_genes)} genes x {len(terms)} terms -> {output_path.name}")
        return str(output_path)

    def _create_description_file(
        self,
        term_names: Dict[str, str],
        term_hierarchies: Dict[str, str],
        db_name: str,
        outdir: Path
    ) -> str:
        """Write stable term identifiers, names, and hierarchy metadata."""
        output_path = outdir / f"{db_name}2disc.gz"

        with gzip.open(output_path, 'wt', encoding='utf-8') as f:
            for term_id in sorted(term_names.keys()):
                term_name = term_names[term_id]
                hierarchy = term_hierarchies.get(term_id, term_name)
                f.write(f"{term_id}\t{term_name}\t{hierarchy}\n")

        print(f"    Description table: {len(term_names)} terms -> {output_path.name}")
        return str(output_path)

    def _create_gmt_file(
        self,
        term_genes: Dict[str, List[str]],
        term_names: Dict[str, str],
        species: str,
        db_name: str,
        outdir: Path
    ) -> str:
        """Write a compressed GMT gene-set file."""
        output_path = outdir / f"{species}.{db_name}.gmt.gz"

        count = 0
        with gzip.open(output_path, 'wt', encoding='utf-8') as f:
            for term_id in sorted(term_genes.keys()):
                genes = sorted(set(term_genes[term_id]))
                if not genes:
                    continue
                term_name = term_names.get(term_id, "")
                line = f"{term_id}\t{term_name}\t" + "\t".join(genes) + "\n"
                f.write(line)
                count += 1

        print(f"    GMT file: {count} gene sets -> {output_path.name}")
        return str(output_path)
