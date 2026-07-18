"""Generate compressed GMT gene-set files from built database artifacts."""

import gzip
from pathlib import Path
from typing import Dict, List, Optional, Tuple


class GMTGenerator:
    """Generate compressed GMT files from database matrices and descriptions."""

    def __init__(self, organism_dir: str):
        """Initialization GMT Generator

        Args:
organism_dir: Path to the species database catalogue, 
Usually database/organism/v{date}/{species}/
        """
        self.organism_dir = Path(organism_dir)

    # ============================
    # Internal tools methodology
    # ============================
    @staticmethod
    def _open_gz_or_text(filepath: str):
        """Open plain text or gzip-compressed input transparently."""
        if filepath.endswith('.gz'):
            return gzip.open(filepath, 'rt', encoding='utf-8')
        else:
            return open(filepath, 'r', encoding='utf-8')

    def _read_tab_matrix(self, tab_path: str) -> Tuple[List[str], Dict[str, List[str]]]:
        """Read a binary gene-by-term matrix into gene sets."""
        terms = []
        term_to_genes: Dict[str, List[str]] = {}

        with self._open_gz_or_text(tab_path) as f:
            header_line = f.readline().strip()
            if not header_line:
                return terms, term_to_genes

            parts = header_line.split('\t')
            terms = parts[1:]  # Skipping the first column "Gene"

            # Initialize
            for t in terms:
                term_to_genes[t] = []

            for line in f:
                line = line.strip()
                if not line:
                    continue
                cols = line.split('\t')
                gene = cols[0]
                for i, val in enumerate(cols[1:]):
                    if i < len(terms) and val == '1':
                        term_to_genes[terms[i]].append(gene)

        return terms, term_to_genes

    def _read_description(self, disc_path: str) -> Dict[str, str]:
        """Read term identifiers and display names."""
        descriptions: Dict[str, str] = {}

        with self._open_gz_or_text(disc_path) as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                # The description document may also contain parent entries/tier columns; the second column of GMT can only write descriptions, but the second column is not the same as the other
                # Otherwise, the following list will be misread as a gene.
                parts = line.split('\t')
                if len(parts) >= 2:
                    descriptions[parts[0]] = parts[1]
                elif len(parts) == 1 and parts[0]:
                    descriptions[parts[0]] = ""

        return descriptions

    def _write_gmt(self, term_to_genes: Dict[str, List[str]],
                   descriptions: Dict[str, str],
                   output_path: str) -> str:
        """Write gene sets in compressed GMT format."""
        out_path = Path(output_path)
        out_path.parent.mkdir(parents=True, exist_ok=True)

        count = 0
        with gzip.open(out_path, 'wt', encoding='utf-8') as f:
            for term_id, genes in sorted(term_to_genes.items()):
                if not genes:
                    continue
                desc = descriptions.get(term_id, "")
                line = f"{term_id}\t{desc}\t" + "\t".join(genes) + "\n"
                f.write(line)
                count += 1

        print(f"|--- GMTGenerator: Write {count} A gene set. -> {out_path}")
        return str(out_path)

    # ============================
    # GMT generation in databases
    # ============================
    def generate_go_gmt(self, species: str) -> str:
        """Generate GMT output for Gene Ontology."""
        tab_path = self.organism_dir / f"{species}.GO2gene.tab.gz"
        disc_path = self.organism_dir / "GO2disc.gz"

        if not tab_path.exists():
            raise FileNotFoundError(f"GO Genome File does not exist: {tab_path}")
        if not disc_path.exists():
            raise FileNotFoundError(f"GO description file does not exist: {disc_path}")

        print(f"|---GMTGenerator: Generate GO GMT (species=){species})")

        terms, term_to_genes = self._read_tab_matrix(str(tab_path))
        descriptions = self._read_description(str(disc_path))

        output_path = str(self.organism_dir / f"{species}.GO.gmt.gz")
        return self._write_gmt(term_to_genes, descriptions, output_path)

    def generate_kegg_gmt(self, species: str) -> str:
        """Generate GMT output for KEGG pathways."""
        tab_path = self.organism_dir / f"{species}.kegg2gene.tab.gz"
        disc_path = self.organism_dir / f"{species}.kegg2disc.gz"

        if not tab_path.exists():
            raise FileNotFoundError(f"KEGG gene matrix file does not exist: {tab_path}")
        if not disc_path.exists():
            raise FileNotFoundError(f"The KEGG description file does not exist: {disc_path}")

        print(f"|---GMTGenerator: Generate KEG GMT (species=){species})")

        terms, term_to_genes = self._read_tab_matrix(str(tab_path))
        descriptions = self._read_description(str(disc_path))

        output_path = str(self.organism_dir / f"{species}.KEGG.gmt.gz")
        return self._write_gmt(term_to_genes, descriptions, output_path)

    def generate_reactome_gmt(self, species: str) -> str:
        """Generate GMT output for Reactome pathways."""
        tab_path = self.organism_dir / f"{species}.Reactome2gene.tab.gz"
        disc_path = self.organism_dir / f"{species}.Reactome2disc.gz"

        if not tab_path.exists():
            raise FileNotFoundError(f"Reactome gene matrix file does not exist: {tab_path}")
        if not disc_path.exists():
            raise FileNotFoundError(f"The Reactome profile does not exist: {disc_path}")

        print(f"|---GMTGenerator: Generate{species})")

        terms, term_to_genes = self._read_tab_matrix(str(tab_path))
        descriptions = self._read_description(str(disc_path))

        output_path = str(self.organism_dir / f"{species}.Reactome.gmt.gz")
        return self._write_gmt(term_to_genes, descriptions, output_path)

    def generate_do_gmt(self, species: str = "hsa") -> str:
        """Generate GMT output for Disease Ontology."""
        tab_path = self.organism_dir / "hsa.DO2gene.tab.gz"
        disc_path = self.organism_dir / "hsa.DO2disc.gz"

        if not tab_path.exists():
            raise FileNotFoundError(f"DO gene-term matrix does not exist: {tab_path}")
        if not disc_path.exists():
            raise FileNotFoundError(f"The DO profile does not exist: {disc_path}")

        print(f"|---GMTGenerator: Generate DO GMT (SPECies=){species})")

        terms, term_to_genes = self._read_tab_matrix(str(tab_path))
        descriptions = self._read_description(str(disc_path))

        output_path = str(self.organism_dir / f"{species}.DO.gmt.gz")
        return self._write_gmt(term_to_genes, descriptions, output_path)

    def generate_disgenet_gmt(self, species: str = "hsa") -> str:
        """Generate GMT output for DisGeNET."""
        tab_path = self.organism_dir / "hsa.CUI2gene.tab.gz"
        disc_path = self.organism_dir / "hsa.CUI2disc.gz"

        if not tab_path.exists():
            raise FileNotFoundError(f"The DisGeNET gene matrix file does not exist: {tab_path}")
        if not disc_path.exists():
            raise FileNotFoundError(f"DisGeNET Description file does not exist: {disc_path}")

        print(f"--- GMTGENATOR: Generate DisGeNET GMT (SPECies={species})")

        terms, term_to_genes = self._read_tab_matrix(str(tab_path))
        descriptions = self._read_description(str(disc_path))

        output_path = str(self.organism_dir / f"{species}.DisGeNET.gmt.gz")
        return self._write_gmt(term_to_genes, descriptions, output_path)

    def generate_wikipathways_gmt(self, species: str) -> str:
        """Generate GMT output for WikiPathways."""
        tab_path = self.organism_dir / f"{species}.WikiPathways2gene.tab.gz"
        disc_path = self.organism_dir / f"{species}.WikiPathways2disc.gz"

        if not tab_path.exists():
            raise FileNotFoundError(f"The WikiPathways gene matrix file does not exist: {tab_path}")
        if not disc_path.exists():
            raise FileNotFoundError(f"The WikiPathways profile does not exist: {disc_path}")

        print(f"--- GMTGENATOR: Generate WikiPathways GMT (SPECies={species})")

        terms, term_to_genes = self._read_tab_matrix(str(tab_path))
        descriptions = self._read_description(str(disc_path))

        output_path = str(self.organism_dir / f"{species}.WikiPathways.gmt.gz")
        return self._write_gmt(term_to_genes, descriptions, output_path)

    def generate_all_gmt(self, species: str) -> Dict[str, str]:
        """Generate GMT files for every available database."""
        results: Dict[str, str] = {}

        generators = [
            ("GO", lambda: self.generate_go_gmt(species)),
            ("KEGG", lambda: self.generate_kegg_gmt(species)),
            ("Reactome", lambda: self.generate_reactome_gmt(species)),
            ("DO", lambda: self.generate_do_gmt(species)),
            ("DisGeNET", lambda: self.generate_disgenet_gmt(species)),
            ("WikiPathways", lambda: self.generate_wikipathways_gmt(species)),
        ]

        for db_name, gen_func in generators:
            try:
                output_path = gen_func()
                results[db_name] = output_path
            except FileNotFoundError:
                print(f"|---GMTGenerator: Skipped{db_name}(Data file does not exist)")

        print(f"|---GMTGenerator: Generate{len(results)}A GMT file")
        return results
