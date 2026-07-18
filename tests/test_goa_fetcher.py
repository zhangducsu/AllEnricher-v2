import pytest
import gzip
from pathlib import Path
from allenricher.database.goa_fetcher import GOAFetcher


class TestGOAFetcher:
    """GOAFetcher Test"""

    @pytest.fixture
    def sample_gaf_file(self, tmp_path):
        """Creates a 17 column simulation file that corresponds to GAF2.2, with Taxon at the 12 below."""
        gaf_path = tmp_path / "test.gaf"
        rows = [
            ["UniProtKB", "P12345", "GOT2", "involved_in", "GO:0006457", "PMID:1", "IEA", "", "P", "GOT2", "mAspAT", "protein", "taxon:9606", "20260101", "UniProtKB", "", ""],
            ["UniProtKB", "P12345", "GOT2", "located_in", "GO:0005739", "PMID:1", "IEA", "", "C", "GOT2", "mAspAT", "protein", "taxon:9606", "20260101", "UniProtKB", "", ""],
            ["UniProtKB", "P67890", "TP53", "enables", "GO:0003677", "PMID:2", "EXP", "", "F", "TP53", "p53", "protein", "taxon:9606", "20260101", "HGNC", "", ""],
            ["UniProtKB", "P67890", "TP53", "NOT|enables", "GO:9999999", "PMID:2", "EXP", "", "F", "TP53", "p53", "protein", "taxon:9606", "20260101", "HGNC", "", ""],
            ["UniProtKB", "Q9Y6K1", "BRCA1", "involved_in", "GO:0006974", "PMID:3", "IDA", "", "P", "BRCA1", "", "protein", "taxon:9606", "20260101", "HGNC", "", ""],
            ["UniProtKB", "P62988", "RPL10A", "enables", "GO:0003735", "GOREF:3", "IEA", "", "F", "RPL10A", "", "protein", "taxon:9606", "20260101", "UniProtKB", "", ""],
            ["UniProtKB", "P62988", "RPL10A", "involved_in", "GO:0006412", "GOREF:3", "IEA", "", "P", "RPL10A", "", "protein", "taxon:9606", "20260101", "UniProtKB", "", ""],
            ["UniProtKB", "P99999", "UNKNOWN", "enables", "GO:0003674", "GOREF:3", "IEA", "", "F", "Unknown", "", "protein", "taxon:10090", "20260101", "UniProtKB", "", ""],
            ["UniProtKB", "P11111", "P11111", "enables", "GO:0003674", "GOREF:3", "IEA", "", "F", "", "", "protein", "taxon:9606", "20260101", "UniProtKB", "", ""],
            ["UniProtKB", "Q00001", "7SK", "enables", "GO:0003674", "GOREF:3", "IEA", "", "F", "7SK", "", "ncRNA", "taxon:9606|taxon:10090", "20260101", "UniProtKB", "", ""],
        ]
        gaf_path.write_text(
            "!gaf-version: 2.2\n" + "\n".join("\t".join(row) for row in rows) + "\n",
            encoding="utf-8",
        )
        return gaf_path

    def test_parse_goa_file_filters_by_taxid(self, sample_gaf_file):
        """Test filter with TaxID"""
        fetcher = GOAFetcher(cache_dir="/tmp/test_goa")
        gene_to_go, all_genes = fetcher.parse_goa_file(sample_gaf_file, taxid=9606)

        assert len(all_genes) == 6
        assert "GOT2" in all_genes
        assert "TP53" in all_genes
        assert "BRCA1" in all_genes
        assert "RPL10A" in all_genes
        assert "UNKNOWN" not in all_genes

    def test_parse_goa_file_preserves_valid_symbols(self, sample_gaf_file):
        """The GAF allows Symbol to equal the object ID and also allows valid symbols at the beginning of the number."""
        fetcher = GOAFetcher(cache_dir="/tmp/test_goa")
        gene_to_go, all_genes = fetcher.parse_goa_file(sample_gaf_file, taxid=9606)

        assert "P11111" in all_genes
        assert "7SK" in all_genes

    def test_parse_goa_file_gene_to_go_mapping(self, sample_gaf_file):
        """Test the gene to the GO map."""
        fetcher = GOAFetcher(cache_dir="/tmp/test_goa")
        gene_to_go, all_genes = fetcher.parse_goa_file(sample_gaf_file, taxid=9606)

        # GoT2 should have two Go term
        assert len(gene_to_go.get("GOT2", set())) == 2
        assert "GO:0006457" in gene_to_go["GOT2"]
        assert "GO:0005739" in gene_to_go["GOT2"]

        # RPL10A should have two Go term
        assert len(gene_to_go.get("RPL10A", set())) == 2
        assert "GO:9999999" not in gene_to_go["TP53"]

    def test_build_go2gene_matrix(self, tmp_path):
        """Test Go2ge Matrix Generation"""
        gene_to_go = {
            "GOT2": {"GO:0006457", "GO:0005739"},
            "TP53": {"GO:0003677"},
            "BRCA1": {"GO:0006974"},
        }
        all_genes = {"GOT2", "TP53", "BRCA1"}
        all_go_terms = {"GO:0006457", "GO:0003677", "GO:0005739", "GO:0006974"}
        output_path = tmp_path / "GO2gene.tab.gz"

        GOAFetcher.build_go2gene_matrix(gene_to_go, all_genes, all_go_terms, output_path)

        assert output_path.exists()
        # Read and authenticate
        with gzip.open(output_path, 'rt') as f:
            header = f.readline().strip().split('\t')
            assert header[0] == "Gene"
            assert "GO:0006457" in header

            lines = f.readlines()
            assert len(lines) == 3  # Three genes.

            # GOT2 lines
            got2_line = [l for l in lines if l.startswith("GOT2")][0]
            fields = got2_line.strip().split('\t')
            got2_idx = header.index("GO:0006457")
            assert fields[got2_idx] == "1"

    def test_build_gene2go_list(self, tmp_path):
        """Testing gene2go list generation"""
        gene_to_go = {
            "GOT2": {"GO:0006457", "GO:0005739"},
            "TP53": {"GO:0003677"},
        }
        go_names = {
            "GO:0006457": "protein binding",
            "GO:0005739": "mitochondrion",
            "GO:0003677": "DNA binding",
        }
        output_path = tmp_path / "gene2go.txt"

        GOAFetcher.build_gene2go_list(gene_to_go, go_names, output_path)

        assert output_path.exists()
        content = output_path.read_text()
        lines = content.strip().split('\n')
        assert len(lines) == 3  # 3 notes

        # GoT2 should have two lines.
        got2_lines = [l for l in lines if l.startswith("GOT2")]
        assert len(got2_lines) == 2
