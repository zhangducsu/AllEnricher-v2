"""
Level-level annotation file resolver unit test
"""

import gzip

import pytest
from pathlib import Path

from allenricher.database.parsers.annotation_parser import (
    AnnotationParser,
    AnnotationRecord,
)


# ============================================================
# AnnotationRecord Test
# ============================================================

class TestAnnotationRecord:
    def test_hierarchy_levels_with_hierarchy(self):
        rec = AnnotationRecord(
            gene="TP53", term_id="GO:0006915", term_name="apoptotic process",
            hierarchy="Biological Process|Cellular Process"
        )
        assert rec.hierarchy_levels == ["Biological Process", "Cellular Process"]

    def test_hierarchy_levels_without_hierarchy(self):
        rec = AnnotationRecord(
            gene="TP53", term_id="GO:0006915", term_name="apoptotic process"
        )
        assert rec.hierarchy_levels == []

    def test_hierarchy_levels_empty_string(self):
        rec = AnnotationRecord(
            gene="TP53", term_id="GO:0006915", term_name="apoptotic process",
            hierarchy=""
        )
        assert rec.hierarchy_levels == []

    def test_repr(self):
        rec = AnnotationRecord(
            gene="TP53", term_id="GO:0006915", term_name="apoptotic process"
        )
        r = repr(rec)
        assert "TP53" in r
        assert "GO:0006915" in r


# ============================================================
# Four Columns Resolution
# ============================================================

class TestParseFourColumn:
    def test_parse_four_column_with_hierarchy(self, tmp_path):
        """Four-column formatting, validation level"""
        f = tmp_path / "annotation.tsv"
        f.write_text(
            "TP53\tGO:0006915\tapoptotic process\tBiological Process|Cellular Process\n"
            "BRCA1\tGO:0006915\tapoptotic process\tBiological Process|Cellular Process\n"
            "EGFR\tGO:0007155\tcell adhesion\tBiological Process|Cell adhesion\n",
            encoding='utf-8'
        )

        parser = AnnotationParser(str(f), format_type='four_column')
        records = parser.parse()

        assert len(records) == 3
        assert records[0].gene == "TP53"
        assert records[0].term_id == "GO:0006915"
        assert records[0].term_name == "apoptotic process"
        assert records[0].hierarchy == "Biological Process|Cellular Process"
        assert records[0].hierarchy_levels == ["Biological Process", "Cellular Process"]


# ============================================================
# 3-column format resolution
# ============================================================

class TestParseThreeColumn:
    def test_parse_three_column(self, tmp_path):
        """3-column format resolution"""
        f = tmp_path / "annotation.tsv"
        f.write_text(
            "TP53\tGO:0006915\tapoptotic process\n"
            "BRCA1\tGO:0006915\tapoptotic process\n"
            "EGFR\tGO:0007155\tcell adhesion\n",
            encoding='utf-8'
        )

        parser = AnnotationParser(str(f), format_type='three_column')
        records = parser.parse()

        assert len(records) == 3
        assert records[0].gene == "TP53"
        assert records[0].term_id == "GO:0006915"
        assert records[0].term_name == "apoptotic process"
        assert records[0].hierarchy is None


# ============================================================
# Two-column format resolution
# ============================================================

class TestParseTwoColumn:
    def test_parse_two_column(self, tmp_path):
        """Two-column format resolution"""
        f = tmp_path / "annotation.tsv"
        f.write_text(
            "TP53\tapoptotic process\n"
            "BRCA1\tapoptotic process\n"
            "EGFR\tcell adhesion\n",
            encoding='utf-8'
        )

        parser = AnnotationParser(str(f), format_type='two_column')
        records = parser.parse()

        assert len(records) == 3
        assert records[0].gene == "TP53"
        # The term_name under two columns is also used as the term_id
        assert records[0].term_id == "apoptotic process"
        assert records[0].term_name == "apoptotic process"
        assert records[0].hierarchy is None


# ============================================================
# get_term_genes test
# ============================================================

class TestGetTermGenes:
    def test_get_term_genes(self, tmp_path):
        """Term to Genes map"""
        f = tmp_path / "annotation.tsv"
        f.write_text(
            "TP53\tGO:0006915\tapoptotic process\n"
            "BRCA1\tGO:0006915\tapoptotic process\n"
            "EGFR\tGO:0007155\tcell adhesion\n",
            encoding='utf-8'
        )

        parser = AnnotationParser(str(f))
        term_genes = parser.get_term_genes()

        assert len(term_genes) == 2
        assert term_genes["GO:0006915"] == {"TP53", "BRCA1"}
        assert term_genes["GO:0007155"] == {"EGFR"}


# ============================================================
# _hierarchy_tree test
# ============================================================

class TestGetHierarchyTree:
    def test_get_hierarchy_tree(self, tmp_path):
        """Tier Tree Structure"""
        f = tmp_path / "annotation.tsv"
        f.write_text(
            "TP53\tGO:0006915\tapoptotic process\tBiological Process|Cellular Process\n"
            "BRCA1\tGO:0006915\tapoptotic process\tBiological Process|Cellular Process\n"
            "EGFR\tGO:0007155\tcell adhesion\tBiological Process|Cell adhesion\n",
            encoding='utf-8'
        )

        parser = AnnotationParser(str(f))
        tree = parser.get_hierarchy_tree()

        # Top level should be Biologicial Process
        assert "Biological Process" in tree
        bp = tree["Biological Process"]
        assert "Cellular Process" in bp
        assert "Cell adhesion" in bp

        # Verify bottom term data
        cp = bp["Cellular Process"]
        assert "GO:0006915" in cp
        assert cp["GO:0006915"]["genes"] == {"TP53", "BRCA1"}
        assert cp["GO:0006915"]["term_name"] == "apoptotic process"

        ca = bp["Cell adhesion"]
        assert "GO:0007155" in ca
        assert ca["GO:0007155"]["genes"] == {"EGFR"}


# ============================================================
# KEGG Level 3 Test
# ============================================================

class TestKEGGStyleHierarchy:
    def test_kegg_style_hierarchy(self, tmp_path):
        """KEG Level 3: Category Subcategory Pathway"""
        f = tmp_path / "kegg_annotation.tsv"
        f.write_text(
            "TP53\thsa04110\tCell Cycle\tMetabolism|Global and overview|Cell Cycle\n"
            "BRCA1\thsa04110\tCell Cycle\tMetabolism|Global and overview|Cell Cycle\n"
            "EGFR\thsa04010\tMAPK signaling pathway\tSignal transduction|Signaling|MAPK signaling pathway\n",
            encoding='utf-8'
        )

        parser = AnnotationParser(str(f))
        tree = parser.get_hierarchy_tree()

        # Level 3: Heerarchy path with 3 layers, term_id under 3
        assert "Metabolism" in tree
        metabolism = tree["Metabolism"]
        assert "Global and overview" in metabolism
        overview = metabolism["Global and overview"]
        # The last layer of "Cell Cycle" will be the node of the tree.
        assert "Cell Cycle" in overview
        cell_cycle = overview["Cell Cycle"]
        assert "hsa04110" in cell_cycle
        assert cell_cycle["hsa04110"]["genes"] == {"TP53", "BRCA1"}

        assert "Signal transduction" in tree
        st = tree["Signal transduction"]
        assert "Signaling" in st
        signaling = st["Signaling"]
        assert "MAPK signaling pathway" in signaling
        mapk = signaling["MAPK signaling pathway"]
        assert "hsa04010" in mapk
        assert mapk["hsa04010"]["genes"] == {"EGFR"}


# ============================================================
# GO Level 2 Test
# ============================================================

class TestGOStyleHierarchy:
    def test_go_style_hierarchy(self, tmp_path):
        """GO Level 2: OntologicTerm"""
        f = tmp_path / "go_annotation.tsv"
        f.write_text(
            "TP53\tGO:0006915\tapoptotic process\tBiological Process|apoptotic process\n"
            "BRCA1\tGO:0006915\tapoptotic process\tBiological Process|apoptotic process\n"
            "EGFR\tGO:0005576\textracellular region\tCellular Component|extracellular region\n",
            encoding='utf-8'
        )

        parser = AnnotationParser(str(f))
        tree = parser.get_hierarchy_tree()

        # The path to the hierarchy is 2 floor, term_id below 2 floor
        assert "Biological Process" in tree
        bp = tree["Biological Process"]
        assert "apoptotic process" in bp
        ap = bp["apoptotic process"]
        assert "GO:0006915" in ap
        assert ap["GO:0006915"]["genes"] == {"TP53", "BRCA1"}

        assert "Cellular Component" in tree
        cc = tree["Cellular Component"]
        assert "extracellular region" in cc
        er = cc["extracellular region"]
        assert "GO:0005576" in er
        assert er["GO:0005576"]["genes"] == {"EGFR"}


# ============================================================
# AutoFormat Test Test
# ============================================================

class TestAutoDetectFormat:
    def test_auto_detect_four_column(self, tmp_path):
        f = tmp_path / "four_col.tsv"
        f.write_text(
            "TP53\tGO:0006915\tapoptotic process\tBiological Process\n"
        )
        parser = AnnotationParser(str(f))
        assert parser._detect_format() == 'four_column'

    def test_auto_detect_three_column(self, tmp_path):
        f = tmp_path / "three_col.tsv"
        f.write_text(
            "TP53\tGO:0006915\tapoptotic process\n"
        )
        parser = AnnotationParser(str(f))
        assert parser._detect_format() == 'three_column'

    def test_auto_detect_two_column(self, tmp_path):
        f = tmp_path / "two_col.tsv"
        f.write_text(
            "TP53\tapoptotic process\n"
        )
        parser = AnnotationParser(str(f))
        assert parser._detect_format() == 'two_column'

    def test_auto_detect_with_comments(self, tmp_path):
        """Ignore comment lines while detecting the table format."""
        f = tmp_path / "with_comments.tsv"
        f.write_text(
            "# this is a comment\n"
            "# another comment\n"
            "TP53\tGO:0006915\tapoptotic process\n"
        )
        parser = AnnotationParser(str(f))
        assert parser._detect_format() == 'three_column'


# ============================================================
# Empty lines and comments
# ============================================================

class TestEmptyLinesAndComments:
    def test_empty_lines_and_comments(self, tmp_path):
        """Ignore blank lines and comments while parsing records."""
        f = tmp_path / "mixed.tsv"
        f.write_text(
            "# header comment\n"
            "\n"
            "TP53\tGO:0006915\tapoptotic process\n"
            "\n"
            "# inline comment\n"
            "BRCA1\tGO:0006915\tapoptotic process\n"
            "\n"
        )
        parser = AnnotationParser(str(f))
        records = parser.parse()

        assert len(records) == 2
        assert records[0].gene == "TP53"
        assert records[1].gene == "BRCA1"

    def test_empty_file_raises(self, tmp_path):
        """Reject an empty annotation file."""
        f = tmp_path / "empty.tsv"
        f.write_text("", encoding='utf-8')

        parser = AnnotationParser(str(f))
        with pytest.raises(ValueError, match="empty or contains only"):
            parser._detect_format()

    def test_only_comments_raises(self, tmp_path):
        """Reject an annotation file containing only comments."""
        f = tmp_path / "comments_only.tsv"
        f.write_text("# comment1\n# comment2\n", encoding='utf-8')

        parser = AnnotationParser(str(f))
        with pytest.raises(ValueError, match="empty or contains only"):
            parser._detect_format()

    def test_gzip_header_and_custom_hierarchy_separator(self, tmp_path):
        path = tmp_path / "annotation.tsv.gz"
        with gzip.open(path, "wt", encoding="utf-8") as handle:
            handle.write("gene\tterm_id\tterm_name\thierarchy\n")
            handle.write(" TP53 \t GO:0006915 \t Apoptosis \t Biology > Cell death \n")

        records = AnnotationParser(
            str(path), hierarchy_separator=">"
        ).parse()

        assert len(records) == 1
        assert records[0].gene == "TP53"
        assert records[0].hierarchy == "Biology|Cell death"


# ============================================================
# File does not have a test
# ============================================================

class TestFileNotFound:
    def test_file_not_found(self):
        """File does not have an anomaly"""
        parser = AnnotationParser("/nonexistent/path/file.tsv")
        with pytest.raises(FileNotFoundError, match="Annotation file does not exist"):
            parser.parse()

    def test_detect_format_file_not_found(self):
        """_Defect_format file does not appear abnormal"""
        parser = AnnotationParser("/nonexistent/path/file.tsv")
        with pytest.raises(FileNotFoundError, match="Annotation file does not exist"):
            parser._detect_format()


# ============================================================
# Resolve result cache test
# ============================================================

class TestParseCaching:
    def test_parse_caching(self, tmp_path):
        """Repeated call parse() should return the cache result"""
        f = tmp_path / "cache.tsv"
        f.write_text("TP53\tGO:0006915\tapoptotic process\n", encoding='utf-8')

        parser = AnnotationParser(str(f))
        records1 = parser.parse()
        records2 = parser.parse()
        assert records1 is records2


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
