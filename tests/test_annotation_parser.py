"""
层级注释文件解析器单元测试
"""

import pytest
from pathlib import Path

from allenricher.database.parsers.annotation_parser import (
    AnnotationParser,
    AnnotationRecord,
)


# ============================================================
# AnnotationRecord 测试
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
# 四列格式解析
# ============================================================

class TestParseFourColumn:
    def test_parse_four_column_with_hierarchy(self, tmp_path):
        """四列格式解析，验证层级"""
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
# 三列格式解析
# ============================================================

class TestParseThreeColumn:
    def test_parse_three_column(self, tmp_path):
        """三列格式解析"""
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
# 两列格式解析
# ============================================================

class TestParseTwoColumn:
    def test_parse_two_column(self, tmp_path):
        """两列格式解析"""
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
        # 两列格式下 term_name 同时作为 term_id
        assert records[0].term_id == "apoptotic process"
        assert records[0].term_name == "apoptotic process"
        assert records[0].hierarchy is None


# ============================================================
# get_term_genes 测试
# ============================================================

class TestGetTermGenes:
    def test_get_term_genes(self, tmp_path):
        """term 到 genes 映射"""
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
# get_hierarchy_tree 测试
# ============================================================

class TestGetHierarchyTree:
    def test_get_hierarchy_tree(self, tmp_path):
        """层级树结构"""
        f = tmp_path / "annotation.tsv"
        f.write_text(
            "TP53\tGO:0006915\tapoptotic process\tBiological Process|Cellular Process\n"
            "BRCA1\tGO:0006915\tapoptotic process\tBiological Process|Cellular Process\n"
            "EGFR\tGO:0007155\tcell adhesion\tBiological Process|Cell adhesion\n",
            encoding='utf-8'
        )

        parser = AnnotationParser(str(f))
        tree = parser.get_hierarchy_tree()

        # 顶层应有 Biological Process
        assert "Biological Process" in tree
        bp = tree["Biological Process"]
        assert "Cellular Process" in bp
        assert "Cell adhesion" in bp

        # 验证最底层 term 数据
        cp = bp["Cellular Process"]
        assert "GO:0006915" in cp
        assert cp["GO:0006915"]["genes"] == {"TP53", "BRCA1"}
        assert cp["GO:0006915"]["term_name"] == "apoptotic process"

        ca = bp["Cell adhesion"]
        assert "GO:0007155" in ca
        assert ca["GO:0007155"]["genes"] == {"EGFR"}


# ============================================================
# KEGG 三级层级测试
# ============================================================

class TestKEGGStyleHierarchy:
    def test_kegg_style_hierarchy(self, tmp_path):
        """KEGG 三级层级: Category|Subcategory|Pathway"""
        f = tmp_path / "kegg_annotation.tsv"
        f.write_text(
            "TP53\thsa04110\tCell Cycle\tMetabolism|Global and overview|Cell Cycle\n"
            "BRCA1\thsa04110\tCell Cycle\tMetabolism|Global and overview|Cell Cycle\n"
            "EGFR\thsa04010\tMAPK signaling pathway\tSignal transduction|Signaling|MAPK signaling pathway\n",
            encoding='utf-8'
        )

        parser = AnnotationParser(str(f))
        tree = parser.get_hierarchy_tree()

        # 三级层级：hierarchy 路径有 3 层，term_id 在第 3 层之下
        assert "Metabolism" in tree
        metabolism = tree["Metabolism"]
        assert "Global and overview" in metabolism
        overview = metabolism["Global and overview"]
        # hierarchy 最后一层 "Cell Cycle" 也会作为树的节点
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
# GO 两级层级测试
# ============================================================

class TestGOStyleHierarchy:
    def test_go_style_hierarchy(self, tmp_path):
        """GO 两级层级: Ontology|Term"""
        f = tmp_path / "go_annotation.tsv"
        f.write_text(
            "TP53\tGO:0006915\tapoptotic process\tBiological Process|apoptotic process\n"
            "BRCA1\tGO:0006915\tapoptotic process\tBiological Process|apoptotic process\n"
            "EGFR\tGO:0005576\textracellular region\tCellular Component|extracellular region\n",
            encoding='utf-8'
        )

        parser = AnnotationParser(str(f))
        tree = parser.get_hierarchy_tree()

        # hierarchy 路径有 2 层，term_id 在第 2 层之下
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
# 自动格式检测测试
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
        """自动检测应跳过注释行"""
        f = tmp_path / "with_comments.tsv"
        f.write_text(
            "# this is a comment\n"
            "# another comment\n"
            "TP53\tGO:0006915\tapoptotic process\n"
        )
        parser = AnnotationParser(str(f))
        assert parser._detect_format() == 'three_column'


# ============================================================
# 空行和注释处理测试
# ============================================================

class TestEmptyLinesAndComments:
    def test_empty_lines_and_comments(self, tmp_path):
        """空行和注释处理"""
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
        """空文件应抛出 ValueError"""
        f = tmp_path / "empty.tsv"
        f.write_text("", encoding='utf-8')

        parser = AnnotationParser(str(f))
        with pytest.raises(ValueError, match="为空或仅包含注释行"):
            parser._detect_format()

    def test_only_comments_raises(self, tmp_path):
        """仅包含注释行的文件应抛出 ValueError"""
        f = tmp_path / "comments_only.tsv"
        f.write_text("# comment1\n# comment2\n", encoding='utf-8')

        parser = AnnotationParser(str(f))
        with pytest.raises(ValueError, match="为空或仅包含注释行"):
            parser._detect_format()


# ============================================================
# 文件不存在测试
# ============================================================

class TestFileNotFound:
    def test_file_not_found(self):
        """文件不存在异常"""
        parser = AnnotationParser("/nonexistent/path/file.tsv")
        with pytest.raises(FileNotFoundError, match="注释文件不存在"):
            parser.parse()

    def test_detect_format_file_not_found(self):
        """_detect_format 文件不存在异常"""
        parser = AnnotationParser("/nonexistent/path/file.tsv")
        with pytest.raises(FileNotFoundError, match="注释文件不存在"):
            parser._detect_format()


# ============================================================
# 解析结果缓存测试
# ============================================================

class TestParseCaching:
    def test_parse_caching(self, tmp_path):
        """多次调用 parse() 应返回缓存结果"""
        f = tmp_path / "cache.tsv"
        f.write_text("TP53\tGO:0006915\tapoptotic process\n", encoding='utf-8')

        parser = AnnotationParser(str(f))
        records1 = parser.parse()
        records2 = parser.parse()
        assert records1 is records2


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
