"""
数据库构建模块单元测试（v2 新架构）

新架构：
  download → database/basic/{type}/{ver}/ (全体物种通用数据)
  build    → database/organism/v{date}/{species}/ (指定物种格式化数据)
"""

import pytest
import gzip
import tempfile
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).parent.parent))

from allenricher.database.parsers.go import GOParser
from allenricher.database.parsers.kegg import KEGGParser
from allenricher.database.parsers.reactome import ReactomeParser
from allenricher.database.parsers.do import DOParser
from allenricher.database.parsers.disgenet import DisGeNETParser
from allenricher.database.downloader import DataDownloader
from allenricher.database.builder import DatabaseBuilder


# ============================================================
# Test Helpers
# ============================================================

def _create_mock_go_basic(root: Path, version: str = "GO20250101"):
    """创建模拟的 GO 全体物种基础数据目录"""
    go_dir = root / "basic" / "go" / version
    go_dir.mkdir(parents=True)

    with gzip.open(go_dir / "gene2go.gz", 'wt') as f:
        f.write("9606\t1\tGO:0005576\t\t\textracellular region\t\t\tcellular_component\n")
        f.write("9606\t2\tGO:0005576\t\t\textracellular region\t\t\tcellular_component\n")
        f.write("9606\t1\tGO:0051301\t\t\tcell division\t\t\tbiological_process\n")
        f.write("9606\t3\tGO:0051301\t\t\tcell division\t\t\tbiological_process\n")
        f.write("10090\t4\tGO:0005576\t\t\textracellular region\t\t\tcellular_component\n")

    with gzip.open(go_dir / "gene_info.gz", 'wt') as f:
        f.write("9606\t1\tGENE_A\t\t\t\t\t\t\t\t\n")
        f.write("9606\t2\tGENE_B\t\t\t\t\t\t\t\t\n")
        f.write("9606\t3\tGENE_C\t\t\t\t\t\t\t\t\n")

    with open(go_dir / "go-basic.obo", 'w') as f:
        f.write("format-version: 1.2\n")
        f.write("[Term]\nid: GO:0005576\nname: extracellular region\nnamespace: cellular_component\n")
        f.write("is_a: GO:0005615 ! extracellular space\n\n")
        f.write("[Term]\nid: GO:0051301\nname: cell division\nnamespace: biological_process\n")
    return go_dir


def _create_mock_reactome_basic(root: Path, version: str = "Reactome20250101"):
    """创建模拟的 Reactome 全体物种基础数据目录"""
    re_dir = root / "basic" / "reactome" / version
    re_dir.mkdir(parents=True)

    with gzip.open(re_dir / "NCBI2Reactome_All_Levels.txt.gz", 'wt') as f:
        f.write("1\tR-HSA-12345\tPathway Name 1\turl1\n")
        f.write("2\tR-HSA-12345\tPathway Name 1\turl2\n")
        f.write("1\tR-HSA-67890\tPathway Name 2\turl3\n")
        f.write("3\tR-MMU-12345\tPathway Name 1\turl4\n")

    with gzip.open(re_dir / "gene_info.gz", 'wt') as f:
        f.write("9606\t1\tGENE_A\t\t\t\t\t\t\t\t\n")
        f.write("9606\t2\tGENE_B\t\t\t\t\t\t\t\t\n")
    return re_dir


# ============================================================
# 解析器测试（不变）
# ============================================================

class TestGOParser:
    def test_parse_gene2go(self, tmp_path):
        gene2go_path = tmp_path / "gene2go.gz"
        with gzip.open(gene2go_path, 'wt') as f:
            f.write("9606\t1\tGO:0005576\t\t\textracellular region\t\t\tcellular_component\n")
            f.write("9606\t2\tGO:0005576\t\t\textracellular region\t\t\tcellular_component\n")
            f.write("9606\t1\tGO:0051301\t\t\tcell division\t\t\tbiological_process\n")
            f.write("10090\t4\tGO:0005576\t\t\textracellular region\t\t\tcellular_component\n")

        gene_info_path = tmp_path / "gene_info.gz"
        with gzip.open(gene_info_path, 'wt') as f:
            f.write("9606\t1\tGENE_A\t\t\t\t\t\t\t\t\n")
            f.write("9606\t2\tGENE_B\t\t\t\t\t\t\t\t\n")

        parser = GOParser()
        parser.parse_gene2go(str(gene2go_path), str(gene_info_path), 9606, "hsa", str(tmp_path))

        output_tab = tmp_path / "hsa.GO2gene.tab.gz"
        assert output_tab.exists()
        with gzip.open(output_tab, 'rt') as f:
            lines = f.readlines()
            header = lines[0].strip().split('\t')
            assert header[0] == 'Gene'
            assert 'GO:0005576' in header
            assert 'GO:0051301' in header

    def test_parse_obo(self, tmp_path):
        obo_path = tmp_path / "go-basic.obo"
        with open(obo_path, 'w') as f:
            f.write("format-version: 1.2\n")
            f.write("[Term]\nid: GO:0005576\nname: extracellular region\nnamespace: cellular_component\n")
            f.write("is_a: GO:0005615 ! extracellular space\n")

        parser = GOParser()
        parser.parse_obo(str(obo_path), str(tmp_path))
        output_disc = tmp_path / "GO2disc.gz"
        assert output_disc.exists()
        with gzip.open(output_disc, 'rt') as f:
            lines = f.readlines()
            assert 'GO:0005576' in lines[0]


class TestKEGGParser:
    def test_build_database(self, tmp_path):
        gene_info_path = tmp_path / "gene_info.gz"
        with gzip.open(gene_info_path, 'wt') as f:
            f.write("9606\t1\tGENE_A\t\t\t\t\t\t\t\t\n")
            f.write("9606\t2\tGENE_B\t\t\t\t\t\t\t\t\n")

        gene2pathway_path = tmp_path / "gene2pathway.txt"
        with open(gene2pathway_path, 'w') as f:
            f.write("GENE_A\t1\thsa04110\tCell Cycle\n")
            f.write("GENE_B\t2\thsa04110\tCell Cycle\n")

        pathway_summary_path = tmp_path / "pathway_summary.txt"
        with open(pathway_summary_path, 'w') as f:
            f.write("Metabolism\tGlobal and overview maps\t04110\tCell Cycle\n")

        parser = KEGGParser()
        parser.build_database(
            species="hsa",
            gene_info_path=str(gene_info_path),
            gene2pathway_path=str(gene2pathway_path),
            outdir=str(tmp_path),
            pathway_summary_path=str(pathway_summary_path)
        )

        output_tab = tmp_path / "hsa.kegg2gene.tab.gz"
        assert output_tab.exists()
        with gzip.open(output_tab, 'rt') as f:
            header = f.readline().strip().split('\t')
            assert 'hsa04110' in header


class TestReactomeParser:
    def test_parse_ncbi2reactome(self, tmp_path):
        ncbi2reactome_path = tmp_path / "NCBI2Reactome.txt.gz"
        with gzip.open(ncbi2reactome_path, 'wt') as f:
            f.write("1\tR-HSA-12345\tPathway Name 1\turl1\n")
            f.write("2\tR-HSA-12345\tPathway Name 1\turl2\n")
            f.write("3\tR-MMU-12345\tPathway Name 1\turl4\n")

        gene_info_path = tmp_path / "gene_info.gz"
        with gzip.open(gene_info_path, 'wt') as f:
            f.write("9606\t1\tGENE_A\t\t\t\t\t\t\t\t\n")
            f.write("9606\t2\tGENE_B\t\t\t\t\t\t\t\t\n")

        parser = ReactomeParser()
        parser.parse_ncbi2reactome(str(ncbi2reactome_path), str(gene_info_path), 9606, "hsa", str(tmp_path))

        output_tab = tmp_path / "hsa.Reactome2gene.tab.gz"
        assert output_tab.exists()


class TestDOParser:
    def test_parse_disease_files(self, tmp_path):
        disease_file = tmp_path / "human_disease_knowledge.tsv"
        with open(disease_file, 'w') as f:
            f.write("col0\tGENE_A\tDOID:1234\tBreast Cancer\tcol4\n")
            f.write("col0\tGENE_B\tDOID:5678\tLung Cancer\tcol4\n")

        gene_info_path = tmp_path / "gene_info.gz"
        with gzip.open(gene_info_path, 'wt') as f:
            f.write("9606\t1\tGENE_A\t\t\t\t\t\t\t\t\n")
            f.write("9606\t2\tGENE_B\t\t\t\t\t\t\t\t\n")

        parser = DOParser()
        parser.parse_disease_files([str(disease_file)], str(gene_info_path), 9606, str(tmp_path))
        assert (tmp_path / "hsa.DO2gene.tab.gz").exists()
        assert (tmp_path / "hsa.DO2disc.gz").exists()


class TestDisGeNETParser:
    def test_parse_associations(self, tmp_path):
        assoc_path = tmp_path / "associations.tsv.gz"
        with gzip.open(assoc_path, 'wt') as f:
            f.write("GENE_A\t1\t\t\tCUI:1234\tBreast Cancer\n")
            f.write("GENE_B\t2\t\t\tCUI:5678\tLung Cancer\n")

        gene_info_path = tmp_path / "gene_info.gz"
        with gzip.open(gene_info_path, 'wt') as f:
            f.write("9606\t1\tGENE_A\t\t\t\t\t\t\t\t\n")
            f.write("9606\t2\tGENE_B\t\t\t\t\t\t\t\t\n")

        parser = DisGeNETParser()
        parser.parse_associations(str(assoc_path), str(gene_info_path), 9606, str(tmp_path))
        assert (tmp_path / "hsa.CUI2gene.tab.gz").exists()


# ============================================================
# 新架构测试：download → build 完整流程
# ============================================================

class TestDatabaseBuilderNew:
    """测试新的 DatabaseBuilder（v2 架构）"""

    def test_build_go_from_basic(self, tmp_path):
        """从 database/basic/ 构建 GO 数据库"""
        root = Path(tmp_path) / "database"
        _create_mock_go_basic(root)

        builder = DatabaseBuilder(root_dir=str(root))
        outdir = builder.build_go(species="hsa", taxid=9606)

        # 验证输出位置：database/organism/v{date}/hsa/
        assert outdir.startswith(str(root / "organism"))
        assert "hsa" in outdir

        out_path = Path(outdir)
        assert (out_path / "hsa.GO2gene.tab.gz").exists()
        assert (out_path / "GO2disc.gz").exists()
        assert (out_path / "hsa.gene2go.txt").exists()

    def test_build_reactome_from_basic(self, tmp_path):
        """从 database/basic/ 构建 Reactome 数据库"""
        root = Path(tmp_path) / "database"
        _create_mock_reactome_basic(root)

        builder = DatabaseBuilder(root_dir=str(root))
        outdir = builder.build_reactome(species="hsa", taxid=9606)

        out_path = Path(outdir)
        assert (out_path / "hsa.Reactome2gene.tab.gz").exists()
        assert (out_path / "hsa.Reactome2disc.gz").exists()

    def test_build_species_db(self, tmp_path):
        """一键构建物种数据库"""
        root = Path(tmp_path) / "database"
        _create_mock_go_basic(root)
        _create_mock_reactome_basic(root)

        builder = DatabaseBuilder(root_dir=str(root))
        outdir = builder.build_species_db(
            species="hsa", taxid=9606,
            databases=["GO", "Reactome"]
        )

        out_path = Path(outdir)
        assert (out_path / "hsa.GO2gene.tab.gz").exists()
        assert (out_path / "GO2disc.gz").exists()
        assert (out_path / "hsa.Reactome2gene.tab.gz").exists()

    def test_build_no_basic_data(self, tmp_path):
        """测试没有基础数据时的错误提示"""
        root = Path(tmp_path) / "database"
        root.mkdir(parents=True)
        (root / "basic").mkdir()

        builder = DatabaseBuilder(root_dir=str(root))
        with pytest.raises(FileNotFoundError):
            builder.build_go(species="hsa", taxid=9606)

    def test_skip_disgenet_for_non_human(self, tmp_path):
        """测试非人类物种自动跳过 DisGeNET"""
        root = Path(tmp_path) / "database"
        _create_mock_go_basic(root)

        builder = DatabaseBuilder(root_dir=str(root))
        # mmu（小鼠）不会触发 DisGeNET 构建，不应报错
        outdir = builder.build_species_db(
            species="mmu", taxid=10090,
            databases=["GO", "DisGeNET"]
        )
        # 应该正常完成（DisGeNET 被跳过）
        assert Path(outdir).exists()


class TestDataDownloaderNew:
    """测试新的 DataDownloader"""

    def test_init_creates_dirs(self):
        downloader = DataDownloader(root_dir="/tmp/test_downloader")
        assert downloader.basic_dir == Path("/tmp/test_downloader/basic")

    def test_version_listing(self, tmp_path):
        """测试版本列表功能"""
        root = Path(tmp_path) / "database"
        _create_mock_go_basic(root)

        downloader = DataDownloader(root_dir=str(root))
        versions = downloader.list_go_versions()
        assert "GO20250101" in versions
        assert downloader.get_latest_go_version() == "GO20250101"

    def test_no_versions(self, tmp_path):
        """测试无版本时返回空"""
        downloader = DataDownloader(root_dir=str(tmp_path))
        assert downloader.get_latest_go_version() is None
        assert downloader.list_go_versions() == []


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
