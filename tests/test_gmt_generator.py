"""
GMT 基因集文件生成器单元测试

测试 GMTGenerator 的核心功能：
- GMT 文件格式正确性（至少3列，第一列名称，第二列描述，后续为基因）
- 压缩文件可正确读取
- 空数据处理
- 与 DatabaseBuilder 的集成
"""

import pytest
import gzip
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).parent.parent))

from allenricher.database.gmt_generator import GMTGenerator
from allenricher.database.builder import DatabaseBuilder


# ============================================================
# 测试辅助函数
# ============================================================

def _create_mock_go_db(org_dir: Path, species: str = "hsa"):
    """创建模拟的 GO 数据库产物文件"""
    org_dir.mkdir(parents=True, exist_ok=True)

    # {species}.GO2gene.tab.gz: Gene\tGO_ID1\tGO_ID2\t...
    tab_path = org_dir / f"{species}.GO2gene.tab.gz"
    with gzip.open(tab_path, 'wt', encoding='utf-8') as f:
        f.write("Gene\tGO:0005576\tGO:0051301\n")
        f.write("GENE_A\t1\t0\n")
        f.write("GENE_B\t1\t1\n")
        f.write("GENE_C\t0\t1\n")

    # GO2disc.gz: GO_ID\tnamespace:name\tfather1;father2;...
    disc_path = org_dir / "GO2disc.gz"
    with gzip.open(disc_path, 'wt', encoding='utf-8') as f:
        f.write("GO:0005576\tcellular_component:extracellular_region\tGO:0005615\n")
        f.write("GO:0051301\tbiological_process:cell_division\t\n")


def _create_mock_kegg_db(org_dir: Path, species: str = "hsa"):
    """创建模拟的 KEGG 数据库产物文件"""
    org_dir.mkdir(parents=True, exist_ok=True)

    tab_path = org_dir / f"{species}.kegg2gene.tab.gz"
    with gzip.open(tab_path, 'wt', encoding='utf-8') as f:
        f.write("Gene\thsa04110\thsa04150\n")
        f.write("GENE_A\t1\t0\n")
        f.write("GENE_B\t1\t1\n")

    disc_path = org_dir / f"{species}.kegg2disc.gz"
    with gzip.open(disc_path, 'wt', encoding='utf-8') as f:
        f.write("hsa04110\tMetabolism|Global_and_overview_maps|Cell_Cycle\n")
        f.write("hsa04150\tMetabolism|Lipid_metabolism|PPAR_signaling_pathway\n")


def _create_mock_reactome_db(org_dir: Path, species: str = "hsa"):
    """创建模拟的 Reactome 数据库产物文件"""
    org_dir.mkdir(parents=True, exist_ok=True)

    tab_path = org_dir / f"{species}.Reactome2gene.tab.gz"
    with gzip.open(tab_path, 'wt', encoding='utf-8') as f:
        f.write("Gene\tR-HSA-12345\tR-HSA-67890\n")
        f.write("GENE_A\t1\t0\n")
        f.write("GENE_B\t1\t1\n")

    disc_path = org_dir / f"{species}.Reactome2disc.gz"
    with gzip.open(disc_path, 'wt', encoding='utf-8') as f:
        f.write("R-HSA-12345\tSignal_Transduction\n")
        f.write("R-HSA-67890\tImmune_System\n")


def _create_mock_do_db(org_dir: Path):
    """创建模拟的 DO 数据库产物文件"""
    org_dir.mkdir(parents=True, exist_ok=True)

    tab_path = org_dir / "hsa.DO2gene.tab.gz"
    with gzip.open(tab_path, 'wt', encoding='utf-8') as f:
        f.write("Gene\tDOID:1234\tDOID:5678\n")
        f.write("GENE_A\t1\t0\n")
        f.write("GENE_B\t0\t1\n")

    disc_path = org_dir / "hsa.DO2disc.gz"
    with gzip.open(disc_path, 'wt', encoding='utf-8') as f:
        f.write("DOID:1234\tBreast_Cancer\n")
        f.write("DOID:5678\tLung_Cancer\n")


def _create_mock_disgenet_db(org_dir: Path):
    """创建模拟的 DisGeNET 数据库产物文件"""
    org_dir.mkdir(parents=True, exist_ok=True)

    tab_path = org_dir / "hsa.CUI2gene.tab.gz"
    with gzip.open(tab_path, 'wt', encoding='utf-8') as f:
        f.write("Gene\tCUI:0001\tCUI:0002\n")
        f.write("GENE_A\t1\t0\n")
        f.write("GENE_B\t1\t1\n")

    disc_path = org_dir / "hsa.CUI2disc.gz"
    with gzip.open(disc_path, 'wt', encoding='utf-8') as f:
        f.write("CUI:0001\tBreast_Cancer\n")
        f.write("CUI:0002\tLung_Cancer\n")


def _create_mock_full_db(org_dir: Path, species: str = "hsa"):
    """创建模拟的完整物种数据库（所有类型）"""
    _create_mock_go_db(org_dir, species)
    _create_mock_kegg_db(org_dir, species)
    _create_mock_reactome_db(org_dir, species)
    if species == "hsa":
        _create_mock_do_db(org_dir)
        _create_mock_disgenet_db(org_dir)


# ============================================================
# GMT 格式验证辅助
# ============================================================

def _validate_gmt_format(filepath: str):
    """验证 GMT 文件格式正确性

    Args:
        filepath: .gmt.gz 文件路径

    Returns:
        list: 解析后的所有行数据 [[name, desc, gene1, gene2, ...], ...]
    """
    rows = []
    with gzip.open(filepath, 'rt', encoding='utf-8') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split('\t')
            # 至少3列：名称、描述、至少一个基因
            assert len(parts) >= 3, f"GMT 行少于3列: {parts}"
            # 第一列为名称（非空）
            assert parts[0], "GMT 第一列（名称）为空"
            # 第二列为描述（可以为空字符串）
            # 后续列为基因（至少一个）
            assert len(parts) >= 3, f"GMT 行无基因: {line}"
            rows.append(parts)
    return rows


# ============================================================
# GMTGenerator 单元测试
# ============================================================

class TestGMTGeneratorGO:
    """测试 GO GMT 生成"""

    def test_generate_go_gmt(self, tmp_path):
        """测试从 GO 数据库产物生成 GMT"""
        org_dir = tmp_path / "hsa"
        _create_mock_go_db(org_dir)

        gen = GMTGenerator(organism_dir=str(org_dir))
        output = gen.generate_go_gmt("hsa")

        assert output.endswith("hsa.GO.gmt.gz")
        assert Path(output).exists()

        rows = _validate_gmt_format(output)
        # 应有 2 个 GO term
        assert len(rows) == 2
        # GO:0005576 应有 GENE_A, GENE_B
        go5576 = [r for r in rows if r[0] == "GO:0005576"][0]
        assert "GENE_A" in go5576[2:]
        assert "GENE_B" in go5576[2:]
        assert "GENE_C" not in go5576[2:]
        # 描述列
        assert "extracellular_region" in go5576[1]

    def test_go_gmt_missing_files(self, tmp_path):
        """测试 GO 数据文件缺失时抛出异常"""
        org_dir = tmp_path / "hsa"
        org_dir.mkdir(parents=True)

        gen = GMTGenerator(organism_dir=str(org_dir))
        with pytest.raises(FileNotFoundError, match="GO"):
            gen.generate_go_gmt("hsa")


class TestGMTGeneratorKEGG:
    """测试 KEGG GMT 生成"""

    def test_generate_kegg_gmt(self, tmp_path):
        """测试从 KEGG 数据库产物生成 GMT"""
        org_dir = tmp_path / "hsa"
        _create_mock_kegg_db(org_dir)

        gen = GMTGenerator(organism_dir=str(org_dir))
        output = gen.generate_kegg_gmt("hsa")

        assert output.endswith("hsa.KEGG.gmt.gz")
        assert Path(output).exists()

        rows = _validate_gmt_format(output)
        assert len(rows) == 2
        # hsa04110 应有 GENE_A, GENE_B
        pathway = [r for r in rows if r[0] == "hsa04110"][0]
        assert "GENE_A" in pathway[2:]
        assert "GENE_B" in pathway[2:]
        # 描述列包含分类信息
        assert "Cell_Cycle" in pathway[1]


class TestGMTGeneratorReactome:
    """测试 Reactome GMT 生成"""

    def test_generate_reactome_gmt(self, tmp_path):
        """测试从 Reactome 数据库产物生成 GMT"""
        org_dir = tmp_path / "hsa"
        _create_mock_reactome_db(org_dir)

        gen = GMTGenerator(organism_dir=str(org_dir))
        output = gen.generate_reactome_gmt("hsa")

        assert output.endswith("hsa.Reactome.gmt.gz")
        assert Path(output).exists()

        rows = _validate_gmt_format(output)
        assert len(rows) == 2
        pathway = [r for r in rows if r[0] == "R-HSA-12345"][0]
        assert "GENE_A" in pathway[2:]
        assert "Signal_Transduction" in pathway[1]


class TestGMTGeneratorDO:
    """测试 DO GMT 生成"""

    def test_generate_do_gmt(self, tmp_path):
        """测试从 DO 数据库产物生成 GMT"""
        org_dir = tmp_path / "hsa"
        _create_mock_do_db(org_dir)

        gen = GMTGenerator(organism_dir=str(org_dir))
        output = gen.generate_do_gmt("hsa")

        assert output.endswith("hsa.DO.gmt.gz")
        assert Path(output).exists()

        rows = _validate_gmt_format(output)
        assert len(rows) == 2
        doid = [r for r in rows if r[0] == "DOID:1234"][0]
        assert "GENE_A" in doid[2:]
        assert "Breast_Cancer" in doid[1]


class TestGMTGeneratorDisGeNET:
    """测试 DisGeNET GMT 生成"""

    def test_generate_disgenet_gmt(self, tmp_path):
        """测试从 DisGeNET 数据库产物生成 GMT"""
        org_dir = tmp_path / "hsa"
        _create_mock_disgenet_db(org_dir)

        gen = GMTGenerator(organism_dir=str(org_dir))
        output = gen.generate_disgenet_gmt("hsa")

        assert output.endswith("hsa.DisGeNET.gmt.gz")
        assert Path(output).exists()

        rows = _validate_gmt_format(output)
        assert len(rows) == 2
        cui = [r for r in rows if r[0] == "CUI:0001"][0]
        assert "GENE_A" in cui[2:]
        assert "Breast_Cancer" in cui[1]


class TestGMTGeneratorAll:
    """测试 generate_all_gmt 批量生成"""

    def test_generate_all_gmt_full(self, tmp_path):
        """测试完整数据库时生成所有 GMT"""
        org_dir = tmp_path / "hsa"
        _create_mock_full_db(org_dir)

        gen = GMTGenerator(organism_dir=str(org_dir))
        results = gen.generate_all_gmt("hsa")

        assert len(results) == 5
        assert "GO" in results
        assert "KEGG" in results
        assert "Reactome" in results
        assert "DO" in results
        assert "DisGeNET" in results

        # 验证所有文件存在
        for db_name, path in results.items():
            assert Path(path).exists(), f"{db_name} GMT 文件不存在: {path}"

    def test_generate_all_gmt_partial(self, tmp_path):
        """测试仅有部分数据库时只生成可用的 GMT"""
        org_dir = tmp_path / "hsa"
        _create_mock_go_db(org_dir)
        _create_mock_kegg_db(org_dir)
        # 不创建 Reactome/DO/DisGeNET

        gen = GMTGenerator(organism_dir=str(org_dir))
        results = gen.generate_all_gmt("hsa")

        assert len(results) == 2
        assert "GO" in results
        assert "KEGG" in results
        assert "Reactome" not in results
        assert "DO" not in results
        assert "DisGeNET" not in results

    def test_generate_all_gmt_non_human(self, tmp_path):
        """测试非人类物种不生成 DO/DisGeNET"""
        org_dir = tmp_path / "mmu"
        _create_mock_go_db(org_dir, "mmu")
        _create_mock_kegg_db(org_dir, "mmu")
        _create_mock_reactome_db(org_dir, "mmu")

        gen = GMTGenerator(organism_dir=str(org_dir))
        results = gen.generate_all_gmt("mmu")

        assert len(results) == 3
        assert "GO" in results
        assert "KEGG" in results
        assert "Reactome" in results
        assert "DO" not in results
        assert "DisGeNET" not in results


class TestGMTGeneratorEdgeCases:
    """测试边界情况"""

    def test_empty_matrix(self, tmp_path):
        """测试空矩阵数据（仅有表头无数据行）"""
        org_dir = tmp_path / "hsa"
        org_dir.mkdir(parents=True)

        tab_path = org_dir / "hsa.GO2gene.tab.gz"
        with gzip.open(tab_path, 'wt', encoding='utf-8') as f:
            f.write("Gene\tGO:0005576\n")

        disc_path = org_dir / "GO2disc.gz"
        with gzip.open(disc_path, 'wt', encoding='utf-8') as f:
            f.write("GO:0005576\tcellular_component:extracellular_region\n")

        gen = GMTGenerator(organism_dir=str(org_dir))
        output = gen.generate_go_gmt("hsa")

        # 应生成文件，但内容为空（无基因的 term 被跳过）
        assert Path(output).exists()
        with gzip.open(output, 'rt', encoding='utf-8') as f:
            content = f.read()
            assert content == ""

    def test_no_description(self, tmp_path):
        """测试描述文件中缺少某个 term 的描述"""
        org_dir = tmp_path / "hsa"
        org_dir.mkdir(parents=True)

        tab_path = org_dir / "hsa.GO2gene.tab.gz"
        with gzip.open(tab_path, 'wt', encoding='utf-8') as f:
            f.write("Gene\tGO:0005576\tGO:9999999\n")
            f.write("GENE_A\t1\t1\n")

        disc_path = org_dir / "GO2disc.gz"
        with gzip.open(disc_path, 'wt', encoding='utf-8') as f:
            # GO:9999999 没有描述
            f.write("GO:0005576\tcellular_component:extracellular_region\n")

        gen = GMTGenerator(organism_dir=str(org_dir))
        output = gen.generate_go_gmt("hsa")

        rows = _validate_gmt_format(output)
        assert len(rows) == 2
        # GO:9999999 的描述应为空字符串
        missing_desc = [r for r in rows if r[0] == "GO:9999999"][0]
        assert missing_desc[1] == ""

    def test_gmt_compression_readable(self, tmp_path):
        """测试生成的 .gmt.gz 压缩文件可正确读取"""
        org_dir = tmp_path / "hsa"
        _create_mock_go_db(org_dir)

        gen = GMTGenerator(organism_dir=str(org_dir))
        output = gen.generate_go_gmt("hsa")

        # 使用 gzip 标准库读取
        with gzip.open(output, 'rt', encoding='utf-8') as f:
            lines = f.readlines()
            assert len(lines) == 2
            for line in lines:
                parts = line.strip().split('\t')
                assert len(parts) >= 3


# ============================================================
# DatabaseBuilder 集成测试
# ============================================================

class TestDatabaseBuilderGMTIntegration:
    """测试 DatabaseBuilder 与 GMTGenerator 的集成"""

    def test_generate_gmt_files_with_output_dir(self, tmp_path):
        """测试指定 output_dir 时生成 GMT"""
        org_dir = tmp_path / "hsa"
        _create_mock_full_db(org_dir)

        builder = DatabaseBuilder(root_dir=str(tmp_path / "database"))
        results = builder.generate_gmt_files("hsa", output_dir=str(org_dir))

        assert len(results) == 5
        assert "GO" in results

    def test_generate_gmt_files_auto_detect(self, tmp_path):
        """测试自动检测最新物种数据库目录"""
        root = tmp_path / "database"
        organism_dir = root / "organism" / "v20260101" / "hsa"
        _create_mock_full_db(organism_dir)

        builder = DatabaseBuilder(root_dir=str(root))
        results = builder.generate_gmt_files("hsa")

        assert len(results) == 5

    def test_generate_gmt_files_no_data(self, tmp_path):
        """测试无数据时返回空字典"""
        root = tmp_path / "database"
        root.mkdir(parents=True)
        (root / "organism").mkdir()

        builder = DatabaseBuilder(root_dir=str(root))
        results = builder.generate_gmt_files("hsa")

        assert results == {}

    def test_build_go_generates_gmt(self, tmp_path):
        """测试 build_go 完成后自动生成 GMT"""
        root = Path(tmp_path) / "database"

        # 创建 GO 基础数据
        go_dir = root / "basic" / "go" / "GO20250101"
        go_dir.mkdir(parents=True)

        with gzip.open(go_dir / "gene2go.gz", 'wt') as f:
            f.write("9606\t1\tGO:0005576\t\t\textracellular region\t\t\tcellular_component\n")
            f.write("9606\t2\tGO:0005576\t\t\textracellular region\t\t\tcellular_component\n")

        with gzip.open(go_dir / "gene_info.gz", 'wt') as f:
            f.write("9606\t1\tGENE_A\t\t\t\t\t\t\t\t\n")
            f.write("9606\t2\tGENE_B\t\t\t\t\t\t\t\t\n")

        with open(go_dir / "go-basic.obo", 'w') as f:
            f.write("format-version: 1.2\n")
            f.write("[Term]\nid: GO:0005576\nname: extracellular region\nnamespace: cellular_component\n")

        builder = DatabaseBuilder(root_dir=str(root))
        outdir = builder.build_go(species="hsa", taxid=9606)

        out_path = Path(outdir)
        # 验证原始数据库产物
        assert (out_path / "hsa.GO2gene.tab.gz").exists()
        assert (out_path / "GO2disc.gz").exists()
        # 验证 GMT 文件自动生成
        assert (out_path / "hsa.GO.gmt.gz").exists()


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
