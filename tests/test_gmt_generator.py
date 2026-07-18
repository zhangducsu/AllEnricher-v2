"""Tests for generating standard GMT gene-set files.

Coverage includes GMT columns, compressed inputs, empty databases, and DatabaseBuilder
integration.
"""

import pytest
import gzip
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).parent.parent))

from allenricher.database.gmt_generator import GMTGenerator
from allenricher.database.builder import DatabaseBuilder


# ============================================================
# Test fixtures
# ============================================================

def _create_mock_go_db(org_dir: Path, species: str = "hsa"):
    """Create a synthetic GO database in the standard on-disk format."""
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
    """Create simulated KEGG database product file"""
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
    """Create a simulated Reactome database product file"""
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
    """Create a mock DO database product file"""
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
    """Create a mock DisGeNET database product file"""
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
    """Create simulated complete species database (all types)"""
    _create_mock_go_db(org_dir, species)
    _create_mock_kegg_db(org_dir, species)
    _create_mock_reactome_db(org_dir, species)
    if species == "hsa":
        _create_mock_do_db(org_dir)
        _create_mock_disgenet_db(org_dir)


def _create_mock_wikipathways_db(org_dir: Path, species: str = "hsa"):
    org_dir.mkdir(parents=True, exist_ok=True)
    with gzip.open(org_dir / f"{species}.WikiPathways2gene.tab.gz", "wt", encoding="utf-8") as f:
        f.write("Gene\tWP1\nGENE_A\t1\nGENE_B\t1\n")
    with gzip.open(org_dir / f"{species}.WikiPathways2disc.gz", "wt", encoding="utf-8") as f:
        f.write("WP1\tExample_pathway\n")


# ============================================================
# GMT Format Validation
# ============================================================

def _validate_gmt_format(filepath: str):
    """Validation GMT File Format Correct

    Args:
filepath: .gmt.gz File Path

    Returns:
list: Resarse all row data [[name, desc, gene1, gene2, ...], ...]
    """
    rows = []
    with gzip.open(filepath, 'rt', encoding='utf-8') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split('\t')
            # At least three columns: name, description, at least one gene
            assert len(parts) >= 3, f"GMT rows are less than three columns: {parts}"
            # First named (non-empty)
            assert parts[0], "GMT First Column (Name) is empty"
            # Second one is a description (for empty string)
            # Follow-up as gene (at least one)
            assert len(parts) >= 3, f"GMT Gene Free: {line}"
            rows.append(parts)
    return rows


# ============================================================
# GMTGenerator Unit Test
# ============================================================

class TestGMTGeneratorGO:
    """Test GO GMT generation"""

    def test_generate_go_gmt(self, tmp_path):
        """Test to generate GMT from GO database product"""
        org_dir = tmp_path / "hsa"
        _create_mock_go_db(org_dir)

        gen = GMTGenerator(organism_dir=str(org_dir))
        output = gen.generate_go_gmt("hsa")

        assert output.endswith("hsa.GO.gmt.gz")
        assert Path(output).exists()

        rows = _validate_gmt_format(output)
        # 2 GO term
        assert len(rows) == 2
        extracellular = [row for row in rows if row[0] == "GO:0005576"][0]
        assert set(extracellular[2:]) == {"GENE_A", "GENE_B"}
        assert "GO:0005615" not in extracellular[2:]
        # GO: 0005576 with GENE_A, GENE_B
        go5576 = [r for r in rows if r[0] == "GO:0005576"][0]
        assert "GENE_A" in go5576[2:]
        assert "GENE_B" in go5576[2:]
        assert "GENE_C" not in go5576[2:]
        # Description Row
        assert "extracellular_region" in go5576[1]

    def test_go_gmt_missing_files(self, tmp_path):
        """Test GO data file missing to throw an anomaly"""
        org_dir = tmp_path / "hsa"
        org_dir.mkdir(parents=True)

        gen = GMTGenerator(organism_dir=str(org_dir))
        with pytest.raises(FileNotFoundError, match="GO"):
            gen.generate_go_gmt("hsa")


class TestGMTGeneratorKEGG:
    """Test KEG GMT generation"""

    def test_generate_kegg_gmt(self, tmp_path):
        """Test to generate GMT from KEG database product"""
        org_dir = tmp_path / "hsa"
        _create_mock_kegg_db(org_dir)

        gen = GMTGenerator(organism_dir=str(org_dir))
        output = gen.generate_kegg_gmt("hsa")

        assert output.endswith("hsa.KEGG.gmt.gz")
        assert Path(output).exists()

        rows = _validate_gmt_format(output)
        assert len(rows) == 2
        # Hsa04110 should be GENE_A, GENE_B
        pathway = [r for r in rows if r[0] == "hsa04110"][0]
        assert "GENE_A" in pathway[2:]
        assert "GENE_B" in pathway[2:]
        # Description column contains classified information
        assert "Cell_Cycle" in pathway[1]


class TestGMTGeneratorReactome:
    """Test Reactome GMT Generation"""

    def test_generate_reactome_gmt(self, tmp_path):
        """Test to generate GMT from Reactome database product"""
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
    """Test Do GMT generation"""

    def test_generate_do_gmt(self, tmp_path):
        """Test to generate GMT from DO database product"""
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
    """Test DisGeNET GMT generation"""

    def test_generate_disgenet_gmt(self, tmp_path):
        """Test to generate GMT from DisGeNET database product"""
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
    """Test volume generation of generate_all_gmt"""

    def test_generate_all_gmt_full(self, tmp_path):
        """Generate all GMTs when testing the complete database"""
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

        # Verify that all files exist
        for db_name, path in results.items():
            assert Path(path).exists(), f"{db_name}The GMT file does not exist: {path}"

    def test_generate_all_gmt_partial(self, tmp_path):
        """Only available GMT is generated when only some databases are tested"""
        org_dir = tmp_path / "hsa"
        _create_mock_go_db(org_dir)
        _create_mock_kegg_db(org_dir)
        # Do Not Create Reactome/DO/DisGeNET

        gen = GMTGenerator(organism_dir=str(org_dir))
        results = gen.generate_all_gmt("hsa")

        assert len(results) == 2
        assert "GO" in results
        assert "KEGG" in results
        assert "Reactome" not in results
        assert "DO" not in results
        assert "DisGeNET" not in results

    def test_generate_all_gmt_non_human(self, tmp_path):
        """Test non-human species for NO/DisGeNET"""
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

    def test_generate_all_includes_wikipathways(self, tmp_path):
        org_dir = tmp_path / "hsa"
        _create_mock_wikipathways_db(org_dir)

        results = GMTGenerator(str(org_dir)).generate_all_gmt("hsa")

        assert set(results) == {"WikiPathways"}
        rows = _validate_gmt_format(results["WikiPathways"])
        assert rows == [["WP1", "Example_pathway", "GENE_A", "GENE_B"]]


class TestGMTGeneratorEdgeCases:
    """Test the boundary."""

    def test_empty_matrix(self, tmp_path):
        """Test empty matrix data (no data lines for tophead only)"""
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

        # The file is created but remains empty because the invalid term is skipped.
        assert Path(output).exists()
        with gzip.open(output, 'rt', encoding='utf-8') as f:
            content = f.read()
            assert content == ""

    def test_no_description(self, tmp_path):
        """A description of the term missing in the test description file"""
        org_dir = tmp_path / "hsa"
        org_dir.mkdir(parents=True)

        tab_path = org_dir / "hsa.GO2gene.tab.gz"
        with gzip.open(tab_path, 'wt', encoding='utf-8') as f:
            f.write("Gene\tGO:0005576\tGO:9999999\n")
            f.write("GENE_A\t1\t1\n")

        disc_path = org_dir / "GO2disc.gz"
        with gzip.open(disc_path, 'wt', encoding='utf-8') as f:
            # GO: 99999999 Not described
            f.write("GO:0005576\tcellular_component:extracellular_region\n")

        gen = GMTGenerator(organism_dir=str(org_dir))
        output = gen.generate_go_gmt("hsa")

        rows = _validate_gmt_format(output)
        assert len(rows) == 2
        # GO: 99999999 shall be described as an empty string
        missing_desc = [r for r in rows if r[0] == "GO:9999999"][0]
        assert missing_desc[1] == ""

    def test_gmt_compression_readable(self, tmp_path):
        """Test generated.gmt.gz compression file to read correctly"""
        org_dir = tmp_path / "hsa"
        _create_mock_go_db(org_dir)

        gen = GMTGenerator(organism_dir=str(org_dir))
        output = gen.generate_go_gmt("hsa")

        # Read using gzip standard library
        with gzip.open(output, 'rt', encoding='utf-8') as f:
            lines = f.readlines()
            assert len(lines) == 2
            for line in lines:
                parts = line.strip().split('\t')
                assert len(parts) >= 3


# ============================================================
# DatabaseBuilder Integration Test
# ============================================================

class TestDatabaseBuilderGMTIntegration:
    """Test the integration of the DatabaseBuilder with GMTGenerator"""

    def test_generate_gmt_files_with_output_dir(self, tmp_path):
        """GMT generated when testing output_dir"""
        org_dir = tmp_path / "hsa"
        _create_mock_full_db(org_dir)

        builder = DatabaseBuilder(root_dir=str(tmp_path / "database"))
        results = builder.generate_gmt_files("hsa", output_dir=str(org_dir))

        assert len(results) == 5
        assert "GO" in results

    def test_generate_gmt_files_auto_detect(self, tmp_path):
        """Test the most recent species database catalogue for automatic testing"""
        root = tmp_path / "database"
        organism_dir = root / "organism" / "v20260101" / "hsa"
        _create_mock_full_db(organism_dir)

        builder = DatabaseBuilder(root_dir=str(root))
        results = builder.generate_gmt_files("hsa")

        assert len(results) == 5

    def test_generate_gmt_files_no_data(self, tmp_path):
        """Returns empty dictionary when testing non-data"""
        root = tmp_path / "database"
        root.mkdir(parents=True)
        (root / "organism").mkdir()

        builder = DatabaseBuilder(root_dir=str(root))
        results = builder.generate_gmt_files("hsa")

        assert results == {}

    def test_build_go_generates_gmt(self, tmp_path):
        """Auto-generated GMT after testing build_go"""
        root = Path(tmp_path) / "database"

        # Create GO Basic Data
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
        # Validate original database product
        assert (out_path / "hsa.GO2gene.tab.gz").exists()
        assert (out_path / "GO2disc.gz").exists()
        # Validate GMT file automatically generated
        assert (out_path / "hsa.GO.gmt.gz").exists()


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
