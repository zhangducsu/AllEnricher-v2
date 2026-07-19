"""
Database Build Module Module Testing (v2 New architecture)

New architecture:
download → database/basic/{type}/{ver}/ (Common data for all species)
build    → database/organism/v{date}/{species}/(Formatized data for specified species)
"""

import pytest
import gzip
import json
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
from allenricher.database.kegg_fetcher import KEGGFetcher
from allenricher.database.species_registry import SpeciesRegistry


# ============================================================
# Test Helpers
# ============================================================

def _create_mock_go_basic(root: Path, version: str = "GO20250101"):
    """Create simulated GO basic data catalogue for all species"""
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
    """Create a simulated Reactome Basic Data Directory for All Species"""
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

    (re_dir / "ReactomePathways.txt").write_text(
        "R-HSA-10000\tMetabolism\tHomo sapiens\n"
        "R-HSA-12345\tPathway Name 1\tHomo sapiens\n"
        "R-HSA-67890\tPathway Name 2\tHomo sapiens\n",
        encoding="utf-8",
    )
    (re_dir / "ReactomePathwaysRelation.txt").write_text(
        "R-HSA-10000\tR-HSA-12345\n"
        "R-HSA-10000\tR-HSA-67890\n",
        encoding="utf-8",
    )
    return re_dir


# ============================================================
# Parser testing (no change)
# ============================================================

class TestGOParser:
    def test_parse_gene2go(self, tmp_path):
        gene2go_path = tmp_path / "gene2go.gz"
        with gzip.open(gene2go_path, 'wt') as f:
            f.write("9606\t1\tGO:0005576\tIEA\tlocated_in\textracellular region\t-\tComponent\n")
            f.write("9606\t2\tGO:0005576\tIEA\t\textracellular region\t-\tComponent\n")
            f.write("9606\t1\tGO:0051301\tEXP\tinvolved_in\tcell division\t1\tProcess\n")
            f.write("9606\t1\tGO:9999999\tEXP\tNOT|involved_in\tnegated term\t1\tProcess\n")
            f.write("10090\t4\tGO:0005576\tIEA\tlocated_in\textracellular region\t-\tComponent\n")

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
            assert 'GO:9999999' not in header

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

    def test_gene_info_is_filtered_by_taxid(self, tmp_path, capsys):
        gene_info_path = tmp_path / "gene_info.gz"
        with gzip.open(gene_info_path, "wt") as f:
            f.write("9606\t1\tHUMAN\n")
            f.write("10090\t2\tMOUSE\n")
        gene2pathway_path = tmp_path / "gene2pathway.txt"
        gene2pathway_path.write_text("HUMAN\t1\thsa00010\tPathway\n", encoding="utf-8")

        KEGGParser.build_database(
            species="hsa",
            taxid=9606,
            gene_info_path=str(gene_info_path),
            gene2pathway_path=str(gene2pathway_path),
            outdir=str(tmp_path),
        )

        assert "Valid genes for the requested TaxID: 1" in capsys.readouterr().out


class TestKEGGFetcher:
    def test_ncbi_symbol_mapping_is_filtered_by_taxid(self, tmp_path):
        gene_info_path = tmp_path / "gene_info.gz"
        with gzip.open(gene_info_path, "wt") as f:
            f.write("9606\t1\tHUMAN\n")
            f.write("10090\t2\tMOUSE\n")

        mapping = KEGGFetcher(str(tmp_path))._ncbi_id_to_symbol(
            str(gene_info_path), taxid=9606
        )

        assert mapping == {"1": "HUMAN"}

    def test_api_get_retries_interrupted_requests_response(self, tmp_path, monkeypatch):
        import requests

        attempts = []

        class Response:
            text = "path:bta00010\tbta:1\n"

            @staticmethod
            def raise_for_status():
                return None

        def fake_get(*_args, **_kwargs):
            attempts.append(1)
            if len(attempts) < 3:
                raise requests.exceptions.ChunkedEncodingError("premature")
            return Response()

        monkeypatch.setattr(requests, "get", fake_get)
        monkeypatch.setattr("allenricher.database.kegg_fetcher.time.sleep", lambda *_args: None)

        assert KEGGFetcher(str(tmp_path))._api_get("link/bta/pathway").startswith("path:bta")
        assert len(attempts) == 3

    def test_pathway_class_request_uses_species_qualified_kegg_ids(self, tmp_path, monkeypatch):
        fetcher = KEGGFetcher(str(tmp_path))
        endpoints = []

        def fake_api_get(endpoint):
            endpoints.append(endpoint)
            return (
                "ENTRY       sce01100                    Pathway\n"
                "CLASS       Metabolism; Global and overview maps\n"
                "///\n"
                "ENTRY       sce00010                    Pathway\n"
                "CLASS       Metabolism; Carbohydrate metabolism\n"
                "///\n"
            )

        monkeypatch.setattr(fetcher, "_api_get", fake_api_get)
        monkeypatch.setattr("allenricher.database.kegg_fetcher.time.sleep", lambda *_args: None)

        categories = fetcher._get_brite_categories(
            "sce", [("01100", "Metabolic pathways"), ("00010", "Glycolysis")]
        )

        assert endpoints == ["get/sce01100+sce00010"]
        assert categories["01100"] == ("Metabolism", "Global_and_overview_maps")

    def test_cached_pathway_ids_are_normalized_to_bare_ids(self, tmp_path):
        (tmp_path / "sce_pathways.txt").write_text(
            "sce01100\tMetabolic pathways\n00010\tGlycolysis\n", encoding="utf-8"
        )

        pathways = KEGGFetcher(str(tmp_path))._list_pathways("sce")

        assert pathways == [("01100", "Metabolic pathways"), ("00010", "Glycolysis")]

    def test_cached_pathway_classes_ignore_stale_species_prefix(self, tmp_path):
        (tmp_path / "dme_pathway_classes.txt").write_text(
            "hsa00010\tMetabolism\tCarbohydrate_Metabolism\n"
            "00010\tUncategorized\tUncategorized\n",
            encoding="utf-8",
        )

        categories = KEGGFetcher(str(tmp_path))._get_brite_categories(
            "dme", [("00010", "Glycolysis")]
        )

        assert categories == {"00010": ("Metabolism", "Carbohydrate_Metabolism")}
        assert KEGGFetcher._clean_pathway_name(
            "Glycolysis - Drosophila melanogaster (fruit fly)"
        ) == "Glycolysis"


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

    def test_parse_ncbi2reactome_preserves_name_and_hierarchy(self, tmp_path):
        basic = _create_mock_reactome_basic(tmp_path)
        outdir = tmp_path / "output"

        ReactomeParser.parse_ncbi2reactome(
            str(basic / "NCBI2Reactome_All_Levels.txt.gz"),
            str(basic / "gene_info.gz"),
            9606,
            "hsa",
            str(outdir),
            str(basic / "ReactomePathways.txt"),
            str(basic / "ReactomePathwaysRelation.txt"),
        )

        with gzip.open(outdir / "hsa.Reactome2disc.gz", "rt", encoding="utf-8") as handle:
            rows = {parts[0]: parts[1:] for parts in map(lambda line: line.rstrip().split("\t"), handle)}
        assert rows["R-HSA-12345"] == ["Pathway Name 1", "Metabolism|Pathway Name 1"]


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

    def test_parse_disease_files_uses_official_name_hierarchy_and_obsolete_flag(self, tmp_path):
        disease_file = tmp_path / "human_disease_knowledge.tsv"
        disease_file.write_text(
            "col0\tGENE_A\tDOID:2\tLegacy Child\tcol4\n"
            "col0\tGENE_A\tDOID:3\tLegacy Obsolete\tcol4\n",
            encoding="utf-8",
        )
        gene_info_path = tmp_path / "gene_info.gz"
        with gzip.open(gene_info_path, "wt", encoding="utf-8") as handle:
            handle.write("9606\t1\tGENE_A\t\t\t\t\t\t\t\t\n")
        ontology_path = tmp_path / "doid.obo"
        ontology_path.write_text(
            "[Term]\nid: DOID:1\nname: disease\n\n"
            "[Term]\nid: DOID:2\nname: child disease\nis_a: DOID:1 ! disease\n\n"
            "[Term]\nid: DOID:3\nname: obsolete disease\nis_obsolete: true\n",
            encoding="utf-8",
        )

        DOParser.parse_disease_files(
            [str(disease_file)],
            str(gene_info_path),
            9606,
            str(tmp_path / "output"),
            str(ontology_path),
        )

        with gzip.open(tmp_path / "output/hsa.DO2disc.gz", "rt", encoding="utf-8") as handle:
            rows = [line.rstrip().split("\t") for line in handle]
        assert rows == [["DOID:2", "child disease", "disease|child disease"]]


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
# New architecture test: Download  build complete process
# ============================================================

class TestDatabaseBuilderNew:
    """Test new DatasBuilder (v2 architecture)"""

    def test_build_go_from_basic(self, tmp_path):
        """From database/basic/Build GO databases"""
        root = Path(tmp_path) / "database"
        _create_mock_go_basic(root)

        builder = DatabaseBuilder(root_dir=str(root))
        outdir = builder.build_go(species="hsa", taxid=9606)

        # Validate the output position: database/organism/v{date}/hsa/
        assert outdir.startswith(str(root / "organism"))
        assert "hsa" in outdir

        out_path = Path(outdir)
        assert (out_path / "hsa.GO2gene.tab.gz").exists()
        assert (out_path / "GO2disc.gz").exists()
        assert (out_path / "hsa.gene2go.txt").exists()

    def test_build_reactome_from_basic(self, tmp_path):
        """From database/basic/Build Reactome Database"""
        root = Path(tmp_path) / "database"
        _create_mock_reactome_basic(root)

        builder = DatabaseBuilder(root_dir=str(root))
        outdir = builder.build_reactome(species="hsa", taxid=9606)

        out_path = Path(outdir)
        assert (out_path / "hsa.Reactome2gene.tab.gz").exists()
        assert (out_path / "hsa.Reactome2disc.gz").exists()

    def test_build_do_uses_ontology_metadata(self, tmp_path):
        root = Path(tmp_path) / "database"
        _create_mock_go_basic(root)
        do_dir = root / "basic/do"
        do_dir.mkdir(parents=True)
        (do_dir / "human_disease_knowledge_filtered.tsv").write_text(
            "col0\tGENE_A\tDOID:2\tLegacy Child\tcol4\n",
            encoding="utf-8",
        )
        (do_dir / "doid.obo").write_text(
            "[Term]\nid: DOID:1\nname: disease\n\n"
            "[Term]\nid: DOID:2\nname: child disease\nis_a: DOID:1 ! disease\n",
            encoding="utf-8",
        )

        outdir = DatabaseBuilder(root_dir=str(root)).build_do(9606, "GO20250101")

        with gzip.open(Path(outdir) / "hsa.DO2disc.gz", "rt", encoding="utf-8") as handle:
            assert handle.read().strip() == "DOID:2\tchild disease\tdisease|child disease"

    def test_build_species_db(self, tmp_path):
        """One-click to build a species database"""
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
        """Error hint when no basic data are available for testing"""
        root = Path(tmp_path) / "database"
        root.mkdir(parents=True)
        (root / "basic").mkdir()

        builder = DatabaseBuilder(root_dir=str(root))
        with pytest.raises(FileNotFoundError):
            builder.build_go(species="hsa", taxid=9606)

    def test_skip_disgenet_for_non_human(self, tmp_path):
        """Test non-human species for automatic skipping DisGeNET"""
        root = Path(tmp_path) / "database"
        _create_mock_go_basic(root)

        builder = DatabaseBuilder(root_dir=str(root))
        # Mouse is not eligible for the human-only DisGeNET build and must be reported as skipped.
        outdir = builder.build_species_db(
            species="mmu", taxid=10090,
            databases=["GO", "DisGeNET"]
        )
        # Should be done properly (DisGeNET skipped)
        assert Path(outdir).exists()

    def test_only_unsupported_database_still_writes_manifest(self, tmp_path):
        builder = DatabaseBuilder(root_dir=str(tmp_path / "database"))

        outdir = Path(builder.build_species_db("mmu", 10090, ["DO"]))

        manifest = __import__("json").loads(
            (outdir / "build_manifest.json").read_text(encoding="utf-8")
        )
        assert manifest["built_databases"] == []
        assert manifest["skipped_databases"] == ["DO"]
        assert "DO" not in manifest["dependencies"]

    def test_build_disgenet_reuses_v1_database(self, tmp_path):
        root = tmp_path / "database"
        legacy = root / "organism" / "v20190612" / "hsa"
        legacy.mkdir(parents=True)
        for name in ("hsa.CUI2gene.tab.gz", "hsa.CUI2disc.gz", "hsa.DisGeNET.gmt.gz"):
            (legacy / name).write_bytes(b"legacy")
        outdir = Path(DatabaseBuilder(root_dir=str(root)).build_disgenet(9606))
        assert (outdir / "hsa.CUI2gene.tab.gz").read_bytes() == b"legacy"
        assert (outdir / "hsa.CUI2disc.gz").read_bytes() == b"legacy"

    def test_build_species_db_propagates_database_failure(self, tmp_path, monkeypatch):
        builder = DatabaseBuilder(root_dir=str(tmp_path / "database"))
        monkeypatch.setattr(builder, "build_go_with_fallback", lambda *args, **kwargs: (_ for _ in ()).throw(RuntimeError("broken")))
        with pytest.raises(RuntimeError, match="GO: broken"):
            builder.build_species_db("hsa", 9606, ["GO"])


class TestDataDownloaderNew:
    """Test the new DataDownloader"""

    def test_init_creates_dirs(self):
        downloader = DataDownloader(root_dir="/tmp/test_downloader")
        assert downloader.basic_dir == Path("/tmp/test_downloader/basic")

    def test_version_listing(self, tmp_path):
        """Test version list function"""
        root = Path(tmp_path) / "database"
        _create_mock_go_basic(root)

        downloader = DataDownloader(root_dir=str(root))
        versions = downloader.list_go_versions()
        assert "GO20250101" in versions
        assert downloader.get_latest_go_version() == "GO20250101"

    def test_no_versions(self, tmp_path):
        """Return empty when no version is tested"""
        downloader = DataDownloader(root_dir=str(tmp_path))
        assert downloader.get_latest_go_version() is None
        assert downloader.list_go_versions() == []

    def test_go_download_records_aggregate_version(self, tmp_path, monkeypatch):
        downloader = DataDownloader(root_dir=str(tmp_path))

        def fake_download(_mirrors, filename, destination, **_kwargs):
            destination.parent.mkdir(parents=True, exist_ok=True)
            if filename.endswith(".gz"):
                with gzip.open(destination, "wt", encoding="utf-8") as stream:
                    stream.write("9606\t1\n")
            else:
                destination.write_text("format-version: 1.2\n", encoding="utf-8")
            return destination

        monkeypatch.setattr(downloader.manager, "download_with_mirror_fallback", fake_download)
        monkeypatch.setattr(
            "allenricher.database.version.RemoteVersionChecker.check_head",
            lambda *_args, **_kwargs: {"last_modified": "today"},
        )
        monkeypatch.setattr(
            "allenricher.database.version.RemoteVersionChecker.check_go_obo_version",
            lambda *_args, **_kwargs: {"remote_version": "release-current", "last_modified": "today"},
        )

        downloader.download_go_basic(version="GO20260713")
        manifest = json.loads((tmp_path / "versions.json").read_text(encoding="utf-8"))

        assert manifest["versions"]["go"]["local_version"] == "GO20260713"
        assert manifest["versions"]["go"]["remote_version"] == "release-current"

    def test_disgenet_download_reuses_v1_database(self, tmp_path):
        downloader = DataDownloader(root_dir=str(tmp_path))
        legacy = tmp_path / "organism" / "v20190612" / "hsa" / "hsa.DisGeNET.gmt.gz"
        legacy.parent.mkdir(parents=True)
        legacy.write_bytes(b"legacy")
        assert downloader.download_disgenet() == str(legacy)

    def test_disgenet_download_records_human_in_unified_registry(self, tmp_path):
        downloader = DataDownloader(root_dir=str(tmp_path))
        legacy = tmp_path / "organism" / "v20190612" / "hsa" / "hsa.DisGeNET.gmt.gz"
        legacy.parent.mkdir(parents=True)
        legacy.write_bytes(b"legacy")

        downloader.download_all(["disgenet"])
        registry = SpeciesRegistry(tmp_path / "basic" / "supported_species.tsv")
        registry.load()

        assert registry.get_summary()["disgenet"]["count"] == 1
        assert registry.entries[9606].has_disgenet is True

    def test_disgenet_download_fails_without_v1_database(self, tmp_path):
        downloader = DataDownloader(root_dir=str(tmp_path))
        with pytest.raises(RuntimeError, match="Current DisGeNET releases require authorization"):
            downloader.download_all(["disgenet"])

    def test_taxonomy_download_is_dispatched_without_registry_rebuild(self, tmp_path, monkeypatch):
        downloader = DataDownloader(root_dir=str(tmp_path))
        taxonomy_file = tmp_path / "basic" / "taxonomy" / "taxid_to_name.tsv"

        def fake_download(_output_dir):
            taxonomy_file.parent.mkdir(parents=True)
            taxonomy_file.write_text("taxid\tlatin_name\n", encoding="utf-8")
            return taxonomy_file

        monkeypatch.setattr(downloader, "_download_taxonomy_names", fake_download)
        result = downloader.download_all(["taxonomy"])

        assert result == {"taxonomy": str(taxonomy_file.parent)}
        assert not (tmp_path / "basic" / "supported_species.tsv").exists()

    def test_kegg_registry_matches_current_api_format_by_taxonomy(self, tmp_path, monkeypatch):
        downloader = DataDownloader(root_dir=str(tmp_path))
        taxonomy = tmp_path / "basic" / "taxonomy" / "taxid_to_name.tsv"
        taxonomy.parent.mkdir(parents=True)
        taxonomy.write_text(
            "taxid\tlatin_name\n9606\tHomo sapiens\n10090\tMus musculus\n",
            encoding="utf-8",
        )
        monkeypatch.setattr(
            downloader,
            "_api_get_with_retry",
            lambda *args, **kwargs: (
                "T01001\thsa; Homo sapiens (human)\n"
                "T01002\tmmu; Mus musculus (mouse)\n"
            ),
        )

        path = downloader._build_kegg_registry(tmp_path / "basic" / "kegg")
        rows = path.read_text(encoding="utf-8").splitlines()

        assert rows[1].startswith("9606\tHomo sapiens\thsa\t")
        assert rows[2].startswith("10090\tMus musculus\tmmu\t")

    def test_kegg_registry_reuses_cache_on_network_failure(self, tmp_path, monkeypatch):
        downloader = DataDownloader(root_dir=str(tmp_path))
        cached = tmp_path / "basic" / "kegg" / "kegg_species_registry.tsv"
        cached.parent.mkdir(parents=True)
        cached.write_text("taxid\tlatin_name\n9606\tHomo sapiens\n", encoding="utf-8")
        monkeypatch.setattr(
            downloader,
            "_api_get_with_retry",
            lambda *args, **kwargs: (_ for _ in ()).throw(RuntimeError("offline")),
        )

        assert downloader._build_kegg_registry(cached.parent) == cached

    def test_reactome_registry_reads_pathway_id_column(self, tmp_path):
        downloader = DataDownloader(root_dir=str(tmp_path))
        source = tmp_path / "NCBI2Reactome.txt.gz"
        with gzip.open(source, "wt", encoding="utf-8") as stream:
            stream.write("1\tR-HSA-109582\thttps://reactome.org/x\tHemostasis\n")

        path = downloader._build_reactome_registry(source, tmp_path)

        assert "9606\tHSA" in path.read_text(encoding="utf-8")

    def test_auto_kegg_abbreviation_is_not_database_support(self, tmp_path):
        downloader = DataDownloader(root_dir=str(tmp_path))
        go = tmp_path / "go.tsv"
        go.write_text(
            "taxid\tlatin_name\tsource\tgene_count\tterm_count\n"
            "1\troot\tncbi_gene2go\t1\t1\n",
            encoding="utf-8",
        )
        missing = tmp_path / "missing.tsv"
        output = tmp_path / "supported_species.tsv"

        downloader._merge_all_registries(go, missing, missing, missing, missing, output)
        lines = output.read_text(encoding="utf-8").splitlines()
        row = dict(zip(lines[0].split("\t"), lines[1].split("\t")))

        assert row["kegg_code"] == "roo"
        assert row["kegg_code_source"] == "auto"
        assert row["has_kegg"] == "False"

    def test_registry_refresh_preserves_tf_coverage_when_core_registry_changes(self, tmp_path):
        downloader = DataDownloader(root_dir=str(tmp_path))
        downloader.record_database_species(
            "TRRUST", [(9606, "Homo sapiens"), (10090, "Mus musculus")]
        )
        downloader.refresh_supported_species_registry()

        go = tmp_path / "basic" / "go" / "GO1" / "go_species_registry.tsv"
        go.parent.mkdir(parents=True)
        go.write_text(
            "taxid\tlatin_name\tsource\tgene_count\tterm_count\n"
            "9606\tHomo sapiens\tncbi_gene2go\t1\t1\n",
            encoding="utf-8",
        )
        downloader.refresh_supported_species_registry()

        registry = SpeciesRegistry(tmp_path / "basic" / "supported_species.tsv")
        registry.load()
        assert registry.entries[9606].has_go is True
        assert registry.entries[9606].has_trrust is True
        assert registry.entries[10090].has_trrust is True

    def test_wikipathways_prefers_species_name_over_colliding_code(self, tmp_path):
        downloader = DataDownloader(root_dir=str(tmp_path))
        go = tmp_path / "go.tsv"
        go.write_text(
            "taxid\tlatin_name\tsource\tgene_count\tterm_count\n"
            "9796\tEquus caballus\tuniprot_goa\t1\t1\n"
            "4577\tZea mays\tuniprot_goa\t1\t1\n"
            "9615\tCanis lupus familiaris\tuniprot_goa\t1\t1\n",
            encoding="utf-8",
        )
        kegg = tmp_path / "kegg.tsv"
        kegg.write_text(
            "taxid\tlatin_name\tkegg_code\tkegg_code_source\tgene_count\n"
            "9796\tEquus caballus\tecb\tkegg\t1\n"
            "218491\tPectobacterium atrosepticum SCRI1043\teca\tkegg\t1\n"
            "9615\tCanis lupus familiaris\tcfa\tkegg\t1\n",
            encoding="utf-8",
        )
        wiki = tmp_path / "wiki.tsv"
        wiki.write_text(
            "species_latin_name\tspecies_code\tgene_count\tpathway_count\n"
            "Equus caballus\teca\t1\t1\n"
            "Zea mays\t-\t1\t1\n"
            "Canis familiaris\tcfa\t1\t1\n",
            encoding="utf-8",
        )
        names = tmp_path / "basic" / "taxonomy" / "names.dmp"
        names.parent.mkdir(parents=True)
        names.write_text("9615\t|\tCanis familiaris\t|\t\t|\tsynonym\t|\n", encoding="utf-8")
        missing = tmp_path / "missing.tsv"
        output = tmp_path / "supported_species.tsv"

        downloader._merge_all_registries(go, kegg, missing, missing, wiki, output)
        registry = SpeciesRegistry(output)
        registry.load()
        assert registry.entries[9796].has_wikipathways is True
        assert registry.entries[4577].has_wikipathways is True
        assert registry.entries[9615].has_wikipathways is True
        assert registry.entries[218491].has_wikipathways is False


if __name__ == "__main__":
    pytest.main([__file__, "-v"])


def test_species_summary_excludes_human_from_animaltfdb_coverage(tmp_path):
    from allenricher.database.species_registry import SpeciesEntry, SpeciesRegistry

    registry = SpeciesRegistry(tmp_path / "supported_species.tsv")
    registry.add_entry(SpeciesEntry(taxid=9606, latin_name="Homo sapiens", kegg_code="hsa", has_animaltfdb=True))
    registry.add_entry(SpeciesEntry(taxid=7227, latin_name="Drosophila melanogaster", kegg_code="dme", has_animaltfdb=True))

    summary = registry.get_summary()

    assert summary["animaltfdb"]["count"] == 1
