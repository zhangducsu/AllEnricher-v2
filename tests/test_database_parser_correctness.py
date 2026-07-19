from pathlib import Path

import pytest

from allenricher.database.builder import DatabaseBuilder
from allenricher.database.parsers.animaltfdb import AnimalTFDBParser
from allenricher.database.parsers.wikipathways import WikiPathwaysParser


def test_wikipathways_release_gmt_reads_all_gene_columns(tmp_path):
    gmt = tmp_path / "wikipathways.gmt"
    gmt.write_text(
        "Glutathione metabolism%WikiPathways_20260510%WP100%Homo sapiens\t"
        "https://www.wikipathways.org/instance/WP100\t2687\t2678\t2876\n",
        encoding="utf-8",
    )
    gene_sets, descriptions = WikiPathwaysParser.parse_gmt(str(gmt))
    assert gene_sets == {"WP100": {"2687", "2678", "2876"}}
    assert descriptions == {"WP100": "Glutathione metabolism"}


def test_wikipathways_parser_keeps_standard_and_legacy_gmt_compatibility(tmp_path):
    gmt = tmp_path / "mixed.gmt"
    gmt.write_text(
        "WP1\tPathway one\tG1\tG2\n"
        "WP2\tPathway two\tG3/G4\n",
        encoding="utf-8",
    )
    gene_sets, descriptions = WikiPathwaysParser.parse_gmt(str(gmt))
    assert gene_sets["WP1"] == {"G1", "G2"}
    assert gene_sets["WP2"] == {"G3", "G4"}
    assert descriptions["WP1"] == "Pathway one"


def test_wikipathways_numeric_ncbi_ids_are_mapped_to_symbols():
    converted = WikiPathwaysParser.convert_ncbi_to_symbol(
        {"WP100": {"2687", "ncbigene:2678", "UNCHANGED"}},
        {"2687": "GGT1", "2678": "GGT5"},
    )
    assert converted["WP100"] == {"GGT1", "GGT5", "UNCHANGED"}


def test_htftarget_builder_matches_parser_contract(tmp_path, monkeypatch):
    source = tmp_path / "basic" / "htftarget" / "tf-target-information.txt"
    source.parent.mkdir(parents=True)
    source.write_text(
        "TF\ttarget\ttissue\nTF1\tGENE1\tliver\nTF1\tGENE2\tlung\n",
        encoding="utf-8",
    )
    builder = DatabaseBuilder(str(tmp_path))
    monkeypatch.setattr(builder, "_get_gene_info_path", lambda: None)

    outdir = Path(builder.build_htftarget("hsa", 9606))

    assert (outdir / "hsa.hTF_2gene.tab.gz").is_file()
    assert (outdir / "hsa.hTF_2disc.gz").is_file()


@pytest.mark.parametrize(
    ("method", "database"),
    [("build_chea3", "ChEA3"), ("build_htftarget", "hTFtarget")],
)
def test_human_only_tf_builders_reject_nonhuman_species(tmp_path, method, database):
    with pytest.raises(ValueError, match=database):
        getattr(DatabaseBuilder(str(tmp_path)), method)("mmu", 10090)


def test_animaltfdb_builder_forwards_gene_info_for_official_external_ids(tmp_path, monkeypatch):
    cache = tmp_path / "basic" / "animaltfdb" / "AnimalTFDBv4.0"
    cache.mkdir(parents=True)
    (cache / "Drosophila_melanogaster_TF").write_text("Symbol\nTF1\n", encoding="utf-8")
    (cache / "Drosophila_melanogaster_ortholog_to_human").write_text(
        "Ensembl ID\tOrtholog ID\nFBgn1\tENSG1\n", encoding="utf-8"
    )
    htftarget = tmp_path / "basic" / "htftarget" / "tf-target-information.txt"
    htftarget.parent.mkdir(parents=True)
    htftarget.write_text("TF\ttarget\nTF1\tG1\n", encoding="utf-8")
    gene_info = tmp_path / "gene_info.gz"
    gene_info.write_bytes(b"fixture")
    captured = {}

    def capture_build_database(**kwargs):
        captured.update(kwargs)
        raise RuntimeError("captured")

    monkeypatch.setattr(AnimalTFDBParser, "build_database", capture_build_database)
    builder = DatabaseBuilder(str(tmp_path))
    monkeypatch.setattr(builder, "_load_valid_genes", lambda *_args: {"TF1"})

    with pytest.raises(RuntimeError, match="captured"):
        builder.build_animaltfdb(
            "dme", 7227,
            species_latin="Drosophila_melanogaster",
            gene_info_path=str(gene_info),
        )

    assert captured["gene_info_path"] == str(gene_info)
    assert captured["species_taxid"] == 7227


def test_builder_finds_go_version_from_overridden_basic_directory(tmp_path):
    source_basic = tmp_path / "frozen_sources"
    (source_basic / "go" / "GO20260715").mkdir(parents=True)
    builder = DatabaseBuilder(str(tmp_path / "snapshot"))
    builder.basic_dir = source_basic

    assert builder._get_go_version() == "GO20260715"
