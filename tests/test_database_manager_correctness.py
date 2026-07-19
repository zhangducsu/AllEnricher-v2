import gzip

import pytest

from allenricher.core.config import Config
from allenricher.database.custom_builder import CustomDatabaseBuilder
from allenricher.database.manager import DatabaseManager


def _write_matrix(path, term, genes):
    path.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(path, "wt", encoding="utf-8") as handle:
        handle.write(f"Gene\t{term}\n")
        for gene in genes:
            handle.write(f"{gene}\t1\n")


def _write_text_gzip(path, text):
    path.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(path, "wt", encoding="utf-8") as handle:
        handle.write(text)


def test_each_database_uses_latest_version_that_contains_it(tmp_path):
    root = tmp_path / "database"
    go_dir = root / "organism" / "v20260101" / "hsa"
    kegg_dir = root / "organism" / "v20260201" / "hsa"
    unrelated_dir = root / "organism" / "v20260301" / "hsa"

    _write_matrix(go_dir / "hsa.GO2gene.tab.gz", "GO:1", ["G1", "G2"])
    _write_text_gzip(go_dir / "GO2disc.gz", "GO:1\tbiological_process:test\n")
    _write_matrix(kegg_dir / "hsa.kegg2gene.tab.gz", "hsa00010", ["G2", "G3"])
    _write_text_gzip(kegg_dir / "hsa.kegg2disc.gz", "hsa00010\tTest_pathway\n")
    _write_matrix(
        unrelated_dir / "hsa.WikiPathways2gene.tab.gz", "WP1", ["G4"]
    )

    manager = DatabaseManager(str(root), "hsa")
    manager.load_databases(["GO", "KEGG"])

    assert manager.database_versions == {
        "GO": "v20260101",
        "KEGG": "v20260201",
    }
    assert manager.active_version == "mixed"
    assert manager.databases["GO"]["GO:1"]["genes"] == ["G1", "G2"]
    assert manager.databases["KEGG"]["hsa00010"]["genes"] == ["G2", "G3"]


def test_genome_background_finds_gene_info_under_database_basic(tmp_path):
    root = tmp_path / "database"
    gene_info = root / "basic" / "go" / "GO20260101" / "gene_info.gz"
    _write_text_gzip(
        gene_info,
        "#tax_id\tGeneID\tSymbol\tLocusTag\tSynonyms\tdbXrefs\tchromosome\tmap_location\tdescription\ttype_of_gene\n"
        "9606\t1\tGENE1\t-\t-\t-\t1\t-\tGene 1\tprotein-coding\n"
        "10090\t2\tMouseGene\t-\t-\t-\t1\t-\tMouse gene\tprotein-coding\n",
    )

    manager = DatabaseManager(str(root), "hsa")

    assert manager.get_genome_genes(taxid=9606) == {"GENE1"}


def test_tf_loader_ignores_newer_unrelated_version(tmp_path):
    root = tmp_path / "database"
    tf_dir = root / "organism" / "v20260101" / "hsa"
    newer = root / "organism" / "v20260201" / "hsa"
    _write_matrix(tf_dir / "hsa.ChEA3_2gene.tab.gz", "ENCODE|TF1", ["G1"])
    _write_text_gzip(
        tf_dir / "hsa.ChEA3_2disc.gz",
        "Term_ID\tTerm_Name\tTF\tLibrary\tContext\tEvidence_Type\tTarget_Set_Size\tSource\n"
        "ENCODE|TF1\tTF1 [ENCODE]\tTF1\tENCODE\tENCODE\tChIP-seq\t1\tChEA3\n",
    )
    _write_matrix(newer / "hsa.GO2gene.tab.gz", "GO:1", ["G1"])

    loaded = DatabaseManager(str(root), "hsa").load_chea3()

    assert loaded is not None
    assert loaded["gene2tf"].columns.tolist() == ["Gene", "ENCODE|TF1"]


@pytest.mark.parametrize("database", ["ChEA3", "hTFtarget"])
def test_human_only_tf_databases_reject_nonhuman_species(tmp_path, database):
    manager = DatabaseManager(str(tmp_path), "mmu")
    with pytest.raises(ValueError, match="does not support species"):
        manager.load_database(database)


def test_human_animaltfdb_does_not_alias_htftarget(tmp_path):
    with pytest.raises(ValueError, match="hTFtarget"):
        DatabaseManager(str(tmp_path), "hsa").load_animaltfdb()


def test_standard_analysis_prefers_sparse_gmt_with_matrix_fallback(tmp_path):
    root = tmp_path / "database"
    db_dir = root / "organism" / "v20260101" / "hsa"
    _write_matrix(db_dir / "hsa.GO2gene.tab.gz", "GO:1", ["MATRIX_ONLY"])
    _write_text_gzip(db_dir / "GO2disc.gz", "GO:1\tbiological_process:test\n")
    _write_text_gzip(
        db_dir / "hsa.GO.gmt.gz",
        "GO:1\tbiological_process:test\tGMT_A\tGMT_B\n",
    )

    manager = DatabaseManager(str(root), "hsa")
    manager.load_database("GO")

    assert manager.databases["GO"]["GO:1"]["genes"] == ["GMT_A", "GMT_B"]


def test_custom_database_name_builds_and_loads_end_to_end(tmp_path):
    annotation = tmp_path / "annotation.tsv"
    annotation.write_text(
        "gene\tterm_id\tterm_name\n"
        "G1\tTERM_A\tPathway A\n"
        "G2\tTERM_A\tPathway A\n"
        "G3\tTERM_B\tPathway B\n",
        encoding="utf-8",
    )
    database_root = tmp_path / "database"
    CustomDatabaseBuilder(str(database_root)).build_from_annotation(
        str(annotation), "hsa", 9606, "MyPathways"
    )

    assert Config(species="hsa", databases=["MyPathways"]).validate() == []
    manager = DatabaseManager(str(database_root), "hsa")
    manager.load_databases(["MyPathways"])

    assert set(manager.databases["MyPathways"]) == {"TERM_A", "TERM_B"}
    assert manager.databases["MyPathways"]["TERM_A"]["genes"] == ["G1", "G2"]


def test_custom_database_preserves_hierarchy_for_plotting(tmp_path):
    annotation = tmp_path / "annotation.tsv"
    annotation.write_text(
        "gene\tterm_id\tterm_name\thierarchy\n"
        "G1\tTERM_A\tPathway A\tMetabolism|Carbohydrate|Pathway A\n"
        "G2\tTERM_B\tPathway B\tDisease|Cancer|Pathway B\n",
        encoding="utf-8",
    )
    database_root = tmp_path / "database"
    CustomDatabaseBuilder(str(database_root)).build_from_annotation(
        str(annotation), "hsa", 9606, "MyPathways", format_type="four_column"
    )

    manager = DatabaseManager(str(database_root), "hsa")
    manager.load_database("MyPathways")

    assert manager.databases["MyPathways"]["TERM_A"]["hierarchy"] == (
        "Metabolism|Carbohydrate|Pathway A"
    )


def test_sparse_gmt_does_not_reformat_mapped_term_name(tmp_path):
    species_dir = tmp_path / "organism" / "v20260716" / "hsa"
    species_dir.mkdir(parents=True)
    with gzip.open(species_dir / "hsa.GO.gmt.gz", "wt", encoding="utf-8") as handle:
        handle.write("GO:0005515\tunused\tGENE1\tGENE2\n")
    with gzip.open(species_dir / "GO2disc.gz", "wt", encoding="utf-8") as handle:
        handle.write("GO:0005515\tmolecular_function:protein binding\n")

    manager = DatabaseManager(str(tmp_path), "hsa")
    manager.load_database("GO")

    assert manager.databases["GO"]["GO:0005515"]["name"] == (
        "Molecular Function|Protein Binding"
    )


def test_non_tf_terms_without_real_names_are_excluded_but_tf_symbols_are_valid():
    terms = {
        "DOID:1": {"name": "DOID:1", "genes": ["G1"]},
        "DOID:2": {"name": "Named disease", "genes": ["G2"]},
    }

    assert set(DatabaseManager._drop_terms_without_names("DO", terms)) == {"DOID:2"}
    assert set(DatabaseManager._drop_terms_without_names(
        "TRRUST", {"NFKB1": {"name": "NFKB1", "genes": ["G1"]}}
    )) == {"NFKB1"}


def test_custom_database_name_rejects_path_characters():
    errors = Config(species="hsa", databases=["../outside"]).validate()
    assert any("Invalid database name" in error for error in errors)


def test_custom_species_code_is_validated_by_database_manager():
    assert Config(species="mfor", databases=["MFORKEGG"]).validate() == []


def test_custom_species_code_rejects_path_traversal():
    errors = Config(species="../outside", databases=["MFORKEGG"]).validate()
    assert any("Invalid species code" in error for error in errors)
def test_custom_term_names_preserve_user_capitalization():
    manager = DatabaseManager("unused", "hsa")

    assert manager._format_term_name(
        "MY_CUSTOM_DB", "Molecular Function|Single-stranded DNA Endonuclease Activity"
    ) == "Molecular Function|Single-stranded DNA Endonuclease Activity"
