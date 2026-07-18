from __future__ import annotations

import gzip
from pathlib import Path

import pytest

from allenricher.database.manager import DatabaseManager
from allenricher.database.custom_builder import CustomDatabaseBuilder


@pytest.mark.parametrize(
    ("database", "species", "matrix_name", "disc_name", "disc_content"),
    [
        ("TRRUST", "hsa", "hsa.gene2TF.tab.gz", "hsa.TF2disc.gz",
         "Term_ID\tTerm_Name\tTF\tMode\tTarget_Set_Size\tSource\n"
         "TF_A\tTF_A\tTF_A\tActivation\t2\tTRRUST\n"
         "TF_B\tTF_B\tTF_B\tRepression\t2\tTRRUST\n"),
        ("ChEA3", "hsa", "hsa.ChEA3_2gene.tab.gz", "hsa.ChEA3_2disc.gz",
         "Term_ID\tTerm_Name\tTF\tLibrary\tContext\tEvidence_Type\tTarget_Set_Size\tSource\n"
         "TF_A\tTF_A\tTF_A\tENCODE\tENCODE\tChIP-seq\t2\tChEA3\n"
         "TF_B\tTF_B\tTF_B\tReMap\tReMap\tChIP-seq\t2\tChEA3\n"),
        ("AnimalTFDB", "dme", "dme.AnimalTFDB_2gene.tab.gz", "dme.AnimalTFDB_mapped_2disc.gz",
         "Term_ID\tTerm_Name\tTF\tLibrary\tContext\tEvidence_Type\tInference_Type\tHuman_TF\tFamily\tTarget_Set_Size\tSource\n"
         "TF_A\tTF_A\tTF_A\thTFtarget\tliver\tTF-target\tortholog_inferred\tHUMAN_A\tHomeobox\t2\tAnimalTFDB_ortholog_inferred_hTFtarget\n"
         "TF_B\tTF_B\tTF_B\thTFtarget\tlung\tTF-target\tortholog_inferred\tHUMAN_B\tbHLH\t2\tAnimalTFDB_ortholog_inferred_hTFtarget\n"),
        ("hTFtarget", "hsa", "hsa.hTF_2gene.tab.gz", "hsa.hTF_2disc.gz",
         "Term_ID\tTerm_Name\tTF\tLibrary\tContext\tEvidence_Type\tInference_Type\tTarget_Set_Size\tSource\n"
         "TF_A\tTF_A\tTF_A\thTFtarget\tliver\tTF-target\tdirect\t2\thTFtarget\n"
         "TF_B\tTF_B\tTF_B\thTFtarget\tlung\tTF-target\tdirect\t2\thTFtarget\n"),
    ],
)
def test_generic_analyze_loader_uses_real_tf_builder_filenames(
    tmp_path, database, species, matrix_name, disc_name, disc_content
):
    species_dir = tmp_path / "organism" / "v20260715" / species
    species_dir.mkdir(parents=True)
    with gzip.open(species_dir / matrix_name, "wt", encoding="utf-8") as handle:
        handle.write("Gene\tTF_A\tTF_B\n")
        handle.write("GENE1\t1\t0\n")
        handle.write("GENE2\t1\t1\n")
        handle.write("GENE3\t0\t1\n")
    with gzip.open(species_dir / disc_name, "wt", encoding="utf-8") as handle:
        handle.write(disc_content)

    manager = DatabaseManager(str(tmp_path), species)
    manager.load_database(database)

    terms = manager.get_all_term_data()[database]
    assert set(terms["TF_A"]["genes"]) == {"GENE1", "GENE2"}
    assert terms["TF_A"]["name"] == f"TF_A targets [{database}]"
    assert set(terms["TF_B"]["genes"]) == {"GENE2", "GENE3"}


def test_tf_sparse_gmt_keeps_term_name_from_metadata(tmp_path):
    species_dir = tmp_path / "organism" / "v20260716" / "hsa"
    species_dir.mkdir(parents=True)
    with gzip.open(species_dir / "hsa.ChEA3.gmt.gz", "wt", encoding="utf-8") as handle:
        handle.write("ENCODE|TF1\tENCODE|TF1\tGENE1\tGENE2\n")
    with gzip.open(species_dir / "hsa.ChEA3_2disc.gz", "wt", encoding="utf-8") as handle:
        handle.write(
            "Term_ID\tTerm_Name\tTF\tLibrary\tContext\tEvidence_Type\tTarget_Set_Size\tSource\n"
            "ENCODE|TF1\tTF1 [ENCODE]\tTF1\tENCODE\tENCODE\tChIP-seq\t2\tChEA3\n"
        )

    manager = DatabaseManager(str(tmp_path), "hsa")
    manager.load_database("ChEA3")

    assert manager.databases["ChEA3"]["ENCODE|TF1"]["name"] == "TF1 [ENCODE]"


def test_custom_gmt_build_stays_sparse_and_loadable(tmp_path):
    source = tmp_path / "public_go.gmt"
    source.write_text(
        "GO:1\tBiological Process|Signal\tGENE1\tGENE2\n"
        "GO:2\tMolecular Function|Binding\tGENE2\tGENE3\n",
        encoding="utf-8",
    )

    outdir = Path(
        CustomDatabaseBuilder(str(tmp_path)).build_from_annotation(
            annotation_file=str(source),
            species="hsa",
            taxid=9606,
            db_name="PUBLIC_GO_CUSTOM",
        )
    )

    assert (outdir / "hsa.PUBLIC_GO_CUSTOM.gmt.gz").is_file()
    assert (outdir / "PUBLIC_GO_CUSTOM2disc.gz").is_file()
    assert not (outdir / "hsa.PUBLIC_GO_CUSTOM2gene.tab.gz").exists()
    manager = DatabaseManager(str(tmp_path), "hsa")
    manager.load_database("PUBLIC_GO_CUSTOM")
    assert set(manager.get_all_term_data()["PUBLIC_GO_CUSTOM"]["GO:1"]["genes"]) == {"GENE1", "GENE2"}
