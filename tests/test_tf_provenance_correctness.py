from __future__ import annotations

import gzip

import pandas as pd
import pytest
from scipy.stats import hypergeom
from statsmodels.stats.multitest import multipletests

from allenricher.analysis.tf_enrichment import TFEnrichmentAnalyzer
from allenricher.analysis.tf_meta_analyzer import TFMetaAnalyzer
from allenricher.database.ortholog_mapper import OrthologMapper
from allenricher.database.parsers.chea3 import ChEA3Parser
from allenricher.database.parsers.htftarget import HTFtargetParser
from allenricher.database.parsers.trrust import TRRUSTParser


def _read_gzip_tsv(path):
    return pd.read_csv(path, sep="\t", compression="gzip", dtype=str)


def test_chea3_keeps_first_target_and_separates_source_libraries(tmp_path):
    encode = tmp_path / "ENCODE_tf.gmt"
    remap = tmp_path / "ReMap_tf.gmt"
    encode.write_text("TFX\tGENE_A\tGENE_B\n", encoding="utf-8")
    remap.write_text("TFX\tGENE_C\tGENE_D\n", encoding="utf-8")

    output = tmp_path / "output"
    ChEA3Parser.build_database(
        [str(encode), str(remap)], str(output), "hsa"
    )

    matrix = _read_gzip_tsv(output / "hsa.ChEA3_2gene.tab.gz")
    metadata = _read_gzip_tsv(output / "hsa.ChEA3_2disc.gz")

    assert set(matrix.columns) == {"Gene", "ENCODE|TFX", "ReMap|TFX"}
    assert matrix.loc[matrix["Gene"] == "GENE_A", "ENCODE|TFX"].item() == "1"
    assert matrix.loc[matrix["Gene"] == "GENE_A", "ReMap|TFX"].item() == "0"
    assert matrix.loc[matrix["Gene"] == "GENE_C", "ReMap|TFX"].item() == "1"
    assert metadata.set_index("Term_ID").loc["ENCODE|TFX", "Library"] == "ENCODE"
    assert metadata.set_index("Term_ID").loc["ReMap|TFX", "Library"] == "ReMap"

    with pytest.raises(ValueError, match="separate"):
        ChEA3Parser.merge_libraries(
            {"ENCODE": {"TFX": {"GENE_A"}}}, method="union"
        )


def test_htftarget_keeps_tissue_specific_target_sets(tmp_path):
    source = tmp_path / "tf-target-information.txt"
    source.write_text(
        "TF\tTarget\tTissue\n"
        "TFX\tGENE_A\tLiver\n"
        "TFX\tGENE_B\tLung\n"
        "TFX\tGENE_C\tLiver,Lung\n",
        encoding="utf-8",
    )

    output = tmp_path / "output"
    HTFtargetParser.build_database(str(source), str(output), "hsa")
    matrix = _read_gzip_tsv(output / "hsa.hTF_2gene.tab.gz")
    metadata = _read_gzip_tsv(output / "hsa.hTF_2disc.gz").set_index("Term_ID")

    assert set(matrix.columns) == {"Gene", "TFX|Liver", "TFX|Lung"}
    liver = set(matrix.loc[matrix["TFX|Liver"] == "1", "Gene"])
    lung = set(matrix.loc[matrix["TFX|Lung"] == "1", "Gene"])
    assert liver == {"GENE_A", "GENE_C"}
    assert lung == {"GENE_B", "GENE_C"}
    assert metadata.loc["TFX|Liver", "Context"] == "Liver"
    assert metadata.loc["TFX|Lung", "Context"] == "Lung"


def test_trrust_preserves_edges_and_regulation_filter_changes_gene_set(tmp_path):
    source = tmp_path / "trrust.tsv"
    source.write_text(
        "TFX\tGENE_A\tActivation\t111\n"
        "TFX\tGENE_B\tRepression\t222\n"
        "TFX\tGENE_C\tUnknown\t333\n",
        encoding="utf-8",
    )
    output = tmp_path / "output"
    TRRUSTParser.build_database(str(source), str(output), "hsa")

    edges = _read_gzip_tsv(output / "hsa.TRRUST_edges.tsv.gz")
    matrix = _read_gzip_tsv(output / "hsa.gene2TF.tab.gz")
    metadata = _read_gzip_tsv(output / "hsa.TF2disc.gz")
    assert set(edges["Mode"]) == {"Activation", "Repression", "Unknown"}
    assert set(edges["PMID"]) == {"111", "222", "333"}
    assert set(matrix["Gene"]) == {"GENE_A", "GENE_B", "GENE_C"}
    assert "TFX" not in set(matrix["Gene"])
    assert metadata.loc[metadata["Term_ID"] == "TFX", "Term_Name"].item() == (
        "TFX targets [TRRUST]"
    )

    database = {
        "gene2tf": pd.read_csv(
            output / "hsa.gene2TF.tab.gz", sep="\t", compression="gzip"
        ),
        "tf_info": pd.read_csv(
            output / "hsa.TF2disc.gz", sep="\t", compression="gzip"
        ),
        "edges": edges,
    }
    result = TFEnrichmentAnalyzer(database).ora(
        ["GENE_A"], min_overlap=1, regulation="activation", min_size=1
    )

    assert result.loc[0, "TF"] == "TFX"
    assert result.loc[0, "TF_Targets"] == 1
    assert result.loc[0, "Background_Size"] == 1
    assert result.loc[0, "Overlap_Genes"] == "GENE_A"
    assert result.loc[0, "Mode"] == "activation"

    all_result = TFEnrichmentAnalyzer(database).ora(
        ["GENE_A"], min_overlap=1, regulation=None, min_size=1
    )
    assert all_result.loc[0, "Background_Size"] == 3


def test_tf_ora_bh_uses_all_tested_terms_and_corrects_within_library():
    genes = [f"G{i}" for i in range(1, 9)]
    matrix = pd.DataFrame({
        "Gene": genes,
        "A|TF1": [1, 1, 0, 0, 0, 0, 0, 0],
        "A|TF2": [0, 0, 1, 1, 0, 0, 0, 0],
        "B|TF3": [1, 0, 0, 0, 1, 0, 0, 0],
    })
    metadata = pd.DataFrame([
        {"Term_ID": "A|TF1", "Term_Name": "TF1 [A]", "TF": "TF1", "Library": "A"},
        {"Term_ID": "A|TF2", "Term_Name": "TF2 [A]", "TF": "TF2", "Library": "A"},
        {"Term_ID": "B|TF3", "Term_Name": "TF3 [B]", "TF": "TF3", "Library": "B"},
    ])
    result = TFEnrichmentAnalyzer(
        {"gene2tf": matrix, "tf_info": metadata}
    ).ora(["G1", "G2"], min_overlap=1, min_size=1)
    indexed = result.set_index("Term_ID")

    # A Library background is G1-G4; B Library background is G1/G5, cannot share the full matrix of 8 lines.
    p_a1 = hypergeom.sf(2 - 1, 4, 2, 2)
    p_a2 = hypergeom.sf(0 - 1, 4, 2, 2)
    expected_a1_fdr = multipletests([p_a1, p_a2], method="fdr_bh")[1][0]
    p_b = hypergeom.sf(1 - 1, 2, 2, 1)

    assert indexed.loc["A|TF1", "FDR"] == pytest.approx(expected_a1_fdr)
    assert indexed.loc["A|TF1", "Background_Size"] == 4
    assert indexed.loc["A|TF1", "GeneSet_Size"] == 2
    assert indexed.loc["B|TF3", "FDR"] == pytest.approx(p_b)
    assert indexed.loc["B|TF3", "Background_Size"] == 2
    assert indexed.loc["B|TF3", "GeneSet_Size"] == 1
    assert "A|TF2" not in indexed.index


def test_tf_ora_custom_background_keeps_unannotated_measurable_genes():
    matrix = pd.DataFrame({
        "Gene": ["G1", "G2", "G3"],
        "TF1": [1, 1, 0],
        "TF2": [0, 0, 1],
    })
    result = TFEnrichmentAnalyzer({"gene2tf": matrix}).ora(
        ["G1", "MEASURED_ONLY", "OUTSIDE"],
        background_genes={"G1", "G2", "G3", "G4", "MEASURED_ONLY", "G6"},
        min_overlap=1,
        min_size=1,
    )
    row = result.set_index("Term_ID").loc["TF1"]

    assert row["Input_GeneSet_Size"] == 3
    assert row["GeneSet_Size"] == 2
    assert row["Background_Size"] == 6
    assert row["Mapping_Rate"] == pytest.approx(2 / 3)
    assert row["Pvalue"] == pytest.approx(hypergeom.sf(0, 6, 2, 2))


def test_tf_ora_tissue_filter_redefines_the_annotated_background():
    matrix = pd.DataFrame({
        "Gene": ["G1", "G2", "G3", "G4"],
        "TFX|Lung": [1, 1, 0, 0],
        "TFX|Liver": [0, 0, 1, 1],
    })
    metadata = pd.DataFrame([
        {"Term_ID": "TFX|Lung", "TF": "TFX", "Library": "hTFtarget", "Context": "Lung"},
        {"Term_ID": "TFX|Liver", "TF": "TFX", "Library": "hTFtarget", "Context": "Liver"},
    ])
    result = TFEnrichmentAnalyzer(
        {"gene2tf": matrix, "tf_info": metadata}
    ).ora(["G1", "G3"], tissue="lung", min_overlap=1, min_size=1)

    assert result["Term_ID"].tolist() == ["TFX|Lung"]
    assert result.loc[0, "Background_Size"] == 2
    assert result.loc[0, "GeneSet_Size"] == 1


def test_animaltfdb_mapping_is_one_to_one_and_keeps_inference_provenance():
    mapper = OrthologMapper(
        human_tf_to_targets={
            "HUMAN_TF|Liver": {"HUMAN_A", "HUMAN_DUP", "HUMAN_B"}
        },
        species_to_human={
            "species_tf": "HUMAN_TF",
            "species_a": "HUMAN_A",
            "species_dup1": "HUMAN_DUP",
            "species_dup2": "HUMAN_DUP",
            "species_b": "HUMAN_B",
        },
        species_tf_set={"species_tf"},
        human_term_metadata={
            "HUMAN_TF|Liver": {
                "TF": "HUMAN_TF",
                "Context": "Liver",
                "Evidence_Type": "ChIP-seq",
            }
        },
    )

    mapped, reverse = mapper.map_tf_targets()

    assert mapped == {"species_tf|Liver": {"species_a", "species_b"}}
    assert "species_dup1" not in reverse
    assert "species_dup2" not in reverse
    assert mapper.mapped_term_metadata["species_tf|Liver"] == {
        "Term_Name": "species_tf [hTFtarget inferred; Liver]",
        "TF": "species_tf",
        "Library": "AnimalTFDB_hTFtarget",
        "Context": "Liver",
        "Evidence_Type": "ChIP-seq",
        "Inference_Type": "ortholog-inferred",
        "Human_TF": "HUMAN_TF",
    }


def test_tf_consensus_is_rank_based_and_stouffer_meta_is_disabled():
    results = {
        "ENCODE": pd.DataFrame({"TF": ["TF1", "TF2"], "Pvalue": [0.01, 0.02]}),
        "ReMap": pd.DataFrame({"TF": ["TF2", "TF1"], "Pvalue": [0.001, 0.5]}),
    }

    consensus = TFMetaAnalyzer.rank_consensus(results, method="meanrank")
    assert consensus["TF"].tolist() == ["TF1", "TF2"]
    assert consensus["Consensus_Score"].nunique() == 1
    assert set(consensus["Consensus_Method"]) == {"MeanRank"}
    with pytest.raises(ValueError, match="Stouffer"):
        TFMetaAnalyzer({}).combine_results(results, method="meta")
