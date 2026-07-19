import gzip

from allenricher.database.parsers.animaltfdb import AnimalTFDBParser


def test_official_ortholog_ids_are_mapped_to_symbols_and_best_hit_is_kept(tmp_path):
    gene_info = tmp_path / "gene_info.gz"
    with gzip.open(gene_info, "wt", encoding="utf-8") as handle:
        handle.write("#tax_id\tGeneID\tSymbol\tLocusTag\tSynonyms\tdbXrefs\n")
        handle.write("7227\t1\tflyA\t-\t-\tFlyBase:FBgn0001\n")
        handle.write("9606\t2\tHUMAN_LOW\t-\t-\tEnsembl:ENSG_LOW\n")
        handle.write("9606\t3\tHUMAN_BEST\t-\t-\tEnsembl:ENSG_BEST\n")
    ortholog = tmp_path / "Drosophila_melanogaster_ortholog_to_human"
    ortholog.write_text(
        "Species\tEnsembl ID\tCoverage\tIdentity\tOrtholog ID\tCoverage\tIdentity\tOrtholog Species\n"
        "Drosophila_melanogaster\tFBgn0001\t30\t20\tENSG_LOW\t30\t20\tHomo sapiens\n"
        "Drosophila_melanogaster\tFBgn0001\t90\t80\tENSG_BEST\t90\t80\tHomo sapiens\n",
        encoding="utf-8",
    )

    mappings = AnimalTFDBParser.load_external_id_symbol_maps(str(gene_info), {7227, 9606})
    species_map = mappings[7227]
    human_map = mappings[9606]
    result = AnimalTFDBParser.parse_ortholog_to_human(
        str(ortholog), species_map, human_map
    )

    assert result == {"flyA": "HUMAN_BEST"}


def test_historical_two_column_symbol_file_remains_supported(tmp_path):
    ortholog = tmp_path / "legacy.tsv"
    ortholog.write_text("flyA\tHUMAN_A\nflyB\tHUMAN_B\n", encoding="utf-8")

    assert AnimalTFDBParser.parse_ortholog_to_human(str(ortholog)) == {
        "flyA": "HUMAN_A",
        "flyB": "HUMAN_B",
    }
