from unittest.mock import Mock, patch

import pytest

from allenricher.database.chea3_fetcher import ChEA3Fetcher


def test_chea3_download_urls_match_current_official_asset_names():
    assert ChEA3Fetcher.CHEA3_GMT_LIBS == {
        "ENCODE": "https://maayanlab.cloud/chea3/assets/tflibs/ENCODE_ChIP-seq.gmt",
        "ReMap": "https://maayanlab.cloud/chea3/assets/tflibs/ReMap_ChIP-seq.gmt",
        "LiteratureChIP": "https://maayanlab.cloud/chea3/assets/tflibs/Literature_ChIP-seq.gmt",
        "GTExCoexpression": "https://maayanlab.cloud/chea3/assets/tflibs/GTEx_Coexpression.gmt",
        "ARCHS4Coexpression": "https://maayanlab.cloud/chea3/assets/tflibs/ARCHS4_Coexpression.gmt",
        "EnrichrQueries": "https://maayanlab.cloud/chea3/assets/tflibs/Enrichr_Queries.gmt",
    }


def test_chea3_species_list_is_human_only():
    assert ChEA3Fetcher.get_supported_species() == ["Homo sapiens"]


def test_chea3_fetcher_rejects_html_saved_as_gmt(tmp_path):
    response = Mock(content=b"<html>not a GMT</html>")
    response.raise_for_status.return_value = None
    with patch("allenricher.database.chea3_fetcher.requests.get", return_value=response):
        with pytest.raises(ValueError, match="GMT"):
            ChEA3Fetcher(str(tmp_path)).download_gmt_library("ENCODE")

    assert not (tmp_path / "chea3" / "ChEA3v2024" / "ENCODE_tf.gmt").exists()
