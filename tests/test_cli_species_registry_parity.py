from types import SimpleNamespace

from allenricher.cli import _visible_database_support


def test_list_species_default_outputs_include_tf_database_fields():
    from pathlib import Path

    source = (Path(__file__).resolve().parents[1] / "allenricher" / "cli.py").read_text(encoding="utf-8-sig")

    for label in ("TRRUST", "ChEA3", "AnimalTFDB", "hTFtarget"):
        assert label in source
    for field in ("has_trrust", "has_chea3", "has_animaltfdb", "has_htftarget"):
        assert field in source
    assert "has_go\\thas_kegg\\thas_reactome\\thas_do\\thas_disgenet\\thas_wikipathways\\thas_trrust\\thas_chea3\\thas_animaltfdb\\thas_htftarget" in source


def test_species_display_uses_runtime_tf_species_rules():
    human = SimpleNamespace(
        taxid=9606,
        kegg_code="hsa",
        has_trrust=True,
        has_chea3=True,
        has_animaltfdb=True,
        has_htftarget=True,
    )
    fly = SimpleNamespace(taxid=7227, kegg_code="dme", has_animaltfdb=True)

    assert _visible_database_support(human, "TRRUST", "has_trrust") is True
    assert _visible_database_support(human, "ChEA3", "has_chea3") is True
    assert _visible_database_support(human, "hTFtarget", "has_htftarget") is True
    assert _visible_database_support(human, "AnimalTFDB", "has_animaltfdb") is False
    assert _visible_database_support(fly, "AnimalTFDB", "has_animaltfdb") is True