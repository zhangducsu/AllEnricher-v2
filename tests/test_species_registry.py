import pytest
import tempfile
from pathlib import Path
from allenricher.database.species_registry import SpeciesRegistry, SpeciesEntry


class TestSpeciesEntry:
    """SpeciesEntry Data class test"""

    def test_create_minimal_entry(self):
        """Test to minimize creation"""
        entry = SpeciesEntry(taxid=9606, latin_name="Homo sapiens")
        assert entry.taxid == 9606
        assert entry.latin_name == "Homo sapiens"
        assert entry.has_go is False
        assert entry.has_kegg is False
        assert entry.go_source is None

    def test_create_full_entry(self):
        """Test full field creation"""
        entry = SpeciesEntry(
            taxid=9606, latin_name="Homo sapiens", common_name="Human",
            has_go=True, go_source="ncbi_gene2go", go_gene_count=19500, go_term_count=18500,
            has_kegg=True, kegg_code="hsa", kegg_code_source="kegg", kegg_gene_count=22345, kegg_pathway_count=350,
            has_reactome=True, reactome_code="HSA", reactome_gene_count=10500, reactome_pathway_count=1500,
            has_do=True, do_gene_count=12000, do_term_count=8000
        )
        assert entry.taxid == 9606
        assert entry.has_go is True
        assert entry.kegg_code == "hsa"
        assert entry.kegg_code_source == "kegg"


class TestSpeciesRegistry:
    """SpeciesRegistry Test"""

    @pytest.fixture
    def temp_registry(self, tmp_path):
        """Create a temporary species registry."""
        registry = SpeciesRegistry(registry_path=tmp_path / "supported_species.tsv")
        # Add Test Data
        entries = [
            SpeciesEntry(taxid=9606, latin_name="Homo sapiens", common_name="Human",
                has_go=True, go_source="ncbi_gene2go", go_gene_count=19500, go_term_count=18500,
                has_kegg=True, kegg_code="hsa", kegg_code_source="kegg", kegg_gene_count=22345, kegg_pathway_count=350,
                has_reactome=True, reactome_code="HSA", reactome_gene_count=10500, reactome_pathway_count=1500,
                has_do=True, do_gene_count=12000, do_term_count=8000,
                has_disgenet=True),
            SpeciesEntry(taxid=10090, latin_name="Mus musculus", common_name="Mouse",
                has_go=True, go_source="ncbi_gene2go", go_gene_count=18200, go_term_count=16800,
                has_kegg=True, kegg_code="mmu", kegg_code_source="kegg", kegg_gene_count=25245, kegg_pathway_count=320,
                has_reactome=True, reactome_code="MMU", reactome_gene_count=9800, reactome_pathway_count=1400,
                has_do=False),
            SpeciesEntry(taxid=3702, latin_name="Arabidopsis thaliana", common_name="Thale cress",
                has_go=True, go_source="uniprot_goa", go_gene_count=15000, go_term_count=12000,
                has_kegg=False,
                has_reactome=False,
                has_do=False),
        ]
        for e in entries:
            registry.add_entry(e)
        registry.save()
        return registry

    def test_save_and_load(self, temp_registry):
        """Test Save and Load"""
        # Reload
        loaded = SpeciesRegistry(registry_path=temp_registry.registry_path)
        loaded.load()
        assert len(loaded.entries) == 3
        assert loaded.query_by_taxid(9606) is not None
        assert loaded.query_by_taxid(10090) is not None
        assert loaded.query_by_taxid(3702) is not None

    def test_load_default_prefers_current_basic_registry(self, tmp_path):
        legacy = SpeciesRegistry(tmp_path / "supported_species.tsv")
        legacy.add_entry(SpeciesEntry(taxid=1, latin_name="Legacy species"))
        legacy.save()
        current = SpeciesRegistry(tmp_path / "basic" / "supported_species.tsv")
        current.add_entry(SpeciesEntry(taxid=2, latin_name="Current species"))
        current.save()

        loaded = SpeciesRegistry.load_default(tmp_path)

        assert loaded.registry_path == current.registry_path
        assert set(loaded.entries) == {2}

    def test_default_registry_adds_authoritative_tf_database_coverage(self, tmp_path):
        from allenricher.database.animaltfdb_fetcher import AnimalTFDBFetcher
        from allenricher.database.downloader import DataDownloader

        downloader = DataDownloader(root_dir=str(tmp_path))
        downloader.record_database_species(
            "TRRUST", [(9606, "Homo sapiens"), (10090, "Mus musculus")]
        )
        downloader.record_database_species("ChEA3", [(9606, "Homo sapiens")])
        downloader.record_database_species(
            "AnimalTFDB", AnimalTFDBFetcher.get_supported_species_records()
        )
        downloader.record_database_species("hTFtarget", [(9606, "Homo sapiens")])
        downloader.record_database_species("DisGeNET", [(9606, "Homo sapiens")])
        downloader.refresh_supported_species_registry()

        loaded = SpeciesRegistry.load_default(tmp_path)
        summary = loaded.get_summary()

        assert summary["disgenet"]["count"] == 1
        assert summary["trrust"]["count"] == 2
        assert summary["chea3"]["count"] == 1
        assert summary["animaltfdb"]["count"] == 183
        assert summary["htftarget"]["count"] == 1
        assert len(loaded.filter_by_databases(animaltfdb=True)) == 183
        assert [item.taxid for item in loaded.filter_by_databases(chea3=True)] == [9606]
        assert [item.taxid for item in loaded.filter_by_databases(htftarget=True)] == [9606]

    def test_query_by_taxid(self, temp_registry):
        """Test for TaxID"""
        entry = temp_registry.query_by_taxid(9606)
        assert entry is not None
        assert entry.latin_name == "Homo sapiens"
        assert entry.has_go is True
        assert entry.kegg_code == "hsa"

        # Query does not exist
        assert temp_registry.query_by_taxid(99999) is None

    def test_query_by_latin_name(self, temp_registry):
        """Test for Latin names"""
        results = temp_registry.query_by_latin_name("Homo sapiens")
        assert len(results) == 1
        assert results[0].taxid == 9606

        # Fuzzy Query (case insensitive)
        results = temp_registry.query_by_latin_name("mus")
        assert len(results) == 1
        assert results[0].taxid == 10090

        # Substring Query
        results = temp_registry.query_by_latin_name("Arabidopsis")
        assert len(results) == 1

        # No match
        results = temp_registry.query_by_latin_name("NotExists")
        assert len(results) == 0

    def test_query_by_kegg_code(self, temp_registry):
        """Test for Kegg code Query"""
        entry = temp_registry.query_by_kegg_code("hsa")
        assert entry is not None
        assert entry.taxid == 9606

        entry = temp_registry.query_by_kegg_code("mmu")
        assert entry is not None
        assert entry.taxid == 10090

        # Cannot initialise Evolution's mail component.
        assert temp_registry.query_by_kegg_code("xxx") is None

    def test_filter_by_databases(self, temp_registry):
        """Test Filter on Database"""
        # Filter only Go
        results = temp_registry.filter_by_databases(go=True)
        assert len(results) == 3

        # Filter Kegg only
        results = temp_registry.filter_by_databases(kegg=True)
        assert len(results) == 2

        # GO + KEGG
        results = temp_registry.filter_by_databases(go=True, kegg=True)
        assert len(results) == 2

        # DO
        results = temp_registry.filter_by_databases(do=True)
        assert len(results) == 1

        results = temp_registry.filter_by_databases(disgenet=True)
        assert len(results) == 1

        # Reactome
        results = temp_registry.filter_by_databases(reactome=True)
        assert len(results) == 2

    def test_get_summary(self, temp_registry):
        """Statistical summary of tests"""
        summary = temp_registry.get_summary()
        assert summary['total_species'] == 3
        assert summary['go']['count'] == 3
        assert summary['kegg']['count'] == 2
        assert summary['reactome']['count'] == 2
        assert summary['do']['count'] == 1
        assert summary['disgenet']['count'] == 1

    def test_get_species_detail(self, temp_registry):
        """Test species details"""
        detail = temp_registry.get_species_detail(9606)
        assert detail is not None
        assert detail['taxid'] == 9606
        assert detail['latin_name'] == 'Homo sapiens'
        assert detail['has_go'] is True
        assert detail['kegg_code'] == 'hsa'

        # Cannot initialise Evolution's mail component.
        assert temp_registry.get_species_detail(99999) is None

    def test_generate_kegg_abbreviation(self):
        """Test KEG Abbreviation"""
        assert SpeciesRegistry.generate_kegg_abbreviation("Homo sapiens") == "hsa"
        assert SpeciesRegistry.generate_kegg_abbreviation("Mus musculus") == "mmu"
        assert SpeciesRegistry.generate_kegg_abbreviation("Arabidopsis thaliana") == "ath"
        assert SpeciesRegistry.generate_kegg_abbreviation("Rattus norvegicus") == "rno"
        assert SpeciesRegistry.generate_kegg_abbreviation("Danio rerio") == "dre"
        assert SpeciesRegistry.generate_kegg_abbreviation("Saccharomyces cerevisiae") == "sce"
        assert SpeciesRegistry.generate_kegg_abbreviation("Caenorhabditis elegans") == "cel"

    def test_add_entry_overwrite(self, temp_registry):
        """Test Add Entry Overwrite"""
        new_entry = SpeciesEntry(taxid=9606, latin_name="Homo sapiens Updated",
            has_go=True, go_source="uniprot_goa")
        temp_registry.add_entry(new_entry)
        assert len(temp_registry.entries) == 3  # No increase
        assert temp_registry.query_by_taxid(9606).latin_name == "Homo sapiens Updated"
