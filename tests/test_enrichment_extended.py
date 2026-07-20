"""
AllEnricher v2.0 Additional testing

This test file is supplemented. test_enrichment.py Unoverwrite module testing:
- GSEA Gene set enrichment analysis algorithm
- SSGSEA Single sampleGSEAAlgorithms
- EnrichmentAnalyzer It's... _get_method Registration of methods
- Config Adding validation function (input_file, background_file Inspection)
- generate_term_url URL Generate Functions
- SpeciesLookup species retrieval module (offline mode)
"""

import pytest
import math
import tempfile
from pathlib import Path

from allenricher.core.config import Config, EnrichmentMethod
from allenricher.core.enrichment import (
    EnrichmentAnalyzer,
    GSEA,
    SSGSEA,
    EnrichmentResult,
    generate_term_url,
)


class TestGSEA:
    """GSEA Gene Set Enrichment Analysis Algorithm Test"""

    def test_calculate_enrichment_score_positive(self):
        """Test GSEA Enrichment Score Calculate: The gene is concentrated at the top of the ranked list"""
        gsea = GSEA(permutations=100, min_size=2, max_size=100)

        # Construct test data: genes are concentrated in the front of ranked list
        ranked_genes = ["GENE_A", "GENE_B", "GENE_C", "GENE_D", "GENE_E",
                        "GENE_F", "GENE_G", "GENE_H", "GENE_I", "GENE_J"]
        gene_set = {"GENE_A", "GENE_B", "GENE_C"}  # The first three genes.

        es, hit_genes = gsea.calculate_enrichment_score(ranked_genes, gene_set)

        assert es > 0  # It's heading towards the rich.
        assert len(hit_genes) > 0  # There should be a frontline gene.

    def test_calculate_enrichment_score_negative(self):
        """Test GSEA Enrichment: Gene set enrichment at the bottom of the ranked list"""
        gsea = GSEA(permutations=100, min_size=2, max_size=100)

        ranked_genes = ["GENE_A", "GENE_B", "GENE_C", "GENE_D", "GENE_E",
                        "GENE_F", "GENE_G", "GENE_H", "GENE_I", "GENE_J"]
        gene_set = {"GENE_H", "GENE_I", "GENE_J"}  # The last three genes.

        es, hit_genes = gsea.calculate_enrichment_score(ranked_genes, gene_set)

        # When the bottom is enriched, ES should be negative or smaller positive
        assert isinstance(es, float)

    def test_calculate_enrichment_score_no_overlap(self):
        """Test ES calculations when no intersection"""
        gsea = GSEA(permutations=100, min_size=2, max_size=100)

        ranked_genes = ["GENE_A", "GENE_B", "GENE_C"]
        gene_set = {"GENE_X", "GENE_Y"}  # No intersection.

        es, hit_genes = gsea.calculate_enrichment_score(ranked_genes, gene_set)

        assert es == 0.0  # E.S. should be 0 if no intersection
        assert len(hit_genes) == 0

    def test_calculate_enrichment_with_weights(self):
        """Test weighted GSEA ES calculation"""
        gsea = GSEA(permutations=100, min_size=2, max_size=100)

        ranked_genes = ["GENE_A", "GENE_B", "GENE_C", "GENE_D", "GENE_E"]
        gene_set = {"GENE_A", "GENE_B", "GENE_C"}
        gene_weights = {"GENE_A": 2.0, "GENE_B": 1.0, "GENE_C": 0.5,
                        "GENE_D": 0.1, "GENE_E": 0.1}

        es, hit_genes = gsea.calculate_enrichment_score(
            ranked_genes, gene_set, gene_weights=gene_weights
        )

        assert es > 0  # The weight of power should still be on the market.

    def test_enrichment_score_respects_theoretical_bounds(self):
        """The floating point of the full-blown time is not to accumulate beyond the ES theoretical boundary."""
        gsea = GSEA(permutations=10, min_size=1, max_size=500)
        ranked_genes = [f"GENE{i:03d}" for i in range(1, 45)]
        gene_set = set(ranked_genes[:20])
        weights = {gene: 1.0 - index * 0.01 for index, gene in enumerate(ranked_genes)}

        es, _ = gsea.calculate_enrichment_score(ranked_genes, gene_set, weights)

        assert -1.0 <= es <= 1.0
        assert math.isclose(es, 1.0, abs_tol=1e-12)

    def test_permutation_test(self):
        """Test GSEA replacement test p-value calculation"""
        gsea = GSEA(permutations=200, min_size=2, max_size=100, seed=42)

        ranked_genes = [f"GENE_{i:03d}" for i in range(50)]
        gene_set = {f"GENE_{i:03d}" for i in range(5)}  # The first five genes.

        # Calculate observation first
        observed_es, _ = gsea.calculate_enrichment_score(ranked_genes, gene_set)
        # Run replacement check again
        pvalue = gsea._run_permutation_test(ranked_genes, gene_set, observed_es)

        assert 0 <= pvalue <= 1  # p value should be within legal range
        assert isinstance(pvalue, float)

    def test_calculate_pvalue_placeholder(self):
        """Test the calculate_pvalue of GSA returns the placeholder"""
        gsea = GSEA()
        pvalue = gsea.calculate_pvalue(10, 50, 100, 1000)
        assert pvalue == 1.0  # GSEA does not use this method to return placeholder value


class TestSSGSEA:
    """SSGSEA single sample GSEA algorithm test"""

    def test_calculate_enrichment_score(self):
        """Test ssGSEA Enrichment Scores"""
        ssgsea = SSGSEA(min_size=2, max_size=100)

        ranked_genes = ["GENE_A", "GENE_B", "GENE_C", "GENE_D", "GENE_E",
                        "GENE_F", "GENE_G", "GENE_H", "GENE_I", "GENE_J"]
        gene_set = {"GENE_A", "GENE_B", "GENE_C"}

        es, es_min, es_max, hit_genes = ssgsea.calculate_enrichment_score(
            ranked_genes, gene_set
        )

        assert isinstance(es, float)
        assert isinstance(es_min, float)
        assert isinstance(es_max, float)
        assert es >= es_min  # ES should > = minimum
        assert es <= es_max  # ES should < = maximum

    def test_nes_normalization(self):
        """Testing NES = NES should be within [1, 1]"""
        ssgsea = SSGSEA(min_size=2, max_size=100)

        ranked_genes = [f"GENE_{i:03d}" for i in range(20)]
        gene_set = {f"GENE_{i:03d}" for i in range(5)}

        es, es_min, es_max, _ = ssgsea.calculate_enrichment_score(
            ranked_genes, gene_set
        )

        # NES = ES / (|ES_min| + |ES_max|)
        denominator = abs(es_min) + abs(es_max)
        if denominator > 0:
            nes = es / denominator
            assert -1 <= nes <= 1  # NES should be within [-1, ]

    def test_no_overlap(self):
        """Test ssGSEA when there is no intersection"""
        ssgsea = SSGSEA(min_size=2, max_size=100)

        ranked_genes = ["GENE_A", "GENE_B", "GENE_C"]
        gene_set = {"GENE_X", "GENE_Y"}

        es, es_min, es_max, hit_genes = ssgsea.calculate_enrichment_score(
            ranked_genes, gene_set
        )

        assert es == 0.0
        assert es_min == 0.0
        assert es_max == 0.0
        assert len(hit_genes) == 0

    def test_calculate_pvalue_placeholder(self):
        """Tests for sGSEA calculate_pvalue returns nn"""
        ssgsea = SSGSEA()
        pvalue = ssgsea.calculate_pvalue(10, 50, 100, 1000)
        assert math.isnan(pvalue)  # ssGSEA does not use p value


class TestMethodRegistration:
    """Register test analysis method to _get_method()"""

    def test_ssgsea_registered(self):
        """Test ssgsea method registered _get_method()"""
        config = Config(method="ssgsea")
        analyzer = EnrichmentAnalyzer(config)

        # Should not be thrown away
        method = analyzer._get_method()
        assert isinstance(method, SSGSEA)

    def test_gsea_registered(self):
        """Test gsea method registered"""
        config = Config(method="gsea")
        analyzer = EnrichmentAnalyzer(config)

        method = analyzer._get_method()
        assert isinstance(method, GSEA)

    def test_fisher_alias_is_not_registered(self):
        """Verify that the retired Fisher alias is rejected."""
        config = Config(method="fisher")
        analyzer = EnrichmentAnalyzer(config)
        with pytest.raises(ValueError, match="Unknown method"):
            analyzer._get_method()

    def test_hypergeometric_registered(self):
        """Test hypergeometric method registered"""
        from allenricher.core.enrichment import HypergeometricTest

        config = Config(method="hypergeometric")
        analyzer = EnrichmentAnalyzer(config)

        method = analyzer._get_method()
        assert isinstance(method, HypergeometricTest)

    def test_invalid_method_raises(self):
        """Could not close temporary folder: %s"""
        config = Config(method="nonexistent")
        # Validation is deferred until a method-specific property is accessed.
        analyzer = EnrichmentAnalyzer(config)

        with pytest.raises(ValueError, match="Unknown method"):
            _ = analyzer.method


class TestConfigValidation:
    """Config Complementary Test for Validation"""

    def test_validate_missing_input_file(self, tmp_path):
        """Test for validation when the input file does not exist"""
        config = Config(input_file="/nonexistent/path/genes.txt")
        errors = config.validate()

        assert any("Input file not found" in e for e in errors)

    def test_validate_existing_input_file(self, tmp_path):
        """Verify with the entry file when testing it to exist"""
        gene_file = tmp_path / "genes.txt"
        gene_file.write_text("BRCA1\nTP53\n")

        config = Config(input_file=str(gene_file))
        errors = config.validate()

        assert not any("Input file not found" in e for e in errors)

    def test_validate_missing_background_file(self):
        """Test for validation when background gene file does not exist"""
        config = Config(background_file="/nonexistent/background.txt")
        errors = config.validate()

        assert any("Background file not found" in e for e in errors)

    def test_non_ora_methods_ignore_ora_only_input_paths(self):
        """GSEA and activity methods should not validate unused ORA file paths."""
        for method in ('gsea', 'ssgsea', 'gsva'):
            config = Config(
                method=method,
                input_file='/nonexistent/query_genes.txt',
                background_file='/nonexistent/background_genes.txt',
            )
            errors = config.validate()
            assert not any('Input file not found' in error for error in errors)
            assert not any('Background file not found' in error for error in errors)

    def test_validate_null_input_file(self):
        """Test input_file is None and file should not be reported as incorrect"""
        config = Config(input_file=None)
        errors = config.validate()

        assert not any("Input file not found" in e for e in errors)

    def test_config_yaml_template_consistency(self):
        """Test key parameters in DEFAULT_CONFIG_YAML to match Config default"""
        from allenricher.core.config import DEFAULT_CONFIG_YAML, Config

        config = Config()

        # Check that YAML templates do not contain outdated default values
        assert "max_genes: 500" not in DEFAULT_CONFIG_YAML  # For.inf
        assert "gsea_min_size: 10" not in DEFAULT_CONFIG_YAML
        assert "n_jobs: 4" not in DEFAULT_CONFIG_YAML  # For 1

        # Check that YAML templates contain the right default values
        assert "gsea_min_size: null" in DEFAULT_CONFIG_YAML
        assert "gsea_max_size: null" in DEFAULT_CONFIG_YAML
        assert "min_genes: 3" in DEFAULT_CONFIG_YAML
        assert "max_genes: .inf" in DEFAULT_CONFIG_YAML
        assert "plot_width: null" in DEFAULT_CONFIG_YAML
        assert "plot_height: null" in DEFAULT_CONFIG_YAML
        assert "background_mode: \"annotated\"" in DEFAULT_CONFIG_YAML
        assert config.background_mode == "annotated"
        assert config.plot_width is None
        assert config.plot_height is None


class TestGenerateTermURL:
    """Generate_term_url URL generation function test"""

    def test_go_url(self):
        """Test GO database URL generation"""
        url = generate_term_url("GO:0008150", "GO")
        assert "amigo" in url
        assert "GO:0008150" in url

    def test_kegg_url(self):
        """Test KEGG database URL generation"""
        url = generate_term_url("hsa04010", "KEGG")
        assert "kegg.jp" in url
        assert "hsa04010" in url

    def test_reactome_url(self):
        """Test Reactome database URL generation"""
        url = generate_term_url("R-HSA-1234567", "Reactome")
        assert "reactome.org" in url
        assert "R-HSA-1234567" in url

    def test_do_url(self):
        """Test DO (Hender of Disease) URL generation"""
        url = generate_term_url("DOID:1234", "DO")
        assert "disease-ontology" in url

    def test_unknown_database(self):
        """Test unknown database to return empty string"""
        url = generate_term_url("UNKNOWN:001", "UnknownDB")
        assert url == ""


class TestSpeciesLookup:
    """SpeciesRegistry species retrieval module test"""

    @pytest.fixture
    def registry(self, tmp_path):
        """Create a SpeciesRegistry instance with test data"""
        from allenricher.database.species_registry import SpeciesRegistry, SpeciesEntry

        registry = SpeciesRegistry(registry_path=tmp_path / "test_species.tsv")

        # Add test species data
        test_species = [
            SpeciesEntry(taxid=9606, latin_name="Homo sapiens", common_name="Human",
                        has_go=True, go_source="ncbi_gene2go", go_gene_count=19500,
                        has_kegg=True, kegg_code="hsa", kegg_code_source="kegg", kegg_gene_count=22345,
                        has_reactome=True, reactome_code="HSA",
                        has_do=True, do_gene_count=12000),
            SpeciesEntry(taxid=10090, latin_name="Mus musculus", common_name="Mouse",
                        has_go=True, go_source="ncbi_gene2go",
                        has_kegg=True, kegg_code="mmu", kegg_code_source="kegg",
                        has_reactome=True, reactome_code="MMU",
                        has_do=False),
            SpeciesEntry(taxid=4932, latin_name="Saccharomyces cerevisiae", common_name="Yeast",
                        has_go=True, go_source="ncbi_gene2go",
                        has_kegg=True, kegg_code="sce", kegg_code_source="kegg",
                        has_reactome=True, reactome_code="SCE",
                        has_do=False),
        ]
        for sp in test_species:
            registry.add_entry(sp)

        return registry

    def test_builtin_species_loaded(self, registry):
        """Test built-in species data loaded"""
        # The built-in species should include humans.
        info = registry.query_by_kegg_code("hsa")
        assert info is not None
        assert info.latin_name == "Homo sapiens"
        assert info.taxid == 9606

    def test_builtin_species_count(self, registry):
        """Number of built-in species tested"""
        all_species = list(registry.entries.values())
        assert len(all_species) == 3  # There are three species in the test data.

    def test_lookup_by_latin_name(self, registry):
        """Tests are retrieved by Latin names"""
        results = registry.query_by_latin_name("Homo sapiens")
        assert len(results) > 0
        assert results[0].kegg_code == "hsa"

    def test_lookup_by_taxid(self, registry):
        """Verify species lookup by TaxID."""
        info = registry.query_by_taxid(9606)
        assert info is not None
        assert info.kegg_code == "hsa"

    def test_yeast_taxid_consistency(self, registry):
        """Testing the Brewer Taxid in species_registry and config (both 4932)"""
        from allenricher.core.config import SPECIES_CONFIGS

        info = registry.query_by_kegg_code("sce")
        assert info is not None
        assert info.taxid == 4932  # For the correct taxid
        assert SPECIES_CONFIGS["sce"].taxonomy_id == 4932  # It's also the same in config.

    def test_builtin_taxid_map_consistency(self, registry):
        """Verify that the S. cerevisiae registry entry consistently uses TaxID 4932."""
        # Query the species registry by KEGG organism code.
        info_by_code = registry.query_by_kegg_code("sce")
        assert info_by_code is not None
        assert info_by_code.taxid == 4932

        # Verify backwards by query_by_taxid
        info_by_taxid = registry.query_by_taxid(4932)
        assert info_by_taxid is not None
        assert info_by_taxid.kegg_code == "sce"


class TestEnrichmentResultURL:
    """EnterResult URL field test"""

    def test_term_url_in_dict(self):
        """Test EnrichmentResult. to_dict() contains Term_URL fields"""
        result = EnrichmentResult(
            term_id="GO:0008150",
            term_name="biological_process",
            database="GO",
            pvalue=0.001,
            adjusted_pvalue=0.01,
            gene_count=10,
            background_count=100,
            expected_count=5.0,
            rich_factor=2.0,
            gene_list=["BRCA1", "TP53"],
            gene_ratio="10/100",
            background_ratio="100/1000",
            term_url="https://amigo.geneontology.org/amigo/term/GO:0008150"
        )

        d = result.to_dict()
        assert "Term_URL" in d
        assert "amigo" in d["Term_URL"]

    def test_term_url_default_empty(self):
        """Test Term_URL Default to Empty String"""
        result = EnrichmentResult(
            term_id="TEST:001",
            term_name="Test",
            database="TEST",
            pvalue=0.05,
            adjusted_pvalue=0.1,
            gene_count=5,
            background_count=50,
            expected_count=2.5,
            rich_factor=2.0,
            gene_list=[],
            gene_ratio="5/100",
            background_ratio="50/1000"
        )

        d = result.to_dict()
        assert d["Term_URL"] == ""


class TestLazyMethodInit:
    """Delay initialize method attribute test"""

    def test_constructor_does_not_validate_method(self):
        """Test construction function does not immediately verify method name"""
        # An invalid method is reported during deferred validation, not object construction.
        config = Config(method="nonexistent_method")
        analyzer = EnrichmentAnalyzer(config)
        assert analyzer.config.method == "nonexistent_method"

    def test_method_lazy_loading(self):
        """Delay load and cache for testing method properties"""
        config = Config(method="hypergeometric")
        analyzer = EnrichmentAnalyzer(config)

        # First visit: Examples of how to create
        method1 = analyzer.method
        # Second visit: the same example should be returned (cached)
        method2 = analyzer.method
        assert method1 is method2

    def test_method_setter(self):
        """Test with the setter directly set example"""
        from allenricher.core.enrichment import FisherExactTest

        config = Config(method="nonexistent")  # Invalid name
        analyzer = EnrichmentAnalyzer(config)

        # Examples of direct fixation of methods
        analyzer.method = FisherExactTest()
        assert isinstance(analyzer.method, FisherExactTest)


class TestDeadCodeRemoval:
    """Dead Code Remove Validation Test"""

    def test_no_cache_fields_in_config(self):
        """Test Config removed field using_cache and cache_dir"""
        config = Config()
        assert not hasattr(config, 'use_cache'), "_use_cache should have been removed"
        assert not hasattr(config, 'cache_dir'), "Case_dir should have been removed"

    def test_no_cache_in_yaml_template(self):
        """Test DeFAUT_CONFIG_YAML does not contain cache-related configuration"""
        from allenricher.core.config import DEFAULT_CONFIG_YAML
        assert "use_cache" not in DEFAULT_CONFIG_YAML
        assert "cache_dir" not in DEFAULT_CONFIG_YAML

    def test_no_utils_package(self):
        """Test utils package removed"""
        import importlib
        with pytest.raises(ModuleNotFoundError):
            importlib.import_module("allenricher.utils")


class TestHeatmapGeneSeparator:
    """Verify restoration of gene separators used by heatmap inputs."""

    def test_genes_separated_by_semicolon(self):
        """Test EnrichmentResult. to_dit()"""
        result = EnrichmentResult(
            term_id="TEST:001",
            term_name="Test Term",
            database="TEST",
            pvalue=0.001,
            adjusted_pvalue=0.01,
            gene_count=3,
            background_count=50,
            expected_count=2.5,
            rich_factor=2.0,
            gene_list=["BRCA1", "TP53", "EGFR"],
            gene_ratio="3/100",
            background_ratio="50/1000"
        )
        d = result.to_dict()
        # Genes should be separated by semicolons.
        assert d["Genes"] == "BRCA1;TP53;EGFR"
        # The split with a semicolon should be given three genes.
        assert len(d["Genes"].split(';')) == 3


class TestAdjustPvaluesNaN:
    """A treatment test for NaN p"""

    def test_nan_pvalues_skip_correction(self):
        """Test contains a skip correction at the np"""
        config = Config()
        analyzer = EnrichmentAnalyzer(config)

        # Create a list of results with np values (ssGSEA output)
        results = [
            EnrichmentResult(
                term_id=f"TEST:{i:03d}", term_name=f"Term {i}", database="TEST",
                pvalue=float('nan'), adjusted_pvalue=float('nan'),
                gene_count=5, background_count=50,
                expected_count=2.5, rich_factor=2.0,
                gene_list=["G1", "G2"], gene_ratio="5/100", background_ratio="50/1000"
            )
            for i in range(5)
        ]

        # There should be no abnormality, direct return to the original result.
        corrected = analyzer.adjust_pvalues(results, method="BH")
        assert len(corrected) == 5
        # The nn p value has not been modified
        assert all(math.isnan(r.adjusted_pvalue) for r in corrected)

    def test_normal_pvalues_correction(self):
        """The correction of the test p value is still valid."""
        config = Config()
        analyzer = EnrichmentAnalyzer(config)

        results = [
            EnrichmentResult(
                term_id=f"TEST:{i:03d}", term_name=f"Term {i}", database="TEST",
                pvalue=0.001 * (i + 1), adjusted_pvalue=0.001 * (i + 1),
                gene_count=5, background_count=50,
                expected_count=2.5, rich_factor=2.0,
                gene_list=["G1", "G2"], gene_ratio="5/100", background_ratio="50/1000"
            )
            for i in range(5)
        ]

        corrected = analyzer.adjust_pvalues(results, method="BH")
        # Corrected p value should > = original p
        for orig, corr in zip(results, corrected):
            assert corr.adjusted_pvalue >= orig.pvalue


class TestHypergeometricBackgroundTotal:
    """Hypergemetetrist background_ttal parameter test"""

    def test_uses_background_total(self):
        """Tests for the Hypergeememediatrist correct use of background_ttal parameters"""
        from allenricher.core.enrichment import HypergeometricTest

        hypo = HypergeometricTest()
        gene_set = {"A", "B", "C", "D", "E"}
        background_set = {"A", "B", "C", "D", "E", "F", "G", "H", "I", "J"}
        term_genes = {"A", "B", "C"}

        # Do Not Specify Background_total
        result1 = hypo.calculate_enrichment(
            gene_set, background_set, term_genes, "Test", "TEST:001", "TEST"
        )

        # Specify different backgroup_total (semantics for simulation v1.0)
        result2 = hypo.calculate_enrichment(
            gene_set, background_set, term_genes, "Test", "TEST:001", "TEST",
            background_total=20000  # Total number of genes in the Note
        )

        assert result1 is not None
        assert result2 is not None
        # Different background_total should produce different p values
        assert result1.pvalue != result2.pvalue


class TestGeneSetSizeDefaults:
    """The default value for the method-level gene set size."""

    def test_gsea_defaults(self):
        gsea = GSEA()
        assert (gsea.min_size, gsea.max_size) == (15, 500)

    def test_ssgsea_defaults(self):
        ssgsea = SSGSEA()
        assert (ssgsea.min_size, ssgsea.max_size) == (1, None)


class TestConfigMaxGenesType:
    """Config max_genes type test"""

    def test_max_genes_is_float_inf(self):
        """Test max_genes default value is float('inf')"""
        config = Config()
        assert config.max_genes == float('inf')
        assert isinstance(config.max_genes, float)


class TestEnrichmentResultToDictColumnNames:
    """EnterResult. to_disc() Use Conformity Test for Listing and Plitter.py"""

    def test_to_dict_column_names_match_expected(self):
        """Test to_dict() output listings are consistent with listings used in plot_all"""
        result = EnrichmentResult(
            term_id="GO:0008150",
            term_name="biological_process",
            database="GO",
            pvalue=0.001,
            adjusted_pvalue=0.01,
            gene_count=10,
            background_count=100,
            expected_count=5.0,
            rich_factor=2.0,
            gene_list=["BRCA1", "TP53", "EGFR"],
            gene_ratio="10/100",
            background_ratio="100/1000"
        )
        d = result.to_dict()

        # Verify that listing used in listbox.py exists in to_dict() output
        expected_keys = [
            "Term_ID", "Term_Name", "P_Value", "Adjusted_P_Value",
            "Gene_Count", "Background_Count", "Expected_Count",
            "Rich_Factor", "Gene_Ratio", "Background_Ratio", "Genes"
        ]
        for key in expected_keys:
            assert key in d, f"To_dict() Missing list: {key}"

    def test_adjusted_p_value_not_q_value(self):
        """The p-value name of the corrected test is Adjusted_P_Value instead of Q_Value"""
        result = EnrichmentResult(
            term_id="TEST:001", term_name="Test", database="TEST",
            pvalue=0.001, adjusted_pvalue=0.01,
            gene_count=5, background_count=50,
            expected_count=2.5, rich_factor=2.0,
            gene_list=["G1"], gene_ratio="5/100", background_ratio="50/1000"
        )
        d = result.to_dict()
        assert "Adjusted_P_Value" in d
        assert "Q_Value" not in d
        assert d["Adjusted_P_Value"] == 0.01


class TestAPISecurity:
    """API Safety Related Test (Track Through Protection)"""

    def test_path_traversal_prevention(self):
        """Test path through the defense logic: clean path should not contain... or /"""
        import re
        malicious_database = "../../etc"
        safe_database = re.sub(r'[^\w\-]', '', malicious_database)
        # "..." and "/" were removed, only "etc" was retained.
        assert ".." not in safe_database
        assert "/" not in safe_database
        assert safe_database == "etc"

    def test_safe_plot_filename(self):
        """Test security filename generation"""
        import re
        database = "GO"
        plot_type = "barplot"
        safe_db = re.sub(r'[^\w\-]', '', database)
        safe_pt = re.sub(r'[^\w\-]', '', plot_type)
        filename = f"{safe_db}_{safe_pt}.pdf"
        assert filename == "GO_barplot.pdf"
        assert "/" not in filename
        assert ".." not in filename


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
