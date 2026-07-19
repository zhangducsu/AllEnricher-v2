"""
AllEnricher v2.0 Unit Test

Test coverage:
- Config Default for the configuration class, Certification and species configuration
- FisherExactTest Fisher's p-calculation and enrichment analysis
- HypergeometricTest It's a geometric test.pValue calculation andFisherTest consistency
- EnrichmentResult Results in Data Class Dictionary Conversion
- EnrichmentAnalyzer Loading of the gene list for the analysis engine, pValue correction and result filtering
- Enumeration class (EnrichmentMethod, CorrectionMethod)Complete check
"""

import pytest
import pandas as pd
import numpy as np
from pathlib import Path
import tempfile

from allenricher.core.config import Config, EnrichmentMethod, CorrectionMethod
from allenricher.core.enrichment import (
    EnrichmentAnalyzer,
    FisherExactTest,
    HypergeometricTest,
    EnrichmentResult
)


class TestConfig:
    """Configuration class (Config) test"""

    def test_default_config(self):
        """Test for default configuration values to be correct"""
        config = Config()

        assert config.species == "hsa"          # Default species: humans
        assert config.method == "hypergeometric"  # Default method: Hypergeometric check
        assert config.correction == "BH"        # The video is available on the website of the website.
        assert config.pvalue_cutoff == 0.05     # Default p threshold
        assert config.qvalue_cutoff == 0.05     # Default q threshold

    def test_config_validation(self):
        """Verify default configuration (no error)"""
        config = Config()
        errors = config.validate()

        assert len(errors) == 0  # The default configuration should be all valid

    def test_invalid_method(self):
        """Test for invalid method name"""
        config = Config(method="invalid")  # Unintentional setting of invalid method
        errors = config.validate()

        assert any("method" in e.lower() for e in errors)  # Error in methods to be reported

    def test_get_species_config(self):
        """Test for species configuration acquisition (hsa -> Human)"""
        config = Config(species="hsa")
        species_config = config.get_species_config()

        assert species_config.name == "Homo sapiens"   # Academic name
        assert species_config.taxonomy_id == 9606       # ID Classification


class TestFisherExactTest:
    """Fisher's precision test."""

    def test_calculate_pvalue(self):
        """Test p calculation: significant enrichment"""
        fisher = FisherExactTest()

        # Example of test: 50 out of 1, 000 background genes are in this entry and 10 out of 100 input genes are hit
        # Expected: should be significant (p < 0.05)
        pvalue = fisher.calculate_pvalue(
            gene_count=10,           # Enter the number of hits from the gene for the purpose
            background_count=50,     # Quantity of background genes falling under the Article
            gene_total=100,          # Enter the total number of genes
            background_total=1000    # Total number of background genes
        )

        assert 0 <= pvalue <= 1          # The p-value should be within [0, 1]
        assert pvalue < 0.05             # The scene should be quite rich.

    def test_calculate_enrichment(self):
        """Test complete enrichment analysis calculation"""
        fisher = FisherExactTest()

        # Construct Test Data
        gene_set = {"BRCA1", "TP53", "EGFR", "MYC"}                                    # Enter a gene set
        background_set = {"BRCA1", "TP53", "EGFR", "MYC", "KRAS", "PIK3CA", "PTEN", "RB1"}  # Background gene set
        term_genes = {"BRCA1", "TP53", "EGFR", "KRAS", "PIK3CA"}  # Genes annotated to the term.

        result = fisher.calculate_enrichment(
            gene_set=gene_set,
            background_set=background_set,
            term_genes=term_genes,
            term_name="Test Term",
            term_id="TEST:001",
            database="TEST"
        )

        assert result is not None              # Should return result
        assert result.term_id == "TEST:001"    # Entry ID is correct
        assert result.gene_count == 3  # BRCA1, TP53, and EGFR are query hits.
        assert result.rich_factor > 1          # The rich factor should be greater than 1 (meaning rich)


class TestHypergeometricTest:
    """Super-geometric test"""

    def test_calculate_pvalue(self):
        """Test for extra geometric p-value calculation"""
        hyper = HypergeometricTest()

        pvalue = hyper.calculate_pvalue(
            gene_count=10,
            background_count=50,
            gene_total=100,
            background_total=1000
        )

        assert 0 <= pvalue <= 1  # p is within legal range.

    def test_comparison_with_fisher(self):
        """Compared with hypergeometric testsFisherConsistency of results precisely tested

Mathal, Hypergeometric check withFisherThe exact test is equal.,
So for the same parameter, Both.pIt should be very close..
        """
        fisher = FisherExactTest()
        hyper = HypergeometricTest()

        params = {
            "gene_count": 10,
            "background_count": 50,
            "gene_total": 100,
            "background_total": 1000
        }

        p_fisher = fisher.calculate_pvalue(**params)
        p_hyper = hyper.calculate_pvalue(**params)

        # The p-value for both methods should be very close (minimal float error allowed)
        assert abs(p_fisher - p_hyper) < 0.01


class TestEnrichmentResult:
    """Data class testing for enrichment analysis"""

    def test_to_dict(self):
        """Test result object to dictionaries conversion"""
        result = EnrichmentResult(
            term_id="GO:0008150",            # GO entry ID
            term_name="biological_process",  # Entry Name
            database="GO",                    # Database sources
            pvalue=0.001,                     # Original p value
            adjusted_pvalue=0.01,             # Correct p value after correction
            gene_count=10,                    # The number of genes.
            background_count=100,             # Background genes
            expected_count=5.0,               # The number of genes expected.
            rich_factor=2.0,                  # Futrim.
            gene_list=["BRCA1", "TP53"],      # Hit Gene List
            gene_ratio="10/100",              # Query-gene ratio.
            background_ratio="100/1000"       # Background ratio
        )

        d = result.to_dict()

        assert d["Term_ID"] == "GO:0008150"  # Dictionary name should be in upper case underlined format
        assert d["P_Value"] == 0.001
        assert d["Gene_Count"] == 10


class TestEnrichmentAnalyzer:
    """The Eutra Analysis Engine."""

    @pytest.fixture
    def config(self):
        """Creates a test profile object"""
        return Config(
            species="hsa",
            databases=["GO"],
            method="hypergeometric",
            qvalue_cutoff=0.05
        )

    @pytest.fixture
    def analyzer(self, config):
        """Create an example of test analyser"""
        return EnrichmentAnalyzer(config)

    def test_load_gene_list(self, analyzer, tmp_path):
        """Test gene-list file loading"""
        # Create a temporary gene-list file with one identifier per row.
        gene_file = tmp_path / "genes.txt"
        gene_file.write_text("BRCA1\nTP53\nEGFR\n")

        genes = analyzer.load_gene_list(str(gene_file))

        assert len(genes) == 3        # There's three genes to be loaded.
        assert "BRCA1" in genes       # The name should be loaded correctly.

    def test_adjust_pvalues(self, analyzer):
        """Test multiple testspValue Correction (BHMethodology)

        BHCorrect, CorrectpValue >= OriginalpValue.
        """
        # Construct 5 test results, p increment
        results = [
            EnrichmentResult(
                term_id=f"TERM{i}",
                term_name=f"Term {i}",
                database="TEST",
                pvalue=0.01 * (i + 1),      # p value: 0.01, 0.02, 0.03, 0.04, 0.05
                adjusted_pvalue=0.01 * (i + 1),
                gene_count=5,
                background_count=50,
                expected_count=2.5,
                rich_factor=2.0,
                gene_list=[],
                gene_ratio="5/100",
                background_ratio="50/1000"
            )
            for i in range(5)
        ]

        adjusted = analyzer.adjust_pvalues(results, "BH")

        assert len(adjusted) == 5
        # Correct p value should > = original p (minus floating point error allowed)
        for r in adjusted:
            assert r.adjusted_pvalue >= r.pvalue * 0.99

    def test_filter_results(self, analyzer):
        """Test Results Filter

Settings output_all=False, Only keep meets q Visible entry for value threshold.
        """
        analyzer.config.output_all = False
        results = [
            EnrichmentResult(
                term_id="TERM1",
                term_name="Term 1",
                database="TEST",
                pvalue=0.001,
                adjusted_pvalue=0.01,    # Below threshold, should be retained
                gene_count=5,
                background_count=50,
                expected_count=2.5,
                rich_factor=2.0,
                gene_list=[],
                gene_ratio="5/100",
                background_ratio="50/1000"
            ),
            EnrichmentResult(
                term_id="TERM2",
                term_name="Term 2",
                database="TEST",
                pvalue=0.1,
                adjusted_pvalue=0.5,     # Over threshold 0.05, filtered.
                gene_count=5,
                background_count=50,
                expected_count=2.5,
                rich_factor=2.0,
                gene_list=[],
                gene_ratio="5/100",
                background_ratio="50/1000"
            )
        ]

        filtered = analyzer.filter_results(results)

        assert len(filtered) == 1            # Only 1 result retained
        assert filtered[0].term_id == "TERM1"  # It's a TRRM1.

    def test_filter_results_output_all(self, analyzer):
        """Test all results of the default output (withv1Unanimous)

output_all=True Filter only is not satisfactory min_genes entry, Keep All p Value.
        """
        # Default output_all=True
        assert analyzer.config.output_all == True
        results = [
            EnrichmentResult(
                term_id="TERM1",
                term_name="Term 1",
                database="TEST",
                pvalue=0.001,
                adjusted_pvalue=0.01,
                gene_count=5,
                background_count=50,
                expected_count=2.5,
                rich_factor=2.0,
                gene_list=[],
                gene_ratio="5/100",
                background_ratio="50/1000"
            ),
            EnrichmentResult(
                term_id="TERM2",
                term_name="Term 2",
                database="TEST",
                pvalue=0.1,
                adjusted_pvalue=0.5,     # Not significant, but also when output_all=True
                gene_count=5,
                background_count=50,
                expected_count=2.5,
                rich_factor=2.0,
                gene_list=[],
                gene_ratio="5/100",
                background_ratio="50/1000"
            ),
            EnrichmentResult(
                term_id="TERM3",
                term_name="Term 3",
                database="TEST",
                pvalue=0.001,
                adjusted_pvalue=0.01,
                gene_count=1,              # Below Min_genes, should be filtered
                background_count=50,
                expected_count=2.5,
                rich_factor=2.0,
                gene_list=[],
                gene_ratio="1/100",
                background_ratio="50/1000"
            )
        ]

        filtered = analyzer.filter_results(results)

        assert len(filtered) == 2            # Keeps TERM1 and TERM2 and TERM3 is filtered for lack of genes
        assert filtered[0].term_id == "TERM1"
        assert filtered[1].term_id == "TERM2"


class TestEnrichmentMethodEnum:
    """Enumeration test for enrichment analysis"""

    def test_all_methods_exist(self):
        """All expected enrichment analysis methods tested have been defined"""
        methods = [m.value for m in EnrichmentMethod]

        assert "hypergeometric" in methods   # Super-geometric test
        assert "gsea" in methods             # GSEA gene set enrichment analysis
        assert "ssgsea" in methods           # SGSEA single sample GSEA


class TestCorrectionMethodEnum:
    """Multiple Tests Correct Method Enumeration Test"""

    def test_all_corrections_exist(self):
        """Test all expected correction methods are defined"""
        corrections = [c.value for c in CorrectionMethod]

        assert "BH" in corrections           # Benjamini-Hochberg
        assert "BY" in corrections           # Benjamini-Yekutieli
        assert "bonferroni" in corrections   # Bonferroni Correction
        assert "holm" in corrections         # Holm Correction
        assert "none" in corrections         # No correction


if __name__ == "__main__":
    pytest.main([__file__, "-v"])  # Direct Run: python test_enrichment.py
