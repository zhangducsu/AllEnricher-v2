"""
CLI Integrated Test - GSEA/ssGSEA/GSVA Related Functions

Test content:
1. CLI Support hypergeometric/gsea/ssgsea/gsva Methodology
2. --expression-matrix Parameters
3. --ranked-genes Parameters
4. Config Organisation GSVA Relevant Fields
5. EnrichmentMethod Organisation GSVA
"""

import pytest
import sys
import tempfile
import os
from pathlib import Path
from unittest.mock import patch, MagicMock
from io import StringIO

sys.path.insert(0, str(Path(__file__).parent.parent))

from allenricher.cli import create_parser
from allenricher.core.config import Config, EnrichmentMethod


class TestCLIMethodChoices:
    """Test CLI method options to support all enrichment analysis methods"""

    def test_cli_method_choices_include_gsva(self):
        """Validation CLI--methodOptions for gsva"""
        parser = create_parser()
        # gsva should be accepted as a valid method value
        args = parser.parse_args(['analyze', '-i', 'genes.txt', '-m', 'gsva'])
        assert args.method == 'gsva'

    def test_cli_method_choices_all_methods(self):
        """Validation CLI--methodOptions support all formal methods"""
        parser = create_parser()
        valid_methods = ['hypergeometric', 'gsea', 'ssgsea', 'gsva']
        for method in valid_methods:
            args = parser.parse_args(['analyze', '-i', 'genes.txt', '-m', method])
            assert args.method == method, f"Method '{method}\"should be supported"

    def test_cli_method_invalid_rejected(self):
        """Validation CLI--methodOption rejects invalid name"""
        parser = create_parser()
        with pytest.raises(SystemExit):
            parser.parse_args(['analyze', '-i', 'genes.txt', '-m', 'invalid_method'])


class TestCLIExpressionMatrixOption:
    """Test--expression-matrixParameters"""

    def test_cli_expression_matrix_option(self):
        """Validate expression-matrixParameters can be correctly parsed"""
        parser = create_parser()
        args = parser.parse_args([
            'analyze', '-i', 'genes.txt',
            '-e', 'expression_matrix.tsv'
        ])
        assert args.expression_matrix == 'expression_matrix.tsv'

    def test_cli_expression_matrix_long_option(self):
        """Validate expression-matrixLong Options Format"""
        parser = create_parser()
        args = parser.parse_args([
            'analyze', '-i', 'genes.txt',
            '--expression-matrix', 'expr.csv'
        ])
        assert args.expression_matrix == 'expr.csv'

    def test_cli_expression_matrix_default_none(self):
        """Validate expression-matrixDefault value is None"""
        parser = create_parser()
        args = parser.parse_args(['analyze', '-i', 'genes.txt'])
        assert args.expression_matrix is None


class TestCLIRankedGenesOption:
    """Test--ranked-genesParameters"""

    def test_cli_ranked_genes_option(self):
        """Validate ranked-genesParameters can be correctly parsed"""
        parser = create_parser()
        args = parser.parse_args([
            'analyze', '-i', 'genes.txt',
            '-r', 'ranked_genes.tsv'
        ])
        assert args.ranked_genes == 'ranked_genes.tsv'

    def test_cli_ranked_genes_long_option(self):
        """Validate ranked-genesLong Options Format"""
        parser = create_parser()
        args = parser.parse_args([
            'analyze', '-i', 'genes.txt',
            '--ranked-genes', 'ranked.txt'
        ])
        assert args.ranked_genes == 'ranked.txt'

    def test_cli_ranked_genes_default_none(self):
        """Validate --ranked-genes default value of None"""
        parser = create_parser()
        args = parser.parse_args(['analyze', '-i', 'genes.txt'])
        assert args.ranked_genes is None

    def test_cli_all_new_options_together(self):
        """Verify that all new parameters can be used simultaneously"""
        parser = create_parser()
        args = parser.parse_args([
            'analyze', '-i', 'genes.txt',
            '-m', 'gsva',
            '-e', 'expression_matrix.tsv',
            '-r', 'ranked_genes.tsv',
            '-s', 'hsa',
            '-d', 'GO,KEGG'
        ])
        assert args.method == 'gsva'
        assert args.expression_matrix == 'expression_matrix.tsv'
        assert args.ranked_genes == 'ranked_genes.tsv'
        assert args.species == 'hsa'
        assert args.databases == 'GO,KEGG'


class TestConfigGSVAFields:
    """Verify Config Class contains GSVA-related fields"""

    def test_config_gsva_method_field(self):
        """Validation Config with gsva_method field with default value 'gsva'"""
        config = Config()
        assert hasattr(config, 'gsva_method')
        assert config.gsva_method == 'gsva'

    def test_config_gsva_kcdf_field(self):
        """Validation Config contains gsva_kcdf field with default value 'Gaussian'"""
        config = Config()
        assert hasattr(config, 'gsva_kcdf')
        assert config.gsva_kcdf == 'Gaussian'

    def test_config_gsva_tau_field(self):
        """Verify that Config contains gsva_tau field with default value of 1.0"""
        config = Config()
        assert hasattr(config, 'gsva_tau')
        assert config.gsva_tau == 1.0

    def test_config_gsva_fields_custom_values(self):
        """Verify that GSVA fields can be customised for a value"""
        config = Config(
            gsva_method="plage",
            gsva_kcdf="Poisson",
            gsva_tau=0.5
        )
        assert config.gsva_method == "plage"
        assert config.gsva_kcdf == "Poisson"
        assert config.gsva_tau == 0.5

    def test_config_gsva_fields_do_not_break_existing_defaults(self):
        """Validation of new GSVA field does not affect default values for existing configurations"""
        config = Config()
        # Check that the default values for the existing field have not been changed
        assert config.method == 'hypergeometric'
        assert config.species == 'hsa'
        assert config.correction == 'BH'
        assert config.pvalue_cutoff == 0.05
        assert config.qvalue_cutoff == 0.05
        assert config.gsea_permutations == 1000
        assert config.gsea_min_size is None
        assert config.gsea_max_size is None


class TestConfigMethodEnum:
    """Verify that Enrichment Method contains GSVA"""

    def test_enrichment_method_enum_has_gsva(self):
        """Verify that Enrichment Method contains GSVA"""
        assert hasattr(EnrichmentMethod, 'GSVA')

    def test_enrichment_method_gsva_value(self):
        """Verify GSVA emulation value 'gsva'"""
        assert EnrichmentMethod.GSVA.value == 'gsva'

    def test_enrichment_method_all_values(self):
        """Verify the correctness of all the enumerations"""
        expected = {
            EnrichmentMethod.HYPERGEOMETRIC: 'hypergeometric',
            EnrichmentMethod.GSEA: 'gsea',
            EnrichmentMethod.SSGSEA: 'ssgsea',
            EnrichmentMethod.GSVA: 'gsva',
        }
        for enum_member, expected_value in expected.items():
            assert enum_member.value == expected_value

    def test_config_validate_accepts_gsva(self):
        """Validation Config. validate() Accept gsva method"""
        config = Config(method='gsva')
        errors = config.validate()
        # gsva should be a legitimate method and should not cause errors in the method.
        method_errors = [e for e in errors if 'method' in e.lower()]
        assert len(method_errors) == 0, f"There should be no methodological error but get: {method_errors}"


class TestGSVAImport:
    """Validation GSVA module can be correctly imported"""

    def test_gsva_module_import(self):
        """Verify allenricher.core.gsva modules to import"""
        from allenricher.core.gsva import GSVA
        assert GSVA is not None

    def test_gsva_class_instantiation(self):
        """Validate GSVA class can be properly executable"""
        from allenricher.core.gsva import GSVA
        gsva = GSVA()
        assert gsva.method == 'gsva'
        assert gsva.kcdf == 'Gaussian'
        assert gsva.tau == 1.0

    def test_gsva_class_custom_params(self):
        """Instantiate GSVA with explicit custom parameters."""
        from allenricher.core.gsva import GSVA
        gsva = GSVA(method='plage', kcdf='Poisson', tau=0.5, min_size=5, max_size=1000)
        assert gsva.method == 'plage'
        assert gsva.kcdf == 'Poisson'
        assert gsva.tau == 0.5
        assert gsva.min_size == 5
        assert gsva.max_size == 1000

    def test_gsva_in_enrichment_method_mapping(self):
        """Validate GSVA in Enrichment Analyzer._get_method"""
        from allenricher.core.enrichment import EnrichmentAnalyzer
        config = Config(method='gsva')
        analyzer = EnrichmentAnalyzer(config)
        method = analyzer.method
        assert method is not None
        assert method.method == 'gsva'
