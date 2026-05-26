"""
CLI 集成测试 - GSEA/ssGSEA/GSVA 相关功能

测试内容：
1. CLI 支持 fisher/hypergeometric/gsea/ssgsea/gsva 方法
2. --expression-matrix 参数
3. --ranked-genes 参数
4. Config 类包含 GSVA 相关字段
5. EnrichmentMethod 枚举包含 GSVA
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
    """测试 CLI 方法选项支持所有富集分析方法"""

    def test_cli_method_choices_include_gsva(self):
        """验证 CLI --method 选项支持 gsva"""
        parser = create_parser()
        # gsva 应该被接受为合法的 method 值
        args = parser.parse_args(['analyze', '-i', 'genes.txt', '-m', 'gsva'])
        assert args.method == 'gsva'

    def test_cli_method_choices_all_methods(self):
        """验证 CLI --method 选项支持所有五种方法"""
        parser = create_parser()
        valid_methods = ['fisher', 'hypergeometric', 'gsea', 'ssgsea', 'gsva']
        for method in valid_methods:
            args = parser.parse_args(['analyze', '-i', 'genes.txt', '-m', method])
            assert args.method == method, f"方法 '{method}' 应该被支持"

    def test_cli_method_invalid_rejected(self):
        """验证 CLI --method 选项拒绝无效的方法名"""
        parser = create_parser()
        with pytest.raises(SystemExit):
            parser.parse_args(['analyze', '-i', 'genes.txt', '-m', 'invalid_method'])


class TestCLIExpressionMatrixOption:
    """测试 --expression-matrix 参数"""

    def test_cli_expression_matrix_option(self):
        """验证 --expression-matrix 参数可以被正确解析"""
        parser = create_parser()
        args = parser.parse_args([
            'analyze', '-i', 'genes.txt',
            '-e', 'expression_matrix.tsv'
        ])
        assert args.expression_matrix == 'expression_matrix.tsv'

    def test_cli_expression_matrix_long_option(self):
        """验证 --expression-matrix 长选项形式"""
        parser = create_parser()
        args = parser.parse_args([
            'analyze', '-i', 'genes.txt',
            '--expression-matrix', 'expr.csv'
        ])
        assert args.expression_matrix == 'expr.csv'

    def test_cli_expression_matrix_default_none(self):
        """验证 --expression-matrix 默认值为 None"""
        parser = create_parser()
        args = parser.parse_args(['analyze', '-i', 'genes.txt'])
        assert args.expression_matrix is None


class TestCLIRankedGenesOption:
    """测试 --ranked-genes 参数"""

    def test_cli_ranked_genes_option(self):
        """验证 --ranked-genes 参数可以被正确解析"""
        parser = create_parser()
        args = parser.parse_args([
            'analyze', '-i', 'genes.txt',
            '-r', 'ranked_genes.tsv'
        ])
        assert args.ranked_genes == 'ranked_genes.tsv'

    def test_cli_ranked_genes_long_option(self):
        """验证 --ranked-genes 长选项形式"""
        parser = create_parser()
        args = parser.parse_args([
            'analyze', '-i', 'genes.txt',
            '--ranked-genes', 'ranked.txt'
        ])
        assert args.ranked_genes == 'ranked.txt'

    def test_cli_ranked_genes_default_none(self):
        """验证 --ranked-genes 默认值为 None"""
        parser = create_parser()
        args = parser.parse_args(['analyze', '-i', 'genes.txt'])
        assert args.ranked_genes is None

    def test_cli_all_new_options_together(self):
        """验证所有新增参数可以同时使用"""
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
    """验证 Config 类包含 GSVA 相关字段"""

    def test_config_gsva_method_field(self):
        """验证 Config 包含 gsva_method 字段且默认值为 'gsva'"""
        config = Config()
        assert hasattr(config, 'gsva_method')
        assert config.gsva_method == 'gsva'

    def test_config_gsva_kcdf_field(self):
        """验证 Config 包含 gsva_kcdf 字段且默认值为 'Gaussian'"""
        config = Config()
        assert hasattr(config, 'gsva_kcdf')
        assert config.gsva_kcdf == 'Gaussian'

    def test_config_gsva_tau_field(self):
        """验证 Config 包含 gsva_tau 字段且默认值为 1.0"""
        config = Config()
        assert hasattr(config, 'gsva_tau')
        assert config.gsva_tau == 1.0

    def test_config_gsva_fields_custom_values(self):
        """验证 GSVA 字段可以被自定义赋值"""
        config = Config(
            gsva_method="plage",
            gsva_kcdf="Poisson",
            gsva_tau=0.5
        )
        assert config.gsva_method == "plage"
        assert config.gsva_kcdf == "Poisson"
        assert config.gsva_tau == 0.5

    def test_config_gsva_fields_do_not_break_existing_defaults(self):
        """验证新增 GSVA 字段不影响现有配置的默认值"""
        config = Config()
        # 检查原有字段默认值未被改变
        assert config.method == 'fisher'
        assert config.species == 'hsa'
        assert config.correction == 'BH'
        assert config.pvalue_cutoff == 0.05
        assert config.qvalue_cutoff == 0.05
        assert config.gsea_permutations == 1000
        assert config.gsea_min_size == 10
        assert config.gsea_max_size == 500


class TestConfigMethodEnum:
    """验证 EnrichmentMethod 枚举包含 GSVA"""

    def test_enrichment_method_enum_has_gsva(self):
        """验证 EnrichmentMethod 枚举包含 GSVA"""
        assert hasattr(EnrichmentMethod, 'GSVA')

    def test_enrichment_method_gsva_value(self):
        """验证 GSVA 枚举值为 'gsva'"""
        assert EnrichmentMethod.GSVA.value == 'gsva'

    def test_enrichment_method_all_values(self):
        """验证所有枚举值的正确性"""
        expected = {
            EnrichmentMethod.FISHER: 'fisher',
            EnrichmentMethod.HYPERGEOMETRIC: 'hypergeometric',
            EnrichmentMethod.GSEA: 'gsea',
            EnrichmentMethod.SSGSEA: 'ssgsea',
            EnrichmentMethod.GSVA: 'gsva',
        }
        for enum_member, expected_value in expected.items():
            assert enum_member.value == expected_value

    def test_config_validate_accepts_gsva(self):
        """验证 Config.validate() 接受 gsva 方法"""
        config = Config(method='gsva')
        errors = config.validate()
        # gsva 应该是合法方法，不应产生方法相关的错误
        method_errors = [e for e in errors if 'method' in e.lower()]
        assert len(method_errors) == 0, f"不应有方法相关错误，但得到: {method_errors}"


class TestGSVAImport:
    """验证 GSVA 模块可以被正确导入"""

    def test_gsva_module_import(self):
        """验证 allenricher.core.gsva 模块可以导入"""
        from allenricher.core.gsva import GSVA
        assert GSVA is not None

    def test_gsva_class_instantiation(self):
        """验证 GSVA 类可以正常实例化"""
        from allenricher.core.gsva import GSVA
        gsva = GSVA()
        assert gsva.method == 'gsva'
        assert gsva.kcdf == 'Gaussian'
        assert gsva.tau == 1.0

    def test_gsva_class_custom_params(self):
        """验证 GSVA 类可以使用自定义参数实例化"""
        from allenricher.core.gsva import GSVA
        gsva = GSVA(method='plage', kcdf='Poisson', tau=0.5, min_size=5, max_size=1000)
        assert gsva.method == 'plage'
        assert gsva.kcdf == 'Poisson'
        assert gsva.tau == 0.5
        assert gsva.min_size == 5
        assert gsva.max_size == 1000

    def test_gsva_in_enrichment_method_mapping(self):
        """验证 GSVA 在 EnrichmentAnalyzer._get_method 方法映射中"""
        from allenricher.core.enrichment import EnrichmentAnalyzer
        config = Config(method='gsva')
        analyzer = EnrichmentAnalyzer(config)
        method = analyzer.method
        assert method is not None
        assert method.method == 'gsva'
