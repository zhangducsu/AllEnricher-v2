"""
CLI build 子命令自定义注释文件参数 单元测试
"""

import pytest
import sys
from pathlib import Path
from unittest.mock import patch, MagicMock

sys.path.insert(0, str(Path(__file__).parent.parent))

from allenricher.cli import create_parser


class TestBuildParserCustomAnnotArgs:
    """测试 build 子命令自定义注释文件参数"""

    def test_build_parser_accepts_go_annot(self):
        """验证 --go-annot 参数被接受"""
        parser = create_parser()
        args = parser.parse_args([
            'build', '-s', 'hsa', '-t', '9606',
            '--go-annot', 'go_annotations.tsv'
        ])
        assert args.go_annot == 'go_annotations.tsv'

    def test_build_parser_accepts_kegg_annot(self):
        """验证 --kegg-annot 参数被接受"""
        parser = create_parser()
        args = parser.parse_args([
            'build', '-s', 'hsa', '-t', '9606',
            '--kegg-annot', 'kegg_annotations.tsv'
        ])
        assert args.kegg_annot == 'kegg_annotations.tsv'

    def test_build_parser_accepts_custom_annot(self):
        """验证 --custom-annot 参数被接受"""
        parser = create_parser()
        args = parser.parse_args([
            'build', '-s', 'hsa', '-t', '9606',
            '--custom-annot', 'custom_annotations.tsv'
        ])
        assert args.custom_annot == 'custom_annotations.tsv'

    def test_build_parser_custom_db_name_default(self):
        """验证 --custom-db-name 默认值为 CUSTOM"""
        parser = create_parser()
        args = parser.parse_args(['build', '-s', 'hsa', '-t', '9606'])
        assert args.custom_db_name == 'CUSTOM'

    def test_build_parser_custom_db_name_explicit(self):
        """验证 --custom-db-name 可显式指定"""
        parser = create_parser()
        args = parser.parse_args([
            'build', '-s', 'hsa', '-t', '9606',
            '--custom-db-name', 'MYDB'
        ])
        assert args.custom_db_name == 'MYDB'

    def test_build_parser_annot_format_choices(self):
        """验证 --annot-format 可选值"""
        parser = create_parser()
        for fmt in ['three_column', 'four_column', 'two_column', 'auto']:
            args = parser.parse_args([
                'build', '-s', 'hsa', '-t', '9606',
                '--annot-format', fmt
            ])
            assert args.annot_format == fmt

    def test_build_parser_annot_format_invalid(self):
        """验证 --annot-format 无效值被拒绝"""
        parser = create_parser()
        with pytest.raises(SystemExit):
            parser.parse_args([
                'build', '-s', 'hsa', '-t', '9606',
                '--annot-format', 'invalid_format'
            ])

    def test_build_parser_annot_format_default(self):
        """验证 --annot-format 默认值为 auto"""
        parser = create_parser()
        args = parser.parse_args(['build', '-s', 'hsa', '-t', '9606'])
        assert args.annot_format == 'auto'

    def test_build_parser_hierarchy_sep_default(self):
        """验证 --hierarchy-sep 默认值为 |"""
        parser = create_parser()
        args = parser.parse_args(['build', '-s', 'hsa', '-t', '9606'])
        assert args.hierarchy_sep == '|'

    def test_build_parser_hierarchy_sep_explicit(self):
        """验证 --hierarchy-sep 可显式指定"""
        parser = create_parser()
        args = parser.parse_args([
            'build', '-s', 'hsa', '-t', '9606',
            '--hierarchy-sep', '/'
        ])
        assert args.hierarchy_sep == '/'

    def test_existing_build_params_unchanged(self):
        """验证现有参数不受影响"""
        parser = create_parser()
        args = parser.parse_args([
            'build', '-s', 'hsa', '-t', '9606', '-d', 'GO,KEGG',
            '--database-dir', '/data/db', '--gene-info', 'gene_info.gz'
        ])
        assert args.command == 'build'
        assert args.species == 'hsa'
        assert args.taxonomy == 9606
        assert args.databases == 'GO,KEGG'
        assert args.database_dir == '/data/db'
        assert args.gene_info == 'gene_info.gz'

    def test_all_custom_annot_params_together(self):
        """验证所有自定义注释参数可同时使用"""
        parser = create_parser()
        args = parser.parse_args([
            'build', '-s', 'mmu', '-t', '10090',
            '--go-annot', 'go.tsv',
            '--kegg-annot', 'kegg.tsv',
            '--custom-annot', 'custom.tsv',
            '--custom-db-name', 'MYDB',
            '--annot-format', 'four_column',
            '--hierarchy-sep', '/'
        ])
        assert args.go_annot == 'go.tsv'
        assert args.kegg_annot == 'kegg.tsv'
        assert args.custom_annot == 'custom.tsv'
        assert args.custom_db_name == 'MYDB'
        assert args.annot_format == 'four_column'
        assert args.hierarchy_sep == '/'


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
