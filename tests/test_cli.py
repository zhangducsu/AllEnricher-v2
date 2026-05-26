"""
CLI 模块单元测试
"""

import pytest
import sys
import tempfile
import os
from pathlib import Path
from unittest.mock import patch, MagicMock
from io import StringIO

sys.path.insert(0, str(Path(__file__).parent.parent))

from allenricher.cli import create_parser, main, cmd_list, cmd_config


class TestCreateParser:
    """测试参数解析器"""

    def test_no_args(self):
        """测试无参数"""
        parser = create_parser()
        args = parser.parse_args([])
        assert args.command is None

    def test_analyze_required_args(self):
        """测试 analyze 必需参数"""
        parser = create_parser()
        args = parser.parse_args(['analyze', '-i', 'genes.txt'])
        assert args.command == 'analyze'
        assert args.input == 'genes.txt'
        assert args.species == 'hsa'
        assert args.databases == 'GO,KEGG'
        assert args.method == 'fisher'

    def test_analyze_all_args(self):
        """测试 analyze 全部参数"""
        parser = create_parser()
        args = parser.parse_args([
            'analyze', '-i', 'genes.txt', '-s', 'mmu',
            '-d', 'GO,KEGG,Reactome', '-o', 'out/',
            '-m', 'hypergeometric', '-c', 'bonferroni',
            '-p', '0.01', '-q', '0.01', '-n', '5',
            '-j', '4', '--no-plot', '--no-report', '--verbose'
        ])
        assert args.command == 'analyze'
        assert args.species == 'mmu'
        assert args.databases == 'GO,KEGG,Reactome'
        assert args.method == 'hypergeometric'
        assert args.correction == 'bonferroni'
        assert args.pvalue == 0.01
        assert args.qvalue == 0.01
        assert args.min_genes == 5
        assert args.jobs == 4
        assert args.no_plot is True
        assert args.no_report is True
        assert args.verbose is True

    def test_download_args(self):
        """测试 download 参数"""
        parser = create_parser()
        args = parser.parse_args(['download', '-d', 'GO,KEGG', '-s', 'hsa'])
        assert args.command == 'download'
        assert args.databases == 'GO,KEGG'
        assert args.species == 'hsa'

    def test_build_args(self):
        """测试 build 参数"""
        parser = create_parser()
        args = parser.parse_args(['build', '-s', 'hsa', '-t', '9606', '-d', 'GO'])
        assert args.command == 'build'
        assert args.species == 'hsa'
        assert args.taxonomy == 9606
        assert args.databases == 'GO'

    def test_list_species(self):
        """测试 list species"""
        parser = create_parser()
        args = parser.parse_args(['list', 'species'])
        assert args.command == 'list'
        assert args.resource == 'species'

    def test_list_databases(self):
        """测试 list databases"""
        parser = create_parser()
        args = parser.parse_args(['list', 'databases'])
        assert args.command == 'list'
        assert args.resource == 'databases'

    def test_config_args(self):
        """测试 config 参数"""
        parser = create_parser()
        args = parser.parse_args(['config', '-o', 'my_config.yaml'])
        assert args.command == 'config'
        assert args.output == 'my_config.yaml'

    def test_serve_args(self):
        """测试 serve 参数"""
        parser = create_parser()
        args = parser.parse_args(['serve', '--port', '9000'])
        assert args.command == 'serve'
        assert args.port == 9000
        assert args.host == '0.0.0.0'

    def test_invalid_resource(self):
        """测试无效 resource"""
        parser = create_parser()
        with pytest.raises(SystemExit):
            parser.parse_args(['list', 'invalid'])


class TestCmdList:
    """测试 list 命令"""

    def test_list_species(self, capsys):
        """测试列出物种"""
        args = MagicMock(resource='species')
        ret = cmd_list(args)
        assert ret == 0
        captured = capsys.readouterr()
        assert 'hsa' in captured.out
        assert 'mmu' in captured.out

    def test_list_databases(self, capsys):
        """测试列出数据库"""
        args = MagicMock(resource='databases')
        ret = cmd_list(args)
        assert ret == 0
        captured = capsys.readouterr()
        assert 'GO' in captured.out
        assert 'KEGG' in captured.out


class TestCmdConfig:
    """测试 config 命令"""

    def test_generate_config(self, tmp_path):
        """测试生成配置文件"""
        config_path = str(tmp_path / "test_config.yaml")
        args = MagicMock(output=config_path)
        ret = cmd_config(args)
        assert ret == 0
        assert Path(config_path).exists()

        # 验证内容
        with open(config_path) as f:
            content = f.read()
            assert 'species' in content
            assert 'databases' in content


class TestMain:
    """测试主入口"""

    def test_no_command(self, capsys):
        """测试无子命令"""
        with patch('sys.argv', ['allenricher']):
            ret = main()
            assert ret == 0
            captured = capsys.readouterr()
            assert 'usage' in captured.out.lower() or 'AllEnricher' in captured.out

    def test_version(self):
        """测试版本号"""
        with pytest.raises(SystemExit) as exc_info:
            with patch('sys.argv', ['allenricher', '-v']):
                main()
        assert exc_info.value.code == 0


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
