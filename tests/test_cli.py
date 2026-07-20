"""
CLI Module Module Test
"""

import pytest
import sys
import tempfile
import os
import io
import pandas as pd
from pathlib import Path
from argparse import Namespace
from unittest.mock import patch, MagicMock
from io import StringIO

sys.path.insert(0, str(Path(__file__).parent.parent))

from allenricher.cli import (
    _cli_option_was_provided,
    _cmd_tf_enrich,
    _load_tf_gene_file,
    _parse_gmt_term_data,
    _run_tf_analysis,
    create_parser,
    main,
    cmd_analyze,
    cmd_list,
    cmd_config,
    cmd_download,
    cmd_serve,
    cmd_check_update,
)


def test_gmt_parser_preserves_ids_names_and_hierarchy(tmp_path):
    gmt = tmp_path / "custom.gmt"
    gmt.write_text(
        "P1\tMetabolism|Carbohydrate|Pathway one\tG1\tG2\n",
        encoding="utf-8",
    )

    terms = _parse_gmt_term_data(str(gmt))

    assert terms["P1"]["name"] == "Metabolism|Carbohydrate|Pathway one"
    assert terms["P1"]["hierarchy"] == "Metabolism|Carbohydrate|Pathway one"
    assert terms["P1"]["genes"] == {"G1", "G2"}


class TestCreateParser:
    """Test Parameter Solver"""

    def test_no_args(self):
        """Test No Parameters"""
        parser = create_parser()
        args = parser.parse_args([])
        assert args.command is None

    def test_analyze_default_args_without_parser_required_input(self):
        """Parser accepts analyze without a universal --input file."""
        parser = create_parser()
        args = parser.parse_args(['analyze'])
        assert args.command == 'analyze'
        assert args.input is None
        assert args.species == 'hsa'
        assert args.databases == 'GO,KEGG'
        assert args.method == 'hypergeometric'
        assert args.min_genes == 3

    def test_analyze_all_args(self):
        """Test all allyze parameters"""
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

    def test_explicit_default_cli_option_is_treated_as_override(self):
        parser = create_parser()
        args = parser.parse_args([
            'analyze', '-i', 'genes.txt', '-s', 'hsa', '--config', 'config.yaml'
        ])
        args._provided_options = {'-s', '--config'}

        assert _cli_option_was_provided(
            args, '-s', '--species', fallback=args.species != 'hsa'
        )

    def test_download_args(self):
        """Test download parameters"""
        parser = create_parser()
        args = parser.parse_args(['download', '-d', 'GO,KEGG', '-s', 'hsa'])
        assert args.command == 'download'
        assert args.databases == 'GO,KEGG'
        assert args.species == 'hsa'

    def test_build_args(self):
        """Test Build Parameters"""
        parser = create_parser()
        args = parser.parse_args(['build', '-s', 'hsa', '-t', '9606', '-d', 'GO'])
        assert args.command == 'build'
        assert args.species == 'hsa'
        assert args.taxonomy == 9606
        assert args.databases == 'GO'

    def test_list_species(self):
        """Test list species"""
        parser = create_parser()
        args = parser.parse_args(['list', 'species'])
        assert args.command == 'list'
        assert args.resource == 'species'

    def test_list_databases(self):
        """Test list datases"""
        parser = create_parser()
        args = parser.parse_args(['list', 'databases'])
        assert args.command == 'list'
        assert args.resource == 'databases'

    def test_config_args(self):
        """Test config parameters"""
        parser = create_parser()
        args = parser.parse_args(['config', '-o', 'my_config.yaml'])
        assert args.command == 'config'
        assert args.output == 'my_config.yaml'

    def test_serve_args(self):
        """Testserv parameters"""
        parser = create_parser()
        args = parser.parse_args(['serve', '--port', '9000'])
        assert args.command == 'serve'
        assert args.port == 9000
        assert args.host == '127.0.0.1'

    def test_invalid_resource(self):
        """Test invalid"""
        parser = create_parser()
        with pytest.raises(SystemExit):
            parser.parse_args(['list', 'invalid'])


class TestCmdAnalyzeMethodInputs:
    """Runtime checks for method-specific analyze inputs."""

    def test_gsea_uses_ranked_genes_without_query_input_or_background(self, tmp_path):
        parser = create_parser()
        args = parser.parse_args([
            'analyze', '-m', 'gsea', '-r', 'ranked.tsv', '-d', 'GO',
            '--background-mode', 'custom', '--no-plot', '--no-report',
            '-o', str(tmp_path),
        ])

        result = pd.DataFrame({
            'pathway': ['GO:0001'],
            'pval': [0.001],
            'padj': [0.01],
            'ES': [0.6],
            'NES': [1.8],
            'size': [2],
            'leadingEdge': ['G1;G2'],
        })

        with patch('allenricher.cli.EnrichmentAnalyzer') as analyzer_cls, \
             patch('allenricher.cli.DatabaseManager') as manager_cls:
            analyzer = analyzer_cls.return_value
            analyzer.load_ranked_gene_list.return_value = [('G1', 2.0), ('G2', 1.0)]
            analyzer.run_analysis.return_value = {'GO': result}

            manager = manager_cls.return_value
            manager.active_version = 'mock'
            manager.database_versions = {}
            manager.get_all_term_data.return_value = {
                'GO': {'GO:0001': {'name': 'Mock pathway', 'genes': {'G1', 'G2'}}}
            }
            manager.get_build_metadata.return_value = {}

            assert cmd_analyze(args) == 0

        analyzer.load_gene_list.assert_not_called()
        analyzer.load_ranked_gene_list.assert_called_once_with('ranked.tsv')
        call_args = analyzer.run_analysis.call_args
        assert call_args.args[0] == set()
        assert call_args.args[1] == set()
        assert call_args.kwargs['ranked_gene_list'] == [('G1', 2.0), ('G2', 1.0)]

class TestCmdList:
    """Test list command"""

    def test_list_species(self, capsys):
        """Tests listed species"""
        args = MagicMock(resource='species')
        ret = cmd_list(args)
        assert ret == 0
        captured = capsys.readouterr()
        assert 'hsa' in captured.out
        assert 'mmu' in captured.out

    def test_list_databases(self, capsys):
        """Test list database"""
        args = MagicMock(resource='databases')
        ret = cmd_list(args)
        assert ret == 0
        captured = capsys.readouterr()
        for database in (
            'GO', 'KEGG', 'Reactome', 'WikiPathways', 'DO', 'DisGeNET',
            'TRRUST', 'ChEA3', 'AnimalTFDB', 'hTFtarget', 'CUSTOM',
        ):
            assert database in captured.out
        assert 'MSigDB' not in captured.out

    def test_tf_database_parser_exposes_all_local_backends(self):
        parser = create_parser()
        for database in ('trrust', 'chea3', 'animaltfdb', 'htftarget'):
            args = parser.parse_args(['tf-enrich', '-i', 'genes.txt', '-d', database])
            assert args.database == database

    def test_species_registry_parser_exposes_disease_and_tf_filters(self):
        args = create_parser().parse_args([
            'list-species', '--disgenet', '--trrust', '--chea3',
            '--animaltfdb', '--htftarget',
        ])
        assert args.disgenet is True
        assert args.trrust is True
        assert args.chea3 is True
        assert args.animaltfdb is True
        assert args.htftarget is True

    def test_tf_filter_defaults_match_method_contract(self):
        args = create_parser().parse_args(['tf-enrich', '-i', 'genes.txt'])
        assert args.background is None
        assert args.tf_library is None
        assert args.tf_tissue is None
        assert args.tf_regulation == 'all'
        assert args.tf_min_size is None
        assert args.tf_max_size is None
        assert args.tf_combine == 'none'

    def test_tf_ora_background_option_and_gene_file_contract(self, tmp_path):
        genes = tmp_path / 'genes.txt'
        genes.write_text('# comment\nG1\nG1\nG2\n', encoding='utf-8')
        assert _load_tf_gene_file(genes) == ['G1', 'G2']

        args = create_parser().parse_args([
            'tf-enrich', '-i', str(genes), '--background', str(genes)
        ])
        assert args.background == str(genes)

        genes.write_text('G1\t2\n', encoding='utf-8')
        with pytest.raises(ValueError, match='TF ORA requires one gene ID per line'):
            _load_tf_gene_file(genes)

    def test_analyze_tf_overlay_rejects_non_ora_method(self, tmp_path):
        genes = tmp_path / 'genes.txt'
        ranks = tmp_path / 'ranks.tsv'
        genes.write_text('G1\nG2\n', encoding='utf-8')
        ranks.write_text('gene\tweight\nG1\t1\nG2\t-1\n', encoding='utf-8')

        argv = [
            'allenricher', 'analyze', '-i', str(genes), '-m', 'gsea',
            '-r', str(ranks), '--tf-database', 'trrust', '--no-plot',
            '--no-report',
        ]
        with patch.object(sys, 'argv', argv):
            assert main() == 1


class TestCmdConfig:
    """Test config command"""

    def test_generate_config(self, tmp_path):
        """Test Generate Profile"""
        config_path = str(tmp_path / "test_config.yaml")
        args = MagicMock(output=config_path)
        ret = cmd_config(args)
        assert ret == 0
        assert Path(config_path).exists()

        # Validation
        with open(config_path) as f:
            content = f.read()
            assert 'species' in content
            assert 'databases' in content


def test_download_forwards_runtime_options():
    args = Namespace(
        databases='do', species='hsa', database_dir='test-db', workers=8,
        no_multi_thread=True, no_verify=True, force=True,
        trrust=False, chea3=False, animaltfdb=False,
    )
    with patch('allenricher.database.downloader.DataDownloader') as downloader:
        downloader.return_value.download_all.return_value = {}
        assert cmd_download(args) == 0
    downloader.assert_called_once_with(
        root_dir='test-db', overwrite=True, max_workers=8,
        use_multi_thread=False, verify_integrity=False,
    )


def test_trrust_download_updates_the_unified_species_registry():
    from allenricher.cli import _cmd_download_trrust

    args = Namespace(database_dir="test-db", force=False)
    with patch("allenricher.database.trrust_fetcher.TRRUSTFetcher") as fetcher_cls, patch(
        "allenricher.database.downloader.DataDownloader"
    ) as downloader_cls:
        fetcher = fetcher_cls.return_value
        fetcher.download_all.return_value = {"Homo sapiens": Path("trrust.tsv")}
        fetcher.get_supported_species_records.return_value = [
            (9606, "Homo sapiens"),
            (10090, "Mus musculus"),
        ]

        assert _cmd_download_trrust(args) == 0

    downloader = downloader_cls.return_value
    downloader.record_database_species.assert_called_once_with(
        "TRRUST", [(9606, "Homo sapiens"), (10090, "Mus musculus")]
    )
    downloader.refresh_supported_species_registry.assert_called_once_with()


def test_chea3_download_updates_the_unified_species_registry():
    from allenricher.cli import _cmd_download_chea3

    args = Namespace(database_dir="test-db", force=False)
    with patch("allenricher.database.chea3_fetcher.ChEA3Fetcher") as fetcher_cls, patch(
        "allenricher.database.downloader.DataDownloader"
    ) as downloader_cls:
        fetcher = fetcher_cls.return_value
        fetcher.download_all_gmt_libraries.return_value = {"ENCODE": Path("encode.gmt")}
        fetcher.get_supported_species_records.return_value = [(9606, "Homo sapiens")]

        assert _cmd_download_chea3(args) == 0

    downloader = downloader_cls.return_value
    downloader.record_database_species.assert_called_once_with(
        "ChEA3", [(9606, "Homo sapiens")]
    )
    downloader.refresh_supported_species_registry.assert_called_once_with()


def test_animaltfdb_download_records_animal_and_human_tf_coverage():
    from allenricher.cli import _cmd_download_animaltfdb

    args = Namespace(database_dir="test-db", force=False, species=None)
    animal_records = [(9606, "Homo sapiens"), (10090, "Mus musculus")]
    with patch("allenricher.database.animaltfdb_fetcher.AnimalTFDBFetcher") as fetcher_cls, patch(
        "allenricher.database.htftarget_fetcher.HTFtargetFetcher"
    ) as htftarget_cls, patch(
        "allenricher.database.downloader.DataDownloader"
    ) as downloader_cls:
        fetcher = fetcher_cls.return_value
        fetcher.get_supported_species_records.return_value = animal_records
        htftarget_cls.get_supported_species_records.return_value = [(9606, "Homo sapiens")]

        assert _cmd_download_animaltfdb(args) == 0

    downloader = downloader_cls.return_value
    calls = downloader.record_database_species.call_args_list
    assert calls[0].args == ("hTFtarget", [(9606, "Homo sapiens")])
    assert calls[1].args == ("AnimalTFDB", animal_records)
    downloader.refresh_supported_species_registry.assert_called_once_with()


def test_serve_forwards_reload():
    with patch('allenricher.api.server.start_api') as start_api:
        assert cmd_serve(Namespace(host='127.0.0.1', port=8765, reload=True)) == 0
    start_api.assert_called_once_with(host='127.0.0.1', port=8765, reload=True)


def test_tf_gsea_uses_user_weights_and_exports_tf_names_with_official_values(tmp_path):
    import pandas as pd

    ranked = tmp_path / "ranked.tsv"
    ranked.write_text("gene\tweight\nG2\t-2\nG1\t3\nG3\t1\n", encoding="utf-8")
    output = tmp_path / "output"
    official = pd.DataFrame({
        "pathway": ["TF1"],
        "pval": [0.01],
        "padj": [0.02],
        "log2err": [0.2],
        "ES": [0.8],
        "NES": [1.9],
        "size": [3],
        "leadingEdge": ["G1;G3"],
    })
    args = Namespace(
        input=str(ranked), species="hsa", database="trrust",
        output=str(output), top_n=20, database_dir=str(tmp_path / "database"),
        online=False, method="gsea", report=False,
    )

    with patch("allenricher.database.manager.DatabaseManager") as manager_cls, \
         patch("allenricher.analysis.tf_enrichment.TFEnrichmentAnalyzer") as analyzer_cls, \
         patch("allenricher.report.visualizer.Visualizer") as visualizer_cls:
        manager_cls.return_value.load_trrust.return_value = {"gene2tf": pd.DataFrame()}
        analyzer_cls.return_value.gsea.return_value = official.copy()
        analyzer_cls.return_value.metadata_frame.return_value = pd.DataFrame({
            "Term_ID": ["TF1"],
            "Term_Name": ["TF one"],
            "TF": ["TF1"],
            "Library": ["TRRUST"],
            "Context": [""],
            "Evidence_Type": ["curated"],
            "Inference_Type": ["direct"],
        })
        visualizer_cls.return_value.plot_tf_enrichment_bar.return_value = MagicMock()

        assert _cmd_tf_enrich(args) == 0

    manager_cls.assert_called_once_with(str(tmp_path / "database"), "hsa")
    analyzer_cls.return_value.gsea.assert_called_once_with(
        ranked_genes=[("G1", 3.0), ("G3", 1.0), ("G2", -2.0)],
        library=None,
        tissue=None,
        regulation='all',
        min_size=15,
        max_size=5000,
    )
    saved = pd.read_csv(output / "tf_enrichment_trrust_hsa.csv")
    assert saved.loc[0, "Term_ID"] == "TF1"
    assert saved.loc[0, "Term_Name"] == "TF one"
    assert saved.loc[0, "TF"] == "TF1"
    pd.testing.assert_frame_equal(saved[official.columns], official)
    assert not (output / "tf_analysis_metadata_trrust_hsa.json").exists()


def test_tf_gsea_rejects_ora_background_before_loading_database(tmp_path):
    ranked = tmp_path / "ranked.tsv"
    background = tmp_path / "background.txt"
    ranked.write_text("gene\tweight\nG1\t2\nG2\t-1\n", encoding="utf-8")
    background.write_text("G1\nG2\n", encoding="utf-8")
    args = Namespace(
        input=str(ranked), species="hsa", database="trrust",
        output=str(tmp_path / "output"), top_n=20, database_dir=None,
        online=False, method="gsea", report=False, background=str(background),
    )

    with patch("allenricher.database.manager.DatabaseManager") as manager_cls:
        assert _cmd_tf_enrich(args) == 1
    manager_cls.assert_not_called()


def test_analyze_tf_overlay_forwards_real_background_to_ora(tmp_path):
    result = pd.DataFrame({"TF": ["TF1"], "Pvalue": [0.01]})
    args = Namespace(
        tf_database="trrust", database_dir=str(tmp_path), tf_library=None,
        tf_tissue=None, tf_regulation="all", tf_min_size=None, tf_max_size=None,
    )

    with patch("allenricher.database.manager.DatabaseManager") as manager_cls, \
         patch("allenricher.analysis.tf_enrichment.TFEnrichmentAnalyzer") as analyzer_cls:
        manager_cls.return_value.load_trrust.return_value = {"gene2tf": pd.DataFrame()}
        analyzer_cls.return_value.ora.return_value = result
        actual = _run_tf_analysis(
            args, ["G1"], "hsa", background_genes={"G1", "G2", "G3"}
        )

    assert actual is not None
    analyzer_cls.return_value.ora.assert_called_once_with(
        gene_set=["G1"],
        library=None,
        tissue=None,
        regulation="all",
        min_size=3,
        max_size=None,
        background_genes={"G1", "G2", "G3"},
    )


def test_tf_online_gsea_is_rejected_before_api_call(tmp_path):
    ranked = tmp_path / "ranked.tsv"
    ranked.write_text("gene\tweight\nG1\t2\nG2\t-1\n", encoding="utf-8")
    args = Namespace(
        input=str(ranked), species="hsa", database="chea3",
        output=str(tmp_path / "output"), top_n=20, database_dir=None,
        online=True, method="gsea", report=False,
    )

    with patch("allenricher.database.chea3_fetcher.ChEA3Fetcher") as fetcher:
        assert _cmd_tf_enrich(args) == 1
    fetcher.assert_not_called()


def test_check_update_is_gbk_safe():
    stream = io.TextIOWrapper(io.BytesIO(), encoding='gbk')
    status = {'go': {'has_update': False, 'local': {'version': 'v1'}, 'remote': {'remote_version': 'v1'}}}
    with patch('allenricher.database.version.RemoteVersionChecker.check_updates', return_value=status), \
         patch('sys.stdout', stream):
        assert cmd_check_update(Namespace(database_dir='test-db', json=False)) == 0
        stream.flush()


class TestMain:
    """Test Main Entry"""

    def test_no_command(self, capsys):
        """Test No Child Command"""
        with patch('sys.argv', ['allenricher']):
            ret = main()
            assert ret == 0
            captured = capsys.readouterr()
            assert 'usage' in captured.out.lower() or 'AllEnricher' in captured.out

    def test_version(self):
        """Test Version Number"""
        with pytest.raises(SystemExit) as exc_info:
            with patch('sys.argv', ['allenricher', '-v']):
                main()
        assert exc_info.value.code == 0


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
