import gzip
import importlib.util
import json
from types import SimpleNamespace
from pathlib import Path

from allenricher import cli
from allenricher.core.config import Config


def test_analyze_cli_args_override_config(monkeypatch, tmp_path):
    gene_file = tmp_path / "genes.txt"
    gene_file.write_text("TP53\n", encoding="utf-8")
    config_file = tmp_path / "config.yaml"
    config_file.write_text(
        "\n".join([
            f"input_file: {gene_file.as_posix()}",
            "species: mmu",
            "databases: [KEGG]",
            "database_dir: ./from_config",
            f"output_dir: {(tmp_path / 'from_config_output').as_posix()}",
        ]),
        encoding="utf-8",
    )

    captured = {}

    class FakeAnalyzer:
        def __init__(self, config):
            captured["config_species"] = config.species
            captured["config_databases"] = list(config.databases)
            captured["config_database_dir"] = config.database_dir

        def load_gene_list(self, _path):
            return {"TP53"}

        def run_analysis(self, *_args, **_kwargs):
            return {}

    class FakeDatabaseManager:
        active_version = None

        def __init__(self, database_dir, species):
            captured["manager_database_dir"] = database_dir
            captured["manager_species"] = species

        def load_databases(self, databases, version=None):
            captured["loaded_databases"] = list(databases)
            captured["loaded_version"] = version

        def get_background_genes(self):
            return {"TP53", "BRCA1"}

        def get_all_term_data(self):
            return {}

        def get_build_metadata(self):
            return {}

    monkeypatch.setattr(cli, "EnrichmentAnalyzer", FakeAnalyzer)
    monkeypatch.setattr(cli, "DatabaseManager", FakeDatabaseManager)

    args = SimpleNamespace(
        verbose=False,
        config=str(config_file),
        input=str(gene_file),
        species="hsa",
        databases="GO",
        output=str(tmp_path / "cli_output"),
        background=None,
        background_mode="annotated",
        method="hypergeometric",
        correction="BH",
        pvalue=0.05,
        qvalue=0.05,
        min_genes=2,
        jobs=1,
        no_plot=True,
        no_report=True,
        only_significant=False,
        ai=None,
        ai_key=None,
        ai_model=None,
        database_dir=str(tmp_path / "cli_database"),
        use_version=None,
        expression_matrix=None,
        ranked_genes=None,
        gmt=None,
        plot_types=None,
        groups=None,
        plot_format="png",
        plot_dpi=300,
        style="nature",
        palette=None,
        tf_database=None,
        tf_only=False,
        use_r_plots=False,
    )

    assert cli.cmd_analyze(args) == 0
    assert captured["config_species"] == "hsa"
    assert captured["config_databases"] == ["GO"]
    assert captured["config_database_dir"] == str(tmp_path / "cli_database")
    assert captured["manager_species"] == "hsa"
    assert captured["manager_database_dir"] == str(tmp_path / "cli_database")
    assert captured["loaded_databases"] == ["GO"]


def test_config_validate_accepts_local_custom_database(tmp_path):
    gene_file = tmp_path / "genes.txt"
    gene_file.write_text("TP53\n", encoding="utf-8")
    species_dir = tmp_path / "database" / "organism" / "vtest" / "e2e"
    species_dir.mkdir(parents=True)
    with gzip.open(species_dir / "e2e.CUSTOM2gene.tab.gz", "wt", encoding="utf-8") as fh:
        fh.write("Gene\tCUSTOM0001\nTP53\t1\n")

    config = Config(
        input_file=str(gene_file),
        species="e2e",
        databases=["CUSTOM"],
        database_dir=str(tmp_path / "database"),
    )

    assert config.validate() == []


def test_e2e_summary_splits_expected_failure_types(tmp_path):
    runner_path = Path(__file__).parents[1] / "test_e2e_2026" / "run_all_e2e.py"
    spec = importlib.util.spec_from_file_location("run_all_e2e", runner_path)
    runner = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(runner)

    results = [
        {
            "id": "BUG001",
            "title": "real bug",
            "group": "g",
            "status": "EXPECTED_FAIL",
            "expected_failure": True,
            "expected_failure_kind": "",
            "exit_code": 1,
            "output_file_count": 0,
            "duration_seconds": 0.1,
            "failed_checks": [],
        },
        {
            "id": "NEG001",
            "title": "negative case",
            "group": "g",
            "status": "EXPECTED_FAIL",
            "expected_failure": True,
            "expected_failure_kind": "intentional_negative",
            "exit_code": 1,
            "output_file_count": 0,
            "duration_seconds": 0.1,
            "failed_checks": [],
        },
    ]

    runner.write_summary(tmp_path, "local", results, {"coverage_policy": {}})

    summary = json.loads((tmp_path / "E2E_SUMMARY.json").read_text(encoding="utf-8"))
    report = (tmp_path / "E2E_SUMMARY.md").read_text(encoding="utf-8")
    assert summary["expected_failure_counts"] == {
        "ACTIONABLE_EXPECTED_FAIL": 1,
        "INTENTIONAL_NEGATIVE": 1,
    }
    assert "## ACTIONABLE_EXPECTED_FAIL" in report
    assert "## INTENTIONAL_NEGATIVE" in report
