"""
Visualize IC-to-endpoint testing

Test all visualization functions to return correctlyFigureObject
"""

import pytest
import pandas as pd
import numpy as np
from pathlib import Path
import matplotlib.pyplot as plt
import matplotlib.figure

# Import Visualization Module
from allenricher.visualization.gsea_plots import (
    plot_gsea_enrichment, plot_gsea_lollipop
)
from allenricher.visualization.gsva_plots import (
    plot_pathway_heatmap, plot_group_comparison, plot_sample_correlation
)
# Test data path
TEST_DATA_DIR = Path(__file__).parent.parent / "test_data"
RESULTS_DIR = TEST_DATA_DIR / "e2e_results"
OUTPUT_DIR = RESULTS_DIR / "plots"


@pytest.fixture
def gsea_results():
    """Load GSEA results"""
    return pd.read_csv(RESULTS_DIR / "gsea_results.csv")


@pytest.fixture
def ranked_genes():
    """Load Ranked Gene List"""
    return pd.read_csv(TEST_DATA_DIR / "ranked_genes.tsv", sep='\t')


@pytest.fixture
def ssgsea_results():
    """Loading of ssGSEA results"""
    return pd.read_csv(RESULTS_DIR / "ssgsea_results.csv", index_col=0)


@pytest.fixture
def gsva_results():
    """Load GSVA results"""
    return pd.read_csv(TEST_DATA_DIR / "gsva_results.csv", index_col=0)


@pytest.fixture
def gene_sets():
    """Loading of gene sets"""
    gene_sets = {}
    gmt_file = TEST_DATA_DIR / "test_pathways_from_gmt.gmt"
    with open(gmt_file, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 3:
                gene_sets[parts[0]] = set(parts[2:])
    return gene_sets


class TestGSEAPlots:
    """Test GSEA visualization function"""

    def test_plot_gsea_enrichment_returns_figure(self, ranked_genes, gsea_results, gene_sets):
        """Test prot_gsea_enrichment to return to Figure objects"""
        row = gsea_results.iloc[0]
        gene_set = list(gene_sets.values())[0] if gene_sets else set()
        gene_weights = dict(zip(ranked_genes['gene'], ranked_genes['weight']))

        fig = plot_gsea_enrichment(
            ranked_genes=ranked_genes['gene'].tolist()[:500],
            gene_weights=gene_weights,
            gene_set=gene_set,
            es=row['es'],
            nes=row['nes'],
            pvalue=row['pvalue'],
            title="Test Enrichment"
        )

        assert isinstance(fig, matplotlib.figure.Figure)
        plt.close(fig)

class TestGSVAPlots:
    """Test GSVA/ssGSEAVisualisation Functions"""

    def test_plot_pathway_heatmap_returns_figure(self, ssgsea_results):
        """Tests prot_pathway_headmap to return to Figure objects"""
        fig = plot_pathway_heatmap(ssgsea_results, cluster_rows=False, cluster_cols=False)

        assert isinstance(fig, matplotlib.figure.Figure)
        plt.close(fig)

    def test_plot_sample_correlation_returns_figure(self, ssgsea_results):
        """Test plot_sample_colrelation to return to Figure objects"""
        fig = plot_sample_correlation(ssgsea_results, method="pearson")

        assert isinstance(fig, matplotlib.figure.Figure)
        plt.close(fig)

    def test_plot_group_comparison_returns_figure(self, ssgsea_results):
        """Tests prot_group_comparison to return to Figure objects"""
        groups = {
            "Group_A": list(ssgsea_results.columns[:3]),
            "Group_B": list(ssgsea_results.columns[3:])
        }
        fig = plot_group_comparison(
            ssgsea_results,
            groups=groups,
            pathways=list(ssgsea_results.index[:4]),
            plot_type="box"
        )

        assert isinstance(fig, matplotlib.figure.Figure)
        plt.close(fig)

    def test_plot_pathway_heatmap_with_output(self, ssgsea_results, tmp_path):
        """Tests the plot_pathway_heatmap saved files"""
        output_file = tmp_path / "test_heatmap.png"
        fig = plot_pathway_heatmap(
            ssgsea_results,
            cluster_rows=False,
            cluster_cols=False,
            output_file=str(output_file),
            dpi=150
        )

        assert output_file.exists()
        plt.close(fig)

    def test_plot_sample_correlation_spearman(self, ssgsea_results):
        """Tests plot_sample_correlation with pearman method"""
        fig = plot_sample_correlation(ssgsea_results, method="spearman")

        assert isinstance(fig, matplotlib.figure.Figure)
        plt.close(fig)

    def test_plot_group_comparison_violin(self, ssgsea_results):
        """Test plot_group_comparison with virolin type"""
        groups = {
            "Group_A": list(ssgsea_results.columns[:3]),
            "Group_B": list(ssgsea_results.columns[3:])
        }
        fig = plot_group_comparison(
            ssgsea_results,
            groups=groups,
            pathways=list(ssgsea_results.index[:4]),
            plot_type="violin"
        )

        assert isinstance(fig, matplotlib.figure.Figure)
        plt.close(fig)

    def test_plot_group_comparison_bar(self, ssgsea_results):
        """Tests for bar type for plot_group_comparison"""
        groups = {
            "Group_A": list(ssgsea_results.columns[:3]),
            "Group_B": list(ssgsea_results.columns[3:])
        }
        fig = plot_group_comparison(
            ssgsea_results,
            groups=groups,
            pathways=list(ssgsea_results.index[:4]),
            plot_type="bar"
        )

        assert isinstance(fig, matplotlib.figure.Figure)
        plt.close(fig)


class TestIntegration:
    """Integrated Test - Test Full Workstream"""

    def test_full_visualization_pipeline(self, gsea_results, ssgsea_results, gsva_results, gene_sets, tmp_path):
        """Test full visualization pipe"""
        output_dir = tmp_path / "viz_output"
        output_dir.mkdir(exist_ok=True)

        fig1 = plot_gsea_lollipop(
            gsea_results,
            output_file=str(output_dir / "gsea_lollipop.png"),
        )
        fig2 = plot_pathway_heatmap(
            ssgsea_results,
            cluster_rows=False,
            cluster_cols=False,
            output_file=str(output_dir / "ssgsea_heatmap.png"),
            dpi=150
        )
        fig3 = plot_sample_correlation(
            ssgsea_results,
            output_file=str(output_dir / "correlation.png"),
            dpi=150,
        )

        # Validation all files created
        expected_files = [
            "gsea_lollipop.png",
            "ssgsea_heatmap.png",
            "correlation.png"
        ]

        for filename in expected_files:
            assert (output_dir / filename).exists(), f"Documentation{filename}Not created"

        # Clean
        for fig in [fig1, fig2, fig3]:
            plt.close(fig)
