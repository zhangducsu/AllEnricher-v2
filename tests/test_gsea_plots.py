"""
GSEA Visualization module test
========================

Test coverage: 
- plot_gsea_enrichment generates the correct Figure object
- plot_gsea_lollipop / plot_gsea_ridgeplot to generate the correct Figure object
- Use E2E Test Data (test_data/ranked_genes.tsv and test_data/gene_sets.gmt)
- Verify chart saved to temporary file
"""

import os
import pytest
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from pathlib import Path

# Use non-interactive backends to avoid pop windows
matplotlib.use("Agg")

from allenricher.visualization.gsea_plots import (
    plot_gsea_enrichment,
    plot_gsea_lollipop,
    plot_gsea_multi_enrichment,
    plot_gsea_ridgeplot,
)


# ============================================================
# Test Data Load Tool
# ============================================================

TEST_DATA_DIR = Path(__file__).parent.parent / "test_data"


def _load_ranked_genes() -> tuple:
    """
From test_data/ranked_genes.tsv Loading of the ranked gene list and weight

    Returns:
        (ranked_genes: list[str], gene_weights: dict[str, float])
    """
    tsv_path = TEST_DATA_DIR / "ranked_genes.tsv"
    df = pd.read_csv(tsv_path, sep="\t")
    ranked_genes = df["gene"].tolist()
    gene_weights = dict(zip(df["gene"], df["weight"]))
    return ranked_genes, gene_weights


def _load_gene_sets() -> dict:
    """
From test_data/gene_sets.gmt Loading gene set

    Returns:
        {pathway_name: set(genes)}
    """
    gmt_path = TEST_DATA_DIR / "gene_sets.gmt"
    gene_sets = {}
    with open(gmt_path, "r") as f:
        for line in f:
            parts = line.strip().split("\t")
            name = parts[0]
            genes = set(parts[2:])
            gene_sets[name] = genes
    return gene_sets


def _make_gsea_results_df(n_pathways: int = 10, seed: int = 42) -> pd.DataFrame:
    """
Construct GSEA Outcome DataFrame (For bar and bubble chart testing)

    Returns:
Organisation pathway, nes, pvalue, gene_count Columns DataFrame
    """
    rng = np.random.default_rng(seed)
    pathways = [f"Pathway_{i}" for i in range(n_pathways)]
    nes = rng.uniform(-3.0, 3.0, size=n_pathways)
    pvalue = rng.uniform(1e-10, 0.05, size=n_pathways)
    gene_count = rng.integers(10, 100, size=n_pathways)

    return pd.DataFrame({
        "pathway": pathways,
        "nes": nes,
        "pvalue": pvalue,
        "gene_count": gene_count,
    })


# ============================================================
# Test Class
# ============================================================

class TestPlotGseaEnrichment:
    """Test"""

    def test_returns_figure_object(self):
        """Test returns the matlotlib Figure object"""
        ranked_genes, gene_weights = _load_ranked_genes()
        gene_sets = _load_gene_sets()
        # Take the first gene set.
        pathway_name = list(gene_sets.keys())[0]
        gene_set = gene_sets[pathway_name]

        fig = plot_gsea_enrichment(
            ranked_genes=ranked_genes,
            gene_weights=gene_weights,
            gene_set=gene_set,
            es=0.65,
            nes=1.82,
            pvalue=0.003,
        )

        assert isinstance(fig, matplotlib.figure.Figure)
        plt.close(fig)

    def test_figure_has_three_axes(self):
        """Test chart with three subcharts (three panels)"""
        ranked_genes, gene_weights = _load_ranked_genes()
        gene_sets = _load_gene_sets()
        pathway_name = list(gene_sets.keys())[0]
        gene_set = gene_sets[pathway_name]

        fig = plot_gsea_enrichment(
            ranked_genes=ranked_genes,
            gene_weights=gene_weights,
            gene_set=gene_set,
            es=0.65,
            nes=1.82,
            pvalue=0.003,
        )

        assert len(fig.axes) == 3
        plt.close(fig)

    def test_custom_title(self):
        """Test Custom Titles"""
        ranked_genes, gene_weights = _load_ranked_genes()
        gene_sets = _load_gene_sets()
        pathway_name = list(gene_sets.keys())[0]
        gene_set = gene_sets[pathway_name]

        fig = plot_gsea_enrichment(
            ranked_genes=ranked_genes,
            gene_weights=gene_weights,
            gene_set=gene_set,
            es=0.65,
            nes=1.82,
            pvalue=0.003,
            title="Custom Title Test",
        )

        # suptitle should contain custom titles
        assert fig._suptitle is not None
        assert "Custom Title Test" in fig._suptitle.get_text()
        plt.close(fig)

    def test_save_to_file(self, tmp_path):
        """Test Chart Save to Temporary File"""
        ranked_genes, gene_weights = _load_ranked_genes()
        gene_sets = _load_gene_sets()
        pathway_name = list(gene_sets.keys())[0]
        gene_set = gene_sets[pathway_name]

        output_file = str(tmp_path / "gsea_enrichment.png")
        fig = plot_gsea_enrichment(
            ranked_genes=ranked_genes,
            gene_weights=gene_weights,
            gene_set=gene_set,
            es=0.65,
            nes=1.82,
            pvalue=0.003,
            output_file=output_file,
        )

        assert os.path.exists(output_file)
        assert os.path.getsize(output_file) > 0
        assert fig.number not in plt.get_fignums()
        plt.close(fig)

    def test_empty_gene_set(self):
        """Verify that an empty gene-set collection is handled without error."""
        ranked_genes, gene_weights = _load_ranked_genes()

        fig = plot_gsea_enrichment(
            ranked_genes=ranked_genes[:50],
            gene_weights=gene_weights,
            gene_set=set(),
            es=0.0,
            nes=0.0,
            pvalue=1.0,
        )

        assert isinstance(fig, matplotlib.figure.Figure)
        plt.close(fig)

    def test_all_gene_sets(self):
        """Test all E2E genomes for normal drawing"""
        ranked_genes, gene_weights = _load_ranked_genes()
        gene_sets = _load_gene_sets()

        for pathway_name, gene_set in gene_sets.items():
            fig = plot_gsea_enrichment(
                ranked_genes=ranked_genes,
                gene_weights=gene_weights,
                gene_set=gene_set,
                es=0.5,
                nes=1.5,
                pvalue=0.01,
            )
            assert isinstance(fig, matplotlib.figure.Figure)
            plt.close(fig)


class TestPlotGseaMultiEnrichment:
    def test_returns_one_curve_panel_and_one_hit_panel_per_pathway(self, tmp_path):
        ranked_genes, gene_weights = _load_ranked_genes()
        all_gene_sets = _load_gene_sets()
        selected_ids = list(all_gene_sets)[:2]
        gene_sets = {term_id: all_gene_sets[term_id] for term_id in selected_ids}
        results = pd.DataFrame({
            "Term_ID": selected_ids,
            "Term_Name": selected_ids,
            "NES": [2.0, -2.1],
            "p_value": [0.001, 0.002],
            "Adjusted_P_Value": [0.01, 0.02],
        })
        output = tmp_path / "multi.png"

        fig = plot_gsea_multi_enrichment(
            results,
            selected_ids,
            ranked_genes,
            gene_weights,
            gene_sets,
            output_file=str(output),
        )

        assert len(fig.axes) == len(selected_ids) + 1
        assert tuple(fig.get_size_inches()) == pytest.approx((7.2, 3.47))
        assert output.exists() and output.stat().st_size > 0
        assert fig.number not in plt.get_fignums()
        plt.close(fig)


class TestPlotGseaRidgeplot:
    def test_uses_pathway_gene_scores_and_saves_plot(self, tmp_path):
        ranked_genes, gene_weights = _load_ranked_genes()
        all_gene_sets = _load_gene_sets()
        selected_ids = list(all_gene_sets)[:3]
        results = pd.DataFrame({
            "Term_ID": selected_ids,
            "Term_Name": [f"GSEA|{term_id}" for term_id in selected_ids],
            "NES": [2.1, -1.9, 1.6],
            "p_value": [0.001, 0.003, 0.01],
        })
        output = tmp_path / "ridgeplot.png"

        fig = plot_gsea_ridgeplot(
            results,
            ranked_genes,
            gene_weights,
            all_gene_sets,
            output_file=str(output),
        )

        assert isinstance(fig, matplotlib.figure.Figure)
        assert output.exists() and output.stat().st_size > 0
        assert len(fig.axes[0].collections) >= 3
        assert fig.number not in plt.get_fignums()
        plt.close(fig)


class TestPlotGseaLollipop:
    """Tests"""

    def test_returns_figure_object_with_nes_fallback(self):
        """No Enrich Factor should return to NES."""
        df = _make_gsea_results_df()

        fig = plot_gsea_lollipop(df)

        assert isinstance(fig, matplotlib.figure.Figure)
        assert fig.axes[0].collections[-1].get_alpha() == pytest.approx(1.0)
        plt.close(fig)

    def test_save_to_file_with_enrichfactor_alias(self, tmp_path):
        """Rich_Factor should be directly available for lollipops."""
        df = pd.DataFrame({
            "Term_Name": ["Cell Cycle", "DNA Repair", "Immune Response"],
            "Rich_Factor": [4.2, 3.5, 2.8],
            "Adjusted_P_Value": [1e-4, 2e-3, 5e-3],
            "Gene_Count": [22, 18, 14],
        })
        output_file = str(tmp_path / "gsea_lollipop.png")

        fig = plot_gsea_lollipop(df, output_file=output_file)

        assert os.path.exists(output_file)
        assert os.path.getsize(output_file) > 0
        plt.close(fig)

    def test_enrichfactor_top_n_uses_fdr_before_metric(self):
        df = pd.DataFrame({
            "Term_Name": ["Weak extreme", "Strong one", "Strong two"],
            "Rich_Factor": [100.0, 4.0, 3.0],
            "Adjusted_P_Value": [0.9, 1e-4, 2e-4],
            "Gene_Count": [2, 20, 18],
        })

        fig = plot_gsea_lollipop(df, top_n=2)
        labels = {tick.get_text() for tick in fig.axes[0].get_yticklabels()}

        assert labels == {"Strong one", "Strong two"}
        plt.close(fig)

    def test_constant_fdr_colorbar_shows_the_real_single_value(self):
        df = pd.DataFrame({
            "Term_Name": ["Pathway A", "Pathway B", "Pathway C"],
            "Rich_Factor": [4.2, 3.5, 2.8],
            "Adjusted_P_Value": [0.02, 0.02, 0.02],
            "Gene_Count": [22, 18, 14],
        })

        fig = plot_gsea_lollipop(df)
        colorbar = next(axis for axis in fig.axes if axis.get_title() == "FDR")

        assert len(colorbar.get_yticks()) == 1
        assert colorbar.get_yticklabels()[0].get_text() == "2.0e-2"
        plt.close(fig)


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
