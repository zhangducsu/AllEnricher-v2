import matplotlib
matplotlib.use("Agg", force=True)
import matplotlib.pyplot as plt
import pandas as pd

from allenricher.visualization.gsva_plots import plot_pathway_heatmap, plot_sample_correlation


def test_pathway_heatmap_disables_clustering_for_single_pathway(tmp_path):
    scores = pd.DataFrame([[0.1, 0.3, -0.2]], index=["GO|E2E|Signal transduction"], columns=["S1", "S2", "S3"])
    output = tmp_path / "heatmap.png"

    fig = plot_pathway_heatmap(scores, output_file=str(output))

    assert fig is not None
    assert output.exists()
    plt.close(fig)


def test_sample_correlation_skips_single_pathway(tmp_path):
    scores = pd.DataFrame([[0.1, 0.3, -0.2]], index=["GO|E2E|Signal transduction"], columns=["S1", "S2", "S3"])
    output = tmp_path / "sample_correlation.png"

    fig = plot_sample_correlation(scores, output_file=str(output))

    assert fig is None
    assert not output.exists()


def test_sample_correlation_uses_sample_columns(tmp_path):
    scores = pd.DataFrame(
        [[0.1, 0.3, -0.2], [0.5, 0.2, -0.4], [-0.1, 0.0, 0.4]],
        index=["Pathway A", "Pathway B", "Pathway C"],
        columns=["S1", "S2", "S3"],
    )
    output = tmp_path / "sample_correlation.png"

    fig = plot_sample_correlation(scores, output_file=str(output))

    assert fig is not None
    assert output.exists()
    sample_labels = ["S1", "S2", "S3"]
    assert any(
        set(tick.get_text() for tick in ax.get_xticklabels()) == set(sample_labels)
        for ax in fig.axes
    )
    plt.close(fig)
