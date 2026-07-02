import matplotlib
matplotlib.use("Agg", force=True)
import matplotlib.pyplot as plt
import pandas as pd

from allenricher.visualization.gsea_plots import plot_gsea_dotplot, plot_gsea_nes_barplot


def _single_term_gsea_df():
    return pd.DataFrame({
        "Term_ID": ["GO:0000001"],
        "Term_Name": ["GO|E2E|Signal transduction"],
        "NES": [1.73],
        "P_Value": [0.004],
        "Adjusted_P_Value": [0.02],
        "Gene_Count": [12],
    })


def test_gsea_nes_single_term_uses_clean_compact_label():
    fig = plot_gsea_nes_barplot(_single_term_gsea_df())
    ax = fig.axes[0]

    assert fig.get_size_inches()[1] < 3.5
    assert [tick.get_text() for tick in ax.get_yticklabels()] == ["Signal transduction"]
    assert any("NES=1.73" in text.get_text() for text in ax.texts)
    plt.close(fig)


def test_gsea_dotplot_single_term_uses_clean_compact_label():
    fig = plot_gsea_dotplot(_single_term_gsea_df())
    ax = fig.axes[0]

    assert fig.get_size_inches()[1] < 3.5
    assert [tick.get_text() for tick in ax.get_yticklabels()] == ["Signal transduction"]
    assert any("Genes=12" in text.get_text() for text in ax.texts)
    plt.close(fig)
