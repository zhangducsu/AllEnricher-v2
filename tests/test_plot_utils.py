import matplotlib
matplotlib.use("Agg", force=True)
import matplotlib.pyplot as plt

from allenricher.visualization.plot_theme import save_figure_dual
from allenricher.visualization.plot_utils import clean_pathway_label, safe_plot_stem, term_figure_size


def test_clean_pathway_label_removes_internal_prefixes():
    assert clean_pathway_label("GO|E2E|Signal transduction") == "Signal transduction"
    assert clean_pathway_label("Go|e2e|signal Transduction") == "Signal transduction"
    assert clean_pathway_label("biological_process:cell Cycle") == "Cell cycle"


def test_safe_plot_stem_is_ascii_and_colon_safe():
    assert safe_plot_stem("GO:0000001") == "GO_0000001"
    assert safe_plot_stem("GO\uf03a0000001") == "GO_0000001"
    assert safe_plot_stem("path/with spaces") == "path_with_spaces"


def test_single_term_figure_size_is_compact():
    width, height = term_figure_size(1)
    assert width == 8.0
    assert height < 3.5


def test_save_figure_dual_writes_requested_svg(tmp_path):
    fig, ax = plt.subplots()
    ax.plot([0, 1], [0, 1])
    output = tmp_path / "plot.svg"

    save_figure_dual(fig, str(output), dpi=120)

    assert output.exists()
    assert (tmp_path / "plot.png").exists()
    assert (tmp_path / "plot.pdf").exists()
    plt.close(fig)
