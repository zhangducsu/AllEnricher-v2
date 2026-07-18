"""
Visualization of module testing
"""

import pytest
import sys
from pathlib import Path
import pandas as pd
import numpy as np

# Add Item Path
sys.path.insert(0, str(Path(__file__).parent.parent))

from allenricher.visualization.plotter import Plotter
from allenricher.core.config import Config


class TestPlotter:
    """Test Platter Class"""

    @pytest.fixture
    def sample_data(self):
        """Create an example of enrichment analysis data"""
        return pd.DataFrame({
            'Term_ID': ['GO:0005576', 'GO:0051301', 'GO:0062023', 'GO:0005515', 'GO:0000070'],
            'Term_Name': [
                'cellular_component:extracellular region',
                'biological_process:cell division',
                'cellular_component:collagen-containing extracellular matrix',
                'molecular_function:protein binding',
                'biological_process:mitotic sister chromatid segregation'
            ],
            'Gene_Count': [172, 55, 54, 605, 14],
            'Background_Count': [172, 55, 54, 605, 14],
            'P_Value': [1e-10, 1e-8, 1e-7, 1e-6, 1e-5],
            'Adjusted_P_Value': [1e-9, 1e-7, 1e-6, 1e-5, 1e-4],
            'Rich_Factor': [2.5, 3.2, 2.8, 1.9, 4.1],
            'Expected_Count': [68.8, 17.2, 19.3, 318.4, 3.4],
            'Genes': ['gene1;gene2;gene3', 'gene4;gene5', 'gene6;gene7', 'gene8;gene9', 'gene10']
        })

    @pytest.fixture
    def plotter(self, tmp_path):
        """Create Plutter instance"""
        return Plotter(output_dir=str(tmp_path))

    def test_init(self, tmp_path):
        """Test Initialisation"""
        plotter = Plotter(output_dir=str(tmp_path / "plots"))
        assert plotter.output_dir.exists()
        assert plotter.output_dir.name == "plots"

    def test_plot_barplot(self, plotter, sample_data):
        """Test column diagram generation"""
        output_file = "test_barplot.pdf"
        result = plotter.plot_barplot(sample_data, "GO", output_file, top_n=5)

        # Check Return Path
        assert result is not None
        assert "test_barplot.pdf" in result

        # Check if file is generated (if R is available)
        output_path = Path(result)
        if output_path.exists():
            assert output_path.stat().st_size > 0

    def test_plot_lollipop(self, plotter, sample_data):
        """Test ORA lollipop production"""
        output_file = "test_lollipop.pdf"
        result = plotter.plot_lollipop(sample_data, "GO", output_file, top_n=5)

        assert result is not None
        assert "test_lollipop.pdf" in result

        output_path = Path(result)
        if output_path.exists():
            assert output_path.stat().st_size > 0

    def test_plot_all(self, plotter, sample_data):
        """Test Batch Generate All Charts"""
        plots = plotter.plot_all(sample_data, "GO", top_n=5)

        # Check for returned dictionary
        assert "barplot" in plots
        assert "lollipop" in plots
        assert "bubble" not in plots

        # Check Path (now using png format as default)
        assert "GO_barplot.png" in plots["barplot"]
        assert "GO_lollipop.png" in plots["lollipop"]

    def test_plot_all_honors_config_formats_dpi_and_size(self, tmp_path, sample_data):
        config = Config(
            plot_formats=["svg"], plot_dpi=150,
            plot_width=7.0, plot_height=5.0,
        )
        plots = Plotter(str(tmp_path), config).plot_all(sample_data, "GO", top_n=3)
        assert plots["barplot"].endswith("GO_barplot.svg")
        assert plots["lollipop"].endswith("GO_lollipop.svg")
        assert (tmp_path / "GO_barplot.svg").stat().st_size > 0
        assert (tmp_path / "GO_lollipop.svg").stat().st_size > 0
        assert not (tmp_path / "GO_barplot.png").exists()

    def test_top_n_filtering(self, plotter, sample_data):
        """Test top_n parameter filter"""
        # Use top_n=3
        result = plotter.plot_barplot(sample_data, "GO", "test_top3.pdf", top_n=3)
        assert result is not None

    def test_different_databases(self, plotter, sample_data):
        """Test different databases"""
        for db in ["GO", "KEGG", "DO", "Reactome", "DisGeNET"]:
            result = plotter.plot_barplot(sample_data, db, f"test_{db}.pdf", top_n=3)
            assert result is not None


class TestPythonPlots:
    """Test Python Drawing Module"""

    def test_barplot_py_exists(self):
        """Test Barplot.py Python Module Exists"""
        script_path = Path(__file__).parent.parent / "allenricher" / "visualization" / "barplot.py"
        assert script_path.exists()

    def test_hierarchy_level_selection_ignores_terminal_terms(self):
        from allenricher.visualization.barplot import _select_hierarchy_level

        level, counts = _select_hierarchy_level([
            "Root|Class A|Term 1",
            "Root|Class B|Term 2",
            "Root|Class C|Term 3",
        ])

        assert counts == {0: 1, 1: 3}
        assert level == 1

    def test_barplot_uses_hierarchy_legend_or_significance_colorbar(self):
        import matplotlib.pyplot as plt
        from allenricher.visualization.barplot import plot_barplot

        data = pd.DataFrame({
            "term_id": ["T1", "T2", "T3"],
            "term": ["Term 1", "Term 2", "Term 3"],
            "qvalue": [0.001, 0.01, 0.04],
            "gene_count": [9, 7, 5],
            "rich_factor": [0.5, 0.4, 0.3],
        })
        hierarchy_map = {
            "T1": "Metabolism|Carbohydrate|Term 1",
            "T2": "Disease|Cancer|Term 2",
            "T3": "Metabolism|Lipid|Term 3",
        }

        hierarchical = plot_barplot(data, database="Custom", hierarchy_map=hierarchy_map)
        legend = hierarchical.axes[0].get_legend()
        assert legend is not None
        assert legend.get_title().get_text() == "Hierarchy level 1"
        assert len({patch.get_facecolor() for patch in hierarchical.axes[0].patches}) == 2
        plt.close(hierarchical)

        flat = plot_barplot(data, database="Custom")
        assert flat.axes[0].get_legend() is None
        assert flat.axes[1].get_ylabel() == "-log10(Q-value)"
        assert flat.axes[1].get_position().height <= flat.axes[0].get_position().height * 0.35
        assert len({patch.get_facecolor() for patch in flat.axes[0].patches}) == 3
        plt.close(flat)

    def test_bubble_py_removed(self):
        """The old bubble.py is no longer retained."""
        script_path = Path(__file__).parent.parent / "allenricher" / "visualization" / "bubble.py"
        assert not script_path.exists()

    def test_plot_theme_exists(self):
        """Test plot_theme.py module exists"""
        script_path = Path(__file__).parent.parent / "allenricher" / "visualization" / "plot_theme.py"
        assert script_path.exists()


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
