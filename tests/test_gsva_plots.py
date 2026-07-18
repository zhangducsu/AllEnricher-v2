"""Tests for ssGSEA and GSVA visualizations.

Coverage includes pathway activity heatmaps, group comparisons, sample correlations,
group annotations, correlation methods, and file output.
"""

import os
from pathlib import Path

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pytest

from allenricher.visualization.gsva_plots import (
    DEFAULT_ACTIVITY_HEATMAP_TOP_N,
    _activity_heatmap_layout,
    _clean_pathway_label,
    _wrap_panel_label,
    plot_group_comparison,
    plot_pathway_heatmap,
    plot_sample_correlation,
    select_activity_heatmap_scores,
)


def test_activity_heatmap_layout_is_data_aware_and_bounded():
    small = _activity_heatmap_layout(5, 6, 24, 1)
    medium = _activity_heatmap_layout(20, 30, 60, 1)
    large = _activity_heatmap_layout(500, 600, 100, 3)

    assert small["heat_width"] < medium["heat_width"] <= large["heat_width"]
    assert small["heat_height"] < medium["heat_height"] <= large["heat_height"]
    assert large["heat_width"] == pytest.approx(9.2)
    assert large["heat_height"] == pytest.approx(9.0)
    assert large["figsize"][0] <= 15.0
    assert large["figsize"][1] <= 12.0
    assert large["sample_stride"] > 1
    assert large["pathway_stride"] > 1


def test_pathway_display_label_removes_hierarchy_without_changing_source_id():
    source = "Metabolism|Carbohydrate metabolism|Glycolysis_Gluconeogenesis"
    assert _clean_pathway_label(source) == "Glycolysis Gluconeogenesis"
    assert source.startswith("Metabolism|")


def test_group_comparison_wraps_long_panel_labels_to_two_lines():
    label = _wrap_panel_label(
        "Age-rage Signaling Pathway In Diabetic Complications And Related Processes"
    )

    assert "\n" in label
    assert len(label.splitlines()) == 2
    assert all(len(line) <= 32 for line in label.splitlines())


def test_activity_heatmap_defaults_to_40_high_variance_rows_without_mutating_source():
    scores = pd.DataFrame(
        {"S1": np.zeros(60), "S2": np.arange(60), "S3": -np.arange(60)},
        index=[f"P{index:02d}" for index in range(60)],
    )

    selected = select_activity_heatmap_scores(scores)

    assert DEFAULT_ACTIVITY_HEATMAP_TOP_N == 40
    assert selected.shape == (40, 3)
    assert "P59" in selected.index
    assert "P00" not in selected.index
    assert scores.shape == (60, 3)


def test_activity_heatmap_prefers_group_difference_when_groups_are_available():
    scores = pd.DataFrame(
        {
            "A1": [0, -20, 0], "A2": [0, 20, 0],
            "B1": [10, -20, 1], "B2": [10, 20, -1],
        },
        index=["group_difference", "within_group_variance", "small_difference"],
    )
    annotation = pd.DataFrame(
        {"Group": ["A", "A", "B", "B"]}, index=scores.columns
    )

    selected = select_activity_heatmap_scores(scores, annotation, top_n=1)

    assert selected.index.tolist() == ["group_difference"]

# Use non-interactive backends to avoid popup windows during testing
matplotlib.use("Agg")


class TestPlotPathwayHeatmap:
    """plot_pathway_headmap test"""

    @pytest.fixture
    def scores_df(self):
        """Create a simulation circuit activity score matrix (5 routes x 6 samples)"""
        np.random.seed(42)
        pathways = [
            "HSA_Cell_Cycle",
            "HSA_DNA_Repair",
            "HSA_PI3K_AKT",
            "HSA_MAPK_Signaling",
            "HSA_Apoptosis",
        ]
        samples = [f"Normal_{i}" for i in range(1, 4)] + [f"Disease_{i}" for i in range(1, 4)]
        data = np.random.randn(len(pathways), len(samples))
        return pd.DataFrame(data, index=pathways, columns=samples)

    @pytest.fixture
    def annotation_col(self):
        """Create Sample Group Comment"""
        samples = [f"Normal_{i}" for i in range(1, 4)] + [f"Disease_{i}" for i in range(1, 4)]
        return pd.DataFrame({
            "Group": ["Normal"] * 3 + ["Disease"] * 3,
        }, index=samples)

    def test_returns_figure(self, scores_df):
        """Test returns the matlotlib Figure object"""
        fig = plot_pathway_heatmap(scores_df)
        assert fig is not None
        assert isinstance(fig, matplotlib.figure.Figure)
        plt.close(fig)

    def test_with_annotation(self, scores_df, annotation_col):
        """Test the thermal chart with group comments"""
        fig = plot_pathway_heatmap(scores_df, annotation_col=annotation_col)
        assert fig is not None
        assert isinstance(fig, matplotlib.figure.Figure)
        assert any(text.get_text() == "Group" for axis in fig.axes for text in axis.texts)
        assert any(text.get_text().startswith("Activity score") for axis in fig.axes for text in axis.texts)
        assert fig.get_size_inches()[0] <= 8.5
        heat_axis = next(
            axis for axis in fig.axes
            if len([tick for tick in axis.get_xticklabels() if tick.get_text()]) == scores_df.shape[1]
        )
        assert all(spine.get_visible() for spine in heat_axis.spines.values())
        assert len(heat_axis.patches) == scores_df.size
        assert heat_axis.patches[0].get_width() == pytest.approx(0.92)
        legend_axis = next(axis for axis in fig.axes if axis.get_label() == "activity_legend")
        category_positions = {
            text.get_text(): text.get_position()
            for text in legend_axis.texts
            if text.get_text() in {"Normal", "Disease"}
        }
        assert category_positions["Normal"][0] == category_positions["Disease"][0]
        assert category_positions["Normal"][1] > category_positions["Disease"][1]
        color_axis = next(axis for axis in legend_axis.child_axes if axis.get_label() == "activity_colorbar")
        assert color_axis.get_position().width < 0.03
        assert color_axis.get_position().x0 > legend_axis.get_position().x0
        assert color_axis.get_yticks().tolist() == pytest.approx([-color_axis.get_ylim()[1], 0, color_axis.get_ylim()[1]])
        fig.canvas.draw()
        renderer = fig.canvas.get_renderer()
        group_title = next(text for text in legend_axis.texts if text.get_text() == "Group")
        color_title = next(text for text in legend_axis.texts if text.get_text().startswith("Activity score"))
        assert min(position[1] for position in category_positions.values()) > color_title.get_position()[1]
        within_group_gap = (
            category_positions["Normal"][1] - category_positions["Disease"][1]
        ) * fig.get_size_inches()[1]
        between_group_gap = (
            min(position[1] for position in category_positions.values())
            - color_title.get_position()[1]
        ) * fig.get_size_inches()[1]
        assert within_group_gap < between_group_gap < 0.5
        assert color_title.get_position()[0] == 0
        assert all(
            not group_title.get_window_extent(renderer).overlaps(patch.get_window_extent(renderer))
            for patch in legend_axis.patches
        )
        color_gap = (
            color_title.get_window_extent(renderer).y0
            - color_axis.get_window_extent(renderer).y1
        ) / fig.dpi
        assert 0 < color_gap < 0.25
        assert all(
            not color_title.get_window_extent(renderer).overlaps(tick.get_window_extent(renderer))
            for tick in color_axis.get_yticklabels()
        )
        plt.close(fig)

    def test_no_clustering(self, scores_df):
        """Tests are not for groups."""
        fig = plot_pathway_heatmap(scores_df, cluster_rows=False, cluster_cols=False)
        assert fig is not None
        plt.close(fig)

    def test_save_to_file(self, scores_df, tmp_path):
        """Could not close temporary folder: %s"""
        output_file = str(tmp_path / "heatmap.png")
        fig = plot_pathway_heatmap(scores_df, output_file=output_file)
        assert os.path.exists(output_file)
        assert os.path.getsize(output_file) > 0
        plt.close(fig)

    def test_save_pdf(self, scores_df, tmp_path):
        """Test Save as PDF"""
        output_file = str(tmp_path / "heatmap.pdf")
        fig = plot_pathway_heatmap(scores_df, output_file=output_file)
        assert os.path.exists(output_file)
        assert os.path.getsize(output_file) > 0
        plt.close(fig)

    def test_custom_figsize(self, scores_df):
        """Test Custom Figsize"""
        fig = plot_pathway_heatmap(scores_df, figsize=(12, 10))
        assert fig is not None
        plt.close(fig)


class TestPlotGroupComparison:
    """plot_group_comparison test"""

    @pytest.fixture
    def scores_df(self):
        """Create a synthetic pathway activity matrix."""
        np.random.seed(42)
        pathways = [f"Pathway_{i}" for i in range(1, 16)]
        samples = [f"Normal_{i}" for i in range(1, 4)] + [f"Disease_{i}" for i in range(1, 4)]
        data = np.random.randn(len(pathways), len(samples))
        return pd.DataFrame(data, index=pathways, columns=samples)

    @pytest.fixture
    def groups(self):
        """Create Group Dictionary"""
        return {
            "Normal": [f"Normal_{i}" for i in range(1, 4)],
            "Disease": [f"Disease_{i}" for i in range(1, 4)],
        }

    def test_box_plot(self, scores_df, groups):
        """Test box chart type"""
        fig = plot_group_comparison(scores_df, groups, plot_type="box")
        assert fig is not None
        assert isinstance(fig, matplotlib.figure.Figure)
        plt.close(fig)

    def test_violin_plot(self, scores_df, groups):
        """Test violin chart type"""
        fig = plot_group_comparison(scores_df, groups, plot_type="violin")
        assert fig is not None
        assert isinstance(fig, matplotlib.figure.Figure)
        plt.close(fig)

    def test_bar_plot(self, scores_df, groups):
        """Test column chart type"""
        fig = plot_group_comparison(scores_df, groups, plot_type="bar")
        assert fig is not None
        assert isinstance(fig, matplotlib.figure.Figure)
        plt.close(fig)

    def test_specific_pathways(self, scores_df, groups):
        """Test the specified route list"""
        pathways = ["Pathway_1", "Pathway_2", "Pathway_3"]
        fig = plot_group_comparison(scores_df, groups, pathways=pathways)
        assert fig is not None
        plt.close(fig)

    def test_invalid_plot_type(self, scores_df, groups):
        """Could not close temporary folder: %s"""
        with pytest.raises(ValueError, match="Unsupported plot_type"):
            plot_group_comparison(scores_df, groups, plot_type="scatter")

    def test_save_to_file(self, scores_df, groups, tmp_path):
        """Test and save the combo charts and the statistical tables"""
        output_file = str(tmp_path / "group_comparison.png")
        fig = plot_group_comparison(scores_df, groups, output_file=output_file)
        statistics_file = tmp_path / "group_comparison.statistics.tsv"
        assert os.path.exists(output_file)
        assert os.path.getsize(output_file) > 0
        statistics = pd.read_csv(statistics_file, sep="\t")
        assert set(statistics["test_type"]) == {"pairwise"}
        assert statistics["adjusted_pvalue"].isna().all()
        labels = [text.get_text() for axis in fig.axes for text in axis.texts]
        assert not any("Kruskal-Wallis, P=" in label or "BH P=" in label for label in labels)
        assert any(label.startswith("Wilcoxon, P=") for label in labels)
        density_lines = [
            line for axis in fig.axes for line in axis.lines
            if len(line.get_xdata()) == 256
        ]
        assert density_lines and any(np.isnan(line.get_xdata()).any() for line in density_lines)
        assert all(spine.get_visible() for axis in fig.axes for spine in axis.spines.values())
        assert set(statistics.columns) == {
            "Pathway", "test_type", "group1", "group2", "raw_pvalue",
            "adjusted_pvalue", "p_adjust_method",
        }
        plt.close(fig)

    def test_three_groups_keep_global_and_adjusted_pvalues(self, scores_df, tmp_path):
        """The three groups and above retain the global tests and multiple corrections."""
        groups = {
            "A": ["Normal_1", "Normal_2"],
            "B": ["Normal_3", "Disease_1"],
            "C": ["Disease_2", "Disease_3"],
        }
        statistics_file = tmp_path / "three_groups.tsv"
        fig = plot_group_comparison(
            scores_df,
            groups,
            pathways=["Pathway_1", "Pathway_2"],
            statistics_file=str(statistics_file),
        )
        statistics = pd.read_csv(statistics_file, sep="\t")
        assert set(statistics["test_type"]) == {"global", "pairwise"}
        assert set(statistics.loc[statistics["test_type"] == "pairwise", "p_adjust_method"]) == {"BH"}
        labels = [text.get_text() for axis in fig.axes for text in axis.texts]
        assert any("Kruskal-Wallis, P=" in label for label in labels)
        assert any("Wilcoxon, BH P=" in label for label in labels)
        plt.close(fig)


class TestPlotSampleCorrelation:
    """plot_sample_colrelation test"""

    @pytest.fixture
    def scores_df(self):
        """Create a synthetic pathway activity matrix."""
        np.random.seed(42)
        pathways = [f"Pathway_{i}" for i in range(1, 11)]
        samples = [f"Normal_{i}" for i in range(1, 4)] + [f"Disease_{i}" for i in range(1, 4)]
        data = np.random.randn(len(pathways), len(samples))
        return pd.DataFrame(data, index=pathways, columns=samples)

    @pytest.fixture
    def annotation_col(self):
        """Create Sample Group Comment"""
        samples = [f"Normal_{i}" for i in range(1, 4)] + [f"Disease_{i}" for i in range(1, 4)]
        return pd.DataFrame({
            "Group": ["Normal"] * 3 + ["Disease"] * 3,
        }, index=samples)

    def test_pearson(self, scores_df):
        """Test Pearson Relevance"""
        fig = plot_sample_correlation(scores_df, method="pearson")
        assert fig is not None
        assert isinstance(fig, matplotlib.figure.Figure)
        labels = {
            tick.get_text()
            for axis in fig.axes
            for tick in [*axis.get_xticklabels(), *axis.get_yticklabels()]
            if tick.get_text()
        }
        assert set(scores_df.columns).issubset(labels)
        assert any(axis.patches for axis in fig.axes)
        assert fig._suptitle.get_text() == "Sample Correlation"
        assert fig._suptitle.get_position()[0] < 0.5
        assert fig._suptitle.get_position()[1] < 0.95
        assert sum(text.get_text() == "1" for axis in fig.axes for text in axis.texts) == scores_df.shape[1]
        heat_axis = next(
            axis for axis in fig.axes
            if sum(text.get_text() == "1" for text in axis.texts) == scores_df.shape[1]
        )
        upper_values = [text for text in heat_axis.texts if text.get_text() != "1"]
        expected_triangle = scores_df.shape[1] * (scores_df.shape[1] - 1) // 2
        assert len(upper_values) == expected_triangle
        assert {text.get_color() for text in heat_axis.texts} == {"#404040"}
        assert len(heat_axis.patches) == expected_triangle
        assert all(spine.get_visible() for spine in heat_axis.spines.values())
        plt.close(fig)

    def test_spearman(self, scores_df):
        """Test Spearman Relevance"""
        fig = plot_sample_correlation(scores_df, method="spearman")
        assert fig is not None
        plt.close(fig)

    def test_with_annotation(self, scores_df, annotation_col):
        """Test with Group Comment"""
        annotation_col["Group"] = ["Control"] * 3 + ["Treatment"] * 3
        fig = plot_sample_correlation(scores_df, annotation_col=annotation_col)
        assert fig is not None
        legend_text = {text.get_text() for axis in fig.axes for text in axis.texts}
        assert {"Group", "Control", "Treatment", "Correlation"}.issubset(legend_text)
        legend_y = {
            text.get_text(): text.get_position()[1]
            for axis in fig.axes for text in axis.texts
            if text.get_text() in {"Group", "Correlation"}
        }
        assert legend_y["Group"] > legend_y["Correlation"]
        legend_axis = next(
            axis for axis in fig.axes
            if axis.get_label() == "sample_correlation_legend"
        )
        assert legend_axis.patches[0].get_facecolor() != legend_axis.patches[1].get_facecolor()
        category_y = {
            text.get_text(): text.get_position()[1]
            for text in legend_axis.texts
            if text.get_text() in {"Control", "Treatment"}
        }
        assert min(category_y.values()) > legend_y["Correlation"]
        within_group_gap = (
            category_y["Control"] - category_y["Treatment"]
        ) * fig.get_size_inches()[1]
        between_group_gap = (
            min(category_y.values()) - legend_y["Correlation"]
        ) * fig.get_size_inches()[1]
        assert within_group_gap < between_group_gap < 0.5
        color_axis = next(
            axis for axis in legend_axis.child_axes
            if axis.get_label() == "sample_correlation_colorbar"
        )
        assert color_axis.get_position().x0 >= legend_axis.get_position().x0
        fig.canvas.draw()
        renderer = fig.canvas.get_renderer()
        group_title = next(text for text in legend_axis.texts if text.get_text() == "Group")
        assert all(
            not group_title.get_window_extent(renderer).overlaps(patch.get_window_extent(renderer))
            for patch in legend_axis.patches
        )
        plt.close(fig)

    def test_invalid_method(self, scores_df):
        """Reject an unsupported correlation method."""
        with pytest.raises(ValueError, match="Unsupported correlation method"):
            plot_sample_correlation(scores_df, method="kendall")

    def test_save_to_file(self, scores_df, tmp_path):
        """Could not close temporary folder: %s"""
        output_file = str(tmp_path / "correlation.png")
        fig = plot_sample_correlation(scores_df, output_file=output_file)
        assert os.path.exists(output_file)
        assert os.path.getsize(output_file) > 0
        plt.close(fig)


class TestE2EWithRealData:
    """End-to-end testing using E2E test data"""

    @pytest.fixture
    def ssgsea_scores(self):
        """Loading of ssGSEA test results"""
        test_data_dir = Path(__file__).parent.parent / "test_data"
        csv_path = test_data_dir / "ssgsea_results.csv"
        if not csv_path.exists():
            pytest.skip(f"E2E Test data not available: {csv_path}")
        df = pd.read_csv(csv_path, index_col=0)
        return df

    @pytest.fixture
    def gsva_scores(self):
        """Load GSVA test results"""
        test_data_dir = Path(__file__).parent.parent / "test_data"
        csv_path = test_data_dir / "gsva_results.csv"
        if not csv_path.exists():
            pytest.skip(f"E2E Test data not available: {csv_path}")
        df = pd.read_csv(csv_path, index_col=0)
        return df

    @pytest.fixture
    def groups(self):
        """E2E Group of Test Data"""
        return {
            "Normal": ["Normal_1", "Normal_2", "Normal_3"],
            "Disease": ["Disease_1", "Disease_2", "Disease_3"],
        }

    @pytest.fixture
    def annotation_col(self):
        """E2E Group Comment for Test Data"""
        samples = [f"Normal_{i}" for i in range(1, 4)] + [f"Disease_{i}" for i in range(1, 4)]
        return pd.DataFrame({
            "Group": ["Normal"] * 3 + ["Disease"] * 3,
        }, index=samples)

    def test_ssgsea_heatmap(self, ssgsea_scores, annotation_col, tmp_path):
        """SSGSEA Data Hot Chart"""
        output_file = str(tmp_path / "ssgsea_heatmap.png")
        fig = plot_pathway_heatmap(
            ssgsea_scores, annotation_col=annotation_col, output_file=output_file
        )
        assert isinstance(fig, matplotlib.figure.Figure)
        assert os.path.exists(output_file)
        plt.close(fig)

    def test_gsva_heatmap(self, gsva_scores, annotation_col, tmp_path):
        """GSVA Data Hot Chart"""
        output_file = str(tmp_path / "gsva_heatmap.png")
        fig = plot_pathway_heatmap(
            gsva_scores, annotation_col=annotation_col, output_file=output_file
        )
        assert isinstance(fig, matplotlib.figure.Figure)
        assert os.path.exists(output_file)
        plt.close(fig)

    def test_ssgsea_group_comparison(self, ssgsea_scores, groups, tmp_path):
        """SGSEA inter-group comparison"""
        output_file = str(tmp_path / "ssgsea_group_comparison.png")
        fig = plot_group_comparison(
            ssgsea_scores, groups, plot_type="box", output_file=output_file
        )
        assert isinstance(fig, matplotlib.figure.Figure)
        assert os.path.exists(output_file)
        plt.close(fig)

    def test_ssgsea_correlation(self, ssgsea_scores, annotation_col, tmp_path):
        """sSGSEA Data Relevance Hotchart"""
        output_file = str(tmp_path / "ssgsea_correlation.png")
        fig = plot_sample_correlation(
            ssgsea_scores, annotation_col=annotation_col, output_file=output_file
        )
        assert isinstance(fig, matplotlib.figure.Figure)
        assert os.path.exists(output_file)
        plt.close(fig)


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
