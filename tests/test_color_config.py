"""
Colour Configuration System Test - TDD

Test principles:
1. No color for a picture can be hard coded.
2. Colours must be obtained from colour settings parameters
3. Colour settings parameters provide for multiple colour selections that have been designed
"""

import hashlib

import pytest
import matplotlib.pyplot as plt
from unittest.mock import patch, MagicMock
import pandas as pd
import numpy as np


class TestColorConfiguration:
    """Test colour configuration system"""
    
    def test_color_config_exists(self):
        """Test colour configuration classes exist and can be imported"""
        from allenricher.visualization.color_config import ColorConfig
        config = ColorConfig()
        assert config is not None
    
    def test_color_config_has_preset_palettes(self):
        """Test colour configuration provides multiple preset colour options"""
        from allenricher.visualization.color_config import ColorConfig
        config = ColorConfig()
        
        # At least 10 preset colours must be provided
        palettes = config.get_available_palettes()
        assert len(palettes) >= 10, f"Less than 10 preset colour schemes, only{len(palettes)}H-H-H-H-H-H-H-H-H-H-H-H"
        
        # Must contain common colour schemes
        required_palettes = [
            'tol_bright', 'nature', 'science', 'cell', 'lancet', 'nejm', 'jama',
            'colorbrewer_purd', 'echarts_v4',
        ]
        for palette in required_palettes:
            assert palette in palettes, f"Lack of required colour options: {palette}"
    
    def test_color_config_returns_color_list(self):
        """Test colour profile returns the color list"""
        from allenricher.visualization.color_config import ColorConfig
        config = ColorConfig()
        
        colors = config.get_colors('default', n=5)
        assert isinstance(colors, list)
        assert len(colors) == 5
        
        # Every color must be a valid hex format
        for color in colors:
            assert isinstance(color, str)
            assert color.startswith('#')
            assert len(color) in [4, 7, 9]  # #RGB, #RRGGBB, #RRGGBBAA
    
    def test_color_config_different_palettes_return_different_colors(self):
        """Test different color schemes to return different colours"""
        from allenricher.visualization.color_config import ColorConfig
        config = ColorConfig()
        
        colors_default = config.get_colors('default', n=3)
        colors_nature = config.get_colors('nature', n=3)
        colors_science = config.get_colors('science', n=3)
        
        # The first color of different colour options should be different.
        assert colors_default[0] != colors_nature[0], "The color of the default and nature should not be the same."
        assert colors_nature[0] != colors_science[0], "The color of nature and science should not be the same."
    
    def test_color_config_categorical_colors(self):
        """Test Class Colour Configuration"""
        from allenricher.visualization.color_config import ColorConfig
        config = ColorConfig()
        
        # GO class colour
        go_colors = config.get_categorical_colors('go')
        assert 'biological_process' in go_colors
        assert 'cellular_component' in go_colors
        assert 'molecular_function' in go_colors
        
        # KEGC Class Colour
        kegg_colors = config.get_categorical_colors('kegg')
        assert len(kegg_colors) >= 6  # Kegg has six major categories.
    

class TestNoHardcodedColors:
    """No hard code color in the test chart"""
    
    def test_barplot_uses_color_config(self):
        """Test barplot with colour configuration instead of hard encoding"""
        from allenricher.visualization.barplot import _get_category_colors
        from allenricher.visualization.color_config import ColorConfig
        
        # Mock ClorConfig.get_categorical_colors to verify that it is called
        with patch.object(ColorConfig, 'get_categorical_colors') as mock_get_colors:
            mock_get_colors.return_value = {
                'biological_process': '#FF0000',
                'cellular_component': '#00FF00',
                'molecular_function': '#0000FF',
            }
            
            # Call _get_category_colors
            result = _get_category_colors('GO')
            
            # Validation colour configuration called
            mock_get_colors.assert_called_once_with('go', palette=None)
            
            # Verify that returned is mock's color
            assert result['biological_process'] == '#FF0000'
    
class TestPlotThemeIntegration:
    """Test the plot Theme in combination with the colour configuration"""
    
    def test_plot_theme_returns_colors_from_palette(self):
        """Test ProtTheme to return color from the colorboard"""
        from allenricher.visualization.plot_theme import PlotTheme
        
        # Validate the Prot Theme.get_palette returns the colour list
        colors = PlotTheme.get_palette('default', n=3)
        assert len(colors) == 3
        assert all(c.startswith('#') for c in colors)
        
        # Verify different colour panels to return different colours
        colors_nature = PlotTheme.get_palette('nature', n=3)
        colors_science = PlotTheme.get_palette('science', n=3)
        assert colors_nature[0] != colors_science[0]
    
    def test_plot_theme_palette_parameter(self):
        """Test the Prot Theme support paraette parameters"""
        from allenricher.visualization.plot_theme import PlotTheme
        
        # The colour scheme should be supported by the paraette parameter
        colors = PlotTheme.get_palette(n=3, palette='nature')
        assert len(colors) == 3
        
        colors2 = PlotTheme.get_palette(n=3, palette='science')
        assert len(colors2) == 3
        
        # Different colour schemes should return different colors
        assert colors[0] != colors2[0]


class TestAllPlotsUseConfigurableColors:
    """Test all chart types with configurable colours"""
    
    @pytest.mark.parametrize("plot_func_name,plot_module", [
        ('plot_barplot', 'allenricher.visualization.barplot'),
        ('plot_gsea_enrichment', 'allenricher.visualization.gsea_plots'),
        ('plot_gsea_multi_enrichment', 'allenricher.visualization.gsea_plots'),
        ('plot_gsea_ridgeplot', 'allenricher.visualization.gsea_plots'),
        ('plot_gsea_lollipop', 'allenricher.visualization.gsea_plots'),
        ('plot_pathway_heatmap', 'allenricher.visualization.gsva_plots'),
        ('plot_group_comparison', 'allenricher.visualization.gsva_plots'),
        ('plot_sample_correlation', 'allenricher.visualization.gsva_plots'),
    ])
    def test_all_plots_accept_palette_parameter(self, plot_func_name, plot_module):
        """Test all drawing functions to accept the palette parameter"""
        import importlib
        
        module = importlib.import_module(plot_module)
        plot_func = getattr(module, plot_func_name)
        
        # Check whether the function signature contains a palette parameter
        import inspect
        sig = inspect.signature(plot_func)
        params = list(sig.parameters.keys())
        
        assert 'palette' in params, f"{plot_func_name}Missing paraette parameters"
        assert 'style' in params, f"{plot_func_name}Missing style parameters"


def _figure_digest(fig):
    fig.canvas.draw()
    digest = hashlib.sha256(memoryview(fig.canvas.buffer_rgba())).hexdigest()
    plt.close(fig)
    return digest


class TestColorSelectionBehavior:
    def test_continuous_palettes_use_role_appropriate_light_anchors(self):
        from matplotlib.colors import to_rgb

        from allenricher.visualization.color_config import (
            DIVERGING_WHITE_MIDPOINT,
            DIVERGING_PALETTES,
            MIN_GRADIENT_DISTANCE_FROM_WHITE,
            PaletteSelection,
            SEQUENTIAL_GRADIENT_ANCHORS,
            SEQUENTIAL_PALETTES,
            visible_gradient_colors,
        )
        from allenricher.visualization.plot_theme import PlotTheme

        def distance_from_white(color):
            return sum((1.0 - channel) ** 2 for channel in to_rgb(color)) ** 0.5

        for name in SEQUENTIAL_PALETTES:
            anchors = visible_gradient_colors(
                SEQUENTIAL_PALETTES[name], "sequential", name
            )
            assert anchors == SEQUENTIAL_GRADIENT_ANCHORS[name]
            assert len(anchors) == 2
            assert all(anchor in SEQUENTIAL_PALETTES[name] for anchor in anchors)
            selection = PaletteSelection(sequential=name)
            colors = PlotTheme.get_plot_colors(
                palette=selection, role="sequential", n=257
            )
            assert min(map(distance_from_white, colors)) >= MIN_GRADIENT_DISTANCE_FROM_WHITE

        for name in DIVERGING_PALETTES:
            selection = PaletteSelection(diverging=name)
            colors = PlotTheme.get_plot_colors(
                palette=selection, role="diverging", n=257
            )
            midpoint = colors[len(colors) // 2].upper()
            midpoint_rgb = to_rgb(midpoint)
            midpoint_chroma = max(midpoint_rgb) - min(midpoint_rgb)
            assert (
                midpoint == DIVERGING_WHITE_MIDPOINT
                or distance_from_white(midpoint) <= 0.02
                or midpoint_chroma >= 0.02
            )
            assert "#D9D9D9" not in {color.upper() for color in colors}

    def test_registry_is_shared(self):
        from allenricher.visualization.color_config import (
            CATEGORICAL_PALETTES,
            DIVERGING_PALETTES,
            PALETTES as color_palettes,
            PALETTE_ROLES,
            PUBLIC_PALETTES,
            SEQUENTIAL_PALETTES,
        )
        from allenricher.visualization.plot_theme import PALETTES as theme_palettes

        assert theme_palettes is color_palettes
        assert len(CATEGORICAL_PALETTES) == 15
        assert len(SEQUENTIAL_PALETTES) == 4
        assert len(DIVERGING_PALETTES) == 4
        assert len(color_palettes) == 24
        assert len(PUBLIC_PALETTES) == 23
        assert set(PALETTE_ROLES) == set(PUBLIC_PALETTES)
        assert set(CATEGORICAL_PALETTES).isdisjoint(SEQUENTIAL_PALETTES)
        assert set(CATEGORICAL_PALETTES).isdisjoint(DIVERGING_PALETTES)
        assert set(SEQUENTIAL_PALETTES).isdisjoint(DIVERGING_PALETTES)
        assert "default" not in PUBLIC_PALETTES
        assert color_palettes["default"] == color_palettes["tol_bright"]
        assert "colorbrewer_purd" in SEQUENTIAL_PALETTES
        assert "tol_sunset" in DIVERGING_PALETTES
        assert "echarts_v4" in CATEGORICAL_PALETTES
        assert {
            "go_bp", "go_cc", "go_mf", "kegg_pathway", "gsea",
            "tol_burga", "china_style",
        }.isdisjoint(color_palettes)

    def test_invalid_style_and_palette_are_rejected_before_plotting(self):
        from allenricher.core.config import Config
        from allenricher.visualization.plot_theme import PlotTheme

        assert any("plot_style" in error for error in Config(plot_style="invalid").validate())
        assert any("plot_palette" in error for error in Config(plot_palette="invalid").validate())
        assert any(
            "categorical_palette" in error
            for error in Config(categorical_palette="viridis").validate()
        )
        assert any(
            "sequential_palette" in error
            for error in Config(sequential_palette="nature").validate()
        )
        assert any(
            "diverging_palette" in error
            for error in Config(diverging_palette="colorbrewer_purd").validate()
        )
        with pytest.raises(ValueError, match="Unknown palette"):
            with PlotTheme.context("nature", "invalid"):
                pass
        with pytest.raises(ValueError, match="Unknown sequential palette"):
            PlotTheme.get_sequential_cmap("invalid")
        with pytest.raises(ValueError, match="Unknown diverging palette"):
            PlotTheme.get_diverging_cmap("invalid")

    @pytest.mark.parametrize(
        "palette", ["invalid", "gsea", "tol_burga", "china_style", "default"]
    )
    def test_cli_rejects_unknown_palette(self, palette):
        from allenricher.cli import create_parser

        with pytest.raises(SystemExit):
            create_parser().parse_args([
                "analyze", "-i", "genes.txt", "--palette", palette,
            ])

    @pytest.mark.parametrize("palette", ["colorbrewer_purd", "echarts_v4"])
    def test_cli_accepts_source_named_palettes(self, palette):
        from allenricher.cli import create_parser

        args = create_parser().parse_args([
            "analyze", "-i", "genes.txt", "--palette", palette,
        ])
        assert args.palette == palette

    def test_cli_accepts_role_specific_palettes(self):
        from allenricher.cli import create_parser

        args = create_parser().parse_args([
            "analyze", "-i", "genes.txt",
            "--categorical-palette", "nature",
            "--sequential-palette", "viridis",
            "--diverging-palette", "colorbrewer_brbg",
        ])
        assert args.categorical_palette == "nature"
        assert args.sequential_palette == "viridis"
        assert args.diverging_palette == "colorbrewer_brbg"

    def test_legacy_palette_only_overrides_its_compatible_role(self):
        from allenricher.visualization.color_config import resolve_palette_selection

        sequential = resolve_palette_selection(legacy_palette="colorbrewer_purd")
        assert sequential.categorical == "tol_bright"
        assert sequential.sequential == "colorbrewer_purd"
        assert sequential.diverging == "colorbrewer_rdbu"

        diverging = resolve_palette_selection(legacy_palette="tol_sunset")
        assert diverging.categorical == "tol_bright"
        assert diverging.sequential == "colorbrewer_blues"
        assert diverging.diverging == "tol_sunset"

    def test_explicit_role_palette_has_priority_over_legacy_palette(self):
        from allenricher.visualization.color_config import resolve_palette_selection

        selection = resolve_palette_selection(
            legacy_palette="colorbrewer_purd",
            sequential_palette="viridis",
        )
        assert selection.sequential == "viridis"

    def test_categorical_colors_never_cycle(self, caplog):
        from allenricher.visualization.color_config import categorical_colors

        colors = categorical_colors("tol_high_contrast", 6)
        assert len(colors) == len(set(colors)) == 6
        assert "fallback" in caplog.text
        with pytest.raises(ValueError, match="maximum supported"):
            categorical_colors("tol_high_contrast", 21)

    def test_style_fonts_use_portable_families(self):
        from allenricher.visualization.plot_theme import PRESETS

        assert {preset.font_family for preset in PRESETS.values()} <= {
            "sans-serif", "serif",
        }

    def test_public_styles_are_distinct_and_legacy_aliases_resolve(self):
        from allenricher.visualization.plot_theme import PlotTheme, resolve_style

        assert PlotTheme.available_styles() == ["nature", "science", "presentation"]
        assert resolve_style("cell") == "nature"
        assert resolve_style("omicshare") == "science"
        with pytest.raises(ValueError, match="Unknown style"):
            resolve_style("colorblind")

    def test_figure_style_applies_typography_border_grid_and_linewidth(self):
        from allenricher.visualization.plot_theme import PlotTheme, apply_figure_style

        snapshots = {}
        for style in PlotTheme.available_styles():
            with PlotTheme.context(style):
                fig, ax = plt.subplots()
                line, = ax.plot([0, 1], [0, 1])
                ax.set_title("Style")
                apply_figure_style(fig, style, axes=[ax], grid_axis="x")
                snapshots[style] = {
                    "line_width": line.get_linewidth(),
                    "top_border": ax.spines["top"].get_visible(),
                    "grid": any(item.get_visible() for item in ax.get_xgridlines()),
                    "family": ax.title.get_fontfamily()[0],
                    "title_size": ax.title.get_fontsize(),
                }
                plt.close(fig)

        assert snapshots["nature"]["top_border"] is False
        assert snapshots["nature"]["grid"] is False
        assert snapshots["science"]["top_border"] is True
        assert snapshots["science"]["family"] in {"serif", "DejaVu Serif"}
        assert snapshots["presentation"]["grid"] is True
        assert snapshots["presentation"]["line_width"] > snapshots["science"]["line_width"]
        assert snapshots["presentation"]["title_size"] > snapshots["nature"]["title_size"]

    def test_barplot_default_style_is_valid(self):
        from allenricher.visualization.barplot import plot_barplot

        data = pd.DataFrame({
            "term": ["biological_process|Cell cycle"],
            "qvalue": [0.01],
            "gene_count": [4],
            "rich_factor": [0.25],
        })
        fig = plot_barplot(data)
        assert fig is not None
        plt.close(fig)

    def test_palette_roles_only_change_compatible_plots(self):
        from allenricher.visualization.barplot import plot_barplot
        from allenricher.visualization.color_config import PaletteSelection
        from allenricher.visualization.gsea_plots import plot_gsea_enrichment, plot_gsea_lollipop
        from allenricher.visualization.gsva_plots import plot_group_comparison

        categorical_a = PaletteSelection(categorical="nature")
        categorical_b = PaletteSelection(categorical="science")
        sequential_a = PaletteSelection(sequential="colorbrewer_purd")
        sequential_b = PaletteSelection(sequential="viridis")
        diverging_a = PaletteSelection(diverging="colorbrewer_rdbu")
        diverging_b = PaletteSelection(diverging="tol_sunset")

        ora = pd.DataFrame({
            "term": ["biological_process|Cell cycle", "cellular_component|Nucleus"],
            "qvalue": [0.001, 0.01],
            "gene_count": [8, 6],
            "rich_factor": [0.4, 0.3],
        })
        assert _figure_digest(plot_barplot(ora, palette=categorical_a)) != _figure_digest(
            plot_barplot(ora, palette=categorical_b)
        )
        assert _figure_digest(plot_barplot(ora, palette=sequential_a)) == _figure_digest(
            plot_barplot(ora, palette=sequential_b)
        )

        flat_ora = ora.assign(term=["Cell cycle", "Nucleus"])
        assert _figure_digest(plot_barplot(flat_ora, palette=sequential_a)) != _figure_digest(
            plot_barplot(flat_ora, palette=sequential_b)
        )
        assert _figure_digest(plot_barplot(flat_ora, palette=categorical_a)) == _figure_digest(
            plot_barplot(flat_ora, palette=categorical_b)
        )

        gsea = pd.DataFrame({
            "Term_Name": ["Pathway A", "Pathway B", "Pathway C"],
            "NES": [2.1, -1.8, 1.3],
            "Adjusted_P_Value": [0.001, 0.01, 0.04],
            "Gene_Count": [12, 9, 5],
        })
        assert _figure_digest(plot_gsea_lollipop(gsea, palette=sequential_a)) != _figure_digest(
            plot_gsea_lollipop(gsea, palette=sequential_b)
        )
        assert _figure_digest(plot_gsea_lollipop(gsea, palette=categorical_a)) == _figure_digest(
            plot_gsea_lollipop(gsea, palette=categorical_b)
        )

        ranked = [f"G{index}" for index in range(30)]
        weights = {gene: 2.0 - index * 0.13 for index, gene in enumerate(ranked)}
        gene_set = set(ranked[2:16:2])
        enrichment_args = (ranked, weights, gene_set, 0.4, 1.8, 0.01)
        assert _figure_digest(
            plot_gsea_enrichment(*enrichment_args, palette=diverging_a)
        ) != _figure_digest(
            plot_gsea_enrichment(*enrichment_args, palette=diverging_b)
        )
        assert _figure_digest(
            plot_gsea_enrichment(*enrichment_args, palette=sequential_a)
        ) == _figure_digest(
            plot_gsea_enrichment(*enrichment_args, palette=sequential_b)
        )

        scores = pd.DataFrame(
            [[-1.2, -0.8, 0.7, 1.1], [-0.4, 0.2, 1.3, 0.8], [1.0, 0.5, -0.6, -1.1]],
            index=["P1", "P2", "P3"],
            columns=["S1", "S2", "S3", "S4"],
        )
        groups = {"Control": ["S1", "S2"], "Disease": ["S3", "S4"]}
        assert _figure_digest(
            plot_group_comparison(scores, groups, palette=categorical_a, top_n=2)
        ) != _figure_digest(
            plot_group_comparison(scores, groups, palette=categorical_b, top_n=2)
        )

    def test_style_does_not_select_colors(self):
        from allenricher.visualization.plot_theme import PlotTheme

        assert PlotTheme.get_plot_colors("nature", role="categorical") == (
            PlotTheme.get_plot_colors("presentation", role="categorical")
        )
        assert PlotTheme.get_plot_colors("nature", role="sequential") == (
            PlotTheme.get_plot_colors("science", role="sequential")
        )
        assert PlotTheme.get_plot_colors("nature", role="diverging") == (
            PlotTheme.get_plot_colors("cell", role="diverging")
        )
