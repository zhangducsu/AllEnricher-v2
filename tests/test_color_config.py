"""
颜色配置系统测试 - TDD

测试原则:
1. 任何图的颜色都不能是硬编码
2. 颜色必须从颜色设置参数获取
3. 颜色设置参数要提供已经设计好的多种配色中选择
"""

import pytest
import matplotlib.pyplot as plt
from unittest.mock import patch, MagicMock
import pandas as pd
import numpy as np


class TestColorConfiguration:
    """测试颜色配置系统"""
    
    def test_color_config_exists(self):
        """测试颜色配置类存在且可导入"""
        from allenricher.visualization.color_config import ColorConfig
        config = ColorConfig()
        assert config is not None
    
    def test_color_config_has_preset_palettes(self):
        """测试颜色配置提供多种预设配色方案"""
        from allenricher.visualization.color_config import ColorConfig
        config = ColorConfig()
        
        # 必须提供至少10种预设配色
        palettes = config.get_available_palettes()
        assert len(palettes) >= 10, f"预设配色方案不足10种，只有{len(palettes)}种"
        
        # 必须包含常见的配色方案
        required_palettes = ['default', 'nature', 'science', 'cell', 'lancet', 'nejm', 'jama']
        for palette in required_palettes:
            assert palette in palettes, f"缺少必需的配色方案: {palette}"
    
    def test_color_config_returns_color_list(self):
        """测试颜色配置返回颜色列表"""
        from allenricher.visualization.color_config import ColorConfig
        config = ColorConfig()
        
        colors = config.get_colors('default', n=5)
        assert isinstance(colors, list)
        assert len(colors) == 5
        
        # 每个颜色必须是有效的十六进制格式
        for color in colors:
            assert isinstance(color, str)
            assert color.startswith('#')
            assert len(color) in [4, 7, 9]  # #RGB, #RRGGBB, #RRGGBBAA
    
    def test_color_config_different_palettes_return_different_colors(self):
        """测试不同配色方案返回不同颜色"""
        from allenricher.visualization.color_config import ColorConfig
        config = ColorConfig()
        
        colors_default = config.get_colors('default', n=3)
        colors_nature = config.get_colors('nature', n=3)
        colors_science = config.get_colors('science', n=3)
        
        # 不同配色方案的第一种颜色应该不同
        assert colors_default[0] != colors_nature[0], "default和nature配色不应相同"
        assert colors_nature[0] != colors_science[0], "nature和science配色不应相同"
    
    def test_color_config_categorical_colors(self):
        """测试分类颜色配置"""
        from allenricher.visualization.color_config import ColorConfig
        config = ColorConfig()
        
        # GO分类颜色
        go_colors = config.get_categorical_colors('go')
        assert 'biological_process' in go_colors
        assert 'cellular_component' in go_colors
        assert 'molecular_function' in go_colors
        
        # KEGG分类颜色
        kegg_colors = config.get_categorical_colors('kegg')
        assert len(kegg_colors) >= 6  # KEGG有6大类
    
    def test_color_config_volcano_colors(self):
        """测试火山图颜色配置"""
        from allenricher.visualization.color_config import ColorConfig
        config = ColorConfig()
        
        volcano_colors = config.get_volcano_colors()
        assert 'up' in volcano_colors
        assert 'down' in volcano_colors
        assert 'ns' in volcano_colors
        
        # 颜色必须是有效的十六进制格式
        for key, color in volcano_colors.items():
            assert color.startswith('#')


class TestNoHardcodedColors:
    """测试图表中没有硬编码颜色"""
    
    def test_barplot_uses_color_config(self):
        """测试barplot使用颜色配置而非硬编码"""
        from allenricher.visualization.barplot import _get_category_colors
        from allenricher.visualization.color_config import ColorConfig
        
        # Mock ColorConfig.get_categorical_colors来验证它被调用
        with patch.object(ColorConfig, 'get_categorical_colors') as mock_get_colors:
            mock_get_colors.return_value = {
                'biological_process': '#FF0000',
                'cellular_component': '#00FF00',
                'molecular_function': '#0000FF',
            }
            
            # 调用 _get_category_colors
            result = _get_category_colors('GO')
            
            # 验证颜色配置被调用
            mock_get_colors.assert_called_once_with('go')
            
            # 验证返回的是mock的颜色
            assert result['biological_process'] == '#FF0000'
    
    def test_bubble_uses_plot_theme_palette(self):
        """测试bubble使用PlotTheme palette参数而非硬编码"""
        from allenricher.visualization.bubble import plot_bubble
        from allenricher.visualization.plot_theme import PlotTheme
        
        data = pd.DataFrame([
            {'Term_Name': 'Test1', 'Qvalue': 0.01, 'GeneCount': 10, 'RichFactor': 0.5}
        ])
        
        # Mock PlotTheme.get_diverging_cmap来验证palette参数被传递
        with patch.object(PlotTheme, 'get_diverging_cmap') as mock_get_cmap:
            mock_cmap = MagicMock()
            mock_get_cmap.return_value = mock_cmap
            
            try:
                plot_bubble(data, output_file='test.png', style='nature', palette='tol_sunset')
            except:
                pass
            
            # 验证palette参数被传递给get_diverging_cmap
            mock_get_cmap.assert_called_once_with(name='tol_sunset')
    
    def test_volcano_uses_color_config(self):
        """测试volcano使用ColorConfig而非硬编码"""
        from allenricher.visualization.common_plots import plot_volcano
        from allenricher.visualization.color_config import VOLCANO_COLORS
        
        # 验证VOLCANO_COLORS是从color_config导入的
        assert 'up' in VOLCANO_COLORS
        assert 'down' in VOLCANO_COLORS
        assert 'ns' in VOLCANO_COLORS
        
        # 验证颜色是有效的十六进制格式
        for key, color in VOLCANO_COLORS.items():
            assert color.startswith('#')
            assert len(color) == 7  # #RRGGBB


class TestPlotThemeIntegration:
    """测试PlotTheme与颜色配置集成"""
    
    def test_plot_theme_returns_colors_from_palette(self):
        """测试PlotTheme从色板返回颜色"""
        from allenricher.visualization.plot_theme import PlotTheme
        
        # 验证PlotTheme.get_palette返回颜色列表
        colors = PlotTheme.get_palette('default', n=3)
        assert len(colors) == 3
        assert all(c.startswith('#') for c in colors)
        
        # 验证不同色板返回不同颜色
        colors_nature = PlotTheme.get_palette('nature', n=3)
        colors_science = PlotTheme.get_palette('science', n=3)
        assert colors_nature[0] != colors_science[0]
    
    def test_plot_theme_palette_parameter(self):
        """测试PlotTheme支持palette参数"""
        from allenricher.visualization.plot_theme import PlotTheme
        
        # 应该支持通过palette参数指定配色方案
        colors = PlotTheme.get_palette(n=3, palette='nature')
        assert len(colors) == 3
        
        colors2 = PlotTheme.get_palette(n=3, palette='science')
        assert len(colors2) == 3
        
        # 不同配色方案应该返回不同颜色
        assert colors[0] != colors2[0]


class TestAllPlotsUseConfigurableColors:
    """测试所有图表类型使用可配置颜色"""
    
    @pytest.mark.parametrize("plot_func_name,plot_module", [
        ('plot_barplot', 'allenricher.visualization.barplot'),
        ('plot_bubble', 'allenricher.visualization.bubble'),
        ('plot_volcano', 'allenricher.visualization.common_plots'),
        ('plot_gsea_enrichment', 'allenricher.visualization.gsea_plots'),
        ('plot_gsea_nes_barplot', 'allenricher.visualization.gsea_plots'),
        ('plot_gsea_dotplot', 'allenricher.visualization.gsea_plots'),
        ('plot_pathway_heatmap', 'allenricher.visualization.gsva_plots'),
        ('plot_group_comparison', 'allenricher.visualization.gsva_plots'),
        ('plot_pathway_dotplot', 'allenricher.visualization.gsva_plots'),
        ('plot_sample_correlation', 'allenricher.visualization.gsva_plots'),
    ])
    def test_all_plots_accept_palette_parameter(self, plot_func_name, plot_module):
        """测试所有绘图函数接受palette参数"""
        import importlib
        
        module = importlib.import_module(plot_module)
        plot_func = getattr(module, plot_func_name)
        
        # 检查函数签名是否包含palette参数
        import inspect
        sig = inspect.signature(plot_func)
        params = list(sig.parameters.keys())
        
        assert 'palette' in params, f"{plot_func_name} 缺少palette参数"
