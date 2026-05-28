"""
可视化模块单元测试
"""

import pytest
import sys
from pathlib import Path
import pandas as pd
import numpy as np

# 添加项目路径
sys.path.insert(0, str(Path(__file__).parent.parent))

from allenricher.visualization.plotter import Plotter


class TestPlotter:
    """测试 Plotter 类"""

    @pytest.fixture
    def sample_data(self):
        """创建示例富集分析数据"""
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
        """创建 Plotter 实例"""
        return Plotter(output_dir=str(tmp_path))

    def test_init(self, tmp_path):
        """测试初始化"""
        plotter = Plotter(output_dir=str(tmp_path / "plots"))
        assert plotter.output_dir.exists()
        assert plotter.output_dir.name == "plots"

    def test_plot_barplot(self, plotter, sample_data):
        """测试柱状图生成"""
        output_file = "test_barplot.pdf"
        result = plotter.plot_barplot(sample_data, "GO", output_file, top_n=5)

        # 检查返回路径
        assert result is not None
        assert "test_barplot.pdf" in result

        # 检查文件是否生成（如果 R 可用）
        output_path = Path(result)
        if output_path.exists():
            assert output_path.stat().st_size > 0

    def test_plot_bubble(self, plotter, sample_data):
        """测试气泡图生成"""
        output_file = "test_bubble.pdf"
        result = plotter.plot_bubble(sample_data, output_file, top_n=5)

        # 检查返回路径
        assert result is not None
        assert "test_bubble.pdf" in result

        # 检查文件是否生成（如果 R 可用）
        output_path = Path(result)
        if output_path.exists():
            assert output_path.stat().st_size > 0

    def test_plot_all(self, plotter, sample_data):
        """测试批量生成所有图表"""
        plots = plotter.plot_all(sample_data, "GO", top_n=5)

        # 检查返回的字典
        assert "barplot" in plots
        assert "bubble" in plots

        # 检查路径（现在使用 png 格式作为默认）
        assert "GO_barplot.png" in plots["barplot"]
        assert "GO_bubble.png" in plots["bubble"]

    def test_top_n_filtering(self, plotter, sample_data):
        """测试 top_n 参数过滤"""
        # 使用 top_n=3
        result = plotter.plot_barplot(sample_data, "GO", "test_top3.pdf", top_n=3)
        assert result is not None

    def test_different_databases(self, plotter, sample_data):
        """测试不同数据库"""
        for db in ["GO", "KEGG", "DO", "Reactome", "DisGeNET"]:
            result = plotter.plot_barplot(sample_data, db, f"test_{db}.pdf", top_n=3)
            assert result is not None


class TestPythonPlots:
    """测试 Python 绘图模块"""

    def test_barplot_py_exists(self):
        """测试 barplot.py Python 模块存在"""
        script_path = Path(__file__).parent.parent / "allenricher" / "visualization" / "barplot.py"
        assert script_path.exists()

    def test_bubble_py_exists(self):
        """测试 bubble.py Python 模块存在"""
        script_path = Path(__file__).parent.parent / "allenricher" / "visualization" / "bubble.py"
        assert script_path.exists()

    def test_plot_theme_exists(self):
        """测试 plot_theme.py 风格模块存在"""
        script_path = Path(__file__).parent.parent / "allenricher" / "visualization" / "plot_theme.py"
        assert script_path.exists()


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
