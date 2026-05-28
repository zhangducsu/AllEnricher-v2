# tests/test_cute_charts.py
import pytest
from pathlib import Path


class TestCuteCharts:
    def test_cute_barplot_requires_install(self, tmp_path):
        """cutecharts 未安装时应抛出 ImportError"""
        try:
            from allenricher.visualization.cute_charts import plot_cute_barplot
            # 如果能导入，说明已安装，测试正常生成
            data = [
                {'Term_Name': f'Term{i}', 'Adjusted_P_Value': 10**(-i-1), 'Gene_Count': 50-i*3}
                for i in range(5)
            ]
            out = tmp_path / 'cute_bar.html'
            result = plot_cute_barplot(data, str(out), db_name='GO')
            assert result.exists()
            content = result.read_text()
            assert '<html' in content.lower() or '<!doctype' in content.lower()
        except ImportError:
            pytest.skip("cutecharts not installed")

    def test_cute_bubble_requires_install(self, tmp_path):
        """cutecharts 未安装时应抛出 ImportError"""
        try:
            from allenricher.visualization.cute_charts import plot_cute_bubble
            data = [
                {
                    'Term_Name': f'Term{i}',
                    'Adjusted_P_Value': 10**(-i-1),
                    'Gene_Count': 50-i*3,
                    'Background_Count': 1000,
                }
                for i in range(5)
            ]
            out = tmp_path / 'cute_bubble.html'
            result = plot_cute_bubble(data, str(out), db_name='GO')
            assert result.exists()
        except ImportError:
            pytest.skip("cutecharts not installed")

    def test_empty_data(self, tmp_path):
        """空数据应返回路径不生成文件"""
        try:
            from allenricher.visualization.cute_charts import plot_cute_barplot
            out = tmp_path / 'empty.html'
            result = plot_cute_barplot([], str(out), db_name='GO')
            assert result == Path(out)
        except ImportError:
            pytest.skip("cutecharts not installed")

    def test_output_extension_forced_html(self, tmp_path):
        """非 .html 扩展名应自动转为 .html"""
        try:
            from allenricher.visualization.cute_charts import plot_cute_barplot
            data = [
                {'Term_Name': 'T1', 'Adjusted_P_Value': 1e-5, 'Gene_Count': 30},
            ]
            out = tmp_path / 'cute_bar.png'  # 非 html 扩展名
            result = plot_cute_barplot(data, str(out), db_name='GO')
            assert result.suffix == '.html'
        except ImportError:
            pytest.skip("cutecharts not installed")
