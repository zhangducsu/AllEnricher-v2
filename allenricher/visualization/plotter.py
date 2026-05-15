"""
可视化模块 - AllEnricher v2.0
=============================

本模块负责富集分析结果的可视化展示。
- barplot.R: R base原生绘图
- bubble.R: ggplot2绘图 (v1原版)

依赖：
- R 语言环境
"""

import os
import subprocess
import logging
from pathlib import Path
from typing import Dict, List, Optional
import pandas as pd

logger = logging.getLogger(__name__)


class Plotter:
    """可视化绘图器类"""

    def __init__(self, output_dir: str, config=None):
        """初始化绘图器"""
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.config = config
        self.r_script_dir = Path(__file__).parent

    def _run_r_script(self, script_name: str, args: List[str]) -> tuple:
        """运行 R 脚本"""
        script_path = self.r_script_dir / script_name

        if not script_path.exists():
            logger.error(f"R 脚本不存在: {script_path}")
            return 1, "", f"Script not found: {script_path}"

        cmd = ["Rscript", str(script_path)] + args

        try:
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                timeout=300
            )
            return result.returncode, result.stdout, result.stderr
        except subprocess.TimeoutExpired:
            logger.warning("R 脚本执行超时（300秒）")
            return 1, "", "Timeout"
        except Exception as e:
            logger.warning(f"R 脚本执行异常: {e}")
            return 1, "", str(e)

    def _prepare_barplot_data(self, data: pd.DataFrame) -> pd.DataFrame:
        """准备 barplot.R 所需的数据格式
        
        v1 barplot.R 读取的列：
        - rt[,c(3,5,7)]: 第3列(ObservedGeneNum), 第5列(RichFactor), 第7列(adjP)
        - rt[,8]: 第8列(TermName)
        
        期望的列顺序：
        1: (任意), 2: (任意), 3: ObservedGeneNum, 4: (任意), 
        5: RichFactor, 6: (任意), 7: adjP, 8: TermName
        """
        df = data.copy()
        
        # 创建符合 v1 期望的列顺序
        result = pd.DataFrame()
        result['col1'] = df.get('Term_ID', '')  # 第1列：TermID（占位）
        result['col2'] = df.get('Background_Count', 0)  # 第2列：TermGeneNum
        result['col3'] = df.get('Gene_Count', 0)  # 第3列：ObservedGeneNum（实际使用）
        result['col4'] = df.get('Expected_Count', 0)  # 第4列：ExpectedGeneNum
        result['col5'] = df.get('Rich_Factor', 0)  # 第5列：RichFactor（实际使用）
        result['col6'] = df.get('P_Value', 1)  # 第6列：rawP
        result['col7'] = df.get('Adjusted_P_Value', 1)  # 第7列：adjP（实际使用）
        result['col8'] = df.get('Term_Name', '')  # 第8列：TermName（实际使用）
        
        return result

    def _prepare_bubble_data(self, data: pd.DataFrame) -> pd.DataFrame:
        """准备 bubble.R 所需的数据格式
        
        v1 bubble.R 期望的列：RichFactor, TermName, ObservedGeneNum, adjP
        注意：只输出这4列，避免 Genes 列中的分号导致 R read.table 解析失败
        """
        column_mapping = {
            'Rich_Factor': 'RichFactor',
            'Term_Name': 'TermName',
            'Gene_Count': 'ObservedGeneNum',
            'Adjusted_P_Value': 'adjP'
        }
        
        df = data.copy()
        result = pd.DataFrame()
        for old_col, new_col in column_mapping.items():
            if old_col in df.columns:
                result[new_col] = df[old_col]
        
        return result

    def plot_barplot(
        self,
        data: pd.DataFrame,
        database: str,
        output_file: str,
        top_n: int = 20
    ) -> str:
        """生成富集分析柱状图"""
        # 取前 top_n 条数据
        plot_data = data.head(top_n).copy()

        # 准备 v1 格式的数据
        barplot_data = self._prepare_barplot_data(plot_data)

        # 保存为 TSV 文件 (无header，空格分隔，v1格式)
        tsv_file = self.output_dir / f"temp_{database}_barplot.txt"
        barplot_data.to_csv(tsv_file, sep='\t', index=False, header=False)

        # 构建输出路径
        output_path = self.output_dir / output_file

        # 调用 R 脚本
        returncode, stdout, stderr = self._run_r_script(
            "barplot.R",
            [database, str(tsv_file), str(output_path.with_suffix(''))]
        )

        # 清理临时文件
        if tsv_file.exists():
            tsv_file.unlink()

        if returncode != 0:
            logger.warning(f"柱状图生成失败: {stderr}")

        return str(output_path)

    def plot_bubble(
        self,
        data: pd.DataFrame,
        output_file: str,
        top_n: int = 20
    ) -> str:
        """生成富集分析气泡图（按Q值排序，与条形图一致）"""
        # 按校正P值排序后取前 top_n 条数据
        plot_data = data.sort_values(
            by='Adjusted_P_Value', ascending=True
        ).head(top_n).copy()

        # 准备 v1 格式的数据
        bubble_data = self._prepare_bubble_data(plot_data)

        # 保存为 TSV 文件 (有header，v1格式)
        tsv_file = self.output_dir / "temp_bubble.txt"
        bubble_data.to_csv(tsv_file, sep='\t', index=False, header=True)

        # 构建输出路径
        output_path = self.output_dir / output_file

        # 调用 R 脚本
        returncode, stdout, stderr = self._run_r_script(
            "bubble.R",
            [str(tsv_file), str(output_path.with_suffix(''))]
        )

        # 清理临时文件
        if tsv_file.exists():
            tsv_file.unlink()

        if returncode != 0:
            logger.warning(f"气泡图生成失败: {stderr}")

        return str(output_path)

    def plot_all(
        self,
        data: pd.DataFrame,
        database: str,
        top_n: int = 20
    ) -> Dict[str, str]:
        """生成所有标准图表（R脚本同时输出PDF和PNG）"""
        plots = {}

        # 生成柱状图（R脚本会同时生成PDF和PNG）
        bar_file = f"{database}_barplot.pdf"
        self.plot_barplot(data, database, bar_file, top_n)
        plots["barplot"] = str(self.output_dir / bar_file)
        plots["barplot_png"] = str(self.output_dir / bar_file.replace('.pdf', '.png'))

        # 生成气泡图（R脚本会同时生成PDF和PNG）
        bubble_file = f"{database}_bubble.pdf"
        self.plot_bubble(data, bubble_file, top_n)
        plots["bubble"] = str(self.output_dir / bubble_file)
        plots["bubble_png"] = str(self.output_dir / bubble_file.replace('.pdf', '.png'))

        return plots
