"""
可视化模块 - AllEnricher v2.0
=============================

本模块负责富集分析结果的可视化展示。核心思路是：在Python端准备数据，
动态生成R脚本，然后通过系统调用 Rscript 执行脚本，利用 R 语言的
ggplot2、igraph、ComplexHeatmap 等包生成高质量图表。

支持的图表类型：
- 柱状图 (Bar Plot)：展示富集条目的 -log10(P-value)
- 气泡图 (Bubble Plot)：以 Rich Factor 为横轴、气泡大小表示基因数量
- 点图 (Dot Plot)：以 Gene Ratio 为横轴、展示基因比例与显著性
- 富集图谱 (Enrichment Map)：基于 Jaccard 相似性的术语关系网络图
- 热图 (Heatmap)：跨数据库的富集结果比较热图
- UpSet 图 (UpSet Plot)：展示基因在不同富集条目间的重叠关系
- 基因-条目网络图 (CNet Plot)：展示基因与富集条目之间的关联网络

依赖：
- R 语言环境（需安装 Rscript）
- R 包：ggplot2, igraph, ComplexHeatmap, UpSetR, reshape2
"""

import os
import subprocess
import logging
from pathlib import Path
from typing import Dict, List, Optional
import pandas as pd
import numpy as np

logger = logging.getLogger(__name__)


class Plotter:
    """
    可视化绘图器类
    ==============

    本类是 AllEnricher 的核心可视化组件，通过生成 R 脚本并调用 ggplot2 等
    R 包来生成各种富集分析图表。

    工作流程：
        1. 接收 Python 端的 pandas DataFrame 数据
        2. 对数据进行预处理（排序、计算派生指标等）
        3. 将数据嵌入到动态生成的 R 脚本中
        4. 将 R 脚本写入临时文件并调用 Rscript 执行
        5. 执行完成后删除临时 R 脚本文件
        6. 返回生成的图表文件路径

    使用示例：
        >>> plotter = Plotter(output_dir="./results/plots")
        >>> plotter.plot_barplot(data, "GO Enrichment", "go_barplot.pdf")

    注意事项：
        - 需要系统已安装 R 语言环境
        - 需要安装相应的 R 包（ggplot2, igraph 等）
        - 输出目录会自动创建
    """
    
    def __init__(self, output_dir: str, config=None):
        """
        初始化绘图器

        参数:
            output_dir (str): 图表输出目录路径，所有生成的图表文件将保存到此目录
            config: 可选的配置对象，用于传递全局配置参数（如自定义颜色方案等）

        属性:
            output_dir (Path): 输出目录的 Path 对象
            config: 全局配置对象
            go_colors (dict): GO 富集分析的分类颜色方案（生物过程/细胞组分/分子功能）
            kegg_colors (dict): KEGG 通路分析的分类颜色方案（代谢/遗传信息处理等）
        """
        # 将输出目录转换为 Path 对象，并自动创建（包括父目录）
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.config = config
        
        # GO（Gene Ontology）分类颜色方案
        # 分别对应三大本体：生物过程(BP)、细胞组分(CC)、分子功能(MF)
        self.go_colors = {
            "biological_process": "#E64B35",    # 生物过程 - 红色系
            "cellular_component": "#4DBBD5",     # 细胞组分 - 蓝色系
            "molecular_function": "#00A087"      # 分子功能 - 绿色系
        }
        
        # KEGG 通路分类颜色方案
        # 分别对应 KEGG 的六大功能分类
        self.kegg_colors = {
            "Metabolism": "#66C2A5",                              # 代谢 - 青绿色
            "Genetic_Information_Processing": "#FC8D62",          # 遗传信息处理 - 橙色
            "Environmental_Information_Processing": "#8DA0CB",    # 环境信息处理 - 灰蓝色
            "Cellular_Processes": "#E78AC3",                      # 细胞过程 - 粉色
            "Organismal_Systems": "#A6D854",                      # 生物体系统 - 黄绿色
            "Human_Diseases": "#FFD92F"                           # 人类疾病 - 黄色
        }
    
    def plot_barplot(
        self,
        data: pd.DataFrame,
        title: str,
        output_file: str,
        top_n: int = 20,
        color: str = "#3498db"
    ) -> str:
        """
        生成富集分析柱状图

        根据富集分析结果生成水平柱状图，纵轴为富集条目名称，横轴为
        -log10(校正后P值)，柱子按显著性从低到高排列。

        参数:
            data (pd.DataFrame): 富集分析结果数据框，需包含以下列：
                - Term_Name: 富集条目名称
                - Adjusted_P_Value: 校正后的P值
                - Gene_Count: 基因数量
                - Rich_Factor: 富集因子
            title (str): 图表标题
            output_file (str): 输出文件路径（通常为 .pdf 格式）
            top_n (int): 显示前N个最显著的富集条目，默认为20
            color (str): 柱状图填充颜色（十六进制颜色码），默认为 "#3498db"（蓝色）

        返回:
            str: 生成的图表文件路径

        R脚本说明:
            - 使用 ggplot2 的 geom_bar 绘制水平柱状图
            - coord_flip() 将坐标轴翻转，使条目名称显示在纵轴
            - reorder() 按显著性对条目进行排序
            - ggsave() 以 300dpi 保存为 PDF 格式
        """
        # 取前 top_n 条数据并复制，避免修改原始数据
        plot_data = data.head(top_n).copy()
        # 计算 -log10(校正后P值)，值越大表示越显著
        plot_data['neg_log_pval'] = -np.log10(plot_data['Adjusted_P_Value'])
        # 按显著性升序排列（绘图时翻转后最显著的在最上方）
        plot_data = plot_data.sort_values('neg_log_pval', ascending=True)
        
        # 动态生成 R 脚本：将 Python 数据嵌入 R 代码中
        r_script = f'''
library(ggplot2)

# 将 Python 端的数据嵌入 R 数据框
data <- data.frame(
    Term = c({",".join([f'"{t}"' for t in plot_data["Term_Name"]])}),
    NegLogPval = c({",".join([str(v) for v in plot_data["neg_log_pval"]])}),
    GeneCount = c({",".join([str(v) for v in plot_data["Gene_Count"]])}),
    RichFactor = c({",".join([str(v) for v in plot_data["Rich_Factor"]])})
)

# 使用 ggplot2 绘制水平柱状图
p <- ggplot(data, aes(x = reorder(Term, NegLogPval), y = NegLogPval)) +
    geom_bar(stat = "identity", fill = "{color}") +    # 静态柱状图，使用指定颜色填充
    coord_flip() +                                       # 翻转坐标轴，使条目名称竖排显示
    labs(
        title = "{title}",
        x = "",
        y = "-log10(Adjusted P-value)"
    ) +
    theme_bw() +                                         # 使用黑白主题
    theme(
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),  # 标题居中加粗
        axis.text.y = element_text(size = 10),            # 纵轴文字大小
        axis.text.x = element_text(size = 10),            # 横轴文字大小
        panel.grid.minor = element_blank()                # 隐藏次要网格线
    )

# 保存图表，宽度10英寸，高度8英寸，分辨率300dpi
ggsave("{output_file}", p, width = 10, height = 8, dpi = 300)
'''
        
        # 将 R 脚本写入临时文件
        script_file = self.output_dir / "temp_barplot.R"
        with open(script_file, 'w') as f:
            f.write(r_script)
        
        # 调用系统 Rscript 命令执行脚本
        try:
            result = subprocess.run(
                ["Rscript", script_file],
                capture_output=True, text=True, timeout=300
            )
            if result.returncode != 0:
                logger.warning(f"R 脚本执行失败 (返回码: {result.returncode})，请检查是否已安装所需的 R 包")
                if result.stderr:
                    logger.debug(f"R stderr: {result.stderr[:500]}")
        except subprocess.TimeoutExpired:
            logger.warning(f"R 脚本执行超时（300秒），请检查数据量是否过大")
        except Exception as e:
            logger.warning(f"R 脚本执行异常: {e}")
        finally:
            if os.path.exists(script_file):
                os.remove(script_file)
        
        return output_file
    
    def plot_bubble(
        self,
        data: pd.DataFrame,
        title: str,
        output_file: str,
        top_n: int = 20
    ) -> str:
        """
        生成富集分析气泡图

        以 Rich Factor（富集因子）为横轴，富集条目名称为纵轴，生成气泡图。
        气泡大小表示基因数量，颜色表示 -log10(P-value) 的显著性程度。

        参数:
            data (pd.DataFrame): 富集分析结果数据框，需包含以下列：
                - Term_Name: 富集条目名称
                - Adjusted_P_Value: 校正后的P值
                - Gene_Count: 基因数量
                - Rich_Factor: 富集因子（富集基因数/背景基因数）
            title (str): 图表标题
            output_file (str): 输出文件路径（通常为 .pdf 格式）
            top_n (int): 显示前N个最显著的富集条目，默认为20

        返回:
            str: 生成的图表文件路径

        R脚本说明:
            - 使用 geom_point 绘制气泡，aes(size=GeneCount, color=NegLogPval) 映射大小和颜色
            - scale_color_gradient 设置颜色渐变：蓝色(#56B4E9)到橙红色(#D55E00)
            - scale_size_continuous 限制气泡大小范围为 3~10
        """
        # 取前 top_n 条数据并复制
        plot_data = data.head(top_n).copy()
        # 计算 -log10(校正后P值)，用于颜色映射
        plot_data['neg_log_pval'] = -np.log10(plot_data['Adjusted_P_Value'])
        
        # 动态生成 R 脚本：使用 ggplot2 绘制气泡图
        r_script = f'''
library(ggplot2)

# 将 Python 端的数据嵌入 R 数据框
data <- data.frame(
    Term = c({",".join([f'"{t}"' for t in plot_data["Term_Name"]])}),
    RichFactor = c({",".join([str(v) for v in plot_data["Rich_Factor"]])}),
    GeneCount = c({",".join([str(v) for v in plot_data["Gene_Count"]])}),
    NegLogPval = c({",".join([str(v) for v in plot_data["neg_log_pval"]])})
)

# 使用 ggplot2 绘制气泡图
p <- ggplot(data, aes(x = RichFactor, y = reorder(Term, RichFactor))) +
    geom_point(aes(size = GeneCount, color = NegLogPval)) +   # 气泡图：大小=基因数，颜色=显著性
    scale_color_gradient(low = "#56B4E9", high = "#D55E00") + # 颜色渐变：蓝色（低显著性）到橙红色（高显著性）
    scale_size_continuous(range = c(3, 10)) +                  # 气泡大小范围限制
    labs(
        title = "{title}",
        x = "Rich Factor",
        y = "",
        size = "Gene Count",
        color = "-log10(P-value)"
    ) +
    theme_bw() +
    theme(
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        axis.text.y = element_text(size = 10),
        legend.position = "right"                               # 图例放在右侧
    )

# 保存图表
ggsave("{output_file}", p, width = 10, height = 8, dpi = 300)
'''
        
        # 将 R 脚本写入临时文件
        script_file = self.output_dir / "temp_bubble.R"
        with open(script_file, 'w') as f:
            f.write(r_script)
        
        # 调用系统 Rscript 命令执行脚本，然后删除临时文件
        try:
            result = subprocess.run(
                ["Rscript", script_file],
                capture_output=True, text=True, timeout=300
            )
            if result.returncode != 0:
                logger.warning(f"R 脚本执行失败 (返回码: {result.returncode})，请检查是否已安装所需的 R 包")
                if result.stderr:
                    logger.debug(f"R stderr: {result.stderr[:500]}")
        except subprocess.TimeoutExpired:
            logger.warning(f"R 脚本执行超时（300秒），请检查数据量是否过大")
        except Exception as e:
            logger.warning(f"R 脚本执行异常: {e}")
        finally:
            if os.path.exists(script_file):
                os.remove(script_file)
        
        return output_file
    
    def plot_dotplot(
        self,
        data: pd.DataFrame,
        title: str,
        output_file: str,
        top_n: int = 20
    ) -> str:
        """
        生成富集分析点图（Dot Plot）

        以 Gene Ratio（基因比例）为横轴，富集条目名称为纵轴，生成点图。
        点的大小表示基因数量，颜色表示 -log10(P-value) 的显著性程度。
        与气泡图的区别在于横轴使用 Gene Ratio 而非 Rich Factor。

        参数:
            data (pd.DataFrame): 富集分析结果数据框，需包含以下列：
                - Term_Name: 富集条目名称
                - Adjusted_P_Value: 校正后的P值
                - Gene_Count: 基因数量
                - Gene_Ratio: 基因比例（格式为 "分子/分母"，如 "15/200"）
            title (str): 图表标题
            output_file (str): 输出文件路径（通常为 .pdf 格式）
            top_n (int): 显示前N个最显著的富集条目，默认为20

        返回:
            str: 生成的图表文件路径

        R脚本说明:
            - 使用 geom_point 绘制点图，映射大小和颜色
            - 颜色渐变：蓝色(#3498DB)到红色(#E74C3C)
        """
        # 取前 top_n 条数据并复制
        plot_data = data.head(top_n).copy()
        
        # 解析 Gene Ratio 字符串（格式如 "15/200"），计算为浮点数比例
        def parse_ratio(ratio_str):
            try:
                num, denom = ratio_str.split('/')
                return float(num) / float(denom)
            except:
                return 0
        
        # 对每行数据应用比例解析
        plot_data['gene_ratio'] = plot_data['Gene_Ratio'].apply(parse_ratio)
        # 计算 -log10(校正后P值)，用于颜色映射
        plot_data['neg_log_pval'] = -np.log10(plot_data['Adjusted_P_Value'])
        
        # 动态生成 R 脚本：使用 ggplot2 绘制点图
        r_script = f'''
library(ggplot2)

# 将 Python 端的数据嵌入 R 数据框
data <- data.frame(
    Term = c({",".join([f'"{t}"' for t in plot_data["Term_Name"]])}),
    GeneRatio = c({",".join([str(v) for v in plot_data["gene_ratio"]])}),
    GeneCount = c({",".join([str(v) for v in plot_data["Gene_Count"]])}),
    NegLogPval = c({",".join([str(v) for v in plot_data["neg_log_pval"]])})
)

# 使用 ggplot2 绘制点图
p <- ggplot(data, aes(x = GeneRatio, y = reorder(Term, GeneRatio))) +
    geom_point(aes(size = GeneCount, color = NegLogPval)) +    # 点图：大小=基因数，颜色=显著性
    scale_color_gradient(low = "#3498DB", high = "#E74C3C") +  # 颜色渐变：蓝色（低显著性）到红色（高显著性）
    scale_size_continuous(range = c(3, 10)) +                   # 点大小范围限制
    labs(
        title = "{title}",
        x = "Gene Ratio",
        y = "",
        size = "Gene Count",
        color = "-log10(P-value)"
    ) +
    theme_bw() +
    theme(
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        axis.text.y = element_text(size = 10)
    )

# 保存图表
ggsave("{output_file}", p, width = 10, height = 8, dpi = 300)
'''
        
        # 将 R 脚本写入临时文件
        script_file = self.output_dir / "temp_dotplot.R"
        with open(script_file, 'w') as f:
            f.write(r_script)
        
        # 调用系统 Rscript 命令执行脚本，然后删除临时文件
        try:
            result = subprocess.run(
                ["Rscript", script_file],
                capture_output=True, text=True, timeout=300
            )
            if result.returncode != 0:
                logger.warning(f"R 脚本执行失败 (返回码: {result.returncode})，请检查是否已安装所需的 R 包")
                if result.stderr:
                    logger.debug(f"R stderr: {result.stderr[:500]}")
        except subprocess.TimeoutExpired:
            logger.warning(f"R 脚本执行超时（300秒），请检查数据量是否过大")
        except Exception as e:
            logger.warning(f"R 脚本执行异常: {e}")
        finally:
            if os.path.exists(script_file):
                os.remove(script_file)
        
        return output_file
    
    def plot_enrichment_map(
        self,
        data: pd.DataFrame,
        output_file: str,
        similarity_threshold: float = 0.3
    ) -> str:
        """
        生成富集图谱（Enrichment Map）

        基于富集条目之间的基因集合 Jaccard 相似性构建网络图，
        展示不同富集条目之间的功能关联关系。相似性超过阈值的条目对
        之间会绘制连线，边的粗细表示相似性大小。

        参数:
            data (pd.DataFrame): 富集分析结果数据框，需包含以下列：
                - Term_ID: 富集条目ID
                - Genes: 关联基因列表（分号分隔的字符串，如 "GENE1;GENE2;GENE3"）
            output_file (str): 输出文件路径（通常为 .pdf 格式）
            similarity_threshold (float): Jaccard 相似性阈值，仅当两个条目的
                相似性 >= 此阈值时才绘制连线，默认为 0.3

        返回:
            str: 生成的图表文件路径

        算法说明:
            1. 解析每个富集条目的基因集合
            2. 计算所有条目对之间的 Jaccard 相似性
               Jaccard(A, B) = |A ∩ B| / |A ∪ B|
            3. 筛选相似性 >= threshold 的边
            4. 使用 R 的 igraph 包构建无向图并可视化
            5. 使用 Fruchterman-Reingold 力导向布局算法

        R脚本说明:
            - 使用 igraph 的 graph_from_data_frame 从边列表构建图
            - layout_with_fr 使用 Fruchterman-Reingold 力导向布局
            - edge.width 按相似性权重设置边的粗细（权重 * 5）
        """
        # 定义 Jaccard 相似性计算函数
        # Jaccard 相似性 = 两个集合的交集大小 / 并集大小
        # 取值范围为 [0, 1]，值越大表示两个集合越相似
        def jaccard_similarity(set1: set, set2: set) -> float:
            intersection = len(set1 & set2)  # 交集：两个条目共有的基因数
            union = len(set1 | set2)         # 并集：两个条目涉及的所有基因数
            return intersection / union if union > 0 else 0  # 避免除以零
        
        # 解析每个富集条目关联的基因集合
        # Genes 列为分号分隔的基因名称字符串
        gene_sets = {}
        for _, row in data.iterrows():
            genes = set(row['Genes'].split(';')) if 'Genes' in row else set()
            gene_sets[row['Term_ID']] = genes
        
        # 计算所有富集条目对之间的 Jaccard 相似性，构建边列表
        edges = []
        terms = list(gene_sets.keys())
        for i, term1 in enumerate(terms):
            for term2 in terms[i+1:]:  # 仅计算上三角，避免重复
                sim = jaccard_similarity(gene_sets[term1], gene_sets[term2])
                if sim >= similarity_threshold:  # 仅保留相似性达到阈值的边
                    edges.append((term1, term2, sim))
        
        # 动态生成 R 脚本：使用 igraph 构建并可视化富集图谱
        r_script = f'''
library(igraph)
library(ggplot2)

# 从边列表数据框创建无向图
edges <- data.frame(
    from = c({",".join([f'"{e[0]}"' for e in edges]) if edges else ""}),
    to = c({",".join([f'"{e[1]}"' for e in edges]) if edges else ""}),
    weight = c({",".join([str(e[2]) for e in edges]) if edges else ""})
)

if(nrow(edges) > 0) {{
    # graph_from_data_frame: 从数据框创建图对象，directed=FALSE 表示无向图
    g <- graph_from_data_frame(edges, directed = FALSE)
    
    # 使用 PDF 设备输出图表
    pdf("{output_file}", width = 12, height = 10)
    plot(g,
         vertex.size = 15,                                 # 节点大小
         vertex.label.cex = 0.8,                           # 节点标签文字大小
         vertex.label.color = "black",                     # 节点标签颜色
         edge.width = E(g)$weight * 5,                     # 边的粗细 = 相似性权重 × 5
         layout = layout_with_fr(g),                       # Fruchterman-Reingold 力导向布局
         main = "Enrichment Map"
    )
    dev.off()
}} else {{
    # 如果没有满足阈值的边，则输出提示信息
    pdf("{output_file}", width = 10, height = 8)
    plot.new()
    text(0.5, 0.5, "No significant term overlaps", cex = 1.5)
    dev.off()
}}
'''
        
        # 将 R 脚本写入临时文件
        script_file = self.output_dir / "temp_network.R"
        with open(script_file, 'w') as f:
            f.write(r_script)
        
        # 调用系统 Rscript 命令执行脚本，然后删除临时文件
        try:
            result = subprocess.run(
                ["Rscript", script_file],
                capture_output=True, text=True, timeout=300
            )
            if result.returncode != 0:
                logger.warning(f"R 脚本执行失败 (返回码: {result.returncode})，请检查是否已安装所需的 R 包")
                if result.stderr:
                    logger.debug(f"R stderr: {result.stderr[:500]}")
        except subprocess.TimeoutExpired:
            logger.warning(f"R 脚本执行超时（300秒），请检查数据量是否过大")
        except Exception as e:
            logger.warning(f"R 脚本执行异常: {e}")
        finally:
            if os.path.exists(script_file):
                os.remove(script_file)
        
        return output_file
    
    def plot_heatmap(
        self,
        results: List,
        output_file: str,
        database: str = None,
        top_n: int = 20,
        title: str = "Enrichment Heatmap"
    ) -> str:
        """
        生成富集分析热图（Heatmap）

        展示 top N 个显著富集条目与输入基因的关联关系。
        行 = 基因，列 = 条目，颜色 = 是否关联（1=关联，0=不关联）。
        使用 R 的 pheatmap 包绘制，包含行/列聚类树状图。

        参数:
            results (List): 富集分析结果列表，每个元素需包含以下属性：
                - term_id: 条目ID
                - term_name: 条目名称
                - gene_list: 命中的基因名称列表
            output_file (str): 输出文件路径（通常为 .pdf 格式）
            database (str): 数据库名称（用于日志记录），默认为 None
            top_n (int): 显示前N个最显著的富集条目，默认为20
            title (str): 图表标题，默认为 "Enrichment Heatmap"

        返回:
            str: 生成的图表文件路径

        R脚本说明:
            - 使用 pheatmap 包绘制热图
            - 构建基因 x 条目的二值矩阵（1=关联，0=不关联）
            - 行名 = 基因，列名 = 条目ID
            - 颜色方案：白色（0）到红色（1）
            - 启用行和列的层次聚类树状图
            - 如果条目太多，只取 top_n 个最显著的
        """
        # 如果结果列表为空，输出占位图表并返回
        if not results:
            logger.warning("plot_heatmap: 结果列表为空，输出占位图表")
            r_script = f'''
pdf("{output_file}", width = 10, height = 8)
plot.new()
text(0.5, 0.5, "No data for heatmap", cex = 1.5)
dev.off()
'''
            script_file = self.output_dir / "temp_heatmap.R"
            with open(script_file, 'w') as f:
                f.write(r_script)
            try:
                subprocess.run(
                    ["Rscript", script_file],
                    capture_output=True, text=True, timeout=300
                )
            except subprocess.TimeoutExpired:
                logger.warning(f"R 脚本执行超时（300秒），请检查数据量是否过大")
            except Exception as e:
                logger.warning(f"R 脚本执行异常: {e}")
            finally:
                if os.path.exists(script_file):
                    os.remove(script_file)
            return output_file

        # 仅取前 top_n 个富集条目
        top_results = results[:top_n]

        # 收集所有涉及的唯一基因名称（用于构建矩阵行）
        all_genes = set()
        for r in top_results:
            all_genes.update(r.gene_list)
        all_genes = sorted(all_genes)  # 排序以确保可重复性

        # 如果没有基因数据，输出占位图表
        if not all_genes:
            logger.warning("plot_heatmap: 没有基因数据，输出占位图表")
            r_script = f'''
pdf("{output_file}", width = 10, height = 8)
plot.new()
text(0.5, 0.5, "No gene data for heatmap", cex = 1.5)
dev.off()
'''
            script_file = self.output_dir / "temp_heatmap.R"
            with open(script_file, 'w') as f:
                f.write(r_script)
            try:
                subprocess.run(
                    ["Rscript", script_file],
                    capture_output=True, text=True, timeout=300
                )
            except subprocess.TimeoutExpired:
                logger.warning(f"R 脚本执行超时（300秒），请检查数据量是否过大")
            except Exception as e:
                logger.warning(f"R 脚本执行异常: {e}")
            finally:
                if os.path.exists(script_file):
                    os.remove(script_file)
            return output_file

        # 构建基因 x 条目的二值关联矩阵
        # mat_dict: {基因名: {条目索引: 0或1}}
        # 行 = 基因，列 = 条目
        mat_rows = []
        for gene in all_genes:
            row = []
            for r in top_results:
                # 1 表示该基因与该条目关联，0 表示不关联
                row.append(1 if gene in r.gene_list else 0)
            mat_rows.append(row)

        # 为 R 脚本准备数据：行名（基因）和列名（条目ID）
        gene_names = [g.replace('"', '\\"').replace("'", "\\'") for g in all_genes]
        term_ids = []
        for r in top_results:
            tid = getattr(r, 'term_id', getattr(r, 'Term_ID', str(r)))
            term_ids.append(str(tid).replace('"', '\\"').replace("'", "\\'"))

        # 将矩阵数据序列化为 R 向量格式（按列优先顺序，即 R 的 matrix 默认行为）
        # R 的 matrix() 默认按列填充，所以需要转置：将 Python 的行优先转为列优先
        n_genes = len(all_genes)
        n_terms = len(top_results)
        # 按列展开矩阵数据（R matrix bycol=FALSE 默认按列填充）
        mat_data_by_col = []
        for col_idx in range(n_terms):
            for row_idx in range(n_genes):
                mat_data_by_col.append(mat_rows[row_idx][col_idx])

        mat_data_str = ",".join(str(v) for v in mat_data_by_col)

        # 动态生成 R 脚本：使用 pheatmap 绘制热图
        r_script = f'''
# 尝试加载 pheatmap 包，如果未安装则给出提示并退出
if (!requireNamespace("pheatmap", quietly = TRUE)) {{
    message("错误: pheatmap 包未安装。请运行: install.packages('pheatmap')")
    # 输出占位图表
    pdf("{output_file}", width = 10, height = 8)
    plot.new()
    text(0.5, 0.5, "pheatmap package not available\\nPlease install: install.packages('pheatmap')", cex = 1.2)
    dev.off()
    quit(status = 1)
}}

library(pheatmap)

# 构建基因 x 条目的二值关联矩阵
# 行名 = 基因名称，列名 = 条目ID
# 值 = 1（基因与条目关联）或 0（不关联）
mat <- matrix(
    c({mat_data_str}),
    nrow = {n_genes},
    ncol = {n_terms},
    byrow = FALSE
)
rownames(mat) <- c({",".join([f'"{g}"' for g in gene_names])})
colnames(mat) <- c({",".join([f'"{t}"' for t in term_ids])})

# 过滤掉全为零的行（即没有任何关联的基因），避免热图过于稀疏
mat <- mat[rowSums(mat) > 0, , drop = FALSE]

# 如果过滤后矩阵为空，输出提示信息
if (nrow(mat) == 0 || ncol(mat) == 0) {{
    pdf("{output_file}", width = 10, height = 8)
    plot.new()
    text(0.5, 0.5, "No significant associations for heatmap", cex = 1.2)
    dev.off()
    quit(status = 0)
}}

# 根据基因数量动态调整图表高度，避免基因过多时图表过于拥挤
# 每个基因约分配 0.15 英寸的高度，最小 8 英寸，最大 30 英寸
plot_height <- max(8, min(30, nrow(mat) * 0.15 + 2))

# 使用 pheatmap 绘制热图
# - cluster_rows=TRUE: 对基因（行）进行层次聚类
# - cluster_cols=TRUE: 对条目（列）进行层次聚类
# - color: 白色(0)到红色(1)的渐变色
# - fontsize: 字体大小设为8，适应较多基因/条目的场景
# - main: 图表标题
p <- pheatmap(
    mat,
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    color = colorRampPalette(c("white", "red"))(100),
    main = "{title}",
    fontsize = 8,
    fontsize_row = 6,
    fontsize_col = 8,
    border_color = NA,
    silent = TRUE
)

# 保存图表为 PDF 文件
pdf("{output_file}", width = 12, height = plot_height)
grid::grid.draw(p$gtable)
dev.off()
'''

        # 将 R 脚本写入临时文件
        script_file = self.output_dir / "temp_heatmap.R"
        with open(script_file, 'w') as f:
            f.write(r_script)

        # 调用系统 Rscript 命令执行脚本
        try:
            result = subprocess.run(
                ["Rscript", script_file],
                capture_output=True, text=True, timeout=300
            )
            if result.returncode != 0:
                logger.error("plot_heatmap: R 脚本执行失败，请检查是否已安装 pheatmap 包")
                if result.stderr:
                    logger.debug(f"R stderr: {result.stderr[:500]}")
        except subprocess.TimeoutExpired:
            logger.warning(f"R 脚本执行超时（300秒），请检查数据量是否过大")
        except Exception as e:
            logger.warning(f"R 脚本执行异常: {e}")
        finally:
            if os.path.exists(script_file):
                os.remove(script_file)

        return output_file
    
    def plot_upset(
        self,
        results_dict: Dict[str, List],
        output_file: str,
        top_n: int = 20,
        title: str = "Gene Set Overlap"
    ) -> str:
        """
        生成 UpSet 图（UpSet Plot）

        展示不同数据库/条目之间的基因重叠关系。UpSet 图是韦恩图的高级替代方案，
        当集合数量较多时比韦恩图更清晰直观。使用 R 的 UpSetR 包绘制。

        参数:
            results_dict (Dict[str, List]): 按数据库分组的富集结果字典，
                键为数据库名称（如 "GO_BP", "KEGG"），值为该数据库的富集分析结果列表。
                每个结果元素需包含 gene_list 属性（命中的基因名称列表）。
            output_file (str): 输出文件路径（通常为 .pdf 格式）
            top_n (int): 每个数据库取前N个条目参与 UpSet 图计算，默认为20
            title (str): 图表标题，默认为 "Gene Set Overlap"

        返回:
            str: 生成的图表文件路径

        R脚本说明:
            - 使用 UpSetR 包绘制 UpSet 图
            - 构建基因集列表（每个条目对应的基因集合）
            - fromList() 将命名列表转换为 UpSetR 所需的输入格式
            - order.by="freq" 按交集频率排序
            - main.bar.color 设置主条形图颜色
            - sets.bar.color 设置集合大小条形图颜色
        """
        # 如果输入字典为空，输出占位图表并返回
        if not results_dict:
            logger.warning("plot_upset: 结果字典为空，输出占位图表")
            r_script = f'''
pdf("{output_file}", width = 12, height = 8)
plot.new()
text(0.5, 0.5, "No data for UpSet plot", cex = 1.5)
dev.off()
'''
            script_file = self.output_dir / "temp_upset.R"
            with open(script_file, 'w') as f:
                f.write(r_script)
            try:
                subprocess.run(
                    ["Rscript", script_file],
                    capture_output=True, text=True, timeout=300
                )
            except subprocess.TimeoutExpired:
                logger.warning(f"R 脚本执行超时（300秒），请检查数据量是否过大")
            except Exception as e:
                logger.warning(f"R 脚本执行异常: {e}")
            finally:
                if os.path.exists(script_file):
                    os.remove(script_file)
            return output_file

        # 从每个数据库中收集前 top_n 个条目的基因集合
        # 构建命名列表：键 = "数据库名_条目ID"，值 = 基因列表
        gene_sets = {}
        for db_name, results in results_dict.items():
            if not results:
                continue
            # 仅取前 top_n 个条目
            for r in results[:top_n]:
                # 获取条目ID，兼容 EnrichmentResult 对象和字典两种格式
                if hasattr(r, 'term_id'):
                    tid = r.term_id
                elif isinstance(r, dict):
                    tid = r.get('Term_ID', r.get('term_id', 'unknown'))
                else:
                    tid = str(r)

                # 获取基因列表，兼容 EnrichmentResult 对象和字典两种格式
                if hasattr(r, 'gene_list'):
                    genes = r.gene_list
                elif isinstance(r, dict):
                    genes_str = r.get('Genes', '')
                    genes = genes_str.split(';') if genes_str else []
                else:
                    genes = []

                # 使用 "数据库名_条目ID" 作为集合名称，避免不同数据库间 ID 冲突
                set_name = f"{db_name}_{tid}"
                # 去重并过滤空基因名
                gene_sets[set_name] = list(set(g for g in genes if g))

        # 如果没有有效的基因集合数据，输出占位图表
        if not gene_sets:
            logger.warning("plot_upset: 没有有效的基因集合数据，输出占位图表")
            r_script = f'''
pdf("{output_file}", width = 12, height = 8)
plot.new()
text(0.5, 0.5, "No valid gene sets for UpSet plot", cex = 1.2)
dev.off()
'''
            script_file = self.output_dir / "temp_upset.R"
            with open(script_file, 'w') as f:
                f.write(r_script)
            try:
                subprocess.run(
                    ["Rscript", script_file],
                    capture_output=True, text=True, timeout=300
                )
            except subprocess.TimeoutExpired:
                logger.warning(f"R 脚本执行超时（300秒），请检查数据量是否过大")
            except Exception as e:
                logger.warning(f"R 脚本执行异常: {e}")
            finally:
                if os.path.exists(script_file):
                    os.remove(script_file)
            return output_file

        # 为 R 脚本准备数据：将基因集列表序列化为 R 的 list() 格式
        # R 的 list(name1 = c("gene1", "gene2"), name2 = c("gene3", "gene4"), ...)
        list_elements = []
        for set_name, genes in gene_sets.items():
            # 转义集合名称中的特殊字符
            safe_name = set_name.replace('"', '\\"').replace("'", "\\'")
            # 转义基因名称中的特殊字符
            safe_genes = [g.replace('"', '\\"').replace("'", "\\'") for g in genes]
            genes_str = ",".join([f'"{g}"' for g in safe_genes])
            list_elements.append(f'`{safe_name}` = c({genes_str})')

        list_str = ",\n    ".join(list_elements)

        # 计算实际参与 UpSet 图的集合数量（UpSetR 的 nsets 参数）
        n_sets = len(gene_sets)

        # 动态生成 R 脚本：使用 UpSetR 绘制 UpSet 图
        r_script = f'''
# 尝试加载 UpSetR 包，如果未安装则给出提示并退出
if (!requireNamespace("UpSetR", quietly = TRUE)) {{
    message("错误: UpSetR 包未安装。请运行: install.packages('UpSetR')")
    # 输出占位图表
    pdf("{output_file}", width = 12, height = 8)
    plot.new()
    text(0.5, 0.5, "UpSetR package not available\\nPlease install: install.packages('UpSetR')", cex = 1.2)
    dev.off()
    quit(status = 1)
}}

library(UpSetR)

# 构建基因集列表（命名列表格式）
# 每个元素代表一个富集条目，值为该条目关联的基因集合
listInput <- list(
    {list_str}
)

# 检查列表是否有效
if (length(listInput) == 0 || all(sapply(listInput, length) == 0)) {{
    pdf("{output_file}", width = 12, height = 8)
    plot.new()
    text(0.5, 0.5, "No valid gene sets for UpSet plot", cex = 1.2)
    dev.off()
    quit(status = 0)
}}

# 使用 UpSetR 绘制 UpSet 图
# - fromList(): 将命名列表转换为 UpSetR 输入格式
# - order.by="freq": 按交集大小（频率）降序排列
# - nsets: 显示的集合数量（即条目数量）
# - main.bar.color: 主条形图（交集大小）的颜色
# - sets.bar.color: 集合大小条形图的颜色
# - number.angles: 交集大小数字标签的角度（0=水平）
# - point.size: 交集矩阵中点的大小
# - line.size: 交集连线粗细
pdf("{output_file}", width = 14, height = 10)
upset(
    fromList(listInput),
    order.by = "freq",
    nsets = {n_sets},
    main.bar.color = "steelblue",
    sets.bar.color = "gray",
    number.angles = 0,
    point.size = 3.5,
    line.size = 1,
    mainbar.y.label = "Gene Intersection Size",
    sets.x.label = "Gene Set Size",
    text.scale = c(1.2, 1.2, 1.2, 1.2, 1.5, 0.8)
)
title(main = "{title}", cex.main = 1.5)
dev.off()
'''

        # 将 R 脚本写入临时文件
        script_file = self.output_dir / "temp_upset.R"
        with open(script_file, 'w') as f:
            f.write(r_script)

        # 调用系统 Rscript 命令执行脚本
        try:
            result = subprocess.run(
                ["Rscript", script_file],
                capture_output=True, text=True, timeout=300
            )
            if result.returncode != 0:
                logger.error("plot_upset: R 脚本执行失败，请检查是否已安装 UpSetR 包")
                if result.stderr:
                    logger.debug(f"R stderr: {result.stderr[:500]}")
        except subprocess.TimeoutExpired:
            logger.warning(f"R 脚本执行超时（300秒），请检查数据量是否过大")
        except Exception as e:
            logger.warning(f"R 脚本执行异常: {e}")
        finally:
            if os.path.exists(script_file):
                os.remove(script_file)

        return output_file
    
    def plot_cnet(
        self,
        data: pd.DataFrame,
        output_file: str,
        top_n: int = 10
    ) -> str:
        """
        生成基因-条目网络图（Concept Network Plot）

        将富集条目与其关联的基因构建为二分网络图（bipartite graph），
        直观展示基因与富集条目之间的归属关系。条目节点和基因节点
        使用不同颜色和大小进行区分。

        参数:
            data (pd.DataFrame): 富集分析结果数据框，需包含以下列：
                - Term_Name: 富集条目名称
                - Genes: 关联基因列表（分号分隔的字符串，如 "GENE1;GENE2;GENE3"）
            output_file (str): 输出文件路径（通常为 .pdf 格式）
            top_n (int): 取前N个最显著的富集条目构建网络，默认为10

        返回:
            str: 生成的图表文件路径

        R脚本说明:
            - 使用 igraph 构建二分网络图
            - 条目节点（Term）：红色(#E74C3C)，大小为20
            - 基因节点（Gene）：蓝色(#3498DB)，大小为8
            - 使用 Fruchterman-Reingold 力导向布局
            - 边列表最多保留500条（防止网络过于密集）
        """
        # 取前 top_n 条富集结果并复制
        plot_data = data.head(top_n).copy()
        
        # 构建边列表：每个富集条目与其关联的每个基因之间建立一条边
        edges = []
        for _, row in plot_data.iterrows():
            term = row['Term_Name']
            genes = row['Genes'].split(';') if 'Genes' in row else []
            for gene in genes:
                edges.append((term, gene))  # (条目名称, 基因名称) 构成一条边
        
        # 动态生成 R 脚本：使用 igraph 构建并可视化基因-条目网络
        r_script = f'''
library(igraph)
library(ggplot2)

# 从边列表数据框创建网络图
# 限制最多500条边，避免网络过于密集导致图面混乱
edges <- data.frame(
    from = c({",".join([f'"{e[0]}"' for e in edges[:500]]) if edges else ""}),
    to = c({",".join([f'"{e[1]}"' for e in edges[:500]]) if edges else ""})
)

if(nrow(edges) > 0) {{
    # 从边列表创建无向图
    g <- graph_from_data_frame(edges, directed = FALSE)
    
    # 根据节点名称区分节点类型（条目 vs 基因）
    # 出现在 edges$from 中的节点为条目节点，其余为基因节点
    node_types <- ifelse(V(g)$name %in% unique(edges$from), "Term", "Gene")
    V(g)$color <- ifelse(node_types == "Term", "#E74C3C", "#3498DB")  # 条目=红色，基因=蓝色
    V(g)$size <- ifelse(node_types == "Term", 20, 8)                 # 条目节点更大
    
    # 使用 PDF 设备输出图表
    pdf("{output_file}", width = 14, height = 12)
    plot(g,
         vertex.size = V(g)$size,                            # 节点大小
         vertex.color = V(g)$color,                          # 节点颜色
         vertex.label.cex = 0.6,                             # 节点标签文字大小
         vertex.label.color = "black",                       # 节点标签颜色
         layout = layout_with_fr(g),                         # Fruchterman-Reingold 力导向布局
         main = "Gene-Term Network"
    )
    # 添加图例，说明节点颜色含义
    legend("topright", 
           legend = c("Term", "Gene"),
           col = c("#E74C3C", "#3498DB"),
           pch = 19,
           pt.cex = c(2, 1))
    dev.off()
}} else {{
    # 如果没有边数据，则输出提示信息
    pdf("{output_file}", width = 10, height = 8)
    plot.new()
    text(0.5, 0.5, "No network data", cex = 1.5)
    dev.off()
}}
'''
        
        # 将 R 脚本写入临时文件
        script_file = self.output_dir / "temp_cnet.R"
        with open(script_file, 'w') as f:
            f.write(r_script)
        
        # 调用系统 Rscript 命令执行脚本，然后删除临时文件
        try:
            result = subprocess.run(
                ["Rscript", script_file],
                capture_output=True, text=True, timeout=300
            )
            if result.returncode != 0:
                logger.warning(f"R 脚本执行失败 (返回码: {result.returncode})，请检查是否已安装所需的 R 包")
                if result.stderr:
                    logger.debug(f"R stderr: {result.stderr[:500]}")
        except subprocess.TimeoutExpired:
            logger.warning(f"R 脚本执行超时（300秒），请检查数据量是否过大")
        except Exception as e:
            logger.warning(f"R 脚本执行异常: {e}")
        finally:
            if os.path.exists(script_file):
                os.remove(script_file)
        
        return output_file
    
    def plot_all(
        self,
        data: pd.DataFrame,
        database: str,
        top_n: int = 20
    ) -> Dict[str, str]:
        """
        生成所有标准图表

        为指定数据库的富集分析结果批量生成所有标准图表类型，
        包括柱状图、气泡图、点图和基因-条目网络图。

        参数:
            data (pd.DataFrame): 富集分析结果数据框
            database (str): 数据库名称（如 "GO", "KEGG", "Reactome" 等），
                用于确定图表标题和柱状图颜色
            top_n (int): 每种图表显示前N个最显著的富集条目，默认为20

        返回:
            Dict[str, str]: 字典，键为图表类型名称，值为对应的输出文件路径。
                包含以下键：
                - "barplot": 柱状图文件路径
                - "bubble": 气泡图文件路径
                - "dotplot": 点图文件路径
                - "cnet": 基因-条目网络图文件路径
        """
        plots = {}
        
        # 根据数据库名称选择对应的柱状图颜色方案
        # 每个数据库使用不同的主题色，便于区分不同来源的富集结果
        colors = {
            "GO": "#E64B35",           # GO - 红色系（与 GO BP 颜色一致）
            "KEGG": "#4DBBD5",         # KEGG - 蓝色系（与 GO CC 颜色一致）
            "Reactome": "#3C5488",     # Reactome - 深蓝色
            "WikiPathways": "#00A087", # WikiPathways - 绿色系（与 GO MF 颜色一致）
            "MSigDB": "#F39B7F",       # MSigDB - 浅橙色
            "DO": "#8491B4",           # DO（疾病本体）- 灰紫色
            "DisGeNET": "#91D1C2"      # DisGeNET - 浅青色
        }
        color = colors.get(database, "#3498DB")  # 未匹配的数据库使用默认蓝色
        
        # 构建图表标题
        title = f"{database} Enrichment Analysis"
        
        # 生成柱状图（Bar Plot）
        bar_file = str(self.output_dir / f"{database}_barplot.pdf")
        self.plot_barplot(data, title, bar_file, top_n, color)
        plots["barplot"] = bar_file
        
        # 生成气泡图（Bubble Plot）
        bubble_file = str(self.output_dir / f"{database}_bubble.pdf")
        self.plot_bubble(data, title, bubble_file, top_n)
        plots["bubble"] = bubble_file
        
        # 生成点图（Dot Plot）
        dot_file = str(self.output_dir / f"{database}_dotplot.pdf")
        self.plot_dotplot(data, title, dot_file, top_n)
        plots["dotplot"] = dot_file
        
        # 生成基因-条目网络图（Concept Network Plot）
        cnet_file = str(self.output_dir / f"{database}_cnet.pdf")
        self.plot_cnet(data, cnet_file, top_n)
        plots["cnet"] = cnet_file
        
        # 生成热图（Heatmap）— 仅当数据量充足时生成
        # 将 DataFrame 转换为 EnrichmentResult 列表格式，供 plot_heatmap 使用
        try:
            from allenricher.core.enrichment import EnrichmentResult
            results_list = []
            for _, row in data.iterrows():
                # 解析基因列表（EnrichmentResult.to_dict() 中基因以分号分隔）
                gene_str = str(row.get('Genes', ''))
                gene_list = [g.strip() for g in gene_str.split(';') if g.strip()] if gene_str else []
                results_list.append(EnrichmentResult(
                    term_id=str(row.get('Term_ID', '')),
                    term_name=str(row.get('Term_Name', '')),
                    database=database,
                    pvalue=float(row.get('P_Value', 1)),
                    adjusted_pvalue=float(row.get('Adjusted_P_Value', 1)),
                    gene_count=int(row.get('Gene_Count', 0)),
                    background_count=int(row.get('Background_Count', 0)) if 'Background_Count' in row.index else 0,
                    expected_count=float(row.get('Expected_Count', 0)) if 'Expected_Count' in row.index else 0.0,
                    rich_factor=float(row.get('Rich_Factor', 0)) if 'Rich_Factor' in row.index else 0.0,
                    gene_list=gene_list,
                    gene_ratio=str(row.get('Gene_Ratio', '')),
                    background_ratio=str(row.get('Background_Ratio', ''))
                ))
            if results_list and len(results_list) >= 2:
                heatmap_file = str(self.output_dir / f"{database}_heatmap.pdf")
                self.plot_heatmap(results_list, heatmap_file, database=database)
                plots["heatmap"] = heatmap_file
        except Exception as e:
            logger.warning(f"跳过热图生成: {e}")
        
        return plots
