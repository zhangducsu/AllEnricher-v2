"""
Interactive HTML Report Generator for AllEnricher v2.0

生成交互式HTML报告模块

本模块负责将富集分析的结果数据生成为一个完整的、可交互的HTML报告页面。
主要功能包括：
- 统计摘要（Summary Statistics）：展示输入基因数量、分析数据库数量、富集条目总数等概览信息
- 交互式数据表格（Interactive Tables）：支持排序、分页、搜索的数据表格，并提供TSV下载和剪贴板复制功能
- 嵌入式图表展示（Embedded Plots）：展示柱状图和气泡图等可视化结果
- AI智能解读（AI Interpretation）：可选的AI生成分析解读内容

依赖项：
- pandas：用于数据处理
- base64：用于图片编码（如有需要）
- 标准库：os, json, logging, datetime, pathlib
"""

import os
import json
import logging
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional, Any
import pandas as pd
import base64

logger = logging.getLogger(__name__)


class ReportGenerator:
    """
    交互式HTML报告生成器

    本类负责将富集分析的结果数据整合并生成为一个功能完整的交互式HTML报告。
    报告包含以下核心部分：
    1. 统计摘要 - 展示分析的整体概览信息
    2. 数据表格 - 每个数据库对应的富集结果表格，支持排序、下载和复制
    3. 图表区域 - 嵌入或链接到可视化图表（柱状图、气泡图等）
    4. AI解读 - 可选的AI生成分析结果解读

    使用方式：
        generator = ReportGenerator(output_dir="./reports")
        report_path = generator.generate(
            results=results_dict,
            output_file="./reports/report.html",
            gene_list=["GENE1", "GENE2", ...],
            ai_interpretation={"GO": "...", "KEGG": "..."}
        )
    """
    
    def __init__(self, output_dir: str, config=None):
        """
        初始化报告生成器

        Args:
            output_dir: 报告输出目录路径，用于存放生成的HTML文件和图表文件
            config: 可选的配置对象，用于自定义报告生成行为
        """
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)  # 确保输出目录存在
        self.config = config
    
    def generate(
        self,
        results: Dict[str, pd.DataFrame],
        output_file: str,
        gene_list: List[str] = None,
        ai_interpretation: Dict[str, str] = None
    ) -> str:
        """
        生成完整的HTML报告

        这是报告生成的主入口方法，负责协调各子模块生成报告的各个部分，
        并将最终结果组装为完整的HTML文档写入文件。

        Args:
            results: 富集分析结果字典，键为数据库名称（如 "GO", "KEGG"），
                     值为对应的 pandas DataFrame 结果表格
            output_file: 输出HTML文件的完整路径
            gene_list: 可选，用户输入的基因列表，用于在摘要中展示输入基因数量
            ai_interpretation: 可选，AI生成的分析解读字典，键为数据库名称，
                               值为对应的解读文本内容

        Returns:
            str: 生成的报告文件的路径
        """
        # 检查是否有富集结果
        # 如果结果为空或所有 DataFrame 都为空，生成无结果提示页面
        has_results = results and any(len(df) > 0 for df in results.values())
        
        if not has_results:
            # 生成无结果提示页面
            html = self._generate_no_results_page(gene_list)
            with open(output_file, 'w', encoding='utf-8') as f:
                f.write(html)
            logger.info(f"No enrichment results found. Generated no-results page: {output_file}")
            return output_file
        
        # 生成报告的各个部分
        summary = self._generate_summary(results, gene_list)   # 统计摘要部分
        tables = self._generate_tables(results)                 # 交互式数据表格部分
        plots = self._generate_plot_section(results)             # 图表展示部分
        ai_section = self._generate_ai_section(ai_interpretation) if ai_interpretation else ""  # AI解读部分（可选）
        
        # 构建完整的HTML文档
        # 传入实际有结果的数据库名称，用于动态生成导航栏
        active_db_names = [db for db, df in results.items() if len(df) > 0]
        html = self._build_html(
            summary=summary,
            tables=tables,
            plots=plots,
            ai_section=ai_section,
            db_names=active_db_names
        )
        
        # 将HTML内容写入文件
        with open(output_file, 'w', encoding='utf-8') as f:
            f.write(html)
        
        logger.info(f"Report generated: {output_file}")
        return output_file
    
    def _generate_no_results_page(self, gene_list: List[str] = None) -> str:
        """
        生成无富集结果提示页面

        当富集分析没有找到任何显著结果时，生成一个友好的提示页面，
        包含可能的原因分析和改进建议，帮助用户调整分析参数。

        Args:
            gene_list: 可选，用户输入的基因列表，用于显示输入基因数量

        Returns:
            str: 无结果提示页面的完整HTML字符串
        """
        html = f'''<!DOCTYPE html>
<html lang="zh-CN">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>AllEnricher v2.0 - 无富集结果</title>
    
    <!-- Font Awesome 图标库 -->
    <link rel="stylesheet" href="https://cdnjs.cloudflare.com/ajax/libs/font-awesome/6.4.0/css/all.min.css">
    
    <style>
        /* === CSS变量定义 === */
        :root {{
            --primary-color: #3498db;
            --secondary-color: #2ecc71;
            --accent-color: #e74c3c;
            --warning-color: #f39c12;
            --dark-color: #2c3e50;
            --light-color: #ecf0f1;
        }}
        
        /* === 全局重置 === */
        * {{
            margin: 0;
            padding: 0;
            box-sizing: border-box;
        }}
        
        body {{
            font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif;
            line-height: 1.6;
            color: var(--dark-color);
            background-color: #f5f6fa;
            min-height: 100vh;
            display: flex;
            flex-direction: column;
        }}
        
        /* === 页面头部 === */
        .header {{
            background: linear-gradient(135deg, var(--primary-color), var(--secondary-color));
            color: white;
            padding: 2rem;
            text-align: center;
            box-shadow: 0 2px 10px rgba(0,0,0,0.1);
        }}
        
        .header h1 {{
            font-size: 2.5rem;
            margin-bottom: 0.5rem;
        }}
        
        .header .version {{
            font-size: 1rem;
            opacity: 0.9;
        }}
        
        /* === 主内容区域 === */
        .main {{
            flex: 1;
            max-width: 900px;
            margin: 2rem auto;
            padding: 0 2rem;
            width: 100%;
        }}
        
        /* === 无结果提示卡片 === */
        .no-results-card {{
            background: white;
            border-radius: 10px;
            padding: 2rem;
            box-shadow: 0 2px 10px rgba(0,0,0,0.05);
            text-align: center;
        }}
        
        .no-results-icon {{
            font-size: 5rem;
            color: var(--warning-color);
            margin-bottom: 1rem;
        }}
        
        .no-results-title {{
            font-size: 1.8rem;
            color: var(--dark-color);
            margin-bottom: 1rem;
        }}
        
        .no-results-message {{
            font-size: 1.1rem;
            color: #666;
            margin-bottom: 2rem;
        }}
        
        /* === 原因和建议列表 === */
        .info-section {{
            text-align: left;
            margin-top: 2rem;
        }}
        
        .info-section h3 {{
            color: var(--primary-color);
            margin-bottom: 1rem;
            display: flex;
            align-items: center;
            gap: 0.5rem;
        }}
        
        .info-list {{
            list-style: none;
            padding: 0;
        }}
        
        .info-list li {{
            padding: 0.75rem 1rem;
            margin-bottom: 0.5rem;
            background: var(--light-color);
            border-radius: 5px;
            border-left: 4px solid var(--primary-color);
        }}
        
        .info-list.suggestions li {{
            border-left-color: var(--secondary-color);
        }}
        
        .info-list li i {{
            margin-right: 0.5rem;
            color: var(--primary-color);
        }}
        
        .info-list.suggestions li i {{
            color: var(--secondary-color);
        }}
        
        /* === 统计信息卡片 === */
        .stats-grid {{
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(150px, 1fr));
            gap: 1rem;
            margin-bottom: 2rem;
        }}
        
        .stat-card {{
            background: var(--light-color);
            padding: 1rem;
            border-radius: 8px;
            text-align: center;
        }}
        
        .stat-value {{
            font-size: 1.5rem;
            font-weight: bold;
            color: var(--primary-color);
        }}
        
        .stat-label {{
            font-size: 0.9rem;
            color: #666;
        }}
        
        /* === 页脚 === */
        .footer {{
            text-align: center;
            padding: 2rem;
            color: #666;
            font-size: 0.9rem;
        }}
        
        /* === 响应式布局 === */
        @media (max-width: 768px) {{
            .header h1 {{
                font-size: 1.8rem;
            }}
            
            .no-results-icon {{
                font-size: 4rem;
            }}
            
            .no-results-title {{
                font-size: 1.5rem;
            }}
        }}
    </style>
</head>
<body>
    <!-- 页面头部 -->
    <header class="header">
        <h1><i class="fas fa-dna"></i> AllEnricher Report</h1>
        <p class="version">Version 2.0 | Generated on {datetime.now().strftime("%Y-%m-%d %H:%M")}</p>
    </header>
    
    <!-- 主内容区域 -->
    <main class="main">
        <div class="no-results-card">
            <!-- 无结果图标 -->
            <div class="no-results-icon">
                <i class="fas fa-search"></i>
            </div>
            
            <!-- 无结果标题 -->
            <h2 class="no-results-title">未找到显著富集的结果</h2>
            
            <!-- 无结果说明 -->
            <p class="no-results-message">
                根据当前的分析参数，未能找到任何显著富集的功能条目。
                请参考以下可能的原因和改进建议。
            </p>
            
            <!-- 统计信息 -->
            <div class="stats-grid">
                <div class="stat-card">
                    <div class="stat-value">{len(gene_list) if gene_list else "N/A"}</div>
                    <div class="stat-label">输入基因数</div>
                </div>
                <div class="stat-card">
                    <div class="stat-value">{datetime.now().strftime("%Y-%m-%d")}</div>
                    <div class="stat-label">分析日期</div>
                </div>
            </div>
            
            <!-- 可能的原因 -->
            <div class="info-section">
                <h3><i class="fas fa-question-circle"></i> 可能的原因</h3>
                <ul class="info-list">
                    <li><i class="fas fa-circle"></i> 输入基因列表过小或与数据库无交集</li>
                    <li><i class="fas fa-circle"></i> p 值/q 值阈值过于严格</li>
                    <li><i class="fas fa-circle"></i> 背景基因集设置不当</li>
                    <li><i class="fas fa-circle"></i> 基因 ID 格式与数据库不匹配</li>
                </ul>
            </div>
            
            <!-- 改进建议 -->
            <div class="info-section">
                <h3><i class="fas fa-lightbulb"></i> 改进建议</h3>
                <ul class="info-list suggestions">
                    <li><i class="fas fa-check-circle"></i> 增加输入基因数量（建议至少 10 个以上）</li>
                    <li><i class="fas fa-check-circle"></i> 放宽 p 值/q 值阈值（如 <code>-p 0.1 -q 0.1</code>）</li>
                    <li><i class="fas fa-check-circle"></i> 检查基因 ID 格式是否正确（如 Ensembl ID、Gene Symbol 等）</li>
                    <li><i class="fas fa-check-circle"></i> 尝试使用其他数据库进行分析</li>
                    <li><i class="fas fa-check-circle"></i> 检查背景基因集是否包含输入基因</li>
                </ul>
            </div>
        </div>
    </main>
    
    <!-- 页脚 -->
    <footer class="footer">
        <p>Generated by <a href="https://github.com/zd105/AllEnricher" target="_blank">AllEnricher v2.0</a></p>
        <p>&copy; {datetime.now().year} AllEnricher Team. Licensed under MIT.</p>
    </footer>
</body>
</html>
'''
        return html
    
    def _generate_summary(
        self,
        results: Dict[str, pd.DataFrame],
        gene_list: List[str] = None
    ) -> str:
        """
        生成统计摘要部分

        根据富集分析结果计算并生成报告顶部的统计概览区域，包括：
        - 输入基因数量卡片
        - 分析数据库数量卡片
        - 富集条目总数卡片
        - 分析日期卡片
        - 各数据库的详细统计表格（条目数、最小P值、最小校正P值）

        Args:
            results: 富集分析结果字典，键为数据库名称，值为DataFrame
            gene_list: 可选，用户输入的基因列表

        Returns:
            str: 统计摘要部分的HTML字符串
        """
        # 计算总体统计信息
        total_terms = sum(len(df) for df in results.values())  # 所有数据库的富集条目总数
        databases = list(results.keys())  # 获取所有数据库名称列表
        
        # 逐个数据库收集统计信息
        summary_stats = []
        for db_name, df in results.items():
            if len(df) > 0:
                summary_stats.append({
                    "database": db_name,
                    "terms": len(df),
                    "min_pval": df['P_Value'].min() if 'P_Value' in df.columns else 0,
                    "min_adj_pval": df['Adjusted_P_Value'].min() if 'Adjusted_P_Value' in df.columns else 0
                })
        
        html = f'''
        <div class="section" id="summary">
            <h2><i class="fas fa-chart-pie"></i> Analysis Summary</h2>
            <div class="summary-grid">
                <div class="summary-card">
                    <div class="summary-value">{len(gene_list) if gene_list else "N/A"}</div>
                    <div class="summary-label">Input Genes</div>
                </div>
                <div class="summary-card">
                    <div class="summary-value">{len(databases)}</div>
                    <div class="summary-label">Databases Analyzed</div>
                </div>
                <div class="summary-card">
                    <div class="summary-value">{total_terms}</div>
                    <div class="summary-label">Enriched Terms</div>
                </div>
                <div class="summary-card">
                    <div class="summary-value">{datetime.now().strftime("%Y-%m-%d")}</div>
                    <div class="summary-label">Analysis Date</div>
                </div>
            </div>

            <h3>Database Statistics</h3>
            <table class="data-table">
                <thead>
                    <tr>
                        <th>Database</th>
                        <th>Enriched Terms</th>
                        <th>Min P-value</th>
                        <th>Min Adj. P-value</th>
                    </tr>
                </thead>
                <tbody>
                    {self._render_summary_rows(summary_stats)}
                </tbody>
            </table>
        </div>
        '''

        return html

    @staticmethod
    def _render_summary_rows(summary_stats):
        """
        渲染摘要表格行

        将各数据库的统计信息列表转换为HTML表格行字符串。
        每行包含：数据库名称（带锚点链接）、富集条目数、最小P值、最小校正P值。

        Args:
            summary_stats: 统计信息字典列表，每个字典包含 database, terms, min_pval, min_adj_pval

        Returns:
            str: 拼接后的HTML表格行字符串
        """
        rows = []
        for s in summary_stats:
            db = s["database"]
            terms = s["terms"]
            min_pval = f"{s['min_pval']:.2e}"          # 格式化为科学计数法
            min_adj_pval = f"{s['min_adj_pval']:.2e}"  # 格式化为科学计数法
            rows.append(f'''
                    <tr>
                        <td><a href="#{db}-table">{db}</a></td>
                        <td>{terms}</td>
                        <td>{min_pval}</td>
                        <td>{min_adj_pval}</td>
                    </tr>
                    ''')
        return " ".join(rows)
    
    def _generate_tables(self, results: Dict[str, pd.DataFrame]) -> str:
        """
        生成交互式数据表格

        为每个数据库的富集分析结果生成一个独立的交互式HTML表格。
        表格特性：
        - 支持列排序（通过 DataTables 插件）
        - 显示 Term ID（带超链接）、Term Name、Gene Count、Rich Factor、P-value、Adj. P-value、Genes 列
        - Term ID 列包含指向原始数据库的超链接（如 GO AmiGO、KEGG 等）
        - 提供 TSV 下载和复制到剪贴板功能
        - 基因列超长时自动截断，鼠标悬停可查看完整内容

        Args:
            results: 富集分析结果字典，键为数据库名称，值为DataFrame

        Returns:
            str: 所有数据表格的HTML字符串拼接结果
        """
        tables_html = []
        
        for db_name, df in results.items():
            if len(df) == 0:
                continue  # 跳过空结果
            
            # 逐行生成表格内容
            rows = []
            for idx, row in df.iterrows():
                tid = row.get('Term_ID', 'N/A')           # 条目ID
                tname = row.get('Term_Name', 'N/A')        # 条目名称
                gcount = row.get('Gene_Count', 0)          # 基因数量
                rf = f"{row.get('Rich_Factor', 0):.4f}"   # 富集因子，保留4位小数
                pv = f"{row.get('P_Value', 1):.2e}"       # P值，科学计数法
                adjpv = f"{row.get('Adjusted_P_Value', 1):.2e}"  # 校正后P值，科学计数法
                term_url = row.get('Term_URL', '')         # 条目的数据库链接 URL
                genes_str = str(row.get('Genes', ''))      # 关联基因列表字符串
                genes_short = genes_str[:50] + ('...' if len(genes_str) > 50 else '')  # 基因列表截断显示
                
                # 如果有 URL，将 Term ID 转换为超链接
                if term_url:
                    tid_cell = f'<a href="{term_url}" target="_blank" title="点击查看数据库详情">{tid}</a>'
                else:
                    tid_cell = tid
                
                rows.append(f'''
                <tr>
                    <td>{tid_cell}</td>
                    <td class="term-name">{tname}</td>
                    <td>{gcount}</td>
                    <td>{rf}</td>
                    <td>{pv}</td>
                    <td class="pval">{adjpv}</td>
                    <td class="genes" title="{genes_str}">{genes_short}</td>
                </tr>
                ''')
            
            table_html = f'''
            <div class="section" id="{db_name}-table">
                <h2><i class="fas fa-table"></i> {db_name} Enrichment Results</h2>
                <p class="result-count">Found <strong>{len(df)}</strong> enriched terms</p>
                
                <div class="table-container">
                    <table class="data-table sortable" id="table-{db_name}">
                        <thead>
                            <tr>
                                <th data-sort="string">Term ID</th>
                                <th data-sort="string">Term Name</th>
                                <th data-sort="int">Gene Count</th>
                                <th data-sort="float">Rich Factor</th>
                                <th data-sort="float">P-value</th>
                                <th data-sort="float">Adj. P-value</th>
                                <th>Genes</th>
                            </tr>
                        </thead>
                        <tbody>
                            {" ".join(rows)}
                        </tbody>
                    </table>
                </div>
                
                <div class="table-actions">
                    <button onclick="downloadTable('{db_name}')" class="btn">
                        <i class="fas fa-download"></i> Download TSV
                    </button>
                    <button onclick="copyTable('{db_name}')" class="btn">
                        <i class="fas fa-copy"></i> Copy to Clipboard
                    </button>
                </div>
            </div>
            '''
            
            tables_html.append(table_html)
        
        return "\n".join(tables_html)
    
    def _generate_plot_section(self, results: Dict[str, pd.DataFrame]) -> str:
        """
        生成图表展示区域

        检查输出目录中是否存在各数据库对应的图表文件（PDF格式的柱状图和气泡图），
        如果存在则生成可点击的链接卡片，用户可以在新标签页中打开查看。

        图表文件命名规则：
        - 柱状图：{数据库名}_barplot.pdf
        - 气泡图：{数据库名}_bubble.pdf

        Args:
            results: 富集分析结果字典，键为数据库名称，值为DataFrame

        Returns:
            str: 图表展示区域的HTML字符串
        """
        plots_html = ['''
        <div class="section" id="plots">
            <h2><i class="fas fa-chart-bar"></i> Visualization</h2>
        ''']
        
        for db_name in results.keys():
            plot_dir = self.output_dir / "plots"  # 图表文件存放子目录
            
            # 检查是否存在对应的图表文件
            barplot_file = plot_dir / f"{db_name}_barplot.pdf"   # 柱状图文件路径
            bubble_file = plot_dir / f"{db_name}_bubble.pdf"     # 气泡图文件路径
            
            if barplot_file.exists() or bubble_file.exists():
                plots_html.append(f'''
                <div class="plot-section">
                    <h3>{db_name}</h3>
                    <div class="plot-grid">
                ''')
                
                # 遍历两种图表类型，如果文件存在则生成链接卡片
                for plot_type, plot_file in [("Bar Plot", barplot_file), ("Bubble Plot", bubble_file)]:
                    if plot_file.exists():
                        # PDF文件以链接形式嵌入，点击后在新标签页打开
                        plots_html.append(f'''
                        <div class="plot-card">
                            <h4>{plot_type}</h4>
                            <div class="plot-placeholder">
                                <a href="plots/{plot_file.name}" target="_blank">
                                    <i class="fas fa-file-pdf"></i>
                                    <span>Open {plot_type}</span>
                                </a>
                            </div>
                        </div>
                        ''')
                
                plots_html.append('</div></div>')
        
        plots_html.append('</div>')
        return "\n".join(plots_html)
    
    def _generate_ai_section(self, ai_interpretation: Dict[str, str]) -> str:
        """
        生成AI解读部分

        将AI生成的富集分析解读内容渲染为HTML段落。
        每个数据库的解读内容会以独立卡片的形式展示，
        并附带免责声明，提醒用户AI生成的内容需要领域专家审核。

        Args:
            ai_interpretation: AI解读字典，键为数据库名称，值为解读文本内容

        Returns:
            str: AI解读部分的HTML字符串；如果输入为空则返回空字符串
        """
        if not ai_interpretation:
            return ""  # 无AI解读内容时返回空字符串
        
        # 为每个数据库生成独立的解读卡片
        interpretations = []
        for db_name, interpretation in ai_interpretation.items():
            interpretations.append(f'''
            <div class="ai-interpretation">
                <h3><i class="fas fa-robot"></i> {db_name} Interpretation</h3>
                <div class="interpretation-content">
                    {interpretation}
                </div>
            </div>
            ''')
        
        html = f'''
        <div class="section" id="ai-interpretation">
            <h2><i class="fas fa-brain"></i> AI-Powered Interpretation</h2>
            <p class="ai-disclaimer">
                <i class="fas fa-info-circle"></i>
                The following interpretations are generated by AI and should be reviewed by domain experts.
            </p>
            {" ".join(interpretations)}
        </div>
        '''
        
        return html
    
    def _build_html(self, summary: str, tables: str, plots: str, ai_section: str, db_names: List[str] = None) -> str:
        """
        构建完整HTML文档

        将报告的各个部分（摘要、表格、图表、AI解读）组装到一个完整的HTML页面中。
        页面包含以下组成部分：
        - HTML头部（meta标签、外部CSS/JS引用）
        - 内嵌CSS样式（响应式布局、表格样式、卡片样式等）
        - 页面头部（标题、版本号、生成时间）
        - 导航栏（各部分锚点链接，根据实际分析的数据库动态生成）
        - 主内容区域（摘要 -> 图表 -> 表格 -> AI解读）
        - 页脚（版权信息）
        - 内嵌JavaScript（DataTables初始化、下载和复制功能）

        Args:
            summary: 统计摘要部分的HTML字符串
            tables: 数据表格部分的HTML字符串
            plots: 图表展示部分的HTML字符串
            ai_section: AI解读部分的HTML字符串
            db_names: 实际分析的数据库名称列表，用于动态生成导航栏链接

        Returns:
            str: 完整的HTML文档字符串
        """
        html = f'''<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>AllEnricher v2.0 - Enrichment Analysis Report</title>
    
    <!-- Font Awesome 图标库 -->
    <link rel="stylesheet" href="https://cdnjs.cloudflare.com/ajax/libs/font-awesome/6.4.0/css/all.min.css">
    
    <!-- DataTables 表格插件样式 -->
    <link rel="stylesheet" href="https://cdn.datatables.net/1.13.4/css/jquery.dataTables.min.css">
    
    <style>
        /* === CSS变量定义：统一管理全局颜色主题 === */
        :root {{
            --primary-color: #3498db;     /* 主色调：蓝色 */
            --secondary-color: #2ecc71;   /* 辅助色：绿色 */
            --accent-color: #e74c3c;      /* 强调色：红色 */
            --dark-color: #2c3e50;        /* 深色文字 */
            --light-color: #ecf0f1;       /* 浅色背景 */
            --border-color: #bdc3c7;      /* 边框颜色 */
        }}
        
        /* === 全局重置样式 === */
        * {{
            margin: 0;
            padding: 0;
            box-sizing: border-box;
        }}
        
        /* === 页面主体样式 === */
        body {{
            font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif;
            line-height: 1.6;
            color: var(--dark-color);
            background-color: #f5f6fa;
        }}
        
        /* === 页面头部：渐变背景的标题栏 === */
        .header {{
            background: linear-gradient(135deg, var(--primary-color), var(--secondary-color));
            color: white;
            padding: 2rem;
            text-align: center;
            box-shadow: 0 2px 10px rgba(0,0,0,0.1);
        }}
        
        .header h1 {{
            font-size: 2.5rem;
            margin-bottom: 0.5rem;
        }}
        
        .header .version {{
            font-size: 1rem;
            opacity: 0.9;
        }}
        
        /* === 导航栏：粘性定位，始终显示在页面顶部 === */
        .nav {{
            background: white;
            padding: 1rem;
            position: sticky;
            top: 0;
            z-index: 100;
            box-shadow: 0 2px 5px rgba(0,0,0,0.1);
        }}
        
        .nav ul {{
            list-style: none;
            display: flex;
            justify-content: center;
            gap: 2rem;
            flex-wrap: wrap;
        }}
        
        .nav a {{
            color: var(--dark-color);
            text-decoration: none;
            padding: 0.5rem 1rem;
            border-radius: 5px;
            transition: all 0.3s;
        }}
        
        .nav a:hover {{
            background: var(--primary-color);
            color: white;
        }}
        
        /* === 主内容区域：居中布局，最大宽度1400px === */
        .main {{
            max-width: 1400px;
            margin: 0 auto;
            padding: 2rem;
        }}
        
        /* === 各内容板块：白色圆角卡片，带阴影 === */
        .section {{
            background: white;
            border-radius: 10px;
            padding: 2rem;
            margin-bottom: 2rem;
            box-shadow: 0 2px 10px rgba(0,0,0,0.05);
        }}
        
        .section h2 {{
            color: var(--primary-color);
            margin-bottom: 1.5rem;
            padding-bottom: 0.5rem;
            border-bottom: 2px solid var(--light-color);
        }}
        
        .section h2 i {{
            margin-right: 0.5rem;
        }}
        
        /* === 摘要卡片网格：自适应布局 === */
        .summary-grid {{
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(200px, 1fr));
            gap: 1.5rem;
            margin-bottom: 2rem;
        }}
        
        .summary-card {{
            background: linear-gradient(135deg, var(--primary-color), var(--secondary-color));
            color: white;
            padding: 1.5rem;
            border-radius: 10px;
            text-align: center;
        }}
        
        .summary-value {{
            font-size: 2.5rem;
            font-weight: bold;
        }}
        
        .summary-label {{
            font-size: 0.9rem;
            opacity: 0.9;
            margin-top: 0.5rem;
        }}
        
        /* === 数据表格样式 === */
        .table-container {{
            overflow-x: auto;  /* 水平溢出时显示滚动条 */
        }}
        
        .data-table {{
            width: 100%;
            border-collapse: collapse;
            margin-top: 1rem;
        }}
        
        .data-table th,
        .data-table td {{
            padding: 0.75rem;
            text-align: left;
            border-bottom: 1px solid var(--border-color);
        }}
        
        .data-table th {{
            background: var(--light-color);
            font-weight: 600;
            cursor: pointer;
        }}
        
        .data-table th:hover {{
            background: var(--primary-color);
            color: white;
        }}
        
        .data-table tr:hover {{
            background: #f8f9fa;
        }}
        
        .data-table .pval {{
            font-family: monospace;
        }}
        
        .data-table .genes {{
            font-size: 0.85rem;
            color: #666;
            max-width: 200px;
            overflow: hidden;
            text-overflow: ellipsis;
            white-space: nowrap;
        }}
        
        .data-table .term-name {{
            max-width: 300px;
        }}
        
        /* === 操作按钮样式 === */
        .btn {{
            background: var(--primary-color);
            color: white;
            border: none;
            padding: 0.5rem 1rem;
            border-radius: 5px;
            cursor: pointer;
            font-size: 0.9rem;
            transition: all 0.3s;
        }}
        
        .btn:hover {{
            background: var(--dark-color);
        }}
        
        .table-actions {{
            margin-top: 1rem;
            display: flex;
            gap: 1rem;
        }}
        
        /* === 图表展示区域：网格布局 === */
        .plot-grid {{
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(400px, 1fr));
            gap: 1.5rem;
        }}
        
        .plot-card {{
            border: 1px solid var(--border-color);
            border-radius: 10px;
            padding: 1rem;
            text-align: center;
        }}
        
        .plot-card h4 {{
            margin-bottom: 1rem;
            color: var(--dark-color);
        }}
        
        .plot-placeholder {{
            background: var(--light-color);
            padding: 3rem;
            border-radius: 5px;
        }}
        
        .plot-placeholder a {{
            color: var(--primary-color);
            text-decoration: none;
            display: flex;
            flex-direction: column;
            align-items: center;
            gap: 0.5rem;
        }}
        
        .plot-placeholder i {{
            font-size: 3rem;
        }}
        
        /* === AI解读区域样式：左侧绿色边框标识 === */
        .ai-interpretation {{
            background: #f8f9fa;
            border-left: 4px solid var(--secondary-color);
            padding: 1.5rem;
            margin-top: 1rem;
            border-radius: 0 10px 10px 0;
        }}
        
        .ai-interpretation h3 {{
            color: var(--secondary-color);
            margin-bottom: 1rem;
        }}
        
        .interpretation-content {{
            line-height: 1.8;
        }}
        
        .ai-disclaimer {{
            background: #fff3cd;
            padding: 1rem;
            border-radius: 5px;
            margin-bottom: 1rem;
        }}
        
        /* === 页脚样式 === */
        .footer {{
            text-align: center;
            padding: 2rem;
            color: #666;
            font-size: 0.9rem;
        }}
        
        /* === 响应式布局：适配移动端和小屏幕设备 === */
        @media (max-width: 768px) {{
            .header h1 {{
                font-size: 1.8rem;
            }}
            
            .nav ul {{
                flex-direction: column;
                align-items: center;
            }}
            
            .summary-grid {{
                grid-template-columns: repeat(2, 1fr);
            }}
            
            .plot-grid {{
                grid-template-columns: 1fr;
            }}
        }}
    </style>
</head>
<body>
    <!-- 页面头部：显示报告标题、版本号和生成时间 -->
    <header class="header">
        <h1><i class="fas fa-dna"></i> AllEnricher Report</h1>
        <p class="version">Version 2.0 | Generated on {datetime.now().strftime("%Y-%m-%d %H:%M")}</p>
    </header>
    
    <!-- 导航栏：包含各部分的锚点链接，支持快速跳转 -->
    <nav class="nav">
        <ul>
            <li><a href="#summary"><i class="fas fa-chart-pie"></i> Summary</a></li>
            <li><a href="#plots"><i class="fas fa-chart-bar"></i> Plots</a></li>
            {" ".join([f'<li><a href="#{db}-table"><i class="fas fa-table"></i> {db}</a></li>' for db in (db_names or [])])}
            {f'<li><a href="#ai-interpretation"><i class="fas fa-brain"></i> AI Interpretation</a></li>' if ai_section else ''}
        </ul>
    </nav>
    
    <!-- 主内容区域：依次展示摘要、图表、数据表格和AI解读 -->
    <main class="main">
        {summary}
        {plots}
        {tables}
        {ai_section}
    </main>
    
    <!-- 页脚：版权信息 -->
    <footer class="footer">
        <p>Generated by <a href="https://github.com/zd105/AllEnricher" target="_blank">AllEnricher v2.0</a></p>
        <p>&copy; {datetime.now().year} AllEnricher Team. Licensed under MIT.</p>
    </footer>
    
    <!-- jQuery 库：DataTables 插件的前置依赖 -->
    <script src="https://code.jquery.com/jquery-3.6.0.min.js"></script>
    
    <!-- DataTables 插件：为HTML表格添加排序、分页、搜索等交互功能 -->
    <script src="https://cdn.datatables.net/1.13.4/js/jquery.dataTables.min.js"></script>
    
    <script>
        // 初始化 DataTables 插件，为所有数据表格启用交互功能
        $(document).ready(function() {{
            $('.data-table').DataTable({{
                pageLength: 25,           // 每页显示25条记录
                order: [[5, 'asc']],      // 默认按第6列（校正后P值）升序排列
                responsive: true          // 启用响应式布局
            }});
        }});
        
        // 下载表格数据为TSV文件
        // 从表格DOM中提取表头和行数据，组装为制表符分隔的文本并触发下载
        function downloadTable(dbName) {{
            const table = document.getElementById('table-' + dbName);
            let tsv = [];
            
            // 提取表头文本
            const headers = Array.from(table.querySelectorAll('th')).map(th => th.textContent);
            tsv.push(headers.join('\\t'));
            
            // 提取每行数据
            const rows = table.querySelectorAll('tbody tr');
            rows.forEach(row => {{
                const cells = Array.from(row.querySelectorAll('td')).map(td => td.textContent);
                tsv.push(cells.join('\\t'));
            }});
            
            // 创建Blob对象并触发文件下载
            const blob = new Blob([tsv.join('\\n')], {{ type: 'text/tab-separated-values' }});
            const url = URL.createObjectURL(blob);
            const a = document.createElement('a');
            a.href = url;
            a.download = dbName + '_enrichment.tsv';  // 文件名格式：数据库名_enrichment.tsv
            a.click();
            URL.revokeObjectURL(url);  // 释放临时URL资源
        }}
        
        // 复制表格内容到系统剪贴板
        // 通过选中表格DOM节点并执行复制命令实现
        function copyTable(dbName) {{
            const table = document.getElementById('table-' + dbName);
            const range = document.createRange();
            range.selectNode(table);
            window.getSelection().removeAllRanges();   // 清除已有选区
            window.getSelection().addRange(range);     // 选中表格内容
            document.execCommand('copy');               // 执行复制命令
            window.getSelection().removeAllRanges();   // 清除选区
            alert('Table copied to clipboard!');        // 提示用户复制成功
        }}
    </script>
</body>
</html>
'''
        
        return html
