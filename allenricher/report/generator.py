"""
Interactive HTML Report Generator for AllEnricher v2.0

生成交互式HTML报告模块 - 学术风格设计

采用低调、专业的学术审美风格，符合研究分析报告的视觉标准。
设计参考：Nature/Science 期刊 supplementary data 的排版风格。
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
    交互式HTML报告生成器 - 学术风格

    设计特点：
    - 中性色调：使用灰度色系，避免鲜艳颜色
    - 专业字体：Noto Serif + Source Sans Pro 组合
    - 充足留白：类似学术论文的排版密度
    - 清晰层级：通过字体大小和字重区分信息重要性
    """

    def __init__(self, output_dir: str, config=None):
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.config = config

    def generate(
        self,
        results: Dict[str, pd.DataFrame],
        output_file: str,
        gene_list: List[str] = None,
        ai_interpretation: Dict[str, str] = None
    ) -> str:
        has_results = results and any(len(df) > 0 for df in results.values())

        if not has_results:
            html = self._generate_no_results_page(gene_list)
            with open(output_file, 'w', encoding='utf-8') as f:
                f.write(html)
            return output_file

        summary = self._generate_summary(results, gene_list)
        tables = self._generate_tables(results)
        plots = self._generate_plot_section(results)
        ai_section = self._generate_ai_section(ai_interpretation) if ai_interpretation else ""

        active_db_names = [db for db, df in results.items() if len(df) > 0]
        html = self._build_html(
            summary=summary,
            tables=tables,
            plots=plots,
            ai_section=ai_section,
            db_names=active_db_names
        )

        with open(output_file, 'w', encoding='utf-8') as f:
            f.write(html)

        return output_file

    def _generate_no_results_page(self, gene_list: List[str] = None) -> str:
        """生成无富集结果提示页面"""
        html = f'''<!DOCTYPE html>
<html lang="zh-CN">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>AllEnricher v2.0 - 无富集结果</title>
    <link href="https://fonts.googleapis.com/css2?family=Noto+Serif:wght@400;600&family=Source+Sans+Pro:wght@400;600&display=swap" rel="stylesheet">
    <style>
        :root {{
            --text-primary: #1a1a1a;
            --text-secondary: #4a4a4a;
            --text-muted: #6b6b6b;
            --border-color: #d1d1d1;
            --bg-primary: #ffffff;
            --bg-secondary: #f8f8f8;
            --accent-color: #2c5282;
        }}
        * {{ margin: 0; padding: 0; box-sizing: border-box; }}
        body {{
            font-family: 'Source Sans Pro', -apple-system, sans-serif;
            line-height: 1.6;
            color: var(--text-primary);
            background-color: var(--bg-secondary);
        }}
        .header {{
            background: var(--bg-primary);
            border-bottom: 1px solid var(--border-color);
            padding: 1.5rem 2rem;
        }}
        .header h1 {{
            font-family: 'Noto Serif', Georgia, serif;
            font-size: 1.5rem;
            font-weight: 600;
        }}
        .meta {{
            font-size: 0.875rem;
            color: var(--text-muted);
            margin-top: 0.25rem;
        }}
        .main {{
            max-width: 700px;
            margin: 3rem auto;
            padding: 0 2rem;
        }}
        .no-results-box {{
            background: var(--bg-primary);
            border: 1px solid var(--border-color);
            border-radius: 4px;
            padding: 2.5rem;
        }}
        .no-results-box h2 {{
            font-family: 'Noto Serif', Georgia, serif;
            font-size: 1.25rem;
            margin-bottom: 1rem;
            color: var(--text-primary);
        }}
        .info-box {{
            background: var(--bg-secondary);
            border: 1px solid var(--border-color);
            padding: 1.5rem;
            margin-top: 1.5rem;
        }}
        .info-box h3 {{
            font-size: 0.9rem;
            font-weight: 600;
            margin-bottom: 0.75rem;
            color: var(--text-secondary);
        }}
        .info-box ul {{
            margin-left: 1.25rem;
            color: var(--text-secondary);
            font-size: 0.9rem;
        }}
        .info-box li {{
            margin-bottom: 0.4rem;
        }}
        .footer {{
            text-align: center;
            padding: 2rem;
            font-size: 0.8rem;
            color: var(--text-muted);
        }}
        .footer a {{
            color: var(--accent-color);
        }}
    </style>
</head>
<body>
    <header class="header">
        <h1>AllEnricher Report</h1>
        <p class="meta">Version 2.0 | {datetime.now().strftime("%Y-%m-%d")}</p>
    </header>
    <main class="main">
        <div class="no-results-box">
            <h2>未找到显著富集的结果</h2>
            <p style="color: var(--text-secondary);">
                根据当前的分析参数，未能检测到统计上显著富集的功能条目。
            </p>
            <div class="info-box">
                <h3>可能的原因</h3>
                <ul>
                    <li>输入基因列表过小或与数据库无交集</li>
                    <li>P值/Q值阈值过于严格</li>
                    <li>背景基因集设置不当</li>
                    <li>基因ID格式与数据库不匹配</li>
                </ul>
            </div>
            <div class="info-box">
                <h3>建议</h3>
                <ul>
                    <li>增加输入基因数量（建议至少10个以上）</li>
                    <li>放宽P值/Q值阈值（如 <code>-p 0.1 -q 0.1</code>）</li>
                    <li>检查基因ID格式是否正确</li>
                    <li>尝试使用其他数据库进行分析</li>
                </ul>
            </div>
        </div>
    </main>
    <footer class="footer">
        <p>Generated by <a href="https://github.com/zd105/AllEnricher">AllEnricher v2.0</a></p>
    </footer>
</body>
</html>'''
        return html

    def _generate_summary(self, results: Dict[str, pd.DataFrame], gene_list: List[str] = None) -> str:
        """生成统计摘要部分"""
        total_terms = sum(len(df) for df in results.values())
        databases = list(results.keys())

        summary_stats = []
        for db_name, df in results.items():
            if len(df) > 0:
                summary_stats.append({
                    "database": db_name,
                    "terms": len(df),
                    "min_pval": df['P_Value'].min() if 'P_Value' in df.columns else 0,
                    "min_adj_pval": df['Adjusted_P_Value'].min() if 'Adjusted_P_Value' in df.columns else 0
                })

        rows_html = "".join([
            f'<tr><td><a href="#{s["database"]}-table">{s["database"]}</a></td>'
            f'<td>{s["terms"]}</td>'
            f'<td>{s["min_pval"]:.2e}</td>'
            f'<td>{s["min_adj_pval"]:.2e}</td></tr>'
            for s in summary_stats
        ])

        html = f'''
        <div class="section" id="summary">
            <h2>Analysis Summary</h2>
            <div class="summary-grid">
                <div class="stat-item">
                    <span class="stat-value">{len(gene_list) if gene_list else "—"}</span>
                    <span class="stat-label">Input Genes</span>
                </div>
                <div class="stat-item">
                    <span class="stat-value">{len(databases)}</span>
                    <span class="stat-label">Databases</span>
                </div>
                <div class="stat-item">
                    <span class="stat-value">{total_terms}</span>
                    <span class="stat-label">Enriched Terms</span>
                </div>
            </div>
            <table class="data-table">
                <thead>
                    <tr>
                        <th>Database</th>
                        <th>Terms</th>
                        <th>Min P-value</th>
                        <th>Min Adj. P-value</th>
                    </tr>
                </thead>
                <tbody>{rows_html}</tbody>
            </table>
        </div>'''
        return html

    def _generate_tables(self, results: Dict[str, pd.DataFrame]) -> str:
        """生成交互式数据表格"""
        tables_html = []

        for db_name, df in results.items():
            if len(df) == 0:
                continue

            rows = []
            for idx, row in df.iterrows():
                tid = row.get('Term_ID', 'N/A')
                tname = row.get('Term_Name', 'N/A')
                gcount = row.get('Gene_Count', 0)
                rf = f"{row.get('Rich_Factor', 0):.4f}"
                pv = f"{row.get('P_Value', 1):.2e}"
                adjpv = f"{row.get('Adjusted_P_Value', 1):.2e}"
                term_url = row.get('Term_URL', '')
                genes_str = str(row.get('Genes', ''))
                genes_short = genes_str[:60] + ('...' if len(genes_str) > 60 else '')

                if term_url:
                    tid_cell = f'<a href="{term_url}" target="_blank">{tid}</a>'
                else:
                    tid_cell = tid

                rows.append(f'<tr><td>{tid_cell}</td><td>{tname}</td><td>{gcount}</td>'
                            f'<td>{rf}</td><td>{pv}</td><td>{adjpv}</td>'
                            f'<td class="genes" title="{genes_str}">{genes_short}</td></tr>')

            table_html = f'''
            <div class="section" id="{db_name}-table">
                <h2>{db_name} Enrichment Results <span class="result-count">({len(df)} terms)</span></h2>
                <div class="table-wrapper">
                    <table class="data-table" id="table-{db_name}">
                        <thead>
                            <tr>
                                <th>Term ID</th>
                                <th>Term Name</th>
                                <th>Genes</th>
                                <th>Rich Factor</th>
                                <th>P-value</th>
                                <th>Adj. P-value</th>
                            </tr>
                        </thead>
                        <tbody>{"".join(rows)}</tbody>
                    </table>
                </div>
                <div class="table-actions">
                    <button onclick="downloadTable('{db_name}')">Download TSV</button>
                    <button onclick="copyTable('{db_name}')">Copy</button>
                </div>
            </div>'''
            tables_html.append(table_html)

        return "\n".join(tables_html)

    def _generate_plot_section(self, results: Dict[str, pd.DataFrame]) -> str:
        """生成图表展示区域 - 使用PNG图片直接嵌入HTML"""
        plots_html = ['<div class="section" id="plots"><h2>Visualization</h2>']

        for db_name in results.keys():
            plot_dir = self.output_dir / "plots"
            
            # 优先使用PNG，回退到PDF链接
            barplot_png = plot_dir / f"{db_name}_barplot.png"
            barplot_pdf = plot_dir / f"{db_name}_barplot.pdf"
            bubble_png = plot_dir / f"{db_name}_bubble.png"
            bubble_pdf = plot_dir / f"{db_name}_bubble.pdf"

            has_bar = barplot_png.exists() or barplot_pdf.exists()
            has_bubble = bubble_png.exists() or bubble_pdf.exists()

            if has_bar or has_bubble:
                plots_html.append(f'<div class="plot-group"><h3>{db_name}</h3>')
                
                # Bar Plot
                if barplot_png.exists():
                    # 使用PNG图片
                    img_data = self._encode_image_to_base64(barplot_png)
                    plots_html.append(f'''
                        <div class="plot-container">
                            <img src="data:image/png;base64,{img_data}" alt="{db_name} Bar Plot" class="plot-img">
                            <p class="plot-caption">Bar Plot (Top enriched terms by Q-value)</p>
                        </div>
                    ''')
                elif barplot_pdf.exists():
                    # 回退到PDF链接
                    plots_html.append(f'<a href="plots/{barplot_pdf.name}" target="_blank" class="plot-link">Bar Plot (PDF)</a>')
                
                # Bubble Plot
                if bubble_png.exists():
                    # 使用PNG图片
                    img_data = self._encode_image_to_base64(bubble_png)
                    plots_html.append(f'''
                        <div class="plot-container">
                            <img src="data:image/png;base64,{img_data}" alt="{db_name} Bubble Plot" class="plot-img">
                            <p class="plot-caption">Bubble Plot (Gene count vs Rich factor)</p>
                        </div>
                    ''')
                elif bubble_pdf.exists():
                    # 回退到PDF链接
                    plots_html.append(f'<a href="plots/{bubble_pdf.name}" target="_blank" class="plot-link">Bubble Plot (PDF)</a>')
                
                plots_html.append('</div>')

        plots_html.append('</div>')
        return "\n".join(plots_html)

    def _encode_image_to_base64(self, image_path: Path) -> str:
        """将图片文件编码为base64字符串"""
        try:
            with open(image_path, 'rb') as f:
                return base64.b64encode(f.read()).decode('utf-8')
        except Exception as e:
            logger.warning(f"图片编码失败: {image_path}, {e}")
            return ""

    def _generate_ai_section(self, ai_interpretation: Dict[str, str]) -> str:
        """生成AI解读部分"""
        if not ai_interpretation:
            return ""

        interpretations = []
        for db_name, interpretation in ai_interpretation.items():
            # 将 Markdown 格式转换为 HTML：\n 转为 <br>，**bold** 转为 <strong>bold</strong>
            html_text = interpretation.replace('\n', '<br>')
            import re
            # 使用正则替换成对的 ** 为 <strong>...</strong>
            html_text = re.sub(r'\*\*(.+?)\*\*', r'<strong>\1</strong>', html_text)
            interpretations.append(f'<div class="ai-block"><h3>{db_name}</h3><p>{html_text}</p></div>')

        html = f'''
        <div class="section" id="ai-interpretation">
            <h2>AI Interpretation <span class="ai-note">(Review recommended)</span></h2>
            {"".join(interpretations)}
            <div class="ai-disclaimer">
                <strong>Disclaimer:</strong> This AI-generated interpretation is for reference only and should be reviewed by domain experts. Verify all biological conclusions with literature.
            </div>
        </div>'''
        return html

    def _build_html(self, summary: str, tables: str, plots: str, ai_section: str, db_names: List[str] = None) -> str:
        """构建完整HTML文档"""

        nav_items = '<li><a href="#summary">Summary</a></li>'
        nav_items += '<li><a href="#plots">Plots</a></li>'
        if db_names:
            nav_items += ''.join([f'<li><a href="#{db}-table">{db}</a></li>' for db in db_names])
        if ai_section:
            nav_items += '<li><a href="#ai-interpretation">AI</a></li>'

        html = f'''<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>AllEnricher v2.0 - Enrichment Analysis Report</title>
    <link href="https://fonts.googleapis.com/css2?family=Noto+Serif:wght@400;600&family=Source+Sans+Pro:wght@400;600&display=swap" rel="stylesheet">
    <style>
        :root {{
            --text-primary: #1a1a1a;
            --text-secondary: #4a4a4a;
            --text-muted: #6b6b6b;
            --border-color: #d1d1d1;
            --bg-primary: #ffffff;
            --bg-secondary: #f8f8f8;
            --accent-color: #2c5282;
            --accent-hover: #1a365d;
        }}

        * {{ margin: 0; padding: 0; box-sizing: border-box; }}

        body {{
            font-family: 'Source Sans Pro', -apple-system, BlinkMacSystemFont, sans-serif;
            font-size: 14px;
            line-height: 1.6;
            color: var(--text-primary);
            background-color: var(--bg-secondary);
        }}

        /* Header */
        .header {{
            background: var(--bg-primary);
            border-bottom: 1px solid var(--border-color);
            padding: 1.25rem 2rem;
        }}

        .header-content {{
            max-width: 1200px;
            margin: 0 auto;
            display: flex;
            justify-content: space-between;
            align-items: baseline;
        }}

        .header h1 {{
            font-family: 'Noto Serif', Georgia, serif;
            font-size: 1.25rem;
            font-weight: 600;
            color: var(--text-primary);
        }}

        .header .meta {{
            font-size: 0.8rem;
            color: var(--text-muted);
        }}

        /* Navigation */
        .nav {{
            background: var(--bg-primary);
            border-bottom: 1px solid var(--border-color);
            position: sticky;
            top: 0;
            z-index: 100;
        }}

        .nav ul {{
            max-width: 1200px;
            margin: 0 auto;
            padding: 0 2rem;
            list-style: none;
            display: flex;
            gap: 2rem;
        }}

        .nav a {{
            display: block;
            padding: 0.75rem 0;
            color: var(--text-secondary);
            text-decoration: none;
            font-size: 0.875rem;
            border-bottom: 2px solid transparent;
            transition: border-color 0.15s;
        }}

        .nav a:hover {{
            color: var(--accent-color);
            border-bottom-color: var(--accent-color);
        }}

        /* Main Content */
        .main {{
            max-width: 1200px;
            margin: 0 auto;
            padding: 2rem;
        }}

        /* Section */
        .section {{
            background: var(--bg-primary);
            border: 1px solid var(--border-color);
            border-radius: 4px;
            padding: 1.5rem 2rem;
            margin-bottom: 1.5rem;
        }}

        .section h2 {{
            font-family: 'Noto Serif', Georgia, serif;
            font-size: 1.1rem;
            font-weight: 600;
            color: var(--text-primary);
            margin-bottom: 1.25rem;
            padding-bottom: 0.5rem;
            border-bottom: 1px solid var(--border-color);
        }}

        .section h2 .result-count,
        .section h2 .ai-note {{
            font-family: 'Source Sans Pro', sans-serif;
            font-weight: 400;
            font-size: 0.875rem;
            color: var(--text-muted);
        }}

        /* Summary Grid */
        .summary-grid {{
            display: flex;
            gap: 2rem;
            margin-bottom: 1.5rem;
            padding: 1rem 0;
        }}

        .stat-item {{
            display: flex;
            flex-direction: column;
        }}

        .stat-value {{
            font-size: 1.5rem;
            font-weight: 600;
            color: var(--text-primary);
            line-height: 1.2;
        }}

        .stat-label {{
            font-size: 0.8rem;
            color: var(--text-muted);
            margin-top: 0.25rem;
        }}

        /* Data Table */
        .table-wrapper {{
            overflow-x: auto;
            margin: 0 -2rem;
            padding: 0 2rem;
        }}

        .data-table {{
            width: 100%;
            border-collapse: collapse;
            font-size: 0.875rem;
        }}

        .data-table th {{
            text-align: left;
            padding: 0.6rem 0.75rem;
            background: var(--bg-secondary);
            border-bottom: 2px solid var(--border-color);
            font-weight: 600;
            color: var(--text-secondary);
            white-space: nowrap;
        }}

        .data-table td {{
            padding: 0.6rem 0.75rem;
            border-bottom: 1px solid var(--border-color);
            vertical-align: top;
        }}

        .data-table tr:hover td {{
            background: #fafafa;
        }}

        .data-table a {{
            color: var(--accent-color);
            text-decoration: none;
        }}

        .data-table a:hover {{
            text-decoration: underline;
        }}

        .data-table .genes {{
            font-size: 0.8rem;
            color: var(--text-muted);
            max-width: 200px;
            overflow: hidden;
            text-overflow: ellipsis;
            white-space: nowrap;
        }}

        /* Table Actions */
        .table-actions {{
            display: flex;
            gap: 1rem;
            margin-top: 1rem;
            padding-top: 1rem;
            border-top: 1px solid var(--border-color);
        }}

        .table-actions button {{
            background: var(--bg-secondary);
            border: 1px solid var(--border-color);
            padding: 0.4rem 0.75rem;
            font-size: 0.8rem;
            color: var(--text-secondary);
            cursor: pointer;
            border-radius: 3px;
            transition: all 0.15s;
        }}

        .table-actions button:hover {{
            background: var(--bg-primary);
            border-color: var(--accent-color);
            color: var(--accent-color);
        }}

        /* Plots */
        .plot-group {{
            margin-bottom: 1.25rem;
        }}

        .plot-group h3 {{
            font-size: 0.9rem;
            font-weight: 600;
            color: var(--text-secondary);
            margin-bottom: 0.5rem;
        }}

        .plot-links {{
            display: flex;
            gap: 1rem;
            font-size: 0.875rem;
        }}

        .plot-link {{
            color: var(--accent-color);
            text-decoration: none;
        }}

        .plot-link:hover {{
            text-decoration: underline;
        }}

        /* Plot Images - PNG嵌入样式 */
        .plot-container {{
            margin: 1rem 0;
            padding: 1rem;
            background: var(--bg-secondary);
            border: 1px solid var(--border-color);
            border-radius: 4px;
            text-align: center;
        }}

        .plot-img {{
            max-width: 100%;
            height: auto;
            max-height: 600px;
            border: 1px solid var(--border-color);
            background: white;
        }}

        .plot-caption {{
            font-size: 0.8rem;
            color: var(--text-muted);
            margin-top: 0.5rem;
            font-style: italic;
        }}

        /* AI Interpretation */
        .ai-block {{
            margin-bottom: 1.25rem;
            padding: 1rem;
            background: var(--bg-secondary);
            border-left: 3px solid var(--border-color);
        }}

        .ai-block h3 {{
            font-size: 0.9rem;
            font-weight: 600;
            color: var(--text-secondary);
            margin-bottom: 0.5rem;
        }}

        .ai-block p {{
            font-size: 0.875rem;
            color: var(--text-secondary);
            line-height: 1.7;
        }}

        .ai-disclaimer {{
            margin-top: 1rem;
            padding: 0.75rem;
            background: #fff3cd;
            border: 1px solid #ffc107;
            border-radius: 4px;
            font-size: 0.8rem;
            color: #856404;
        }}

        /* Footer */
        .footer {{
            text-align: center;
            padding: 2rem;
            font-size: 0.8rem;
            color: var(--text-muted);
            border-top: 1px solid var(--border-color);
            margin-top: 2rem;
        }}

        .footer a {{
            color: var(--accent-color);
            text-decoration: none;
        }}

        /* DataTables customization */
        .dataTables_wrapper .dataTables_length,
        .dataTables_wrapper .dataTables_filter {{
            margin-bottom: 1rem;
        }}

        .dataTables_wrapper .dataTables_info {{
            color: var(--text-muted);
            font-size: 0.8rem;
        }}

        .paginate_button {{
            color: var(--text-secondary) !important;
        }}

        .paginate_button.current {{
            color: var(--accent-color) !important;
            font-weight: 600;
        }}

        /* Responsive */
        @media (max-width: 768px) {{
            .header-content {{
                flex-direction: column;
                gap: 0.5rem;
            }}

            .nav ul {{
                overflow-x: auto;
                gap: 1.5rem;
            }}

            .summary-grid {{
                flex-wrap: wrap;
                gap: 1rem;
            }}

            .stat-item {{
                min-width: 80px;
            }}
        }}
    </style>
</head>
<body>
    <header class="header">
        <div class="header-content">
            <h1>AllEnricher Report</h1>
            <span class="meta">Version 2.0 | {datetime.now().strftime("%Y-%m-%d")}</span>
        </div>
    </header>

    <nav class="nav">
        <ul>{nav_items}</ul>
    </nav>

    <main class="main">
        {summary}
        {plots}
        {tables}
        {ai_section}
    </main>

    <footer class="footer">
        <p>Generated by <a href="https://github.com/zd105/AllEnricher">AllEnricher v2.0</a></p>
    </footer>

    <script src="https://code.jquery.com/jquery-3.6.0.min.js"></script>
    <script src="https://cdn.datatables.net/1.13.4/js/jquery.dataTables.min.js"></script>
    <script>
        $(document).ready(function() {{
            $('.data-table').DataTable({{
                pageLength: 25,
                order: [[5, 'asc']],
                columnDefs: [
                    {{ orderable: false, targets: [2, 6] }}
                ]
            }});
        }});

        function downloadTable(dbName) {{
            const table = document.getElementById('table-' + dbName);
            let tsv = [];
            const headers = Array.from(table.querySelectorAll('th')).map(th => th.textContent.trim());
            tsv.push(headers.join('\t'));
            const rows = table.querySelectorAll('tbody tr');
            rows.forEach(row => {{
                const cells = Array.from(row.querySelectorAll('td')).map(td => td.textContent.trim());
                tsv.push(cells.join('\t'));
            }});
            const blob = new Blob([tsv.join('\n')], {{ type: 'text/tab-separated-values' }});
            const url = URL.createObjectURL(blob);
            const a = document.createElement('a');
            a.href = url;
            a.download = dbName + '_enrichment.tsv';
            a.click();
            URL.revokeObjectURL(url);
        }}

        function copyTable(dbName) {{
            const table = document.getElementById('table-' + dbName);
            const range = document.createRange();
            range.selectNode(table);
            window.getSelection().removeAllRanges();
            window.getSelection().addRange(range);
            document.execCommand('copy');
            window.getSelection().removeAllRanges();
            alert('Table copied to clipboard');
        }}
    </script>
</body>
</html>'''
        return html
