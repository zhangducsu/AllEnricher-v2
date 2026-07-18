"""Generate self-contained, publication-oriented HTML reports for AllEnricher analyses."""

import os
import json
import logging
import re
import tempfile
import shutil
from html import escape
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional, Any
import pandas as pd
import base64
import allenricher
from allenricher.core.config import database_display_name
from allenricher.visualization.color_config import PaletteLike, resolve_palette_selection
from allenricher.report.methods_reference import render_methods_reference_html

logger = logging.getLogger(__name__)


def _first_value_for_report(row: pd.Series, names: tuple[str, ...], default: Any = None) -> Any:
    for name in names:
        if name in row.index and pd.notna(row.get(name)) and str(row.get(name)).strip():
            return row.get(name)
    return default


class ReportGenerator:
    """Build searchable result tables, figures, provenance, and Methods text into one HTML report."""

    _TABLE_JS = '''
    <script>
        // Vanilla JS Table Component
        (function() {
            const tables = document.querySelectorAll('.data-table[id^="table-"]');

            tables.forEach(table => {
                const tbody = table.querySelector('tbody');
                const rows = Array.from(tbody.querySelectorAll('tr'));
                const headers = Array.from(table.querySelectorAll('thead th'));
                const dbName = table.id.replace('table-', '');

                // State
                let currentPage = 1;
                let pageSize = 20;
                let sortColumn = Math.max(0, headers.findIndex(th =>
                    ['adj. p-value', 'padj', 'fdr', 'p-value', 'pval'].includes(th.textContent.trim().toLowerCase())
                ));
                let sortDirection = 'asc';
                let searchTerm = '';
                let filteredRows = [...rows];

                // Make headers sortable
                headers.forEach((th, index) => {
                    th.classList.add('sortable');
                    th.addEventListener('click', () => {
                        if (sortColumn === index) {
                            sortDirection = sortDirection === 'asc' ? 'desc' : 'asc';
                        } else {
                            sortColumn = index;
                            sortDirection = 'asc';
                        }
                        updateSortIcons();
                        applyFilters();
                    });
                });

                function updateSortIcons() {
                    headers.forEach((th, index) => {
                        th.classList.remove('sort-asc', 'sort-desc');
                        if (index === sortColumn) {
                            th.classList.add(sortDirection === 'asc' ? 'sort-asc' : 'sort-desc');
                        }
                    });
                }

                // Create controls container
                const controlsDiv = document.createElement('div');
                controlsDiv.className = 'table-controls';
                controlsDiv.innerHTML = `
                    <div class="table-length">
                        <label>Show</label>
                        <select class="page-size-select">
                            <option value="10" ${pageSize === 10 ? 'selected' : ''}>10</option>
                            <option value="20" ${pageSize === 20 ? 'selected' : ''}>20</option>
                            <option value="50" ${pageSize === 50 ? 'selected' : ''}>50</option>
                            <option value="100" ${pageSize === 100 ? 'selected' : ''}>100</option>
                        </select>
                        <label>entries</label>
                    </div>
                    <div class="table-search">
                        <label>Search:</label>
                        <input type="text" class="search-input" placeholder="Filter results..." value="${searchTerm}">
                    </div>
                `;
                table.parentNode.insertBefore(controlsDiv, table);

                // Create info and pagination container
                const footerDiv = document.createElement('div');
                footerDiv.innerHTML = `
                    <div class="table-info"></div>
                    <div class="pagination"></div>
                `;
                table.parentNode.insertBefore(footerDiv, table.nextSibling);

                const pageSizeSelect = controlsDiv.querySelector('.page-size-select');
                const searchInput = controlsDiv.querySelector('.search-input');
                const infoDiv = footerDiv.querySelector('.table-info');
                const paginationDiv = footerDiv.querySelector('.pagination');

                // Event listeners
                pageSizeSelect.addEventListener('change', (e) => {
                    pageSize = parseInt(e.target.value);
                    currentPage = 1;
                    applyFilters();
                });

                searchInput.addEventListener('input', (e) => {
                    searchTerm = e.target.value.toLowerCase();
                    currentPage = 1;
                    applyFilters();
                });

                function filterRows() {
                    if (!searchTerm) return [...rows];
                    return rows.filter(row => {
                        const cells = row.querySelectorAll('td');
                        return Array.from(cells).some(cell =>
                            cell.textContent.toLowerCase().includes(searchTerm)
                        );
                    });
                }

                function sortRows(rowsToSort) {
                    const sorted = [...rowsToSort].sort((a, b) => {
                        const aVal = a.cells[sortColumn]?.textContent.trim() || '';
                        const bVal = b.cells[sortColumn]?.textContent.trim() || '';

                        // Try numeric comparison
                        const aNum = parseFloat(aVal);
                        const bNum = parseFloat(bVal);

                        if (!isNaN(aNum) && !isNaN(bNum)) {
                            return sortDirection === 'asc' ? aNum - bNum : bNum - aNum;
                        }

                        // String comparison
                        return sortDirection === 'asc'
                            ? aVal.localeCompare(bVal)
                            : bVal.localeCompare(aVal);
                    });
                    // Re-append rows in sorted order to update DOM
                    sorted.forEach(row => tbody.appendChild(row));
                    return sorted;
                }

                function renderPagination() {
                    const totalPages = Math.ceil(filteredRows.length / pageSize);
                    let html = '';

                    // Previous button
                    html += `<button ${currentPage === 1 ? 'disabled' : ''} onclick="goToPage('${dbName}', ${currentPage - 1})">Previous</button>`;

                    // Page numbers
                    const maxButtons = 5;
                    let startPage = Math.max(1, currentPage - Math.floor(maxButtons / 2));
                    let endPage = Math.min(totalPages, startPage + maxButtons - 1);

                    if (endPage - startPage < maxButtons - 1) {
                        startPage = Math.max(1, endPage - maxButtons + 1);
                    }

                    if (startPage > 1) {
                        html += `<button onclick="goToPage('${dbName}', 1)">1</button>`;
                        if (startPage > 2) html += `<button disabled>...</button>`;
                    }

                    for (let i = startPage; i <= endPage; i++) {
                        html += `<button class="${i === currentPage ? 'active' : ''}" onclick="goToPage('${dbName}', ${i})">${i}</button>`;
                    }

                    if (endPage < totalPages) {
                        if (endPage < totalPages - 1) html += `<button disabled>...</button>`;
                        html += `<button onclick="goToPage('${dbName}', ${totalPages})">${totalPages}</button>`;
                    }

                    // Next button
                    html += `<button ${currentPage === totalPages || totalPages === 0 ? 'disabled' : ''} onclick="goToPage('${dbName}', ${currentPage + 1})">Next</button>`;

                    paginationDiv.innerHTML = html;
                }

                function renderTable() {
                    const start = (currentPage - 1) * pageSize;
                    const end = start + pageSize;
                    const pageRows = filteredRows.slice(start, end);

                    // Hide all rows first
                    rows.forEach(row => row.style.display = 'none');

                    // Show only page rows
                    pageRows.forEach(row => row.style.display = '');

                    // Update info
                    const total = filteredRows.length;
                    const showingStart = total === 0 ? 0 : start + 1;
                    const showingEnd = Math.min(end, total);
                    infoDiv.textContent = `Showing ${showingStart} to ${showingEnd} of ${total} entries${searchTerm ? ' (filtered from ' + rows.length + ' total entries)' : ''}`;

                    renderPagination();
                }

                function applyFilters() {
                    filteredRows = filterRows();
                    filteredRows = sortRows(filteredRows);
                    renderTable();
                }

                // Store table state globally for pagination callbacks
                if (!window.tableStates) window.tableStates = {};
                window.tableStates[dbName] = {
                    setPage: (page) => {
                        currentPage = page;
                        renderTable();
                    }
                };

                // Initial render
                updateSortIcons();
                applyFilters();
            });

            // Global pagination function
            window.goToPage = function(dbName, page) {
                if (window.tableStates && window.tableStates[dbName]) {
                    window.tableStates[dbName].setPage(page);
                }
            };
        })();

        function exportCellText(cell) {
            const fullValue = cell.querySelector('[data-full]');
            return fullValue ? fullValue.dataset.full : cell.textContent.trim();
        }

        function downloadTable(dbName) {
            const table = document.getElementById('table-' + dbName);
            let tsv = [];
            const headers = Array.from(table.querySelectorAll('th')).map(th => th.textContent.trim());
            tsv.push(headers.join('\\t'));
            const rows = table.querySelectorAll('tbody tr');
            rows.forEach(row => {
                const cells = Array.from(row.querySelectorAll('td')).map(exportCellText);
                tsv.push(cells.join('\\t'));
            });
            const blob = new Blob([tsv.join('\\n')], { type: 'text/tab-separated-values' });
            const url = URL.createObjectURL(blob);
            const a = document.createElement('a');
            a.href = url;
            a.download = dbName + '_enrichment.tsv';
            a.click();
            URL.revokeObjectURL(url);
        }

        function copyTable(dbName) {
            const table = document.getElementById('table-' + dbName);
            const clone = table.cloneNode(true);
            clone.querySelectorAll('tbody tr').forEach(row => row.style.display = '');
            clone.querySelectorAll('[data-full]').forEach(value => {
                value.textContent = value.dataset.full;
            });

            const range = document.createRange();
            range.selectNode(clone);
            window.getSelection().removeAllRanges();
            window.getSelection().addRange(range);
            document.execCommand('copy');
            window.getSelection().removeAllRanges();
            alert('Table copied to clipboard');
        }
    </script>'''

    def __init__(self, output_dir: str, config=None):
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.config = config

    @staticmethod
    def _html_id(value: object) -> str:
        """Return a stable identifier that is safe in HTML and JavaScript."""
        slug = re.sub(r"[^A-Za-z0-9_-]+", "-", str(value)).strip("-")
        return slug or "results"

    @staticmethod
    def _html_text(value: object) -> str:
        if value is None or (isinstance(value, float) and pd.isna(value)):
            return ""
        return escape(str(value), quote=True)

    def generate(
        self,
        results: Dict[str, pd.DataFrame],
        output_file: str,
        gene_list: List[str] = None,
        ai_interpretation: Dict[str, Any] = None,
        pvalue_cutoff: Optional[float] = None,
        qvalue_cutoff: Optional[float] = None,
        gsea_results: pd.DataFrame = None,
        gsea_ranked_genes: List[str] = None,
        gsea_gene_weights: Dict[str, float] = None,
        gsea_gene_sets: Dict[str, set] = None,
        gsva_results: pd.DataFrame = None,
        gsva_groups: Dict[str, List[str]] = None,
        analysis_method: str = None,
        plot_types: List[str] = None,
        metadata: dict = None,
        ai_interpretation_error: Dict[str, Any] = None,
    ) -> str:
        self._ai_interpretation_error = ai_interpretation_error
        pvalue_cutoff = (
            getattr(self.config, 'pvalue_cutoff', 0.05)
            if pvalue_cutoff is None else pvalue_cutoff
        )
        qvalue_cutoff = (
            getattr(self.config, 'qvalue_cutoff', 0.05)
            if qvalue_cutoff is None else qvalue_cutoff
        )
        has_results = results and any(len(df) > 0 for df in results.values())

        if not has_results:
            html = self._generate_no_results_page(
                gene_list, metadata=metadata, ai_interpretation=ai_interpretation
            )
            with open(output_file, 'w', encoding='utf-8') as f:
                f.write(html)
            return output_file

        # ORA and GSEA reports apply the recorded significance thresholds.
        # Activity matrices have no per-row significance columns and remain unfiltered.
        sig_results = {}
        if analysis_method in ('ssgsea', 'gsva'):
            sig_results = {db_name: df.copy() for db_name, df in results.items() if len(df) > 0}
        else:
            for db_name, df in results.items():
                if len(df) == 0:
                    continue
                mask = pd.Series(True, index=df.index)
                p_column = next((column for column in ('pval', 'P_Value', 'p_value') if column in df.columns), None)
                q_column = next((column for column in ('padj', 'Adjusted_P_Value', 'FDR') if column in df.columns), None)
                if p_column:
                    mask &= df[p_column] <= pvalue_cutoff
                if q_column:
                    mask &= df[q_column] <= qvalue_cutoff
                sig_results[db_name] = df.loc[mask]

        # Preserve report structure even when no term passes the recorded filters.
        if not sig_results or all(len(df) == 0 for df in sig_results.values()):
            html = self._generate_no_results_page(
                gene_list, metadata=metadata, ai_interpretation=ai_interpretation
            )
            with open(output_file, 'w', encoding='utf-8') as f:
                f.write(html)
            return output_file

        if gsea_results is not None and len(gsea_results) > 0:
            db_col = 'Database' if 'Database' in gsea_results.columns else None
            if db_col:
                gsea_sig = {}
                for db_name in gsea_results[db_col].unique():
                    sub = gsea_results[gsea_results[db_col] == db_name].copy()
                    mask = pd.Series(True, index=sub.index)
                    q_column = next((column for column in ('padj', 'FDR', 'Adjusted_P_Value') if column in sub.columns), None)
                    p_column = next((column for column in ('pval', 'p_value', 'P_Value') if column in sub.columns), None)
                    if q_column:
                        mask &= sub[q_column] <= qvalue_cutoff
                    if p_column:
                        mask &= sub[p_column] <= pvalue_cutoff
                    filtered = sub.loc[mask]
                    if len(filtered) > 0:
                        gsea_sig[db_name] = filtered
                if gsea_sig:
                    sig_results = gsea_sig

        # ssGSEA and GSVA share the activity-matrix report layout.
        is_ssgsea_branch = analysis_method == 'ssgsea'
        is_gsva_branch = analysis_method == 'gsva' and gsva_results is not None
        if is_ssgsea_branch or is_gsva_branch:
            sections = {}
            saved_plots = self._generate_plot_section(results, analysis_method)
            has_saved_plots = bool(self._plot_file_groups())
            ai_section = self._generate_ai_content(ai_interpretation)
            for db_name, df in results.items():
                if 'NES' in df.columns:
                    df = df.reindex(df['NES'].abs().sort_values(ascending=False).index)
                sections[f'{db_name}_table'] = self._generate_tables(
                    {db_name: df}, analysis_method=analysis_method,
                    ai_interpretation=ai_interpretation,
                )

            # Generate fallback figures only when the analysis did not save plot files.
            gsva_section = ''
            if gsva_results is not None and not has_saved_plots:
                gsva_section = self._generate_gsva_plots_section(
                    gsva_results=gsva_results,
                    groups=gsva_groups,
                    plot_types=plot_types or ['heatmap', 'group_comparison', 'correlation'],
                    analysis_method=analysis_method,
                )
            if gsva_section:
                sections['gsva_visualization'] = gsva_section

            sections['summary'] = self._generate_summary(
                results, analysis_method=analysis_method
            )

            active_db_names = [db for db, df in results.items() if len(df) > 0]
            html = self._build_html(
                summary=sections.get('summary', ''),
                tables='\n'.join(v for k, v in sections.items() if k.endswith('_table')),
                plots=saved_plots,
                ai_section=ai_section,
                db_names=active_db_names,
                gsva_plots_html=sections.get('gsva_visualization', ''),
                analysis_method=analysis_method,
                metadata=metadata
            )

            with open(output_file, 'w', encoding='utf-8') as f:
                f.write(html)

            return output_file

        summary = self._generate_summary(sig_results, gene_list, analysis_method)
        tables = self._generate_tables(
            sig_results, analysis_method, ai_interpretation=ai_interpretation
        )
        plots = self._generate_plot_section(sig_results, analysis_method)
        has_saved_plots = bool(self._plot_file_groups())
        ai_section = self._generate_ai_content(ai_interpretation)

        # Generate fallback method-specific figures only when no saved figures exist.
        gsea_plots_html = ""
        gsva_plots_html = ""
        # Apply the exact style and semantic palette selection recorded by the run.
        plot_style = getattr(self.config, 'plot_style', 'nature') if self.config else 'nature'
        plot_palette = resolve_palette_selection(
            legacy_palette=getattr(self.config, 'plot_palette', None) if self.config else None,
            categorical_palette=(
                getattr(self.config, 'categorical_palette', None) if self.config else None
            ),
            sequential_palette=(
                getattr(self.config, 'sequential_palette', None) if self.config else None
            ),
            diverging_palette=(
                getattr(self.config, 'diverging_palette', None) if self.config else None
            ),
        )
        if not has_saved_plots and gsea_results is not None and len(gsea_results) > 0:
            gsea_plots_html = self._generate_gsea_plots_section(
                gsea_results, gsea_ranked_genes, gsea_gene_weights,
                gsea_gene_sets, plot_types, plot_style, plot_palette
            )
        if not has_saved_plots and gsva_results is not None and len(gsva_results) > 0:
            gsva_plots_html = self._generate_gsva_plots_section(
                gsva_results, gsva_groups, plot_types, plot_style, plot_palette
            )

        active_db_names = [db for db, df in sig_results.items() if len(df) > 0]
        html = self._build_html(
            summary=summary,
            tables=tables,
            plots=plots,
            ai_section=ai_section,
            db_names=active_db_names,
            gsea_plots_html=gsea_plots_html,
            gsva_plots_html=gsva_plots_html,
            analysis_method=analysis_method,
            metadata=metadata
        )

        with open(output_file, 'w', encoding='utf-8') as f:
            f.write(html)

        return output_file

    def _generate_no_results_page(
        self,
        gene_list: List[str] = None,
        metadata: dict = None,
        ai_interpretation: Dict[str, Any] = None,
    ) -> str:
        """Build a complete report when no terms pass the recorded filters."""
        _version = allenricher.__version__
        _db_ver = metadata.get("database_version", "") if metadata else ""
        _version_str = f"Version {_version}"
        if _db_ver:
            _version_str += f" | DB: {_db_ver}"
        methods_reference_html = render_methods_reference_html(metadata)
        ai_section = self._generate_ai_content(ai_interpretation)
        html = f'''<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>AllEnricher v{_version} - No Enrichment Results</title>
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
            font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Arial, sans-serif;
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
            font-family: Georgia, 'Times New Roman', serif;
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
        .methods-reference {{
            background: var(--bg-primary);
            border: 1px solid var(--border-color);
            border-radius: 4px;
            padding: 2rem 2.5rem;
            margin-top: 1.5rem;
        }}
        .methods-reference h2 {{
            font-family: Georgia, 'Times New Roman', serif;
            font-size: 1.25rem;
            margin-bottom: 1rem;
        }}
        .methods-reference h3 {{ font-size: 0.95rem; margin: 1rem 0 0.5rem; }}
        .methods-prose p {{ color: var(--text-secondary); margin-bottom: 0.75rem; }}
        .methods-references {{ margin-left: 1.25rem; color: var(--text-secondary); }}
        .methods-references li {{ margin-bottom: 0.5rem; }}
        .methods-reference a {{ color: var(--accent-color); }}
        .no-results-box h2 {{
            font-family: Georgia, 'Times New Roman', serif;
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
            border-top: 1px solid var(--border-color);
        }}
        .footer a {{
            color: var(--accent-color);
        }}
        .footer .citation {{
            max-width: 900px;
            margin: 0.4rem auto 0;
            line-height: 1.5;
        }}
    </style>
</head>
<body>
    <header class="header">
        <h1>AllEnricher Report</h1>
        <p class="meta">{_version_str} | {datetime.now().strftime("%Y-%m-%d")}</p>
    </header>
    <main class="main">
        <div class="no-results-box">
            <h2>No significant enrichment results found</h2>
            <p style="color: var(--text-secondary);">
                No terms met the significance criteria recorded for this analysis.
            </p>
            <div class="info-box">
                <h3>Possible causes</h3>
                <ul>
                    <li>The query contains too few genes that map to the selected database.</li>
                    <li>The configured P-value or adjusted P-value cutoff excludes all terms.</li>
                    <li>The selected background does not represent the tested gene universe.</li>
                    <li>The query identifier type does not match the database identifier type.</li>
                </ul>
            </div>
            <div class="info-box">
                <h3>Recommendations</h3>
                <ul>
                    <li>Review the input-to-database mapping reported in the analysis log.</li>
                    <li>Confirm that the background contains all genes eligible for selection.</li>
                    <li>Confirm that gene IDs match the identifier type used by the database.</li>
                    <li>Use a different database only when it addresses the analysis question.</li>
                </ul>
            </div>
        </div>
        {ai_section}
        {methods_reference_html}
    </main>
    <footer class="footer">
        <p>Generated by <a href="https://github.com/zhangducsu/AllEnricher-v2" target="_blank" rel="noopener noreferrer">AllEnricher v{_version}</a></p>
        <p class="citation"><strong>Citation:</strong> <cite>Zhang D, Hu Q, Liu X, et al. AllEnricher: a comprehensive gene set function enrichment tool for both model and non-model species. BMC Bioinformatics. 2020;21:106.</cite></p>
    </footer>
</body>
</html>'''
        return html

    def _generate_summary(self, results: Dict[str, pd.DataFrame], gene_list: List[str] = None, analysis_method: str = None) -> str:
        """Summarize the recorded analysis without interpreting its biological meaning."""
        total_terms = sum(len(df) for df in results.values())
        databases = list(results.keys())

        # ssGSEA and GSVA output a pathway-by-sample activity matrix.
        if analysis_method in ('ssgsea', 'gsva'):
            all_scores = []
            sample_cols_set = set()
            summary_stats = []
            for db_name, df in results.items():
                if len(df) > 0:
                    # Only numeric columns represent sample activity scores.
                    sample_cols = [c for c in df.columns if c not in [
                        'Term_ID', 'Term_Name', 'Gene_Count', 'Background_Count',
                        'Term_URL', 'NES', 'ES', 'P_Value', 'Adjusted_P_Value',
                        'FDR', 'Genes', 'Leading_Edge', 'Database',
                        'Rich_Factor', 'Gene_Ratio', 'Background_Ratio',
                        'Expected_Count'
                    ] and pd.api.types.is_numeric_dtype(df[c])]
                    sample_cols_set.update(sample_cols)
                    for col in sample_cols:
                        vals = pd.to_numeric(df[col], errors='coerce').dropna()
                        all_scores.extend(vals.tolist())
                    summary_stats.append({
                        "database": db_name,
                        "terms": len(df),
                    })

            # Ignore missing values when reporting the observed score range.
            numeric_scores = [x for x in all_scores if isinstance(x, (int, float))]
            score_min = min(numeric_scores) if numeric_scores else 0
            score_max = max(numeric_scores) if numeric_scores else 0
            n_samples = len(sample_cols_set)

            rows_html = "".join([
                f'<tr><td><a href="#{self._html_id(s["database"])}-table">'
                f'{self._html_text(database_display_name(s["database"]))}</a></td>'
                f'<td>{s["terms"]}</td></tr>'
                for s in summary_stats
            ])

            html = f'''
            <div class="section" id="summary">
                <h2>{'ssGSEA' if analysis_method == 'ssgsea' else 'GSVA'} Analysis Summary</h2>
                <div class="summary-grid">
                    <div class="stat-item">
                        <span class="stat-value">{total_terms}</span>
                        <span class="stat-label">Pathways</span>
                    </div>
                    <div class="stat-item">
                        <span class="stat-value">{n_samples}</span>
                        <span class="stat-label">Samples</span>
                    </div>
                    <div class="stat-item">
                        <span class="stat-value">{score_min:.2f} ~ {score_max:.2f}</span>
                        <span class="stat-label">Activity Score Range</span>
                    </div>
                </div>
                <table class="data-table">
                    <thead>
                        <tr>
                            <th>Database</th>
                            <th>Pathways</th>
                        </tr>
                    </thead>
                    <tbody>{rows_html}</tbody>
                </table>
            </div>'''
            return html

        summary_stats = []
        for db_name, df in results.items():
            if len(df) > 0:
                summary_stats.append({
                    "database": db_name,
                    "terms": len(df),
                    "min_pval": df['pval'].min() if 'pval' in df.columns
                                else (df['P_Value'].min() if 'P_Value' in df.columns
                                      else (df['p_value'].min() if 'p_value' in df.columns else 0)),
                    "min_adj_pval": df['padj'].min() if 'padj' in df.columns
                                    else (df['Adjusted_P_Value'].min() if 'Adjusted_P_Value' in df.columns
                                          else (df['FDR'].min() if 'FDR' in df.columns else 0))
                })

        rows_html = "".join([
            f'<tr><td><a href="#{self._html_id(s["database"])}-table">'
            f'{self._html_text(database_display_name(s["database"]))}</a></td>'
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
                    <span class="stat-value">{len(gene_list) if gene_list else " - "}</span>
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

    def _generate_tables(
        self,
        results: Dict[str, pd.DataFrame],
        analysis_method: str = None,
        ai_interpretation: Dict[str, Any] = None,
    ) -> str:
        """Build searchable result tables with stable term identifiers and descriptive names."""
        tables_html = []
        is_activity = analysis_method in {'ssgsea', 'gsva'}

        def format_number(value: object, pattern: str) -> str:
            try:
                if pd.isna(value):
                    return ""
                return format(float(value), pattern)
            except (TypeError, ValueError):
                return self._html_text(value)

        for db_name, df in results.items():
            if len(df) == 0:
                continue
            db_id = self._html_id(db_name)
            db_label = self._html_text(database_display_name(db_name))
            rows = []
            is_gsea = 'NES' in df.columns or ('ES' in df.columns and 'setSize' in df.columns)
            has_hierarchy = (
                'Hierarchy' in df.columns
                and df['Hierarchy'].fillna('').astype(str).str.strip().ne('').any()
            )

            if is_activity:
                metadata_columns = {
                    'Term_ID', 'Term_Name', 'Hierarchy', 'Database', 'Gene_Count',
                    'Background_Count', 'Term_URL', 'NES', 'ES', 'P_Value',
                    'Adjusted_P_Value', 'FDR', 'Genes', 'Leading_Edge',
                    'Rich_Factor', 'Gene_Ratio', 'Background_Ratio', 'Expected_Count',
                }
                sample_cols = [
                    column for column in df.columns
                    if column not in metadata_columns and pd.api.types.is_numeric_dtype(df[column])
                ]
                headers = ['Term ID', 'Term Name']
                if has_hierarchy:
                    headers.append('Hierarchy')
                headers.extend(map(str, sample_cols))
                for index, row in df.iterrows():
                    term_id = row.get('Term_ID', index)
                    term_name = row.get('Term_Name', term_id)
                    row_anchor = self._evidence_row_anchor(db_name, row, ai_interpretation)
                    cells = [self._html_text(term_id), self._html_text(term_name)]
                    if has_hierarchy:
                        cells.append(self._html_text(row.get('Hierarchy', '')))
                    cells.extend(format_number(row.get(column, ''), '.6g') for column in sample_cols)
                    rows.append(f'<tr{row_anchor}>' + ''.join(f'<td>{cell}</td>' for cell in cells) + '</tr>')
                method_label = 'ssGSEA' if analysis_method == 'ssgsea' else 'GSVA'
                count_label = 'pathways'
            elif is_gsea:
                headers = ['Term ID', 'Term Name']
                if has_hierarchy:
                    headers.append('Hierarchy')
                headers.extend(['pathway', 'pval', 'padj', 'log2err', 'ES', 'NES', 'size', 'leadingEdge'])
                for index, row in df.iterrows():
                    term_id = row.get('Term_ID', row.get('pathway', index))
                    term_name = row.get('Term_Name', term_id)
                    row_anchor = self._evidence_row_anchor(db_name, row, ai_interpretation)
                    leading_edge = str(row.get('leadingEdge', row.get('Lead_genes', '')))
                    cells = [self._html_text(term_id), self._html_text(term_name)]
                    if has_hierarchy:
                        cells.append(self._html_text(row.get('Hierarchy', '')))
                    cells.extend([
                        self._html_text(row.get('pathway', term_id)),
                        format_number(row.get('pval', row.get('P_Value', 1)), '.2e'),
                        format_number(row.get('padj', row.get('Adjusted_P_Value', 1)), '.2e'),
                        format_number(row.get('log2err', ''), '.4f'),
                        format_number(row.get('ES', row.get('es', 0)), '.4f'),
                        format_number(row.get('NES', row.get('nes', 0)), '.4f'),
                        self._html_text(row.get('size', row.get('setSize', 0))),
                        f'<span class="genes" data-full="{self._html_text(leading_edge)}">'
                        f'{self._html_text(self._truncate_genes(leading_edge))}</span>',
                    ])
                    rows.append(f'<tr{row_anchor}>' + ''.join(f'<td>{cell}</td>' for cell in cells) + '</tr>')
                method_label = 'GSEA'
                count_label = 'pathways'
            else:
                headers = ['Term ID', 'Term Name']
                if has_hierarchy:
                    headers.append('Hierarchy')
                headers.extend(['Gene Count', 'Rich Factor', 'P-value', 'Adj. P-value', 'Gene List'])
                for _, row in df.iterrows():
                    term_id = row.get('Term_ID', 'N/A')
                    term_name = row.get('Term_Name', term_id)
                    row_anchor = self._evidence_row_anchor(db_name, row, ai_interpretation)
                    term_url = self._html_text(row.get('Term_URL', ''))
                    term_id_text = self._html_text(term_id)
                    term_id_cell = (
                        f'<a href="{term_url}" target="_blank" rel="noopener">{term_id_text}</a>'
                        if term_url else term_id_text
                    )
                    genes = str(row.get('Genes', ''))
                    cells = [term_id_cell, self._html_text(term_name)]
                    if has_hierarchy:
                        cells.append(self._html_text(row.get('Hierarchy', '')))
                    cells.extend([
                        self._html_text(row.get('Gene_Count', 0)),
                        format_number(row.get('Rich_Factor', 0), '.4f'),
                        format_number(row.get('P_Value', 1), '.2e'),
                        format_number(row.get('Adjusted_P_Value', 1), '.2e'),
                        f'<span class="genes" data-full="{self._html_text(genes)}">'
                        f'{self._html_text(self._truncate_genes(genes))}</span>',
                    ])
                    rows.append(f'<tr{row_anchor}>' + ''.join(f'<td>{cell}</td>' for cell in cells) + '</tr>')
                method_label = 'ORA'
                count_label = 'terms'

            header_html = "".join(f'<th>{self._html_text(header)}</th>' for header in headers)

            table_html = f'''
            <div class="section" id="{db_id}-table">
                <h2>{db_label} {method_label} Results <span class="result-count">({len(df)} {count_label})</span></h2>
                <div class="table-wrapper">
                    <table class="data-table" id="table-{db_id}">
                        <thead>
                            <tr>{header_html}</tr>
                        </thead>
                        <tbody>{"".join(rows)}</tbody>
                    </table>
                </div>
                <div class="table-actions">
                    <button onclick="downloadTable('{db_id}')">Download TSV</button>
                    <button onclick="copyTable('{db_id}')">Copy</button>
                </div>
            </div>'''
            tables_html.append(table_html)

        return "\n".join(tables_html)

    def _evidence_row_anchor(
        self,
        database: str,
        row: pd.Series,
        ai_interpretation: Optional[Dict[str, Any]],
    ) -> str:
        """Return the stable HTML row id for code-selected evidence."""
        if not isinstance(ai_interpretation, dict):
            return ""
        evidence = ai_interpretation.get("evidence", {})
        if not isinstance(evidence, dict):
            return ""
        row_id = str(_first_value_for_report(row, ("Term_ID", "ID", "term_id", "pathway"), ""))
        row_name = str(_first_value_for_report(row, ("Term_Name", "Description", "term_name", "pathway"), ""))
        for evidence_id, record in evidence.items():
            if not isinstance(record, dict) or str(record.get("database")) != str(database):
                continue
            if row_id and str(record.get("term_id", "")) == row_id:
                return f' id="evidence-{self._html_id(evidence_id)}"'
            if row_name and str(record.get("term_name", "")) == row_name:
                return f' id="evidence-{self._html_id(evidence_id)}"'
        return ""

    @staticmethod
    def _truncate_genes(genes_str: str, max_len: int = 80) -> str:
        """Truncate long gene lists for table display while preserving the full result file."""
        if not genes_str or genes_str == 'nan':
            return ''
        if len(genes_str) > max_len:
            return genes_str[:max_len] + '...'
        return genes_str

    def _plot_file_groups(self) -> List[tuple[str, List[Path]]]:
        """Group alternate plot formats by figure stem."""
        supported = {'.png', '.jpg', '.jpeg', '.svg', '.pdf'}
        groups: Dict[str, List[Path]] = {}
        for path in self.output_dir.rglob('*'):
            if path.is_file() and path.suffix.lower() in supported:
                key = path.relative_to(self.output_dir).with_suffix('').as_posix()
                groups.setdefault(key, []).append(path)
        return sorted(groups.items())

    def _plot_identity(
        self,
        stem: str,
        results: Dict[str, pd.DataFrame],
        analysis_method: Optional[str],
    ) -> tuple[str, str]:
        term_names: Dict[str, tuple[str, str, str]] = {}
        for database, frame in results.items():
            id_column = next(
                (column for column in ('Term_ID', 'pathway', 'term_id') if column in frame.columns),
                None,
            )
            if id_column is None:
                continue
            for _, row in frame.iterrows():
                term_id = str(row.get(id_column, '')).strip()
                term_name = str(row.get('Term_Name', term_id)).strip() or term_id
                if term_id:
                    identity = (database, term_name, term_id)
                    term_names[term_id] = identity
                    term_names[re.sub(r'[^A-Za-z0-9_.-]+', '_', term_id)] = identity

        for safe_term_id, (database, term_name, term_id) in term_names.items():
            if stem == f'{safe_term_id}_enrichment':
                return database, f'{term_name} ({term_id}) - GSEA Enrichment Plot'

        database = next(
            (name for name in results if stem.casefold().startswith(f'{name}_'.casefold())),
            next(iter(results), 'Analysis'),
        )
        base = stem[len(database) + 1:] if stem.casefold().startswith(f'{database}_'.casefold()) else stem
        method_label = {
            'hypergeometric': 'ORA', 'gsea': 'GSEA',
            'ssgsea': 'ssGSEA', 'gsva': 'GSVA',
        }.get(analysis_method, 'Enrichment')
        labels = (
            ('activity_heatmap', f'{method_label} Pathway Activity Heatmap'),
            ('sample_correlation', f'{method_label} Sample Correlation'),
            ('group_comparison', f'{method_label} Group Comparison'),
            ('enrichment2_up', 'GSEA Multi-pathway Enrichment - Up'),
            ('enrichment2_down', 'GSEA Multi-pathway Enrichment - Down'),
            ('lollipop', f'{method_label} Lollipop Plot'),
            ('emapplot', f'{method_label} Pathway Network'),
            ('ridgeplot', 'GSEA Ridge Plot'),
            ('barplot', f'{method_label} Bar Plot'),
        )
        title = next((label for token, label in labels if token in base.casefold()), base.replace('_', ' '))
        return database, title

    def _generate_plot_section(
        self,
        results: Dict[str, pd.DataFrame],
        analysis_method: Optional[str] = None,
    ) -> str:
        """Embed every generated figure type with one preview per format group."""
        groups = self._plot_file_groups()
        if not groups:
            return (
                '<div class="section" id="plots"><h2>Figures</h2>'
                '<div class="plot-placeholder"><p>No figures were generated.</p></div></div>'
            )

        by_database: Dict[str, List[str]] = {}
        priority = {'.png': 0, '.jpg': 1, '.jpeg': 1, '.svg': 2, '.pdf': 3}
        mime_types = {
            '.png': 'image/png', '.jpg': 'image/jpeg', '.jpeg': 'image/jpeg',
            '.svg': 'image/svg+xml', '.pdf': 'application/pdf',
        }
        for _, files in groups:
            files = sorted(files, key=lambda path: priority[path.suffix.lower()])
            preview = files[0]
            database, title = self._plot_identity(preview.stem, results, analysis_method)
            encoded = self._encode_image_to_base64(preview)
            suffix = preview.suffix.lower()
            title_html = self._html_text(title)
            if suffix == '.pdf':
                media = (
                    f'<object data="data:application/pdf;base64,{encoded}" '
                    f'type="application/pdf" class="plot-pdf" title="{title_html}">'
                    f'<p>PDF preview unavailable.</p></object>'
                )
            else:
                media = (
                    f'<img src="data:{mime_types[suffix]};base64,{encoded}" '
                    f'alt="{title_html}" class="plot-img">'
                )
            format_links = ''.join(
                f'<a href="{self._html_text(path.relative_to(self.output_dir).as_posix())}" '
                f'target="_blank" rel="noopener">{path.suffix[1:].upper()}</a>'
                for path in files
            )
            relative = self._html_text(preview.relative_to(self.output_dir).as_posix())
            card = (
                '<div class="plot-container">'
                f'<p class="plot-title">{title_html}</p>{media}'
                f'<p class="plot-caption">{relative}</p>'
                f'<div class="plot-formats">{format_links}</div></div>'
            )
            by_database.setdefault(database, []).append(card)

        parts = [
            '<div class="section" id="plots">',
            f'<h2>Figures <span class="result-count">({len(groups)} plot types)</span></h2>',
        ]
        for database, cards in by_database.items():
            parts.extend([
                f'<div class="plot-group"><h3>{self._html_text(database_display_name(database))} Figures</h3>',
                '<div class="plot-grid">', *cards, '</div></div>',
            ])
        parts.append('</div>')
        return '\n'.join(parts)

    def _encode_image_to_base64(self, image_path: Path) -> str:
        """Encode an image for embedding in a self-contained HTML report."""
        try:
            with open(image_path, 'rb') as f:
                return base64.b64encode(f.read()).decode('utf-8')
        except Exception as e:
            logger.warning("Failed to embed image %s: %s", image_path, e)
            return ""

    def _generate_gsea_plots_section(
        self,
        gsea_results: pd.DataFrame,
        ranked_genes: List[str] = None,
        gene_weights: Dict[str, float] = None,
        gene_sets: Dict[str, set] = None,
        plot_types: List[str] = None,
        plot_style: str = 'nature',
        plot_palette: PaletteLike = None,
    ) -> str:
        """Build the GSEA figure section from recorded analysis inputs."""
        if plot_types is None:
            plot_types = ["enrichment"]

        # Temporary figures are embedded before the directory is removed.
        temp_dir = Path(tempfile.mkdtemp(prefix="gsea_plots_"))
        html_parts = ['<div class="section" id="gsea-plots">',
                       '<h2>GSEA Visualization</h2>',
                       '<div class="plot-grid">']

        has_any_plot = False

        try:
            # Render at most five leading pathways when running-ES inputs are available.
            if "enrichment" in plot_types and ranked_genes and gene_weights and gene_sets:
                top_pathways = gsea_results.head(5)
                for _, row in top_pathways.iterrows():
                    pw_name = row.get("pathway", "")
                    if pw_name not in gene_sets:
                        continue
                    try:
                        from allenricher.visualization.gsea_plots import plot_gsea_enrichment
                        output_file = str(temp_dir / f"enrichment_{pw_name[:50]}.png")
                        fig = plot_gsea_enrichment(
                            ranked_genes=ranked_genes,
                            gene_weights=gene_weights,
                            gene_set=gene_sets[pw_name],
                            es=row.get("ES", row.get("es", 0.0)),
                            nes=row.get("NES", row.get("nes", 0.0)),
                            pvalue=row.get("pval", row.get("pvalue", 1.0)),
                            padj=row.get("padj"),
                            title=pw_name,
                            output_file=output_file,
                            dpi=150,
                            style=plot_style,
                            palette=plot_palette,
                        )
                        import matplotlib.pyplot as plt
                        plt.close(fig)

                        img_data = self._encode_image_to_base64(Path(output_file))
                        if img_data:
                            has_any_plot = True
                            html_parts.append(
                                f'<div class="gsea-plot-container">'
                                f'<p class="plot-title">{pw_name}</p>'
                                f'<img src="data:image/png;base64,{img_data}" '
                                f'alt="GSEA Enrichment: {pw_name}" class="plot-img">'
                                f'<p class="plot-caption">NES = {row.get("NES", row.get("nes", 0)):.2f}, '
                                f'P-value = {row.get("pval", row.get("pvalue", 1)):.2e}</p>'
                                f'</div>'
                            )
                    except Exception as e:
                        logger.warning(f"Failed to generate GSEA enrichment plot ({pw_name}): {e}")

        finally:
            # Embedded images no longer depend on their temporary files.
            try:
                shutil.rmtree(temp_dir, ignore_errors=True)
            except Exception:
                pass

        if not has_any_plot:
            html_parts.append(
                '<div class="plot-placeholder">'
                '<p>No GSEA figure could be generated from the recorded plotting inputs.</p>'
                '</div>'
            )

        html_parts.append('</div></div>')
        return "\n".join(html_parts)

    def _generate_gsva_plots_section(
        self,
        gsva_results: pd.DataFrame,
        groups: Dict[str, List[str]] = None,
        plot_types: List[str] = None,
        plot_style: str = 'nature',
        plot_palette: PaletteLike = None,
        analysis_method: str = 'gsva',
    ) -> str:
        """Build the ssGSEA or GSVA activity figure section."""
        if plot_types is None:
            plot_types = ["heatmap", "group_comparison", "correlation"]
        method_label = 'ssGSEA' if analysis_method == 'ssgsea' else 'GSVA'

        # Temporary figures are embedded before the directory is removed.
        temp_dir = Path(tempfile.mkdtemp(prefix="gsva_plots_"))
        html_parts = ['<div class="section" id="gsva-plots">',
                       f'<h2>{method_label} Visualization</h2>',
                       '<div class="plot-grid">']

        has_any_plot = False

        try:
            # Pathway-by-sample activity heatmap.
            if "heatmap" in plot_types:
                try:
                    from allenricher.visualization.gsva_plots import plot_pathway_heatmap
                    output_file = str(temp_dir / "pathway_heatmap.png")
                    fig = plot_pathway_heatmap(
                        scores_df=gsva_results,
                        title=f'{method_label} Pathway Activity',
                        output_file=output_file,
                        dpi=150,
                        style=plot_style,
                        palette=plot_palette,
                    )
                    import matplotlib.pyplot as plt
                    plt.close(fig)

                    img_data = self._encode_image_to_base64(Path(output_file))
                    if img_data:
                        has_any_plot = True
                        html_parts.append(
                            f'<div class="gsva-plot-container">'
                            f'<p class="plot-title">{method_label} Pathway Activity</p>'
                            f'<img src="data:image/png;base64,{img_data}" '
                            f'alt="Pathway Heatmap" class="plot-img">'
                            f'<p class="plot-caption">Sample-pathway activity scores with hierarchical clustering</p>'
                            f'</div>'
                        )
                except Exception as e:
                    logger.warning("Failed to generate the pathway activity heatmap: %s", e)

            # Group comparison requires at least two recorded sample groups.
            if "group_comparison" in plot_types and groups and len(groups) >= 2:
                try:
                    from allenricher.visualization.gsva_plots import plot_group_comparison
                    output_file = str(temp_dir / "group_comparison.png")
                    fig = plot_group_comparison(
                        scores_df=gsva_results,
                        groups=groups,
                        output_file=output_file,
                        dpi=150,
                        style=plot_style,
                        palette=plot_palette,
                    )
                    import matplotlib.pyplot as plt
                    plt.close(fig)

                    img_data = self._encode_image_to_base64(Path(output_file))
                    if img_data:
                        has_any_plot = True
                        html_parts.append(
                            f'<div class="gsva-plot-container">'
                            f'<p class="plot-title">Group Comparison</p>'
                            f'<img src="data:image/png;base64,{img_data}" '
                            f'alt="Group Comparison" class="plot-img">'
                            f'<p class="plot-caption">Pathway activity comparison between sample groups</p>'
                            f'</div>'
                        )
                except Exception as e:
                    logger.warning("Failed to generate the group comparison figure: %s", e)

            # Sample correlation based on pathway activity profiles.
            if "correlation" in plot_types:
                try:
                    from allenricher.visualization.gsva_plots import plot_sample_correlation
                    output_file = str(temp_dir / "sample_correlation.png")
                    fig = plot_sample_correlation(
                        scores_df=gsva_results,
                        output_file=output_file,
                        dpi=150,
                        style=plot_style,
                        palette=plot_palette,
                    )
                    import matplotlib.pyplot as plt
                    plt.close(fig)

                    img_data = self._encode_image_to_base64(Path(output_file))
                    if img_data:
                        has_any_plot = True
                        html_parts.append(
                            f'<div class="gsva-plot-container">'
                            f'<p class="plot-title">Sample Correlation</p>'
                            f'<img src="data:image/png;base64,{img_data}" '
                            f'alt="Sample Correlation" class="plot-img">'
                            f'<p class="plot-caption">Pearson correlation between samples based on pathway activity</p>'
                            f'</div>'
                        )
                except Exception as e:
                    logger.warning(f"Failed to generate sample-correlation plot: {e}")

        finally:
            # Embedded images no longer depend on their temporary files.
            try:
                shutil.rmtree(temp_dir, ignore_errors=True)
            except Exception:
                pass

        if not has_any_plot:
            html_parts.append(
                '<div class="plot-placeholder">'
                f'<p>No {method_label} figure could be generated from the recorded plotting inputs.</p>'
                '</div>'
            )

        html_parts.append('</div></div>')
        return "\n".join(html_parts)

    def _generate_ai_content(self, ai_interpretation: Dict[str, Any] = None) -> str:
        if ai_interpretation:
            return self._generate_ai_section(ai_interpretation)
        error = getattr(self, "_ai_interpretation_error", None)
        if not error:
            return ""
        code = self._html_text(error.get("error_code", "AI_INTERPRETATION_FAILED"))
        backend = self._html_text(error.get("backend", "unknown"))
        mode = self._html_text(error.get("mode", "summary"))
        message = self._html_text(error.get("message", "Unknown AI interpretation error"))
        return f'''
        <div class="section" id="ai-interpretation">
            <h2>AI Interpretation</h2>
            <div class="ai-error">
                <strong>AI interpretation unavailable.</strong>
                The enrichment analysis completed successfully; only the optional AI layer failed.
                <p><strong>Error code:</strong> {code}<br>
                <strong>Backend:</strong> {backend}<br>
                <strong>Profile:</strong> {mode}<br>
                <strong>Details:</strong> {message}</p>
            </div>
        </div>'''

    def _generate_ai_section(self, ai_interpretation: Dict[str, Any]) -> str:
        """Build the optional AI interpretation section."""
        if not ai_interpretation:
            return ""

        if ai_interpretation.get("schema_version") == 1:
            def evidence_links(evidence_ids: List[str]) -> str:
                return " ".join(
                    f'<a class="ai-evidence-link" href="#evidence-{self._html_id(evidence_id)}">'
                    f'{self._html_text(evidence_id)}</a>'
                    for evidence_id in evidence_ids
                )

            section_labels = {
                "core_themes": "Core themes",
                "key_evidence": "Key evidence",
                "limitations": "Contradictions and limitations",
                "computational_checks": "Computational checks",
            }
            database_blocks = []
            for database, content in ai_interpretation.get("databases", {}).items():
                groups = []
                for section, label in section_labels.items():
                    items = content.get(section, []) if isinstance(content, dict) else []
                    entries = []
                    for item in items:
                        text = self._html_text(item.get("text", ""))
                        ids = item.get("evidence_ids", [])
                        links = evidence_links(ids)
                        evidence_markup = (
                            f' <span class="ai-evidence">[{links}]</span>' if links else ""
                        )
                        entries.append(
                            f'<li><span>{text}</span>'
                            f'{evidence_markup}</li>'
                        )
                    if not entries:
                        entries.append('<li class="ai-empty">None recorded.</li>')
                    groups.append(
                        f'<div class="ai-subsection"><h4>{label}</h4><ul>{"".join(entries)}</ul></div>'
                    )
                database_blocks.append(
                    f'<div class="ai-block"><h3>{self._html_text(database_display_name(database))}</h3>'
                    f'{"".join(groups)}</div>'
                )
            mode = self._html_text(ai_interpretation.get("profile", "summary"))
            backend_labels = {
                "openai": "OpenAI",
                "claude": "Claude",
                "deepseek": "DeepSeek",
                "glm": "GLM",
                "minimax": "MiniMax",
                "ollama": "Ollama",
                "mock": "Mock",
            }
            backend_value = str(ai_interpretation.get("backend", "unknown"))
            backend = self._html_text(
                backend_labels.get(backend_value.lower(), backend_value)
            )
            return f'''
        <div class="section" id="ai-interpretation">
            <h2>AI Interpretation <span class="ai-note">(Model: {backend} | Profile: {mode}; requires expert review)</span></h2>
            {"".join(database_blocks)}
            <div class="ai-disclaimer">
                <strong>Important:</strong> AI text is restricted to the linked evidence and must be checked against the source tables before use.
            </div>
        </div>'''

        interpretations = []
        for db_name, interpretation in ai_interpretation.items():
            if not isinstance(interpretation, str):
                continue
            # Escape model output before supporting line breaks and bold Markdown.
            html_text = escape(str(interpretation), quote=False).replace('\n', '<br>')
            html_text = re.sub(r'\*\*(.+?)\*\*', r'<strong>\1</strong>', html_text)
            interpretations.append(
                f'<div class="ai-block"><h3>{self._html_text(database_display_name(db_name))}</h3>'
                f'<p>{html_text}</p></div>'
            )

        html = f'''
        <div class="section" id="ai-interpretation">
            <h2>AI Interpretation <span class="ai-note">(Requires expert review)</span></h2>
            {"".join(interpretations)}
            <div class="ai-disclaimer">
                <strong>Important:</strong> This text was generated from the displayed enrichment tables. It may contain errors and must be checked against the source results, the study design, and verified literature before use.
            </div>
        </div>'''
        return html

    def _build_html(self, summary: str, tables: str, plots: str, ai_section: str,
                    db_names: List[str] = None,
                    gsea_plots_html: str = "",
                    gsva_plots_html: str = "",
                    analysis_method: str = None,
                    metadata: dict = None) -> str:
        """Assemble the complete self-contained HTML document."""

        _version = allenricher.__version__
        _db_ver = metadata.get("database_version", "") if metadata else ""
        _version_str = f"Version {_version}"
        if _db_ver:
            _version_str += f" | DB: {_db_ver}"

        # Match the report title to the method that actually ran.
        _title_map = {
            'hypergeometric': 'ORA Enrichment Analysis Report',
            'ssgsea': 'ssGSEA Pathway Activity Report',
            'gsva': 'GSVA Pathway Activity Report',
            'gsea': 'GSEA Enrichment Analysis Report',
        }
        _report_title = _title_map.get(analysis_method, 'Enrichment Analysis Report')
        methods_reference_html = render_methods_reference_html(metadata)

        nav_items = '<li><a href="#summary">Summary</a></li>'
        nav_items += '<li><a href="#methods-reference">Materials and Methods</a></li>'
        if plots:
            nav_items += '<li><a href="#plots">Figures</a></li>'
        # Add method-specific figure links only when those sections exist.
        if gsea_plots_html:
            nav_items += '<li><a href="#gsea-plots">GSEA Figures</a></li>'
        if gsva_plots_html:
            activity_label = 'ssGSEA' if analysis_method == 'ssgsea' else 'GSVA'
            nav_items += f'<li><a href="#gsva-plots">{activity_label} Figures</a></li>'
        if db_names:
            nav_items += ''.join(
                f'<li><a href="#{self._html_id(db)}-table">'
                f'{self._html_text(database_display_name(db))} Results</a></li>'
                for db in db_names
            )
        if ai_section:
            nav_items += '<li><a href="#ai-interpretation">AI Interpretation</a></li>'

        # Keep method-specific figures together before the result tables.
        gsea_gsva_sections = ""
        if gsea_plots_html:
            gsea_gsva_sections += gsea_plots_html
        if gsva_plots_html:
            gsea_gsva_sections += gsva_plots_html

        html = f'''<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>AllEnricher v{_version} - {_report_title}</title>
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
            font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Arial, sans-serif;
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
            display: block;
        }}

        .header h1 {{
            font-family: Georgia, 'Times New Roman', serif;
            font-size: 1.25rem;
            font-weight: 600;
            color: var(--text-primary);
        }}

        .header .meta {{
            display: block;
            font-size: 0.8rem;
            color: var(--text-muted);
            margin-top: 0.35rem;
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
            font-family: Georgia, 'Times New Roman', serif;
            font-size: 1.1rem;
            font-weight: 600;
            color: var(--text-primary);
            margin-bottom: 1.25rem;
            padding-bottom: 0.5rem;
            border-bottom: 1px solid var(--border-color);
        }}

        .section h2 .result-count,
        .section h2 .ai-note {{
            font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Arial, sans-serif;
            font-weight: 400;
            font-size: 0.875rem;
            color: var(--text-muted);
        }}

        .methods-reference h3 {{
            font-size: 0.95rem;
            margin: 1rem 0 0.5rem;
            color: var(--text-secondary);
        }}

        .methods-prose p {{
            margin-bottom: 0.75rem;
            color: var(--text-secondary);
        }}

        .methods-references {{
            margin-left: 1.25rem;
            color: var(--text-secondary);
        }}

        .methods-references li {{ margin-bottom: 0.55rem; }}

        .methods-reference a {{ color: var(--accent-color); }}

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
            max-width: 350px;
            overflow: hidden;
            text-overflow: ellipsis;
            white-space: nowrap;
        }}

        .data-table .genes:hover {{
            white-space: normal;
            word-break: break-all;
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

        /* Prot Images - PNG Embedded Styles
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
            border: 1px solid var(--border-color);
            background: white;
        }}

        .plot-pdf {{
            width: 100%;
            min-height: 700px;
            border: 1px solid var(--border-color);
            background: white;
        }}

        .plot-formats {{
            display: flex;
            justify-content: center;
            gap: 0.75rem;
            margin-top: 0.35rem;
            font-size: 0.8rem;
        }}

        .plot-formats a {{
            color: var(--accent-color);
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

        .ai-subsection {{ margin-top: 0.75rem; }}
        .ai-subsection h4 {{
            font-size: 0.82rem;
            color: var(--text-secondary);
            margin-bottom: 0.25rem;
        }}
        .ai-subsection ul {{ margin: 0 0 0 1.1rem; color: var(--text-secondary); }}
        .ai-subsection li {{ margin: 0.25rem 0; line-height: 1.55; }}
        .ai-evidence-link {{ color: var(--accent-color); white-space: nowrap; }}
        .ai-empty {{ color: var(--text-muted); font-style: italic; }}

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

        .footer .citation {{
            max-width: 900px;
            margin: 0.4rem auto 0;
            line-height: 1.5;
        }}

        /* Table Controls */
        .table-controls {{
            display: flex;
            justify-content: space-between;
            align-items: center;
            margin-bottom: 1rem;
            flex-wrap: wrap;
            gap: 0.75rem;
        }}

        .table-search {{
            display: flex;
            align-items: center;
            gap: 0.5rem;
        }}

        .table-search input {{
            padding: 0.4rem 0.75rem;
            border: 1px solid var(--border-color);
            border-radius: 4px;
            font-size: 0.85rem;
            min-width: 200px;
        }}

        .table-search input:focus {{
            outline: none;
            border-color: var(--accent-color);
        }}

        .table-length {{
            display: flex;
            align-items: center;
            gap: 0.5rem;
            font-size: 0.85rem;
            color: var(--text-secondary);
        }}

        .table-length select {{
            padding: 0.4rem 0.5rem;
            border: 1px solid var(--border-color);
            border-radius: 4px;
            font-size: 0.85rem;
            background: var(--bg-primary);
        }}

        .table-info {{
            color: var(--text-muted);
            font-size: 0.8rem;
            margin-top: 0.75rem;
        }}

        /* Pagination */
        .pagination {{
            display: flex;
            justify-content: center;
            align-items: center;
            gap: 0.25rem;
            margin-top: 1rem;
            flex-wrap: wrap;
        }}

        .pagination button {{
            padding: 0.4rem 0.75rem;
            border: 1px solid var(--border-color);
            background: var(--bg-primary);
            color: var(--text-secondary);
            border-radius: 4px;
            cursor: pointer;
            font-size: 0.85rem;
            transition: all 0.2s;
        }}

        .pagination button:hover:not(:disabled) {{
            background: var(--bg-secondary);
            border-color: var(--accent-color);
        }}

        .pagination button:disabled {{
            opacity: 0.5;
            cursor: not-allowed;
        }}

        .pagination button.active {{
            background: var(--accent-color);
            color: white;
            border-color: var(--accent-color);
        }}

        /* Sortable Headers */
        th.sortable {{
            cursor: pointer;
            user-select: none;
            position: relative;
            padding-right: 1.5rem;
        }}

        th.sortable:hover {{
            background: var(--bg-secondary);
        }}

        th.sortable::after {{
            content: '↕';
            position: absolute;
            right: 0.5rem;
            opacity: 0.4;
            font-size: 0.75rem;
        }}

        th.sortable.sort-asc::after {{
            content: '↑';
            opacity: 1;
            color: var(--accent-color);
        }}

        th.sortable.sort-desc::after {{
            content: '↓';
            opacity: 1;
            color: var(--accent-color);
        }}

        th.sortable.sort-disabled {{
            cursor: default;
        }}

        th.sortable.sort-disabled::after {{
            content: none;
        }}

        /* Plot Placeholder */
        .plot-placeholder {{
            padding: 2rem;
            text-align: center;
            color: var(--text-muted);
            background: var(--bg-secondary);
            border: 1px dashed var(--border-color);
            border-radius: 4px;
        }}

        .plot-placeholder code {{
            background: var(--bg-primary);
            padding: 0.15rem 0.4rem;
            border-radius: 3px;
            font-size: 0.8rem;
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

            .plot-grid {{
                grid-template-columns: 1fr !important;
            }}
        }}

        /* GSEA/GSVA Plot Styles */
        .gsea-plot-container,
        .gsva-plot-container {{
            background: var(--bg-primary);
            border: 1px solid var(--border-color);
            border-radius: 4px;
            padding: 1rem;
            text-align: center;
        }}

        .plot-grid {{
            display: grid;
            grid-template-columns: minmax(0, 1fr);
            gap: 1.5rem;
        }}

        .plot-title {{
            font-family: Georgia, 'Times New Roman', serif;
            font-size: 0.95rem;
            font-weight: 600;
            color: var(--text-primary);
            margin-bottom: 0.75rem;
        }}
    </style>
</head>
<body>
    <header class="header">
        <div class="header-content">
            <h1>AllEnricher Report - {_report_title}</h1>
            <span class="meta">{_version_str} | {datetime.now().strftime("%Y-%m-%d")}</span>
        </div>
    </header>

    <nav class="nav">
        <ul>{nav_items}</ul>
    </nav>

    <main class="main">
        {summary}
        {methods_reference_html}
        {plots}
        {gsea_gsva_sections}
        {tables}
        {ai_section}
    </main>

    <footer class="footer">
        <p>Generated by <a href="https://github.com/zhangducsu/AllEnricher-v2" target="_blank" rel="noopener noreferrer">AllEnricher v{_version}</a></p>
        <p class="citation"><strong>Citation:</strong> <cite>Zhang D, Hu Q, Liu X, et al. AllEnricher: a comprehensive gene set function enrichment tool for both model and non-model species. BMC Bioinformatics. 2020;21:106.</cite></p>
    </footer>

    {self._TABLE_JS}
</body>
</html>'''
        return html
