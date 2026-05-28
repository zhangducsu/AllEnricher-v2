#!/usr/bin/env python3
"""
AllEnricher v2.0 命令行接口模块 (Command Line Interface)

本模块是 AllEnricher 工具的命令行入口，提供以下核心功能：
  - analyze:  运行基因集功能富集分析（主工作流）
  - download: 下载指定的富集分析数据库（如 GO、KEGG 等）
  - build:    为指定物种构建本地数据库
  - serve:    启动 RESTful API 服务器，提供在线分析服务
  - list:     列出支持的物种或数据库资源
  - config:   生成默认配置文件（YAML/JSON 格式）

使用示例：
    allenricher analyze -i genes.txt -s hsa -d GO,KEGG
    allenricher download -d GO,KEGG -s hsa
    allenricher serve --port 8000
"""

import argparse
import logging
import sys
from pathlib import Path
from typing import Dict, List, Optional, Set

from allenricher import __version__
from allenricher.core.config import Config
from allenricher.core.enrichment import EnrichmentAnalyzer
from allenricher.database.manager import DatabaseManager
from allenricher.visualization.plotter import Plotter
from allenricher.report.generator import ReportGenerator
from allenricher.ai.interpreter import create_interpreter, create_interpreter_from_config

# 配置日志输出格式
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# GSEA/GSVA/ssGSEA 可视化辅助函数
# ---------------------------------------------------------------------------

# 各方法支持的图表类型白名单
_METHOD_PLOT_TYPES = {
    'gsea': {'enrichment', 'nes_barplot', 'dotplot'},
    'ssgsea': {'heatmap', 'group_comparison', 'dotplot', 'correlation'},
    'gsva': {'heatmap', 'group_comparison', 'dotplot', 'correlation'},
}

# 通用图表类型（不依赖于特定分析方法）
_COMMON_PLOT_TYPES = {'network', 'upset', 'volcano'}


def _parse_gmt_file(gmt_file: str) -> Dict[str, Set[str]]:
    """读取GMT格式基因集文件

    GMT文件格式：每行以TAB分隔，第1列为基因集名称，第2列为描述（可选），
    第3列起为基因名称列表。支持 .gmt 和 .gmt.gz 两种格式。

    Args:
        gmt_file: GMT文件路径

    Returns:
        Dict[str, Set[str]]: 基因集名称到基因集合的映射

    Raises:
        FileNotFoundError: 文件不存在
        ValueError: 文件格式不正确
    """
    import gzip

    gmt_path = Path(gmt_file)
    if not gmt_path.exists():
        raise FileNotFoundError(f"GMT文件不存在: {gmt_file}")

    # 根据扩展名选择打开方式
    if gmt_path.suffix.lower() == '.gz':
        opener = gzip.open
        mode = 'rt'
    else:
        opener = open
        mode = 'r'

    gene_sets: Dict[str, Set[str]] = {}
    with opener(gmt_path, mode, encoding='utf-8') as f:
        for line_num, line in enumerate(f, 1):
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            parts = line.split('\t')
            if len(parts) < 3:
                logger.warning(f"GMT文件第{line_num}行格式不正确（至少需要3列），已跳过")
                continue
            set_name = parts[0]
            # 第3列起为基因名称
            genes = {g.strip() for g in parts[2:] if g.strip()}
            if genes:
                gene_sets[set_name] = genes

    logger.info(f"GMT文件加载完成: {len(gene_sets)} 个基因集, 来源: {gmt_file}")
    return gene_sets


def _parse_groups(groups_str: str) -> Dict[str, List[str]]:
    """解析分组字符串

    将 'Group1:sample1,sample2;Group2:sample3,sample4' 格式的字符串
    解析为 {组名: [样本名列表]} 字典。

    Args:
        groups_str: 分组定义字符串

    Returns:
        Dict[str, List[str]]: 分组名称到样本名列表的映射

    Raises:
        ValueError: 格式不正确
    """
    if not groups_str:
        return {}

    groups: Dict[str, List[str]] = {}
    for group_def in groups_str.split(';'):
        group_def = group_def.strip()
        if not group_def:
            continue
        if ':' not in group_def:
            raise ValueError(
                f"分组定义格式错误: '{group_def}'，"
                f"期望格式为 'GroupName:sample1,sample2'"
            )
        group_name, samples_str = group_def.split(':', 1)
        group_name = group_name.strip()
        samples = [s.strip() for s in samples_str.split(',') if s.strip()]
        if not group_name:
            raise ValueError("分组名称不能为空")
        if not samples:
            raise ValueError(f"分组 '{group_name}' 中没有样本")
        groups[group_name] = samples

    logger.info(f"样本分组解析完成: {len(groups)} 个分组 - {list(groups.keys())}")
    return groups


def _generate_plots(
    method: str,
    results: dict,
    ranked_genes: Optional[list],
    gene_weights: Optional[dict],
    gene_sets: Optional[Dict[str, Set[str]]],
    expr_matrix,
    groups: Optional[Dict[str, List[str]]],
    plot_types: List[str],
    output_dir: str,
    plot_format: str = 'png',
    plot_dpi: int = 300,
    plot_style: str = 'nature',
    plot_palette: Optional[str] = None,
) -> List[str]:
    """根据方法类型生成可视化图表

    根据分析方法和用户指定的图表类型，调用对应的可视化函数生成图表。

    Args:
        method: 分析方法 (gsea/ssgsea/gsva)
        results: 分析结果字典 {db_name: DataFrame}
        ranked_genes: 排序基因列表（GSEA需要）
        gene_weights: 基因权重字典（GSEA需要）
        gene_sets: 基因集字典 {set_name: gene_set}（GSEA富集曲线需要）
        expr_matrix: 表达矩阵 DataFrame（ssGSEA/GSVA需要）
        groups: 样本分组字典（组间比较图需要）
        plot_types: 要生成的图表类型列表
        output_dir: 输出目录
        plot_format: 图表格式 (png/pdf/svg)
        plot_dpi: 图表DPI
        plot_style: 图表风格主题 (nature, science, colorblind, presentation, omicshare)
        plot_palette: 自定义配色方案名称（可选）

    Returns:
        List[str]: 生成的图表文件路径列表
    """
    generated_files: List[str] = []

    if not plot_types or not results:
        return generated_files

    # 校验图表类型是否支持
    supported = _METHOD_PLOT_TYPES.get(method, set())
    for pt in plot_types:
        if pt not in supported:
            logger.warning(f"方法 '{method}' 不支持图表类型 '{pt}'，跳过。"
                           f"支持的类型: {sorted(supported)}")

    valid_types = [pt for pt in plot_types if pt in supported]
    if not valid_types:
        return generated_files

    # 确保输出目录存在
    plot_dir = Path(output_dir) / "gsea_plots"
    plot_dir.mkdir(parents=True, exist_ok=True)

    # ---- GSEA 图表 ----
    if method == 'gsea':
        from allenricher.visualization.gsea_plots import (
            plot_gsea_enrichment,
            plot_gsea_nes_barplot,
            plot_gsea_dotplot,
        )

        for db_name, df in results.items():
            if df is None or len(df) == 0:
                continue

            # NES 条形图
            if 'nes_barplot' in valid_types and 'nes' in df.columns:
                out_file = str(plot_dir / f"{db_name}_nes_barplot.{plot_format}")
                try:
                    plot_gsea_nes_barplot(df, output_file=out_file, dpi=plot_dpi,
                                          style=plot_style, palette=plot_palette)
                    generated_files.append(out_file)
                    logger.info(f"NES条形图已生成: {out_file}")
                except Exception as e:
                    logger.error(f"NES条形图生成失败 ({db_name}): {e}")

            # GSEA 气泡图
            if 'dotplot' in valid_types and 'nes' in df.columns:
                out_file = str(plot_dir / f"{db_name}_dotplot.{plot_format}")
                try:
                    plot_gsea_dotplot(df, output_file=out_file, dpi=plot_dpi,
                                      style=plot_style, palette=plot_palette)
                    generated_files.append(out_file)
                    logger.info(f"GSEA气泡图已生成: {out_file}")
                except Exception as e:
                    logger.error(f"GSEA气泡图生成失败 ({db_name}): {e}")

            # GSEA 富集曲线图（需要 ranked_genes + gene_sets）
            if 'enrichment' in valid_types and ranked_genes and gene_sets:
                for set_name, gene_set in gene_sets.items():
                    # 从结果中查找该基因集的统计数据
                    match = None
                    if 'pathway' in df.columns:
                        match_rows = df[df['pathway'] == set_name]
                        if len(match_rows) > 0:
                            match = match_rows.iloc[0]

                    if match is None:
                        continue

                    es_val = match.get('es', 0.0)
                    nes_val = match.get('nes', 0.0)
                    pval = match.get('pvalue', match.get('P_Value', 1.0))

                    out_file = str(plot_dir / f"{set_name}_enrichment.{plot_format}")
                    try:
                        plot_gsea_enrichment(
                            ranked_genes=ranked_genes,
                            gene_weights=gene_weights or {},
                            gene_set=gene_set,
                            es=es_val,
                            nes=nes_val,
                            pvalue=pval,
                            title=set_name,
                            output_file=out_file,
                            dpi=plot_dpi,
                            style=plot_style,
                            palette=plot_palette,
                        )
                        generated_files.append(out_file)
                    except Exception as e:
                        logger.error(f"富集曲线图生成失败 ({set_name}): {e}")

    # ---- 通用图表类型（network, upset, volcano）----
    # 通用图表不依赖于特定的分析方法，可以在任何结果上生成
    common_types_requested = [pt for pt in plot_types if pt in _COMMON_PLOT_TYPES]
    if common_types_requested:
        from allenricher.visualization.common_plots import (
            plot_enrichment_network,
            plot_upset,
            plot_volcano,
        )

        # 确保通用图表输出目录存在
        common_plot_dir = Path(output_dir) / "common_plots"
        common_plot_dir.mkdir(parents=True, exist_ok=True)

        for db_name, df in results.items():
            if df is None or len(df) == 0:
                continue

            # Network 图（需要 gene_sets）
            if 'network' in common_types_requested and gene_sets:
                out_file = str(common_plot_dir / f"{db_name}_network.{plot_format}")
                try:
                    plot_enrichment_network(
                        gene_sets=gene_sets,
                        results_df=df,
                        output_file=out_file,
                        dpi=plot_dpi,
                        style=plot_style,
                        palette=plot_palette,
                    )
                    generated_files.append(out_file)
                    logger.info(f"通路网络图已生成: {out_file}")
                except Exception as e:
                    logger.error(f"通路网络图生成失败 ({db_name}): {e}")

            # Upset 图（需要 gene_sets）
            if 'upset' in common_types_requested and gene_sets:
                out_file = str(common_plot_dir / f"{db_name}_upset.{plot_format}")
                try:
                    plot_upset(
                        gene_sets=gene_sets,
                        output_file=out_file,
                        dpi=plot_dpi,
                        style=plot_style,
                        palette=plot_palette,
                    )
                    generated_files.append(out_file)
                    logger.info(f"UpSet图已生成: {out_file}")
                except Exception as e:
                    logger.error(f"UpSet图生成失败 ({db_name}): {e}")

            # Volcano 图（需要 nes 和 pvalue 列）
            if 'volcano' in common_types_requested:
                has_nes = 'nes' in df.columns or 'NES' in df.columns
                has_pval = 'pvalue' in df.columns or 'pvalue' in df.columns or 'P_Value' in df.columns
                if has_nes and has_pval:
                    out_file = str(common_plot_dir / f"{db_name}_volcano.{plot_format}")
                    try:
                        plot_volcano(
                            results_df=df,
                            output_file=out_file,
                            dpi=plot_dpi,
                            style=plot_style,
                            palette=plot_palette,
                        )
                        generated_files.append(out_file)
                        logger.info(f"火山图已生成: {out_file}")
                    except Exception as e:
                        logger.error(f"火山图生成失败 ({db_name}): {e}")

    # ---- ssGSEA / GSVA 图表 ----
    elif method in ('ssgsea', 'gsva'):
        from allenricher.visualization.gsva_plots import (
            plot_pathway_heatmap,
            plot_group_comparison,
            plot_pathway_dotplot,
            plot_sample_correlation,
        )

        # 从结果中提取活性得分矩阵
        # 结果 DataFrame 应为 行=通路, 列=样本 的活性得分
        scores_df = None
        for db_name, df in results.items():
            if df is not None and len(df) > 0:
                # 尝试将结果转为通路x样本的活性矩阵
                # 如果结果本身就是活性矩阵格式，直接使用
                if isinstance(df, pd.DataFrame):
                    # 检查是否包含数值列（样本列）
                    numeric_cols = df.select_dtypes(include='number').columns
                    non_metric_cols = {'pvalue', 'P_Value', 'Adjusted_P_Value',
                                       'nes', 'es', 'fdr', 'gene_count', 'Gene_Count'}
                    sample_cols = [c for c in numeric_cols if c not in non_metric_cols]
                    if sample_cols:
                        # 尝试使用第一列作为通路名
                        name_col = None
                        for col in ['pathway', 'Term_Name', 'Term_ID', df.index.name]:
                            if col and col in df.columns:
                                name_col = col
                                break
                        if name_col:
                            scores_df = df.set_index(name_col)[sample_cols]
                        else:
                            scores_df = df[sample_cols]
                            scores_df.index.name = 'pathway'
                    elif expr_matrix is not None:
                        # 结果不是活性矩阵，使用表达矩阵
                        scores_df = expr_matrix
                break

        if scores_df is None:
            logger.warning("无法从分析结果中提取活性得分矩阵，跳过可视化")
            return generated_files

        # 构建分组注释 DataFrame（用于热图注释）
        annotation_df = None
        if groups:
            sample_to_group = {}
            for group_name, samples in groups.items():
                for s in samples:
                    sample_to_group[s] = group_name
            common_samples = [s for s in scores_df.columns if s in sample_to_group]
            if common_samples:
                annotation_df = pd.DataFrame({
                    'Group': [sample_to_group[s] for s in common_samples]
                }, index=common_samples)

        # 热图
        if 'heatmap' in valid_types:
            out_file = str(plot_dir / f"activity_heatmap.{plot_format}")
            try:
                plot_pathway_heatmap(
                    scores_df,
                    annotation_col=annotation_df,
                    output_file=out_file,
                    dpi=plot_dpi,
                    style=plot_style,
                    palette=plot_palette,
                )
                generated_files.append(out_file)
                logger.info(f"活性热图已生成: {out_file}")
            except Exception as e:
                logger.error(f"活性热图生成失败: {e}")

        # 组间比较图
        if 'group_comparison' in valid_types and groups:
            out_file = str(plot_dir / f"group_comparison.{plot_format}")
            try:
                plot_group_comparison(
                    scores_df,
                    groups=groups,
                    output_file=out_file,
                    dpi=plot_dpi,
                    style=plot_style,
                    palette=plot_palette,
                )
                generated_files.append(out_file)
                logger.info(f"组间比较图已生成: {out_file}")
            except Exception as e:
                logger.error(f"组间比较图生成失败: {e}")

        # 气泡图
        if 'dotplot' in valid_types:
            out_file = str(plot_dir / f"activity_dotplot.{plot_format}")
            try:
                plot_pathway_dotplot(
                    scores_df,
                    groups=groups,
                    output_file=out_file,
                    dpi=plot_dpi,
                    style=plot_style,
                    palette=plot_palette,
                )
                generated_files.append(out_file)
                logger.info(f"活性气泡图已生成: {out_file}")
            except Exception as e:
                logger.error(f"活性气泡图生成失败: {e}")

        # 样本相关性热图
        if 'correlation' in valid_types:
            out_file = str(plot_dir / f"sample_correlation.{plot_format}")
            try:
                plot_sample_correlation(
                    scores_df,
                    annotation_col=annotation_df,
                    output_file=out_file,
                    dpi=plot_dpi,
                    style=plot_style,
                    palette=plot_palette,
                )
                generated_files.append(out_file)
                logger.info(f"样本相关性热图已生成: {out_file}")
            except Exception as e:
                logger.error(f"样本相关性热图生成失败: {e}")

    return generated_files


def create_parser() -> argparse.ArgumentParser:
    """创建命令行参数解析器

    构建包含以下子命令的参数解析器：
      - analyze:  富集分析子命令（支持多种参数配置）
      - download: 数据库下载子命令
      - build:    物种数据库构建子命令
      - serve:    API 服务器启动子命令
      - list:     资源列表查看子命令
      - config:   配置文件生成子命令

    Returns:
        argparse.ArgumentParser: 配置完成的参数解析器实例
    """
    # 创建顶层解析器
    parser = argparse.ArgumentParser(
        prog='allenricher',
        description='AllEnricher v2.0 - Gene Set Enrichment Analysis Tool',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog='''
Examples:
  # Basic analysis
  allenricher analyze -i genes.txt -s hsa -d GO,KEGG -o results/

  # With AI interpretation
  allenricher analyze -i genes.txt -s hsa --ai openai --ai-key YOUR_KEY

  # Download databases
  allenricher download -d GO,KEGG -s hsa

  # Start API server
  allenricher serve --port 8000
        '''
    )
    
    # 添加版本号参数，使用 -v 或 --version 可查看当前版本
    parser.add_argument('-v', '--version', action='version', version=f'%(prog)s {__version__}')
    
    # 创建子命令解析器组
    subparsers = parser.add_subparsers(dest='command', help='Available commands')
    
    # ==================== analyze 子命令 ====================
    # 运行基因集功能富集分析（主工作流），支持输入基因列表、选择物种和数据库、
    # 设置统计方法和多重检验校正方式，并可生成可视化图表和 AI 解读报告
    analyze_parser = subparsers.add_parser('analyze', help='Run enrichment analysis')
    analyze_parser.add_argument('-i', '--input', required=True, help='Input gene list file')           # 输入基因列表文件路径（必需）
    analyze_parser.add_argument('-s', '--species', default='hsa', help='Species code (default: hsa)')  # 物种代码，默认为人类(hsa)
    analyze_parser.add_argument('-d', '--databases', default='GO,KEGG', help='Comma-separated databases')  # 逗号分隔的数据库名称列表
    analyze_parser.add_argument('-o', '--output', default='./results', help='Output directory')        # 输出目录，默认为 ./results
    analyze_parser.add_argument('-b', '--background', help='Background gene list file')                # 背景基因列表文件（可选）
    analyze_parser.add_argument('--background-mode', dest='background_mode',
                                choices=['annotated', 'genome', 'custom'], default='annotated',
                                help='Background gene set mode: annotated (default), genome, custom')
    analyze_parser.add_argument('-m', '--method', default='fisher', choices=['fisher', 'hypergeometric', 'gsea', 'ssgsea', 'gsva'], help='Enrichment method')  # 富集分析方法：Fisher精确检验/超几何检验/GSEA/ssGSEA/GSVA
    analyze_parser.add_argument('-c', '--correction', default='BH', choices=['BH', 'BY', 'bonferroni', 'holm', 'none'], help='Multiple testing correction')  # 多重检验校正方法
    analyze_parser.add_argument('-p', '--pvalue', type=float, default=0.05, help='P-value cutoff')    # P 值阈值，默认 0.05
    analyze_parser.add_argument('-q', '--qvalue', type=float, default=0.05, help='Q-value cutoff')    # Q 值（校正后 P 值）阈值，默认 0.05
    analyze_parser.add_argument('-n', '--min-genes', type=int, default=2, help='Minimum genes per term')  # 每个功能条目最少包含的基因数
    analyze_parser.add_argument('-j', '--jobs', type=int, default=1, help='Number of parallel jobs')   # 并行任务数
    analyze_parser.add_argument('--no-plot', action='store_true', help='Skip plot generation')          # 跳过可视化图表生成
    analyze_parser.add_argument('--no-report', action='store_true', help='Skip report generation')      # 跳过 HTML 报告生成
    analyze_parser.add_argument('--only-significant', action='store_true', help='Only output significant terms (filter by p/q cutoff)')  # 仅输出显著条目（按 p/q 阈值过滤，默认不启用，输出全部条目）
    analyze_parser.add_argument('--ai', choices=['openai', 'claude', 'deepseek', 'glm', 'minimax', 'ollama', 'mock'],
                                help='AI backend for interpretation (override YAML config)')  # AI 解读后端选择
    analyze_parser.add_argument('--ai-key', help='AI API key (override YAML config, optional if set in YAML)')  # AI 服务 API 密钥
    analyze_parser.add_argument('--ai-model', help='AI model name (override YAML config)')  # AI 模型名称
    analyze_parser.add_argument('--config', help='Configuration file (YAML/JSON)')                      # 外部配置文件路径
    analyze_parser.add_argument('--database-dir', help='Database directory')                       # 数据库目录路径
    analyze_parser.add_argument('--use-version', type=str, default=None,
                                help='指定使用的数据库版本（如 v20260515），默认使用最新版本')
    analyze_parser.add_argument('-e', '--expression-matrix', default=None, help='Expression matrix file (TSV/CSV, rows=genes, cols=samples) for GSEA/ssGSEA/GSVA')  # 表达矩阵文件路径（行=基因，列=样本），用于GSEA/ssGSEA/GSVA
    analyze_parser.add_argument('-r', '--ranked-genes', default=None, help='Ranked gene list file (two columns: gene_name weight) for GSEA')  # 排序基因列表文件路径（两列: 基因名 权重），用于GSEA
    analyze_parser.add_argument('-g', '--gmt', default=None, help='GMT format gene set file path (supports .gmt and .gmt.gz)')  # GMT格式基因集文件路径，用于GSEA/ssGSEA/GSVA的基因集定义
    analyze_parser.add_argument('-pt', '--plot-types', default=None, help='Comma-separated plot types to generate. GSEA: enrichment,nes_barplot,dotplot; GSVA/ssGSEA: heatmap,group_comparison,dotplot,correlation')  # 要生成的图表类型，逗号分隔
    analyze_parser.add_argument('--groups', default=None, help='Sample group definition, format: Group1:sample1,sample2;Group2:sample3,sample4')  # 样本分组定义，用于组间比较图
    analyze_parser.add_argument('--plot-format', default='png', choices=['png', 'pdf', 'svg'], help='Plot output format (default: png)')  # 图表输出格式
    analyze_parser.add_argument('--plot-dpi', type=int, default=300, help='Plot resolution/DPI (default: 300)')  # 图表分辨率(DPI)
    analyze_parser.add_argument('--style', default='nature', choices=['nature', 'science', 'colorblind', 'presentation', 'omicshare'], help='Plot style theme (default: nature)')  # 图表风格主题
    analyze_parser.add_argument('--palette', default=None, help='Custom color palette name (optional)')  # 自定义配色方案
    analyze_parser.add_argument('--verbose', action='store_true', help='Enable verbose (DEBUG) logging')  # 启用详细日志输出
    
    # ==================== download 子命令 ====================
    # 从远程数据源下载指定的富集分析数据库到本地
    download_parser = subparsers.add_parser('download', help='Download databases')
    download_parser.add_argument('-d', '--databases', required=True, help='Comma-separated databases to download')  # 要下载的数据库名称（必需）
    download_parser.add_argument('-s', '--species', default='hsa', help='Species code')                # 物种代码
    download_parser.add_argument('--database-dir', default='./database', help='Database directory')     # 数据库存储目录
    download_parser.add_argument('--workers', type=int, default=4, help='Multi-thread download workers (default: 4)')  # 多线程下载数
    download_parser.add_argument('--no-multi-thread', action='store_true', help='Disable multi-thread download')       # 禁用多线程
    download_parser.add_argument('--no-verify', action='store_true', help='Skip post-download integrity check')        # 跳过完整性校验
    download_parser.add_argument('--force', action='store_true', help='强制重新下载，即使本地已是最新版本')
    
    # ==================== build 子命令 ====================
    # 为指定物种构建本地富集分析数据库，需要提供物种代码和分类学 ID
    build_parser = subparsers.add_parser('build', help='Build species database')
    build_parser.add_argument('-s', '--species', required=True, help='Species code')                    # 物种代码（必需）
    build_parser.add_argument('-t', '--taxonomy', required=True, type=int, help='Taxonomy ID')          # NCBI 分类学 ID（必需）
    build_parser.add_argument('-d', '--databases', default='GO,KEGG,Reactome', help='Comma-separated databases to build')  # 要构建的数据库列表
    build_parser.add_argument('--database-dir', default='./database', help='Database directory')        # 数据库存储目录
    build_parser.add_argument('--gene-info', help='Path to NCBI gene_info.gz file')                    # NCBI gene_info.gz 文件路径（GO和Reactome构建需要）
    
    # 自定义注释文件参数
    build_parser.add_argument('--go-annot',
        help='Path to GO annotation file (TSV: gene<TAB>go_id<TAB>go_name[<TAB>hierarchy])')
    build_parser.add_argument('--kegg-annot',
        help='Path to KEGG annotation file (TSV: gene<TAB>pathway_id<TAB>pathway_name[<TAB>hierarchy])')
    build_parser.add_argument('--custom-annot',
        help='Path to custom annotation file (TSV format with hierarchy support)')
    build_parser.add_argument('--custom-db-name', default='CUSTOM',
        help='Database name for custom annotation (default: CUSTOM)')
    build_parser.add_argument('--annot-format',
        choices=['three_column', 'four_column', 'two_column', 'auto'],
        default='auto', help='Annotation file format (default: auto-detect)')
    build_parser.add_argument('--hierarchy-sep', default='|',
        help='Hierarchy level separator (default: |)')
    
    # ==================== serve 子命令 ====================
    # 启动 RESTful API 服务器，提供在线富集分析服务
    serve_parser = subparsers.add_parser('serve', help='Start API server')
    serve_parser.add_argument('--host', default='0.0.0.0', help='Server host')                          # 服务器监听地址，默认 0.0.0.0
    serve_parser.add_argument('--port', type=int, default=8000, help='Server port')                     # 服务器监听端口，默认 8000
    serve_parser.add_argument('--reload', action='store_true', help='Enable auto-reload')                # 启用热重载（开发模式）
    
    # ==================== list 子命令 ====================
    # 列出支持的物种列表或可用的数据库资源
    list_parser = subparsers.add_parser('list', help='List available resources')
    list_parser.add_argument('resource', choices=['species', 'databases'], help='Resource to list')      # 要查看的资源类型：species（物种）或 databases（数据库）
    
    # ==================== config 子命令 ====================
    # 生成默认的 YAML 配置文件，用户可在此基础上修改
    config_parser = subparsers.add_parser('config', help='Generate configuration file')
    config_parser.add_argument('-o', '--output', default='allenricher.yaml', help='Output config file')  # 输出配置文件路径

    # ==================== check-update 子命令 ====================
    # 检查远程数据源是否有更新
    check_update_parser = subparsers.add_parser('check-update', help='检查远程数据源是否有更新')
    check_update_parser.add_argument('--database-dir', default=None, help='数据库目录路径')
    check_update_parser.add_argument('--json', action='store_true', help='以 JSON 格式输出结果')

    # ==================== cleanup 子命令 ====================
    # 清理旧版本的数据库文件
    cleanup_parser = subparsers.add_parser('cleanup', help='清理旧版本的数据库文件')
    cleanup_parser.add_argument('--keep', type=int, default=2, help='保留的最新版本数量（默认: 2）')
    cleanup_parser.add_argument('--dry-run', action='store_true', help='仅预览，不实际删除')
    cleanup_parser.add_argument('--database-dir', default=None, help='数据库目录路径')

    # ==================== list-versions 子命令 ====================
    # 查看本地已安装的数据库版本
    list_versions_parser = subparsers.add_parser('list-versions', help='查看本地已安装的数据库版本')
    list_versions_parser.add_argument('--database-dir', default=None, help='数据库目录路径')
    list_versions_parser.add_argument('--json', action='store_true', help='以 JSON 格式输出')
    list_versions_parser.add_argument('--lineage', action='store_true', help='显示构建血缘追踪')

    # ==================== list-species 子命令 ====================
    # 列出物种注册表中支持的物种信息
    list_species_parser = subparsers.add_parser('list-species', help='List supported species from registry')
    list_species_parser.add_argument('--go', action='store_true', default=False, help='Filter by GO support')
    list_species_parser.add_argument('--kegg', action='store_true', default=False, help='Filter by KEGG support')
    list_species_parser.add_argument('--reactome', action='store_true', default=False, help='Filter by Reactome support')
    list_species_parser.add_argument('--do', action='store_true', default=False, help='Filter by DO support')
    list_species_parser.add_argument('--format', choices=['table', 'tsv', 'json'], default='table', help='Output format (default: table)')
    list_species_parser.add_argument('--summary', action='store_true', default=False, help='Show summary statistics')

    # ==================== query-species 子命令 ====================
    # 查询特定物种的详细信息
    query_species_parser = subparsers.add_parser('query-species', help='Query species detail from registry')
    query_species_parser.add_argument('--taxid', type=int, default=None, help='Query by NCBI Taxonomy ID')
    query_species_parser.add_argument('--name', type=str, default=None, help='Query by Latin name')
    query_species_parser.add_argument('--kegg', type=str, default=None, help='Query by KEGG organism code')

    return parser


def _resolve_db_dir(args) -> str:
    """统一解析数据库目录路径"""
    if hasattr(args, 'database_dir') and args.database_dir:
        return args.database_dir
    return "./database"


def cmd_analyze(args) -> int:
    """运行富集分析（主工作流）

    这是 AllEnricher 的核心命令处理函数，执行完整的富集分析流程：
      1. 加载或创建配置对象
      2. 验证配置参数的合法性
      3. 创建输出目录
      4. 读取输入基因列表
      5. 加载指定的富集分析数据库
      6. 确定背景基因集
      7. 执行富集分析计算
      8. 保存分析结果
      9. 生成可视化图表（可选）
     10. 生成 AI 解读报告（可选）
     11. 生成 HTML 综合报告（可选）

    Args:
        args: 命令行参数命名空间，包含 analyze 子命令的所有参数

    Returns:
        int: 0 表示成功，1 表示失败
    """
    # ---- 设置日志详细程度 ----
    # 如果用户指定了 --verbose，将日志级别提升为 DEBUG，输出更详细的信息
    if args.verbose:
        logging.getLogger('allenricher').setLevel(logging.DEBUG)
        logger.debug("已启用详细日志模式（DEBUG 级别）")
    
    try:
        logger.info(f"AllEnricher v{__version__} - Starting analysis")
        
        # ---- 第1步：加载或创建配置 ----
        # 如果用户提供了外部配置文件，则从文件加载；否则根据命令行参数创建配置对象
        if args.config:
            config = Config.from_file(args.config)  # 从 YAML/JSON 配置文件加载
            # 命令行参数优先级高于配置文件（如果命令行显式指定了参数则覆盖配置文件中的值）
            if args.input:
                config.input_file = args.input
            if args.species != 'hsa':
                config.species = args.species
            if args.databases != 'GO,KEGG':
                config.databases = args.databases.split(',')
            if args.method != 'fisher':
                config.method = args.method
            if args.correction != 'BH':
                config.correction = args.correction
            if args.pvalue != 0.05:
                config.pvalue_cutoff = args.pvalue
            if args.qvalue != 0.05:
                config.qvalue_cutoff = args.qvalue
            if args.min_genes != 2:
                config.min_genes = args.min_genes
            if args.jobs != 1:
                config.n_jobs = args.jobs
            if args.output != './results':
                config.output_dir = args.output
            if args.background:
                config.background_file = args.background
            # CLI 标志覆盖 output_all（默认 True，--only-significant 时设为 False）
            if args.only_significant:
                config.output_all = False
            # 图表风格参数
            if hasattr(args, 'style') and args.style:
                config.plot_style = args.style
            if hasattr(args, 'palette') and args.palette:
                config.plot_palette = args.palette
        else:
            config = Config(
                input_file=args.input,
                output_dir=args.output,
                species=args.species,
                databases=args.databases.split(','),  # 将逗号分隔的字符串拆分为列表
                method=args.method,
                correction=args.correction,
                pvalue_cutoff=args.pvalue,
                qvalue_cutoff=args.qvalue,
                min_genes=args.min_genes,
                n_jobs=args.jobs,
                background_file=args.background,
                output_all=not args.only_significant,  # 默认输出全部条目（与v1一致）
                plot_style=getattr(args, 'style', 'nature'),
                plot_palette=getattr(args, 'palette', None),
            )
        
        # ---- 第2步：验证配置 ----
        # 检查配置参数是否合法，如输入文件是否存在、物种代码是否有效等
        errors = config.validate()
        if errors:
            for error in errors:
                logger.error(f"Configuration error: {error}")
            return 1  # 配置验证失败，返回错误码 1
        
        # ---- 第3步：创建输出目录 ----
        # 如果输出目录不存在则自动创建（包括父目录）
        output_dir = Path(config.output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        
        # ---- 第4步：加载基因列表 ----
        # 从输入文件中读取待分析的基因列表（优先使用配置文件中的 input_file，其次使用命令行参数）
        input_path = config.input_file or args.input
        if not input_path:
            logger.error("未指定输入基因列表文件！请通过 -i/--input 参数或配置文件指定。")
            return 1
        logger.info(f"Loading gene list from {input_path}")
        analyzer = EnrichmentAnalyzer(config)
        gene_set = analyzer.load_gene_list(input_path)
        
        # ---- 第4.5步：加载表达矩阵或排序基因列表（GSEA/ssGSEA/GSVA专用） ----
        # 如果使用 GSEA/ssGSEA/GSVA 方法，可以提供表达矩阵或排序基因列表
        expression_matrix = None
        ranked_gene_list = None
        if hasattr(args, 'expression_matrix') and args.expression_matrix:
            logger.info(f"Loading expression matrix from {args.expression_matrix}")
            import pandas as pd
            expr_path = Path(args.expression_matrix)
            if expr_path.suffix.lower() in ('.csv',):
                expression_matrix = pd.read_csv(expr_path, index_col=0)
            else:
                # 默认按TSV格式读取
                expression_matrix = pd.read_csv(expr_path, sep='\t', index_col=0)
            logger.info(f"Expression matrix loaded: {expression_matrix.shape[0]} genes x {expression_matrix.shape[1]} samples")
        if hasattr(args, 'ranked_genes') and args.ranked_genes:
            logger.info(f"Loading ranked gene list from {args.ranked_genes}")
            ranked_gene_list = analyzer.load_ranked_gene_list(args.ranked_genes)
            logger.info(f"Ranked gene list loaded: {len(ranked_gene_list)} genes")
        
        # ---- 第4.6步：加载GMT基因集文件（GSEA/ssGSEA/GSVA可视化专用） ----
        gene_sets = None
        if hasattr(args, 'gmt') and args.gmt:
            gene_sets = _parse_gmt_file(args.gmt)
        
        # ---- 第4.7步：解析可视化参数 ----
        # 解析 --plot-types 逗号分隔字符串为列表
        plot_types_list = None
        if hasattr(args, 'plot_types') and args.plot_types:
            plot_types_list = [pt.strip() for pt in args.plot_types.split(',') if pt.strip()]
            logger.info(f"请求生成的图表类型: {plot_types_list}")
        
        # 解析 --groups 分组字符串为字典
        groups_dict = None
        if hasattr(args, 'groups') and args.groups:
            groups_dict = _parse_groups(args.groups)
        
        # 获取图表输出格式和DPI
        plot_format = getattr(args, 'plot_format', 'png')
        plot_dpi = getattr(args, 'plot_dpi', 300)
        
        # ---- 第5步：加载富集分析数据库 ----
        # 根据配置中指定的数据库名称加载对应的数据库数据
        db_dir = args.database_dir if args.database_dir else config.database_dir
        logger.info(f"Loading databases: {config.databases} from {db_dir}")
        db_manager = DatabaseManager(db_dir, config.species)
        # 版本锁定：CLI --use-version 优先于 YAML config.use_version
        use_ver = getattr(args, 'use_version', None) or getattr(config, 'use_version', None)
        if use_ver:
            logger.info(f"使用指定数据库版本: {use_ver}")
        db_manager.load_databases(config.databases, version=use_ver)
        
        # ---- 第6步：确定背景基因集 ----
        # 根据 --background 和 --background-mode 参数确定背景基因集
        background_mode = getattr(args, 'background_mode', 'annotated')

        if config.background_file:
            # 用户提供了 --background 参数，直接使用（忽略 background_mode）
            logger.info(f"Loading background genes from {config.background_file}")
            background_set = analyzer.load_gene_list(config.background_file)
        elif background_mode == 'annotated':
            # 使用注释基因作为背景
            logger.info("Using annotated genes as background (background_mode='annotated')")
            background_set = db_manager.get_background_genes()
        elif background_mode == 'genome':
            # 使用全基因组基因作为背景（来自 gene_info.gz）
            logger.info("Using genome genes as background (background_mode='genome')")
            try:
                # 直接使用物种代码获取全基因组基因
                # manager.py 中的 KEGG_CODE_TO_TAXID 映射会自动处理
                background_set = db_manager.get_genome_genes(species_code=config.species)
                if background_set:
                    logger.info(f"Loaded {len(background_set)} genome genes from gene_info.gz")
                else:
                    logger.warning("No genome genes found in gene_info.gz, falling back to annotated genes")
                    background_set = db_manager.get_background_genes()
            except Exception as e:
                logger.warning(f"Failed to load genome genes: {e}, falling back to annotated genes")
                background_set = db_manager.get_background_genes()
        elif background_mode == 'custom':
            # custom 模式要求必须提供 --background
            logger.error("Background mode 'custom' requires --background parameter")
            return 1
        else:
            # 默认回退
            logger.info("Using all database genes as background")
            background_set = db_manager.get_background_genes()
        
        # ---- 第7步：执行富集分析 ----
        # 对每个数据库运行富集分析，支持并行计算（当 n_jobs > 1 时）
        logger.info("Running enrichment analysis...")
        database_data = db_manager.get_all_term_data()
        results = analyzer.run_analysis(gene_set, background_set, database_data, parallel=config.n_jobs > 1)
        
        # ---- 检查是否有富集结果 ----
        # 如果所有数据库都没有找到显著富集的结果，输出友好的提示信息并正常退出
        if not results or len(results) == 0:
            logger.warning("=" * 60)
            logger.warning("未找到显著富集的结果！")
            logger.warning("可能的原因：")
            logger.warning("  1. 输入基因列表过小或与数据库无交集")
            logger.warning("  2. p 值/q 值阈值过于严格")
            logger.warning("  3. 背景基因集设置不当")
            logger.warning("建议：")
            logger.warning("  - 增加输入基因数量")
            logger.warning("  - 放宽 p 值/q 值阈值（如 -p 0.1 -q 0.1）")
            logger.warning("  - 检查基因 ID 格式是否与数据库匹配")
            logger.warning("=" * 60)
            return 0  # 正常退出，不报错
        
        # ---- 第8步：保存分析结果 ----
        # 将富集分析结果保存到输出目录（保存全部条目）
        logger.info("Saving results...")

        # 构建元数据字典，记录版本和分析信息
        from datetime import datetime, timezone
        metadata = {
            "allenricher_version": __version__,
            "analysis_date": datetime.now(timezone.utc).isoformat(),
            "database_version": db_manager.active_version or "unknown",
            "species": config.species,
            "databases": config.databases,
        }
        build_meta = db_manager.get_build_metadata()
        if build_meta:
            metadata["source_versions"] = build_meta.get("source_versions", {})
            metadata["built_at"] = build_meta.get("built_at", "")

        analyzer.save_results(str(output_dir), metadata=metadata)
        
        # ---- 第8.5步：筛选显著结果用于绘图和报告 ----
        # 按 p/q 阈值过滤，只保留显著富集的条目用于后续可视化和报告
        significant_results = {}
        for db_name, df in results.items():
            if len(df) == 0:
                continue
            filtered = df[
                (df['P_Value'] <= config.pvalue_cutoff) &
                (df['Adjusted_P_Value'] <= config.qvalue_cutoff)
            ].copy()
            if len(filtered) > 0:
                significant_results[db_name] = filtered
        
        sig_total = sum(len(df) for df in significant_results.values())
        all_total = sum(len(df) for df in results.values())
        logger.info(f"Significant results: {sig_total}/{all_total} terms (p<={config.pvalue_cutoff}, q<={config.qvalue_cutoff})")
        
        # ---- 第9步：生成可视化图表 ----
        # 除非用户指定 --no-plot，否则为每个数据库的显著结果生成图表
        if not args.no_plot and significant_results:
            logger.info("Generating plots (significant results only)...")
            plotter = Plotter(str(output_dir / "plots"), config)
            for db_name, df in significant_results.items():
                if len(df) > 0:
                    plotter.plot_all(df, db_name, top_n=config.top_terms,
                                     style=config.plot_style, palette=config.plot_palette)
        
        # ---- 第9.5步：生成GSEA/GSVA/ssGSEA专用可视化图表 ----
        # 如果用户指定了 --plot-types 且方法为 GSEA/ssGSEA/GSVA，则调用专用可视化函数
        gsea_plot_files = []
        if (plot_types_list and config.method in _METHOD_PLOT_TYPES
                and not args.no_plot and results):
            logger.info(f"Generating {config.method} specific plots...")
            # 构建基因权重字典（GSEA富集曲线需要）
            gene_weights = None
            if ranked_gene_list and hasattr(ranked_gene_list, 'items'):
                gene_weights = dict(ranked_gene_list)
            elif ranked_gene_list and isinstance(ranked_gene_list, list):
                # ranked_gene_list 为 [(gene, weight), ...] 列表
                gene_weights = {g: w for g, w in ranked_gene_list}
            
            gsea_plot_files = _generate_plots(
                method=config.method,
                results=results,
                ranked_genes=list(ranked_gene_list) if ranked_gene_list else None,
                gene_weights=gene_weights,
                gene_sets=gene_sets,
                expr_matrix=expression_matrix,
                groups=groups_dict,
                plot_types=plot_types_list,
                output_dir=str(output_dir),
                plot_format=plot_format,
                plot_dpi=plot_dpi,
                plot_style=config.plot_style,
                plot_palette=config.plot_palette,
            )
            if gsea_plot_files:
                logger.info(f"GSEA/GSVA/ssGSEA 可视化完成: {len(gsea_plot_files)} 个图表")
        
        # ---- 第10步：生成 AI 解读报告 ----
        # 如果用户指定了 AI 后端（命令行或YAML配置），则调用 AI 模型对分析结果进行智能解读
        ai_interpretation = None

        # 确定是否启用AI解读：命令行 --ai > YAML ai_interpretation
        ai_enabled = args.ai or (config.ai_interpretation and config.ai_backend)

        if ai_enabled and significant_results:
            # 确定使用的后端：命令行 --ai > YAML ai_backend
            ai_backend = args.ai or config.ai_backend

            logger.info(f"Generating AI interpretation using {ai_backend}...")

            # 如果命令行提供了 --ai-key 或 --ai-model，使用传统方式（命令行参数优先）
            if args.ai_key or (args.ai_model and not config.ai_backends):
                interpreter = create_interpreter(
                    backend=ai_backend,
                    api_key=args.ai_key,
                    model=args.ai_model
                )
            else:
                # 从Config对象创建解释器（支持YAML ai_backends配置）
                interpreter = create_interpreter_from_config(config, backend=ai_backend)

            ai_interpretation = interpreter.interpret_results(significant_results)

            # 将 AI 解读结果保存为 JSON 文件
            import json
            with open(output_dir / "ai_interpretation.json", 'w') as f:
                json.dump(ai_interpretation, f, indent=2)
        
        # ---- 第11步：生成 HTML 综合报告 ----
        # 除非用户指定 --no-report，否则基于显著结果生成 HTML 报告
        if not args.no_report and significant_results:
            logger.info("Generating HTML report (significant results only)...")
            report_gen = ReportGenerator(str(output_dir), config)
            report_gen.generate(
                significant_results,
                str(output_dir / "report.html"),
                gene_list=list(gene_set),
                ai_interpretation=ai_interpretation,
                metadata=metadata
            )
        
        # ---- 打印分析摘要 ----
        logger.info("=" * 50)
        logger.info("Analysis Complete!")
        logger.info("=" * 50)
        for db_name, df in results.items():
            logger.info(f"  {db_name}: {len(df)} enriched terms")
        logger.info(f"Results saved to: {output_dir}")
        
        return 0
    
    except FileNotFoundError as e:
        # 文件未找到错误：输入文件、背景文件或数据库文件不存在
        logger.error(f"找不到文件，请检查路径: {e}")
        return 1
    
    except ValueError as e:
        # 参数错误：配置参数不合法、文件内容为空等
        logger.error(f"参数错误，请检查输入: {e}")
        return 1
    
    except KeyboardInterrupt:
        # 用户中断（Ctrl+C）
        logger.warning("用户中断了分析流程（Ctrl+C）")
        return 130  # Unix 惯例：128 + SIGINT(2) = 130
    
    except ImportError as e:
        # 依赖库缺失错误
        logger.error(f"缺少必要的依赖库: {e}")
        logger.error("请尝试执行: pip install allenricher[all] 安装所有依赖")
        return 1
    
    except Exception as e:
        # 通用异常捕获：记录完整的错误信息（包括堆栈跟踪）
        logger.error(f"分析过程中发生未预期的错误: {e}", exc_info=True)
        return 1


def cmd_download(args) -> int:
    """下载全体物种通用数据库

    下载步骤（对应 v1 的 update_GOdb / update_ReactomeDB）：
    - go:    下载 gene2go.gz + gene_info.gz + go-basic.obo → database/basic/go/GO{date}/
    - reactome: 下载 NCBI2Reactome + gene_info.gz → database/basic/reactome/Reactome{date}/
    - do:    下载 Jensen Lab disease TSV → database/basic/do/
    - disgenet: 下载 DisGeNET TSV → database/basic/disgenet/

    下载的是全体物种的通用原始数据，不区分物种。
    后续用 build 命令从中提取指定物种的数据。

    Args:
        args: 命令行参数命名空间

    Returns:
        int: 0 表示下载成功
    """
    from allenricher.database.downloader import DataDownloader

    databases = [d.strip().lower() for d in args.databases.split(',')]
    download_dir = args.database_dir or "./database"

    # ---- 非强制模式下先检查更新 ----
    if not getattr(args, 'force', False):
        try:
            from allenricher.database.version import RemoteVersionChecker, DatabaseVersionManager
            checker = RemoteVersionChecker()
            vm = DatabaseVersionManager(download_dir)
            update_status = checker.check_updates(vm)

            # 检查是否有任何数据源需要更新
            sources_with_update = [s for s, info in update_status.items() if info["has_update"]]
            sources_checked = list(update_status.keys())

            if sources_checked and not sources_with_update:
                logger.info("所有数据源均为最新版本，无需重新下载。")
                logger.info("如需强制重新下载，请使用 --force 参数。")
                return 0
            elif sources_with_update:
                logger.info(f"以下数据源有更新: {', '.join(sources_with_update)}")
                logger.info("继续下载...")
        except Exception as e:
            logger.warning(f"更新检查失败，继续执行下载: {e}")

    logger.info(f"下载全体物种通用数据 → {download_dir}")
    logger.info(f"数据库类型: {', '.join(databases)}")

    downloader = DataDownloader(
        root_dir=download_dir,
        max_workers=getattr(args, 'workers', 4),
        use_multi_thread=not getattr(args, 'no_multi_thread', False),
        verify_integrity=not getattr(args, 'no_verify', False),
    )

    try:
        downloaded = downloader.download_all(databases)
        for db_type, path in downloaded.items():
            logger.info(f"  ✅ {db_type}: {path}")
        logger.info("全部下载完成！")
        logger.info("")
        logger.info("下一步：构建指定物种的数据库")
        logger.info("  allenricher build -s hsa -t 9606 -d GO,Reactome")
        return 0
    except Exception as e:
        logger.error(f"下载失败: {e}")
        return 1


def cmd_build(args) -> int:
    """构建指定物种的数据库

    构建步骤（对应 v1 的 make_speciesDB）：
    从 database/basic/ 中的全体物种通用数据中，
    提取指定物种的数据，格式化输出到 database/organism/v{date}/{species}/。

    Args:
        args: 命令行参数命名空间

    Returns:
        int: 0 表示构建成功
    """
    from allenricher.database.builder import DatabaseBuilder

    databases = [d.strip().upper() for d in args.databases.split(',')]
    build_dir = args.database_dir or "./database"

    species = args.species
    taxid = args.taxonomy

    # ---- 自定义注释文件构建（优先于标准构建流程） ----
    has_custom_annotation = getattr(args, 'go_annot', None) or \
                           getattr(args, 'kegg_annot', None) or \
                           getattr(args, 'custom_annot', None)

    if has_custom_annotation:
        try:
            from allenricher.database.custom_builder import CustomDatabaseBuilder
            custom_builder = CustomDatabaseBuilder(root_dir=build_dir)

            if getattr(args, 'go_annot', None):
                fmt = None if getattr(args, 'annot_format', 'auto') == 'auto' else args.annot_format
                custom_builder.build_from_annotation(
                    annotation_file=args.go_annot,
                    species=species,
                    taxid=taxid,
                    db_name='GO',
                    format_type=fmt,
                    hierarchy_separator=getattr(args, 'hierarchy_sep', '|')
                )

            if getattr(args, 'kegg_annot', None):
                fmt = None if getattr(args, 'annot_format', 'auto') == 'auto' else args.annot_format
                custom_builder.build_from_annotation(
                    annotation_file=args.kegg_annot,
                    species=species,
                    taxid=taxid,
                    db_name='KEGG',
                    format_type=fmt,
                    hierarchy_separator=getattr(args, 'hierarchy_sep', '|')
                )

            if getattr(args, 'custom_annot', None):
                db_name = getattr(args, 'custom_db_name', 'CUSTOM')
                fmt = None if getattr(args, 'annot_format', 'auto') == 'auto' else args.annot_format
                custom_builder.build_from_annotation(
                    annotation_file=args.custom_annot,
                    species=species,
                    taxid=taxid,
                    db_name=db_name,
                    format_type=fmt,
                    hierarchy_separator=getattr(args, 'hierarchy_sep', '|')
                )
        except ImportError:
            print("Warning: CustomDatabaseBuilder not available. Skipping custom annotation build.")

    # ---- 标准构建流程 ----
    # 尝试从 SpeciesRegistry 查询物种信息
    try:
        from .database.species_registry import SpeciesRegistry
        registry = SpeciesRegistry.load_default()
        entry = registry.query_by_taxid(taxid)
        if entry:
            logger.info(f"Species found in registry: {entry.latin_name} (TaxID: {entry.taxid})")
            logger.info(f"  Database support - GO: {'Yes' if entry.has_go else 'No'}, "
                       f"KEGG: {'Yes' if entry.has_kegg else 'No'}, "
                       f"DO: {'Yes' if entry.has_do else 'No'}, "
                       f"Reactome: {'Yes' if entry.has_reactome else 'No'}")
            # 使用注册表中的 latin_name
            species_display_name = entry.latin_name
        else:
            species_display_name = species
    except Exception:
        # 静默回退，使用命令行提供的物种代码
        species_display_name = species

    logger.info(f"构建物种专属数据库: {species_display_name} (TaxID: {taxid})")
    logger.info(f"数据库根目录: {build_dir}")
    logger.info(f"要构建的数据库: {', '.join(databases)}")

    builder = DatabaseBuilder(root_dir=build_dir)

    try:
        outdir = builder.build_species_db(
            species=species,
            taxid=taxid,
            databases=databases
        )
        logger.info(f"构建完成！输出目录: {outdir}")
        logger.info("")
        logger.info("下一步：运行富集分析")
        logger.info(f"  allenricher analyze -i genes.txt -s {args.species} --database-dir {outdir}")
        return 0
    except FileNotFoundError as e:
        logger.error(f"基础数据未找到: {e}")
        logger.info("请先运行 download 下载基础数据：")
        logger.info("  allenricher download -d go,reactome")
        return 1
    except Exception as e:
        logger.error(f"构建失败: {e}")
        return 1


def cmd_serve(args) -> int:
    """启动 API 服务器

    启动 RESTful API 服务器，提供在线富集分析服务。
    用户可通过 HTTP 接口提交基因列表并获取分析结果。

    Args:
        args: 命令行参数命名空间，包含 host、port、reload 等参数

    Returns:
        int: 0 表示服务器正常启动（注意：服务器运行期间此函数不会返回）
    """
    logger.info(f"Starting API server on {args.host}:{args.port}")
    
    # 延迟导入 API 服务器模块，避免在不需要时加载依赖
    from allenricher.api.server import start_api
    start_api(host=args.host, port=args.port)  # 启动服务器（阻塞调用）
    
    return 0


def cmd_list(args) -> int:
    """列出可用资源

    根据用户指定的资源类型，列出系统支持的物种列表或可用的数据库资源。
      - species:   显示所有支持的物种代码、名称和分类学 ID
      - databases: 显示所有支持的数据库类型

    Args:
        args: 命令行参数命名空间，包含 resource 参数（'species' 或 'databases'）

    Returns:
        int: 0 表示执行成功
    """
    if args.resource == 'species':
        # 列出支持的物种列表
        # 优先尝试从 SpeciesRegistry 加载，失败时回退到 SPECIES_CONFIGS
        registry_loaded = False
        try:
            from .database.species_registry import SpeciesRegistry
            registry = SpeciesRegistry.load_default()
            if registry and registry.entries:
                registry_loaded = True
                print("\nSupported Species (from SpeciesRegistry):")
                print("-" * 70)
                print(f"{'Code':<10} {'Name':<30} {'TaxID':<10} {'GO':<5} {'KEGG':<5}")
                print("-" * 70)
                for entry in registry.entries.values():
                    # 使用 kegg_code 作为代码，如果没有则显示 taxid
                    code = entry.kegg_code if entry.has_kegg else str(entry.taxid)
                    go_status = 'Yes' if entry.has_go else 'No'
                    kegg_status = 'Yes' if entry.has_kegg else 'No'
                    name = entry.latin_name if entry.latin_name else f"TaxID {entry.taxid}"
                    print(f"{code:<10} {name:<30} {entry.taxid:<10} {go_status:<5} {kegg_status:<5}")
        except Exception:
            # 静默回退到 SPECIES_CONFIGS
            pass

        if not registry_loaded:
            from allenricher.core.config import SPECIES_CONFIGS

            print("\nSupported Species:")
            print("-" * 50)
            print(f"{'Code':<10} {'Name':<25} {'Taxonomy ID':<12}")
            print("-" * 50)
            for code, config in SPECIES_CONFIGS.items():
                print(f"{code:<10} {config.display_name:<25} {config.taxonomy_id:<12}")
        
    elif args.resource == 'databases':
        # 列出支持的数据库类型
        print("\nSupported Databases:")
        print("-" * 50)
        for db_name in ['GO', 'KEGG', 'DO', 'Reactome', 'DisGeNET']:
            print(f"  - {db_name}")
    
    return 0


def cmd_config(args) -> int:
    """生成配置文件

    生成一个默认的 YAML 格式配置文件，用户可以在生成的配置文件基础上
    修改参数，然后通过 analyze 命令的 --config 参数加载使用。

    Args:
        args: 命令行参数命名空间，包含 output 参数（输出文件路径）

    Returns:
        int: 0 表示生成成功
    """
    # 创建默认配置对象并写入文件
    config = Config()
    config.to_file(args.output)  # 将配置序列化为 YAML 文件
    logger.info(f"Configuration file generated: {args.output}")
    return 0


def cmd_check_update(args) -> int:
    """检查远程数据源是否有更新

    检测所有远程数据源（gene2go, gene_info, go_obo, kegg, reactome 等）的版本，
    与本地 versions.json 中记录的版本比较，输出更新状态表格。

    Args:
        args: 命令行参数命名空间

    Returns:
        int: 0 表示执行成功
    """
    import json
    from allenricher.database.version import RemoteVersionChecker, DatabaseVersionManager

    db_dir = _resolve_db_dir(args)
    checker = RemoteVersionChecker()
    vm = DatabaseVersionManager(db_dir)

    logger.info(f"正在检查远程数据源更新 (数据库目录: {db_dir}) ...")
    update_status = checker.check_updates(vm)

    # 打印格式化表格
    print("\n远程数据源更新检查")
    print("=" * 80)
    print(f"  {'数据源':<20} {'状态':<10} {'本地版本':<25} {'远程版本'}")
    print(f"  {'-'*20} {'-'*10} {'-'*25} {'-'*25}")

    has_any_update = False
    for source, info in sorted(update_status.items()):
        if info["has_update"]:
            status = "🔄 有更新"
            has_any_update = True
        else:
            status = "✅ 已最新"
        local_ver = info["local"].get("remote_version") or info["local"].get("version") or "-"
        remote_ver = info["remote"].get("remote_version") or "-"
        print(f"  {source:<20} {status:<10} {local_ver:<25} {remote_ver}")

    print("=" * 80)

    if has_any_update:
        print("提示: 有数据源存在更新，可运行 `allenricher download --force` 重新下载。")
    else:
        print("所有数据源均为最新版本。")

    # 如果指定 --json，额外输出 JSON
    if args.json:
        print("\n--- JSON Output ---")
        print(json.dumps(update_status, indent=2, ensure_ascii=False))

    return 0


def cmd_cleanup(args) -> int:
    """清理旧版本的数据库文件

    扫描 database/basic/ 和 database/organism/ 目录，找出旧版本并删除。
    默认保留每个数据源最新的 2 个版本。

    Args:
        args: 命令行参数命名空间

    Returns:
        int: 0 表示执行成功
    """
    from allenricher.database.version import DatabaseVersionManager

    db_dir = _resolve_db_dir(args)
    vm = DatabaseVersionManager(db_dir)

    if args.dry_run:
        logger.info(f"[dry-run] 预览清理操作 (数据库目录: {db_dir}, 保留最新 {args.keep} 个版本)")
    else:
        logger.info(f"执行清理操作 (数据库目录: {db_dir}, 保留最新 {args.keep} 个版本)")

    removed = vm.remove_stale_versions(keep_count=args.keep, dry_run=args.dry_run)

    if not removed:
        print("没有需要清理的旧版本。")
        return 0

    print("\n清理结果:")
    print("-" * 60)
    total_count = 0
    for source, versions in removed.items():
        if versions:
            action = "将删除" if args.dry_run else "已删除"
            print(f"  [{source}] {action}: {', '.join(versions)}")
            total_count += len(versions)

    print("-" * 60)
    action = "将删除" if args.dry_run else "已删除"
    print(f"共 {action} {total_count} 个旧版本目录。")

    if args.dry_run:
        print("\n提示: 这是预览模式。去掉 --dry-run 参数以实际执行清理。")

    return 0


def cmd_list_versions(args) -> int:
    """查看本地已安装的数据库版本

    显示本地数据库目录中已安装的基础数据和物种数据库版本信息。

    Args:
        args: 命令行参数命名空间

    Returns:
        int: 0 表示执行成功
    """
    import json
    from allenricher.database.version import DatabaseVersionManager

    db_dir = _resolve_db_dir(args)
    vm = DatabaseVersionManager(db_dir)

    if args.json:
        print(json.dumps(vm.get_summary_json(), indent=2, ensure_ascii=False))
    elif args.lineage:
        print(vm.get_full_lineage_report())
    else:
        print(vm.get_summary_table())

    return 0


def main():
    """程序主入口函数（命令分发逻辑）

    解析命令行参数，根据用户输入的子命令将执行分发到对应的处理函数：
      - analyze  -> cmd_analyze()
      - download -> cmd_download()
      - build    -> cmd_build()
      - serve    -> cmd_serve()
      - list     -> cmd_list()
      - config   -> cmd_config()

    如果用户未指定任何子命令，则打印帮助信息并退出。

    Returns:
        int: 0 表示正常退出，1 表示执行出错
    """
    # 创建参数解析器并解析命令行参数
    parser = create_parser()
    args = parser.parse_args()
    
    # 如果用户未指定子命令，打印帮助信息后退出
    if args.command is None:
        parser.print_help()
        return 0
    
    # 构建子命令名称到处理函数的映射表
    commands = {
        'analyze': cmd_analyze,
        'download': cmd_download,
        'build': cmd_build,
        'serve': cmd_serve,
        'list': cmd_list,
        'config': cmd_config,
        'list-species': _cmd_list_species,
        'query-species': _cmd_query_species,
        'check-update': cmd_check_update,
        'cleanup': cmd_cleanup,
        'list-versions': cmd_list_versions,
    }
    
    # 根据用户输入的子命令查找并调用对应的处理函数
    handler = commands.get(args.command)
    if handler:
        return handler(args)  # 调用处理函数并返回其退出码
    else:
        parser.print_help()   # 未知命令，打印帮助信息
        return 1


def _cmd_list_species(args) -> int:
    """列出支持的物种"""
    from .database.species_registry import SpeciesRegistry
    import json

    registry = SpeciesRegistry.load_default()

    if args.summary:
        summary = registry.get_summary()
        print(f"\n{'='*50}")
        print("Species Registry Summary")
        print(f"{'='*50}")
        print(f"Total species: {summary['total_species']:,}")
        print(f"\nBy Database:")
        for db, stats in summary.items():
            if db != 'total_species':
                print(f"  - {db.upper()}: {stats['count']:,} species")
        print(f"{'='*50}\n")
        return 0

    entries = registry.filter_by_databases(
        go=args.go or None,
        kegg=args.kegg or None,
        reactome=args.reactome or None,
        do=args.do or None
    )

    if args.format == "table":
        print(f"{'TaxID':<10} {'Latin Name':<30} {'GO':<5} {'KEGG':<6} {'Reactome':<9} {'DO':<4}")
        print("-" * 70)
        for e in entries[:100]:
            print(f"{e.taxid:<10} {e.latin_name:<30} {'Y' if e.has_go else 'N':<5} {'Y' if e.has_kegg else 'N':<6} {'Y' if e.has_reactome else 'N':<9} {'Y' if e.has_do else 'N':<4}")
        if len(entries) > 100:
            print(f"... and {len(entries) - 100} more species")
    elif args.format == "tsv":
        print("taxid\tlatin_name\thas_go\thas_kegg\thas_reactome\thas_do")
        for e in entries:
            print(f"{e.taxid}\t{e.latin_name}\t{e.has_go}\t{e.has_kegg}\t{e.has_reactome}\t{e.has_do}")
    elif args.format == "json":
        data = [{"taxid": e.taxid, "latin_name": e.latin_name, "has_go": e.has_go,
                 "has_kegg": e.has_kegg, "has_reactome": e.has_reactome, "has_do": e.has_do}
                for e in entries]
        print(json.dumps(data, indent=2))

    return 0


def _cmd_query_species(args) -> int:
    """查询物种详细信息"""
    from .database.species_registry import SpeciesRegistry

    registry = SpeciesRegistry.load_default()

    entry = None
    if args.taxid:
        entry = registry.query_by_taxid(args.taxid)
    elif args.name:
        entries = registry.query_by_latin_name(args.name)
        if entries:
            entry = entries[0]
            if len(entries) > 1:
                print(f"Found {len(entries)} matches, showing first:")
    elif args.kegg:
        entry = registry.query_by_kegg_code(args.kegg)

    if not entry:
        print("Species not found in registry.")
        return 1

    detail = registry.get_species_detail(entry.taxid)

    print(f"\n{'='*60}")
    print(f"Species Information")
    print(f"{'='*60}")
    print(f"Taxonomy ID: {detail['taxid']}")
    print(f"Latin Name:  {detail['latin_name']}")
    if detail.get('common_name') and detail['common_name'] != '-':
        print(f"Common Name: {detail['common_name']}")
    print(f"\nDatabase Support:")
    print(f"-" * 60)

    for db_name, db_key in [("GO", "has_go"), ("KEGG", "has_kegg"), ("Reactome", "has_reactome"), ("DO", "has_do")]:
        if detail.get(db_key):
            print(f"\n{db_name}: Supported")
            if db_name == "GO" and detail.get('go_source') and detail['go_source'] != '-':
                print(f"  Source: {detail['go_source']}")
            if db_name == "KEGG" and detail.get('kegg_code') and detail['kegg_code'] != '-':
                print(f"  Code: {detail['kegg_code']} (source: {detail.get('kegg_code_source', '-')})")

    print(f"\nBuild Command:")
    print(f"  allenricher build --taxonomy {detail['taxid']}")
    print(f"{'='*60}\n")

    return 0


if __name__ == '__main__':
    sys.exit(main())
