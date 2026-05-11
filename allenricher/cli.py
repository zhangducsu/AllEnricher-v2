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
from typing import List, Optional

from allenricher import __version__
from allenricher.core.config import Config
from allenricher.core.enrichment import EnrichmentAnalyzer
from allenricher.database.manager import DatabaseManager
from allenricher.visualization.plotter import Plotter
from allenricher.report.generator import ReportGenerator
from allenricher.ai.interpreter import create_interpreter

# 配置日志输出格式
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


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
    analyze_parser.add_argument('-m', '--method', default='fisher', choices=['fisher', 'hypergeometric', 'gsea', 'ssgsea'], help='Enrichment method')  # 富集分析方法：Fisher精确检验/超几何检验/GSEA/ssGSEA
    analyze_parser.add_argument('-c', '--correction', default='BH', choices=['BH', 'BY', 'bonferroni', 'holm', 'none'], help='Multiple testing correction')  # 多重检验校正方法
    analyze_parser.add_argument('-p', '--pvalue', type=float, default=0.05, help='P-value cutoff')    # P 值阈值，默认 0.05
    analyze_parser.add_argument('-q', '--qvalue', type=float, default=0.05, help='Q-value cutoff')    # Q 值（校正后 P 值）阈值，默认 0.05
    analyze_parser.add_argument('-n', '--min-genes', type=int, default=2, help='Minimum genes per term')  # 每个功能条目最少包含的基因数
    analyze_parser.add_argument('-j', '--jobs', type=int, default=1, help='Number of parallel jobs')   # 并行任务数
    analyze_parser.add_argument('--no-plot', action='store_true', help='Skip plot generation')          # 跳过可视化图表生成
    analyze_parser.add_argument('--no-report', action='store_true', help='Skip report generation')      # 跳过 HTML 报告生成
    analyze_parser.add_argument('--ai', choices=['openai', 'claude', 'ollama', 'mock'], help='AI backend for interpretation')  # AI 解读后端选择
    analyze_parser.add_argument('--ai-key', help='AI API key')                                          # AI 服务 API 密钥
    analyze_parser.add_argument('--ai-model', help='AI model name')                                     # AI 模型名称
    analyze_parser.add_argument('--config', help='Configuration file (YAML/JSON)')                      # 外部配置文件路径
    analyze_parser.add_argument('--verbose', action='store_true', help='Enable verbose (DEBUG) logging')  # 启用详细日志输出
    
    # ==================== download 子命令 ====================
    # 从远程数据源下载指定的富集分析数据库到本地
    download_parser = subparsers.add_parser('download', help='Download databases')
    download_parser.add_argument('-d', '--databases', required=True, help='Comma-separated databases to download')  # 要下载的数据库名称（必需）
    download_parser.add_argument('-s', '--species', default='hsa', help='Species code')                # 物种代码
    download_parser.add_argument('--database-dir', default='./database', help='Database directory')     # 数据库存储目录
    
    # ==================== build 子命令 ====================
    # 为指定物种构建本地富集分析数据库，需要提供物种代码和分类学 ID
    build_parser = subparsers.add_parser('build', help='Build species database')
    build_parser.add_argument('-s', '--species', required=True, help='Species code')                    # 物种代码（必需）
    build_parser.add_argument('-t', '--taxonomy', required=True, type=int, help='Taxonomy ID')          # NCBI 分类学 ID（必需）
    build_parser.add_argument('-d', '--databases', default='GO,KEGG,Reactome', help='Comma-separated databases to build')  # 要构建的数据库列表
    build_parser.add_argument('--database-dir', default='./database', help='Database directory')        # 数据库存储目录
    build_parser.add_argument('--gene-info', help='Path to NCBI gene_info.gz file')                    # NCBI gene_info.gz 文件路径（GO和Reactome构建需要）
    
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
    
    return parser


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
                background_file=args.background
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
        
        # ---- 第5步：加载富集分析数据库 ----
        # 根据配置中指定的数据库名称加载对应的数据库数据
        logger.info(f"Loading databases: {config.databases}")
        db_manager = DatabaseManager(config.database_dir, config.species)
        db_manager.load_databases(config.databases)
        
        # ---- 第6步：确定背景基因集 ----
        # 如果用户提供了背景基因文件则使用之，否则使用数据库中所有基因作为背景
        if config.background_file:
            logger.info(f"Loading background genes from {config.background_file}")
            background_set = analyzer.load_gene_list(config.background_file)
        else:
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
        # 将富集分析结果保存到输出目录
        logger.info("Saving results...")
        analyzer.save_results(str(output_dir))
        
        # ---- 第9步：生成可视化图表 ----
        # 除非用户指定 --no-plot，否则为每个数据库的结果生成图表
        if not args.no_plot:
            logger.info("Generating plots...")
            plotter = Plotter(str(output_dir / "plots"), config)
            for db_name, df in results.items():
                if len(df) > 0:
                    plotter.plot_all(df, db_name, top_n=config.top_terms)
        
        # ---- 第10步：生成 AI 解读报告 ----
        # 如果用户指定了 AI 后端，则调用 AI 模型对分析结果进行智能解读
        ai_interpretation = None
        if args.ai:
            logger.info(f"Generating AI interpretation using {args.ai}...")
            interpreter = create_interpreter(
                backend=args.ai,
                api_key=args.ai_key,
                model=args.ai_model
            )
            ai_interpretation = interpreter.interpret_results(results)
            
            # 将 AI 解读结果保存为 JSON 文件
            import json
            with open(output_dir / "ai_interpretation.json", 'w') as f:
                json.dump(ai_interpretation, f, indent=2)
        
        # ---- 第11步：生成 HTML 综合报告 ----
        # 除非用户指定 --no-report，否则生成包含结果表格、图表和 AI 解读的 HTML 报告
        if not args.no_report:
            logger.info("Generating HTML report...")
            report_gen = ReportGenerator(str(output_dir), config)
            report_gen.generate(
                results,
                str(output_dir / "report.html"),
                gene_list=list(gene_set),
                ai_interpretation=ai_interpretation
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
    """下载数据库

    从远程数据源下载指定的富集分析数据库（如 GO、KEGG 等）到本地目录。
    支持一次下载多个数据库（通过逗号分隔指定）。

    Args:
        args: 命令行参数命名空间，包含 databases、species、database_dir 等参数

    Returns:
        int: 0 表示下载成功
    """
    logger.info(f"Downloading databases: {args.databases}")
    
    # 初始化数据库管理器，指定存储目录和目标物种
    db_manager = DatabaseManager(args.database_dir, args.species)
    databases = args.databases.split(',')  # 将逗号分隔的数据库名称拆分为列表
    
    # 逐个下载每个数据库
    for db_name in databases:
        logger.info(f"Downloading {db_name}...")
        db_manager.load_database(db_name)      # 加载数据库配置
        db = db_manager.databases[db_name]
        db.download()                           # 执行实际的下载操作
    
    logger.info("Download complete!")
    return 0


def cmd_build(args) -> int:
    """构建物种数据库

    为指定物种构建本地富集分析数据库（复刻v1.0的make_speciesDB功能）。
    需要提供物种代码和对应的 NCBI 分类学 ID（Taxonomy ID），
    以及 NCBI gene_info.gz 原始数据文件。

    构建流程：
      1. 解析用户指定的数据库列表
      2. 初始化 DatabaseManager
      3. 调用 build_databases 方法，依次构建各数据库
      4. 构建完成后数据库可直接用于富集分析

    Args:
        args: 命令行参数命名空间，包含 species、taxonomy、database_dir、databases、gene_info 等参数

    Returns:
        int: 0 表示构建成功
    """
    from allenricher.database.manager import DatabaseManager

    # 第1步：解析数据库列表
    databases = [d.strip().upper() for d in args.databases.split(',')]
    logger.info(f"Building databases for {args.species} (TaxID: {args.taxonomy})")
    logger.info(f"Databases to build: {', '.join(databases)}")

    # 第2步：初始化数据库管理器
    db_manager = DatabaseManager(args.database_dir, args.species)

    # 第3步：执行一键构建
    try:
        db_manager.build_databases(
            taxid=args.taxonomy,
            databases=databases,
            gene_info_file=args.gene_info
        )
        logger.info("Database build complete!")
        return 0
    except FileNotFoundError as e:
        logger.error(f"Required data file not found: {e}")
        logger.error("Please ensure NCBI gene_info.gz and other source files are available.")
        logger.error("You can download them from: https://ftp.ncbi.nlm.nih.gov/gene/DATA/")
        return 1
    except Exception as e:
        logger.error(f"Database build failed: {e}")
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
        from allenricher.core.config import SPECIES_CONFIGS
        
        print("\nSupported Species:")
        print("-" * 50)
        print(f"{'Code':<10} {'Name':<25} {'Taxonomy ID':<12}")
        print("-" * 50)
        for code, config in SPECIES_CONFIGS.items():
            print(f"{code:<10} {config.display_name:<25} {config.taxonomy_id:<12}")
        
    elif args.resource == 'databases':
        # 列出支持的数据库类型
        from allenricher.database.manager import DatabaseManager
        
        print("\nSupported Databases:")
        print("-" * 50)
        for db_name in DatabaseManager.DATABASE_CLASSES.keys():
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
        'config': cmd_config
    }
    
    # 根据用户输入的子命令查找并调用对应的处理函数
    handler = commands.get(args.command)
    if handler:
        return handler(args)  # 调用处理函数并返回其退出码
    else:
        parser.print_help()   # 未知命令，打印帮助信息
        return 1


if __name__ == '__main__':
    sys.exit(main())
