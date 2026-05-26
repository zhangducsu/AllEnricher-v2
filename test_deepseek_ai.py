#!/usr/bin/env python3
"""
DeepSeek AI 解释功能测试脚本

使用 DeepSeek 模型对 GO 和 KEGG 富集分析结果进行 AI 解读，生成 HTML 报告。
"""

import os
import sys
from pathlib import Path

# 添加项目路径
project_root = Path(__file__).parent
sys.path.insert(0, str(project_root))

from allenricher.core.config import Config
from allenricher.core.enrichment import EnrichmentAnalyzer
from allenricher.database.manager import DatabaseManager
from allenricher.report.generator import ReportGenerator
from allenricher.ai.interpreter import create_interpreter

# 设置 DeepSeek API 密钥
os.environ['DEEPSEEK_API_KEY'] = 'sk-5857ffb7000c42c99f5b5c88ee1d1c51'


def main():
    """主函数"""
    print("=" * 70)
    print("DeepSeek AI 解释功能测试")
    print("=" * 70)
    print()

    # 配置参数
    input_file = project_root / "example_genes.txt"
    output_dir = project_root / "deepseek_test_output"
    species = "hsa"
    databases = ["GO", "KEGG"]
    
    # DeepSeek 配置
    ai_backend = "deepseek"
    ai_model = "deepseek-v4-flash"  # 使用用户指定的模型
    ai_api_key = "sk-5857ffb7000c42c99f5b5c88ee1d1c51"

    print(f"输入文件: {input_file}")
    print(f"输出目录: {output_dir}")
    print(f"物种: {species}")
    print(f"数据库: {databases}")
    print(f"AI 后端: {ai_backend}")
    print(f"AI 模型: {ai_model}")
    print()

    # 创建输出目录
    output_dir.mkdir(exist_ok=True)

    # 步骤 1: 加载基因列表
    print("步骤 1: 加载基因列表...")
    with open(input_file, 'r') as f:
        genes = [line.strip() for line in f if line.strip()]
    gene_set = set(genes)
    print(f"  ✓ 加载了 {len(gene_set)} 个基因")
    print()

    # 步骤 2: 配置和加载数据库
    print("步骤 2: 配置和加载数据库...")
    
    # 使用 v1 兼容的数据库路径
    db_path = project_root / "database" / "organism" / "v20190612" / "hsa"
    
    if not db_path.exists():
        print(f"  ✗ 数据库路径不存在: {db_path}")
        print("  尝试查找其他数据库位置...")
        
        # 尝试查找数据库
        possible_paths = [
            project_root / "database" / "organism" / "v20190612" / "hsa",
            Path("F:/OneDrive/Documents/TraeSOLO/AllEnricher/AllEnricher-v1/database/organism/v20190612/hsa"),
        ]
        
        for path in possible_paths:
            if path.exists():
                db_path = path
                print(f"  ✓ 找到数据库: {db_path}")
                break
    else:
        print(f"  ✓ 数据库路径: {db_path}")

    try:
        db_manager = DatabaseManager(str(db_path), species)
        print("  ✓ DatabaseManager 创建成功")
        
        # 加载数据库
        print("  加载数据库中...")
        db_manager.load_databases(databases)
        print(f"  ✓ 加载了数据库: {databases}")
        
        # 获取数据
        background_set = db_manager.get_background_genes()
        database_data = db_manager.get_all_term_data()
        print(f"  ✓ 背景基因数: {len(background_set)}")
        print()
        
    except Exception as e:
        print(f"  ✗ 数据库加载失败: {e}")
        print("\n将使用 Mock 模式继续测试...")
        database_data = {}
        background_set = set()
        db_manager = None
        print()

    # 步骤 3: 运行富集分析
    print("步骤 3: 运行富集分析...")
    try:
        config = Config(
            species=species,
            databases=databases,
            method="fisher",
            pvalue_cutoff=0.05,
            qvalue_cutoff=0.05,
            min_genes=2
        )
        
        analyzer = EnrichmentAnalyzer(config)
        print("  ✓ EnrichmentAnalyzer 创建成功")
        
        # 运行分析
        if database_data:
            results = analyzer.run_analysis(
                gene_set=gene_set,
                background_set=background_set,
                database_data=database_data,
                parallel=False
            )
            print(f"  ✓ 分析完成，生成了 {len(results)} 个数据库的结果")
            
            # 保存 TSV 结果
            analyzer.save_results(str(output_dir))
            print(f"  ✓ 结果已保存到: {output_dir}")
        else:
            results = {}
            print("  ⚠ 没有数据库数据，跳过富集分析")
            
        print()
        
    except Exception as e:
        print(f"  ✗ 富集分析失败: {e}")
        import traceback
        traceback.print_exc()
        results = {}
        print()

    # 步骤 4: 生成 AI 解读
    print("步骤 4: 生成 AI 解读...")
    print(f"  使用模型: {ai_model}")
    
    try:
        # 创建 DeepSeek 解释器
        interpreter = create_interpreter(
            backend=ai_backend,
            api_key=ai_api_key,
            model=ai_model
        )
        print(f"  ✓ {type(interpreter.interpreter).__name__} 创建成功")
        
        # 生成解读
        print("  正在调用 DeepSeek API...")
        interpretations = interpreter.interpret_results(results)
        print(f"  ✓ 生成了 {len(interpretations)} 个解读")
        
        # 显示解读结果
        for db_name, interpretation in interpretations.items():
            print(f"\n--- {db_name} AI 解读 ---")
            print(interpretation[:500] + "..." if len(interpretation) > 500 else interpretation)
            print()
            
    except Exception as e:
        print(f"  ✗ AI 解读失败: {e}")
        import traceback
        traceback.print_exc()
        interpretations = {}
        print()

    # 步骤 5: 生成 HTML 报告
    print("步骤 5: 生成 HTML 报告...")
    try:
        report_gen = ReportGenerator(str(output_dir), config)
        
        # 生成报告
        report_file = output_dir / "enrichment_report.html"
        html = report_gen.generate(
            results=results,
            output_file=str(report_file),
            gene_list=list(gene_set),
            ai_interpretation=interpretations
        )
        
        print(f"  ✓ HTML 报告已生成: {report_file}")
        print()
        
    except Exception as e:
        print(f"  ✗ HTML 报告生成失败: {e}")
        import traceback
        traceback.print_exc()
        print()

    # 总结
    print("=" * 70)
    print("测试完成！")
    print("=" * 70)
    print()
    print("输出文件:")
    print(f"  - 富集分析结果: {output_dir}/")
    print(f"  - HTML 报告: {output_dir}/enrichment_report.html")
    print()
    
    if interpretations:
        print("AI 解读摘要:")
        for db_name, interpretation in interpretations.items():
            print(f"\n  [{db_name}]")
            print(f"  {interpretation[:200]}...")
    print()


if __name__ == "__main__":
    main()
