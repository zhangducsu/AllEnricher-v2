#!/usr/bin/env python3
"""
AI 解释功能演示脚本

展示 AllEnricher 的 AI 解释功能，支持多种后端：
- MockInterpreter: 测试用模拟后端（无需 API 密钥）
- OpenAI: GPT-4/3.5（需要 API 密钥）
- Claude: Anthropic Claude（需要 API 密钥）
- DeepSeek: 国产大模型（需要 API 密钥）
- GLM: 智谱 AI（需要 API 密钥）
- MiniMax: MiniMax（需要 API 密钥）
- Ollama: 本地部署（需要本地 Ollama 服务）

使用方法:
    python demo_ai_interpreter.py
"""

import json
import pandas as pd
from pathlib import Path

# 添加项目路径
import sys
project_root = Path(__file__).parent
sys.path.insert(0, str(project_root))

from allenricher.ai.interpreter import (
    create_interpreter,
    AIInterpreter,
    get_available_backends
)


def create_sample_results():
    """创建示例富集分析结果"""
    # GO 富集结果
    go_data = {
        "Term_ID": ["GO:0008150", "GO:0009987", "GO:0009653", "GO:0044237", "GO:0005615"],
        "Term_Name": [
            "biological_process",
            "cellular_process",
            "anatomical_structure_morphogenesis",
            "cellular_metabolic_process",
            "extracellular_space"
        ],
        "Gene_Count": [150, 120, 80, 95, 45],
        "P_Value": [1.5e-10, 2.3e-8, 5.6e-6, 8.9e-5, 1.2e-4],
        "Adjusted_P_Value": [3.2e-8, 2.1e-6, 4.5e-4, 0.0032, 0.0089],
    }
    go_df = pd.DataFrame(go_data)

    # KEGG 富集结果
    kegg_data = {
        "Term_ID": ["hsa04110", "hsa04115", "hsa04210", "hsa04010", "hsa04012"],
        "Term_Name": [
            "cell_cycle",
            "p53_signaling_pathway",
            "apoptosis",
            "mapk_signaling_pathway",
            "erbb_signaling_pathway"
        ],
        "Gene_Count": [25, 18, 15, 22, 12],
        "P_Value": [3.2e-12, 5.6e-10, 1.2e-8, 4.5e-7, 8.9e-6],
        "Adjusted_P_Value": [8.9e-10, 5.4e-8, 6.2e-6, 9.8e-5, 0.00045],
    }
    kegg_df = pd.DataFrame(kegg_data)

    return {
        "GO_Biological_Process": go_df,
        "KEGG_Pathway": kegg_df
    }


def demo_mock_interpreter():
    """演示 MockInterpreter（无需 API 密钥）"""
    print("=" * 70)
    print("演示 1: MockInterpreter（测试用模拟后端）")
    print("=" * 70)
    print()

    # 创建解释器
    interpreter = create_interpreter("mock")
    print(f"✓ 创建成功: {interpreter.backend_name}")
    print()

    # 准备示例数据
    results = create_sample_results()
    print(f"✓ 加载示例数据: {len(results)} 个数据库")
    for db_name, df in results.items():
        print(f"  - {db_name}: {len(df)} 条富集条目")
    print()

    # 生成解读
    print("生成 AI 解读...")
    interpretations = interpreter.interpret_results(results)
    print()

    # 显示结果
    for db_name, interpretation in interpretations.items():
        print(f"\n--- {db_name} AI 解读 ---")
        print(interpretation)
        print()

    # 生成 HTML 报告段落
    print("\n--- HTML 报告段落预览 ---")
    html = interpreter.generate_report_section(results)
    print(html[:500] + "..." if len(html) > 500 else html)
    print()

    # 保存为 JSON
    output_dir = project_root / "demo_output"
    output_dir.mkdir(exist_ok=True)

    json_path = output_dir / "mock_interpretation.json"
    with open(json_path, 'w', encoding='utf-8') as f:
        json.dump(interpretations, f, indent=2, ensure_ascii=False)
    print(f"✓ 解读结果已保存: {json_path}")

    # 保存 HTML
    html_path = output_dir / "mock_interpretation.html"
    with open(html_path, 'w', encoding='utf-8') as f:
        f.write(html)
    print(f"✓ HTML 报告已保存: {html_path}")
    print()

    return interpretations


def demo_all_backends():
    """演示所有可用后端"""
    print("\n" + "=" * 70)
    print("演示 2: 所有可用的 AI 后端")
    print("=" * 70)
    print()

    backends = get_available_backends()
    print(f"支持的后端 ({len(backends)} 个):")
    for backend in backends:
        print(f"  - {backend}")
    print()

    # 演示创建不同后端
    results = create_sample_results()

    print("创建各后端解释器实例...")
    for backend in backends:
        try:
            kwargs = {}
            # MiniMax 需要 group_id
            if backend == "minimax":
                kwargs["group_id"] = "demo-group-id"

            interpreter = create_interpreter(backend, **kwargs)
            print(f"  ✓ {backend}: {type(interpreter.interpreter).__name__}")
        except Exception as e:
            print(f"  ✗ {backend}: {e}")


def demo_term_summaries():
    """演示生成单个条目的总结"""
    print("\n" + "=" * 70)
    print("演示 3: 生成单个条目的总结")
    print("=" * 70)
    print()

    interpreter = create_interpreter("mock")
    results = create_sample_results()

    # 生成包含条目总结的解读
    print("生成包含条目总结的解读...")
    interpretations = interpreter.interpret_results(
        results,
        include_term_summaries=True
    )
    print()

    # 显示条目总结
    for db_name, interpretation in interpretations.items():
        if f"{db_name}_term_summaries" in interpretations:
            term_summaries = interpretations[f"{db_name}_term_summaries"]
            print(f"\n--- {db_name} 条目总结 (前 3 个) ---")
            for i, (term_name, summary) in enumerate(term_summaries.items()):
                if i >= 3:
                    break
                print(f"\n[{term_name}]")
                print(summary)
            print()


def demo_multiple_results_scenarios():
    """演示不同结果数量的场景"""
    print("\n" + "=" * 70)
    print("演示 4: 不同结果数量的场景")
    print("=" * 70)
    print()

    interpreter = create_interpreter("mock")

    # 场景 1: 空结果
    print("场景 1: 空结果")
    print("-" * 40)
    empty_results = {"GO": pd.DataFrame()}
    interpretations = interpreter.interpret_results(empty_results)
    print(f"解读结果: {interpretations}")
    print()

    # 场景 2: 少量结果
    print("场景 2: 3 条结果")
    print("-" * 40)
    small_results = {
        "GO": pd.DataFrame({
            "Term_Name": [f"term_{i}" for i in range(3)],
            "P_Value": [1e-5, 1e-4, 1e-3],
            "Gene_Count": [10, 8, 5],
        })
    }
    interpretations = interpreter.interpret_results(small_results)
    print(f"解读包含 {len(small_results['GO'])} 个条目")
    print()

    # 场景 3: 超过 20 条（只展示前 20）
    print("场景 3: 25 条结果（只展示前 20）")
    print("-" * 40)
    many_results = {
        "GO": pd.DataFrame({
            "Term_Name": [f"term_{i}" for i in range(25)],
            "P_Value": [1e-5 * (i+1) for i in range(25)],
            "Gene_Count": [20 - i for i in range(25)],
        })
    }
    interpretations = interpreter.interpret_results(many_results)
    text = interpretations["GO"]
    # 检查是否只展示了前 20 条
    has_term_0 = "term_0" in text
    has_term_19 = "term_19" in text
    has_term_20 = "term_20" in text
    print(f"包含 term_0: {has_term_0}")
    print(f"包含 term_19: {has_term_19}")
    print(f"包含 term_20 (不应该): {has_term_20}")
    print()


def main():
    """主函数"""
    print("\n" + "=" * 70)
    print("AllEnricher AI 解释功能演示")
    print("=" * 70)
    print()

    # 演示 1: MockInterpreter
    demo_mock_interpreter()

    # 演示 2: 所有后端
    demo_all_backends()

    # 演示 3: 条目总结
    demo_term_summaries()

    # 演示 4: 不同结果场景
    demo_multiple_results_scenarios()

    # 总结
    print("\n" + "=" * 70)
    print("演示完成！")
    print("=" * 70)
    print()
    print("下一步:")
    print("1. 使用真实 API 密钥测试其他后端")
    print("2. 查看 demo_output/ 目录下的输出文件")
    print("3. 集成到富集分析流程中")
    print()
    print("环境变量配置:")
    print("  export OPENAI_API_KEY=sk-xxx")
    print("  export ANTHROPIC_API_KEY=sk-ant-xxx")
    print("  export DEEPSEEK_API_KEY=sk-xxx")
    print()


if __name__ == "__main__":
    main()
