#!/usr/bin/env python3
"""
AI Explain function demonstration script

Show AllEnricher It's... AI Explain function, Support multiple backends: 
- MockInterpreter: Test with Simulation Backend (No need API Key)
- OpenAI: GPT-4/3.5 (Yes. API Key)
- Claude: Anthropic Claude (Yes. API Key)
- DeepSeek: Large national model (Yes. API Key)
- GLM: Brain spectrum AI (Yes. API Key)
- MiniMax: MiniMax (Yes. API Key)
- Ollama: Local deployment (Local Ollama Services)

Use of methods:
    python demo_ai_interpreter.py
"""

import json
import pandas as pd
from pathlib import Path

# Add Item Path
import sys
project_root = Path(__file__).parent
sys.path.insert(0, str(project_root))

from allenricher.ai.interpreter import (
    create_interpreter,
    AIInterpreter,
    get_available_backends
)


def create_sample_results():
    """Creates the example of the enrichment analysis"""
    # GO Fuzzy Results
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

    # KEGM Report
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
    """Mocter (no API key required)"""
    print("=" * 70)
    print("Presentation 1: MockInterpreter (test with simulation backend)")
    print("=" * 70)
    print()

    # Create interpreter
    interpreter = create_interpreter("mock")
    print(f"* Created successfully: {interpreter.backend_name}")
    print()

    # Prepare for illustrative data
    results = create_sample_results()
    print(f"* Load example data: {len(results)} database")
    for db_name, df in results.items():
        print(f"- {db_name}: {len(df)}Rich Entry")
    print()

    # Generate interpretation
    print("Generate AI read...")
    interpretations = interpreter.interpret_results(results)
    print()

    # Show results
    for db_name, interpretation in interpretations.items():
        print(f"\n--- {db_name}AI Interpretation---")
        print(interpretation)
        print()

    # Generate HTML Report Paragraph
    print("\n---Preview of HTML report paragraph---")
    html = interpreter.generate_report_section(results)
    print(html[:500] + "..." if len(html) > 500 else html)
    print()

    # Save as JSON
    output_dir = project_root / "demo_output"
    output_dir.mkdir(exist_ok=True)

    json_path = output_dir / "mock_interpretation.json"
    with open(json_path, 'w', encoding='utf-8') as f:
        json.dump(interpretations, f, indent=2, ensure_ascii=False)
    print(f"* The results of the reading have been kept: {json_path}")

    # Save HTML
    html_path = output_dir / "mock_interpretation.html"
    with open(html_path, 'w', encoding='utf-8') as f:
        f.write(html)
    print(f"* HTML report saved: {html_path}")
    print()

    return interpretations


def demo_all_backends():
    """Show all available backends"""
    print("\n" + "=" * 70)
    print("Presentation 2: All available AI backends")
    print("=" * 70)
    print()

    backends = get_available_backends()
    print(f"Backend supported ({len(backends)}):")
    for backend in backends:
        print(f"  - {backend}")
    print()

    # Demonstration to create different backends
    results = create_sample_results()

    print("Create examples of backend interpreters...")
    for backend in backends:
        try:
            kwargs = {}
            # MiniMax needs group_id
            if backend == "minimax":
                kwargs["group_id"] = "demo-group-id"

            interpreter = create_interpreter(backend, **kwargs)
            print(f"  ✓ {backend}: {type(interpreter.interpreter).__name__}")
        except Exception as e:
            print(f"  ✗ {backend}: {e}")


def demo_term_summaries():
    """Shows the summary of the creation of single entries"""
    print("\n" + "=" * 70)
    print("Presentation 3: Generate summary of single entries")
    print("=" * 70)
    print()

    interpreter = create_interpreter("mock")
    results = create_sample_results()

    # Generate an interpretation that includes the summary of entries
    print("Generate an interpretation that includes the summary of entries...")
    interpretations = interpreter.interpret_results(
        results,
        include_term_summaries=True
    )
    print()

    # Show entry summary
    for db_name, interpretation in interpretations.items():
        if f"{db_name}_term_summaries" in interpretations:
            term_summaries = interpretations[f"{db_name}_term_summaries"]
            print(f"\n--- {db_name}Summary of Entry (previous 3)---")
            for i, (term_name, summary) in enumerate(term_summaries.items()):
                if i >= 3:
                    break
                print(f"\n[{term_name}]")
                print(summary)
            print()


def demo_multiple_results_scenarios():
    """Showing different numbers of results"""
    print("\n" + "=" * 70)
    print("Presentation 4: Number of different outcomes")
    print("=" * 70)
    print()

    interpreter = create_interpreter("mock")

    # Scene 1: Empty result
    print("Scene 1: Empty result")
    print("-" * 40)
    empty_results = {"GO": pd.DataFrame()}
    interpretations = interpreter.interpret_results(empty_results)
    print(f"Results of the reading: {interpretations}")
    print()

    # Scenario 2: a small result table.
    print("Scenario 2: 3 results")
    print("-" * 40)
    small_results = {
        "GO": pd.DataFrame({
            "Term_Name": [f"term_{i}" for i in range(3)],
            "P_Value": [1e-5, 1e-4, 1e-3],
            "Gene_Count": [10, 8, 5],
        })
    }
    interpretations = interpreter.interpret_results(small_results)
    print(f"Interpreted {len(small_results['GO'])} terms")
    print()

    # Scenario 3: more than 20 terms; interpretation is capped at 20.
    print("Scenario 3: 25 results (first 20 interpreted)")
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
    # Check if only 20 prior ones are displayed
    has_term_0 = "term_0" in text
    has_term_19 = "term_19" in text
    has_term_20 = "term_20" in text
    print(f"Includes term_0: {has_term_0}")
    print(f"Contains term_19: {has_term_19}")
    print(f"Include term_20 (should not): {has_term_20}")
    print()


def main():
    """Main Function"""
    print("\n" + "=" * 70)
    print("AllEnricher AI Explain function demo")
    print("=" * 70)
    print()

    # Demo 1: MockInterpreter
    demo_mock_interpreter()

    # Presentation 2: All Backends
    demo_all_backends()

    # Presentation 3: Summary of Entry
    demo_term_summaries()

    # Presentation 4: Different Results Scene
    demo_multiple_results_scenarios()

    # Summary
    print("\n" + "=" * 70)
    print("Show's done!")
    print("=" * 70)
    print()
    print("Next:")
    print("1. Test other backends with real API keys")
    print("View output files under demo_output/ directory")
    print("3. Integration into the enrichment analysis process")
    print()
    print("Environment variable configuration:")
    print("  export OPENAI_API_KEY=sk-xxx")
    print("  export ANTHROPIC_API_KEY=sk-ant-xxx")
    print("  export DEEPSEEK_API_KEY=sk-xxx")
    print()


if __name__ == "__main__":
    main()
