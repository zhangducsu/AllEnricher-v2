#!/usr/bin/env python3
"""
DeepSeek AI Explain function test script

Use DeepSeek Models. GO and KEGG The results of the enrichment analysis are being conducted AI Interpretation, Generate HTML Report.
"""

import os
import sys
from pathlib import Path

# Add Item Path
project_root = Path(__file__).parent
sys.path.insert(0, str(project_root))

from allenricher.core.config import Config
from allenricher.core.enrichment import EnrichmentAnalyzer
from allenricher.database.manager import DatabaseManager
from allenricher.report.generator import ReportGenerator
from allenricher.ai.interpreter import create_interpreter

# Set the DeepSeek API key
os.environ['DEEPSEEK_API_KEY'] = 'sk-5857ffb7000c42c99f5b5c88ee1d1c51'


def main():
    """Main Function"""
    print("=" * 70)
    print("DeepSeek AI Explanatorys Test")
    print("=" * 70)
    print()

    # Configure Parameters
    input_file = project_root / "example_genes.txt"
    output_dir = project_root / "deepseek_test_output"
    species = "hsa"
    databases = ["GO", "KEGG"]
    
    # DeepSeek Configuration
    ai_backend = "deepseek"
    ai_model = "deepseek-v4-flash"  # Use user-defined models
    ai_api_key = "sk-5857ffb7000c42c99f5b5c88ee1d1c51"

    print(f"Enter a file: {input_file}")
    print(f"Output directory: {output_dir}")
    print(f"Species: {species}")
    print(f"Database: {databases}")
    print(f"AI backend: {ai_backend}")
    print(f"AI model: {ai_model}")
    print()

    # Create Output Directory
    output_dir.mkdir(exist_ok=True)

    # Step 1: Loading a list of genes
    print("Step 1: Loading a list of genes...")
    with open(input_file, 'r') as f:
        genes = [line.strip() for line in f if line.strip()]
    gene_set = set(genes)
    print(f"It's loaded.{len(gene_set)}Genome.")
    print()

    # Step 2: Configure and load the database
    print("Step 2: Configure and Load Database...")
    
    # Use v1 compatible database path
    db_path = project_root / "database" / "organism" / "v20190612" / "hsa"
    
    if not db_path.exists():
        print(f"Database path does not exist: {db_path}")
        print("Try to find other database locations...")
        
        # Try to find a database
        possible_paths = [
            project_root / "database" / "organism" / "v20190612" / "hsa",
            Path("F:/OneDrive/Documents/TraeSOLO/AllEnricher/AllEnricher-v1/database/organism/v20190612/hsa"),
        ]
        
        for path in possible_paths:
            if path.exists():
                db_path = path
                print(f"Database found: {db_path}")
                break
    else:
        print(f"* Database path: {db_path}")

    try:
        db_manager = DatabaseManager(str(db_path), species)
        print("*DatabaseManager created successfully")
        
        # Loading Database
        print("Loading database...")
        db_manager.load_databases(databases)
        print(f"* Added to database: {databases}")
        
        # Access to data
        background_set = db_manager.get_background_genes()
        database_data = db_manager.get_all_term_data()
        print(f"* Background genes: {len(background_set)}")
        print()
        
    except Exception as e:
        print(f"Could not close temporary folder: %s{e}")
        print("\nContinue testing with Mock mode...")
        database_data = {}
        background_set = set()
        db_manager = None
        print()

    # Step 3: Run enrichment analysis
    print("Step 3: Run Enrichment Analysis...")
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
        print("*EnchmentAnalyzer created successfully")
        
        # Run Analysis
        if database_data:
            results = analyzer.run_analysis(
                gene_set=gene_set,
                background_set=background_set,
                database_data=database_data,
                parallel=False
            )
            print(f"Enrichment analysis completed for {len(results)} databases")
            
            # Save TSV results
            analyzer.save_results(str(output_dir))
            print(f"* Results saved to: {output_dir}")
        else:
            results = {}
            print("No database data, skipping enrichment analysis")
            
        print()
        
    except Exception as e:
        print(f"Enrichment analysis failed: {e}")
        import traceback
        traceback.print_exc()
        results = {}
        print()

    # Step 4: Generate AI Interpretation
    print("Step 4: Generate AI Interpretation...")
    print(f"Use model: {ai_model}")
    
    try:
        # Create DeepSeek interpreter
        interpreter = create_interpreter(
            backend=ai_backend,
            api_key=ai_api_key,
            model=ai_model
        )
        print(f"✓ {type(interpreter.interpreter).__name__}Create Success")
        
        # Generate interpretation
        print("Calling DeepSeek API...")
        interpretations = interpreter.interpret_results(results)
        print(f"* Generates a {len(interpretations)} *")
        
        # Show interpretation results
        for db_name, interpretation in interpretations.items():
            print(f"\n--- {db_name}AI Interpretation---")
            print(interpretation[:500] + "..." if len(interpretation) > 500 else interpretation)
            print()
            
    except Exception as e:
        print(f"AI Interpretation Failed: {e}")
        import traceback
        traceback.print_exc()
        interpretations = {}
        print()

    # Step 5: Generate HTML reports
    print("Step 5: Generate HTML Report...")
    try:
        report_gen = ReportGenerator(str(output_dir), config)
        
        # Generate report
        report_file = output_dir / "enrichment_report.html"
        html = report_gen.generate(
            results=results,
            output_file=str(report_file),
            gene_list=list(gene_set),
            ai_interpretation=interpretations
        )
        
        print(f"* HTML report generated: {report_file}")
        print()
        
    except Exception as e:
        print(f"Could not close temporary folder: %s{e}")
        import traceback
        traceback.print_exc()
        print()

    # Summary
    print("=" * 70)
    print("Test complete!")
    print("=" * 70)
    print()
    print("Output file:")
    print(f"- Results of the analysis of the enrichment: {output_dir}/")
    print(f"- HTML report: {output_dir}/enrichment_report.html")
    print()
    
    if interpretations:
        print("AI Read summary:")
        for db_name, interpretation in interpretations.items():
            print(f"\n  [{db_name}]")
            print(f"  {interpretation[:200]}...")
    print()


if __name__ == "__main__":
    main()
