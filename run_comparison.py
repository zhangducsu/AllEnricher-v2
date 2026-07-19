#!/usr/bin/env python3
"""
Run v2 enrichment analysis and compare with v1 results
"""

import sys
import os
from pathlib import Path
import pandas as pd
import numpy as np

# Add Item Path
project_root = Path(__file__).parent
sys.path.insert(0, str(project_root))

from allenricher.core.config import Config
from allenricher.core.enrichment import EnrichmentAnalyzer
from allenricher.database.manager import DatabaseManager


def load_v1_results(v1_results_dir):
    """Loading of the outcome document of v1"""
    results = {}
    v1_dir = Path(v1_results_dir)
    
    db_files = {
        'GO': 'example.glist.GO.xls',
        'KEGG': 'example.glist.KEGG.xls',
        'Reactome': 'example.glist.Reactome.xls',
        'DO': 'example.glist.DO.xls',
        'DisGeNET': 'example.glist.DisGeNET.xls'
    }
    
    for db_name, filename in db_files.items():
        filepath = v1_dir / filename
        if filepath.exists():
            try:
                df = pd.read_csv(filepath, sep='\t')
                # Make sure TermID is string formatting
                if 'TermID' in df.columns:
                    df['TermID'] = df['TermID'].astype(str)
                results[db_name] = df
                print(f"* Loadingv1{db_name}Outcome: {len(df)}Entry")
            except Exception as e:
                print(f"Failed to load {db_name} from v1: {e}")
    
    return results


def main():
    print("=" * 80)
    print("AllEnricher v1 vs v2 Comparative analysis")
    print("=" * 80)
    
    # Path Settings
    v1_dir = Path(__file__).parent.parent / "AllEnricher-v1"
    v1_results_dir = v1_dir / "example" / "allenricher" / "fisher" / "Q0.05"
    v1_db_dir = v1_dir / "database" / "organism" / "v20190612" / "hsa"
    v2_output_dir = Path(__file__).parent / "comparison_output"
    
    # Create Output Directory
    v2_output_dir.mkdir(exist_ok=True)
    
    # 1. Listing of loaded genes
    print("\n1. Loading of gene lists...")
    gene_list_file = Path(__file__).parent / "example_genes.txt"
    with open(gene_list_file, 'r') as f:
        genes = [line.strip() for line in f if line.strip()]
    print(f"Entering number of genes: {len(genes)}")
    
    # 2. Configuration analysis
    print("\n2. Configuration analysis...")
    config = Config(
        species="hsa",
        databases=["GO", "KEGG", "Reactome", "DO"],
        method="fisher",
        correction="BH",
        pvalue_cutoff=1.0,  # No cut, get the full results for comparison.
        qvalue_cutoff=1.0,
        output_dir=str(v2_output_dir),
        generate_report=False
    )
    
    # 3. Loading of databases
    print("\n3. Loading of databases...")
    try:
        db_manager = DatabaseManager(
            database_dir=str(v1_db_dir),
            species="hsa"
        )
        db_manager.load_databases(["GO", "KEGG", "Reactome", "DO"])
        print("* Database loaded successfully")
    except Exception as e:
        print(f"Could not close temporary folder: %s{e}")
        print(f"Try using test_db...")
        test_db_dir = Path(__file__).parent / "test_db"
        if test_db_dir.exists():
            try:
                db_manager = DatabaseManager(
                    database_dir=str(test_db_dir),
                    species="hsa"
                )
                db_manager.load_databases(["GO", "KEGG"])
                print("* Test_db loaded successfully")
            except Exception as e2:
                print(f"test_db failed to load: {e2}")
                return
        else:
            print(f"test_db does not exist: {test_db_dir}")
            return
    
    # 4. Operational v2 analysis
    print("\n4. Run v2 enrichment analysis...")
    try:
        analyzer = EnrichmentAnalyzer(config)
        database_data = db_manager.get_all_term_data()
        
        # Load background genes
        background_genes = db_manager.get_background_genes()
        print(f"Background genes: {len(background_genes)}")
        
        # Run Analysis
        v2_results = analyzer.run_analysis(
            gene_set=set(genes),
            background_set=background_genes,
            database_data=database_data
        )
        
        print(f"• Analysis completed, result database number: {len(v2_results)}")
        
        # Save v2 results
        for db_name, df in v2_results.items():
            output_file = v2_output_dir / f"v2_{db_name}.tsv"
            df.to_csv(output_file, sep='\t', index=False)
            print(f"Save{db_name}: {len(df)}Entry{output_file.name}")
        
    except Exception as e:
        print(f"Analysis failed: {e}")
        import traceback
        traceback.print_exc()
        return
    
    # 5. Comparison of the results with the loading of v1
    print("\n5. Loading v1 results...")
    v1_results = load_v1_results(v1_results_dir)
    
    # 6. Comparative results
    print("\n" + "=" * 80)
    print("Comparison of results")
    print("=" * 80)
    
    comparison_report = []
    common_dbs = set(v1_results.keys()) & set(v2_results.keys())
    
    for db_name in sorted(common_dbs):
        print(f"\n[{db_name}]")
        v1_df = v1_results[db_name]
        v2_df = v2_results[db_name]
        
        # Basic statistics
        print(f"v1 Entry: {len(v1_df)}")
        print(f"v2 Entry: {len(v2_df)}")
        
        # Try to find a common listing to compare
        v1_cols = set(v1_df.columns)
        v2_cols = set(v2_df.columns)
        
        # Find common entry
        term_id_col_v1 = 'TermID' if 'TermID' in v1_cols else None
        term_id_col_v2 = None
        for col in ['TermID', 'term_id', 'ID']:
            if col in v2_cols:
                term_id_col_v2 = col
                break
        
        if term_id_col_v1 and term_id_col_v2:
            v1_terms = set(v1_df[term_id_col_v1].astype(str))
            v2_terms = set(v2_df[term_id_col_v2].astype(str))
            
            common = v1_terms & v2_terms
            v1_only = v1_terms - v2_terms
            v2_only = v2_terms - v1_terms
            
            print(f"Common entry: {len(common)}")
            print(f"v1 Specialized: {len(v1_only)}")
            print(f"v2 Specialized: {len(v2_only)}")
            
            # Identify salient entries (Q<0.05)
            if 'adjP' in v1_df.columns and 'adjusted_pvalue' in v2_df.columns:
                v1_sig = set(v1_df[v1_df['adjP'] < 0.05][term_id_col_v1].astype(str))
                v2_sig = set(v2_df[v2_df['adjusted_pvalue'] < 0.05][term_id_col_v2].astype(str))
                
                sig_common = v1_sig & v2_sig
                sig_v1_only = v1_sig - v2_sig
                sig_v2_only = v2_sig - v1_sig
                
                print(f"\n  Significant Entry (Q<0.05):")
                print(f"    v1: {len(v1_sig)}")
                print(f"    v2: {len(v2_sig)}")
                print(f"Commonly significant: {len(sig_common)}")
                print(f"V1 is very significant: {len(sig_v1_only)}")
                print(f"V2 is remarkable: {len(sig_v2_only)}")
                
                # Save comparison results
                comparison_report.append({
                    'database': db_name,
                    'v1_total': len(v1_df),
                    'v2_total': len(v2_df),
                    'common': len(common),
                    'v1_sig': len(v1_sig),
                    'v2_sig': len(v2_sig),
                    'sig_common': len(sig_common),
                    'sig_v1_only': len(sig_v1_only),
                    'sig_v2_only': len(sig_v2_only)
                })
                
                # Compare previous entries in detail
                if len(sig_common) > 0:
                    print(f"The \n common cross-references (previous five):")
                    
                    # Get data for common entries
                    common_sig_v1 = v1_df[v1_df[term_id_col_v1].astype(str).isin(sig_common)].head()
                    common_sig_v2 = v2_df[v2_df[term_id_col_v2].astype(str).isin(sig_common)].head()
                    
                    for idx, row in common_sig_v1.head(5).iterrows():
                        term_id = str(row[term_id_col_v1])
                        v2_row = v2_df[v2_df[term_id_col_v2].astype(str) == term_id]
                        
                        if not v2_row.empty:
                            v2_row = v2_row.iloc[0]
                            term_name = row.get('TermName', row.get('term_name', row.get('Description', 'N/A')))
                            v1_p = row.get('rawP', row.get('pvalue', 'N/A'))
                            v1_q = row.get('adjP', 'N/A')
                            v2_p = v2_row.get('pvalue', v2_row.get('rawP', 'N/A'))
                            v2_q = v2_row.get('adjusted_pvalue', v2_row.get('adjP', 'N/A'))
                            
                            print(f"    {term_id} ({term_name[:30]}...)")
                            print(f"      v1: p={v1_p:.6g}, q={v1_q:.6g}")
                            print(f"      v2: p={v2_p:.6g}, q={v2_q:.6g}")
    
    # 7. Preservation of comparative reports
    print("\n" + "=" * 80)
    if comparison_report:
        comparison_df = pd.DataFrame(comparison_report)
        report_file = v2_output_dir / "comparison_report.tsv"
        comparison_df.to_csv(report_file, sep='\t', index=False)
        print(f"* Comparative report saved to: {report_file}")
    
    print("\nDone!")
    print(f"Output directory: {v2_output_dir}")


if __name__ == "__main__":
    main()
