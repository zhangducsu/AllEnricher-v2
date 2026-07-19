#!/usr/bin/env python3
"""GSEA End-to-end testing - Unit testing"""

import sys
import json
import time
from pathlib import Path
import pytest
import pandas as pd

# Add Project Root Directory to Path
sys.path.insert(0, str(Path(__file__).parent.parent))
from allenricher.core.enrichment import GSEA


TEST_DATA_DIR = Path(__file__).parent.parent / "test_data"
RESULTS_DIR = TEST_DATA_DIR / "e2e_results"


def load_test_data():
    """Loading test data"""
    # Read the 500 gene ranked list
    ranked_df = pd.read_csv(TEST_DATA_DIR / "ranked_genes.tsv", sep='\t')
    ranked_df = ranked_df.head(500)  # Only 500.
    ranked_genes = ranked_df['gene'].tolist()
    gene_weights = dict(zip(ranked_df['gene'], ranked_df['weight']))
    
    # Read Test Channel
    gene_sets = {}
    with open(TEST_DATA_DIR / "gene_sets.gmt", 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            pathway = parts[0]
            genes = set(parts[2:])
            gene_sets[pathway] = genes
    
    return ranked_genes, gene_weights, gene_sets


class TestGSEAE2E:
    """GSEA End-to-end testing class"""
    
    @pytest.fixture(scope="class")
    def gsea_results(self):
        """Implementation of GSEA analysis and return of results"""
        ranked_genes, gene_weights, gene_sets = load_test_data()
        
        gsea = GSEA(permutations=100)
        results = []
        
        for pathway_name, pathway_genes in gene_sets.items():
            es, nes, pvalue, leading_edge, _ = gsea.calculate_normalized_es(
                ranked_genes, pathway_genes, gene_weights
            )
            results.append({
                'pathway': pathway_name,
                'es': es,
                'nes': nes,
                'pvalue': pvalue,
                'gene_count': len(pathway_genes),
                'leading_edge_count': len(leading_edge)
            })
        
        return pd.DataFrame(results)
    
    def test_results_dataframe_not_empty(self, gsea_results):
        """Test GSEA results DataFrame is not empty"""
        assert not gsea_results.empty, "GSEA results DataFrame are empty"
        assert len(gsea_results) > 0, "No access to GSEA results"
    
    def test_es_range(self, gsea_results):
        """Test ES in range [-1, 1]"""
        es_min = gsea_results['es'].min()
        es_max = gsea_results['es'].max()
        
        assert es_min >= -1.0, f"ES minimum {es_min} less than 1"
        assert es_max <= 1.0, f"ES Max{es_max}Greater than 1"
    
    def test_pvalue_range(self, gsea_results):
        """Test pvalue in range of [0, 1]"""
        pvalue_min = gsea_results['pvalue'].min()
        pvalue_max = gsea_results['pvalue'].max()
        
        assert pvalue_min >= 0.0, f"pvalue minimum{pvalue_min}Less than 0"
        assert pvalue_max <= 1.0, f"pvalue max{pvalue_max}Greater than 1"
    
    def test_execution_time(self):
        """Test Implementation Time<60sec"""
        ranked_genes, gene_weights, gene_sets = load_test_data()
        
        gsea = GSEA(permutations=100)
        
        start_time = time.time()
        
        for pathway_name, pathway_genes in gene_sets.items():
            gsea.calculate_normalized_es(
                ranked_genes, pathway_genes, gene_weights
            )
        
        elapsed = time.time() - start_time
        
        assert elapsed < 60, f"Implementation time{elapsed: .2f}s over 60 seconds"
    
    def test_report_json_format(self, gsea_results, tmp_path):
        """Test report JSON format correctly"""
        # Generate test report
        report = {
            "test_name": "GSEA Full E2E Test",
            "timestamp": time.strftime("%Y-%m-%d %H:%M:%S"),
            "input_data": {
                "ranked_genes": 500,
                "pathways_tested": len(gsea_results),
                "permutations": 100
            },
            "results": {
                "total_pathways": len(gsea_results),
                "significant_p05": int((gsea_results['pvalue'] < 0.05).sum()),
                "significant_p01": int((gsea_results['pvalue'] < 0.01).sum()),
                "positive_enrichment": int((gsea_results['nes'] > 0).sum()),
                "negative_enrichment": int((gsea_results['nes'] < 0).sum()),
                "nes_range": [float(gsea_results['nes'].min()), float(gsea_results['nes'].max())],
                "execution_time": "0.00s"
            },
            "top_results": gsea_results.nsmallest(5, 'pvalue')[['pathway', 'nes', 'pvalue']].to_dict('records'),
            "status": "passed"
        }
        
        # Save and authenticate JSON format
        report_file = tmp_path / "test_report.json"
        with open(report_file, 'w') as f:
            json.dump(report, f, indent=2)
        
        # Read and authenticate JSON
        with open(report_file, 'r') as f:
            loaded_report = json.load(f)
        
        # Authenticate the necessary fields
        assert "test_name" in loaded_report
        assert "timestamp" in loaded_report
        assert "input_data" in loaded_report
        assert "results" in loaded_report
        assert "top_results" in loaded_report
        assert "status" in loaded_report
        
        # Authenticate results fields
        results = loaded_report["results"]
        assert "total_pathways" in results
        assert "significant_p05" in results
        assert "significant_p01" in results
        assert "positive_enrichment" in results
        assert "negative_enrichment" in results
        assert "nes_range" in results
        assert "execution_time" in results
        
        # Validate Data Type
        assert isinstance(results["total_pathways"], int)
        assert isinstance(results["significant_p05"], int)
        assert isinstance(results["nes_range"], list)
        assert len(results["nes_range"]) == 2
    
    def test_nes_calculation(self, gsea_results):
        """Test NES calculation logic"""
        # NES should be the same as ES (or close to zero)
        for _, row in gsea_results.iterrows():
            es = row['es']
            nes = row['nes']
            
            if es > 0:
                assert nes >= 0, f"ES is positive but negative: ES={es}, NES={nes}O"
            elif es < 0:
                assert nes <= 0, f"ES is negative but NES is positive: NES={es}, NES={nes}"
    
    def test_leading_edge_not_empty(self, gsea_results):
        """Test leaving edge is not empty (when ES is not 0)"""
        # For non-zero-ES, there should be a leaving genes.
        for _, row in gsea_results.iterrows():
            if row['es'] != 0:
                assert row['leading_edge_count'] > 0, \
                    f"Pass.{row['pathway']}ES is not 0 but leaving erge is empty"


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
