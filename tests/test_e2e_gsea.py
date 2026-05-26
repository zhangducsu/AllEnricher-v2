#!/usr/bin/env python3
"""GSEA端到端测试 - 单元测试"""

import sys
import json
import time
from pathlib import Path
import pytest
import pandas as pd

# 添加项目根目录到路径
sys.path.insert(0, str(Path(__file__).parent.parent))
from allenricher.core.enrichment import GSEA


TEST_DATA_DIR = Path(__file__).parent.parent / "test_data"
RESULTS_DIR = TEST_DATA_DIR / "e2e_results"


def load_test_data():
    """加载测试数据"""
    # 读取500基因排序列表
    ranked_df = pd.read_csv(TEST_DATA_DIR / "ranked_genes.tsv", sep='\t')
    ranked_df = ranked_df.head(500)  # 只取前500个
    ranked_genes = ranked_df['gene'].tolist()
    gene_weights = dict(zip(ranked_df['gene'], ranked_df['weight']))
    
    # 读取测试通路
    gene_sets = {}
    with open(TEST_DATA_DIR / "gene_sets.gmt", 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            pathway = parts[0]
            genes = set(parts[2:])
            gene_sets[pathway] = genes
    
    return ranked_genes, gene_weights, gene_sets


class TestGSEAE2E:
    """GSEA端到端测试类"""
    
    @pytest.fixture(scope="class")
    def gsea_results(self):
        """执行GSEA分析并返回结果"""
        ranked_genes, gene_weights, gene_sets = load_test_data()
        
        gsea = GSEA(permutations=100)
        results = []
        
        for pathway_name, pathway_genes in gene_sets.items():
            es, nes, pvalue, leading_edge = gsea.calculate_normalized_es(
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
        """测试GSEA结果DataFrame不为空"""
        assert not gsea_results.empty, "GSEA结果DataFrame为空"
        assert len(gsea_results) > 0, "GSEA结果中没有通路"
    
    def test_es_range(self, gsea_results):
        """测试ES在[-1, 1]范围内"""
        es_min = gsea_results['es'].min()
        es_max = gsea_results['es'].max()
        
        assert es_min >= -1.0, f"ES最小值 {es_min} 小于-1"
        assert es_max <= 1.0, f"ES最大值 {es_max} 大于1"
    
    def test_pvalue_range(self, gsea_results):
        """测试pvalue在[0, 1]范围内"""
        pvalue_min = gsea_results['pvalue'].min()
        pvalue_max = gsea_results['pvalue'].max()
        
        assert pvalue_min >= 0.0, f"pvalue最小值 {pvalue_min} 小于0"
        assert pvalue_max <= 1.0, f"pvalue最大值 {pvalue_max} 大于1"
    
    def test_execution_time(self):
        """测试执行时间<60秒"""
        ranked_genes, gene_weights, gene_sets = load_test_data()
        
        gsea = GSEA(permutations=100)
        
        start_time = time.time()
        
        for pathway_name, pathway_genes in gene_sets.items():
            gsea.calculate_normalized_es(
                ranked_genes, pathway_genes, gene_weights
            )
        
        elapsed = time.time() - start_time
        
        assert elapsed < 60, f"执行时间 {elapsed:.2f}s 超过60秒"
    
    def test_report_json_format(self, gsea_results, tmp_path):
        """测试报告JSON格式正确"""
        # 生成测试报告
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
        
        # 保存并验证JSON格式
        report_file = tmp_path / "test_report.json"
        with open(report_file, 'w') as f:
            json.dump(report, f, indent=2)
        
        # 读取并验证JSON
        with open(report_file, 'r') as f:
            loaded_report = json.load(f)
        
        # 验证必要字段
        assert "test_name" in loaded_report
        assert "timestamp" in loaded_report
        assert "input_data" in loaded_report
        assert "results" in loaded_report
        assert "top_results" in loaded_report
        assert "status" in loaded_report
        
        # 验证results字段
        results = loaded_report["results"]
        assert "total_pathways" in results
        assert "significant_p05" in results
        assert "significant_p01" in results
        assert "positive_enrichment" in results
        assert "negative_enrichment" in results
        assert "nes_range" in results
        assert "execution_time" in results
        
        # 验证数据类型
        assert isinstance(results["total_pathways"], int)
        assert isinstance(results["significant_p05"], int)
        assert isinstance(results["nes_range"], list)
        assert len(results["nes_range"]) == 2
    
    def test_nes_calculation(self, gsea_results):
        """测试NES计算逻辑"""
        # NES应该与ES同号（或接近0）
        for _, row in gsea_results.iterrows():
            es = row['es']
            nes = row['nes']
            
            if es > 0:
                assert nes >= 0, f"ES为正但NES为负: ES={es}, NES={nes}"
            elif es < 0:
                assert nes <= 0, f"ES为负但NES为正: ES={es}, NES={nes}"
    
    def test_leading_edge_not_empty(self, gsea_results):
        """测试leading edge不为空（当ES不为0时）"""
        # 对于非零ES，应该有leading edge基因
        for _, row in gsea_results.iterrows():
            if row['es'] != 0:
                assert row['leading_edge_count'] > 0, \
                    f"通路 {row['pathway']} ES不为0但leading edge为空"


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
