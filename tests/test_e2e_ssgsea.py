#!/usr/bin/env python3
"""ssGSEA E2E测试 - 单元测试"""

import sys
import time
import json
from pathlib import Path
import pandas as pd
import numpy as np
import unittest

# 添加项目根目录到路径
sys.path.insert(0, str(Path(__file__).parent.parent))
from allenricher.core.enrichment import SSGSEA

TEST_DATA_DIR = Path(__file__).parent.parent / "test_data"
RESULTS_DIR = TEST_DATA_DIR / "e2e_results"


class TestSsGSEAEndToEnd(unittest.TestCase):
    """ssGSEA端到端测试类"""
    
    @classmethod
    def setUpClass(cls):
        """测试类初始化，加载测试数据"""
        cls.expr_matrix = pd.read_csv(
            TEST_DATA_DIR / "expression_matrix_6000.tsv", 
            sep='\t', 
            index_col=0
        )
        
        # 读取测试通路
        cls.gene_sets = {}
        with open(TEST_DATA_DIR / "test_pathways_from_gmt.gmt", 'r') as f:
            for line in f:
                parts = line.strip().split('\t')
                pathway = parts[0]
                genes = set(parts[2:])
                cls.gene_sets[pathway] = genes
        
        # 创建ssGSEA分析器
        cls.ssgsea = SSGSEA(min_size=10, max_size=500)
        
        # 运行分析并记录时间
        cls.start_time = time.time()
        cls.results_df = cls.ssgsea.analyze_matrix(cls.expr_matrix, cls.gene_sets)
        cls.elapsed_time = time.time() - cls.start_time
        
        # 确保结果目录存在
        RESULTS_DIR.mkdir(exist_ok=True)
        
        # 保存结果
        cls.results_df.to_csv(RESULTS_DIR / "ssgsea_results.csv")
    
    def test_01_output_shape(self):
        """测试结果矩阵形状正确"""
        expected_pathways = len(self.gene_sets)
        expected_samples = self.expr_matrix.shape[1]
        
        self.assertEqual(
            self.results_df.shape[0], 
            expected_pathways,
            f"通路数不匹配: 期望 {expected_pathways}, 实际 {self.results_df.shape[0]}"
        )
        self.assertEqual(
            self.results_df.shape[1], 
            expected_samples,
            f"样本数不匹配: 期望 {expected_samples}, 实际 {self.results_df.shape[1]}"
        )
        print(f"✓ 输出矩阵形状正确: {self.results_df.shape}")
    
    def test_02_score_range(self):
        """测试得分在[-1, 1]范围内（ssGSEA的NES范围）"""
        min_score = self.results_df.values.min()
        max_score = self.results_df.values.max()
        
        self.assertGreaterEqual(min_score, -1.0, f"最小得分 {min_score} 小于 -1")
        self.assertLessEqual(max_score, 1.0, f"最大得分 {max_score} 大于 1")
        print(f"✓ 得分范围正确: [{min_score:.3f}, {max_score:.3f}]")
    
    def test_03_no_nan_values(self):
        """测试无NaN值"""
        has_nan = self.results_df.isna().any().any()
        self.assertFalse(has_nan, "结果中包含NaN值")
        print("✓ 无NaN值")
    
    def test_04_no_inf_values(self):
        """测试无Inf值"""
        has_inf = np.isinf(self.results_df.values).any()
        self.assertFalse(has_inf, "结果中包含Inf值")
        print("✓ 无Inf值")
    
    def test_05_execution_time(self):
        """测试执行时间<30秒"""
        max_time = 30.0
        self.assertLess(
            self.elapsed_time, 
            max_time,
            f"执行时间 {self.elapsed_time:.2f}s 超过 {max_time}s"
        )
        print(f"✓ 执行时间: {self.elapsed_time:.2f}s (限制: {max_time}s)")
    
    def test_06_report_json_format(self):
        """测试报告JSON格式正确"""
        # 生成测试报告
        sample_means = self.results_df.mean(axis=0)
        pathway_means = self.results_df.mean(axis=1)
        
        report = {
            "test_name": "ssGSEA Full E2E Test",
            "timestamp": time.strftime("%Y-%m-%d %H:%M:%S"),
            "input_data": {
                "expression_matrix": f"{self.expr_matrix.shape[0]}x{self.expr_matrix.shape[1]}",
                "pathways_tested": len(self.gene_sets)
            },
            "results": {
                "output_shape": list(self.results_df.shape),
                "score_range": [
                    float(self.results_df.values.min()), 
                    float(self.results_df.values.max())
                ],
                "score_mean": float(self.results_df.values.mean()),
                "score_std": float(self.results_df.values.std()),
                "sample_means": sample_means.to_dict(),
                "pathway_means": pathway_means.to_dict(),
                "execution_time": f"{self.elapsed_time:.2f}s"
            },
            "validation": {
                "all_scores_in_range": bool(
                    (self.results_df.values >= -1).all() and 
                    (self.results_df.values <= 1).all()
                ),
                "no_nan_values": bool(~self.results_df.isna().any().any()),
                "no_inf_values": bool(~np.isinf(self.results_df.values).any())
            },
            "status": "passed"
        }
        
        # 保存报告
        report_path = RESULTS_DIR / "ssgsea_report.json"
        with open(report_path, 'w') as f:
            json.dump(report, f, indent=2)
        
        # 验证JSON格式
        self.assertTrue(report_path.exists(), "报告文件未生成")
        
        # 验证JSON可以正确加载
        with open(report_path, 'r') as f:
            loaded_report = json.load(f)
        
        # 验证报告结构
        required_keys = ["test_name", "timestamp", "input_data", "results", "validation", "status"]
        for key in required_keys:
            self.assertIn(key, loaded_report, f"报告缺少必要字段: {key}")
        
        print(f"✓ 报告JSON格式正确，保存于: {report_path}")
    
    def test_07_pathway_names_match(self):
        """测试通路名称匹配"""
        result_pathways = set(self.results_df.index)
        expected_pathways = set(self.gene_sets.keys())
        
        self.assertEqual(
            result_pathways, 
            expected_pathways,
            "结果中的通路名称与输入不匹配"
        )
        print("✓ 通路名称匹配")
    
    def test_08_sample_names_match(self):
        """测试样本名称匹配"""
        result_samples = set(self.results_df.columns)
        expected_samples = set(self.expr_matrix.columns)
        
        self.assertEqual(
            result_samples, 
            expected_samples,
            "结果中的样本名称与输入不匹配"
        )
        print("✓ 样本名称匹配")


if __name__ == "__main__":
    # 创建测试套件
    loader = unittest.TestLoader()
    suite = loader.loadTestsFromTestCase(TestSsGSEAEndToEnd)
    
    # 使用TextTestRunner运行测试
    runner = unittest.TextTestRunner(verbosity=2)
    result = runner.run(suite)
    
    # 输出总结
    print("\n" + "=" * 60)
    print("ssGSEA E2E单元测试总结")
    print("=" * 60)
    print(f"测试总数: {result.testsRun}")
    print(f"通过: {result.testsRun - len(result.failures) - len(result.errors)}")
    print(f"失败: {len(result.failures)}")
    print(f"错误: {len(result.errors)}")
    
    if result.wasSuccessful():
        print("\n✓ 所有测试通过!")
    else:
        print("\n✗ 测试未通过，请查看详细输出")
    
    print("=" * 60)
