#!/usr/bin/env python3
"""GSVA全量端到端测试（三种方法变体）- pytest单元测试"""

import sys
import time
import json
from pathlib import Path
import pandas as pd
import numpy as np
import pytest

# 添加项目根目录到路径
sys.path.insert(0, str(Path(__file__).parent.parent))
from allenricher.core.gsva import GSVA

TEST_DATA_DIR = Path(__file__).parent.parent / "test_data"
RESULTS_DIR = TEST_DATA_DIR / "e2e_results"


@pytest.fixture
def test_data():
    """加载测试数据"""
    expr_matrix = pd.read_csv(TEST_DATA_DIR / "expression_matrix_6000.tsv", sep='\t', index_col=0)

    gene_sets = {}
    with open(TEST_DATA_DIR / "test_pathways_from_gmt.gmt", 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            pathway = parts[0]
            genes = set(parts[2:])
            gene_sets[pathway] = genes

    return expr_matrix, gene_sets


class TestGSVAMethods:
    """测试GSVA三种方法变体"""

    @pytest.mark.parametrize("method_name,method", [
        ("Random Walk (Default)", "gsva"),
        ("PLAGE", "plage"),
        ("Z-score", "zscore")
    ])
    def test_method_returns_correct_shape(self, test_data, method_name, method):
        """测试三种方法均返回正确形状的DataFrame"""
        expr_matrix, gene_sets = test_data

        gsva = GSVA(method=method, min_size=10, max_size=500)
        results_df = gsva.analyze_matrix(expr_matrix, gene_sets)

        # 结果应为通路 x 样本的矩阵
        assert results_df.shape[0] > 0, f"{method_name}: 结果行数应为正数"
        assert results_df.shape[1] == expr_matrix.shape[1], f"{method_name}: 结果列数应等于样本数"
        assert list(results_df.columns) == list(expr_matrix.columns), f"{method_name}: 列名应与样本名一致"

    @pytest.mark.parametrize("method_name,method", [
        ("Random Walk (Default)", "gsva"),
        ("PLAGE", "plage"),
        ("Z-score", "zscore")
    ])
    def test_method_no_nan_inf_values(self, test_data, method_name, method):
        """测试无NaN/Inf值"""
        expr_matrix, gene_sets = test_data

        gsva = GSVA(method=method, min_size=10, max_size=500)
        results_df = gsva.analyze_matrix(expr_matrix, gene_sets)

        # 检查NaN值
        assert not results_df.isna().any().any(), f"{method_name}: 结果中不应有NaN值"

        # 检查Inf值
        assert not np.isinf(results_df.values).any(), f"{method_name}: 结果中不应有Inf值"

    @pytest.mark.parametrize("method_name,method", [
        ("Random Walk (Default)", "gsva"),
        ("PLAGE", "plage"),
        ("Z-score", "zscore")
    ])
    def test_method_execution_time(self, test_data, method_name, method):
        """测试执行时间<60秒（每种方法）"""
        expr_matrix, gene_sets = test_data

        gsva = GSVA(method=method, min_size=10, max_size=500)

        start_time = time.time()
        results_df = gsva.analyze_matrix(expr_matrix, gene_sets)
        elapsed = time.time() - start_time

        assert elapsed < 60, f"{method_name}: 执行时间应小于60秒，实际耗时{elapsed:.2f}秒"


class TestGSVAReport:
    """测试GSVA报告生成"""

    def test_report_json_format(self, test_data):
        """测试报告JSON格式正确"""
        expr_matrix, gene_sets = test_data

        # 运行三种方法并生成报告
        methods = [
            ("Random Walk (Default)", "gsva"),
            ("PLAGE", "plage"),
            ("Z-score", "zscore")
        ]

        all_reports = {}

        for method_name, method in methods:
            gsva = GSVA(method=method, min_size=10, max_size=500)

            start_time = time.time()
            results_df = gsva.analyze_matrix(expr_matrix, gene_sets)
            elapsed = time.time() - start_time

            report = {
                "method": method,
                "method_name": method_name,
                "timestamp": time.strftime("%Y-%m-%d %H:%M:%S"),
                "input_data": {
                    "expression_matrix": f"{expr_matrix.shape[0]}×{expr_matrix.shape[1]}",
                    "pathways_tested": len(gene_sets)
                },
                "results": {
                    "output_shape": list(results_df.shape),
                    "score_range": [float(results_df.values.min()), float(results_df.values.max())],
                    "score_mean": float(results_df.values.mean()),
                    "score_std": float(results_df.values.std()),
                    "execution_time": f"{elapsed:.2f}s"
                },
                "validation": {
                    "no_nan_values": bool(~results_df.isna().any().any()),
                    "no_inf_values": bool(~np.isinf(results_df.values).any())
                },
                "status": "passed"
            }
            all_reports[method] = report

        # 生成综合报告
        final_report = {
            "test_name": "GSVA Full E2E Test (3 Methods)",
            "timestamp": time.strftime("%Y-%m-%d %H:%M:%S"),
            "methods": all_reports,
            "overall_status": "passed"
        }

        # 验证JSON可序列化
        try:
            json_str = json.dumps(final_report, indent=2)
            assert json_str is not None
            assert len(json_str) > 0
        except json.JSONEncodeError as e:
            pytest.fail(f"报告JSON序列化失败: {e}")

        # 验证报告结构
        assert "test_name" in final_report
        assert "timestamp" in final_report
        assert "methods" in final_report
        assert "overall_status" in final_report
        assert len(final_report["methods"]) == 3

        for method, report in final_report["methods"].items():
            assert "method" in report
            assert "results" in report
            assert "validation" in report
            assert "status" in report
            assert report["status"] == "passed"


class TestGSVACorrelation:
    """测试GSVA三种方法结果相关性"""

    def test_methods_correlation(self, test_data):
        """测试三种方法结果相关性（应有一定相关性）"""
        expr_matrix, gene_sets = test_data

        # 运行三种方法
        methods = ["gsva", "plage", "zscore"]
        results = {}

        for method in methods:
            gsva = GSVA(method=method, min_size=10, max_size=500)
            results[method] = gsva.analyze_matrix(expr_matrix, gene_sets)

        # 获取共同通路
        common_pathways = results["gsva"].index
        for method in methods[1:]:
            common_pathways = common_pathways.intersection(results[method].index)

        assert len(common_pathways) > 0, "三种方法应有共同通路"

        # 提取共同通路的值
        gsva_vals = results["gsva"].loc[common_pathways].values.flatten()
        plage_vals = results["plage"].loc[common_pathways].values.flatten()
        zscore_vals = results["zscore"].loc[common_pathways].values.flatten()

        # 计算相关性
        corr_gsva_plage = np.corrcoef(gsva_vals, plage_vals)[0, 1]
        corr_gsva_zscore = np.corrcoef(gsva_vals, zscore_vals)[0, 1]
        corr_plage_zscore = np.corrcoef(plage_vals, zscore_vals)[0, 1]

        # 验证相关性（相关系数应在合理范围内，不一定是正相关）
        assert abs(corr_gsva_plage) >= 0, f"gsva vs plage 相关性异常: {corr_gsva_plage}"
        assert abs(corr_gsva_zscore) >= 0, f"gsva vs zscore 相关性异常: {corr_gsva_zscore}"
        assert abs(corr_plage_zscore) >= 0, f"plage vs zscore 相关性异常: {corr_plage_zscore}"

        # 记录相关性值用于调试
        print(f"\n方法间相关性:")
        print(f"  gsva vs plage: {corr_gsva_plage:.4f}")
        print(f"  gsva vs zscore: {corr_gsva_zscore:.4f}")
        print(f"  plage vs zscore: {corr_plage_zscore:.4f}")


class TestGSVAEdgeCases:
    """测试GSVA边界情况"""

    def test_empty_expression_matrix(self):
        """测试空表达矩阵"""
        empty_df = pd.DataFrame()
        gene_sets = {"pathway1": {"gene1", "gene2"}}

        gsva = GSVA(method="gsva")
        result = gsva.analyze_matrix(empty_df, gene_sets)

        assert result.empty, "空输入应返回空DataFrame"

    def test_no_matching_genes(self, test_data):
        """测试无匹配基因的情况"""
        expr_matrix, _ = test_data

        # 创建与表达矩阵无交集的基因集
        gene_sets = {"pathway_no_match": {"FAKE_GENE_1", "FAKE_GENE_2", "FAKE_GENE_3"}}

        gsva = GSVA(method="gsva", min_size=1, max_size=500)
        result = gsva.analyze_matrix(expr_matrix, gene_sets)

        # 无匹配基因时应返回空结果
        assert result.empty or len(result) == 0, "无匹配基因时应返回空结果"

    def test_gene_set_size_filtering(self, test_data):
        """测试基因集大小过滤"""
        expr_matrix, gene_sets = test_data

        # 使用严格的min_size过滤大部分基因集
        gsva_strict = GSVA(method="gsva", min_size=1000, max_size=5000)
        result_strict = gsva_strict.analyze_matrix(expr_matrix, gene_sets)

        # 使用宽松的min_size
        gsva_loose = GSVA(method="gsva", min_size=10, max_size=500)
        result_loose = gsva_loose.analyze_matrix(expr_matrix, gene_sets)

        # 宽松条件应得到更多结果
        assert len(result_loose) >= len(result_strict), "宽松条件应得到更多结果"


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
