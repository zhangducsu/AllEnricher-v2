"""
GSVA（基因集变异分析）模块单元测试

测试覆盖范围：
- GSVA 基本功能：模拟表达矩阵 + 基因集，验证输出形状
- 空输入处理：空表达矩阵应返回空 DataFrame
- 基因集大小过滤：min_size / max_size 边界测试
- 三种方法变体：gsva / plage / zscore
- 无交集场景：基因集与表达矩阵无重叠时应返回 0
- 单样本场景：仅包含一个样本的表达矩阵
- 结果值合理性：非 NaN、非 Inf
"""

import pytest
import pandas as pd
import numpy as np
from typing import Set

from allenricher.core.gsva import GSVA


class TestGSVA:
    """GSVA 基因集变异分析测试"""

    @pytest.fixture
    def expression_matrix(self):
        """
        创建测试用模拟表达矩阵

        100 个基因 x 3 个样本，基因名为 GENE_0001 到 GENE_0100，
        表达值从标准正态分布中随机生成。
        """
        np.random.seed(42)
        n_genes = 100
        n_samples = 3
        gene_names = [f"GENE_{i:04d}" for i in range(1, n_genes + 1)]
        sample_names = ["Sample_A", "Sample_B", "Sample_C"]
        data = np.random.randn(n_genes, n_samples)
        return pd.DataFrame(data, index=gene_names, columns=sample_names)

    @pytest.fixture
    def gene_sets(self):
        """
        创建测试用基因集

        两个基因集，每个包含 15 个基因（满足默认 min_size=10）。
        """
        return {
            "Pathway_1": {f"GENE_{i:04d}" for i in range(1, 16)},      # GENE_0001 ~ GENE_0015
            "Pathway_2": {f"GENE_{i:04d}" for i in range(16, 31)},     # GENE_0016 ~ GENE_0030
        }

    def test_gsva_basic(self, expression_matrix, gene_sets):
        """基本功能测试：100 基因 x 3 样本，2 个基因集，输出形状应为 (2, 3)"""
        gsva = GSVA(method="gsva", min_size=10, max_size=500)
        result = gsva.analyze_matrix(expression_matrix, gene_sets)

        # 验证输出形状：2 个通路 x 3 个样本
        assert result.shape == (2, 3), f"期望形状 (2, 3)，实际为 {result.shape}"
        # 验证通路名称正确
        assert list(result.index) == ["Pathway_1", "Pathway_2"]
        # 验证样本名称正确
        assert list(result.columns) == ["Sample_A", "Sample_B", "Sample_C"]

    def test_gsva_empty_input(self):
        """空输入测试：空表达矩阵应返回空 DataFrame"""
        gsva = GSVA(method="gsva")
        empty_matrix = pd.DataFrame()
        gene_sets = {"Pathway_1": {"GENE_0001"}}

        result = gsva.analyze_matrix(empty_matrix, gene_sets)

        assert result.empty, "空表达矩阵应返回空 DataFrame"

    def test_gsva_small_gene_set(self, expression_matrix):
        """小基因集测试：基因集大小 < min_size 应被跳过"""
        # 创建一个只有 5 个基因的基因集（小于默认 min_size=10）
        small_gene_sets = {
            "Small_Pathway": {f"GENE_{i:04d}" for i in range(1, 6)},  # 仅 5 个基因
        }

        gsva = GSVA(method="gsva", min_size=10, max_size=500)
        result = gsva.analyze_matrix(expression_matrix, small_gene_sets)

        assert result.empty, "基因集大小 < min_size 应被跳过，返回空 DataFrame"

    def test_gsva_large_gene_set(self, expression_matrix):
        """大基因集测试：基因集大小 > max_size 应被跳过"""
        # 创建一个包含 600 个基因的基因集（大于默认 max_size=500）
        # 表达矩阵中只有 100 个基因，所以需要设置较小的 max_size
        large_gene_sets = {
            "Large_Pathway": {f"GENE_{i:04d}" for i in range(1, 101)},  # 100 个基因
        }

        gsva = GSVA(method="gsva", min_size=10, max_size=50)  # max_size=50 < 100
        result = gsva.analyze_matrix(expression_matrix, large_gene_sets)

        assert result.empty, "基因集大小 > max_size 应被跳过，返回空 DataFrame"

    def test_gsva_plage_method(self, expression_matrix, gene_sets):
        """PLAGE 方法测试：验证 PLAGE 方法能正常输出"""
        gsva = GSVA(method="plage", min_size=10, max_size=500)
        result = gsva.analyze_matrix(expression_matrix, gene_sets)

        # 验证输出形状
        assert result.shape == (2, 3), f"PLAGE 输出形状期望 (2, 3)，实际为 {result.shape}"
        # 验证结果为有限数值
        assert np.all(np.isfinite(result.values)), "PLAGE 结果应全部为有限数值"

    def test_gsva_zscore_method(self, expression_matrix, gene_sets):
        """Z-score 方法测试：验证 Z-score 方法能正常输出"""
        gsva = GSVA(method="zscore", min_size=10, max_size=500)
        result = gsva.analyze_matrix(expression_matrix, gene_sets)

        # 验证输出形状
        assert result.shape == (2, 3), f"Z-score 输出形状期望 (2, 3)，实际为 {result.shape}"
        # 验证结果为有限数值
        assert np.all(np.isfinite(result.values)), "Z-score 结果应全部为有限数值"

    def test_gsva_no_overlap(self, expression_matrix):
        """无交集测试：基因集与表达矩阵无交集时应返回空 DataFrame"""
        # 创建一个完全不与表达矩阵重叠的基因集
        no_overlap_sets = {
            "No_Overlap_Pathway": {"UNKNOWN_GENE_1", "UNKNOWN_GENE_2", "UNKNOWN_GENE_3"},
        }

        gsva = GSVA(method="gsva", min_size=1, max_size=500)  # 降低 min_size 以排除大小过滤
        result = gsva.analyze_matrix(expression_matrix, no_overlap_sets)

        assert result.empty, "无交集的基因集应被跳过，返回空 DataFrame"

    def test_gsva_single_sample(self):
        """单样本测试：仅包含一个样本的表达矩阵"""
        np.random.seed(42)
        n_genes = 100
        gene_names = [f"GENE_{i:04d}" for i in range(1, n_genes + 1)]
        data = np.random.randn(n_genes, 1)
        single_sample_matrix = pd.DataFrame(data, index=gene_names, columns=["Only_Sample"])

        gene_sets = {
            "Pathway_1": {f"GENE_{i:04d}" for i in range(1, 16)},
        }

        gsva = GSVA(method="gsva", min_size=10, max_size=500)
        result = gsva.analyze_matrix(single_sample_matrix, gene_sets)

        # 验证输出形状：1 个通路 x 1 个样本
        assert result.shape == (1, 1), f"单样本输出形状期望 (1, 1)，实际为 {result.shape}"
        # 验证结果为有限数值
        assert np.all(np.isfinite(result.values)), "单样本结果应全部为有限数值"

    def test_gsva_result_range(self, expression_matrix, gene_sets):
        """结果值合理性测试：结果值应为有限数值（非 NaN、非 Inf）"""
        gsva = GSVA(method="gsva", min_size=10, max_size=500)
        result = gsva.analyze_matrix(expression_matrix, gene_sets)

        # 验证所有值都是有限数值
        assert np.all(np.isfinite(result.values)), "GSVA 结果不应包含 NaN 或 Inf"

        # 验证结果值在合理范围内（GSVA 得分通常在 [-1, 1] 附近，但具体范围取决于数据）
        # 这里仅验证不是极端值
        assert np.all(np.abs(result.values) < 1e6), "GSVA 结果值不应为极端值"


class TestGSVAInit:
    """GSVA 初始化参数验证测试"""

    def test_invalid_method(self):
        """测试无效方法名称应抛出 ValueError"""
        with pytest.raises(ValueError, match="不支持的 GSVA 方法变体"):
            GSVA(method="invalid_method")

    def test_invalid_kcdf(self):
        """测试无效核函数类型应抛出 ValueError"""
        with pytest.raises(ValueError, match="不支持的核函数类型"):
            GSVA(method="gsva", kcdf="Invalid")

    def test_default_parameters(self):
        """测试默认参数值"""
        gsva = GSVA()
        assert gsva.method == "gsva"
        assert gsva.kcdf == "Gaussian"
        assert gsva.tau == 1.0
        assert gsva.min_size == 10
        assert gsva.max_size == 500


class TestGSVACalculatePvalue:
    """GSVA calculate_pvalue 方法测试"""

    def test_returns_nan(self):
        """GSVA 的 calculate_pvalue 应返回 NaN"""
        gsva = GSVA()
        result = gsva.calculate_pvalue(10, 50, 100, 1000)
        assert np.isnan(result), "GSVA 的 calculate_pvalue 应返回 NaN"


class TestGSVACalculateEnrichment:
    """GSVA calculate_enrichment 方法测试（兼容基类接口）"""

    def test_returns_result_for_valid_input(self):
        """有效输入应返回 EnrichmentResult"""
        gsva = GSVA(min_size=2, max_size=100)
        gene_set = {"GENE_1", "GENE_2", "GENE_3"}
        background_set = {"GENE_1", "GENE_2", "GENE_3", "GENE_4", "GENE_5"}
        term_genes = {"GENE_1", "GENE_2", "GENE_3"}

        result = gsva.calculate_enrichment(
            gene_set=gene_set,
            background_set=background_set,
            term_genes=term_genes,
            term_name="Test_Pathway",
            term_id="TEST:001",
            database="TEST"
        )

        assert result is not None
        assert result.term_id == "TEST:001"
        assert result.gene_count == 3

    def test_returns_none_for_small_gene_set(self):
        """基因集太小应返回 None"""
        gsva = GSVA(min_size=10, max_size=500)
        gene_set = {"GENE_1"}
        background_set = {"GENE_1", "GENE_2"}
        term_genes = {"GENE_1"}

        result = gsva.calculate_enrichment(
            gene_set=gene_set,
            background_set=background_set,
            term_genes=term_genes,
            term_name="Small_Pathway",
            term_id="TEST:002",
            database="TEST"
        )

        assert result is None


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
