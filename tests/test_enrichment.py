"""
AllEnricher v2.0 单元测试

测试覆盖范围：
- Config 配置类的默认值、验证和物种配置
- FisherExactTest Fisher精确检验的p值计算和富集分析
- HypergeometricTest 超几何检验的p值计算及与Fisher检验的一致性
- EnrichmentResult 结果数据类的字典转换
- EnrichmentAnalyzer 分析引擎的基因列表加载、p值校正和结果过滤
- 枚举类（EnrichmentMethod, CorrectionMethod）的完整性检查
"""

import pytest
import pandas as pd
import numpy as np
from pathlib import Path
import tempfile

from allenricher.core.config import Config, EnrichmentMethod, CorrectionMethod
from allenricher.core.enrichment import (
    EnrichmentAnalyzer,
    FisherExactTest,
    HypergeometricTest,
    EnrichmentResult
)


class TestConfig:
    """配置类（Config）测试"""
    
    def test_default_config(self):
        """测试默认配置值是否正确"""
        config = Config()
        
        assert config.species == "hsa"          # 默认物种：人类
        assert config.method == "hypergeometric"  # 默认方法：超几何检验
        assert config.correction == "BH"        # 默认校正：Benjamini-Hochberg
        assert config.pvalue_cutoff == 0.05     # 默认p值阈值
        assert config.qvalue_cutoff == 0.05     # 默认q值阈值
    
    def test_config_validation(self):
        """测试默认配置的验证（应无错误）"""
        config = Config()
        errors = config.validate()
        
        assert len(errors) == 0  # 默认配置应该全部合法
    
    def test_invalid_method(self):
        """测试无效方法名的检测"""
        config = Config(method="invalid")  # 故意设置无效方法
        errors = config.validate()
        
        assert any("method" in e.lower() for e in errors)  # 应报告方法错误
    
    def test_get_species_config(self):
        """测试物种配置的获取（hsa -> 人类）"""
        config = Config(species="hsa")
        species_config = config.get_species_config()
        
        assert species_config.name == "Homo sapiens"   # 学名
        assert species_config.taxonomy_id == 9606       # 分类学ID


class TestFisherExactTest:
    """Fisher精确检验测试"""
    
    def test_calculate_pvalue(self):
        """测试p值计算：显著富集的情况"""
        fisher = FisherExactTest()
        
        # 测试用例：1000个背景基因中50个属于该条目，100个输入基因中有10个命中
        # 预期：应该显著富集（p < 0.05）
        pvalue = fisher.calculate_pvalue(
            gene_count=10,           # 输入基因中命中该条目的数量
            background_count=50,     # 背景基因中属于该条目的数量
            gene_total=100,          # 输入基因总数
            background_total=1000    # 背景基因总数
        )
        
        assert 0 <= pvalue <= 1          # p值应在[0,1]范围内
        assert pvalue < 0.05             # 该场景应显著富集
    
    def test_calculate_enrichment(self):
        """测试完整的富集分析计算"""
        fisher = FisherExactTest()
        
        # 构造测试数据
        gene_set = {"BRCA1", "TP53", "EGFR", "MYC"}                                    # 输入基因集
        background_set = {"BRCA1", "TP53", "EGFR", "MYC", "KRAS", "PIK3CA", "PTEN", "RB1"}  # 背景基因集
        term_genes = {"BRCA1", "TP53", "EGFR", "KRAS", "PIK3CA"}                       # 条目关联基因
        
        result = fisher.calculate_enrichment(
            gene_set=gene_set,
            background_set=background_set,
            term_genes=term_genes,
            term_name="Test Term",
            term_id="TEST:001",
            database="TEST"
        )
        
        assert result is not None              # 应返回结果
        assert result.term_id == "TEST:001"    # 条目ID正确
        assert result.gene_count == 3          # BRCA1, TP53, EGFR 三个基因命中
        assert result.rich_factor > 1          # 富集因子应大于1（表示富集）


class TestHypergeometricTest:
    """超几何检验测试"""
    
    def test_calculate_pvalue(self):
        """测试超几何检验的p值计算"""
        hyper = HypergeometricTest()
        
        pvalue = hyper.calculate_pvalue(
            gene_count=10,
            background_count=50,
            gene_total=100,
            background_total=1000
        )
        
        assert 0 <= pvalue <= 1  # p值应在合法范围内
    
    def test_comparison_with_fisher(self):
        """比较超几何检验与Fisher精确检验的结果一致性
        
        数学上，超几何检验与Fisher精确检验是等价的，
        因此对于相同的参数，两者的p值应该非常接近。
        """
        fisher = FisherExactTest()
        hyper = HypergeometricTest()
        
        params = {
            "gene_count": 10,
            "background_count": 50,
            "gene_total": 100,
            "background_total": 1000
        }
        
        p_fisher = fisher.calculate_pvalue(**params)
        p_hyper = hyper.calculate_pvalue(**params)
        
        # 两种方法的p值应该非常接近（允许微小浮点误差）
        assert abs(p_fisher - p_hyper) < 0.01


class TestEnrichmentResult:
    """富集分析结果数据类测试"""
    
    def test_to_dict(self):
        """测试结果对象到字典的转换"""
        result = EnrichmentResult(
            term_id="GO:0008150",            # GO条目ID
            term_name="biological_process",  # 条目名称
            database="GO",                    # 数据库来源
            pvalue=0.001,                     # 原始p值
            adjusted_pvalue=0.01,             # 校正后p值
            gene_count=10,                    # 命中基因数
            background_count=100,             # 背景基因数
            expected_count=5.0,               # 期望基因数
            rich_factor=2.0,                  # 富集因子
            gene_list=["BRCA1", "TP53"],      # 命中基因列表
            gene_ratio="10/100",              # 基因比率
            background_ratio="100/1000"       # 背景比率
        )
        
        d = result.to_dict()
        
        assert d["Term_ID"] == "GO:0008150"  # 字典键名应使用大写下划线格式
        assert d["P_Value"] == 0.001
        assert d["Gene_Count"] == 10


class TestEnrichmentAnalyzer:
    """富集分析引擎测试"""
    
    @pytest.fixture
    def config(self):
        """创建测试用配置对象"""
        return Config(
            species="hsa",
            databases=["GO"],
            method="hypergeometric",
            qvalue_cutoff=0.05
        )
    
    @pytest.fixture
    def analyzer(self, config):
        """创建测试用分析器实例"""
        return EnrichmentAnalyzer(config)
    
    def test_load_gene_list(self, analyzer, tmp_path):
        """测试基因列表文件加载功能"""
        # 创建临时测试基因文件（每行一个基因名）
        gene_file = tmp_path / "genes.txt"
        gene_file.write_text("BRCA1\nTP53\nEGFR\n")
        
        genes = analyzer.load_gene_list(str(gene_file))
        
        assert len(genes) == 3        # 应加载3个基因
        assert "BRCA1" in genes       # 基因名应正确加载
    
    def test_adjust_pvalues(self, analyzer):
        """测试多重检验p值校正（BH方法）
        
        BH校正后，校正后p值应 >= 原始p值。
        """
        # 构造5个测试结果，p值递增
        results = [
            EnrichmentResult(
                term_id=f"TERM{i}",
                term_name=f"Term {i}",
                database="TEST",
                pvalue=0.01 * (i + 1),      # p值: 0.01, 0.02, 0.03, 0.04, 0.05
                adjusted_pvalue=0.01 * (i + 1),
                gene_count=5,
                background_count=50,
                expected_count=2.5,
                rich_factor=2.0,
                gene_list=[],
                gene_ratio="5/100",
                background_ratio="50/1000"
            )
            for i in range(5)
        ]
        
        adjusted = analyzer.adjust_pvalues(results, "BH")
        
        assert len(adjusted) == 5
        # 校正后p值应 >= 原始p值（允许微小浮点误差）
        for r in adjusted:
            assert r.adjusted_pvalue >= r.pvalue * 0.99
    
    def test_filter_results(self, analyzer):
        """测试结果过滤功能

        设置 output_all=False，仅保留满足 q 值阈值的显著条目。
        """
        analyzer.config.output_all = False
        results = [
            EnrichmentResult(
                term_id="TERM1",
                term_name="Term 1",
                database="TEST",
                pvalue=0.001,
                adjusted_pvalue=0.01,    # 低于阈值，应保留
                gene_count=5,
                background_count=50,
                expected_count=2.5,
                rich_factor=2.0,
                gene_list=[],
                gene_ratio="5/100",
                background_ratio="50/1000"
            ),
            EnrichmentResult(
                term_id="TERM2",
                term_name="Term 2",
                database="TEST",
                pvalue=0.1,
                adjusted_pvalue=0.5,     # 高于阈值0.05，应被过滤掉
                gene_count=5,
                background_count=50,
                expected_count=2.5,
                rich_factor=2.0,
                gene_list=[],
                gene_ratio="5/100",
                background_ratio="50/1000"
            )
        ]

        filtered = analyzer.filter_results(results)

        assert len(filtered) == 1            # 只保留1个结果
        assert filtered[0].term_id == "TERM1"  # 保留的是TERM1

    def test_filter_results_output_all(self, analyzer):
        """测试默认输出全部结果（与v1一致）

        output_all=True 时仅过滤不满足 min_genes 的条目，保留全部 p 值。
        """
        # 默认 output_all=True
        assert analyzer.config.output_all == True
        results = [
            EnrichmentResult(
                term_id="TERM1",
                term_name="Term 1",
                database="TEST",
                pvalue=0.001,
                adjusted_pvalue=0.01,
                gene_count=5,
                background_count=50,
                expected_count=2.5,
                rich_factor=2.0,
                gene_list=[],
                gene_ratio="5/100",
                background_ratio="50/1000"
            ),
            EnrichmentResult(
                term_id="TERM2",
                term_name="Term 2",
                database="TEST",
                pvalue=0.1,
                adjusted_pvalue=0.5,     # 不显著，但 output_all=True 时也应保留
                gene_count=5,
                background_count=50,
                expected_count=2.5,
                rich_factor=2.0,
                gene_list=[],
                gene_ratio="5/100",
                background_ratio="50/1000"
            ),
            EnrichmentResult(
                term_id="TERM3",
                term_name="Term 3",
                database="TEST",
                pvalue=0.001,
                adjusted_pvalue=0.01,
                gene_count=1,              # 低于 min_genes=2，应被过滤掉
                background_count=50,
                expected_count=2.5,
                rich_factor=2.0,
                gene_list=[],
                gene_ratio="1/100",
                background_ratio="50/1000"
            )
        ]

        filtered = analyzer.filter_results(results)

        assert len(filtered) == 2            # 保留TERM1和TERM2，TERM3因基因数不足被过滤
        assert filtered[0].term_id == "TERM1"
        assert filtered[1].term_id == "TERM2"


class TestEnrichmentMethodEnum:
    """富集分析方法枚举测试"""
    
    def test_all_methods_exist(self):
        """测试所有预期的富集分析方法都已定义"""
        methods = [m.value for m in EnrichmentMethod]
        
        assert "hypergeometric" in methods   # 超几何检验
        assert "gsea" in methods             # GSEA基因集富集分析
        assert "ssgsea" in methods           # ssGSEA单样本GSEA


class TestCorrectionMethodEnum:
    """多重检验校正方法枚举测试"""
    
    def test_all_corrections_exist(self):
        """测试所有预期的校正方法都已定义"""
        corrections = [c.value for c in CorrectionMethod]
        
        assert "BH" in corrections           # Benjamini-Hochberg
        assert "BY" in corrections           # Benjamini-Yekutieli
        assert "bonferroni" in corrections   # Bonferroni校正
        assert "holm" in corrections         # Holm校正
        assert "none" in corrections         # 不校正


if __name__ == "__main__":
    pytest.main([__file__, "-v"])  # 直接运行：python test_enrichment.py
