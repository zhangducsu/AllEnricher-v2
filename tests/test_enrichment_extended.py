"""
AllEnricher v2.0 补充测试

本测试文件补充了 test_enrichment.py 中未覆盖的模块测试：
- GSEA 基因集富集分析算法
- SSGSEA 单样本GSEA算法
- EnrichmentAnalyzer 的 _get_method 方法注册
- Config 的新增验证功能（input_file、background_file 检查）
- generate_term_url URL 生成函数
- SpeciesLookup 物种检索模块（离线模式）
"""

import pytest
import math
import tempfile
from pathlib import Path

from allenricher.core.config import Config, EnrichmentMethod
from allenricher.core.enrichment import (
    EnrichmentAnalyzer,
    GSEA,
    SSGSEA,
    EnrichmentResult,
    generate_term_url,
)


class TestGSEA:
    """GSEA 基因集富集分析算法测试"""

    def test_calculate_enrichment_score_positive(self):
        """测试 GSEA 富集分数计算：基因集中在排序列表顶部富集"""
        gsea = GSEA(permutations=100, min_size=2, max_size=100)

        # 构造测试数据：基因集中在排序列表前部
        ranked_genes = ["GENE_A", "GENE_B", "GENE_C", "GENE_D", "GENE_E",
                        "GENE_F", "GENE_G", "GENE_H", "GENE_I", "GENE_J"]
        gene_set = {"GENE_A", "GENE_B", "GENE_C"}  # 前三个基因

        es, hit_genes = gsea.calculate_enrichment_score(ranked_genes, gene_set)

        assert es > 0  # 正向富集，ES 应为正值
        assert len(hit_genes) > 0  # 应有前沿基因

    def test_calculate_enrichment_score_negative(self):
        """测试 GSEA 富集分数计算：基因集中在排序列表底部富集"""
        gsea = GSEA(permutations=100, min_size=2, max_size=100)

        ranked_genes = ["GENE_A", "GENE_B", "GENE_C", "GENE_D", "GENE_E",
                        "GENE_F", "GENE_G", "GENE_H", "GENE_I", "GENE_J"]
        gene_set = {"GENE_H", "GENE_I", "GENE_J"}  # 最后三个基因

        es, hit_genes = gsea.calculate_enrichment_score(ranked_genes, gene_set)

        # 底部富集时 ES 应为负值或较小的正值
        assert isinstance(es, float)

    def test_calculate_enrichment_score_no_overlap(self):
        """测试无交集时的 ES 计算"""
        gsea = GSEA(permutations=100, min_size=2, max_size=100)

        ranked_genes = ["GENE_A", "GENE_B", "GENE_C"]
        gene_set = {"GENE_X", "GENE_Y"}  # 完全无交集

        es, hit_genes = gsea.calculate_enrichment_score(ranked_genes, gene_set)

        assert es == 0.0  # 无交集时 ES 应为 0
        assert len(hit_genes) == 0

    def test_calculate_enrichment_with_weights(self):
        """测试带权重的 GSEA ES 计算"""
        gsea = GSEA(permutations=100, min_size=2, max_size=100)

        ranked_genes = ["GENE_A", "GENE_B", "GENE_C", "GENE_D", "GENE_E"]
        gene_set = {"GENE_A", "GENE_B", "GENE_C"}
        gene_weights = {"GENE_A": 2.0, "GENE_B": 1.0, "GENE_C": 0.5,
                        "GENE_D": 0.1, "GENE_E": 0.1}

        es, hit_genes = gsea.calculate_enrichment_score(
            ranked_genes, gene_set, gene_weights=gene_weights
        )

        assert es > 0  # 带权重时仍应正向富集

    def test_permutation_test(self):
        """测试 GSEA 置换检验的 p 值计算"""
        gsea = GSEA(permutations=200, min_size=2, max_size=100, seed=42)

        ranked_genes = [f"GENE_{i:03d}" for i in range(50)]
        gene_set = {f"GENE_{i:03d}" for i in range(5)}  # 前5个基因

        # 先计算观察 ES
        observed_es, _ = gsea.calculate_enrichment_score(ranked_genes, gene_set)
        # 再运行置换检验
        pvalue = gsea._run_permutation_test(ranked_genes, gene_set, observed_es)

        assert 0 <= pvalue <= 1  # p 值应在合法范围内
        assert isinstance(pvalue, float)

    def test_calculate_pvalue_placeholder(self):
        """测试 GSEA 的 calculate_pvalue 返回占位值"""
        gsea = GSEA()
        pvalue = gsea.calculate_pvalue(10, 50, 100, 1000)
        assert pvalue == 1.0  # GSEA 不使用此方法，返回占位值


class TestSSGSEA:
    """ssGSEA 单样本 GSEA 算法测试"""

    def test_calculate_enrichment_score(self):
        """测试 ssGSEA 富集分数计算"""
        ssgsea = SSGSEA(min_size=2, max_size=100)

        ranked_genes = ["GENE_A", "GENE_B", "GENE_C", "GENE_D", "GENE_E",
                        "GENE_F", "GENE_G", "GENE_H", "GENE_I", "GENE_J"]
        gene_set = {"GENE_A", "GENE_B", "GENE_C"}

        es, es_min, es_max, hit_genes = ssgsea.calculate_enrichment_score(
            ranked_genes, gene_set
        )

        assert isinstance(es, float)
        assert isinstance(es_min, float)
        assert isinstance(es_max, float)
        assert es >= es_min  # ES 应 >= 最小值
        assert es <= es_max  # ES 应 <= 最大值

    def test_nes_normalization(self):
        """测试 NES 归一化：NES 应在 [-1, 1] 范围内"""
        ssgsea = SSGSEA(min_size=2, max_size=100)

        ranked_genes = [f"GENE_{i:03d}" for i in range(20)]
        gene_set = {f"GENE_{i:03d}" for i in range(5)}

        es, es_min, es_max, _ = ssgsea.calculate_enrichment_score(
            ranked_genes, gene_set
        )

        # NES = ES / (|ES_min| + |ES_max|)
        denominator = abs(es_min) + abs(es_max)
        if denominator > 0:
            nes = es / denominator
            assert -1 <= nes <= 1  # NES 应在 [-1, 1] 范围内

    def test_no_overlap(self):
        """测试无交集时的 ssGSEA"""
        ssgsea = SSGSEA(min_size=2, max_size=100)

        ranked_genes = ["GENE_A", "GENE_B", "GENE_C"]
        gene_set = {"GENE_X", "GENE_Y"}

        es, es_min, es_max, hit_genes = ssgsea.calculate_enrichment_score(
            ranked_genes, gene_set
        )

        assert es == 0.0
        assert es_min == 0.0
        assert es_max == 0.0
        assert len(hit_genes) == 0

    def test_calculate_pvalue_placeholder(self):
        """测试 ssGSEA 的 calculate_pvalue 返回 NaN"""
        ssgsea = SSGSEA()
        pvalue = ssgsea.calculate_pvalue(10, 50, 100, 1000)
        assert math.isnan(pvalue)  # ssGSEA 不使用 p 值


class TestMethodRegistration:
    """测试分析方法注册到 _get_method()"""

    def test_ssgsea_registered(self):
        """测试 ssgsea 方法已注册到 _get_method()"""
        config = Config(method="ssgsea")
        analyzer = EnrichmentAnalyzer(config)

        # 不应抛出 ValueError
        method = analyzer._get_method()
        assert isinstance(method, SSGSEA)

    def test_gsea_registered(self):
        """测试 gsea 方法已注册"""
        config = Config(method="gsea")
        analyzer = EnrichmentAnalyzer(config)

        method = analyzer._get_method()
        assert isinstance(method, GSEA)

    def test_fisher_registered(self):
        """测试 fisher 方法已注册"""
        from allenricher.core.enrichment import FisherExactTest

        config = Config(method="fisher")
        analyzer = EnrichmentAnalyzer(config)

        method = analyzer._get_method()
        assert isinstance(method, FisherExactTest)

    def test_hypergeometric_registered(self):
        """测试 hypergeometric 方法已注册"""
        from allenricher.core.enrichment import HypergeometricTest

        config = Config(method="hypergeometric")
        analyzer = EnrichmentAnalyzer(config)

        method = analyzer._get_method()
        assert isinstance(method, HypergeometricTest)

    def test_invalid_method_raises(self):
        """测试无效方法名在访问 method 属性时抛出 ValueError"""
        config = Config(method="nonexistent")
        # 延迟初始化：构造不再报错，访问 method 属性时才报错
        analyzer = EnrichmentAnalyzer(config)

        with pytest.raises(ValueError, match="Unknown method"):
            _ = analyzer.method


class TestConfigValidation:
    """Config 验证功能补充测试"""

    def test_validate_missing_input_file(self, tmp_path):
        """测试输入文件不存在时的验证"""
        config = Config(input_file="/nonexistent/path/genes.txt")
        errors = config.validate()

        assert any("Input file not found" in e for e in errors)

    def test_validate_existing_input_file(self, tmp_path):
        """测试输入文件存在时验证通过"""
        gene_file = tmp_path / "genes.txt"
        gene_file.write_text("BRCA1\nTP53\n")

        config = Config(input_file=str(gene_file))
        errors = config.validate()

        assert not any("Input file not found" in e for e in errors)

    def test_validate_missing_background_file(self):
        """测试背景基因文件不存在时的验证"""
        config = Config(background_file="/nonexistent/background.txt")
        errors = config.validate()

        assert any("Background file not found" in e for e in errors)

    def test_validate_null_input_file(self):
        """测试 input_file 为 None 时不应报文件不存在错误"""
        config = Config(input_file=None)
        errors = config.validate()

        assert not any("Input file not found" in e for e in errors)

    def test_config_yaml_template_consistency(self):
        """测试 DEFAULT_CONFIG_YAML 中的关键参数与 Config 默认值一致"""
        from allenricher.core.config import DEFAULT_CONFIG_YAML, Config

        config = Config()

        # 检查 YAML 模板中不包含过时的默认值
        assert "max_genes: 500" not in DEFAULT_CONFIG_YAML  # 应为 .inf
        assert "gsea_min_size: 15" not in DEFAULT_CONFIG_YAML  # 应为 10
        assert "n_jobs: 4" not in DEFAULT_CONFIG_YAML  # 应为 1

        # 检查 YAML 模板中包含正确的默认值
        assert "gsea_min_size: 10" in DEFAULT_CONFIG_YAML
        assert "max_genes: .inf" in DEFAULT_CONFIG_YAML


class TestGenerateTermURL:
    """generate_term_url URL 生成函数测试"""

    def test_go_url(self):
        """测试 GO 数据库 URL 生成"""
        url = generate_term_url("GO:0008150", "GO")
        assert "amigo" in url
        assert "GO:0008150" in url

    def test_kegg_url(self):
        """测试 KEGG 数据库 URL 生成"""
        url = generate_term_url("hsa04010", "KEGG")
        assert "kegg.jp" in url
        assert "hsa04010" in url

    def test_reactome_url(self):
        """测试 Reactome 数据库 URL 生成"""
        url = generate_term_url("R-HSA-1234567", "Reactome")
        assert "reactome.org" in url
        assert "R-HSA-1234567" in url

    def test_do_url(self):
        """测试 DO（疾病本体）URL 生成"""
        url = generate_term_url("DOID:1234", "DO")
        assert "disease-ontology" in url

    def test_unknown_database(self):
        """测试未知数据库返回空字符串"""
        url = generate_term_url("UNKNOWN:001", "UnknownDB")
        assert url == ""


class TestSpeciesLookup:
    """SpeciesLookup 物种检索模块测试（离线模式）"""

    @pytest.fixture
    def lookup(self):
        """创建离线模式的 SpeciesLookup 实例"""
        from allenricher.database.species_lookup import SpeciesLookup
        # auto_load=False 避免网络请求，手动加载内置物种数据
        instance = SpeciesLookup(auto_load=False)
        instance._load_builtin_species()  # 仅加载内置物种，不触发网络请求
        instance._loaded = True  # 标记为已加载，防止 lookup 方法触发网络请求
        return instance

    def test_builtin_species_loaded(self, lookup):
        """测试内置物种数据加载"""
        # 内置物种应包含人类
        info = lookup.lookup_by_kegg_code("hsa")
        assert info is not None
        assert info.latin_name == "Homo sapiens"
        assert info.taxonomy_id == 9606

    def test_builtin_species_count(self, lookup):
        """测试内置物种数量（至少包含 16 个预配置物种）"""
        all_species = lookup.get_all_species()

        assert len(all_species) >= 16  # 至少 16 个预配置物种

    def test_lookup_by_latin_name(self, lookup):
        """测试通过拉丁名检索"""
        info = lookup.lookup_by_latin_name("Homo sapiens")

        assert info is not None
        assert info.kegg_code == "hsa"

    def test_lookup_by_taxid(self, lookup):
        """测试通过 taxid 检索"""
        info = lookup.lookup_by_taxid(9606)

        assert info is not None
        assert info.kegg_code == "hsa"

    def test_yeast_taxid_consistency(self, lookup):
        """测试酿酒酵母 taxid 在 species_lookup 和 config 中一致（均为 4932）"""
        from allenricher.core.config import SPECIES_CONFIGS

        info = lookup.lookup_by_kegg_code("sce")

        assert info is not None
        assert info.taxonomy_id == 4932  # 应为正确的 taxid
        assert SPECIES_CONFIGS["sce"].taxonomy_id == 4932  # config 中也应一致

    def test_builtin_taxid_map_consistency(self, lookup):
        """测试 BUILTIN_TAXID_MAP 中的 sce taxid 与 _load_builtin_species 一致（均为 4932）"""
        # 通过 lookup_by_kegg_code 获取的是 _load_builtin_species 的数据
        info_by_code = lookup.lookup_by_kegg_code("sce")
        assert info_by_code is not None
        assert info_by_code.taxonomy_id == 4932

        # 通过 lookup_by_taxid 反向验证 BUILTIN_TAXID_MAP 中的映射
        info_by_taxid = lookup.lookup_by_taxid(4932)
        assert info_by_taxid is not None
        assert info_by_taxid.kegg_code == "sce"


class TestEnrichmentResultURL:
    """EnrichmentResult URL 字段测试"""

    def test_term_url_in_dict(self):
        """测试 EnrichmentResult.to_dict() 包含 Term_URL 字段"""
        result = EnrichmentResult(
            term_id="GO:0008150",
            term_name="biological_process",
            database="GO",
            pvalue=0.001,
            adjusted_pvalue=0.01,
            gene_count=10,
            background_count=100,
            expected_count=5.0,
            rich_factor=2.0,
            gene_list=["BRCA1", "TP53"],
            gene_ratio="10/100",
            background_ratio="100/1000",
            term_url="https://amigo.geneontology.org/amigo/term/GO:0008150"
        )

        d = result.to_dict()
        assert "Term_URL" in d
        assert "amigo" in d["Term_URL"]

    def test_term_url_default_empty(self):
        """测试 Term_URL 默认为空字符串"""
        result = EnrichmentResult(
            term_id="TEST:001",
            term_name="Test",
            database="TEST",
            pvalue=0.05,
            adjusted_pvalue=0.1,
            gene_count=5,
            background_count=50,
            expected_count=2.5,
            rich_factor=2.0,
            gene_list=[],
            gene_ratio="5/100",
            background_ratio="50/1000"
        )

        d = result.to_dict()
        assert d["Term_URL"] == ""


class TestLazyMethodInit:
    """延迟初始化 method 属性测试"""

    def test_constructor_does_not_validate_method(self):
        """测试构造函数不立即验证方法名"""
        # 即使方法名无效，构造也不应报错（延迟初始化）
        config = Config(method="nonexistent_method")
        analyzer = EnrichmentAnalyzer(config)
        assert analyzer.config.method == "nonexistent_method"

    def test_method_lazy_loading(self):
        """测试 method 属性的延迟加载和缓存"""
        config = Config(method="fisher")
        analyzer = EnrichmentAnalyzer(config)

        # 第一次访问：创建方法实例
        method1 = analyzer.method
        # 第二次访问：应返回同一实例（缓存）
        method2 = analyzer.method
        assert method1 is method2

    def test_method_setter(self):
        """测试通过 setter 直接设置方法实例"""
        from allenricher.core.enrichment import FisherExactTest

        config = Config(method="nonexistent")  # 无效方法名
        analyzer = EnrichmentAnalyzer(config)

        # 直接设置有效的方法实例
        analyzer.method = FisherExactTest()
        assert isinstance(analyzer.method, FisherExactTest)


class TestDeadCodeRemoval:
    """死代码移除验证测试"""

    def test_no_cache_fields_in_config(self):
        """测试 Config 中已移除 use_cache 和 cache_dir 字段"""
        config = Config()
        assert not hasattr(config, 'use_cache'), "use_cache 应已被移除"
        assert not hasattr(config, 'cache_dir'), "cache_dir 应已被移除"

    def test_no_cache_in_yaml_template(self):
        """测试 DEFAULT_CONFIG_YAML 中不包含缓存相关配置"""
        from allenricher.core.config import DEFAULT_CONFIG_YAML
        assert "use_cache" not in DEFAULT_CONFIG_YAML
        assert "cache_dir" not in DEFAULT_CONFIG_YAML

    def test_no_utils_package(self):
        """测试 utils 包已被移除"""
        import importlib
        with pytest.raises(ModuleNotFoundError):
            importlib.import_module("allenricher.utils")


class TestHeatmapGeneSeparator:
    """热图基因分隔符修复测试"""

    def test_genes_separated_by_semicolon(self):
        """测试 EnrichmentResult.to_dict() 中基因以分号分隔"""
        result = EnrichmentResult(
            term_id="TEST:001",
            term_name="Test Term",
            database="TEST",
            pvalue=0.001,
            adjusted_pvalue=0.01,
            gene_count=3,
            background_count=50,
            expected_count=2.5,
            rich_factor=2.0,
            gene_list=["BRCA1", "TP53", "EGFR"],
            gene_ratio="3/100",
            background_ratio="50/1000"
        )
        d = result.to_dict()
        # 基因应以分号分隔
        assert d["Genes"] == "BRCA1;TP53;EGFR"
        # 用分号拆分应得到3个基因
        assert len(d["Genes"].split(';')) == 3


class TestAdjustPvaluesNaN:
    """adjust_pvalues 对 NaN p 值的处理测试"""

    def test_nan_pvalues_skip_correction(self):
        """测试包含 NaN p 值时跳过校正"""
        config = Config()
        analyzer = EnrichmentAnalyzer(config)

        # 创建包含 NaN p 值的结果列表（模拟 ssGSEA 输出）
        results = [
            EnrichmentResult(
                term_id=f"TEST:{i:03d}", term_name=f"Term {i}", database="TEST",
                pvalue=float('nan'), adjusted_pvalue=float('nan'),
                gene_count=5, background_count=50,
                expected_count=2.5, rich_factor=2.0,
                gene_list=["G1", "G2"], gene_ratio="5/100", background_ratio="50/1000"
            )
            for i in range(5)
        ]

        # 不应抛出异常，直接返回原始结果
        corrected = analyzer.adjust_pvalues(results, method="BH")
        assert len(corrected) == 5
        # NaN p 值未被修改
        assert all(math.isnan(r.adjusted_pvalue) for r in corrected)

    def test_normal_pvalues_correction(self):
        """测试正常 p 值的校正仍然有效"""
        config = Config()
        analyzer = EnrichmentAnalyzer(config)

        results = [
            EnrichmentResult(
                term_id=f"TEST:{i:03d}", term_name=f"Term {i}", database="TEST",
                pvalue=0.001 * (i + 1), adjusted_pvalue=0.001 * (i + 1),
                gene_count=5, background_count=50,
                expected_count=2.5, rich_factor=2.0,
                gene_list=["G1", "G2"], gene_ratio="5/100", background_ratio="50/1000"
            )
            for i in range(5)
        ]

        corrected = analyzer.adjust_pvalues(results, method="BH")
        # 校正后 p 值应 >= 原始 p 值
        for orig, corr in zip(results, corrected):
            assert corr.adjusted_pvalue >= orig.pvalue


class TestHypergeometricBackgroundTotal:
    """HypergeometricTest background_total 参数测试"""

    def test_uses_background_total(self):
        """测试 HypergeometricTest 正确使用 background_total 参数"""
        from allenricher.core.enrichment import HypergeometricTest

        hypo = HypergeometricTest()
        gene_set = {"A", "B", "C", "D", "E"}
        background_set = {"A", "B", "C", "D", "E", "F", "G", "H", "I", "J"}
        term_genes = {"A", "B", "C"}

        # 不指定 background_total
        result1 = hypo.calculate_enrichment(
            gene_set, background_set, term_genes, "Test", "TEST:001", "TEST"
        )

        # 指定不同的 background_total（模拟 v1.0 语义）
        result2 = hypo.calculate_enrichment(
            gene_set, background_set, term_genes, "Test", "TEST:001", "TEST",
            background_total=20000  # 注释文件中的基因总数
        )

        assert result1 is not None
        assert result2 is not None
        # 不同的 background_total 应产生不同的 p 值
        assert result1.pvalue != result2.pvalue


class TestGSEADefaultMinSize:
    """GSEA 默认 min_size 一致性测试"""

    def test_gsea_default_min_size_is_10(self):
        """测试 GSEA 默认 min_size 为 10（与 Config 一致）"""
        gsea = GSEA()
        assert gsea.min_size == 10

    def test_ssgsea_default_min_size_is_10(self):
        """测试 SSGSEA 默认 min_size 为 10"""
        ssgsea = SSGSEA()
        assert ssgsea.min_size == 10


class TestConfigMaxGenesType:
    """Config max_genes 类型测试"""

    def test_max_genes_is_float_inf(self):
        """测试 max_genes 默认值为 float('inf')"""
        config = Config()
        assert config.max_genes == float('inf')
        assert isinstance(config.max_genes, float)


class TestEnrichmentResultToDictColumnNames:
    """EnrichmentResult.to_dict() 列名与 plotter.py 使用一致性测试"""

    def test_to_dict_column_names_match_expected(self):
        """测试 to_dict() 输出的列名与 plot_all 中使用的列名一致"""
        result = EnrichmentResult(
            term_id="GO:0008150",
            term_name="biological_process",
            database="GO",
            pvalue=0.001,
            adjusted_pvalue=0.01,
            gene_count=10,
            background_count=100,
            expected_count=5.0,
            rich_factor=2.0,
            gene_list=["BRCA1", "TP53", "EGFR"],
            gene_ratio="10/100",
            background_ratio="100/1000"
        )
        d = result.to_dict()

        # 验证 plotter.py 中使用的列名在 to_dict() 输出中都存在
        expected_keys = [
            "Term_ID", "Term_Name", "P_Value", "Adjusted_P_Value",
            "Gene_Count", "Background_Count", "Expected_Count",
            "Rich_Factor", "Gene_Ratio", "Background_Ratio", "Genes"
        ]
        for key in expected_keys:
            assert key in d, f"to_dict() 缺少列名: {key}"

    def test_adjusted_p_value_not_q_value(self):
        """测试校正后 p 值的键名是 Adjusted_P_Value 而非 Q_Value"""
        result = EnrichmentResult(
            term_id="TEST:001", term_name="Test", database="TEST",
            pvalue=0.001, adjusted_pvalue=0.01,
            gene_count=5, background_count=50,
            expected_count=2.5, rich_factor=2.0,
            gene_list=["G1"], gene_ratio="5/100", background_ratio="50/1000"
        )
        d = result.to_dict()
        assert "Adjusted_P_Value" in d
        assert "Q_Value" not in d
        assert d["Adjusted_P_Value"] == 0.01


class TestAPISecurity:
    """API 安全相关测试（路径遍历防护）"""

    def test_path_traversal_prevention(self):
        """测试路径遍历防护逻辑：清理后的路径不应包含 .. 或 /"""
        import re
        malicious_database = "../../etc"
        safe_database = re.sub(r'[^\w\-]', '', malicious_database)
        # "." 和 "/" 被移除，只保留 "etc"（合法的数据库代码格式）
        assert ".." not in safe_database
        assert "/" not in safe_database
        assert safe_database == "etc"

    def test_safe_plot_filename(self):
        """测试安全文件名生成"""
        import re
        database = "GO"
        plot_type = "barplot"
        safe_db = re.sub(r'[^\w\-]', '', database)
        safe_pt = re.sub(r'[^\w\-]', '', plot_type)
        filename = f"{safe_db}_{safe_pt}.pdf"
        assert filename == "GO_barplot.pdf"
        assert "/" not in filename
        assert ".." not in filename


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
