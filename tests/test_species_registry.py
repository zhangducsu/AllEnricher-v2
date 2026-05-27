import pytest
import tempfile
from pathlib import Path
from allenricher.database.species_registry import SpeciesRegistry, SpeciesEntry


class TestSpeciesEntry:
    """SpeciesEntry 数据类测试"""

    def test_create_minimal_entry(self):
        """测试最小化创建"""
        entry = SpeciesEntry(taxid=9606, latin_name="Homo sapiens")
        assert entry.taxid == 9606
        assert entry.latin_name == "Homo sapiens"
        assert entry.has_go is False
        assert entry.has_kegg is False
        assert entry.go_source is None

    def test_create_full_entry(self):
        """测试完整字段创建"""
        entry = SpeciesEntry(
            taxid=9606, latin_name="Homo sapiens", common_name="Human",
            has_go=True, go_source="ncbi_gene2go", go_gene_count=19500, go_term_count=18500,
            has_kegg=True, kegg_code="hsa", kegg_code_source="kegg", kegg_gene_count=22345, kegg_pathway_count=350,
            has_reactome=True, reactome_code="HSA", reactome_gene_count=10500, reactome_pathway_count=1500,
            has_do=True, do_gene_count=12000, do_term_count=8000
        )
        assert entry.taxid == 9606
        assert entry.has_go is True
        assert entry.kegg_code == "hsa"
        assert entry.kegg_code_source == "kegg"


class TestSpeciesRegistry:
    """SpeciesRegistry 测试"""

    @pytest.fixture
    def temp_registry(self, tmp_path):
        """创建临时注册表"""
        registry = SpeciesRegistry(registry_path=tmp_path / "supported_species.tsv")
        # 添加测试数据
        entries = [
            SpeciesEntry(taxid=9606, latin_name="Homo sapiens", common_name="Human",
                has_go=True, go_source="ncbi_gene2go", go_gene_count=19500, go_term_count=18500,
                has_kegg=True, kegg_code="hsa", kegg_code_source="kegg", kegg_gene_count=22345, kegg_pathway_count=350,
                has_reactome=True, reactome_code="HSA", reactome_gene_count=10500, reactome_pathway_count=1500,
                has_do=True, do_gene_count=12000, do_term_count=8000),
            SpeciesEntry(taxid=10090, latin_name="Mus musculus", common_name="Mouse",
                has_go=True, go_source="ncbi_gene2go", go_gene_count=18200, go_term_count=16800,
                has_kegg=True, kegg_code="mmu", kegg_code_source="kegg", kegg_gene_count=25245, kegg_pathway_count=320,
                has_reactome=True, reactome_code="MMU", reactome_gene_count=9800, reactome_pathway_count=1400,
                has_do=False),
            SpeciesEntry(taxid=3702, latin_name="Arabidopsis thaliana", common_name="Thale cress",
                has_go=True, go_source="uniprot_goa", go_gene_count=15000, go_term_count=12000,
                has_kegg=False,
                has_reactome=False,
                has_do=False),
        ]
        for e in entries:
            registry.add_entry(e)
        registry.save()
        return registry

    def test_save_and_load(self, temp_registry):
        """测试保存和加载"""
        # 重新加载
        loaded = SpeciesRegistry(registry_path=temp_registry.registry_path)
        loaded.load()
        assert len(loaded.entries) == 3
        assert loaded.query_by_taxid(9606) is not None
        assert loaded.query_by_taxid(10090) is not None
        assert loaded.query_by_taxid(3702) is not None

    def test_query_by_taxid(self, temp_registry):
        """测试按TaxID查询"""
        entry = temp_registry.query_by_taxid(9606)
        assert entry is not None
        assert entry.latin_name == "Homo sapiens"
        assert entry.has_go is True
        assert entry.kegg_code == "hsa"

        # 查询不存在的
        assert temp_registry.query_by_taxid(99999) is None

    def test_query_by_latin_name(self, temp_registry):
        """测试按拉丁名查询"""
        results = temp_registry.query_by_latin_name("Homo sapiens")
        assert len(results) == 1
        assert results[0].taxid == 9606

        # 模糊查询（大小写不敏感）
        results = temp_registry.query_by_latin_name("mus")
        assert len(results) == 1
        assert results[0].taxid == 10090

        # 子串查询
        results = temp_registry.query_by_latin_name("Arabidopsis")
        assert len(results) == 1

        # 无匹配
        results = temp_registry.query_by_latin_name("NotExists")
        assert len(results) == 0

    def test_query_by_kegg_code(self, temp_registry):
        """测试按KEGG代码查询"""
        entry = temp_registry.query_by_kegg_code("hsa")
        assert entry is not None
        assert entry.taxid == 9606

        entry = temp_registry.query_by_kegg_code("mmu")
        assert entry is not None
        assert entry.taxid == 10090

        # 不存在
        assert temp_registry.query_by_kegg_code("xxx") is None

    def test_filter_by_databases(self, temp_registry):
        """测试按数据库筛选"""
        # 只筛选GO
        results = temp_registry.filter_by_databases(go=True)
        assert len(results) == 3

        # 只筛选KEGG
        results = temp_registry.filter_by_databases(kegg=True)
        assert len(results) == 2

        # GO + KEGG
        results = temp_registry.filter_by_databases(go=True, kegg=True)
        assert len(results) == 2

        # DO
        results = temp_registry.filter_by_databases(do=True)
        assert len(results) == 1

        # Reactome
        results = temp_registry.filter_by_databases(reactome=True)
        assert len(results) == 2

    def test_get_summary(self, temp_registry):
        """测试统计摘要"""
        summary = temp_registry.get_summary()
        assert summary['total_species'] == 3
        assert summary['go']['count'] == 3
        assert summary['kegg']['count'] == 2
        assert summary['reactome']['count'] == 2
        assert summary['do']['count'] == 1

    def test_get_species_detail(self, temp_registry):
        """测试物种详细信息"""
        detail = temp_registry.get_species_detail(9606)
        assert detail is not None
        assert detail['taxid'] == 9606
        assert detail['latin_name'] == 'Homo sapiens'
        assert detail['has_go'] is True
        assert detail['kegg_code'] == 'hsa'

        # 不存在
        assert temp_registry.get_species_detail(99999) is None

    def test_generate_kegg_abbreviation(self):
        """测试KEGG缩写生成"""
        assert SpeciesRegistry.generate_kegg_abbreviation("Homo sapiens") == "hsa"
        assert SpeciesRegistry.generate_kegg_abbreviation("Mus musculus") == "mmu"
        assert SpeciesRegistry.generate_kegg_abbreviation("Arabidopsis thaliana") == "ath"
        assert SpeciesRegistry.generate_kegg_abbreviation("Rattus norvegicus") == "rno"
        assert SpeciesRegistry.generate_kegg_abbreviation("Danio rerio") == "dre"
        assert SpeciesRegistry.generate_kegg_abbreviation("Saccharomyces cerevisiae") == "sce"
        assert SpeciesRegistry.generate_kegg_abbreviation("Caenorhabditis elegans") == "cel"

    def test_add_entry_overwrite(self, temp_registry):
        """测试添加条目覆盖"""
        new_entry = SpeciesEntry(taxid=9606, latin_name="Homo sapiens Updated",
            has_go=True, go_source="uniprot_goa")
        temp_registry.add_entry(new_entry)
        assert len(temp_registry.entries) == 3  # 不增加
        assert temp_registry.query_by_taxid(9606).latin_name == "Homo sapiens Updated"
