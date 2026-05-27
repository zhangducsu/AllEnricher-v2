"""
自定义物种 A (specA) 全量端对端测试

模拟用户完整工作流：
1. 提供自定义注释文件 → build 构建数据库（自动生成 GMT）
2. 提供 200 个差异基因列表 → ORA 分析
3. 提供排序基因列表 → GSEA 分析
4. 提供全基因表达矩阵 → GSVA 分析（3 种方法）
5. 提供全基因表达矩阵 → ssGSEA 分析
"""

import gzip
import json
import os
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from allenricher.database.custom_builder import CustomDatabaseBuilder
from allenricher.core.enrichment import GSEA, SSGSEA, FisherExactTest
from allenricher.core.gsva import GSVA
from allenricher.report.generator import ReportGenerator

TEST_DATA_DIR = os.path.join(os.path.dirname(__file__), "..", "test_data", "custom_species")


def load_test_data():
    """加载测试数据"""
    annot_file = os.path.join(TEST_DATA_DIR, "specA_annotation.tsv")
    deg_file = os.path.join(TEST_DATA_DIR, "specA_de_genes.txt")
    with open(deg_file) as f:
        de_genes = [line.strip() for line in f if line.strip()]
    ranked_file = os.path.join(TEST_DATA_DIR, "specA_ranked_genes.tsv")
    ranked_df = pd.read_csv(ranked_file, sep='\t')
    ranked_genes = ranked_df['gene'].tolist()
    gene_weights = dict(zip(ranked_df['gene'], ranked_df['weight']))
    expr_file = os.path.join(TEST_DATA_DIR, "specA_expression_matrix.tsv")
    expr_matrix = pd.read_csv(expr_file, sep='\t', index_col=0)
    bg_file = os.path.join(TEST_DATA_DIR, "specA_background.txt")
    with open(bg_file) as f:
        background = [line.strip() for line in f if line.strip()]
    meta_file = os.path.join(TEST_DATA_DIR, "test_data_metadata.json")
    with open(meta_file) as f:
        metadata = json.load(f)
    return {
        'annot_file': annot_file,
        'de_genes': de_genes,
        'ranked_genes': ranked_genes,
        'gene_weights': gene_weights,
        'expr_matrix': expr_matrix,
        'background': background,
        'metadata': metadata,
    }


def load_gmt_from_db(db_dir, species, db_name):
    """从构建的数据库目录加载 GMT 文件"""
    gmt_path = os.path.join(db_dir, f"{species}.{db_name}.gmt.gz")
    gene_sets = {}
    with gzip.open(gmt_path, 'rt', encoding='utf-8') as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 3:
                term_id = parts[0]
                genes = set(parts[2:])
                gene_sets[term_id] = genes
    return gene_sets


def load_gene_matrix(db_dir, species, db_name):
    """从构建的数据库加载基因矩阵"""
    matrix_path = os.path.join(db_dir, f"{species}.{db_name}2gene.tab.gz")
    df = pd.read_csv(matrix_path, sep='\t', compression='gzip')
    return df


class TestCustomSpeciesBuild:
    """Step 1: 自定义数据库构建"""

    @pytest.fixture(scope="class")
    def built_db(self):
        """构建自定义数据库，返回数据库目录"""
        data = load_test_data()
        with tempfile.TemporaryDirectory() as tmpdir:
            builder = CustomDatabaseBuilder(root_dir=os.path.join(tmpdir, "database"))
            outdir = builder.build_from_annotation(
                annotation_file=data['annot_file'],
                species="specA",
                taxid=99999,
                db_name="CustomSpecies"
            )
            yield outdir

    def test_build_creates_all_files(self, built_db):
        """构建生成 3 个必需文件"""
        assert os.path.exists(os.path.join(built_db, "specA.CustomSpecies2gene.tab.gz"))
        assert os.path.exists(os.path.join(built_db, "CustomSpecies2disc.gz"))
        assert os.path.exists(os.path.join(built_db, "specA.CustomSpecies.gmt.gz"))

    def test_gmt_term_count(self, built_db):
        """GMT 文件包含 65 个 term（与元数据一致）"""
        gene_sets = load_gmt_from_db(built_db, "specA", "CustomSpecies")
        assert len(gene_sets) == 65

    def test_gmt_gene_coverage(self, built_db):
        """GMT 基因集覆盖大部分基因"""
        gene_sets = load_gmt_from_db(built_db, "specA", "CustomSpecies")
        all_genes_in_gmt = set()
        for genes in gene_sets.values():
            all_genes_in_gmt.update(genes)
        assert len(all_genes_in_gmt) >= 5400

    def test_gene_matrix_shape(self, built_db):
        """基因矩阵维度正确 (6000 genes x 65 terms + Gene column)"""
        df = load_gene_matrix(built_db, "specA", "CustomSpecies")
        assert df.shape[0] >= 5900
        assert df.shape[1] == 66  # Gene + 65 terms

    def test_description_hierarchy(self, built_db):
        """描述文件包含三层级信息"""
        disc_path = os.path.join(built_db, "CustomSpecies2disc.gz")
        with gzip.open(disc_path, 'rt', encoding='utf-8') as f:
            lines = [l.strip() for l in f if l.strip()]
        assert len(lines) == 65
        for line in lines:
            parts = line.split('\t')
            hierarchy = parts[2] if len(parts) > 2 else ""
            levels = hierarchy.split('|')
            assert len(levels) == 3


class TestCustomSpeciesORA:
    """Step 2: ORA 富集分析"""

    @pytest.fixture(scope="class")
    def ora_results(self):
        """执行 ORA 分析"""
        data = load_test_data()
        with tempfile.TemporaryDirectory() as tmpdir:
            builder = CustomDatabaseBuilder(root_dir=os.path.join(tmpdir, "database"))
            db_dir = builder.build_from_annotation(
                annotation_file=data['annot_file'],
                species="specA",
                taxid=99999,
                db_name="CustomSpecies"
            )
            gene_matrix = load_gene_matrix(db_dir, "specA", "CustomSpecies")
            disc_path = os.path.join(db_dir, "CustomSpecies2disc.gz")
            descriptions = {}
            with gzip.open(disc_path, 'rt', encoding='utf-8') as f:
                for line in f:
                    parts = line.strip().split('\t')
                    descriptions[parts[0]] = {
                        'name': parts[1],
                        'description': parts[2] if len(parts) > 2 else parts[1]
                    }
            database_data = {"CustomSpecies": {}}
            for term in gene_matrix.columns:
                if term == "Gene":
                    continue
                term_genes = set(gene_matrix.loc[gene_matrix[term] == 1, "Gene"].tolist())
                database_data["CustomSpecies"][term] = {
                    'genes': term_genes,
                    'name': descriptions.get(term, {}).get('name', term),
                    'description': descriptions.get(term, {}).get('description', term),
                }
            # 使用 FisherExactTest 直接进行 ORA 分析
            method = FisherExactTest()
            gene_set = set(data['de_genes'])
            background_set = set(data['background'])
            
            results = {"CustomSpecies": []}
            for term_id, term_info in database_data["CustomSpecies"].items():
                term_genes = term_info['genes']
                result = method.calculate_enrichment(
                    gene_set=gene_set,
                    background_set=background_set,
                    term_genes=term_genes,
                    term_name=term_info['name'],
                    term_id=term_id,
                    database="CustomSpecies"
                )
                if result:
                    results["CustomSpecies"].append(result)
            
            # 转换为 DataFrame
            if results["CustomSpecies"]:
                df_data = []
                for r in results["CustomSpecies"]:
                    df_data.append({
                        'Term_ID': r.term_id,
                        'Term_Name': r.term_name,
                        'P_value': r.pvalue,
                        'Adjusted_P_value': r.adjusted_pvalue,
                        'Gene_Count': r.gene_count,
                        'Background_Count': r.background_count,
                        'Expected_Count': r.expected_count,
                        'Rich_Factor': r.rich_factor,
                        'Genes': ','.join(r.gene_list),
                    })
                import pandas as pd
                results["CustomSpecies"] = pd.DataFrame(df_data)
            else:
                results["CustomSpecies"] = pd.DataFrame()
            
            yield results

    def test_ora_returns_results(self, ora_results):
        """ORA 返回结果"""
        assert "CustomSpecies" in ora_results
        df = ora_results["CustomSpecies"]
        assert len(df) > 0

    def test_ora_columns(self, ora_results):
        """ORA 结果包含标准列"""
        df = ora_results["CustomSpecies"]
        expected_cols = ["Term_ID", "P_value", "Adjusted_P_value"]
        for col in expected_cols:
            assert col in df.columns or any(col.lower() in c.lower() for c in df.columns)

    def test_ora_significant_terms(self, ora_results):
        """ORA 有显著富集结果"""
        df = ora_results["CustomSpecies"]
        pval_col = None
        for col in df.columns:
            if 'p_value' in col.lower() or 'pvalue' in col.lower() or 'p.value' in col.lower():
                pval_col = col
                break
        if pval_col is not None:
            significant = df[df[pval_col] < 0.05]
            assert len(significant) > 0


class TestCustomSpeciesGSEA:
    """Step 3: GSEA 分析"""

    @pytest.fixture(scope="class")
    def gsea_results(self):
        """执行 GSEA 分析"""
        data = load_test_data()
        with tempfile.TemporaryDirectory() as tmpdir:
            builder = CustomDatabaseBuilder(root_dir=os.path.join(tmpdir, "database"))
            db_dir = builder.build_from_annotation(
                annotation_file=data['annot_file'],
                species="specA",
                taxid=99999,
                db_name="CustomSpecies"
            )
            gene_sets = load_gmt_from_db(db_dir, "specA", "CustomSpecies")
            gene_sets = {k: v for k, v in gene_sets.items() if 10 <= len(v) <= 500}
            gsea = GSEA(permutations=100, min_size=10, max_size=500)
            results = gsea.analyze_matrix(
                expression_matrix=data['expr_matrix'],
                gene_sets=gene_sets
            )
            yield results

    def test_gsea_returns_dataframe(self, gsea_results):
        """GSEA 返回 DataFrame"""
        assert isinstance(gsea_results, pd.DataFrame)
        assert len(gsea_results) > 0

    def test_gsea_shape(self, gsea_results):
        """GSEA 结果行数大于 0"""
        assert gsea_results.shape[0] > 0

    def test_gsea_has_samples(self, gsea_results):
        """GSEA 结果包含样本列"""
        # GSEA.analyze_matrix 返回行=通路, 列=样本的矩阵
        assert gsea_results.shape[1] == 6  # 6 个样本
        assert list(gsea_results.columns) == ['Sample_1', 'Sample_2', 'Sample_3', 'Sample_4', 'Sample_5', 'Sample_6']


class TestCustomSpeciesGSVA:
    """Step 4: GSVA 分析（3 种方法）"""

    @pytest.fixture(scope="class")
    def gsva_data(self):
        """准备 GSVA 数据"""
        data = load_test_data()
        with tempfile.TemporaryDirectory() as tmpdir:
            builder = CustomDatabaseBuilder(root_dir=os.path.join(tmpdir, "database"))
            db_dir = builder.build_from_annotation(
                annotation_file=data['annot_file'],
                species="specA",
                taxid=99999,
                db_name="CustomSpecies"
            )
            gene_sets = load_gmt_from_db(db_dir, "specA", "CustomSpecies")
            gene_sets = {k: v for k, v in gene_sets.items() if 10 <= len(v) <= 500}
            yield {
                'expr_matrix': data['expr_matrix'],
                'gene_sets': gene_sets,
            }

    @pytest.mark.parametrize("method", ["gsva", "plage", "zscore"])
    def test_gsva_method(self, gsva_data, method):
        """GSVA 三种方法均返回正确结果"""
        gsva = GSVA(method=method, min_size=10, max_size=500)
        results = gsva.analyze_matrix(
            expression_matrix=gsva_data['expr_matrix'],
            gene_sets=gsva_data['gene_sets']
        )
        assert isinstance(results, pd.DataFrame)
        assert results.shape[0] > 0
        assert results.shape[1] == 6


class TestCustomSpeciesSsGSEA:
    """Step 5: ssGSEA 分析"""

    @pytest.fixture(scope="class")
    def ssgsea_results(self):
        """执行 ssGSEA 分析"""
        data = load_test_data()
        with tempfile.TemporaryDirectory() as tmpdir:
            builder = CustomDatabaseBuilder(root_dir=os.path.join(tmpdir, "database"))
            db_dir = builder.build_from_annotation(
                annotation_file=data['annot_file'],
                species="specA",
                taxid=99999,
                db_name="CustomSpecies"
            )
            gene_sets = load_gmt_from_db(db_dir, "specA", "CustomSpecies")
            gene_sets = {k: v for k, v in gene_sets.items() if 10 <= len(v) <= 500}
            ssgsea = SSGSEA(min_size=10, max_size=500)
            results = ssgsea.analyze_matrix(
                expression_matrix=data['expr_matrix'],
                gene_sets=gene_sets
            )
            yield results

    def test_ssgsea_returns_dataframe(self, ssgsea_results):
        """ssGSEA 返回 DataFrame"""
        assert isinstance(ssgsea_results, pd.DataFrame)
        assert len(ssgsea_results) > 0

    def test_ssgsea_shape(self, ssgsea_results):
        """ssGSEA 结果维度正确"""
        assert ssgsea_results.shape[0] > 0
        assert ssgsea_results.shape[1] == 6


class TestCustomSpeciesFullWorkflow:
    """Step 6: 完整工作流"""

    def test_full_workflow_no_errors(self):
        """完整工作流无异常"""
        data = load_test_data()
        with tempfile.TemporaryDirectory() as tmpdir:
            # 1. Build
            builder = CustomDatabaseBuilder(root_dir=os.path.join(tmpdir, "database"))
            db_dir = builder.build_from_annotation(
                annotation_file=data['annot_file'],
                species="specA",
                taxid=99999,
                db_name="CustomSpecies"
            )
            assert os.path.exists(db_dir)
            # 2. 加载 GMT
            gene_sets = load_gmt_from_db(db_dir, "specA", "CustomSpecies")
            gene_sets_filtered = {k: v for k, v in gene_sets.items() if 10 <= len(v) <= 500}
            assert len(gene_sets_filtered) > 0
            # 3. ORA
            gene_matrix = load_gene_matrix(db_dir, "specA", "CustomSpecies")
            disc_path = os.path.join(db_dir, "CustomSpecies2disc.gz")
            descriptions = {}
            with gzip.open(disc_path, 'rt', encoding='utf-8') as f:
                for line in f:
                    parts = line.strip().split('\t')
                    descriptions[parts[0]] = {
                        'name': parts[1],
                        'description': parts[2] if len(parts) > 2 else parts[1],
                    }
            database_data = {"CustomSpecies": {}}
            for term in gene_matrix.columns:
                if term == "Gene":
                    continue
                term_genes = set(gene_matrix.loc[gene_matrix[term] == 1, "Gene"].tolist())
                database_data["CustomSpecies"][term] = {
                    'genes': term_genes,
                    'name': descriptions.get(term, {}).get('name', term),
                    'description': descriptions.get(term, {}).get('description', term),
                }
            # 使用 FisherExactTest 直接进行 ORA 分析
            method = FisherExactTest()
            gene_set = set(data['de_genes'])
            background_set = set(data['background'])
            
            ora_results_list = []
            for term_id, term_info in database_data["CustomSpecies"].items():
                term_genes = term_info['genes']
                result = method.calculate_enrichment(
                    gene_set=gene_set,
                    background_set=background_set,
                    term_genes=term_genes,
                    term_name=term_info['name'],
                    term_id=term_id,
                    database="CustomSpecies"
                )
                if result:
                    ora_results_list.append(result)
            
            # 转换为 DataFrame
            if ora_results_list:
                df_data = []
                for r in ora_results_list:
                    df_data.append({
                        'Term_ID': r.term_id,
                        'Term_Name': r.term_name,
                        'P_value': r.pvalue,
                        'Adjusted_P_value': r.adjusted_pvalue,
                        'Gene_Count': r.gene_count,
                        'Background_Count': r.background_count,
                        'Expected_Count': r.expected_count,
                        'Rich_Factor': r.rich_factor,
                        'Genes': ','.join(r.gene_list),
                    })
                ora_results = {"CustomSpecies": pd.DataFrame(df_data)}
            else:
                ora_results = {"CustomSpecies": pd.DataFrame()}
            assert "CustomSpecies" in ora_results
            # 4. GSEA
            gsea_inst = GSEA(permutations=100, min_size=10, max_size=500)
            gsea_results = gsea_inst.analyze_matrix(
                expression_matrix=data['expr_matrix'],
                gene_sets=gene_sets_filtered
            )
            assert isinstance(gsea_results, pd.DataFrame)
            # 5. GSVA (3 methods)
            for method in ["gsva", "plage", "zscore"]:
                gsva = GSVA(method=method, min_size=10, max_size=500)
                gsva_results = gsva.analyze_matrix(
                    expression_matrix=data['expr_matrix'],
                    gene_sets=gene_sets_filtered
                )
                assert isinstance(gsva_results, pd.DataFrame)
            # 6. ssGSEA
            ssgsea = SSGSEA(min_size=10, max_size=500)
            ssgsea_results = ssgsea.analyze_matrix(
                expression_matrix=data['expr_matrix'],
                gene_sets=gene_sets_filtered
            )
            assert isinstance(ssgsea_results, pd.DataFrame)

            # 7. 生成 HTML 报告（含图表）
            report_output_dir = os.path.join(TEST_DATA_DIR, "e2e_results")
            os.makedirs(report_output_dir, exist_ok=True)
            report_file = os.path.join(report_output_dir, "specA_enrichment_report.html")

            report_generator = ReportGenerator(output_dir=report_output_dir)
            html_path = report_generator.generate(
                results=ora_results,
                output_file=report_file,
                gene_list=data['de_genes'],
                gsea_results=gsea_results,
                gsea_gene_sets=gene_sets_filtered,
                gsva_results=gsva_results,
                analysis_method="fisher",
            )
            assert os.path.exists(html_path)
            assert os.path.getsize(html_path) > 1000

            # 8. 生成可视化图表并保存
            import matplotlib
            matplotlib.use('Agg')
            import matplotlib.pyplot as plt
            from allenricher.visualization.gsea_plots import (
                plot_gsea_dotplot, plot_gsea_nes_barplot
            )
            from allenricher.visualization.gsva_plots import plot_pathway_heatmap

            plots_dir = os.path.join(report_output_dir, "plots")
            os.makedirs(plots_dir, exist_ok=True)

            # 构建 GSEA 可视化所需的 DataFrame（逐通路计算 NES）
            gsea_viz_data = []
            for pathway_name, pathway_genes in gene_sets_filtered.items():
                # 取第一个样本的排序基因做可视化
                sample_expr = data['expr_matrix']['Sample_1']
                ranked_genes = sample_expr.sort_values(ascending=False).index.tolist()
                _, nes, pval, _ = gsea_inst.calculate_normalized_es(
                    ranked_genes, pathway_genes
                )
                gsea_viz_data.append({
                    'pathway': pathway_name,
                    'nes': nes,
                    'pvalue': pval,
                    'gene_count': len(pathway_genes & set(ranked_genes)),
                })
            gsea_viz_df = pd.DataFrame(gsea_viz_data)

            # GSEA NES barplot
            gsea_bar_path = os.path.join(plots_dir, "gsea_nes_barplot.png")
            plot_gsea_nes_barplot(
                results_df=gsea_viz_df,
                top_n=15,
                title='GSEA NES Ranking (specA)',
                output_file=gsea_bar_path,
            )
            plt.close('all')
            assert os.path.exists(gsea_bar_path)

            # GSEA dotplot
            gsea_dotplot_path = os.path.join(plots_dir, "gsea_dotplot.png")
            plot_gsea_dotplot(
                results_df=gsea_viz_df,
                top_n=15,
                title='GSEA Enrichment Dotplot (specA)',
                output_file=gsea_dotplot_path,
            )
            plt.close('all')
            assert os.path.exists(gsea_dotplot_path)

            # GSVA heatmap
            gsva_heatmap_path = os.path.join(plots_dir, "gsva_heatmap.png")
            plot_pathway_heatmap(
                scores_df=gsva_results,
                title='GSVA Activity Heatmap (specA)',
                output_file=gsva_heatmap_path,
            )
            plt.close('all')
            assert os.path.exists(gsva_heatmap_path)
