"""
Custom species A (speca) full-end-endpoint test

Simulate complete workflow:
1. Provide a custom annotation file → build Build Database (Auto Generate GMT)
2. Provision 200 A list of the different genes → ORA Analysis
3. Provides a sorted list of genes → GSEA Analysis
4. Provide a full-gene representation matrix → GSVA Analysis (3 Methods)
5. Provide a full-gene representation matrix → ssGSEA Analysis
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
    """Loading test data"""
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
    """Load GMT files from the built database directory"""
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
    """Loading gene matrices from built databases"""
    matrix_path = os.path.join(db_dir, f"{species}.{db_name}2gene.tab.gz")
    df = pd.read_csv(matrix_path, sep='\t', compression='gzip')
    return df


class TestCustomSpeciesBuild:
    """Step 1: Customize database construction"""

    @pytest.fixture(scope="class")
    def built_db(self):
        """Build a custom database and return the database directory"""
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
        """Build generates 3 required files"""
        assert os.path.exists(os.path.join(built_db, "specA.CustomSpecies2gene.tab.gz"))
        assert os.path.exists(os.path.join(built_db, "CustomSpecies2disc.gz"))
        assert os.path.exists(os.path.join(built_db, "specA.CustomSpecies.gmt.gz"))

    def test_gmt_term_count(self, built_db):
        """GMT files contain 65 term (consistent with metadata)"""
        gene_sets = load_gmt_from_db(built_db, "specA", "CustomSpecies")
        assert len(gene_sets) == 65

    def test_gmt_gene_coverage(self, built_db):
        """Verify that the GMT fixture covers most genes in the custom species."""
        gene_sets = load_gmt_from_db(built_db, "specA", "CustomSpecies")
        all_genes_in_gmt = set()
        for genes in gene_sets.values():
            all_genes_in_gmt.update(genes)
        assert len(all_genes_in_gmt) >= 5400

    def test_gene_matrix_shape(self, built_db):
        """Genome Matrix is correct in dimensions (6, 000 Genes x 65 terms + Gene Colomn)"""
        df = load_gene_matrix(built_db, "specA", "CustomSpecies")
        assert df.shape[0] >= 5900
        assert df.shape[1] == 66  # Gene + 65 terms

    def test_description_hierarchy(self, built_db):
        """The description document contains three levels of information"""
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
    """Step 2: ORA enrichment analysis"""

    @pytest.fixture(scope="class")
    def ora_results(self):
        """Execute ORA analysis"""
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
            # Direct OLA analysis with FisherExactTest
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
            
            # Convert to DataFrame
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
        """ORA returns result"""
        assert "CustomSpecies" in ora_results
        df = ora_results["CustomSpecies"]
        assert len(df) > 0

    def test_ora_columns(self, ora_results):
        """OLA results include standard columns"""
        df = ora_results["CustomSpecies"]
        expected_cols = ["Term_ID", "P_value", "Adjusted_P_value"]
        for col in expected_cols:
            assert col in df.columns or any(col.lower() in c.lower() for c in df.columns)

    def test_ora_significant_terms(self, ora_results):
        """ORA has a remarkable abundance of results."""
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
    """Step 3: GSEA analysis"""

    @pytest.fixture(scope="class")
    def gsea_results(self):
        """Implementation GSEA analysis"""
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
        """GSEA returns DataFrame"""
        assert isinstance(gsea_results, pd.DataFrame)
        assert len(gsea_results) > 0

    def test_gsea_shape(self, gsea_results):
        """GSEA Results Lines > 0"""
        assert gsea_results.shape[0] > 0

    def test_gsea_has_samples(self, gsea_results):
        """GSEA results include sample columns"""
        # GSEA.analyze_matrix returns row = route, column = matrix of samples
        assert gsea_results.shape[1] == 6  # 6 samples
        assert list(gsea_results.columns) == ['Sample_1', 'Sample_2', 'Sample_3', 'Sample_4', 'Sample_5', 'Sample_6']


class TestCustomSpeciesGSVA:
    """Step 4: GSVA analysis (3 methods)"""

    @pytest.fixture(scope="class")
    def gsva_data(self):
        """Preparing GSVA data"""
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
        """GSVA returns the correct result in all three methods"""
        gsva = GSVA(method=method, min_size=10, max_size=500)
        results = gsva.analyze_matrix(
            expression_matrix=gsva_data['expr_matrix'],
            gene_sets=gsva_data['gene_sets']
        )
        assert isinstance(results, pd.DataFrame)
        assert results.shape[0] > 0
        assert results.shape[1] == 6


class TestCustomSpeciesSsGSEA:
    """Step 5: ssGSEA analysis"""

    @pytest.fixture(scope="class")
    def ssgsea_results(self):
        """Execute ssGSEA analysis"""
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
        """sGSEA returns DataFrame"""
        assert isinstance(ssgsea_results, pd.DataFrame)
        assert len(ssgsea_results) > 0

    def test_ssgsea_shape(self, ssgsea_results):
        """SGSEA results are in the right dimension"""
        assert ssgsea_results.shape[0] > 0
        assert ssgsea_results.shape[1] == 6


class TestCustomSpeciesFullWorkflow:
    """Step 6: Full workflow"""

    def test_full_workflow_no_errors(self):
        """The whole job is in order."""
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
            # 2. Load GMT
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
            # Direct OLA analysis with FisherExactTest
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
            
            # Convert to DataFrame
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

            # 7. Generate HTML reports (including graphs)
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
                analysis_method="hypergeometric",
            )
            assert os.path.exists(html_path)
            assert os.path.getsize(html_path) > 1000

            # 8. Generate visualized charts and save
            import matplotlib
            matplotlib.use('Agg')
            import matplotlib.pyplot as plt
            from allenricher.visualization.gsea_plots import plot_gsea_lollipop
            from allenricher.visualization.gsva_plots import plot_pathway_heatmap

            plots_dir = os.path.join(report_output_dir, "plots")
            os.makedirs(plots_dir, exist_ok=True)

            # DataFrame for building GSEA visualize (no route NES)
            gsea_viz_data = []
            for pathway_name, pathway_genes in gene_sets_filtered.items():
                # Use the first sample to exercise the ranked-gene visualization path.
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

            # GSEA lollipop
            gsea_lollipop_path = os.path.join(plots_dir, "gsea_lollipop.png")
            plot_gsea_lollipop(
                results_df=gsea_viz_df,
                top_n=15,
                title='GSEA Enrichment Lollipop (specA)',
                output_file=gsea_lollipop_path,
            )
            plt.close('all')
            assert os.path.exists(gsea_lollipop_path)

            # GSVA heatmap
            gsva_heatmap_path = os.path.join(plots_dir, "gsva_heatmap.png")
            plot_pathway_heatmap(
                scores_df=gsva_results,
                title='GSVA Activity Heatmap (specA)',
                output_file=gsva_heatmap_path,
            )
            plt.close('all')
            assert os.path.exists(gsva_heatmap_path)
