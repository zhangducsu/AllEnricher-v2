# 自定义物种全量端对端测试 Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** 生成一套模拟 A 物种的自定义注释测试数据（100 term × 三层级 × 6000 基因），派生 200 个差异基因列表、排序基因列表、全基因表达矩阵，完成 ORA + GSEA + GSVA 全套端对端测试。

**Architecture:** 使用 Python 脚本生成测试数据（注释文件 → build → 派生数据），然后编写 pytest 测试文件覆盖完整分析流程：自定义 build → ORA → GSEA → GSVA。所有测试数据和结果保存在 `test_data/custom_species/` 目录。

**Tech Stack:** Python 3.10, pytest, pandas, numpy, allenricher (CustomDatabaseBuilder, EnrichmentAnalyzer, GSEA, SSGSEA, GSVA)

---

## 文件结构

| 文件 | 职责 |
|------|------|
| `test_data/custom_species/generate_test_data.py` | 数据生成脚本（注释文件 + 派生数据） |
| `test_data/custom_species/specA_annotation.tsv` | 生成的四列注释文件（100 term × 三层级 × 6000 基因） |
| `test_data/custom_species/specA_de_genes.txt` | 200 个差异表达基因列表 |
| `test_data/custom_species/specA_ranked_genes.tsv` | 排序基因列表（gene, weight） |
| `test_data/custom_species/specA_expression_matrix.tsv` | 6000 基因 × 6 样本表达矩阵 |
| `test_data/custom_species/specA_background.txt` | 6000 个全基因背景列表 |
| `test_data/custom_species/test_data_metadata.json` | 测试数据元信息 |
| `tests/test_e2e_custom_species.py` | 全量端对端测试文件 |

---

## Task 1: 生成自定义物种测试数据脚本

**Files:**
- Create: `test_data/custom_species/generate_test_data.py`

- [ ] **Step 1: 编写数据生成脚本**

```python
"""
自定义物种 A (specA) 测试数据生成器

生成内容：
1. 四列注释文件: 100 term × 三层级 × ~6000 基因
2. 200 个差异表达基因列表
3. 排序基因列表 (gene, weight)
4. 全基因表达矩阵 (6000 × 6)
5. 背景基因列表
6. 元数据 JSON
"""

import json
import os
import random
import numpy as np
import pandas as pd

random.seed(42)
np.random.seed(42)

OUTPUT_DIR = os.path.dirname(os.path.abspath(__file__))

# ===== 1. 定义三层级分类体系 =====
# Level1: 5 个大类
categories = {
    "Metabolism": {
        "Carbohydrate Metabolism": ["Glycolysis", "TCA Cycle", "Pentose Phosphate", "Gluconeogenesis", "Glycogen Metabolism"],
        "Lipid Metabolism": ["Fatty Acid Oxidation", "Lipid Biosynthesis", "Cholesterol Metabolism", "Phospholipid Metabolism"],
        "Amino Acid Metabolism": ["Arginine Biosynthesis", "Glutamate Metabolism", "Tryptophan Degradation", "Histidine Catabolism"],
    },
    "Genetic Information Processing": {
        "Transcription": ["RNA Polymerase II", "Transcription Factor Activity", "Chromatin Remodeling", "mRNA Processing", "Transcription Regulation"],
        "Translation": ["Ribosome Biogenesis", "tRNA Aminoacylation", "Translation Initiation", "Translation Elongation", "Protein Folding"],
        "DNA Replication": ["DNA Polymerase Activity", "DNA Repair", "DNA Recombination", "Telomere Maintenance"],
    },
    "Environmental Information Processing": {
        "Signal Transduction": ["MAPK Cascade", "PI3K-Akt Signaling", "Wnt Signaling", "Notch Signaling", "JAK-STAT Pathway"],
        "Membrane Transport": ["ABC Transporters", "Ion Channel Activity", "Vesicle Transport", "Endocytosis"],
        "Cell Communication": ["Gap Junction", "Synaptic Transmission", "Hormone Signaling", "Cytokine Signaling"],
    },
    "Cellular Processes": {
        "Cell Growth": ["Cell Cycle Checkpoint", "Mitosis", "Meiosis", "Apoptosis Regulation", "Autophagy"],
        "Cell Motility": ["Actin Cytoskeleton", "Microtubule Organization", "Cell Adhesion", "Cell Migration"],
        "Cell Division": ["Cytokinesis", "Spindle Assembly", "Centrosome Duplication"],
    },
    "Organismal Systems": {
        "Immune System": ["Antigen Processing", "T Cell Activation", "B Cell Receptor", "Complement System", "Inflammatory Response"],
        "Nervous System": ["Neurotransmitter Release", "Neuron Differentiation", "Axon Guidance", "Synaptic Plasticity"],
        "Circulatory System": ["Blood Coagulation", "Vasodilation", "Angiogenesis", "Hemoglobin Biosynthesis"],
    },
}

# 收集所有 term 信息: (term_id, term_name, hierarchy_path)
all_terms = []
term_counter = 0
for cat1, subcats in categories.items():
    for cat2, terms in subcats.items():
        for term_name in terms:
            term_id = f"SPEC{term_counter:04d}"
            hierarchy = f"{cat1}|{cat2}|{term_name}"
            all_terms.append((term_id, term_name, hierarchy))
            term_counter += 1

print(f"共生成 {len(all_terms)} 个 term")

# ===== 2. 生成 6000 个基因符号 =====
gene_pool = [f"GENE{i:05d}" for i in range(1, 6001)]
random.shuffle(gene_pool)

# ===== 3. 分配基因到 term（模拟真实注释分布） =====
# 每个 term 30-120 个基因，总注释行数约 6000-8000 行
annotation_rows = []
gene_to_terms = {}  # gene -> set of term_ids

for term_id, term_name, hierarchy in all_terms:
    # 基因数量: 30-120, 偏态分布
    n_genes = min(max(int(np.random.normal(60, 20)), 30), 120)
    # 优先选择已出现过的基因（模拟真实注释重叠）
    available_genes = list(gene_pool)
    random.shuffle(available_genes)

    assigned = 0
    for gene in available_genes:
        if assigned >= n_genes:
            break
        annotation_rows.append(f"{gene}\t{term_id}\t{term_name}\t{hierarchy}")
        if gene not in gene_to_terms:
            gene_to_terms[gene] = set()
        gene_to_terms[gene].add(term_id)
        assigned += 1

# 确保所有 6000 基因至少出现在一个 term 中
annotated_genes = set(gene_to_terms.keys())
unannotated = [g for g in gene_pool if g not in annotated_genes]
for i, gene in enumerate(unannotated):
    term_idx = i % len(all_terms)
    term_id, term_name, hierarchy = all_terms[term_idx]
    annotation_rows.append(f"{gene}\t{term_id}\t{term_name}\t{hierarchy}")
    if gene not in gene_to_terms:
        gene_to_terms[gene] = set()
    gene_to_terms[gene].add(term_id)

random.shuffle(annotation_rows)

print(f"共生成 {len(annotation_rows)} 行注释")
print(f"注释覆盖基因数: {len(gene_to_terms)}")

# ===== 4. 写入注释文件 =====
annot_file = os.path.join(OUTPUT_DIR, "specA_annotation.tsv")
with open(annot_file, 'w', encoding='utf-8') as f:
    f.write("\n".join(annotation_rows) + "\n")
print(f"注释文件: {annot_file}")

# ===== 5. 生成 200 个差异表达基因 =====
# 选择注释覆盖较多的基因作为 DEG
gene_annotation_count = [(g, len(terms)) for g, terms in gene_to_terms.items()]
gene_annotation_count.sort(key=lambda x: x[1], reverse=True)
de_genes = [g for g, _ in gene_annotation_count[:200]]
random.shuffle(de_genes)

deg_file = os.path.join(OUTPUT_DIR, "specA_de_genes.txt")
with open(deg_file, 'w') as f:
    f.write("\n".join(de_genes) + "\n")
print(f"差异基因列表: {deg_file} ({len(de_genes)} genes)")

# ===== 6. 生成排序基因列表 (gene, weight) =====
# 200 个 DEG 权重较高（正负各半），其余基因权重较低
ranked_data = []
# DEG: 正向富集的基因（前100）
for gene in de_genes[:100]:
    weight = round(np.random.uniform(0.5, 1.0), 4)
    ranked_data.append((gene, weight))
# DEG: 负向富集的基因（后100）
for gene in de_genes[100:]:
    weight = round(np.random.uniform(-1.0, -0.5), 4)
    ranked_data.append((gene, weight))
# 非 DEG: 低权重
non_deg = [g for g in gene_pool if g not in set(de_genes)]
random.shuffle(non_deg)
for gene in non_deg:
    weight = round(np.random.uniform(-0.3, 0.3), 4)
    ranked_data.append((gene, weight))

# 按权重降序排列
ranked_data.sort(key=lambda x: x[1], reverse=True)

ranked_file = os.path.join(OUTPUT_DIR, "specA_ranked_genes.tsv")
with open(ranked_file, 'w') as f:
    f.write("gene\tweight\n")
    for gene, weight in ranked_data:
        f.write(f"{gene}\t{weight}\n")
print(f"排序基因列表: {ranked_file}")

# ===== 7. 生成表达矩阵 (6000 × 6) =====
n_samples = 6
sample_names = [f"Sample_{i+1}" for i in range(n_samples)]

# 基础表达: 正态分布 N(5, 2)
expr_matrix = pd.DataFrame(
    np.random.normal(5, 2, size=(len(gene_pool), n_samples)),
    index=gene_pool,
    columns=sample_names
)

# DEG 在前3个样本中上调，后3个样本中下调
de_set = set(de_genes)
for gene in de_set:
    if gene in expr_matrix.index:
        expr_matrix.loc[gene, sample_names[:3]] += np.random.uniform(2, 5)
        expr_matrix.loc[gene, sample_names[3:]] -= np.random.uniform(1, 3)

# 确保无负值（log2 转换前）
expr_matrix = expr_matrix.clip(lower=0.01)

expr_file = os.path.join(OUTPUT_DIR, "specA_expression_matrix.tsv")
expr_matrix.to_csv(expr_file, sep='\t')
print(f"表达矩阵: {expr_file} ({expr_matrix.shape[0]} genes × {expr_matrix.shape[1]} samples)")

# ===== 8. 生成背景基因列表 =====
bg_file = os.path.join(OUTPUT_DIR, "specA_background.txt")
with open(bg_file, 'w') as f:
    f.write("\n".join(gene_pool) + "\n")
print(f"背景基因列表: {bg_file}")

# ===== 9. 生成元数据 =====
metadata = {
    "timestamp": "2026-05-26",
    "species": "specA",
    "taxid": 99999,
    "db_name": "CustomSpecies",
    "gene_pool_size": len(gene_pool),
    "term_count": len(all_terms),
    "hierarchy_levels": 3,
    "hierarchy_structure": {
        "level1_count": len(categories),
        "level2_count": sum(len(v) for v in categories.values()),
        "level3_count": len(all_terms),
    },
    "test_data": {
        "de_genes": len(de_genes),
        "ranked_genes": len(ranked_data),
        "expression_matrix": f"{expr_matrix.shape[0]}x{expr_matrix.shape[1]}",
        "background_genes": len(gene_pool),
    },
    "annotation_stats": {
        "total_rows": len(annotation_rows),
        "annotated_genes": len(gene_to_terms),
        "avg_genes_per_term": round(len(annotation_rows) / len(all_terms), 1),
    },
    "category_breakdown": {
        cat1: sum(len(terms) for terms in subcats.values())
        for cat1, subcats in categories.items()
    }
}

meta_file = os.path.join(OUTPUT_DIR, "test_data_metadata.json")
with open(meta_file, 'w', encoding='utf-8') as f:
    json.dump(metadata, f, indent=2, ensure_ascii=False)
print(f"元数据: {meta_file}")

print("\n===== 数据生成完成 =====")
```

- [ ] **Step 2: 运行数据生成脚本**

Run: `cd f:\OneDrive\Documents\TraeSOLO\AllEnricher\AllEnricher-v2 && python test_data/custom_species/generate_test_data.py`
Expected: 生成 7 个文件，无报错

- [ ] **Step 3: 验证生成的数据**

检查：
- `specA_annotation.tsv`: 行数 > 6000，列数为 4（tab 分隔）
- `specA_de_genes.txt`: 恰好 200 行
- `specA_ranked_genes.tsv`: 6001 行（含表头）
- `specA_expression_matrix.tsv`: 6001 行 × 7 列（含表头）
- `specA_background.txt`: 6000 行
- `test_data_metadata.json`: 可正确解析

---

## Task 2: 编写全量端对端测试

**Files:**
- Create: `tests/test_e2e_custom_species.py`

- [ ] **Step 1: 编写测试文件**

```python
"""
自定义物种 A (specA) 全量端对端测试

模拟用户完整工作流：
1. 提供自定义注释文件 → build 构建数据库（自动生成 GMT）
2. 提供 200 个差异基因列表 → ORA 分析
3. 提供排序基因列表 → GSEA 分析
4. 提供全基因表达矩阵 → GSVA 分析（3 种方法）
5. 提供全基因表达矩阵 → ssGSEA 分析
6. 生成报告
"""

import gzip
import json
import os
import shutil
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from allenricher.database.custom_builder import CustomDatabaseBuilder
from allenricher.core.enrichment import EnrichmentAnalyzer, GSEA, SSGSEA
from allenricher.core.gsva import GSVA

TEST_DATA_DIR = os.path.join(os.path.dirname(__file__), "..", "test_data", "custom_species")


def load_test_data():
    """加载测试数据"""
    # 注释文件
    annot_file = os.path.join(TEST_DATA_DIR, "specA_annotation.tsv")
    # DEG 列表
    deg_file = os.path.join(TEST_DATA_DIR, "specA_de_genes.txt")
    with open(deg_file) as f:
        de_genes = [line.strip() for line in f if line.strip()]
    # 排序基因列表
    ranked_file = os.path.join(TEST_DATA_DIR, "specA_ranked_genes.tsv")
    ranked_df = pd.read_csv(ranked_file, sep='\t')
    ranked_genes = ranked_df['gene'].tolist()
    gene_weights = dict(zip(ranked_df['gene'], ranked_df['weight']))
    # 表达矩阵
    expr_file = os.path.join(TEST_DATA_DIR, "specA_expression_matrix.tsv")
    expr_matrix = pd.read_csv(expr_file, sep='\t', index_col=0)
    # 背景基因
    bg_file = os.path.join(TEST_DATA_DIR, "specA_background.txt")
    with open(bg_file) as f:
        background = [line.strip() for line in f if line.strip()]
    # 元数据
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
        """GMT 文件包含 100 个 term"""
        gene_sets = load_gmt_from_db(built_db, "specA", "CustomSpecies")
        assert len(gene_sets) == 100

    def test_gmt_gene_coverage(self, built_db):
        """GMT 基因集覆盖大部分基因"""
        gene_sets = load_gmt_from_db(built_db, "specA", "CustomSpecies")
        all_genes_in_gmt = set()
        for genes in gene_sets.values():
            all_genes_in_gmt.update(genes)
        # 至少覆盖 90% 的基因
        assert len(all_genes_in_gmt) >= 5400

    def test_gene_matrix_shape(self, built_db):
        """基因矩阵维度正确"""
        df = load_gene_matrix(built_db, "specA", "CustomSpecies")
        assert df.shape[0] >= 5900  # ~6000 基因
        assert df.shape[1] == 101   # Gene 列 + 100 个 term 列

    def test_description_hierarchy(self, built_db):
        """描述文件包含三层级信息"""
        disc_path = os.path.join(built_db, "CustomSpecies2disc.gz")
        with gzip.open(disc_path, 'rt', encoding='utf-8') as f:
            lines = [l.strip() for l in f if l.strip()]
        assert len(lines) == 100
        # 验证层级格式: Level1|Level2|Level3
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

            # 加载基因矩阵构建数据库数据
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

            # 构建 database_data 格式
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

            analyzer = EnrichmentAnalyzer(method='fisher')
            gene_set = set(data['de_genes'])
            background_set = set(data['background'])
            results = analyzer.run_analysis(
                gene_set=gene_set,
                background_set=background_set,
                database_data=database_data,
                parallel=True
            )
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
            assert col in df.columns or any(col.lower() in c.lower() for c in df.columns), f"缺少列: {col}"

    def test_ora_significant_terms(self, ora_results):
        """ORA 有显著富集结果"""
        df = ora_results["CustomSpecies"]
        # 至少有一些显著结果
        pval_col = None
        for col in df.columns:
            if 'p_value' in col.lower() or 'pvalue' in col.lower() or 'p.value' in col.lower():
                pval_col = col
                break
        if pval_col is not None:
            significant = df[df[pval_col] < 0.05]
            assert len(significant) > 0, "ORA 应至少有一个显著富集的 term"


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

            # 过滤基因集大小
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
        """GSEA 结果行数等于有效基因集数"""
        data = load_test_data()
        # 结果行数应 > 0
        assert gsea_results.shape[0] > 0

    def test_gsea_has_nes(self, gsea_results):
        """GSEA 结果包含 NES 列"""
        assert 'NES' in gsea_results.columns or any('NES' in str(c) for c in gsea_results.columns)


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
        # 列数 = 样本数
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
        assert ssgsea_results.shape[1] == 6  # 6 个样本


class TestCustomSpeciesFullWorkflow:
    """Step 6: 完整工作流（build → ORA → GSEA → GSVA → ssGSEA）"""

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
            analyzer = EnrichmentAnalyzer(method='fisher')
            ora_results = analyzer.run_analysis(
                gene_set=set(data['de_genes']),
                background_set=set(data['background']),
                database_data=database_data,
            )
            assert "CustomSpecies" in ora_results
            assert len(ora_results["CustomSpecies"]) > 0

            # 4. GSEA
            gsea = GSEA(permutations=100, min_size=10, max_size=500)
            gsea_results = gsea.analyze_matrix(
                expression_matrix=data['expr_matrix'],
                gene_sets=gene_sets_filtered
            )
            assert isinstance(gsea_results, pd.DataFrame)
            assert len(gsea_results) > 0

            # 5. GSVA (3 methods)
            for method in ["gsva", "plage", "zscore"]:
                gsva = GSVA(method=method, min_size=10, max_size=500)
                gsva_results = gsva.analyze_matrix(
                    expression_matrix=data['expr_matrix'],
                    gene_sets=gene_sets_filtered
                )
                assert isinstance(gsva_results, pd.DataFrame)
                assert gsva_results.shape[1] == 6

            # 6. ssGSEA
            ssgsea = SSGSEA(min_size=10, max_size=500)
            ssgsea_results = ssgsea.analyze_matrix(
                expression_matrix=data['expr_matrix'],
                gene_sets=gene_sets_filtered
            )
            assert isinstance(ssgsea_results, pd.DataFrame)
            assert ssgsea_results.shape[1] == 6
```

- [ ] **Step 2: 运行端对端测试**

Run: `cd f:\OneDrive\Documents\TraeSOLO\AllEnricher\AllEnricher-v2 && python -m pytest tests/test_e2e_custom_species.py -v --tb=short`
Expected: 所有测试通过

- [ ] **Step 3: 运行全量回归测试**

Run: `cd f:\OneDrive\Documents\TraeSOLO\AllEnricher\AllEnricher-v2 && python -m pytest tests/ -v --tb=short`
Expected: 全部通过，无回归

---

## Task 3: 保存测试结果和总结报告

**Files:**
- Create: `test_data/custom_species/e2e_test_summary.json`

- [ ] **Step 1: 运行测试并保存 JSON 结果**

Run: `cd f:\OneDrive\Documents\TraeSOLO\AllEnricher\AllEnricher-v2 && python -m pytest tests/test_e2e_custom_species.py -v --tb=short --json-report --json-report-file=test_data/custom_species/e2e_test_summary.json 2>&1 || python -m pytest tests/test_e2e_custom_species.py -v --tb=short 2>&1`

如果 `--json-report` 不可用，则手动生成 JSON：

```python
import json, subprocess, datetime
result = subprocess.run(
    ["python", "-m", "pytest", "tests/test_e2e_custom_species.py", "-v", "--tb=short"],
    capture_output=True, text=True, cwd="f:\\OneDrive\\Documents\\TraeSOLO\\AllEnricher\\AllEnricher-v2"
)
summary = {
    "timestamp": datetime.datetime.now().isoformat(),
    "test_file": "tests/test_e2e_custom_species.py",
    "stdout": result.stdout,
    "returncode": result.returncode,
    "passed": result.returncode == 0,
}
with open("test_data/custom_species/e2e_test_summary.json", "w") as f:
    json.dump(summary, f, indent=2, ensure_ascii=False)
```

- [ ] **Step 2: 验证所有产物文件存在**

检查 `test_data/custom_species/` 目录包含：
- `specA_annotation.tsv`
- `specA_de_genes.txt`
- `specA_ranked_genes.tsv`
- `specA_expression_matrix.tsv`
- `specA_background.txt`
- `test_data_metadata.json`
- `generate_test_data.py`
- `e2e_test_summary.json`
