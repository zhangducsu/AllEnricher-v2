"""
GSEA/GSVA/ssGSEA 端对端测试数据生成器

生成人类基因的测试数据：
1. 排序基因列表（带权重）用于 GSEA
2. 多样本表达矩阵用于 ssGSEA/GSVA
"""

import gzip
import numpy as np
import pandas as pd
from pathlib import Path
from typing import Set, Dict

# 真实人类基因符号（常见的差异表达基因）
HUMAN_GENES = [
    # 癌症相关
    "TP53", "BRCA1", "BRCA2", "EGFR", "KRAS", "MYC", "BRAF", "AKT1", "PTEN", "RB1",
    "APC", "CTNNB1", "PIK3CA", "ERBB2", "VEGFA", "CDK4", "CDKN2A", "MDM2", "BCL2", "BAX",
    # 细胞周期
    "CCNA2", "CCNB1", "CCND1", "CCNE1", "CDC25A", "CDC25B", "CHEK1", "CHEK2", "AURKA", "AURKB",
    # DNA 修复
    "ATM", "ATR", "RAD51", "BRCA1", "BRCA2", "PALB2", "CHEK1", "MLH1", "MSH2", "PMS2",
    # 免疫相关
    "PDCD1", "PDL1", "CD274", "CTLA4", "LAG3", "HAVCR2", "TIGIT", "IFNG", "IL6", "TNF",
    # 凋亡
    "CASP3", "CASP8", "CASP9", "BCL2", "BCLXL", "MCL1", "BAX", "BAK1", "PMAIP1", "BBC3",
    # 信号通路
    "WNT1", "WNT2", "CTNNB1", "AXIN1", "AXIN2", "APC", "GSK3B", "CTNNB1", "JUN", "FOS",
    # 代谢
    "HK2", "LDHA", "PKM2", "PGK1", "ENO1", "GAPDH", "TPI1", "PFKL", "G6PD", "TKT",
    # 血管生成
    "VEGFA", "VEGFB", "VEGFC", "FLT1", "KDR", "FLT4", "ANGPT1", "ANGPT2", "TEK", "PDGFA",
    # EMT/转移
    "CDH1", "CDH2", "VIM", "SNAI1", "SNAI2", "ZEB1", "ZEB2", "TWIST1", "MMP2", "MMP9",
    # 干细胞
    "NANOG", "OCT4", "SOX2", "KLF4", "MYC", "LIN28A", "LIN28B", "ALDH1A1", "PROM1", "CD44",
    # 炎症
    "IL1B", "IL6", "IL8", "IL10", "TNF", "CXCL1", "CXCL2", "CXCL8", "CCL2", "CCL5",
    # 激素受体
    "ESR1", "ESR2", "PGR", "AR", "NR3C1", "NR3C2", "THRA", "THRB", "RXRA", "PPARA",
    # 其他常见基因
    "GAPDH", "ACTB", "B2M", "PPIA", "RPLP0", "RPS18", "EEF1A1", "TUBB", "STMN1", "TYMS",
    # 激酶
    "MAPK1", "MAPK3", "MAPK8", "MAPK14", "AKT1", "AKT2", "AKT3", "MTOR", "SRC", "ABL1",
    # 磷酸化
    "STAT1", "STAT3", "STAT5A", "STAT5B", "JAK1", "JAK2", "JAK3", "TYK2", "SOCS1", "SOCS3",
    # 泛素化
    "UBE2I", "UBE2D1", "UBE2D2", "UBE2D3", "UBE2E1", "UBE3A", "VHL", "HIF1A", "EP300", "HDAC1",
    # 核受体
    "PPARG", "RXRA", "VDR", "LXR", "FXR", "SREBF1", "SREBF2", "NR1H3", "NR1H4", "NR5A1",
    # 转移相关
    "MMP1", "MMP3", "MMP7", "MMP9", "MMP14", "TIMP1", "TIMP2", "TIMP3", "MMP2", "PLAUR",
]

# 额外基因填充到 2000 个
EXTRA_GENES = [f"GENE{i:04d}" for i in range(1, 1901)]
ALL_GENES = HUMAN_GENES + EXTRA_GENES

def generate_sorted_gene_list(n_genes: int = 2000, seed: int = 42) -> pd.DataFrame:
    """生成排序基因列表（带权重，如 log2FC）
    
    Args:
        n_genes: 总基因数
        seed: 随机种子
    
    Returns:
        DataFrame (gene, weight, rank)
    """
    np.random.seed(seed)
    
    # 选择基因
    genes = np.random.choice(ALL_GENES, size=n_genes, replace=False)
    
    # 生成权重（模拟 log2 fold change）
    # 上调基因（正权重）和下调基因（负权重），其余为无差异
    direction = np.random.choice([-3, 0, 3], size=n_genes, p=[0.25, 0.50, 0.25])
    weights = np.random.randn(n_genes) * 2 + direction
    
    # 创建 DataFrame
    df = pd.DataFrame({
        'gene': genes,
        'weight': weights
    })
    
    # 按权重降序排序（模拟上调基因在顶部）
    df = df.sort_values('weight', ascending=False).reset_index(drop=True)
    df['rank'] = range(1, len(df) + 1)
    
    return df

def generate_expression_matrix(
    n_genes: int = 2000,
    n_samples: int = 6,
    n_samples_group: int = 3,
    seed: int = 42
) -> pd.DataFrame:
    """生成模拟表达矩阵
    
    Args:
        n_genes: 基因数
        n_samples: 样本总数
        n_samples_group: 每组样本数（用于生成两组：疾病 vs 正常）
        seed: 随机种子
    
    Returns:
        DataFrame (基因 x 样本)，index 为基因名，columns 为样本名
    """
    np.random.seed(seed)
    
    # 选择基因
    genes = np.random.choice(ALL_GENES, size=n_genes, replace=False)
    
    # 生成样本名
    sample_names = [f"Normal_{i+1}" for i in range(n_samples_group)] + \
                   [f"Disease_{i+1}" for i in range(n_samples_group)]
    
    # 生成表达矩阵
    # 正常组：基础表达水平
    normal_expr = np.random.randn(n_genes, n_samples_group) * 2 + 8
    
    # 疾病组：在部分基因上差异表达
    disease_expr = np.random.randn(n_genes, n_samples_group) * 2 + 8
    
    # 选择差异表达基因（上调和下调）
    deg_up = np.random.choice(range(n_genes), size=50, replace=False)
    deg_down = np.random.choice(range(n_genes), size=50, replace=False)
    
    disease_expr[deg_up, :] += np.random.uniform(1, 3, (50, n_samples_group))  # 上调
    disease_expr[deg_down, :] -= np.random.uniform(1, 3, (50, n_samples_group))  # 下调
    
    # 合并
    expr_matrix = np.hstack([normal_expr, disease_expr])
    
    # 创建 DataFrame
    df = pd.DataFrame(expr_matrix, index=genes, columns=sample_names)
    
    # 添加一些 housekeeping 基因（高表达）
    for gene in ["GAPDH", "ACTB", "B2M"]:
        if gene in genes:
            idx = list(genes).index(gene)
            df.iloc[idx, :] = np.random.uniform(10, 14, n_samples)
    
    return df

def generate_gene_sets(n_sets: int = 10, genes_per_set: int = 50) -> dict:
    """生成基因集（模拟通路）
    
    Args:
        n_sets: 基因集数量
        genes_per_set: 每个基因集的基因数
    
    Returns:
        dict {通路名: Set[基因]}
    """
    np.random.seed(42)
    
    gene_sets = {}
    
    # 预定义通路名称
    pathway_names = [
        ("Cell_Cycle", ["CCNA2", "CCNB1", "CCND1", "CCNE1", "CDC25A", "CDC25B", "CHEK1", "AURKA", "TP53", "RB1"]),
        ("DNA_Repair", ["BRCA1", "BRCA2", "ATM", "ATR", "RAD51", "PALB2", "CHEK1", "PMS2", "MLH1", "MDM2"]),
        ("PI3K_AKT", ["AKT1", "AKT2", "PTEN", "PIK3CA", "MTOR", "PDK1", "TSC1", "TSC2", "RHEB", "RPS6KB1"]),
        ("MAPK_Signaling", ["MAPK1", "MAPK3", "EGFR", "KRAS", "BRAF", "RAF1", "MAP2K1", "MAP2K2", "MAPK8", "JUN"]),
        ("Apoptosis", ["BCL2", "BCLXL", "BAX", "BAK1", "CASP3", "CASP8", "CASP9", "PMAIP1", "BBC3", "MCL1"]),
        ("Immune_Response", ["PDCD1", "PDL1", "CTLA4", "IFNG", "IL6", "TNF", "CD274", "LAG3", "HAVCR2", "TIGIT"]),
        ("EMT", ["CDH1", "CDH2", "VIM", "SNAI1", "SNAI2", "ZEB1", "ZEB2", "TWIST1", "MMP2", "MMP9"]),
        ("Angiogenesis", ["VEGFA", "VEGFC", "FLT1", "KDR", "ANGPT1", "ANGPT2", "TEK", "PDGFA", "HIF1A", "EP300"]),
        ("Metabolism", ["HK2", "LDHA", "PKM2", "PGK1", "ENO1", "GAPDH", "G6PD", "TPI1", "PFKL", "PDHA1"]),
        ("Stem_Cell", ["NANOG", "OCT4", "SOX2", "KLF4", "MYC", "LIN28A", "ALDH1A1", "PROM1", "CD44", "NES"]),
    ]
    
    for i in range(n_sets):
        if i < len(pathway_names):
            name, core_genes = pathway_names[i]
        else:
            name = f"Pathway_{i+1}"
            core_genes = []
        
        # 从核心基因 + 随机基因中选取
        all_pool = set(core_genes) | set(np.random.choice(ALL_GENES, size=100, replace=False))
        selected = np.random.choice(list(all_pool), size=min(genes_per_set, len(all_pool)), replace=False)
        gene_sets[f"HSA_{name}_{i+1}"] = set(selected)
    
    return gene_sets

def save_test_data(output_dir: str = "test_data"):
    """保存测试数据到文件"""
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)
    
    # 1. 排序基因列表
    ranked_genes = generate_sorted_gene_list(n_genes=2000, seed=42)
    ranked_genes.to_csv(output_path / "ranked_genes.tsv", sep='\t', index=False)
    print(f"✓ 保存排序基因列表: {output_path / 'ranked_genes.tsv'} ({len(ranked_genes)} genes)")
    
    # 2. 表达矩阵
    expr_matrix = generate_expression_matrix(n_genes=2000, n_samples=6, n_samples_group=3, seed=42)
    expr_matrix.to_csv(output_path / "expression_matrix.tsv", sep='\t')
    print(f"✓ 保存表达矩阵: {output_path / 'expression_matrix.tsv'} ({expr_matrix.shape[0]} genes x {expr_matrix.shape[1]} samples)")
    
    # 3. 基因集
    gene_sets = generate_gene_sets(n_sets=10, genes_per_set=50)
    
    # 保存为文本格式（每行一个基因集）
    with open(output_path / "gene_sets.gmt", 'w') as f:
        for pathway, genes in gene_sets.items():
            genes_str = '\t'.join(genes)
            f.write(f"{pathway}\tpathway\t{genes_str}\n")
    print(f"✓ 保存基因集: {output_path / 'gene_sets.gmt'} ({len(gene_sets)} sets)")
    
    # 4. 保存差异表达基因列表（用于验证）
    deg_genes = ranked_genes.head(100)['gene'].tolist()
    with open(output_path / "top_degs.txt", 'w') as f:
        f.write('\n'.join(deg_genes))
    print(f"✓ 保存差异基因: {output_path / 'top_degs.txt'} ({len(deg_genes)} genes)")
    
    return ranked_genes, expr_matrix, gene_sets

def extract_gene_pool_from_gmt(gmt_dir: str) -> Set[str]:
    """从所有GMT文件提取基因池"""
    gene_pool = set()
    gmt_files = [
        "hsa.GO.gmt.gz",
        "hsa.KEGG.gmt.gz",
        "hsa.Reactome.gmt.gz",
        "hsa.DO.gmt.gz"
    ]

    for gmt_file in gmt_files:
        gmt_path = Path(gmt_dir) / gmt_file
        if gmt_path.exists():
            with gzip.open(gmt_path, 'rt', encoding='utf-8') as f:
                for line in f:
                    parts = line.strip().split('\t')
                    if len(parts) >= 3:
                        genes = parts[2:]
                        gene_pool.update(genes)

    return gene_pool

def select_test_pathways_from_gmt(
    gmt_dir: str,
    pathways_per_db: int = 5,
    min_genes: int = 10,
    max_genes: int = 200
) -> Dict[str, Set[str]]:
    """从各数据库GMT中选取适合测试的通路"""
    test_pathways = {}

    db_files = {
        "GO": "hsa.GO.gmt.gz",
        "KEGG": "hsa.KEGG.gmt.gz",
        "Reactome": "hsa.Reactome.gmt.gz",
        "DO": "hsa.DO.gmt.gz"
    }

    for db_name, gmt_file in db_files.items():
        gmt_path = Path(gmt_dir) / gmt_file
        if not gmt_path.exists():
            continue

        candidates = []
        with gzip.open(gmt_path, 'rt', encoding='utf-8') as f:
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) >= 3:
                    pathway_name = parts[0]
                    genes = set(parts[2:])
                    gene_count = len(genes)

                    if min_genes <= gene_count <= max_genes:
                        candidates.append((pathway_name, genes, gene_count))

        # 按基因数排序，选取适中的
        candidates.sort(key=lambda x: abs(x[2] - 50))  # 优先选50个基因的
        selected = candidates[:pathways_per_db]

        for name, genes, count in selected:
            test_pathways[f"{db_name}_{name}"] = genes

    return test_pathways

def generate_sorted_gene_list_from_pool(
    gene_pool: Set[str],
    n_genes: int = 500,
    seed: int = 42
) -> pd.DataFrame:
    """从基因池生成500基因排序列表"""
    np.random.seed(seed)

    # 如果基因池不足，使用全部
    available_genes = list(gene_pool)
    if len(available_genes) < n_genes:
        print(f"Warning: Gene pool ({len(available_genes)}) < requested {n_genes}")
        n_genes = len(available_genes)

    # 随机选取
    selected_genes = np.random.choice(available_genes, size=n_genes, replace=False)

    # 生成权重（模拟log2FC）
    direction = np.random.choice([-2, 0, 2], size=n_genes, p=[0.3, 0.4, 0.3])
    weights = np.random.randn(n_genes) * 1.5 + direction

    df = pd.DataFrame({
        'gene': selected_genes,
        'weight': weights
    })

    df = df.sort_values('weight', ascending=False).reset_index(drop=True)
    df['rank'] = range(1, len(df) + 1)

    return df

def generate_expression_matrix_from_pool(
    gene_pool: Set[str],
    n_genes: int = 6000,
    n_samples: int = 6,
    seed: int = 42
) -> pd.DataFrame:
    """从基因池生成6000×6表达矩阵"""
    np.random.seed(seed)

    available_genes = list(gene_pool)
    if len(available_genes) < n_genes:
        print(f"Warning: Gene pool ({len(available_genes)}) < requested {n_genes}")
        n_genes = len(available_genes)

    # 选取基因
    selected_genes = np.random.choice(available_genes, size=n_genes, replace=False)

    # 样本名
    sample_names = [f"Normal_{i+1}" for i in range(3)] + [f"Disease_{i+1}" for i in range(3)]

    # 生成表达值
    normal_expr = np.random.randn(n_genes, 3) * 1.5 + 8
    disease_expr = np.random.randn(n_genes, 3) * 1.5 + 8

    # 添加差异表达
    deg_indices = np.random.choice(range(n_genes), size=100, replace=False)
    for idx in deg_indices[:50]:
        disease_expr[idx, :] += np.random.uniform(1, 2, 3)
    for idx in deg_indices[50:]:
        disease_expr[idx, :] -= np.random.uniform(1, 2, 3)

    expr_matrix = np.hstack([normal_expr, disease_expr])
    df = pd.DataFrame(expr_matrix, index=selected_genes, columns=sample_names)

    return df

if __name__ == "__main__":
    print("=" * 60)
    print("GSEA/GSVA/ssGSEA 端对端测试数据生成器")
    print("=" * 60)
    print()

    ranked_genes, expr_matrix, gene_sets = save_test_data()

    print()
    print("=" * 60)
    print("数据统计")
    print("=" * 60)
    print(f"排序基因列表: {len(ranked_genes)} genes")
    print(f"  - Top 5 (上调): {ranked_genes.head()['gene'].tolist()}")
    print(f"  - Bottom 5 (下调): {ranked_genes.tail()['gene'].tolist()}")
    print()
    print(f"表达矩阵: {expr_matrix.shape[0]} genes x {expr_matrix.shape[1]} samples")
    print(f"  - 样本: {expr_matrix.columns.tolist()}")
    print(f"  - 表达范围: [{expr_matrix.values.min():.2f}, {expr_matrix.values.max():.2f}]")
    print()
    print(f"基因集: {len(gene_sets)} sets")
    for name, genes in list(gene_sets.items())[:3]:
        print(f"  - {name}: {len(genes)} genes")
