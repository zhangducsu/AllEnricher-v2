# 扩展 GSEA、GSVA、ssGSEA 功能计划

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task.

**Goal:** 在 AllEnricher v2 中扩展 GSEA、GSVA、ssGSEA 三种基因集富集分析方法，支持表达矩阵输入和样本水平分析。

**Architecture:** 
- 基于现有 `GSEA` 类扩展，新增 `GSVA` 和 `ssGSEA` 类
- 统一输入接口：支持基因列表（ORA）和排序基因列表/表达矩阵（GSEA家族）
- 结果输出格式与现有富集分析一致（DataFrame）

**Tech Stack:** Python 3.8+, numpy, pandas, scipy

---

## 功能对比

| 方法 | 输入 | 输出 | 适用场景 |
|------|------|------|----------|
| **ORA** (Fisher/Hyper) | 基因列表 | 条目富集结果 | 差异基因集分析 |
| **GSEA** | 排序基因列表 (带权重) | 条目富集结果 | 全基因组表达趋势分析 |
| **ssGSEA** | 表达矩阵 (单样本) | 样本-条目活性矩阵 | 单样本通路活性评估 |
| **GSVA** | 表达矩阵 (多样本) | 样本-条目活性矩阵 | 样本间通路活性比较 |

---

## 任务分解

### Task 1: 扩展现有 GSEA 类

**Files:**
- Modify: `allenricher/core/enrichment.py`

**目的**: 完善 GSEA 实现，添加 NES 计算和可视化支持

- [ ] **Step 1: 添加 NES (Normalized Enrichment Score) 计算**

```python
def calculate_normalized_es(
    self,
    ranked_genes: List[str],
    gene_set: Set[str],
    gene_weights: Optional[Dict[str, float]] = None
) -> Tuple[float, float, List[str]]:
    """
    计算归一化富集分数 (NES)
    
    步骤:
    1. 计算实际 ES
    2. 进行 permutations 次置换，计算 null ES 分布
    3. NES = ES / mean(|ES_null|)
    4. pvalue = 置换检验中 |ES_null| >= |ES| 的比例
    
    Returns:
        (ES, NES, pvalue, leading_edge_genes)
    """
```

- [ ] **Step 2: 添加 GSEA 结果存储类**

```python
@dataclass
class GSEAResult(EnrichmentResult):
    """GSEA 专用结果类，额外存储 ES 和 NES"""
    es: float  # Enrichment Score
    nes: float  # Normalized Enrichment Score
    leading_edge: List[str]  # 前沿基因列表
    running_sum: List[float]  # 累积分数轨迹（用于可视化）
```

- [ ] **Step 3: 运行测试验证 GSEA 正确性**

Run: `pytest tests/test_enrichment.py::TestGSEA -v`

---

### Task 2: 实现 ssGSEA (单样本 GSEA)

**Files:**
- Create: `allenricher/core/ssgsea.py`

**目的**: 实现单样本基因集富集分析，评估每个样本中各通路的活性

- [ ] **Step 1: 创建 ssGSEA 类**

```python
class ssGSEA:
    """
    单样本基因集富集分析 (single-sample GSEA)
    
    参考: Barbie et al., Nature, 2009
    
    核心思想:
    - 对每个样本独立进行 GSEA
    - 使用样本内基因表达值的秩次作为排序依据
    - 输出每个样本-通路的活性得分
    
    与 GSEA 的区别:
    - GSEA: 比较两组样本间的差异
    - ssGSEA: 评估单个样本中各通路的绝对活性
    """
    
    def __init__(self, alpha: float = 0.25, normalize: bool = True):
        """
        参数:
            alpha: 权重指数，控制基因权重对秩次的影响 (默认 0.25)
            normalize: 是否对结果进行归一化
        """
        self.alpha = alpha
        self.normalize = normalize
    
    def analyze(
        self,
        expression_matrix: pd.DataFrame,  # 行: 基因, 列: 样本
        gene_sets: Dict[str, Set[str]],   # 通路名 -> 基因集合
        min_size: int = 10,
        max_size: int = 500
    ) -> pd.DataFrame:
        """
        执行 ssGSEA 分析
        
        参数:
            expression_matrix: 表达矩阵，index 为基因名
            gene_sets: 基因集字典
            min_size/max_size: 基因集大小过滤
        
        返回:
            DataFrame: 行=通路, 列=样本, 值=通路活性得分
        """
```

- [ ] **Step 2: 实现 ssGSEA 核心算法**

```python
def _calculate_ssgsea_score(
    self,
    gene_ranks: pd.Series,  # 单个样本的基因秩次
    gene_set: Set[str]
) -> float:
    """
    计算单个样本中单个基因集的 ssGSEA 得分
    
    算法:
    1. 计算基因集中基因的秩次位置
    2. 计算加权秩次和 (使用 alpha 指数)
    3. 归一化得到最终得分
    """
    # 获取基因集中的基因
    genes_in_set = gene_ranks.index.intersection(gene_set)
    if len(genes_in_set) < 2:
        return 0.0
    
    # 计算秩次权重
    ranks = gene_ranks[genes_in_set].values
    n_genes = len(gene_ranks)
    
    # 加权秩次和
    weighted_sum = np.sum(np.abs(ranks) ** self.alpha)
    
    # 归一化
    if self.normalize:
        score = weighted_sum / len(genes_in_set) ** self.alpha
    else:
        score = weighted_sum
    
    return score
```

- [ ] **Step 3: 编写 ssGSEA 测试**

Create: `tests/test_ssgsea.py`

```python
def test_ssgsea_basic():
    """测试 ssGSEA 基本功能"""
    # 创建模拟表达矩阵
    np.random.seed(42)
    expr = pd.DataFrame(
        np.random.randn(100, 3),
        index=[f"Gene{i}" for i in range(100)],
        columns=["Sample1", "Sample2", "Sample3"]
    )
    
    # 创建基因集
    gene_sets = {
        "Pathway_A": set([f"Gene{i}" for i in range(20)]),
        "Pathway_B": set([f"Gene{i}" for i in range(20, 40)])
    }
    
    ssgsea = ssGSEA()
    results = ssgsea.analyze(expr, gene_sets)
    
    assert results.shape == (2, 3)  # 2 通路 x 3 样本
    assert not results.isnull().any().any()  # 无缺失值
```

---

### Task 3: 实现 GSVA (基因集变异分析)

**Files:**
- Create: `allenricher/core/gsva.py`

**目的**: 实现 GSVA 算法，用于多样本间的通路活性比较

- [ ] **Step 1: 创建 GSVA 类**

```python
class GSVA:
    """
    基因集变异分析 (Gene Set Variation Analysis)
    
    参考: Hänzelmann et al., BMC Bioinformatics, 2013
    
    核心思想:
    - 不依赖于预定义的样本分组
    - 使用核密度估计评估基因在样本间的表达分布
    - 计算随机游走统计量评估通路活性
    
    与 ssGSEA 的区别:
    - ssGSEA: 对每个样本独立计算
    - GSVA: 利用多样本信息，通过核密度估计更稳健地评估通路活性
    """
    
    def __init__(
        self,
        method: str = "gsva",  # 'gsva' 或 'plage' 或 'zscore'
        kcdf: str = "Gaussian",  # 核函数: 'Gaussian' 或 'Poisson'
        tau: float = 1.0,  # 核密度估计参数
        max_diff: bool = True,  # 是否使用最大差异作为通路得分
        abs_ranking: bool = False  # 是否使用绝对秩次
    ):
        self.method = method
        self.kcdf = kcdf
        self.tau = tau
        self.max_diff = max_diff
        self.abs_ranking = abs_ranking
```

- [ ] **Step 2: 实现 GSVA 核心算法**

```python
def _calculate_gsva_score(
    self,
    expression_matrix: pd.DataFrame,
    gene_set: Set[str]
) -> pd.Series:
    """
    计算单个基因集的 GSVA 得分
    
    算法步骤:
    1. 对表达矩阵进行基因水平的核密度估计
    2. 计算每个基因在样本间的累积分布函数 (ECDF)
    3. 对基因集内的基因，计算随机游走统计量
    4. 返回每个样本的通路活性得分
    """
    from scipy import stats
    from scipy.stats import norm
    
    genes_in_set = expression_matrix.index.intersection(gene_set)
    if len(genes_in_set) < 10:
        return pd.Series([np.nan] * expression_matrix.shape[1], 
                        index=expression_matrix.columns)
    
    # 提取基因集表达子矩阵
    expr_subset = expression_matrix.loc[genes_in_set]
    
    # 计算每个基因在样本间的秩次（核密度估计近似）
    ranks = expr_subset.rank(axis=1, method="average")
    
    # 归一化到 [0, 1]
    n_samples = expression_matrix.shape[1]
    normalized_ranks = (ranks - 1) / (n_samples - 1)
    
    # 计算随机游走统计量 (简化版)
    # 实际 GSVA 使用更复杂的核密度估计
    walk_stat = normalized_ranks.sum(axis=0) / len(genes_in_set)
    
    return walk_stat
```

- [ ] **Step 3: 编写 GSVA 测试**

Create: `tests/test_gsva.py`

---

### Task 4: 集成到 EnrichmentAnalyzer

**Files:**
- Modify: `allenricher/core/enrichment.py`

**目的**: 让主分析引擎支持 GSEA 家族方法

- [ ] **Step 1: 添加分析方法路由**

```python
def run_analysis(
    self,
    input_data: Union[Set[str], pd.DataFrame, Tuple[List[str], Dict[str, float]]],
    background_set: Optional[Set[str]] = None,
    database_data: Optional[Dict[str, Dict[str, Set[str]]]] = None,
    method: Optional[str] = None
) -> Dict[str, pd.DataFrame]:
    """
    执行富集分析（统一入口）
    
    参数:
        input_data: 输入数据，根据方法不同可以是:
            - ORA: Set[str] 基因列表
            - GSEA: Tuple[List[str], Dict[str, float]] (排序基因, 权重)
            - ssGSEA/GSVA: pd.DataFrame 表达矩阵
    """
    method = method or self.config.method
    
    if method in ["fisher", "hypergeometric"]:
        return self._run_ora(input_data, background_set, database_data)
    elif method == "gsea":
        return self._run_gsea(input_data, database_data)
    elif method == "ssgsea":
        from .ssgsea import ssGSEA
        analyzer = ssGSEA()
        return analyzer.analyze(input_data, database_data)
    elif method == "gsva":
        from .gsva import GSVA
        analyzer = GSVA()
        return analyzer.analyze(input_data, database_data)
```

- [ ] **Step 2: 更新配置类**

Modify: `allenricher/core/config.py`

```python
class EnrichmentMethod(Enum):
    """富集分析方法枚举"""
    FISHER = "fisher"
    HYPERGEOMETRIC = "hypergeometric"
    GSEA = "gsea"
    SSGSEA = "ssgsea"
    GSVA = "gsva"
```

---

### Task 5: CLI 支持

**Files:**
- Modify: `allenricher/cli.py`

**目的**: 命令行支持 GSEA 家族方法

- [ ] **Step 1: 添加表达矩阵输入参数**

```python
@click.option(
    "--expression-matrix", "-e",
    type=click.Path(exists=True),
    help="表达矩阵文件 (用于 GSEA/ssGSEA/GSVA)"
)
@click.option(
    "--ranked-genes", "-r",
    type=click.Path(exists=True),
    help="排序基因列表文件 (用于 GSEA)"
)
@click.option(
    "--method", "-m",
    type=click.Choice(["fisher", "hypergeometric", "gsea", "ssgsea", "gsva"]),
    default="fisher",
    help="富集分析方法"
)
```

---

### Task 6: 文档和示例

**Files:**
- Create: `examples/gsea_example.py`
- Create: `examples/ssgsea_example.py`
- Create: `examples/gsva_example.py`

- [ ] **Step 1: 编写 GSEA 示例**

```python
"""GSEA 分析示例"""
import pandas as pd
from allenricher import EnrichmentAnalyzer

# 1. 准备排序基因列表（带权重，如 log2 fold change）
ranked_genes = ["GeneA", "GeneB", "GeneC", ...]  # 按 fold change 降序排列
gene_weights = {"GeneA": 2.5, "GeneB": 2.1, ...}

# 2. 创建分析器
analyzer = EnrichmentAnalyzer(method="gsea")

# 3. 运行分析
results = analyzer.run_analysis(
    input_data=(ranked_genes, gene_weights),
    database_data={"GO": go_gene_sets, "KEGG": kegg_gene_sets}
)

# 4. 查看结果
print(results["GO"].head())
```

---

## 依赖项

```
numpy>=1.20.0
pandas>=1.3.0
scipy>=1.7.0
```

## 时间表

| Task | 预计时间 | 优先级 |
|------|---------|--------|
| Task 1: 完善 GSEA | 2h | High |
| Task 2: ssGSEA | 3h | High |
| Task 3: GSVA | 3h | Medium |
| Task 4: 集成 | 2h | High |
| Task 5: CLI | 1h | Medium |
| Task 6: 文档 | 1h | Low |

---

**Plan complete and saved to `docs/superpowers/plans/2026-05-26-extend-gsea-gsva-ssgsea.md`.**

**执行选项:**

1. **Subagent-Driven (推荐)** - 每个 Task 派发独立子代理
2. **Inline Execution** - 在当前会话顺序执行

**选择哪种方式?**