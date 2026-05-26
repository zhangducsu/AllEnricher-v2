"""
GSVA（Gene Set Variation Analysis）基因集变异分析模块

中文模块说明：
    本模块实现了 GSVA 算法，用于将基因水平的表达矩阵转换为通路水平的活性矩阵。
    GSVA 是一种非参数方法，不需要预先定义样本分组，适用于单样本和多样本分析。

    主要组件：
    - GSVA: GSVA 分析类，支持三种方法变体（gsva / plage / zscore）

    三种方法变体：
    1. gsva:  标准 GSVA（随机游走统计量），基于经验累积分布函数和核密度估计
    2. plage: 主成分分析（PLAGE），取通路基因表达矩阵的第一主成分得分
    3. zscore: Z-score 方法，计算基因集内基因表达均值的标准分

    参考文献：
    - Hänzelmann S, Castelo R, Guinney J. GSVA: gene set variation analysis
      for microarray and RNA-seq data. BMC Bioinformatics. 2013;14:7.
"""

import logging
from abc import ABC, abstractmethod
from dataclasses import dataclass
from typing import TYPE_CHECKING, Dict, List, Optional, Set, Tuple, Any

import numpy as np
import pandas as pd
from scipy.stats import norm

# 使用 getLogger 获取模块级 logger
logger = logging.getLogger(__name__)

# TYPE_CHECKING 下导入类型，避免循环依赖
# enrichment.py 在模块顶层导入了 gsva.py，因此 gsva.py 不能在模块加载时导入 enrichment.py
if TYPE_CHECKING:
    from allenricher.core.enrichment import EnrichmentMethodBase, EnrichmentResult


class GSVA(ABC):
    """
    基因集变异分析（Gene Set Variation Analysis, GSVA）

    GSVA 是一种非参数方法，将基因水平的表达矩阵转换为通路（基因集）水平的
    活性矩阵。与 ORA（过表示分析）和 GSEA 不同，GSVA 不需要预先定义样本分组，
    而是为每个样本独立计算每个基因集的活性得分。

    核心算法（标准 GSVA 方法）：
    1. 对每个基因，计算其在所有样本中的经验累积分布函数（ECDF）
       - ECDF_down: 基因表达值低于当前值的比例（衡量低表达程度）
       - ECDF_up: 基因表达值高于当前值的比例（衡量高表达程度）
    2. 使用核密度估计（KDE）将离散的 ECDF 值转换为连续的分数
       - 支持 Gaussian 和 Poisson 两种核函数
    3. 对基因集内的基因，计算随机游走统计量（Kolmogorov-Smirnov 类似统计量）
       - 遍历基因集中的基因，维护累积分数
       - 命中基因集时增加 ECDF_up 的 KDE 分数，否则减少 ECDF_down 的 KDE 分数
    4. 返回每个样本的通路活性得分（随机游走统计量的最大偏离值）

    属性:
        method: 方法变体名称（"gsva" / "plage" / "zscore"）
        kcdf: 核密度估计的核函数类型（"Gaussian" / "Poisson"）
        tau: 核密度估计的带宽参数（仅对 Gaussian 核有效）
        min_size: 基因集的最小允许大小
        max_size: 基因集的最大允许大小
    """

    def __init__(
        self,
        method: str = "gsva",
        kcdf: str = "Gaussian",
        tau: float = 1.0,
        min_size: int = 10,
        max_size: int = 500
    ):
        """
        初始化 GSVA 分析器

        参数:
            method: 方法变体名称，支持 "gsva"（标准随机游走）、
                    "plage"（主成分分析）、"zscore"（Z-score 方法）
            kcdf: 核密度估计的核函数类型，支持 "Gaussian" 和 "Poisson"
            tau: 核密度估计的带宽参数（仅对 Gaussian 核有效，默认 1.0）
            min_size: 基因集最小大小，小于此值的基因集将被跳过（默认 10）
            max_size: 基因集最大大小，大于此值的基因集将被跳过（默认 500）

        异常:
            ValueError: 当 method 或 kcdf 参数值不合法时抛出
        """
        # 验证方法变体名称
        valid_methods = ("gsva", "plage", "zscore")
        if method not in valid_methods:
            raise ValueError(
                f"不支持的 GSVA 方法变体: '{method}'，"
                f"有效值为: {valid_methods}"
            )
        # 验证核函数类型
        valid_kcdf = ("Gaussian", "Poisson")
        if kcdf not in valid_kcdf:
            raise ValueError(
                f"不支持的核函数类型: '{kcdf}'，"
                f"有效值为: {valid_kcdf}"
            )

        self.method = method  # 方法变体名称
        self.kcdf = kcdf  # 核密度估计核函数类型
        self.tau = tau  # 带宽参数
        self.min_size = min_size  # 基因集最小大小
        self.max_size = max_size  # 基因集最大大小

    def calculate_pvalue(
        self,
        gene_count: int,
        background_count: int,
        gene_total: int,
        background_total: int
    ) -> float:
        """
        GSVA 不使用传统的 p 值计算

        注意：GSVA 是非参数方法，不进行统计检验，因此此方法返回 NaN 作为占位值。
        GSVA 的核心输出是通路活性得分矩阵，而非 p 值。

        参数:
            gene_count: 未使用
            background_count: 未使用
            gene_total: 未使用
            background_total: 未使用

        返回值:
            float: 固定返回 NaN（不适用）
        """
        # GSVA 不进行统计检验，p 值不适用
        return float('nan')

    def calculate_enrichment(
        self,
        gene_set: Set[str],
        background_set: Set[str],
        term_genes: Set[str],
        term_name: str,
        term_id: str,
        database: str,
        ranked_genes: Optional[List[str]] = None,
        gene_weights: Optional[Dict[str, float]] = None
    ) -> Optional["EnrichmentResult"]:
        """
        使用 GSVA 方法对单个条目执行富集分析（兼容基类接口）

        注意：此方法为兼容 EnrichmentMethodBase 基类接口而实现。
        GSVA 的核心功能通过 analyze_matrix 方法使用，该方法直接对表达矩阵
        计算通路活性得分。此方法仅返回一个占位结果。

        参数:
            gene_set: 输入基因集合（查询基因列表）
            background_set: 背景基因集合
            term_genes: 当前条目所包含的基因集合
            term_name: 条目名称
            term_id: 条目ID
            database: 数据库名称
            ranked_genes: 未使用（GSVA 使用表达矩阵而非排序列表）
            gene_weights: 未使用（GSVA 使用表达矩阵而非权重字典）

        返回值:
            Optional[EnrichmentResult]: 占位富集结果对象；如果基因集大小
                                       不满足要求则返回 None
        """
        # 延迟导入以避免循环依赖
        from allenricher.core.enrichment import EnrichmentResult, generate_term_url

        # 检查基因集大小是否在允许范围内
        overlap = gene_set & term_genes
        if len(overlap) < self.min_size or len(overlap) > self.max_size:
            return None

        # 生成条目的数据库链接 URL
        term_url = generate_term_url(term_id, database)

        # 返回占位结果（GSVA 的核心输出是 analyze_matrix 的活性矩阵）
        return EnrichmentResult(
            term_id=term_id,
            term_name=term_name,
            database=database,
            pvalue=float('nan'),  # GSVA 不计算 p 值
            adjusted_pvalue=float('nan'),
            gene_count=len(overlap),
            background_count=len(term_genes),
            expected_count=0,  # GSVA 不使用期望值概念
            rich_factor=0,  # GSVA 不使用富集因子概念
            gene_list=list(overlap),
            gene_ratio=f"{len(overlap)}/{len(gene_set)}",
            background_ratio=f"{len(term_genes)}/{len(background_set)}",
            term_url=term_url,
            nes=float('nan'),
            es=float('nan'),
            fdr=float('nan'),
            leading_edge=None
        )

    def _compute_ecdf_kde(
        self,
        expression_matrix: np.ndarray,
        kcdf: str = "Gaussian",
        tau: float = 1.0
    ) -> Tuple[np.ndarray, np.ndarray]:
        """
        计算每个基因在每个样本中的经验累积分布函数（ECDF）的核密度估计分数

        对每个基因，在所有样本上计算两个方向的 ECDF：
        - ECDF_down: 基因表达值 <= 当前值的比例（衡量低表达程度）
        - ECDF_up: 基因表达值 > 当前值的比例（衡量高表达程度）

        然后使用核密度估计将离散的 ECDF 值转换为连续分数。

        参数:
            expression_matrix: 表达矩阵，形状为 (n_genes, n_samples)
            kcdf: 核函数类型（"Gaussian" 或 "Poisson"）
            tau: 带宽参数（仅对 Gaussian 核有效）

        返回值:
            Tuple[np.ndarray, np.ndarray]: 元组，包含两个元素：
                - ecdf_down: 低表达方向的 KDE 分数矩阵，形状为 (n_genes, n_samples)
                - ecdf_up: 高表达方向的 KDE 分数矩阵，形状为 (n_genes, n_samples)
        """
        n_genes, n_samples = expression_matrix.shape

        # 初始化 ECDF 矩阵
        ecdf_down = np.zeros_like(expression_matrix, dtype=float)
        ecdf_up = np.zeros_like(expression_matrix, dtype=float)

        for i in range(n_genes):
            # 获取当前基因在所有样本中的表达值
            gene_expr = expression_matrix[i, :]

            if kcdf == "Gaussian":
                # Gaussian 核密度估计
                # 计算基因在所有样本中的标准差
                std = np.std(gene_expr, ddof=1)
                if std == 0:
                    # 标准差为 0 时（所有样本表达值相同），ECDF 为阶跃函数
                    # 使用0.5作为中性值，使得随机游走时该基因不贡献净变化
                    ecdf_down[i, :] = 0.5
                    ecdf_up[i, :] = 0.5
                    continue

                # 计算所有样本的z-score
                mean = np.mean(gene_expr)
                z_scores = (gene_expr - mean) / std

                for j in range(n_samples):
                    # 使用 Gaussian 核平滑 ECDF
                    # 对于每个样本，计算其相对于其他样本的核密度估计
                    z_j = z_scores[j]

                    # 计算该样本相对于所有其他样本的核密度权重
                    # 使用高斯核函数: exp(-(z_i - z_j)^2 / (2*tau^2))
                    kernel_weights = np.exp(-((z_scores - z_j) ** 2) / (2 * tau ** 2))
                    kernel_weights = kernel_weights / np.sum(kernel_weights)

                    # ECDF_down: 表达值 <= 当前值的累积概率
                    # ECDF_up: 表达值 > 当前值的累积概率
                    ecdf_down[i, j] = np.sum(kernel_weights * (z_scores <= z_j))
                    ecdf_up[i, j] = np.sum(kernel_weights * (z_scores > z_j))

            elif kcdf == "Poisson":
                # Poisson 核密度估计
                # 使用简单的 ECDF（无平滑），适用于计数数据
                for j in range(n_samples):
                    ecdf_down[i, j] = np.sum(gene_expr <= gene_expr[j]) / n_samples
                    ecdf_up[i, j] = np.sum(gene_expr > gene_expr[j]) / n_samples

        return ecdf_down, ecdf_up

    def _compute_gsva_score(
        self,
        gene_set_genes: List[str],
        gene_to_idx: Dict[str, int],
        ecdf_down: np.ndarray,
        ecdf_up: np.ndarray,
        n_genes_total: int
    ) -> np.ndarray:
        """
        计算单个基因集在所有样本中的 GSVA 活性得分（随机游走统计量）

        算法步骤（基于原始GSVA论文）：
        1. 确定基因集中与表达矩阵有交集的基因
        2. 对每个样本，根据基因表达水平排序（通过ECDF分数间接反映）
        3. 进行随机游走：
           - 遇到属于基因集的基因：累积分数 += ecdf_up - ecdf_down
           - 遇到不属于基因集的基因：累积分数 -= ecdf_up - ecdf_down
        4. 活性得分 = 累积分数的最大偏离值（Kolmogorov-Smirnov统计量）

        参数:
            gene_set_genes: 基因集中的基因列表
            gene_to_idx: 基因名称到表达矩阵行索引的映射
            ecdf_down: 低表达方向的 KDE 分数矩阵，形状为 (n_genes, n_samples)
            ecdf_up: 高表达方向的 KDE 分数矩阵，形状为 (n_genes, n_samples)
            n_genes_total: 表达矩阵中的基因总数

        返回值:
            np.ndarray: 该基因集在所有样本中的 GSVA 活性得分，形状为 (n_samples,)
        """
        n_samples = ecdf_down.shape[1]

        # 获取基因集中存在于表达矩阵中的基因索引
        gene_indices = [gene_to_idx[g] for g in gene_set_genes if g in gene_to_idx]
        gene_index_set = set(gene_indices)

        if len(gene_indices) == 0:
            return np.zeros(n_samples)

        # 初始化结果数组
        scores = np.zeros(n_samples)

        # 计算基因集内和基因集外的基因数量
        n_genes_in_set = len(gene_index_set)
        n_genes_not_in_set = n_genes_total - n_genes_in_set

        if n_genes_in_set == 0 or n_genes_not_in_set == 0:
            return np.zeros(n_samples)

        for j in range(n_samples):
            # 计算每个基因的ECDF差值（反映该基因在样本j中的相对表达水平）
            # 差值越大，表示该基因在该样本中表达越高
            ecdf_diff = ecdf_up[:, j] - ecdf_down[:, j]

            # 根据ECDF差值降序排序基因（高表达在前）
            sorted_indices = np.argsort(-ecdf_diff)

            running_sum = 0.0
            max_deviation = 0.0

            # 随机游走：遍历按表达水平排序的基因
            for idx in sorted_indices:
                if idx in gene_index_set:
                    # 命中：基因属于当前基因集
                    # 增加步长 = 1/基因集大小
                    running_sum += 1.0 / n_genes_in_set
                else:
                    # 未命中：基因不属于当前基因集
                    # 减少步长 = 1/(总基因数 - 基因集大小)
                    running_sum -= 1.0 / n_genes_not_in_set

                # 跟踪最大偏离值
                abs_deviation = abs(running_sum)
                if abs_deviation > max_deviation:
                    max_deviation = abs_deviation

            scores[j] = max_deviation

        return scores

    def _compute_plage_score(
        self,
        gene_set_genes: List[str],
        gene_to_idx: Dict[str, int],
        expression_matrix: np.ndarray
    ) -> np.ndarray:
        """
        计算单个基因集的 PLAGE 活性得分（第一主成分得分）

        PLAGE（Pathway Level Analysis of Gene Expression）方法：
        1. 提取基因集中基因的表达子矩阵（n_genes_in_set x n_samples）
        2. 对每个基因进行 z-score 标准化（减均值、除标准差）
        3. 对标准化后的矩阵执行 SVD，取第一主成分的得分作为通路活性

        参数:
            gene_set_genes: 基因集中的基因列表
            gene_to_idx: 基因名称到表达矩阵行索引的映射
            expression_matrix: 表达矩阵，形状为 (n_genes, n_samples)

        返回值:
            np.ndarray: 该基因集在所有样本中的 PLAGE 活性得分，形状为 (n_samples,)
        """
        # 获取基因集中存在于表达矩阵中的基因索引
        gene_indices = [gene_to_idx[g] for g in gene_set_genes if g in gene_to_idx]

        if len(gene_indices) == 0:
            n_samples = expression_matrix.shape[1]
            return np.zeros(n_samples)

        # 提取基因集的表达子矩阵
        sub_matrix = expression_matrix[gene_indices, :]  # (n_genes_in_set, n_samples)

        # 对每个基因进行 z-score 标准化
        means = np.mean(sub_matrix, axis=1, keepdims=True)
        stds = np.std(sub_matrix, axis=1, ddof=1, keepdims=True)
        # 避免除以 0
        stds[stds == 0] = 1.0
        z_matrix = (sub_matrix - means) / stds

        # 执行 SVD，取第一主成分得分（使用 numpy SVD，避免 sklearn 依赖）
        n_components = min(len(gene_indices), z_matrix.shape[1])
        if n_components == 0:
            return np.zeros(expression_matrix.shape[1])

        try:
            # SVD: z_matrix.T 的形状为 (n_samples, n_genes_in_set)
            U, S, Vt = np.linalg.svd(z_matrix.T, full_matrices=False)
            # 第一主成分得分 = U[:, 0] * S[0]
            scores = U[:, 0] * S[0]
        except np.linalg.LinAlgError:
            # SVD 失败时使用基因集内标准化表达均值作为后备
            scores = np.mean(z_matrix, axis=0)

        return scores

    def _compute_zscore(
        self,
        gene_set_genes: List[str],
        gene_to_idx: Dict[str, int],
        expression_matrix: np.ndarray
    ) -> np.ndarray:
        """
        计算单个基因集的 Z-score 活性得分

        Z-score 方法：
        1. 对每个基因进行 z-score 标准化（减均值、除标准差）
        2. 计算基因集内所有基因标准化表达值的均值
        3. 该均值即为通路活性得分

        参数:
            gene_set_genes: 基因集中的基因列表
            gene_to_idx: 基因名称到表达矩阵行索引的映射
            expression_matrix: 表达矩阵，形状为 (n_genes, n_samples)

        返回值:
            np.ndarray: 该基因集在所有样本中的 Z-score 活性得分，形状为 (n_samples,)
        """
        n_samples = expression_matrix.shape[1]

        # 获取基因集中存在于表达矩阵中的基因索引
        gene_indices = [gene_to_idx[g] for g in gene_set_genes if g in gene_to_idx]

        if len(gene_indices) == 0:
            return np.zeros(n_samples)

        # 提取基因集的表达子矩阵
        sub_matrix = expression_matrix[gene_indices, :]  # (n_genes_in_set, n_samples)

        # 对每个基因进行 z-score 标准化
        means = np.mean(sub_matrix, axis=1, keepdims=True)
        stds = np.std(sub_matrix, axis=1, ddof=1, keepdims=True)
        # 避免除以 0
        stds[stds == 0] = 1.0
        z_matrix = (sub_matrix - means) / stds

        # 计算基因集内所有基因标准化表达值的均值
        scores = np.mean(z_matrix, axis=0)

        return scores

    def analyze_matrix(
        self,
        expression_matrix: pd.DataFrame,
        gene_sets: Dict[str, Set[str]]
    ) -> pd.DataFrame:
        """
        对表达矩阵执行 GSVA 分析，输出样本 x 通路活性矩阵

        这是 GSVA 的核心方法。输入一个基因 x 样本的基因表达矩阵，
        输出一个通路 x 样本的活性得分矩阵。

        分析流程：
        1. 验证输入表达矩阵的有效性
        2. 根据方法变体选择计算策略：
           - gsva: 计算每个基因的 ECDF KDE 分数，然后执行随机游走统计量
           - plage: 对每个基因集的表达子矩阵执行 PCA，取第一主成分得分
           - zscore: 对每个基因标准化后取基因集内均值
        3. 过滤不满足大小要求的基因集
        4. 返回通路 x 样本的活性得分矩阵

        参数:
            expression_matrix: 基因表达矩阵，行为基因（index），列为样本（columns）
            gene_sets: 基因集字典，键为基因集名称，值为基因名称的集合

        返回值:
            pd.DataFrame: 通路活性得分矩阵，行为通路名称，列为样本名称；
                         如果输入为空矩阵则返回空 DataFrame

        示例:
            >>> import pandas as pd
            >>> expr = pd.DataFrame(...)  # 基因 x 样本表达矩阵
            >>> gene_sets = {"pathway1": {"GeneA", "GeneB"}, "pathway2": {"GeneC"}}
            >>> gsva = GSVA(method="gsva")
            >>> activity = gsva.analyze_matrix(expr, gene_sets)
            >>> print(activity.shape)  # (n_pathways, n_samples)
        """
        # --- 输入验证 ---
        if expression_matrix.empty:
            logger.warning("输入表达矩阵为空，返回空 DataFrame")
            return pd.DataFrame()

        # 获取样本名称
        sample_names = expression_matrix.columns.tolist()
        n_samples = len(sample_names)
        n_genes = len(expression_matrix)

        # 构建基因名称到行索引的映射
        gene_to_idx = {gene: idx for idx, gene in enumerate(expression_matrix.index)}

        # --- 根据方法变体选择计算策略 ---
        if self.method == "gsva":
            # 标准 GSVA：计算 ECDF KDE 分数 + 随机游走统计量
            expr_array = expression_matrix.values.astype(float)

            # 计算所有基因的 ECDF KDE 分数
            logger.info("计算 ECDF 核密度估计分数...")
            ecdf_down, ecdf_up = self._compute_ecdf_kde(
                expr_array, kcdf=self.kcdf, tau=self.tau
            )

            # 对每个基因集计算 GSVA 活性得分
            results = {}
            for pathway_name, genes in gene_sets.items():
                # 获取基因集中与表达矩阵有交集的基因
                overlap = [g for g in genes if g in gene_to_idx]

                # 检查基因集大小是否在允许范围内
                if len(overlap) < self.min_size or len(overlap) > self.max_size:
                    logger.debug(
                        f"跳过基因集 '{pathway_name}'："
                        f"交集基因数 {len(overlap)} 不在 [{self.min_size}, {self.max_size}] 范围内"
                    )
                    continue

                # 计算随机游走统计量
                scores = self._compute_gsva_score(
                    gene_set_genes=overlap,
                    gene_to_idx=gene_to_idx,
                    ecdf_down=ecdf_down,
                    ecdf_up=ecdf_up,
                    n_genes_total=n_genes
                )
                results[pathway_name] = scores

        elif self.method == "plage":
            # PLAGE 方法：第一主成分得分
            expr_array = expression_matrix.values.astype(float)

            results = {}
            for pathway_name, genes in gene_sets.items():
                overlap = [g for g in genes if g in gene_to_idx]

                if len(overlap) < self.min_size or len(overlap) > self.max_size:
                    logger.debug(
                        f"跳过基因集 '{pathway_name}'："
                        f"交集基因数 {len(overlap)} 不在 [{self.min_size}, {self.max_size}] 范围内"
                    )
                    continue

                scores = self._compute_plage_score(
                    gene_set_genes=overlap,
                    gene_to_idx=gene_to_idx,
                    expression_matrix=expr_array
                )
                results[pathway_name] = scores

        elif self.method == "zscore":
            # Z-score 方法：标准化后取均值
            expr_array = expression_matrix.values.astype(float)

            results = {}
            for pathway_name, genes in gene_sets.items():
                overlap = [g for g in genes if g in gene_to_idx]

                if len(overlap) < self.min_size or len(overlap) > self.max_size:
                    logger.debug(
                        f"跳过基因集 '{pathway_name}'："
                        f"交集基因数 {len(overlap)} 不在 [{self.min_size}, {self.max_size}] 范围内"
                    )
                    continue

                scores = self._compute_zscore(
                    gene_set_genes=overlap,
                    gene_to_idx=gene_to_idx,
                    expression_matrix=expr_array
                )
                results[pathway_name] = scores

        # --- 构建结果 DataFrame ---
        if not results:
            logger.warning("没有基因集满足大小要求，返回空 DataFrame")
            return pd.DataFrame()

        # 构建通路 x 样本的活性得分矩阵
        activity_df = pd.DataFrame(results, index=sample_names).T
        activity_df.index.name = "Pathway"

        logger.info(f"GSVA 分析完成：{len(results)} 个通路 x {n_samples} 个样本")
        return activity_df
