"""
转录因子富集分析模块

提供基于转录因子-靶基因关系的富集分析功能，支持 ORA（过表示分析）
和 GSEA（基因集富集分析）两种方法。

数据来源：
    - TRRUST: 包含 TF 调控模式（activator/repressor/mixed/unknown）
    - ChEA3: 仅包含 TF-靶基因关系，无调控模式

主要组件：
    - TFEnrichmentAnalyzer: 转录因子富集分析器
"""

import logging
from typing import Dict, List, Optional, Set, Tuple, TYPE_CHECKING

import numpy as np
import pandas as pd
from scipy.stats import hypergeom
from statsmodels.stats.multitest import multipletests

logger = logging.getLogger(__name__)

# TYPE_CHECKING 下导入类型，避免循环依赖
if TYPE_CHECKING:
    from allenricher.core.gsva import GSVA


class TFEnrichmentAnalyzer:
    """转录因子富集分析器

    基于 TF-靶基因关系数据库，对输入基因集执行转录因子富集分析。
    支持 ORA（过表示分析）和 GSEA（基因集富集分析）两种方法。

    数据库格式要求：
        tf_database 来自 DatabaseManager.load_trrust() 或 load_chea3()，
        包含以下键：
        - 'gene2tf': DataFrame，行为基因，列为 TF，值为 0/1 矩阵
        - 'tf_info': DataFrame（可选），包含 TF 的描述信息
          - TRRUST: 列为 [TF, mode, target_count]，mode 取值
            activator/repressor/mixed/unknown
          - ChEA3: 列为 [TF, lib_count, target_count]，无 mode 列

    Attributes:
        tf_to_targets: {TF: set_of_target_genes} 从 gene2tf 矩阵反向构建
        tf_info: TF 描述信息 DataFrame
        background_size: 背景基因总数（默认 20000，与 ChEA3 一致）
    """

    def __init__(self, tf_database: Dict, background_size: int = 20000):
        """初始化转录因子富集分析器

        参数:
            tf_database: 来自 DatabaseManager.load_trrust() 或 load_chea3() 的字典，
                          包含 'gene2tf' DataFrame 和可选的 'tf_info' DataFrame
            background_size: 背景基因总数，默认 20000（与 ChEA3 一致）

        异常:
            ValueError: tf_database 中缺少 'gene2tf' 键时抛出
        """
        if 'gene2tf' not in tf_database:
            raise ValueError(
                "tf_database 必须包含 'gene2tf' 键。"
                "请使用 DatabaseManager.load_trrust() 或 load_chea3() 获取数据。"
            )

        gene2tf_df = tf_database['gene2tf'].copy()
        self.tf_info = tf_database.get('tf_info', pd.DataFrame())
        self.background_size = background_size

        # 从 gene2tf 矩阵反向构建 tf_to_targets: {TF: set_of_targets}
        # gene2tf_df 的列为 TF，行为基因
        tf_columns = [col for col in gene2tf_df.columns if col != 'Gene']

        # 修复：确保数据列为整数类型，避免字符串 '1' 与整数 1 比较失败
        for col in tf_columns:
            gene2tf_df[col] = pd.to_numeric(gene2tf_df[col], errors='coerce').fillna(0).astype(int)

        self.tf_to_targets: Dict[str, Set[str]] = {}

        for tf in tf_columns:
            # 获取该 TF 列中值为 1 的所有基因
            mask = gene2tf_df[tf] == 1
            targets = set(gene2tf_df.loc[mask, 'Gene'].values) if 'Gene' in gene2tf_df.columns else set(gene2tf_df.index[mask].values)
            if targets:
                self.tf_to_targets[tf] = targets

        # 延迟导入 GSVA 类以避免循环依赖
        self._gsva_class = None

        logger.info(
            "TFEnrichmentAnalyzer 初始化完成: %d 个 TF, 背景大小 %d",
            len(self.tf_to_targets), background_size
        )

    def _get_gsva_class(self):
        """延迟获取 GSVA 类以避免循环依赖"""
        if self._gsva_class is None:
            from allenricher.core.gsva import GSVA
            self._gsva_class = GSVA
        return self._gsva_class

    def ora(
        self,
        gene_set: List[str],
        tf_list: Optional[List[str]] = None,
        min_overlap: int = 3
    ) -> pd.DataFrame:
        """执行转录因子过表示分析（ORA）

        使用超几何检验（Fisher 精确检验）评估每个 TF 的靶基因在输入基因集中
        是否显著富集。

        参数:
            gene_set: 输入基因列表（如差异表达基因）
            tf_list: 可选，指定只分析这些 TF；为 None 时分析所有 TF
            min_overlap: 最小重叠基因数，低于此值的 TF 将被跳过（默认 3）

        返回值:
            pd.DataFrame: 富集结果，按 Pvalue 排序，列包括：
                - TF: 转录因子名称
                - Overlap: 重叠基因数
                - TF_Targets: TF 的靶基因总数
                - GeneSet_Size: 输入基因集大小
                - Overlap_Genes: 重叠的基因列表（逗号分隔）
                - Pvalue: 超几何检验 p 值
                - FDR: BH 校正后的 FDR 值
                - Mode: 调控模式（TRRUST 特有，ChEA3 为 'unknown'）
        """
        gene_set_set = set(gene_set)
        n_gene_set = len(gene_set_set)

        # 确定要分析的 TF 列表
        if tf_list is not None:
            tfs_to_analyze = [tf for tf in tf_list if tf in self.tf_to_targets]
        else:
            tfs_to_analyze = list(self.tf_to_targets.keys())

        # 构建 TF mode 映射
        tf_mode_map = self._build_tf_mode_map()

        results = []
        for tf in tfs_to_analyze:
            targets = self.tf_to_targets[tf]
            overlap = gene_set_set & targets
            n_overlap = len(overlap)

            if n_overlap < min_overlap:
                continue

            n_targets = len(targets)

            # 超几何检验：scipy.stats.hypergeom.sf
            # P(X >= k) 其中 X ~ Hypergeom(M=background_size, n=n_targets, N=n_gene_set)
            # M = 背景总数, n = TF 靶基因数（成功标记数）, N = 基因集大小（抽取数）
            # k = 重叠数
            pvalue = hypergeom.sf(
                n_overlap - 1,
                self.background_size,
                n_targets,
                n_gene_set
            )

            mode = tf_mode_map.get(tf, 'unknown')

            results.append({
                'TF': tf,
                'Overlap': n_overlap,
                'TF_Targets': n_targets,
                'GeneSet_Size': n_gene_set,
                'Overlap_Genes': ','.join(sorted(overlap)),
                'Pvalue': pvalue,
                'Mode': mode,
            })

        if not results:
            logger.warning("ORA 分析未找到显著富集的 TF（min_overlap=%d）", min_overlap)
            return pd.DataFrame(
                columns=['TF', 'Overlap', 'TF_Targets', 'GeneSet_Size',
                         'Overlap_Genes', 'Pvalue', 'FDR', 'Mode']
            )

        result_df = pd.DataFrame(results)

        # FDR 校正（Benjamini-Hochberg）
        if len(result_df) > 1:
            reject, fdr, _, _ = multipletests(
                result_df['Pvalue'].values,
                method='fdr_bh'
            )
            result_df['FDR'] = fdr
        else:
            result_df['FDR'] = result_df['Pvalue'].values

        # 按 Pvalue 排序
        result_df = result_df.sort_values('Pvalue').reset_index(drop=True)

        logger.info("ORA 分析完成: %d 个 TF 有显著富集", len(result_df))
        return result_df

    def gsea(
        self,
        ranked_genes: List[Tuple[str, float]],
        tf_list: Optional[List[str]] = None,
        n_permutations: int = 1000,
        seed: int = 42
    ) -> pd.DataFrame:
        """执行转录因子 GSEA 分析

        使用标准 GSEA 算法评估每个 TF 的靶基因在排序列表中的富集程度。

        参数:
            ranked_genes: 排好序的基因列表，每个元素为 (gene_name, score) 元组，
                          按 score 降序排列
            tf_list: 可选，指定只分析这些 TF；为 None 时分析所有 TF
            n_permutations: 置换检验次数（默认 1000）
            seed: 随机种子（默认 42）

        返回值:
            pd.DataFrame: GSEA 结果，按 NES 绝对值降序排列，列包括：
                - TF: 转录因子名称
                - ES: 富集分数（Enrichment Score）
                - NES: 归一化富集分数（Normalized Enrichment Score）
                - Pvalue: 置换检验 p 值
                - FDR: BH 校正后的 FDR 值
        """
        # 提取基因名和分数
        gene_names = [g for g, _ in ranked_genes]
        gene_scores = {g: s for g, s in ranked_genes}
        n_genes = len(gene_names)

        # 确定要分析的 TF 列表
        if tf_list is not None:
            tfs_to_analyze = [tf for tf in tf_list if tf in self.tf_to_targets]
        else:
            tfs_to_analyze = list(self.tf_to_targets.keys())

        # 构建 TF mode 映射
        tf_mode_map = self._build_tf_mode_map()

        results = []
        rng = np.random.default_rng(seed)
        gene_names_array = np.array(gene_names)

        for tf in tfs_to_analyze:
            targets = self.tf_to_targets[tf]
            gene_set = targets & set(gene_names)
            nh = len(gene_set)

            if nh < 3:
                continue

            # 计算 ES（与现有 GSEA 相同的计算方式）
            es = self._calculate_es(gene_names, gene_set, n_genes, nh)

            if es == 0.0:
                continue

            # 置换检验
            null_es_list = []
            count_ge = 0

            for _ in range(n_permutations):
                permuted = rng.permutation(gene_names_array).tolist()
                permuted_es = self._calculate_es(permuted, gene_set, n_genes, nh)
                null_es_list.append(permuted_es)
                if abs(permuted_es) >= abs(es):
                    count_ge += 1

            # 计算 NES
            null_es_arr = np.array(null_es_list)
            null_pos = null_es_arr[null_es_arr >= 0]
            null_neg = np.abs(null_es_arr[null_es_arr < 0])
            mean_pos = np.mean(null_pos) if len(null_pos) > 0 else 1.0
            mean_neg = np.mean(null_neg) if len(null_neg) > 0 else 1.0

            if es > 0:
                nes = es / mean_pos if mean_pos > 0 else 0.0
            else:
                nes = -es / mean_neg if mean_neg > 0 else 0.0
                nes = -nes  # 保持负号

            pvalue = (count_ge + 1) / (n_permutations + 1)
            mode = tf_mode_map.get(tf, 'unknown')

            results.append({
                'TF': tf,
                'ES': es,
                'NES': nes,
                'Pvalue': pvalue,
                'Mode': mode,
            })

        if not results:
            logger.warning("GSEA 分析未找到显著富集的 TF")
            return pd.DataFrame(
                columns=['TF', 'ES', 'NES', 'Pvalue', 'FDR', 'Mode']
            )

        result_df = pd.DataFrame(results)

        # FDR 校正
        if len(result_df) > 1:
            _, fdr, _, _ = multipletests(
                result_df['Pvalue'].values,
                method='fdr_bh'
            )
            result_df['FDR'] = fdr
        else:
            result_df['FDR'] = result_df['Pvalue'].values

        # 按 NES 绝对值降序排列
        result_df = result_df.assign(
            _abs_nes=result_df['NES'].abs()
        ).sort_values('_abs_nes', ascending=False).drop(columns=['_abs_nes']).reset_index(drop=True)

        logger.info("GSEA 分析完成: %d 个 TF", len(result_df))
        return result_df

    def ssgsea(
        self,
        expression_df: pd.DataFrame,
        tf_list: Optional[List[str]] = None,
        min_size: int = 3,
        max_size: int = 500
    ) -> pd.DataFrame:
        """单样本 GSEA 分析

        对每个样本独立计算每个 TF 的 ssGSEA 得分，基于 TF-target 关系作为基因集。
        ssGSEA 不需要排列检验，直接计算归一化富集分数（NES）。

        参数:
            expression_df: 基因 x 样本的表达矩阵 (DataFrame)，行为基因，列为样本
            tf_list: 可选，指定只分析这些 TF；为 None 时分析所有 TF
            min_size: 基因集最小大小，低于此值的 TF 将被跳过（默认 3）
            max_size: 基因集最大大小，高于此值的 TF 将被跳过（默认 500）

        返回值:
            DataFrame: TF x 样本的 ssGSEA 得分矩阵，值为 NES（归一化富集分数）
        """
        if expression_df.empty:
            logger.warning("输入表达矩阵为空，返回空 DataFrame")
            return pd.DataFrame()

        # 确定要分析的 TF 列表
        if tf_list is not None:
            tfs_to_analyze = [tf for tf in tf_list if tf in self.tf_to_targets]
        else:
            tfs_to_analyze = list(self.tf_to_targets.keys())

        samples = expression_df.columns.tolist()
        results = {}

        for tf in tfs_to_analyze:
            targets = self.tf_to_targets[tf]
            # 获取与表达矩阵有交集的靶基因
            overlap = list(targets & set(expression_df.index))

            # 检查基因集大小是否在允许范围内
            if len(overlap) < min_size or len(overlap) > max_size:
                logger.debug(
                    f"跳过 TF '{tf}'：交集基因数 {len(overlap)} 不在 [{min_size}, {max_size}] 范围内"
                )
                continue

            nes_values = []
            for sample in samples:
                # 获取当前样本的表达量并按降序排序基因
                sample_expr = expression_df[sample]
                ranked_genes = sample_expr.sort_values(ascending=False).index.tolist()

                # 将表达量作为权重
                weights = sample_expr.to_dict()

                # 计算富集分数（ES、ES_min、ES_max）
                es, es_min, es_max = self._calculate_ssgsea_score(
                    ranked_genes, set(overlap), weights
                )

                # 归一化: NES = ES / (|ES_min| + |ES_max|)
                denominator = abs(es_min) + abs(es_max)
                nes = es / denominator if denominator > 0 else 0.0
                nes_values.append(nes)

            results[tf] = nes_values

        if not results:
            logger.warning("ssGSEA 分析未找到满足大小要求的 TF")
            return pd.DataFrame()

        # 构建 TF x 样本的得分矩阵
        result_df = pd.DataFrame(results, index=samples).T
        result_df.index.name = "TF"

        logger.info("ssGSEA 分析完成: %d 个 TF x %d 个样本", len(results), len(samples))
        return result_df

    def _calculate_ssgsea_score(
        self,
        ranked_genes: List[str],
        gene_set: Set[str],
        gene_weights: Optional[Dict[str, float]] = None
    ) -> Tuple[float, float, float]:
        """计算 ssGSEA 富集分数

        算法步骤：
        1. 确定基因集中与排序列表有交集的基因数（nh）
        2. 计算命中增量（hit_inc）和未命中增量（miss_inc）
        3. 从排序列表顶部开始遍历，维护累积分数：
           - 命中：running_sum += hit_inc * weight
           - 未命中：running_sum -= miss_inc
        4. 记录 ES_max 和 ES_min 用于归一化

        参数:
            ranked_genes: 按某种指标排序的基因列表
            gene_set: 基因集（TF 的靶基因集合）
            gene_weights: 可选的基因权重字典

        返回值:
            Tuple[float, float, float]: (es, es_min, es_max)
        """
        n = len(ranked_genes)
        nh = len(gene_set & set(ranked_genes))

        if nh == 0:
            return 0.0, 0.0, 0.0

        # 计算权重总和
        if gene_weights:
            nr = sum(abs(gene_weights.get(g, 1.0)) for g in gene_set if g in gene_weights)
        else:
            nr = nh

        hit_inc = 1.0 / nr if nr > 0 else 0
        miss_inc = 1.0 / (n - nh) if (n - nh) > 0 else 0

        running_sum = 0.0
        max_es = 0.0
        min_es = 0.0

        for gene in ranked_genes:
            if gene in gene_set:
                weight = abs(gene_weights.get(gene, 1.0)) if gene_weights else 1.0
                running_sum += hit_inc * weight
                if running_sum > max_es:
                    max_es = running_sum
            else:
                running_sum -= miss_inc
                if running_sum < min_es:
                    min_es = running_sum

        return max_es, min_es, max_es

    def gsva(
        self,
        expression_df: pd.DataFrame,
        tf_list: Optional[List[str]] = None,
        method: str = "gsva",
        kcdf: str = "Gaussian",
        tau: float = 1.0,
        min_size: int = 3,
        max_size: int = 500
    ) -> pd.DataFrame:
        """基因集变异分析（GSVA）

        使用 allenricher/core/gsva.py 中的 GSVA 类，将 TF-target 关系作为基因集输入，
        对整个表达矩阵计算 TF 活性得分。

        参数:
            expression_df: 基因 x 样本的表达矩阵，行为基因，列为样本
            tf_list: 可选，指定只分析这些 TF；为 None 时分析所有 TF
            method: GSVA 方法变体，支持 "gsva"（标准）、"plage"、"zscore"
            kcdf: 核密度估计的核函数类型，支持 "Gaussian" 和 "Poisson"
            tau: 核密度估计的带宽参数（仅对 Gaussian 核有效，默认 1.0）
            min_size: 基因集最小大小，低于此值的 TF 将被跳过（默认 3）
            max_size: 基因集最大大小，高于此值的 TF 将被跳过（默认 500）

        返回值:
            DataFrame: TF x 样本的 GSVA 活性得分矩阵
        """
        if expression_df.empty:
            logger.warning("输入表达矩阵为空，返回空 DataFrame")
            return pd.DataFrame()

        # 确定要分析的 TF 列表
        if tf_list is not None:
            tfs_to_analyze = [tf for tf in tf_list if tf in self.tf_to_targets]
        else:
            tfs_to_analyze = list(self.tf_to_targets.keys())

        # 构建基因集字典 {TF_name: set_of_targets}
        gene_sets = {}
        for tf in tfs_to_analyze:
            targets = self.tf_to_targets[tf]
            # 检查与表达矩阵的交集
            overlap = targets & set(expression_df.index)
            if len(overlap) >= min_size and len(overlap) <= max_size:
                gene_sets[tf] = overlap
            else:
                logger.debug(
                    f"跳过 TF '{tf}'：交集基因数 {len(overlap)} 不在 [{min_size}, {max_size}] 范围内"
                )

        if not gene_sets:
            logger.warning("GSVA 分析未找到满足大小要求的 TF")
            return pd.DataFrame()

        # 使用 GSVA 类进行分析
        GSVAClass = self._get_gsva_class()
        gsva_analyzer = GSVAClass(
            method=method,
            kcdf=kcdf,
            tau=tau,
            min_size=min_size,
            max_size=max_size
        )

        result_df = gsva_analyzer.analyze_matrix(expression_df, gene_sets)

        # 重命名索引为 TF
        if not result_df.empty:
            result_df.index.name = "TF"

        logger.info("GSVA 分析完成: %d 个 TF x %d 个样本", len(gene_sets), len(expression_df.columns))
        return result_df

    def get_activators(self, result_df: pd.DataFrame) -> pd.DataFrame:
        """从 ORA 结果中筛选激活型转录因子

        参数:
            result_df: ora() 方法返回的结果 DataFrame

        返回值:
            pd.DataFrame: 仅包含 Mode=='activator' 的 TF 行；
                         如果没有 Mode 列或无匹配则返回空 DataFrame
        """
        if result_df.empty or 'Mode' not in result_df.columns:
            return pd.DataFrame()

        return result_df[result_df['Mode'] == 'activator'].reset_index(drop=True)

    def get_repressors(self, result_df: pd.DataFrame) -> pd.DataFrame:
        """从 ORA 结果中筛选抑制型转录因子

        参数:
            result_df: ora() 方法返回的结果 DataFrame

        返回值:
            pd.DataFrame: 仅包含 Mode=='repressor' 的 TF 行；
                         如果没有 Mode 列或无匹配则返回空 DataFrame
        """
        if result_df.empty or 'Mode' not in result_df.columns:
            return pd.DataFrame()

        return result_df[result_df['Mode'] == 'repressor'].reset_index(drop=True)

    def _build_tf_mode_map(self) -> Dict[str, str]:
        """从 tf_info 构建 TF 调控模式映射

        返回值:
            Dict[str, str]: {TF_name: mode_string}
            - TRRUST: mode 取值 activator/repressor/mixed/unknown
            - ChEA3: 无 mode 列，所有 TF 返回 'unknown'
        """
        mode_map: Dict[str, str] = {}

        if self.tf_info.empty:
            return mode_map

        # TRRUST 格式: TF, mode, target_count
        if 'mode' in self.tf_info.columns:
            tf_col = self.tf_info.columns[0]  # TF 名称列
            for _, row in self.tf_info.iterrows():
                tf_name = str(row[tf_col])
                mode = str(row['mode']).lower()
                mode_map[tf_name] = mode

        # ChEA3 格式: TF, lib_count, target_count（无 mode 列）
        # 所有 TF 默认为 'unknown'

        return mode_map

    @staticmethod
    def _calculate_es(
        ranked_genes: List[str],
        gene_set: Set[str],
        n_genes: int,
        nh: int
    ) -> float:
        """计算基因集的富集分数（Enrichment Score）

        使用与现有 GSEA 相同的 ES 计算方式：
        - hit_inc = 1 / nh（nh 为基因集中与排序列表有交集的基因数）
        - miss_inc = 1 / (n - nh)
        - 遍历排序列表，命中基因集时 += hit_inc，否则 -= miss_inc
        - ES = 遍历过程中累积分数的最大值

        参数:
            ranked_genes: 排好序的基因名列表
            gene_set: 目标基因集合
            n_genes: 排序列表中的基因总数
            nh: 基因集中与排序列表有交集的基因数

        返回值:
            float: 富集分数（ES），正值表示在列表顶部富集
        """
        if nh == 0 or n_genes == nh:
            return 0.0

        hit_inc = 1.0 / nh
        miss_inc = 1.0 / (n_genes - nh)

        running_sum = 0.0
        max_es = 0.0

        for gene in ranked_genes:
            if gene in gene_set:
                running_sum += hit_inc
                if running_sum > max_es:
                    max_es = running_sum
            else:
                running_sum -= miss_inc

        return max_es
