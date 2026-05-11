"""
Core enrichment analysis engine for AllEnricher v2.0

中文模块说明：
    本模块是 AllEnricher v2.0 的核心富集分析引擎，提供了多种富集分析方法，
    包括 Fisher 精确检验、超几何检验和 GSEA（基因集富集分析）。

    主要组件：
    - EnrichmentResult: 富集分析结果的数据容器，封装单个条目的所有统计指标
    - EnrichmentMethodBase: 富集分析方法的抽象基类，定义统一接口
    - FisherExactTest: 基于 Fisher 精确检验的富集分析方法
    - HypergeometricTest: 基于超几何分布检验的富集分析方法
    - GSEA: 基于基因排名的基因集富集分析方法
    - EnrichmentAnalyzer: 主分析引擎，协调整个富集分析流程

    分析流程：
    1. 加载基因列表和背景基因集
    2. 对每个数据库执行富集分析（支持并行/串行）
    3. 多重检验校正（BH/BY/Bonferroni/Holm）
    4. 结果过滤和排序
    5. 输出为 DataFrame 或 TSV 文件

    v1.0 语义兼容：
    analyze_database 方法实现了与 AllEnricher v1.0 相同的两遍扫描逻辑，
    确保 p 值计算结果与旧版本一致。
"""

import os
import csv
import logging
from abc import ABC, abstractmethod
from dataclasses import dataclass
from typing import Dict, List, Optional, Set, Tuple, Any
from pathlib import Path
import math
import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests
from concurrent.futures import ThreadPoolExecutor, as_completed
from tqdm import tqdm

from allenricher.core.config import Config, EnrichmentMethod, CorrectionMethod

# 使用 getLogger 获取模块级 logger，不调用 basicConfig 以避免与 CLI 的日志配置冲突
logger = logging.getLogger(__name__)


def generate_term_url(term_id: str, database: str) -> str:
    """
    生成条目的数据库链接 URL

    根据条目 ID 和数据库名称生成对应的在线数据库链接 URL。
    支持多种主流功能注释数据库，方便用户点击跳转查看详细信息。

    支持的数据库及 URL 格式：
        - GO (Gene Ontology): https://amigo.geneontology.org/amigo/term/{term_id}
        - KEGG: https://www.kegg.jp/entry/{term_id}
        - Reactome: https://reactome.org/PathwayBrowser/#{term_id}
        - WikiPathways: https://www.wikipathways.org/index.php/Pathway:{term_id}
        - DO (Disease Ontology): https://disease-ontology.org/?id={term_id}
        - DisGeNET: https://www.disgenet.org/browser/0/10/{term_id}/

    Args:
        term_id: 条目 ID（如 "GO:0008150"、"hsa00010"、"R-HSA-109582" 等）
        database: 数据库名称（如 "GO"、"KEGG"、"Reactome" 等，不区分大小写）

    Returns:
        str: 条目的数据库链接 URL；如果数据库不被支持则返回空字符串

    Examples:
        >>> generate_term_url("GO:0008150", "GO")
        'https://amigo.geneontology.org/amigo/term/GO:0008150'
        >>> generate_term_url("hsa00010", "KEGG")
        'https://www.kegg.jp/entry/hsa00010'
        >>> generate_term_url("R-HSA-109582", "Reactome")
        'https://reactome.org/PathwayBrowser/#/R-HSA-109582'
    """
    # 将数据库名称统一转换为大写，便于匹配
    database = database.upper()

    # 根据不同的数据库生成对应的 URL
    if database == "GO":
        # Gene Ontology 官方 AmiGO 浏览器
        return f"https://amigo.geneontology.org/amigo/term/{term_id}"
    elif database == "KEGG":
        # KEGG 通路数据库（支持 pathway ID 如 hsa00010）
        return f"https://www.kegg.jp/entry/{term_id}"
    elif database == "REACTOME":
        # Reactome 通路数据库
        return f"https://reactome.org/PathwayBrowser/#{term_id}"
    elif database == "WIKIPATHWAYS":
        # WikiPathways 通路数据库
        return f"https://www.wikipathways.org/index.php/Pathway:{term_id}"
    elif database == "DO":
        # Disease Ontology 疾病本体数据库
        return f"https://disease-ontology.org/?id={term_id}"
    elif database == "DISGENET":
        # DisGeNET 基因-疾病关联数据库
        return f"https://www.disgenet.org/browser/0/10/{term_id}/"
    else:
        # 不支持的数据库返回空字符串
        return ""


@dataclass
class EnrichmentResult:
    """
    富集分析结果的数据容器

    用于封装单个富集条目（term）的所有统计指标和元信息。
    该类使用 dataclass 装饰器自动生成 __init__、__repr__ 等方法。

    属性说明：
        term_id: 条目ID（如 GO:0008150）
        term_name: 条目名称（如 "biological_process"）
        database: 数据库名称（如 "GO_BP"）
        pvalue: 原始 p 值
        adjusted_pvalue: 多重检验校正后的 p 值（如 FDR q 值）
        gene_count: 在该条目中富集到的基因数量
        background_count: 背景基因集中属于该条目的基因数量
        expected_count: 期望基因数量（基于背景比例计算）
        rich_factor: 富集因子（观察值/期望值）
        gene_list: 富集到的基因列表
        gene_ratio: 基因比例字符串（如 "10/500"，格式为 "命中数/基因列表总数"）
        background_ratio: 背景比例字符串（如 "50/20000"，格式为 "背景命中数/背景总数"）
        term_url: 条目的数据库链接 URL（如 https://amigo.geneontology.org/amigo/term/GO:0008150）
        nes: 归一化富集分数（Normalized Enrichment Score，仅 GSEA）
        es: 富集分数（Enrichment Score，仅 GSEA）
        fdr: FDR q 值（仅 GSEA，基于置换检验计算）
        leading_edge: 前沿基因列表（仅 GSEA，对富集分数贡献最大的基因）
    """
    term_id: str  # 条目唯一标识符
    term_name: str  # 条目名称
    database: str  # 来源数据库名称
    pvalue: float  # 原始 p 值（未校正）
    adjusted_pvalue: float  # 多重检验校正后的 p 值
    gene_count: int  # 在该条目中命中的基因数量
    background_count: int  # 背景基因集中属于该条目的基因数量
    expected_count: float  # 期望命中基因数 = (background_count / background_total) * gene_total
    rich_factor: float  # 富集因子 = gene_count / expected_count
    gene_list: List[str]  # 命中的基因名称列表
    gene_ratio: str  # 基因比例，格式 "命中数/基因列表总数"，如 "10/500"
    background_ratio: str  # 背景比例，格式 "背景命中数/背景总数"，如 "50/20000"
    term_url: str = ""  # 条目的数据库链接 URL，默认为空字符串

    # GSEA 特有字段
    nes: Optional[float] = None  # 归一化富集分数（Normalized Enrichment Score）
    es: Optional[float] = None  # 富集分数（Enrichment Score）
    fdr: Optional[float] = None  # FDR q 值（基于置换检验计算）
    leading_edge: Optional[List[str]] = None  # 前沿基因（对 ES 贡献最大的基因子集）
    
    def to_dict(self) -> Dict[str, Any]:
        """
        将富集结果转换为字典格式，便于序列化和导出

        返回值:
            Dict[str, Any]: 包含所有富集统计指标的字典，键名使用大写下划线命名
                           （如 "Term_ID"、"P_Value"），基因列表以分号分隔。
                           GSEA 特有字段（NES、ES、FDR、Leading_Edge）仅在
                           有值时才包含在字典中。
        """
        result = {
            "Term_ID": self.term_id,
            "Term_Name": self.term_name,
            "Database": self.database,
            "P_Value": self.pvalue,
            "Adjusted_P_Value": self.adjusted_pvalue,
            "Gene_Count": self.gene_count,
            "Background_Count": self.background_count,
            "Expected_Count": round(self.expected_count, 4),  # 保留4位小数
            "Rich_Factor": round(self.rich_factor, 4),  # 保留4位小数
            "Gene_Ratio": self.gene_ratio,
            "Background_Ratio": self.background_ratio,
            "Term_URL": self.term_url,  # 条目的数据库链接 URL
            "Genes": ";".join(self.gene_list),  # 基因列表以分号连接为字符串
        }
        
        # 仅在 GSEA 分析时输出以下字段
        if self.nes is not None:
            result["NES"] = round(self.nes, 4)
        if self.es is not None:
            result["ES"] = round(self.es, 4)
        if self.fdr is not None:
            result["FDR"] = round(self.fdr, 4)
        if self.leading_edge is not None:
            result["Leading_Edge"] = ";".join(self.leading_edge)  # 前沿基因以分号连接
            
        return result


class EnrichmentMethodBase(ABC):
    """
    富集分析方法的抽象基类

    定义了所有富集分析方法必须实现的统一接口。
    子类需要实现以下两个抽象方法：
    - calculate_pvalue: 根据统计检验方法计算 p 值
    - calculate_enrichment: 对单个条目执行完整的富集分析

    设计模式：策略模式（Strategy Pattern），使得不同的统计检验方法
    可以在运行时互换，由 EnrichmentAnalyzer 统一调度。
    """
    
    @abstractmethod
    def calculate_pvalue(
        self,
        gene_count: int,
        background_count: int,
        gene_total: int,
        background_total: int
    ) -> float:
        """
        计算富集分析的 p 值（抽象方法，由子类实现）

        参数:
            gene_count: 在该条目中命中的基因数量（即列联表中的 a）
            background_count: 背景基因集中属于该条目的基因数量（即列联表中的 a+c）
            gene_total: 输入基因列表的总数量（即列联表中的 a+b）
            background_total: 背景基因集的总数量（即列联表中的 a+b+c+d）

        返回值:
            float: 计算得到的 p 值，范围 [0, 1]
        """
        pass
    
    @abstractmethod
    def calculate_enrichment(
        self,
        gene_set: Set[str],
        background_set: Set[str],
        term_genes: Set[str],
        term_name: str,
        term_id: str,
        database: str,
        background_total: Optional[int] = None
    ) -> Optional[EnrichmentResult]:
        """
        对单个条目执行完整的富集分析（抽象方法，由子类实现）

        参数:
            gene_set: 输入基因集合（查询基因列表）
            background_set: 背景基因集合
            term_genes: 当前条目（term）所包含的基因集合
            term_name: 条目名称
            term_id: 条目ID
            database: 数据库名称
            background_total: 背景基因总数（可选，默认使用 background_set 的大小；
                            在 v1.0 语义中，该值应为注释文件中所有基因的总数）

        返回值:
            Optional[EnrichmentResult]: 富集分析结果对象；如果该条目不满足分析条件
                                       （如命中基因数为0），则返回 None
        """
        pass


class FisherExactTest(EnrichmentMethodBase):
    """
    Fisher 精确检验富集分析方法

    使用 Fisher 精确检验来判断输入基因列表中某个功能条目（term）的基因
    是否显著富集。该方法基于 2x2 列联表，通过计算超几何分布的精确概率来
    判断观察到的富集程度是否具有统计显著性。

    适用场景：
    - 基因列表较小或期望频数较低时（不满足卡方检验的条件）
    - 需要精确 p 值而非近似 p 值的场景

    统计假设：
    - 零假设（H0）：输入基因列表与该条目之间没有关联
    - 备择假设（H1）：输入基因列表中该条目的基因显著富集（单尾检验，greater）
    """
    
    def calculate_pvalue(
        self,
        gene_count: int,
        background_count: int,
        gene_total: int,
        background_total: int
    ) -> float:
        """
        使用 Fisher 精确检验计算 p 值

        构建 2x2 列联表如下：

                        在条目中    不在条目中    行合计
        在基因列表中      a           b          gene_total
        不在基因列表中    c           d          background_total - gene_total
        列合计        background_count  ...       background_total

        其中：
            a = gene_count          （在基因列表中且在条目中的基因数）
            b = gene_total - a      （在基因列表中但不在条目中的基因数）
            c = background_count - a（不在基因列表中但在条目中的基因数）
            d = 余数                （既不在基因列表中也不在条目中的基因数）

        参数:
            gene_count: 命中基因数（列联表中的 a）
            background_count: 背景中属于该条目的基因数（列联表中的 a+c）
            gene_total: 输入基因列表总数（列联表中的 a+b）
            background_total: 背景基因总数（列联表中的 a+b+c+d）

        返回值:
            float: Fisher 精确检验的单尾 p 值（备择假设为"greater"）
        """
        # 构建 2x2 列联表
        #                    在条目中    不在条目中
        # 在基因列表中       a          b
        # 不在基因列表中     c          d
        
        a = gene_count  # 在基因列表中且在条目中
        b = gene_total - gene_count  # 在基因列表中但不在条目中
        c = background_count - gene_count  # 不在基因列表中但在条目中
        d = background_total - background_count - gene_total + gene_count  # 既不在基因列表中也不在条目中
        
        # 防御性处理：当 background_total 与实际注释范围不一致时，
        # 可能出现负值，将其截断为 0
        if a < 0: a = 0
        if b < 0: b = 0
        if c < 0: c = 0
        if d < 0: d = 0
        
        # 跳过退化列联表（所有值均为0的情况）
        if a + b + c + d == 0:
            return 1.0
        
        # 执行 Fisher 精确检验，使用单尾检验（alternative='greater'），
        # 即检验输入基因列表中该条目的基因是否显著多于随机期望
        _, pvalue = stats.fisher_exact([[a, b], [c, d]], alternative='greater')
        return pvalue
    
    def calculate_enrichment(
        self,
        gene_set: Set[str],
        background_set: Set[str],
        term_genes: Set[str],
        term_name: str,
        term_id: str,
        database: str,
        background_total: Optional[int] = None
    ) -> Optional[EnrichmentResult]:
        """
        使用当前统计检验方法对单个条目执行富集分析（简化版 API）

        .. warning::
            此方法是简化的单条目分析接口，使用标准的 gene_total = len(gene_set)
            语义。如需与 AllEnricher v1.0 完全一致的结果（两遍扫描语义），
            请使用 EnrichmentAnalyzer.analyze_database() 或 run_analysis() 方法。

        参数:
            gene_set: 输入基因集合（查询基因列表）
            background_set: 背景基因集合
            term_genes: 当前条目所包含的基因集合
            term_name: 条目名称
            term_id: 条目ID
            database: 数据库名称
            background_total: 背景基因总数（可选，默认使用 background_set 的大小）

        返回值:
            Optional[EnrichmentResult]: 富集分析结果对象；如果命中基因数 < 1 则返回 None
        """
        # 计算输入基因列表与当前条目的交集
        genes_in_term = gene_set & term_genes
        gene_count = len(genes_in_term)
        
        # 如果命中基因数为0，跳过该条目
        if gene_count < 1:
            return None
        
        # 计算背景基因集中属于该条目的基因数
        background_in_term = background_set & term_genes
        background_count = len(background_in_term)
        
        # 输入基因列表总数
        gene_total = len(gene_set)
        if background_total is None:
            background_total = len(background_set)
        
        # 计算期望命中基因数：期望值 = (条目在背景中的比例) * 输入基因总数
        expected_count = (background_count / background_total) * gene_total if background_total > 0 else 0
        
        # 计算富集因子（Rich Factor）：观察值 / 期望值
        # 值 > 1 表示富集，值 < 1 表示 depleted
        rich_factor = gene_count / expected_count if expected_count > 0 else 0
        
        # 调用 Fisher 精确检验计算 p 值
        pvalue = self.calculate_pvalue(
            gene_count, background_count, gene_total, background_total
        )
        
        # 生成条目的数据库链接 URL
        term_url = generate_term_url(term_id, database)
        
        # 构建并返回富集结果对象
        # 注意：adjusted_pvalue 初始设为 pvalue，后续由 EnrichmentAnalyzer 统一校正
        return EnrichmentResult(
            term_id=term_id,
            term_name=term_name,
            database=database,
            pvalue=pvalue,
            adjusted_pvalue=pvalue,  # 将在后续的多重检验校正步骤中更新
            gene_count=gene_count,
            background_count=background_count,
            expected_count=expected_count,
            rich_factor=rich_factor,
            gene_list=list(genes_in_term),
            gene_ratio=f"{gene_count}/{gene_total}",
            background_ratio=f"{background_count}/{background_total}",
            term_url=term_url  # 条目的数据库链接 URL
        )


class HypergeometricTest(EnrichmentMethodBase):
    """
    超几何检验富集分析方法

    使用超几何分布检验来判断输入基因列表中某个功能条目的基因是否显著富集。
    超几何分布描述了从有限总体中不放回抽样时，抽到特定数量"成功"项的概率。

    与 Fisher 精确检验的关系：
    超几何检验与 Fisher 精确检验在数学上是等价的（两者都基于超几何分布），
    但实现方式不同。Fisher 精确检验通过构建列联表调用 scipy 的 fisher_exact，
    而超几何检验直接使用 scipy 的 hypergeom 生存函数（sf）计算上尾概率。

    超几何分布模型：
    想象一个包含 M 个球的罐子，其中 n 个是红球（属于该条目的基因），
    N 个是白球（不属于该条目的基因）。从中抽取 gene_total 个球（输入基因列表），
    问恰好抽到 k 个红球的概率是多少？

    参数对应关系：
        M (population size) = background_total    （背景基因总数）
        n (successes in pop) = background_count   （背景中属于该条目的基因数）
        N (sample size)      = gene_total         （输入基因列表总数）
        k (successes in sample) = gene_count      （输入基因列表中属于该条目的基因数）
    """
    
    def calculate_pvalue(
        self,
        gene_count: int,
        background_count: int,
        gene_total: int,
        background_total: int
    ) -> float:
        """
        使用超几何检验计算 p 值

        计算公式：
            P(X >= k) = 1 - P(X < k) = 1 - P(X <= k-1)
            使用超几何分布的生存函数（survival function, sf）直接计算上尾概率

        超几何分布参数：
            M = background_total    （总体大小：背景基因总数）
            n = background_count    （总体中的"成功"数：背景中属于该条目的基因数）
            N = gene_total          （样本大小：输入基因列表总数）
            k = gene_count          （样本中的"成功"数：输入基因列表中属于该条目的基因数）

        参数:
            gene_count: 样本中命中的基因数（k）
            background_count: 背景中属于该条目的基因数（n）
            gene_total: 输入基因列表总数（N）
            background_total: 背景基因总数（M）

        返回值:
            float: 超几何检验的 p 值（上尾概率 P(X >= k)）
        """
        # P(X >= k) = 1 - P(X < k) = 1 - P(X <= k-1)
        # 使用生存函数（survival function）直接计算上尾概率
        
        # 超几何分布参数映射
        M = background_total   # 总体大小（背景基因总数）
        n = background_count   # 总体中的成功数（背景中属于该条目的基因数）
        N = gene_total         # 样本大小（输入基因列表总数）
        k = gene_count         # 样本中的成功数（命中基因数）
        
        # sf(k-1, M, n, N) = P(X >= k)，即至少观察到 k 个命中基因的概率
        pvalue = stats.hypergeom.sf(k - 1, M, n, N)
        return pvalue
    
    def calculate_enrichment(
        self,
        gene_set: Set[str],
        background_set: Set[str],
        term_genes: Set[str],
        term_name: str,
        term_id: str,
        database: str,
        background_total: Optional[int] = None
    ) -> Optional[EnrichmentResult]:
        """
        使用超几何检验对单个条目执行富集分析

        实现策略：复用 FisherExactTest 的 calculate_enrichment 方法完成
        基础统计量的计算（交集、期望值、富集因子等），然后用超几何检验
        重新计算 p 值。这样做避免了代码重复。

        参数:
            gene_set: 输入基因集合（查询基因列表）
            background_set: 背景基因集合
            term_genes: 当前条目所包含的基因集合
            term_name: 条目名称
            term_id: 条目ID
            database: 数据库名称
            background_total: 背景基因总数（可选）

        返回值:
            Optional[EnrichmentResult]: 富集分析结果对象；如果不满足条件则返回 None
        """
        # 复用 FisherExactTest 计算基础统计量（交集、期望值、富集因子等）
        fisher = FisherExactTest()
        result = fisher.calculate_enrichment(
            gene_set, background_set, term_genes, term_name, term_id, database,
            background_total=background_total
        )
        
        if result is None:
            return None
        
        # 使用超几何检验重新计算 p 值（替换 Fisher 精确检验的 p 值）
        # 注意：使用 background_total（如果提供）而非 len(background_set)，
        # 以确保与 v1.0 两遍扫描语义一致
        effective_background = background_total if background_total is not None else len(background_set)
        pvalue = self.calculate_pvalue(
            result.gene_count,
            result.background_count,
            len(gene_set),
            effective_background
        )
        result.pvalue = pvalue
        result.adjusted_pvalue = pvalue  # 初始值，后续由 EnrichmentAnalyzer 统一校正
        
        return result


class GSEA(EnrichmentMethodBase):
    """
    基因集富集分析（Gene Set Enrichment Analysis, GSEA）

    GSEA 是一种基于基因排名的富集分析方法，与 ORA（过表示分析，如 Fisher 检验）
    不同，GSEA 不需要对基因进行阈值筛选，而是利用全部基因的表达变化信息。

    核心思想：
    1. 将所有基因按某种指标（如 log2 fold change）排序
    2. 对于每个基因集（term），计算一个富集分数（Enrichment Score, ES）
       ES 反映了该基因集的基因在排序列表中是集中在顶部还是底部
    3. 通过置换检验评估 ES 的统计显著性
    4. 对多个基因集的 ES 进行归一化，得到 NES（Normalized Enrichment Score）

    富集分数（ES）计算方法：
    - 从排序列表的顶部开始遍历
    - 遇到属于该基因集的基因时，增加一个"命中"增量（hit increment）
    - 遇到不属于该基因集的基因时，减少一个"未命中"增量（miss increment）
    - ES = 遍历过程中累积分数的最大值（正向富集）或最小值（负向富集）

    属性:
        permutations: 置换检验的次数（默认 1000）
        min_size: 基因集的最小允许大小（默认 10，与 clusterProfiler 一致）
        max_size: 基因集的最大允许大小（默认 500）
    """
    
    def __init__(self, permutations: int = 1000, min_size: int = 10, max_size: int = 500, seed: int = 42):
        """
        初始化 GSEA 分析器

        参数:
            permutations: 置换检验次数，次数越多 p 值越精确但计算越慢
            min_size: 基因集最小大小，小于此值的基因集将被跳过
            max_size: 基因集最大大小，大于此值的基因集将被跳过
            seed: 随机种子，用于确保置换检验结果的可重复性
        """
        self.permutations = permutations  # 置换次数
        self.min_size = min_size  # 基因集最小大小
        self.max_size = max_size  # 基因集最大大小
        self.seed = seed  # 随机种子
    
    def calculate_pvalue(
        self,
        gene_count: int,
        background_count: int,
        gene_total: int,
        background_total: int
    ) -> float:
        """
        GSEA 使用基于置换检验的 p 值，而非传统的参数检验

        注意：此方法在 GSEA 中不被实际使用，GSEA 的 p 值在 calculate_enrichment
        中通过置换检验计算。此处返回 1.0 作为占位值。

        参数:
            gene_count: 未使用
            background_count: 未使用
            gene_total: 未使用
            background_total: 未使用

        返回值:
            float: 固定返回 1.0（占位值）
        """
        # GSEA 的 p 值通过置换检验在 calculate_enrichment 中计算
        return 1.0
    
    def calculate_enrichment_score(
        self,
        ranked_genes: List[str],
        gene_set: Set[str],
        gene_weights: Optional[Dict[str, float]] = None
    ) -> Tuple[float, List[str]]:
        """
        计算基因集的富集分数（Enrichment Score, ES）

        算法步骤：
        1. 确定基因集中与排序列表有交集的基因数（nh）
        2. 计算命中增量（hit_inc）和未命中增量（miss_inc）
           - hit_inc = 1 / N_hit，其中 N_hit 是基因集中基因的权重之和
           - miss_inc = 1 / (N - nh)，其中 N 是排序列表的总基因数
        3. 从排序列表顶部开始遍历，维护一个累积分数（running_sum）：
           - 遇到属于基因集的基因：running_sum += hit_inc * weight
           - 遇到不属于基因集的基因：running_sum -= miss_inc
        4. ES = 遍历过程中 running_sum 的最大值
        5. 前沿基因（leading edge）= 达到最大 running_sum 时已遍历的基因中，
           属于该基因集的基因

        参数:
            ranked_genes: 按某种指标排序的基因列表（如按 log2 fold change 降序排列）
            gene_set: 当前条目（term）所包含的基因集合
            gene_weights: 可选的基因权重字典（如 log2 fold change 值），
                         如果未提供则使用等权重 1.0

        返回值:
            Tuple[float, List[str]]: 元组，包含两个元素：
                - enrichment_score: 富集分数（ES），正值表示在排序列表顶部富集
                - leading_edge_genes: 前沿基因列表（对 ES 贡献最大的基因）
        """
        n = len(ranked_genes)  # 排序列表中的基因总数
        nh = len(gene_set & set(ranked_genes))  # 基因集中出现在排序列表中的基因数
        
        # 如果没有交集，ES 为 0
        if nh == 0:
            return 0.0, []
        
        # 初始化累积分数和相关变量
        running_sum = 0.0  # 累积富集分数
        max_es = 0.0  # 最大富集分数（即 ES）
        leading_edge = []  # 前沿基因（临时存储）
        
        # 计算命中增量和未命中增量
        # N_hit：基因集中基因的权重绝对值之和（有权重时）或基因数量（无权重时）
        nr = sum(abs(gene_weights.get(g, 1.0)) for g in gene_set if g in gene_weights) if gene_weights else nh
        hit_inc = 1.0 / nr if nr > 0 else 0  # 命中时的增量
        miss_inc = 1.0 / (n - nh) if (n - nh) > 0 else 0  # 未命中时的减量
        
        # 遍历排序列表，计算累积富集分数
        for i, gene in enumerate(ranked_genes):
            if gene in gene_set:
                # 命中：基因属于当前基因集
                weight = abs(gene_weights.get(gene, 1.0)) if gene_weights else 1.0
                running_sum += hit_inc * weight  # 增加加权命中增量
                if running_sum > max_es:
                    max_es = running_sum  # 更新最大富集分数
                    leading_edge = ranked_genes[:i+1]  # 记录当前前沿基因位置
            else:
                # 未命中：基因不属于当前基因集
                running_sum -= miss_inc  # 减少未命中增量
        
        # 从前沿基因中筛选出属于当前基因集的基因
        return max_es, [g for g in leading_edge if g in gene_set]

    def _run_permutation_test(
        self,
        ranked_genes: List[str],
        gene_set: Set[str],
        observed_es: float,
        gene_weights: Optional[Dict[str, float]] = None,
        n_permutations: int = 1000,
        seed: int = 42
    ) -> float:
        """
        通过置换检验计算经验 p 值

        原理：在零假设（基因集与排名无关）下，随机打乱排序列表中基因标签的顺序，
        重新计算富集分数，统计打乱后的 ES 大于等于观察 ES 的次数，以此估计 p 值。

        步骤：
        1. 将 ranked_genes 中的基因标签随机打乱（保持基因集大小不变）
        2. 对打乱后的列表计算 enrichment score
        3. 重复 n_permutations 次
        4. p 值 = (打乱后 ES >= 观察 ES 的次数 + 1) / (n_permutations + 1)
           加 1 是为了避免 p 值为 0（即保证 p 值有下界）

        参数:
            ranked_genes: 排好序的基因列表
            gene_set: 目标基因集
            observed_es: 观察到的富集分数
            gene_weights: 基因权重（可选）
            n_permutations: 置换次数
            seed: 随机种子

        返回值:
            经验 p 值
        """
        rng = np.random.default_rng(seed)  # 使用指定种子初始化随机数生成器
        ranked_array = np.array(ranked_genes)  # 转换为 numpy 数组以便高效打乱
        count_ge = 0  # 统计打乱后 ES >= 观察 ES 的次数

        for _ in range(n_permutations):
            # 随机打乱基因标签顺序（仅打乱排序列表，不改变基因集本身）
            permuted_genes = rng.permutation(ranked_array).tolist()
            # 计算打乱后的富集分数（只取 ES，忽略 leading_edge）
            permuted_es, _ = self.calculate_enrichment_score(
                permuted_genes, gene_set, gene_weights
            )
            # 统计打乱后 ES 大于等于观察 ES 的次数
            if permuted_es >= observed_es:
                count_ge += 1

        # 经验 p 值公式：(null_ES >= observed_ES 的次数 + 1) / (n_permutations + 1)
        pvalue = (count_ge + 1) / (n_permutations + 1)
        return pvalue

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
    ) -> Optional[EnrichmentResult]:
        """
        使用 GSEA 方法对单个条目执行富集分析

        分析流程：
        1. 检查基因集大小是否在允许范围内 [min_size, max_size]
        2. 计算富集分数（ES）和前沿基因
        3. 计算归一化富集分数（NES）= ES * sqrt(N/nh)
        4. 通过置换检验计算 p 值和 FDR（当前为占位实现）

        参数:
            gene_set: 输入基因集合（查询基因列表）
            background_set: 背景基因集合
            term_genes: 当前条目所包含的基因集合
            term_name: 条目名称
            term_id: 条目ID
            database: 数据库名称
            ranked_genes: 可选的排序基因列表（如按 fold change 排序）；
                         如果未提供，则使用 background_set 作为排序列表
            gene_weights: 可选的基因权重字典

        返回值:
            Optional[EnrichmentResult]: GSEA 富集分析结果对象；如果基因集大小
                                       不满足要求则返回 None
        """
        # 检查基因集大小是否在允许范围内
        overlap = gene_set & term_genes
        if len(overlap) < self.min_size or len(overlap) > self.max_size:
            return None
        
        # 如果未提供排序列表，使用背景基因集作为默认排序列表
        if ranked_genes is None:
            ranked_genes = list(background_set)
        
        # 计算富集分数（ES）和前沿基因
        es, leading_edge = self.calculate_enrichment_score(
            ranked_genes, term_genes, gene_weights
        )
        
        # 计算归一化富集分数（NES）
        # NES = ES * sqrt(N / nh)，其中 N 为排序列表基因总数，nh 为交集基因数
        # 归一化消除了基因集大小对 ES 的影响，使得不同大小的基因集可以比较
        n = len(ranked_genes)
        nh = len(overlap)
        nes = es * np.sqrt(n / nh) if nh > 0 else 0
        
        # 通过置换检验计算经验 p 值
        # 将排序列表中的基因标签随机打乱，重新计算 ES，统计打乱后 ES >= 观察 ES 的次数
        pvalue = self._run_permutation_test(
            ranked_genes=ranked_genes,
            gene_set=term_genes,
            observed_es=es,
            gene_weights=gene_weights,
            n_permutations=self.permutations,
            seed=self.seed
        )
        # FDR 暂时设为 pvalue，后续在 adjust_pvalues 中统一进行 BH 校正
        fdr = pvalue
        
        # 生成条目的数据库链接 URL
        term_url = generate_term_url(term_id, database)
        
        # 构建并返回 GSEA 富集结果对象
        return EnrichmentResult(
            term_id=term_id,
            term_name=term_name,
            database=database,
            pvalue=pvalue,
            adjusted_pvalue=pvalue,
            gene_count=len(overlap),
            background_count=len(term_genes),
            expected_count=0,  # GSEA 不使用期望值概念
            rich_factor=0,  # GSEA 不使用富集因子概念
            gene_list=list(overlap),
            gene_ratio=f"{len(overlap)}/{len(gene_set)}",
            background_ratio=f"{len(term_genes)}/{len(background_set)}",
            term_url=term_url,  # 条目的数据库链接 URL
            nes=nes,  # 归一化富集分数
            es=es,  # 原始富集分数
            fdr=fdr,  # FDR q 值
            leading_edge=leading_edge  # 前沿基因
        )


class SSGSEA(EnrichmentMethodBase):
    """
    单样本 GSEA（Single-Sample GSEA, ssGSEA）算法

    ssGSEA 是 GSEA 的单样本变体，由 Barbie 等人于 2009 年提出。
    与标准 GSEA 的核心区别在于：
    - 不需要排列检验（permutation test），直接计算归一化富集分数（NES）
    - 每个样本独立计算，适用于单样本分析场景（如单细胞测序、肿瘤样本分型等）
    - 归一化方式不同：NES = ES / (|ES_min| + |ES_max|)，范围在 [-1, 1] 之间

    算法步骤：
    1. 将基因按表达量（或其他指标）降序排列
    2. 对每个基因集，计算 running enrichment score（与 GSEA 相同的累积逻辑）
    3. 在计算过程中跟踪累积分数的最大值（ES_max）和最小值（ES_min）
    4. 归一化：NES = ES / (|ES_min| + |ES_max|)

    属性:
        min_size: 基因集的最小允许大小（默认 10，与 clusterProfiler 一致）
        max_size: 基因集的最大允许大小（默认 500）
    """

    def __init__(self, min_size: int = 10, max_size: int = 500):
        """
        初始化 ssGSEA 分析器

        参数:
            min_size: 基因集最小大小，小于此值的基因集将被跳过（默认 10，与 clusterProfiler 一致）
            max_size: 基因集最大大小，大于此值的基因集将被跳过（默认 500）
        """
        self.method_name = "ssgsea"  # 方法名称
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
        ssGSEA 不使用传统的 p 值计算

        注意：ssGSEA 不做统计检验，因此此方法返回 NaN 作为占位值。
        ssGSEA 的核心统计量是 NES（归一化富集分数），而非 p 值。

        参数:
            gene_count: 未使用
            background_count: 未使用
            gene_total: 未使用
            background_total: 未使用

        返回值:
            float: 固定返回 NaN（不适用）
        """
        # ssGSEA 不进行统计检验，p 值不适用
        return float('nan')

    def calculate_enrichment_score(
        self,
        ranked_genes: List[str],
        gene_set: Set[str],
        gene_weights: Optional[Dict[str, float]] = None
    ) -> Tuple[float, float, float, List[str]]:
        """
        计算基因集的富集分数（Enrichment Score, ES）及归一化所需的最小/最大累积分数

        算法步骤：
        1. 确定基因集中与排序列表有交集的基因数（nh）
        2. 计算命中增量（hit_inc）和未命中增量（miss_inc）
           - hit_inc = 1 / N_hit，其中 N_hit 是基因集中基因的权重绝对值之和
           - miss_inc = 1 / (N - nh)，其中 N 是排序列表的总基因数
        3. 从排序列表顶部开始遍历，维护一个累积分数（running_sum）：
           - 遇到属于基因集的基因：running_sum += hit_inc * weight
           - 遇到不属于基因集的基因：running_sum -= miss_inc
        4. ES = 遍历过程中 running_sum 的最大值（正向富集）
        5. 同时记录 ES_min 和 ES_max，用于后续归一化

        参数:
            ranked_genes: 按某种指标排序的基因列表（如按表达量降序排列）
            gene_set: 当前条目（term）所包含的基因集合
            gene_weights: 可选的基因权重字典（如表达量值），
                         如果未提供则使用等权重 1.0

        返回值:
            Tuple[float, float, float, List[str]]: 元组，包含四个元素：
                - es: 富集分数（Enrichment Score），正值表示在排序列表顶部富集
                - es_min: 遍历过程中的最小累积分数
                - es_max: 遍历过程中的最大累积分数
                - leading_edge_genes: 前沿基因列表（对 ES 贡献最大的基因）
        """
        n = len(ranked_genes)  # 排序列表中的基因总数
        nh = len(gene_set & set(ranked_genes))  # 基因集中出现在排序列表中的基因数

        # 如果没有交集，ES 为 0，ES_min 和 ES_max 也为 0
        if nh == 0:
            return 0.0, 0.0, 0.0, []

        # 初始化累积分数和相关变量
        running_sum = 0.0  # 累积富集分数
        max_es = 0.0  # 最大累积分数（即 ES）
        min_es = 0.0  # 最小累积分数（用于归一化）
        leading_edge = []  # 前沿基因（临时存储）

        # 计算命中增量和未命中增量
        # N_hit：基因集中基因的权重绝对值之和（有权重时）或基因数量（无权重时）
        nr = sum(abs(gene_weights.get(g, 1.0)) for g in gene_set if g in gene_weights) if gene_weights else nh
        hit_inc = 1.0 / nr if nr > 0 else 0  # 命中时的增量
        miss_inc = 1.0 / (n - nh) if (n - nh) > 0 else 0  # 未命中时的减量

        # 遍历排序列表，计算累积富集分数
        for i, gene in enumerate(ranked_genes):
            if gene in gene_set:
                # 命中：基因属于当前基因集
                weight = abs(gene_weights.get(gene, 1.0)) if gene_weights else 1.0
                running_sum += hit_inc * weight  # 增加加权命中增量
                if running_sum > max_es:
                    max_es = running_sum  # 更新最大累积分数
                    leading_edge = ranked_genes[:i+1]  # 记录当前前沿基因位置
            else:
                # 未命中：基因不属于当前基因集
                running_sum -= miss_inc  # 减少未命中增量
                if running_sum < min_es:
                    min_es = running_sum  # 更新最小累积分数

        # 从前沿基因中筛选出属于当前基因集的基因
        leading_edge_genes = [g for g in leading_edge if g in gene_set]
        return max_es, min_es, max_es, leading_edge_genes

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
    ) -> Optional[EnrichmentResult]:
        """
        使用 ssGSEA 方法对单个条目执行富集分析

        与标准 GSEA 的区别：
        - 不需要排列检验，直接计算归一化富集分数（NES）
        - NES = ES / (|ES_min| + |ES_max|)，范围在 [-1, 1] 之间
        - p 值设为 NaN（不适用），NES 作为主要统计量

        分析流程：
        1. 检查基因集大小是否在允许范围内 [min_size, max_size]
        2. 计算富集分数（ES）、ES_min、ES_max 和前沿基因
        3. 归一化：NES = ES / (|ES_min| + |ES_max|)
        4. 构建结果对象（p 值为 NaN，FDR 为 NaN）

        参数:
            gene_set: 输入基因集合（查询基因列表）
            background_set: 背景基因集合
            term_genes: 当前条目所包含的基因集合
            term_name: 条目名称
            term_id: 条目ID
            database: 数据库名称
            ranked_genes: 可选的排序基因列表（如按表达量排序）；
                         如果未提供，则使用 background_set 作为排序列表
            gene_weights: 可选的基因权重字典（如基因表达量）

        返回值:
            Optional[EnrichmentResult]: ssGSEA 富集分析结果对象；如果基因集大小
                                       不满足要求则返回 None
        """
        # 检查基因集大小是否在允许范围内
        overlap = gene_set & term_genes
        if len(overlap) < self.min_size or len(overlap) > self.max_size:
            return None

        # 如果未提供排序列表，使用背景基因集作为默认排序列表
        if ranked_genes is None:
            ranked_genes = list(background_set)

        # 计算富集分数（ES）、ES_min、ES_max 和前沿基因
        es, es_min, es_max, leading_edge = self.calculate_enrichment_score(
            ranked_genes, term_genes, gene_weights
        )

        # 计算归一化富集分数（NES）
        # NES = ES / (|ES_min| + |ES_max|)
        # 分母为累积分数的绝对值范围，确保 NES 在 [-1, 1] 之间
        denominator = abs(es_min) + abs(es_max)
        nes = es / denominator if denominator > 0 else 0.0

        # ssGSEA 不做统计检验，p 值和 FDR 均设为 NaN
        pvalue = float('nan')
        fdr = float('nan')

        # 生成条目的数据库链接 URL
        term_url = generate_term_url(term_id, database)

        # 构建并返回 ssGSEA 富集结果对象
        return EnrichmentResult(
            term_id=term_id,
            term_name=term_name,
            database=database,
            pvalue=pvalue,  # NaN：ssGSEA 不做统计检验
            adjusted_pvalue=float('nan'),  # NaN：无校正后 p 值
            gene_count=len(overlap),
            background_count=len(term_genes),
            expected_count=0,  # ssGSEA 不使用期望值概念
            rich_factor=0,  # ssGSEA 不使用富集因子概念
            gene_list=list(overlap),
            gene_ratio=f"{len(overlap)}/{len(gene_set)}",
            background_ratio=f"{len(term_genes)}/{len(background_set)}",
            term_url=term_url,  # 条目的数据库链接 URL
            nes=nes,  # 归一化富集分数（核心统计量，范围 [-1, 1]）
            es=es,  # 原始（未归一化）富集分数
            fdr=fdr,  # NaN：ssGSEA 不计算 FDR
            leading_edge=leading_edge  # 前沿基因
        )


class EnrichmentAnalyzer:
    """
    富集分析主引擎

    协调整个富集分析流程的核心类，负责：
    1. 根据配置选择合适的富集分析方法（Fisher / 超几何 / GSEA）
    2. 加载基因列表和背景基因集
    3. 对每个数据库执行富集分析（支持并行和串行两种模式）
    4. 多重检验校正（BH / BY / Bonferroni / Holm）
    5. 结果过滤、排序和导出

    使用示例：
        config = Config(method="fisher", correction="BH")
        analyzer = EnrichmentAnalyzer(config)
        results = analyzer.run_analysis(gene_set, background_set, database_data)
    """
    
    def __init__(self, config: Config):
        """
        初始化富集分析引擎

        参数:
            config: 配置对象，包含分析方法、校正方法、过滤阈值等参数
        """
        self.config = config  # 全局配置对象
        self._method: Optional[EnrichmentMethodBase] = None  # 延迟初始化的分析方法实例
        self.results: Dict[str, List[EnrichmentResult]] = {}  # 存储各数据库的分析结果
    
    @property
    def method(self) -> EnrichmentMethodBase:
        """
        延迟获取富集分析方法实例（属性描述符）

        首次访问时根据配置创建对应的方法实例并缓存，
        后续访问直接返回缓存值。这样允许在构造后修改 config.method
        而不会导致立即报错。

        返回值:
            EnrichmentMethodBase: 富集分析方法实例
        """
        if self._method is None:
            self._method = self._get_method()
        return self._method
    
    @method.setter
    def method(self, value: EnrichmentMethodBase):
        """允许直接设置方法实例"""
        self._method = value
    
    def _get_method(self) -> EnrichmentMethodBase:
        """
        根据配置创建对应的富集分析方法实例

        返回值:
            EnrichmentMethodBase: 富集分析方法实例（FisherExactTest / HypergeometricTest / GSEA）

        异常:
            ValueError: 当配置中指定了未知的方法名称时抛出
        """
        # 方法名称到实例的映射表
        methods = {
            EnrichmentMethod.FISHER.value: FisherExactTest(),  # Fisher 精确检验
            EnrichmentMethod.HYPERGEOMETRIC.value: HypergeometricTest(),  # 超几何检验
            EnrichmentMethod.GSEA.value: GSEA(  # GSEA 基因集富集分析
                permutations=self.config.gsea_permutations,
                min_size=self.config.gsea_min_size,
                max_size=self.config.gsea_max_size
            ),
            EnrichmentMethod.SSGSEA.value: SSGSEA(  # ssGSEA 单样本基因集富集分析
                min_size=self.config.gsea_min_size,
                max_size=self.config.gsea_max_size
            ),
        }
        
        if self.config.method not in methods:
            raise ValueError(f"Unknown method: {self.config.method}")
        
        return methods[self.config.method]
    
    def load_gene_list(self, file_path: str) -> Set[str]:
        """
        加载基因列表文件，支持多种输入格式

        支持的格式：
        1. 每行一个基因（.txt/.tsv/.gene，原有格式）
        2. CSV 格式（.csv，逗号分隔），读取第一列
        3. Excel 格式（.xlsx/.xls），读取第一列
        4. 空格/逗号/分号分隔的单行格式（自动检测）
        5. 自动检测分隔符：如果第一行包含逗号，按 CSV 解析；否则按行解析

        所有格式均会：
        - 跳过注释行（以 # 开头）
        - 跳过空行
        - 自动去除首尾空白字符
        - 自动去重（使用 Set 存储）

        参数:
            file_path: 基因列表文件的路径

        返回值:
            Set[str]: 去重后的基因名称集合

        异常:
            FileNotFoundError: 文件不存在时抛出
            ValueError: 文件格式无法识别或文件内容为空时抛出
        """
        file_path = Path(file_path)
        
        # 检查文件是否存在
        if not file_path.exists():
            raise FileNotFoundError(f"基因列表文件不存在: {file_path}")
        
        # 获取文件扩展名（统一转为小写）
        suffix = file_path.suffix.lower()
        
        genes = set()
        
        # --- 根据文件扩展名选择解析方式 ---
        if suffix in ('.xlsx', '.xls'):
            # Excel 格式：使用 pandas 读取第一列
            genes = self._load_from_excel(file_path)
        elif suffix == '.csv':
            # CSV 格式：使用 csv 模块，读取第一列
            genes = self._load_from_csv(file_path)
        else:
            # .txt / .tsv / .gene / 无扩展名 / 其他：尝试自动检测格式
            genes = self._load_from_text_auto(file_path)
        
        # 检查最终结果是否为空
        if not genes:
            raise ValueError(
                f"基因列表文件 {file_path} 中未找到有效的基因名称。"
                f"请检查文件内容是否正确（支持 # 开头的注释行）。"
            )
        
        logger.info(f"从 {file_path} 加载了 {len(genes)} 个基因")
        return genes
    
    def _load_from_excel(self, file_path: Path) -> Set[str]:
        """
        从 Excel 文件加载基因列表（读取第一列）

        支持的格式：.xlsx, .xls
        读取第一个工作表的第一列，跳过空值和注释行。

        参数:
            file_path: Excel 文件路径

        返回值:
            Set[str]: 基因名称集合
        """
        genes = set()
        try:
            # 使用 pandas 读取 Excel 文件的第一个工作表
            df = pd.read_excel(file_path, engine='openpyxl', header=None)
            # 读取第一列的所有行
            for item in df.iloc[:, 0].dropna():
                gene = str(item).strip()
                if gene and not gene.startswith('#'):
                    genes.add(gene)
        except ImportError:
            # 如果 openpyxl 未安装，给出明确的安装提示
            raise ImportError(
                "读取 Excel 文件需要 openpyxl 库。请执行: pip install openpyxl"
            )
        except Exception as e:
            raise ValueError(f"读取 Excel 文件 {file_path} 失败: {e}")
        
        return genes
    
    def _load_from_csv(self, file_path: Path) -> Set[str]:
        """
        从 CSV 文件加载基因列表（读取第一列）

        使用 Python csv 模块解析，支持带表头和不带表头的 CSV 文件。
        读取第一列的所有非空值，跳过注释行。

        参数:
            file_path: CSV 文件路径

        返回值:
            Set[str]: 基因名称集合
        """
        genes = set()
        try:
            with open(file_path, 'r', encoding='utf-8') as f:
                reader = csv.reader(f)
                for row in reader:
                    if not row:
                        continue  # 跳过空行
                    gene = row[0].strip()  # 读取第一列
                    if gene and not gene.startswith('#'):
                        genes.add(gene)
        except Exception as e:
            raise ValueError(f"读取 CSV 文件 {file_path} 失败: {e}")
        
        return genes
    
    def _load_from_text_auto(self, file_path: Path) -> Set[str]:
        """
        从文本文件加载基因列表（自动检测格式）

        自动检测策略：
        1. 读取文件全部内容
        2. 如果只有一行且包含逗号/分号/空格分隔符 → 按分隔符拆分
        3. 如果第一行包含逗号且多行 → 按 CSV 格式解析（取第一列）
        4. 否则 → 按每行一个基因的方式解析（原有格式）

        所有格式均跳过注释行（# 开头）和空行。

        参数:
            file_path: 文本文件路径

        返回值:
            Set[str]: 基因名称集合
        """
        genes = set()
        
        with open(file_path, 'r', encoding='utf-8') as f:
            lines = f.readlines()
        
        # 过滤掉注释行和空行，保留有效行
        valid_lines = [line.strip() for line in lines if line.strip() and not line.strip().startswith('#')]
        
        if not valid_lines:
            return genes  # 返回空集合，由调用方处理
        
        # 情况1：只有一行，尝试按分隔符拆分（空格/逗号/分号/制表符）
        if len(valid_lines) == 1:
            line = valid_lines[0]
            # 按优先级尝试不同的分隔符：逗号 > 分号 > 制表符 > 空格
            for sep in [',', ';', '\t', ' ']:
                if sep in line:
                    parts = [p.strip() for p in line.split(sep) if p.strip()]
                    if len(parts) > 1:
                        # 成功按分隔符拆分为多个基因
                        genes.update(parts)
                        return genes
            # 只有一个基因（无分隔符），直接添加
            genes.add(line)
            return genes
        
        # 情况2：多行，检查第一行是否包含逗号（可能是 CSV 格式）
        first_line = valid_lines[0]
        if ',' in first_line:
            # 可能是 CSV 格式，使用 csv 模块解析第一列
            try:
                import io
                content = ''.join(lines)  # 重新拼接原始内容（保留注释行信息由 csv 处理）
                reader = csv.reader(io.StringIO(content))
                for row in reader:
                    if not row:
                        continue
                    gene = row[0].strip()
                    if gene and not gene.startswith('#'):
                        genes.add(gene)
                # 如果解析结果合理（多于1个基因），使用 CSV 解析结果
                if len(genes) > 0:
                    return genes
            except Exception:
                pass  # CSV 解析失败，回退到按行解析
        
        # 情况3：默认按每行一个基因解析（原有格式）
        for line in valid_lines:
            gene = line.strip()
            if gene and not gene.startswith('#'):
                genes.add(gene)
        
        return genes
    
    def adjust_pvalues(
        self,
        results: List[EnrichmentResult],
        method: str = "BH"
    ) -> List[EnrichmentResult]:
        """
        对富集分析结果进行多重检验校正

        由于富集分析通常同时检验数千个条目，需要进行多重检验校正
        以控制假阳性率。本方法支持以下校正方法：

        校正方法说明：
        - BH (Benjamini-Hochberg): 控制 FDR（错误发现率），
          是最常用的校正方法，在控制假阳性的同时保持较高的统计功效
        - BY (Benjamini-Yekutieli): 比 BH 更保守的 FDR 控制方法，
          适用于条目之间存在相关性的情况
        - Bonferroni: 最保守的校正方法，控制家族错误率（FWER），
          将 p 值乘以检验次数，容易产生假阴性
        - Holm (Holm-Bonferroni): Bonferroni 的改进版本，逐步校正，
          比 Bonferroni 略有更高的统计功效
        - NONE: 不进行校正，直接使用原始 p 值

        参数:
            results: 富集分析结果列表
            method: 校正方法名称，默认为 "BH"

        返回值:
            List[EnrichmentResult]: 校正后的结果列表（每个结果的 adjusted_pvalue 已更新）
        """
        if not results:
            return results
        
        # 提取所有结果的原始 p 值
        pvalues = [r.pvalue for r in results]
        
        # 检查是否包含 NaN p 值（如 ssGSEA 结果），如果有则跳过校正
        # 因为多重检验校正算法无法正确处理 NaN 值
        nan_count = sum(1 for p in pvalues if math.isnan(p))
        if nan_count > 0:
            logger.warning(f"检测到 {nan_count}/{len(pvalues)} 个 NaN p 值（可能来自 ssGSEA），跳过多重检验校正")
            return results
        
        # 校正方法名称映射：配置枚举值 -> statsmodels 中的方法标识符
        correction_methods = {
            CorrectionMethod.BH.value: "fdr_bh",  # Benjamini-Hochberg FDR 校正
            CorrectionMethod.BY.value: "fdr_by",  # Benjamini-Yekutieli FDR 校正
            CorrectionMethod.BONFERRONI.value: "bonferroni",  # Bonferroni 校正
            CorrectionMethod.HOLM.value: "holm",  # Holm-Bonferroni 逐步校正
            CorrectionMethod.NONE.value: None,  # 不校正
        }
        
        # 如果选择不校正，直接返回原始结果
        if method == CorrectionMethod.NONE.value:
            return results
        
        # 获取对应的校正方法标识符，默认使用 BH
        corr_method = correction_methods.get(method, "fdr_bh")
        
        # 调用 statsmodels 的 multipletests 进行多重检验校正
        # 返回值：reject (是否拒绝H0), adjusted_pvalues, _, _
        _, adjusted, _, _ = multipletests(pvalues, method=corr_method)
        
        # 将校正后的 p 值写回每个结果对象
        for result, adj_p in zip(results, adjusted):
            result.adjusted_pvalue = adj_p
        
        return results
    
    def filter_results(
        self,
        results: List[EnrichmentResult]
    ) -> List[EnrichmentResult]:
        """
        根据配置参数过滤富集分析结果

        过滤条件：
        1. 校正后 p 值（adjusted_pvalue）<= qvalue_cutoff
        2. 命中基因数（gene_count）>= min_genes
        3. 命中基因数（gene_count）<= max_genes

        参数:
            results: 待过滤的富集分析结果列表

        返回值:
            List[EnrichmentResult]: 过滤后的结果列表
        """
        filtered = []
        
        for result in results:
            # 检查校正后 p 值是否超过阈值
            if result.adjusted_pvalue > self.config.qvalue_cutoff:
                continue
            
            # 检查命中基因数是否低于最小值
            if result.gene_count < self.config.min_genes:
                continue
            
            # 检查命中基因数是否超过最大值（max_genes 为 inf 时表示无限制）
            if self.config.max_genes != float('inf') and result.gene_count > self.config.max_genes:
                continue
            
            filtered.append(result)
        
        return filtered
    
    def analyze_database(
        self,
        gene_set: Set[str],
        background_set: Set[str],
        term_data: Dict[str, Dict[str, Any]],
        database: str
    ) -> List[EnrichmentResult]:
        """
        对单个数据库执行富集分析（v1.0 语义兼容）

        本方法实现了与 AllEnricher v1.0 完全一致的两遍扫描逻辑：

        === v1.0 语义说明 ===
        - background_total: 注释文件中所有基因的总数（即 gene2term 矩阵的总行数），
          而非用户提供的背景基因集大小
        - gene_total: 输入基因列表中至少命中一个条目的基因数（在第一遍扫描完成后确定），
          而非输入基因列表的总数
        - num_in_C (background_count): 条目基因与用户背景基因集的交集大小
        - 仅保留 num_in_O > 1 的条目（即至少有2个输入基因命中）

        === 两遍扫描逻辑 ===
        第一遍扫描：
        - 遍历所有条目，计算每个条目的基础统计量
        - 收集所有至少命中一个条目的输入基因（genes_with_hits）
        - 跳过 num_in_O <= 1 的条目
        - 最终 gene_total = len(genes_with_hits)

        第二遍扫描：
        - 使用第一遍确定的 gene_total 计算每个条目的 p 值
        - 计算期望值和富集因子
        - 构建 EnrichmentResult 对象

        注意：v1.0 中 ExpectedGeneNum/RichFactor 使用动态增长的 gene_total，
        但 p 值计算使用最终的 gene_total。本实现统一使用最终 gene_total
        以保持结果的自洽性。

        参数:
            gene_set: 输入基因集合（查询基因列表）
            background_set: 背景基因集合（用户提供的背景基因）
            term_data: 条目数据字典，格式为 {term_id: {"name": str, "genes": List[str]}}
            database: 数据库名称

        返回值:
            List[EnrichmentResult]: 该数据库的富集分析结果列表
        """
        results = []
        
        # --- 输入验证 ---
        # 检查 gene_set 是否为空
        if not gene_set:
            logger.error("输入基因列表（gene_set）为空，无法执行富集分析。请检查输入文件是否有效。")
            return results
        
        # 检查 background_set 是否为空（给出警告但不中断，因为某些场景下可能不需要背景集）
        if not background_set:
            logger.warning("背景基因集（background_set）为空，分析结果可能不可靠。建议提供背景基因集。")
        
        # 检查 term_data 是否为空
        if not term_data:
            logger.warning(f"数据库 {database} 的条目数据（term_data）为空，跳过该数据库。")
            return results
        
        # --- v1.0 语义：background_total = 注释文件中所有基因的总数 ---
        # 收集注释文件中出现过的所有基因（即 gene2term 矩阵中所有行）
        all_annotated_genes: Set[str] = set()
        for term_info in term_data.values():
            all_annotated_genes.update(term_info.get("genes", []))
        background_total = len(all_annotated_genes)
        
        if background_total == 0:
            logger.warning(f"No annotated genes found for {database}")
            return results
        
        # --- 第一遍扫描：计算每个条目的基础统计量，并确定最终的 gene_total ---
        term_stats = {}  # term_id -> {num_in_O, num_in_C, genes_in_term, term_name}
        genes_with_hits: Set[str] = set()  # 至少命中一个条目的输入基因
        
        for term_id, term_info in term_data.items():
            term_genes = set(term_info.get("genes", []))
            genes_in_term = gene_set & term_genes  # 输入基因与条目基因的交集
            num_in_O = len(genes_in_term)  # 观察到的命中基因数
            
            # v1.0 语义：仅保留 num_in_O > 1 的条目（至少2个输入基因命中）
            if num_in_O <= 1:
                continue
            
            # 记录命中的基因（用于后续计算 gene_total）
            genes_with_hits.update(genes_in_term)
            # 计算背景中属于该条目的基因数（条目基因与用户背景基因集的交集）
            num_in_C = len(term_genes & background_set)
            
            term_stats[term_id] = {
                "num_in_O": num_in_O,  # 观察到的命中基因数
                "num_in_C": num_in_C,  # 背景中属于该条目的基因数
                "genes_in_term": genes_in_term,  # 命中的基因集合
                "term_name": term_info.get("name", term_id),  # 条目名称
            }
        
        # 最终 gene_total（v1.0 语义：等价于 R 代码中 %gene_list1 在完整循环后的值）
        gene_total = len(genes_with_hits)
        
        if gene_total == 0:
            return results
        
        # --- 第二遍扫描：使用最终 gene_total 计算 p 值和其他统计量 ---
        for term_id, ts in tqdm(
            term_stats.items(),
            desc=f"Analyzing {database}",
            leave=False
        ):
            try:
                num_in_O = ts["num_in_O"]  # 观察到的命中基因数
                num_in_C = ts["num_in_C"]  # 背景中属于该条目的基因数
                
                # 使用选定的统计检验方法计算 p 值
                pvalue = self.method.calculate_pvalue(
                    gene_count=num_in_O,
                    background_count=num_in_C,
                    gene_total=gene_total,  # 使用第一遍扫描确定的最终值
                    background_total=background_total  # 使用注释文件中的基因总数
                )
                
                # 计算期望命中基因数 = (条目在背景中的比例) * 输入基因总数
                expected = num_in_C / background_total * gene_total if background_total > 0 else 0
                # 计算富集因子 = 观察值 / 期望值
                rich_factor = num_in_O / expected if expected > 0 else 0
                
                # 生成条目的数据库链接 URL
                term_url = generate_term_url(term_id, database)
                
                # 构建富集结果对象
                result = EnrichmentResult(
                    term_id=term_id,
                    term_name=ts["term_name"],
                    database=database,
                    pvalue=pvalue,
                    adjusted_pvalue=pvalue,  # 初始值，后续由 adjust_pvalues 统一校正
                    gene_count=num_in_O,
                    background_count=num_in_C,
                    expected_count=round(expected, 6),
                    rich_factor=round(rich_factor, 6),
                    gene_list=sorted(ts["genes_in_term"]),  # 基因列表按字母排序
                    gene_ratio=f"{num_in_O}/{gene_total}",
                    background_ratio=f"{num_in_C}/{background_total}",
                    term_url=term_url  # 条目的数据库链接 URL
                )
                results.append(result)
            except Exception as e:
                # 单个条目计算失败不应中断整个分析，记录错误并跳过该条目
                logger.warning(f"计算条目 {term_id}（{database}）时出错，已跳过：{e}")
                continue
        
        return results
    
    def run_analysis(
        self,
        gene_set: Set[str],
        background_set: Set[str],
        database_data: Dict[str, Dict[str, Dict[str, Any]]],
        parallel: bool = True
    ) -> Dict[str, pd.DataFrame]:
        """
        对所有数据库执行富集分析

        分析流程：
        1. 对每个数据库调用 analyze_database 进行富集分析
        2. 对每个数据库的结果进行多重检验校正
        3. 根据配置参数过滤结果
        4. 按校正后 p 值排序
        5. 将结果转换为 DataFrame 格式

        支持两种执行模式：
        - 并行模式（parallel=True 且 n_jobs > 1）：
          使用 ThreadPoolExecutor 并行处理多个数据库。
          注意：由于 EnrichmentAnalyzer 实例持有不可 pickle 的对象（如 Config），
          因此使用线程池而非进程池。对于 CPU 密集型的统计计算，
          线程池在 CPython 中受 GIL 限制，但数据库 I/O 和
          结果处理仍可从并行中受益。
        - 串行模式（parallel=False 或 n_jobs = 1）：
          逐个处理数据库

        参数:
            gene_set: 输入基因集合（查询基因列表）
            background_set: 背景基因集合
            database_data: 数据库数据字典，格式为：
                          {database_name: {term_id: {"name": str, "genes": List[str]}}}
            parallel: 是否使用并行处理，默认为 True

        返回值:
            Dict[str, pd.DataFrame]: 各数据库的富集分析结果，
                                    键为数据库名称，值为结果 DataFrame
        """
        self.results = {}
        
        if parallel and self.config.n_jobs > 1 and len(database_data) > 1:
            # 并行处理模式：使用线程池并行分析多个数据库
            # 注意：使用 ThreadPoolExecutor 而非 ProcessPoolExecutor，
            # 因为 EnrichmentAnalyzer 对象包含不可 pickle 的成员，
            # 无法在进程间序列化传递。
            logger.info(f"使用并行模式分析 {len(database_data)} 个数据库（线程数: {self.config.n_jobs}）")
            with ThreadPoolExecutor(max_workers=self.config.n_jobs) as executor:
                # 提交所有数据库的分析任务
                futures = {
                    executor.submit(
                        self.analyze_database,
                        gene_set,
                        background_set,
                        term_data,
                        database
                    ): database
                    for database, term_data in database_data.items()
                }
                
                # 收集并行任务的结果
                for future in as_completed(futures):
                    database = futures[future]
                    try:
                        results = future.result()
                        results = self.adjust_pvalues(results, self.config.correction)  # 多重检验校正
                        results = self.filter_results(results)  # 过滤结果
                        results.sort(key=lambda x: x.adjusted_pvalue)  # 按校正后 p 值排序
                        self.results[database] = results
                        logger.info(f"Completed {database}: {len(results)} enriched terms")
                    except Exception as e:
                        logger.error(f"Error analyzing {database}: {e}")
        else:
            # 串行处理模式：逐个分析数据库
            for database, term_data in database_data.items():
                try:
                    results = self.analyze_database(
                        gene_set, background_set, term_data, database
                    )
                    results = self.adjust_pvalues(results, self.config.correction)  # 多重检验校正
                    results = self.filter_results(results)  # 过滤结果
                    results.sort(key=lambda x: x.adjusted_pvalue)  # 按校正后 p 值排序
                    self.results[database] = results
                    logger.info(f"Completed {database}: {len(results)} enriched terms")
                except Exception as e:
                    logger.error(f"Error analyzing {database}: {e}")
        
        return self.to_dataframes()
    
    def to_dataframes(self) -> Dict[str, pd.DataFrame]:
        """
        将分析结果转换为 pandas DataFrame 格式

        返回值:
            Dict[str, pd.DataFrame]: 各数据库的结果 DataFrame，
                                    键为数据库名称，值为包含所有富集指标的 DataFrame
        """
        dataframes = {}
        
        for database, results in self.results.items():
            if results:
                # 将每个 EnrichmentResult 转换为字典，然后构建 DataFrame
                data = [r.to_dict() for r in results]
                df = pd.DataFrame(data)
                dataframes[database] = df
        
        return dataframes
    
    def save_results(self, output_dir: str) -> None:
        """
        将分析结果保存为 TSV 文件

        每个数据库的结果保存为单独的 TSV 文件，文件名格式为：
        {database_name}_enrichment.tsv

        参数:
            output_dir: 输出目录路径，如果不存在会自动创建
        """
        output_path = Path(output_dir)
        output_path.mkdir(parents=True, exist_ok=True)  # 自动创建输出目录
        
        for database, results in self.results.items():
            if results:
                # 保存为 TSV 格式（制表符分隔）
                output_file = output_path / f"{database}_enrichment.tsv"
                data = [r.to_dict() for r in results]
                df = pd.DataFrame(data)
                df.to_csv(output_file, sep='\t', index=False)
                logger.info(f"Saved {database} results to {output_file}")
