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

from __future__ import annotations

import os
import csv
import logging
import collections
from abc import ABC, abstractmethod
from dataclasses import dataclass, field
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

    属性说明：
        term_id: 条目ID（如 GO:0008150）
        term_name: 条目名称（如 "biological_process"）
        database: 数据库名称（如 "GO_BP"）
        pvalue: 原始 p 值
        adjusted_pvalue: 多重检验校正后的 p 值（如 FDR q 值）
        gene_count: 在该条目中富集到的基因数量
        background_count: 背景基因集中属于该条目的基因数量（仅 ORA）
        expected_count: 期望基因数量（仅 ORA）
        rich_factor: 富集因子（仅 ORA）
        gene_list: 富集到的基因列表
        gene_ratio: 基因比例字符串（仅 ORA）
        background_ratio: 背景比例字符串（仅 ORA）
        term_url: 条目的数据库链接 URL
        nes: 归一化富集分数（Normalized Enrichment Score，仅 GSEA）
        es: 富集分数（Enrichment Score，仅 GSEA）
        fdr: FDR q 值（仅 GSEA）
        leading_edge: 前沿基因列表（仅 GSEA）
        set_size: 基因集大小（经过 min/max 过滤后，仅 GSEA）
        rank_at_max: 达到最大 ES 时的排位（仅 GSEA）
    """
    term_id: str
    term_name: str
    database: str
    pvalue: float
    adjusted_pvalue: float
    gene_count: int
    background_count: int = 0
    expected_count: float = 0.0
    rich_factor: float = 0.0
    gene_list: List[str] = field(default_factory=list)
    gene_ratio: str = ""
    background_ratio: str = ""
    term_url: str = ""

    # GSEA 特有字段
    nes: Optional[float] = None
    es: Optional[float] = None
    fdr: Optional[float] = None
    leading_edge: Optional[List[str]] = None
    set_size: Optional[int] = None
    rank_at_max: Optional[int] = None
    fwerp: Optional[float] = None          # gseapy: FWER p-val
    tag_pct: str = ""                       # gseapy: Tag %（如 "6/37"）
    gene_pct: str = ""                      # gseapy: Gene %（如 "30.00%"）
    
    def to_dict(self) -> Dict[str, Any]:
        """
        将富集结果转换为字典格式，便于序列化和导出

        返回值:
            Dict[str, Any]: 包含所有富集统计指标的字典。
                           对于 GSEA 结果（nes 不为 None），
                           列名对齐 gseapy.res2d 标准输出：
                           Term, Description, setSize, ES, NES,
                           NOM p-val, FDR q-val, FWER p-val,
                           rank, Tag %, Gene %, Lead_genes,
                           core_enrichment
                           对于 ORA 结果，输出全部 ORA 字段。
        """
        is_gsea = self.nes is not None
        
        result = collections.OrderedDict() if is_gsea else {}
        
        if is_gsea:
            result["Term_ID"] = self.term_id
            result["Term_Name"] = self.term_name
            result["Database"] = self.database
            result["setSize"] = self.set_size if self.set_size is not None else self.gene_count
            result["ES"] = round(self.es, 4) if self.es is not None else 0.0
            result["NES"] = round(self.nes, 4) if self.nes is not None else 0.0
            result["p_value"] = self.pvalue
            result["FDR"] = self.fdr if self.fdr is not None else self.adjusted_pvalue
            result["rank"] = self.rank_at_max if self.rank_at_max is not None else 0
            result["Tag %"] = self.tag_pct if self.tag_pct else ""
            result["Gene %"] = self.gene_pct if self.gene_pct else ""
            result["Lead_genes"] = ";".join(self.leading_edge) if self.leading_edge else ""
            result["matched_genes"] = ";".join(self.gene_list)
            result["Term_URL"] = self.term_url
        else:
            result["Term_ID"] = self.term_id
            result["Term_Name"] = self.term_name
            result["Database"] = self.database
            result["P_Value"] = self.pvalue
            result["Adjusted_P_Value"] = self.adjusted_pvalue
            result["Gene_Count"] = self.gene_count
            result["Background_Count"] = self.background_count
            result["Expected_Count"] = round(self.expected_count, 4)
            result["Rich_Factor"] = round(self.rich_factor, 4)
            result["Gene_Ratio"] = self.gene_ratio
            result["Background_Ratio"] = self.background_ratio
            result["Term_URL"] = self.term_url
            result["Genes"] = ";".join(self.gene_list)
            
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
    ) -> Tuple[float, List[str], int]:
        """
        计算基因集的富集分数（Enrichment Score, ES）

        算法步骤：
        1. 确定基因集中与排序列表有交集的基因数（nh）
        2. 计算命中增量（hit_inc）和未命中增量（miss_inc）
        3. 从排序列表顶部开始遍历，维护一个累积分数（running_sum）
        4. ES = 遍历过程中 running_sum 的最大绝对值（正值）
        5. 前沿基因（leading edge）= 达到最大 ES 时已遍历的基因中，属于该基因集的基因

        参数:
            ranked_genes: 按某种指标排序的基因列表
            gene_set: 当前条目（term）所包含的基因集合
            gene_weights: 可选的基因权重字典

        返回值:
            Tuple[float, List[str], int]: (ES, leading_edge_genes, rank_at_max)
                - enrichment_score: 富集分数（ES）
                - leading_edge_genes: 前沿基因列表
                - rank_at_max: 达到最大 ES 时排序列表中的排位（1-based）
        """
        n = len(ranked_genes)
        nh = len(gene_set & set(ranked_genes))

        if nh == 0:
            return 0.0, [], 0

        running_sum = 0.0
        max_es = 0.0
        rank_at_max = 0
        leading_edge_indices = []

        nr = sum(abs(gene_weights.get(g, 1.0)) for g in gene_set if g in gene_weights) if gene_weights else nh
        hit_inc = 1.0 / nr if nr > 0 else 0
        miss_inc = 1.0 / (n - nh) if (n - nh) > 0 else 0

        for i, gene in enumerate(ranked_genes):
            if gene in gene_set:
                weight = abs(gene_weights.get(gene, 1.0)) if gene_weights else 1.0
                running_sum += hit_inc * weight
                if running_sum > max_es:
                    max_es = running_sum
                    rank_at_max = i + 1  # 1-based rank
                    leading_edge_indices = list(range(i + 1))
            else:
                running_sum -= miss_inc

        leading_edge = [ranked_genes[idx] for idx in leading_edge_indices if ranked_genes[idx] in gene_set]
        return max_es, leading_edge, rank_at_max

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
            permuted_es, _, _ = self.calculate_enrichment_score(
                permuted_genes, gene_set, gene_weights
            )
            # 统计打乱后 ES 大于等于观察 ES 的次数
            if permuted_es >= observed_es:
                count_ge += 1

        # 经验 p 值公式：(null_ES >= observed_ES 的次数 + 1) / (n_permutations + 1)
        pvalue = (count_ge + 1) / (n_permutations + 1)
        return pvalue

    def calculate_normalized_es(
        self,
        ranked_genes: List[str],
        gene_set: Set[str],
        gene_weights: Optional[Dict[str, float]] = None
    ) -> Tuple[float, float, float, List[str], int]:
        """
        计算基于置换检验的归一化富集分数 (NES)

        步骤:
        1. 计算实际 ES、leading edge 和 rank_at_max
        2. 进行 self.permutations 次置换，每次打乱基因标签重新计算 ES
        3. 分别计算正向和负向 null ES 分布的均值
        4. NES = ES / mean(|ES_null|)（正负分别归一化）
        5. pvalue = 置换检验中 |ES_null| >= |ES| 的比例

        参数:
            ranked_genes: 按某种指标排序的基因列表
            gene_set: 目标基因集合
            gene_weights: 可选的基因权重字典

        返回值:
            Tuple[float, float, float, List[str], int]: (ES, NES, pvalue, leading_edge_genes, rank_at_max)
        """
        # 步骤 1: 计算实际 ES、前沿基因和 rank_at_max
        es, leading_edge, rank_at_max = self.calculate_enrichment_score(
            ranked_genes, gene_set, gene_weights
        )

        # 空基因集或无交集的情况
        if es == 0.0:
            return 0.0, 0.0, 1.0, [], 0

        # 步骤 2: 进行置换检验
        rng = np.random.default_rng(self.seed)
        ranked_array = np.array(ranked_genes)
        null_es_positive = []  # 正向 null ES
        null_es_negative = []  # 负向 null ES
        count_ge = 0  # 统计 |ES_null| >= |ES| 的次数

        for _ in range(self.permutations):
            permuted_genes = rng.permutation(ranked_array).tolist()
            permuted_es, _, _ = self.calculate_enrichment_score(
                permuted_genes, gene_set, gene_weights
            )

            # 分别收集正负向 null ES
            if permuted_es >= 0:
                null_es_positive.append(permuted_es)
            else:
                null_es_negative.append(permuted_es)

            # 统计 |ES_null| >= |ES| 的次数
            if abs(permuted_es) >= abs(es):
                count_ge += 1

        # 步骤 3: 计算正负向 null ES 分布的均值
        mean_pos = np.mean(null_es_positive) if null_es_positive else 1.0
        mean_neg = np.mean([abs(e) for e in null_es_negative]) if null_es_negative else 1.0

        # 步骤 4: NES = ES / mean(|ES_null|)（正负分别归一化）
        if es > 0:
            nes = es / mean_pos if mean_pos > 0 else 0.0
        else:
            nes = es / mean_neg if mean_neg > 0 else 0.0

        # 步骤 5: pvalue = 置换检验中 |ES_null| >= |ES| 的比例
        pvalue = (count_ge + 1) / (self.permutations + 1)

        return es, nes, pvalue, leading_edge, rank_at_max

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
        2. 通过基于置换检验的归一化方法计算 ES、NES、pvalue 和前沿基因
        3. 构建富集分析结果对象

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
        set_size = len(overlap)
        if set_size < self.min_size or set_size > self.max_size:
            return None

        # 如果未提供排序列表，使用背景基因集作为默认排序列表
        if ranked_genes is None:
            ranked_genes = list(background_set)

        # 使用基于置换检验的归一化方法计算 ES、NES、pvalue、前沿基因和 rank_at_max
        es, nes, pvalue, leading_edge, rank_at_max = self.calculate_normalized_es(
            ranked_genes, term_genes, gene_weights
        )

        # FDR 初始设为 pvalue，adjust_pvalues() 中会统一进行 BH 校正
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
            gene_count=set_size,
            gene_list=list(overlap),
            term_url=term_url,
            nes=nes,
            es=es,
            fdr=fdr,
            leading_edge=leading_edge,
            set_size=set_size,
            rank_at_max=rank_at_max,
        )

    def analyze_matrix(
        self,
        expression_matrix: pd.DataFrame,
        gene_sets: Dict[str, Set[str]],
        gene_weights_matrix: Optional[pd.DataFrame] = None
    ) -> pd.DataFrame:
        """
        基于表达矩阵运行 GSEA 分析

        对每个样本:
        1. 按表达量排序基因
        2. 对每个基因集计算 ES 和 NES
        3. 汇总为样本 x 通路结果矩阵

        参数:
            expression_matrix: 表达矩阵 (行=基因, 列=样本)
            gene_sets: {通路名: 基因集合}
            gene_weights_matrix: 可选的权重矩阵 (行=基因, 列=样本)

        返回值:
            DataFrame (行=通路, 列=样本, 值=NES)
        """
        samples = expression_matrix.columns.tolist()
        results = {}

        for pathway_name, pathway_genes in gene_sets.items():
            nes_values = []
            for sample in samples:
                # 获取当前样本的表达量并按降序排序基因
                sample_expr = expression_matrix[sample]
                ranked_genes = sample_expr.sort_values(ascending=False).index.tolist()

                # 获取可选的权重
                weights = None
                if gene_weights_matrix is not None and sample in gene_weights_matrix.columns:
                    weights = gene_weights_matrix[sample].to_dict()

                # 计算归一化 ES
                _, nes, _, _ = self.calculate_normalized_es(
                    ranked_genes, pathway_genes, weights
                )
                nes_values.append(nes)

            results[pathway_name] = nes_values

        return pd.DataFrame(results, index=samples).T


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

    def analyze_matrix(
        self,
        expression_matrix: pd.DataFrame,
        gene_sets: Dict[str, Set[str]]
    ) -> pd.DataFrame:
        """
        基于表达矩阵运行 ssGSEA 分析

        对每个样本独立计算每个通路的 ssGSEA 得分。

        参数:
            expression_matrix: 表达矩阵 (行=基因, 列=样本)
            gene_sets: {通路名: 基因集合}

        返回值:
            DataFrame (行=通路, 列=样本, 值=NES)
        """
        samples = expression_matrix.columns.tolist()
        results = {}

        for pathway_name, pathway_genes in gene_sets.items():
            nes_values = []
            for sample in samples:
                # 获取当前样本的表达量并按降序排序基因
                sample_expr = expression_matrix[sample]
                ranked_genes = sample_expr.sort_values(ascending=False).index.tolist()

                # 将表达量作为权重
                weights = sample_expr.to_dict()

                # 计算富集分数（ES、ES_min、ES_max、前沿基因）
                es, es_min, es_max, _ = self.calculate_enrichment_score(
                    ranked_genes, pathway_genes, weights
                )

                # 归一化: NES = ES / (|ES_min| + |ES_max|)
                denominator = abs(es_min) + abs(es_max)
                nes = es / denominator if denominator > 0 else 0.0
                nes_values.append(nes)

            results[pathway_name] = nes_values

        return pd.DataFrame(results, index=samples).T


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
            EnrichmentMethodBase: 富集分析方法实例（HypergeometricTest / GSEA / SSGSEA / GSVA）

        异常:
            ValueError: 当配置中指定了未知的方法名称时抛出
        """
        # 方法名称到实例的映射表
        methods = {
            EnrichmentMethod.HYPERGEOMETRIC.value: HypergeometricTest(),  # 超几何检验（ORA默认方法）
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

        # 延迟导入 GSVA 以避免循环依赖（gsva.py 导入了 enrichment.py 的基类）
        from allenricher.core.gsva import GSVA
        methods[EnrichmentMethod.GSVA.value] = GSVA(  # GSVA 基因集变异分析
            method=self.config.gsva_method,
            kcdf=self.config.gsva_kcdf,
            tau=self.config.gsva_tau,
            min_size=self.config.gsea_min_size,
            max_size=self.config.gsea_max_size
        )
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
    
    def load_ranked_gene_list(self, file_path: str) -> List[Tuple[str, float]]:
        """
        加载排序基因列表文件（GSEA 专用）

        文件格式要求：TSV格式，包含以下列：
        - gene: 基因名称（必需）
        - weight: 排序权重值（必需，如 log2 fold change）
        - rank: 可选，排序位置

        参数:
            file_path: 排序基因列表文件的路径

        返回值:
            List[Tuple[str, float]]: (基因名称, 权重值) 元组列表，按 rank 升序排列

        异常:
            FileNotFoundError: 文件不存在时抛出
            ValueError: 文件格式不正确时抛出
        """
        file_path_obj = Path(file_path)
        if not file_path_obj.exists():
            raise FileNotFoundError(f"排序基因列表文件不存在: {file_path}")
        
        df = pd.read_csv(file_path_obj, sep='\t')
        
        # 检查必需列
        if 'gene' not in df.columns:
            raise ValueError(
                f"排序基因列表文件必须包含 'gene' 列。"
                f"当前列名: {list(df.columns)}"
            )
        
        # 检查 weight 列，如果不存在则使用默认值 1.0
        if 'weight' not in df.columns:
            logger.warning("排序基因列表文件未包含 'weight' 列，默认使用权重 1.0")
            df['weight'] = 1.0
        
        # 按 rank 列排序（如果存在）
        if 'rank' in df.columns:
            df = df.sort_values('rank')
        
        # 构建 (基因, 权重) 元组列表
        result = list(zip(df['gene'].tolist(), df['weight'].tolist()))
        
        logger.info(f"从 {file_path} 加载了 {len(result)} 个排序基因")
        return result
    
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
        # 先检测是否为 TSV 格式（含制表符的多列数据）
        tab_count = sum(1 for line in valid_lines if '\t' in line)
        is_tsv = tab_count > 0 and tab_count >= len(valid_lines) * 0.5  # 超过50%的行含制表符
        
        if is_tsv:
            # TSV 格式：按制表符拆分，取第一列；跳过表头行
            for i, line in enumerate(valid_lines):
                parts = line.split('\t')
                gene = parts[0].strip()
                if not gene:
                    continue
                # 跳过表头行（第0行且值为常见列名）
                if i == 0 and gene.lower() in ('gene', 'symbol', 'id', 'gene_id', 'gene_symbol', 'name', 'genes'):
                    continue
                genes.add(gene)
        else:
            # 普通每行一个基因格式
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
            # GSEA 的 FDR 字段也应同步为校正后值
            if result.fdr is not None:
                result.fdr = adj_p
        
        return results
    
    def filter_results(
        self,
        results: List[EnrichmentResult]
    ) -> List[EnrichmentResult]:
        """
        根据配置参数过滤富集分析结果

        过滤条件（取决于 output_all 配置）：
        - output_all=True（默认，与v1一致）: 仅过滤不满足 min_genes/max_genes 的条目，保留全部 p 值
        - output_all=False: 额外过滤 adjusted_pvalue > qvalue_cutoff 的条目

        参数:
            results: 待过滤的富集分析结果列表

        返回值:
            List[EnrichmentResult]: 过滤后的结果列表
        """
        filtered = []

        for result in results:
            # 仅当 output_all=False 时才按 q 值过滤
            if not self.config.output_all:
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
        database: str,
        ranked_gene_list: Optional[List[Tuple[str, float]]] = None
    ) -> List[EnrichmentResult]:
        """
        对单个数据库执行富集分析

        对于 ORA 方法（Fisher / Hypergeometric），使用 v1.0 语义兼容的两遍扫描逻辑。
        对于 GSEA/ssGSEA/GSVA 方法，使用各方法自身的 calculate_enrichment 计算逻辑。

        v1.0 ORA 语义说明：
        - background_total: 注释文件中所有基因的总数（即 gene2term 矩阵的总行数），
          而非用户提供的背景基因集大小
        - gene_total: 输入基因列表中至少命中一个条目的基因数（在第一遍扫描完成后确定），
          而非输入基因列表的总数
        - 仅保留 num_in_O > 1 的条目

        === 两遍扫描逻辑（ORA 方法）===
        第一遍扫描：遍历所有条目，收集基础统计量，确定 gene_total
        第二遍扫描：使用最终 gene_total 计算 p 值和其他统计量

        参数:
            gene_set: 输入基因集合（查询基因列表）
            background_set: 背景基因集合（用户提供的背景基因）
            term_data: 条目数据字典，格式为 {term_id: {"name": str, "genes": List[str]}}
            database: 数据库名称
            ranked_gene_list: 可选的排序基因列表（GSEA/ssGSEA 使用）

        返回值:
            List[EnrichmentResult]: 该数据库的富集分析结果列表
        """
        # 对于 GSEA/ssGSEA/GSVA 方法，使用各方法自身的计算逻辑
        if self.config.method in (EnrichmentMethod.GSEA.value, EnrichmentMethod.SSGSEA.value):
            return self._analyze_gsea_like(
                gene_set, background_set, term_data, database,
                ranked_gene_list=ranked_gene_list
            )
        elif self.config.method == EnrichmentMethod.GSVA.value:
            return self._analyze_gsva(
                gene_set, background_set, term_data, database
            )
        
        # === ORA 方法（Fisher / Hypergeometric）原有逻辑 ===
        results = []
        
        # --- 输入验证 ---
        if not gene_set:
            logger.error("输入基因列表（gene_set）为空，无法执行富集分析。请检查输入文件是否有效。")
            return results
        
        if not background_set:
            logger.warning("背景基因集（background_set）为空，分析结果可能不可靠。建议提供背景基因集。")
        
        if not term_data:
            logger.warning(f"数据库 {database} 的条目数据（term_data）为空，跳过该数据库。")
            return results
        
        # --- background_total = 背景基因集的总大小 ---
        background_total = len(background_set) if background_set else 0
        
        if background_total == 0:
            logger.warning(f"背景基因集为空，无法执行富集分析。")
            return results
        
        # --- 第一遍扫描：计算每个条目的基础统计量，并确定最终的 gene_total ---
        term_stats = {}
        genes_with_hits: Set[str] = set()
        
        for term_id, term_info in term_data.items():
            term_genes = set(term_info.get("genes", []))
            genes_in_term = gene_set & term_genes
            num_in_O = len(genes_in_term)
            
            if num_in_O <= 1:
                continue
            
            genes_with_hits.update(genes_in_term)
            # 背景中属于该条目的基因数（使用 background_set 过滤）
            background_in_term = background_set & term_genes if background_set else term_genes
            num_in_C = len(background_in_term)
            
            term_stats[term_id] = {
                "num_in_O": num_in_O,
                "num_in_C": num_in_C,
                "genes_in_term": genes_in_term,
                "term_name": term_info.get("name", term_id),
            }
        
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
                num_in_O = ts["num_in_O"]
                num_in_C = ts["num_in_C"]
                
                pvalue = self.method.calculate_pvalue(
                    gene_count=num_in_O,
                    background_count=num_in_C,
                    gene_total=gene_total,
                    background_total=background_total
                )
                
                expected = num_in_C / background_total * gene_total if background_total > 0 else 0
                rich_factor = num_in_O / expected if expected > 0 else 0
                
                term_url = generate_term_url(term_id, database)
                
                result = EnrichmentResult(
                    term_id=term_id,
                    term_name=ts["term_name"],
                    database=database,
                    pvalue=pvalue,
                    adjusted_pvalue=pvalue,
                    gene_count=num_in_O,
                    background_count=num_in_C,
                    expected_count=round(expected, 6),
                    rich_factor=round(rich_factor, 6),
                    gene_list=sorted(ts["genes_in_term"]),
                    gene_ratio=f"{num_in_O}/{gene_total}",
                    background_ratio=f"{num_in_C}/{background_total}",
                    term_url=term_url
                )
                results.append(result)
            except Exception as e:
                logger.warning(f"计算条目 {term_id}（{database}）时出错，已跳过：{e}")
                continue
        
        return results
    
    def _analyze_gsea_like(
        self,
        gene_set: Set[str],
        background_set: Set[str],
        term_data: Dict[str, Dict[str, Any]],
        database: str,
        ranked_gene_list: Optional[List[Tuple[str, float]]] = None
    ) -> List[EnrichmentResult]:
        """
        使用 GSEA/ssGSEA 方法对单个数据库执行富集分析

        - GSEA: 使用 gseapy.prerank（批量分析所有通路，含置换检验 + FDR 校正）
        - ssGSEA: 沿用原有的逐条目循环

        参数:
            gene_set: 输入基因集合
            background_set: 背景基因集合
            term_data: 条目数据字典
            database: 数据库名称
            ranked_gene_list: 可选的 (基因, 权重) 元组列表

        返回值:
            List[EnrichmentResult]: 该数据库的富集分析结果列表
        """
        if not term_data:
            logger.warning(f"数据库 {database} 的条目数据为空，跳过。")
            return []
        
        # 从 ranked_gene_list 中提取基因名称列表和权重字典
        gene_names = None
        gene_weights = None
        if ranked_gene_list:
            gene_names = [g for g, w in ranked_gene_list]
            gene_weights = {g: w for g, w in ranked_gene_list}
        
        # ---------------------------------------------------------------
        # GSEA：使用 gseapy.prerank（准确、经过验证的行业标准实现）
        # ---------------------------------------------------------------
        if self.config.method == EnrichmentMethod.GSEA.value:
            try:
                return self._run_gsea_with_gseapy(
                    term_data, database, gene_names, gene_weights
                )
            except ImportError:
                logger.warning("gseapy 未安装，回退至内置 GSEA 实现。请运行 pip install gseapy")
            except Exception as e:
                logger.warning(f"gseapy 执行失败 ({e})，回退至内置 GSEA 实现。")
        
        # ---------------------------------------------------------------
        # ssGSEA / 回退：原有逐条目循环
        # ---------------------------------------------------------------
        results = []
        for term_id, term_info in tqdm(
            term_data.items(),
            desc=f"Analyzing {database}",
            leave=False
        ):
            try:
                term_genes = set(term_info.get("genes", []))
                
                # 直接使用方法实例的 calculate_enrichment
                result = self.method.calculate_enrichment(
                    gene_set=gene_set,
                    background_set=background_set,
                    term_genes=term_genes,
                    term_name=term_info.get("name", term_id),
                    term_id=term_id,
                    database=database,
                    ranked_genes=gene_names,
                    gene_weights=gene_weights
                )
                
                if result is not None:
                    results.append(result)
            except Exception as e:
                logger.warning(f"计算条目 {term_id}（{database}）时出错，已跳过：{e}")
                continue
        
        return results
    
    def _run_gsea_with_gseapy(
        self,
        term_data: Dict[str, Dict[str, Any]],
        database: str,
        gene_names: Optional[List[str]],
        gene_weights: Optional[Dict[str, float]],
    ) -> List[EnrichmentResult]:
        """
        使用 gseapy.prerank 对单个数据库执行 GSEA 分析

        将数据库的所有基因集一次性传给 gseapy，批量完成 ES/NES/pvalue/FDR 计算，
        然后将结果映射为 EnrichmentResult 列表。

        参数:
            term_data: 条目数据字典 {term_id: {name, genes}}
            database: 数据库名称
            gene_names: 排序基因列表（基因名称）
            gene_weights: 基因权重字典

        返回值:
            List[EnrichmentResult]: 富集分析结果列表
        """
        import gseapy
        import pandas as pd
        import numpy as np

        # 构建 gseapy 所需的基因集字典: {term_id: [gene1, gene2, ...]}
        gene_sets_dict: Dict[str, List[str]] = {}
        term_names: Dict[str, str] = {}
        for term_id, term_info in term_data.items():
            term_genes = list(term_info.get("genes", []))
            if term_genes:
                gene_sets_dict[term_id] = term_genes
                term_names[term_id] = term_info.get("name", term_id)

        if not gene_sets_dict:
            logger.warning(f"数据库 {database} 无有效基因集数据，跳过。")
            return []

        # 构建排序列表为 pandas Series
        if gene_names and gene_weights:
            # 保留原始符号权重
            rnk = pd.Series(
                data=[gene_weights.get(g, 0.0) for g in gene_names],
                index=gene_names,
            )
        elif gene_names:
            rnk = pd.Series(data=np.ones(len(gene_names)), index=gene_names)
        else:
            logger.warning("未提供排序基因列表，跳过 GSEA 分析。")
            return []

        # 获取 GSEA 参数（从 self.method 中读取）
        gsea_min_size = getattr(self.method, 'min_size', 10)
        gsea_max_size = getattr(self.method, 'max_size', 500)
        gsea_permutations = getattr(self.method, 'permutations', 1000)
        gsea_seed = getattr(self.method, 'seed', 42)

        logger.info(
            f"gseapy.prerank: {len(gene_sets_dict)} gene sets, "
            f"{len(rnk)} ranked genes, "
            f"{gsea_permutations} permutations"
        )

        # 执行 gseapy.prerank（行业标准 GSEA 实现）
        gs_res = gseapy.prerank(
            rnk=rnk,
            gene_sets=gene_sets_dict,
            min_size=gsea_min_size,
            max_size=gsea_max_size,
            permutation_num=gsea_permutations,
            no_plot=True,
            verbose=False,
            seed=gsea_seed,
            threads=min(4, getattr(self.config, 'n_jobs', 1)),
            ascending=False,
        )

        # 映射 gseapy 结果 → EnrichmentResult
        # gs_res.results = {term_id: {es, nes, pval, fdr, lead_genes, matched_genes, tag%, gene%, ...}}
        results: List[EnrichmentResult] = []
        n_ranked = len(gene_names) if gene_names else 0
        for term_id, gs_term in gs_res.results.items():
            try:
                # 解析前沿基因
                lead_genes_str = gs_term.get('lead_genes', '')
                lead_genes = lead_genes_str.split(';') if lead_genes_str else []

                # 解析匹配基因
                matched_genes_str = gs_term.get('matched_genes', '')
                matched_genes = matched_genes_str.split(';') if matched_genes_str else []

                # 从 tag% 解析 setSize（格式: "lead_count/setSize"）
                tag = str(gs_term.get('tag %', '0/0'))
                tag_parts = tag.split('/')
                set_size = int(tag_parts[1]) if len(tag_parts) > 1 else len(matched_genes)

                # 从 gene% 计算 rank_at_max
                gene_pct_str = str(gs_term.get('gene %', '0%')).replace('%', '')
                try:
                    gene_pct = float(gene_pct_str) / 100.0
                    rank_at_max = int(gene_pct * n_ranked)
                except (ValueError, ZeroDivisionError):
                    rank_at_max = 0

                # gseapy 的 fdr 已基于置换检验 + BH 校正
                fdr = float(gs_term['fdr'])
                pval = float(gs_term['pval'])
                nes = float(gs_term['nes'])
                es = float(gs_term['es'])
                fwerp = float(gs_term.get('fwerp', 0.0))

                result = EnrichmentResult(
                    term_id=term_id,
                    term_name=term_names.get(term_id, term_id),
                    database=database,
                    pvalue=pval,
                    adjusted_pvalue=fdr,  # gseapy 直接提供 BH 校正后的 FDR
                    gene_count=set_size,
                    gene_list=matched_genes,
                    term_url=generate_term_url(term_id, database),
                    nes=nes,
                    es=es,
                    fdr=fdr,
                    leading_edge=lead_genes,
                    set_size=set_size,
                    rank_at_max=rank_at_max,
                    fwerp=fwerp,
                    tag_pct=tag,
                    gene_pct=str(gs_term.get('gene %', '0%')),
                )
                results.append(result)
            except Exception as e:
                logger.warning(f"解析 gseapy 结果 {term_id} 时出错: {e}")
                continue

        logger.info(
            f"gseapy 完成 {database}: {len(results)} enriched terms"
        )
        return results
    
    def _analyze_gsva(
        self,
        gene_set: Set[str],
        background_set: Set[str],
        term_data: Dict[str, Dict[str, Any]],
        database: str
    ) -> List[EnrichmentResult]:
        """
        使用 GSVA 方法对单个数据库执行富集分析

        GSVA 的核心输出是通路活性矩阵（通过 analyze_matrix），
        但为保持接口统一，此处对每个条目调用 calculate_enrichment 返回占位结果。

        参数:
            gene_set: 输入基因集合
            background_set: 背景基因集合（未使用，为保持接口一致）
            term_data: 条目数据字典
            database: 数据库名称

        返回值:
            List[EnrichmentResult]: 该数据库的富集分析结果列表
        """
        results = []
        
        if not term_data:
            logger.warning(f"数据库 {database} 的条目数据为空，跳过。")
            return results
        
        for term_id, term_info in tqdm(
            term_data.items(),
            desc=f"Analyzing {database}",
            leave=False
        ):
            try:
                term_genes = set(term_info.get("genes", []))
                
                result = self.method.calculate_enrichment(
                    gene_set=gene_set,
                    background_set=background_set,
                    term_genes=term_genes,
                    term_name=term_info.get("name", term_id),
                    term_id=term_id,
                    database=database,
                    ranked_genes=None,
                    gene_weights=None
                )
                
                if result is not None:
                    results.append(result)
            except Exception as e:
                logger.warning(f"计算条目 {term_id}（{database}）时出错，已跳过：{e}")
                continue
        
        return results
    
    def run_analysis(
        self,
        gene_set: Set[str],
        background_set: Set[str],
        database_data: Dict[str, Dict[str, Dict[str, Any]]],
        parallel: bool = True,
        ranked_gene_list: Optional[List[Tuple[str, float]]] = None
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
        - 串行模式（parallel=False 或 n_jobs = 1）：
          逐个处理数据库

        对于 GSEA/ssGSEA 方法，ranked_gene_list 会被传递给 analyze_database，
        进而传递给方法实例的 calculate_enrichment。

        参数:
            gene_set: 输入基因集合（查询基因列表）
            background_set: 背景基因集合
            database_data: 数据库数据字典，格式为：
                          {database_name: {term_id: {"name": str, "genes": List[str]}}}
            parallel: 是否使用并行处理，默认为 True
            ranked_gene_list: 可选的 (基因, 权重) 元组列表（GSEA/ssGSEA 使用）

        返回值:
            Dict[str, pd.DataFrame]: 各数据库的富集分析结果，
                                    键为数据库名称，值为结果 DataFrame
        """
        self.results = {}
        
        if parallel and self.config.n_jobs > 1 and len(database_data) > 1:
            # 并行处理模式：使用线程池并行分析多个数据库
            logger.info(f"使用并行模式分析 {len(database_data)} 个数据库（线程数: {self.config.n_jobs}）")
            with ThreadPoolExecutor(max_workers=self.config.n_jobs) as executor:
                # 提交所有数据库的分析任务
                futures = {
                    executor.submit(
                        self.analyze_database,
                        gene_set,
                        background_set,
                        term_data,
                        database,
                        ranked_gene_list
                    ): database
                    for database, term_data in database_data.items()
                }
                
                # 收集并行任务的结果
                for future in as_completed(futures):
                    database = futures[future]
                    try:
                        results = future.result()
                        # gseapy 已内置 FDR 校正，GSEA 方法跳过二次校正
                        if self.config.method != EnrichmentMethod.GSEA.value:
                            results = self.adjust_pvalues(results, self.config.correction)
                        results = self.filter_results(results)
                        results.sort(key=lambda x: x.adjusted_pvalue)
                        self.results[database] = results
                        logger.info(f"Completed {database}: {len(results)} enriched terms")
                    except Exception as e:
                        logger.error(f"Error analyzing {database}: {e}")
        else:
            # 串行处理模式：逐个分析数据库
            for database, term_data in database_data.items():
                try:
                    results = self.analyze_database(
                        gene_set, background_set, term_data, database,
                        ranked_gene_list=ranked_gene_list
                    )
                    # gseapy 已内置 FDR 校正，GSEA 方法跳过二次校正
                    if self.config.method != EnrichmentMethod.GSEA.value:
                        results = self.adjust_pvalues(results, self.config.correction)
                    results = self.filter_results(results)
                    results.sort(key=lambda x: x.adjusted_pvalue)
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
    
    def save_results(self, output_dir: str, metadata: dict = None) -> None:
        """
        将分析结果保存为 TSV 文件

        每个数据库的结果保存为单独的 TSV 文件，文件名格式为：
        {database_name}_enrichment.tsv

        当提供 metadata 时，会在 TSV 文件头部写入以 # 开头的注释行，
        记录 AllEnricher 版本、分析日期、数据库版本等信息。

        参数:
            output_dir: 输出目录路径，如果不存在会自动创建
            metadata: 可选的元数据字典，包含版本和分析信息
        """
        output_path = Path(output_dir)
        output_path.mkdir(parents=True, exist_ok=True)  # 自动创建输出目录
        
        for database, results in self.results.items():
            if results:
                # 保存为 TSV 格式（制表符分隔）
                output_file = output_path / f"{database}_enrichment.tsv"
                data = [r.to_dict() for r in results]
                df = pd.DataFrame(data)

                if metadata:
                    # 写入元数据注释行
                    header_lines = []
                    header_lines.append(f"# AllEnricher version: {metadata.get('allenricher_version', 'unknown')}")
                    header_lines.append(f"# Analysis date: {metadata.get('analysis_date', 'unknown')}")
                    header_lines.append(f"# Database version: {metadata.get('database_version', 'unknown')}")
                    header_lines.append(f"# Species: {metadata.get('species', 'unknown')}")
                    header_lines.append("# Source data versions:")
                    source_versions = metadata.get("source_versions", {})
                    if source_versions:
                        for src_name, src_ver in source_versions.items():
                            header_lines.append(f"#   {src_name}: {src_ver}")
                    header_lines.append("#")
                    with open(output_file, 'w', encoding='utf-8', newline='') as f:
                        f.write("\n".join(header_lines) + "\n")
                        df.to_csv(f, sep='\t', index=False, lineterminator='\n')
                else:
                    df.to_csv(output_file, sep='\t', index=False, lineterminator='\n')

                logger.info(f"Saved {database} results to {output_file}")

    def get_annotated_genes(self, gene_set_data: Dict[str, Dict[str, Any]]) -> Set[str]:
        """
        从注释数据中提取所有有注释的基因（clusterProfiler方案）
        
        遍历所有 term 的 gene_list，收集所有出现过的基因
        
        参数:
            gene_set_data: 注释数据，格式为 {term_id: {"genes": set, "name": str, ...}}
        
        返回:
            所有有注释的基因集合
        """
        annotated = set()
        for term_data in gene_set_data.values():
            genes = term_data.get("genes", set())
            if isinstance(genes, set):
                annotated.update(genes)
            elif isinstance(genes, (list, tuple)):
                annotated.update(genes)
        return annotated

    def resolve_background(self, 
                           gene_set_data: Dict[str, Dict[str, Any]],
                           user_background: Optional[Set[str]] = None,
                           background_mode: str = "annotated") -> Set[str]:
        """
        解析背景基因集
        
        参数:
            gene_set_data: 注释数据
            user_background: 用户自定义背景基因集
            background_mode: 背景模式
                - "annotated": 使用有注释的基因（默认，clusterProfiler方案）
                - "genome": 使用全基因组基因（需要外部提供，通过 user_background）
                - "custom": 使用用户提供的 user_background
        
        返回:
            背景基因集合
        
        异常:
            ValueError: 当 background_mode="custom" 但未提供 user_background 时
        """
        if background_mode == "custom":
            if not user_background:
                raise ValueError("background_mode='custom' requires user_background to be provided")
            return user_background
        elif background_mode == "annotated":
            return self.get_annotated_genes(gene_set_data)
        elif background_mode == "genome":
            if not user_background:
                # genome 模式下如果没有提供 user_background，
                # 回退到 annotated 模式
                logger.warning("background_mode='genome' but no user_background provided, falling back to 'annotated'")
                return self.get_annotated_genes(gene_set_data)
            return user_background
        else:
            raise ValueError(f"Unknown background_mode: {background_mode}. Expected 'annotated', 'genome', or 'custom'")
