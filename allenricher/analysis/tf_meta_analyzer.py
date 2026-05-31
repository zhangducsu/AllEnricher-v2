"""
多库联合转录因子分析器

支持 TRRUST + ChEA3 + AnimalTFDB + hTFtarget 联合分析，提供：
1. 结果整合（取并集/交集）
2. Meta 分析（Stouffer's Z-score 方法）
3. 一致性评分（多个库都显著的 TF）
"""

import pandas as pd
import numpy as np
from typing import Dict, List, Set, Optional, Tuple
from scipy import stats


class TFMetaAnalyzer:
    """多库联合 TF 分析器
    
    支持的数据库：
        - TRRUST: 包含 TF 调控模式（activator/repressor/mixed/unknown）
        - ChEA3: 仅包含 TF-靶基因关系，无调控模式
        - AnimalTFDB: 动物转录因子数据库（含同源映射）
        - hTFtarget: 人类转录因子靶基因数据库
    
    所有数据库结果需包含 'TF' 和 'Pvalue'（或 'p_value'）列，
    以便进行 Stouffer's Z-score meta 分析。
    """
    
    # 支持的数据库列表
    SUPPORTED_DATABASES = ['trrust', 'chea3', 'animaltfdb', 'htftarget']
    
    def __init__(self, analyzers: Dict[str, 'TFEnrichmentAnalyzer']):
        """
        Args:
            analyzers: {'TRRUST': analyzer1, 'ChEA3': analyzer2, 'AnimalTFDB': analyzer3, 'hTFtarget': analyzer4}
        """
        self.analyzers = analyzers
    
    def get_available_databases(
        self,
        species: str,
        database_dir: str,
        version: Optional[str] = None
    ) -> List[str]:
        """获取指定物种可用的 TF 数据库列表
        
        检查数据库文件是否存在，返回实际可用的数据库列表。
        
        Args:
            species: 物种代码（如 hsa、mmu）
            database_dir: 数据库文件目录路径
            version: 可选，指定数据库版本
            
        Returns:
            List[str]: 可用的数据库列表（如 ['trrust', 'chea3', 'animaltfdb']）
        """
        from ..database.manager import DatabaseManager
        manager = DatabaseManager(database_dir, species)
        available = []
        
        for db in self.SUPPORTED_DATABASES:
            try:
                if db == 'trrust':
                    data = manager.load_trrust()
                    if data and not data.get('tf2target', pd.DataFrame()).empty:
                        available.append(db)
                elif db == 'chea3':
                    data = manager.load_chea3()
                    if data and not data.get('gene2tf', pd.DataFrame()).empty:
                        available.append(db)
                elif db == 'animaltfdb':
                    data = manager.load_animaltfdb()
                    if data and not data.get('gene2tf', pd.DataFrame()).empty:
                        available.append(db)
                elif db == 'htftarget':
                    data = manager.load_htftarget()
                    if data and not data.get('gene2tf', pd.DataFrame()).empty:
                        available.append(db)
            except Exception:
                # 加载失败时跳过该数据库
                continue
        
        return available
    
    def analyze(
        self,
        gene_set: List[str],
        databases: Optional[List[str]] = None,
        method: str = 'ora',
        min_overlap: int = 3,
        background_size: int = 20000,
        combine_method: str = 'meta'
    ) -> pd.DataFrame:
        """对基因集执行多库联合 TF 分析
        
        加载指定数据库，分别执行富集分析，然后整合结果。
        
        Args:
            gene_set: 输入基因列表（如差异表达基因）
            databases: 要使用的数据库列表，为 None 时使用所有可用的数据库
                       支持: ['trrust', 'chea3', 'animaltfdb', 'htftarget']
            method: 富集分析方法，支持 'ora'（过表示分析）或 'gsea'
            min_overlap: 最小重叠基因数（ORA 方法）
            background_size: 背景基因总数
            combine_method: 结果整合方法，支持 'union'、'intersection'、'meta'
            
        Returns:
            pd.DataFrame: 整合后的分析结果
        """
        from ..database.manager import DatabaseManager
        from .tf_enrichment import TFEnrichmentAnalyzer
        
        # 如果未指定数据库，使用已初始化的 analyzers
        if databases is None:
            if not self.analyzers:
                raise ValueError("未指定数据库且 analyzers 为空，请提供 databases 参数")
            databases = list(self.analyzers.keys())
        
        results: Dict[str, pd.DataFrame] = {}
        
        # 使用已初始化的 analyzers 进行分析
        for db_name in databases:
            db_key = db_name.upper()
            if db_key in self.analyzers:
                analyzer = self.analyzers[db_key]
                if method == 'ora':
                    result = analyzer.ora(gene_set, min_overlap=min_overlap)
                    # 标准化列名：Pvalue -> p_value
                    if 'Pvalue' in result.columns and 'p_value' not in result.columns:
                        result['p_value'] = result['Pvalue']
                    if 'FDR' in result.columns and 'fdr' not in result.columns:
                        result['fdr'] = result['FDR']
                    results[db_name.lower()] = result
                elif method == 'gsea':
                    # GSEA 需要排序列表，这里不支持直接使用 gene_set
                    raise NotImplementedError("GSEA 方法需要排序列表，请直接使用 analyzer.gsea()")
        
        # 整合结果
        if not results:
            return pd.DataFrame()
        
        return self.combine_results(results, method=combine_method)
    
    def combine_results(
        self,
        results_dict: Dict[str, pd.DataFrame],
        method: str = "union"
    ) -> pd.DataFrame:
        """整合多个库的分析结果
        
        Args:
            results_dict: {'trrust': df1, 'chea3': df2, 'animaltfdb': df3, 'htftarget': df4}
            method: 'union'(并集), 'intersection'(交集), 'meta'(meta分析)
            
        Returns:
            整合后的 DataFrame
        """
        if method == "union":
            return self._combine_union(results_dict)
        elif method == "intersection":
            return self._combine_intersection(results_dict)
        elif method == "meta":
            return self.meta_analysis(results_dict)
        else:
            raise ValueError(f"Unknown method: {method}. Use 'union', 'intersection', or 'meta'")
    
    def _combine_union(
        self,
        results_dict: Dict[str, pd.DataFrame]
    ) -> pd.DataFrame:
        """取所有 TF 的并集，标记来源"""
        combined_data = []
        
        for db_name, df in results_dict.items():
            if df is None or df.empty:
                continue
            df_copy = df.copy()
            df_copy['source_db'] = db_name
            combined_data.append(df_copy)
        
        if not combined_data:
            return pd.DataFrame()
        
        # 合并所有结果
        combined = pd.concat(combined_data, ignore_index=True)
        
        # 对每个 TF，汇总来源信息
        if 'TF' in combined.columns:
            tf_sources = combined.groupby('TF')['source_db'].apply(
                lambda x: ','.join(sorted(set(x)))
            ).reset_index()
            tf_sources.columns = ['TF', 'sources']
            
            # 选择每个 TF 的最佳结果（最小 p-value）
            best_results = combined.loc[combined.groupby('TF')['p_value'].idxmin()]
            best_results = best_results.drop(columns=['source_db'], errors='ignore')
            
            # 合并来源信息
            result = best_results.merge(tf_sources, on='TF', how='left')
            return result
        
        return combined
    
    def _combine_intersection(
        self,
        results_dict: Dict[str, pd.DataFrame]
    ) -> pd.DataFrame:
        """只保留在多个库中都出现的 TF"""
        # 获取每个库中的 TF 集合
        tf_sets = {}
        for db_name, df in results_dict.items():
            if df is not None and not df.empty and 'TF' in df.columns:
                tf_sets[db_name] = set(df['TF'].unique())
        
        if not tf_sets:
            return pd.DataFrame()
        
        # 计算交集
        common_tfs = set.intersection(*tf_sets.values())
        
        if not common_tfs:
            return pd.DataFrame()
        
        # 筛选共同 TF 的结果
        combined_data = []
        for db_name, df in results_dict.items():
            if df is not None and not df.empty:
                df_copy = df[df['TF'].isin(common_tfs)].copy()
                df_copy['source_db'] = db_name
                combined_data.append(df_copy)
        
        if not combined_data:
            return pd.DataFrame()
        
        combined = pd.concat(combined_data, ignore_index=True)
        
        # 添加来源标记
        tf_sources = combined.groupby('TF')['source_db'].apply(
            lambda x: ','.join(sorted(set(x)))
        ).reset_index()
        tf_sources.columns = ['TF', 'sources']
        
        # 选择每个 TF 的最佳结果
        best_results = combined.loc[combined.groupby('TF')['p_value'].idxmin()]
        best_results = best_results.drop(columns=['source_db'], errors='ignore')
        
        result = best_results.merge(tf_sources, on='TF', how='left')
        return result
    
    def meta_analysis(
        self,
        results_dict: Dict[str, pd.DataFrame]
    ) -> pd.DataFrame:
        """使用 Stouffer's Z-score 方法进行 meta 分析
        
        对每个在多个库中出现的 TF，使用 Stouffer's Z-score 方法合并 p 值。
        公式: Z_meta = sum(Z_i * w_i) / sqrt(sum(w_i^2))
        其中 w_i = sqrt(n_i)，n_i 为第 i 个库的样本量或权重
        
        支持的数据库：TRRUST、ChEA3、AnimalTFDB、hTFtarget
        
        Args:
            results_dict: {'trrust': df1, 'chea3': df2, 'animaltfdb': df3, 'htftarget': df4}
            
        Returns:
            Meta 分析结果 DataFrame
        """
        # 收集每个 TF 在各个库中的 p 值
        tf_pvalues = {}  # {tf: {db: p_value}}
        tf_info = {}     # {tf: {db: row_dict}}
        
        for db_name, df in results_dict.items():
            if df is None or df.empty or 'TF' not in df.columns:
                continue
            
            for _, row in df.iterrows():
                tf = row['TF']
                p_value = row.get('p_value', row.get('Pvalue', 1.0))
                
                if tf not in tf_pvalues:
                    tf_pvalues[tf] = {}
                    tf_info[tf] = {}
                
                tf_pvalues[tf][db_name] = p_value
                tf_info[tf][db_name] = row.to_dict()
        
        # 对每个 TF 进行 meta 分析
        meta_results = []
        
        for tf, pvalues in tf_pvalues.items():
            if len(pvalues) < 1:
                continue
            
            # 计算 Stouffer's Z-score
            z_scores = []
            weights = []
            
            for db_name, pval in pvalues.items():
                # 将 p 值转换为 Z 分数（单侧检验）
                # 使用最小值避免 log(0)
                pval = max(pval, 1e-300)
                z = stats.norm.ppf(1 - pval)
                z_scores.append(z)
                
                # 权重：假设每个库权重相等，或使用样本量
                # 这里使用相等权重
                weights.append(1.0)
            
            # Stouffer's Z-score 公式
            z_array = np.array(z_scores)
            w_array = np.array(weights)
            
            z_meta = np.sum(z_array * w_array) / np.sqrt(np.sum(w_array ** 2))
            
            # 将 Z 分数转换回 p 值
            p_meta = 1 - stats.norm.cdf(z_meta)
            
            # 计算 FDR（这里使用简单的 Bonferroni 校正作为近似）
            n_tests = len(tf_pvalues)
            fdr_meta = min(p_meta * n_tests, 1.0)
            
            # 收集来源信息
            sources = ','.join(sorted(pvalues.keys()))
            n_databases = len(pvalues)
            
            # 获取该 TF 的最佳原始信息
            best_db = min(pvalues.keys(), key=lambda x: pvalues[x])
            best_info = tf_info[tf][best_db]
            
            result_row = {
                'TF': tf,
                'p_value_meta': p_meta,
                'fdr_meta': fdr_meta,
                'z_score_meta': z_meta,
                'n_databases': n_databases,
                'sources': sources,
                'p_values_original': str(pvalues),
            }
            
            # 添加原始结果中的其他有用字段
            for key in ['overlap', 'background_size', 'target_genes', 'mode', 'Overlap', 'TF_Targets', 'Mode']:
                if key in best_info:
                    # 标准化字段名
                    normalized_key = key.lower() if key not in ['Overlap', 'TF_Targets', 'Mode'] else key
                    result_row[normalized_key] = best_info[key]
            
            meta_results.append(result_row)
        
        if not meta_results:
            return pd.DataFrame()
        
        result_df = pd.DataFrame(meta_results)
        
        # 按 meta p 值排序
        result_df = result_df.sort_values('p_value_meta', ascending=True)
        
        return result_df
    
    def calculate_consistency_score(
        self,
        results_dict: Dict[str, pd.DataFrame],
        fdr_threshold: float = 0.05
    ) -> pd.DataFrame:
        """计算 TF 的一致性评分
        
        统计每个 TF 在多少个库中显著（FDR < threshold）。
        一致性高的 TF 更可靠。
        
        Args:
            results_dict: {'trrust': df1, 'chea3': df2, 'animaltfdb': df3, 'htftarget': df4}
            fdr_threshold: FDR 阈值，默认为 0.05
            
        Returns:
            一致性评分 DataFrame
        """
        # 收集每个 TF 在各个库中的显著性
        tf_significance = {}  # {tf: {db: is_significant}}
        tf_info = {}          # {tf: {db: row_dict}}
        
        for db_name, df in results_dict.items():
            if df is None or df.empty or 'TF' not in df.columns:
                continue
            
            for _, row in df.iterrows():
                tf = row['TF']
                
                # 获取 FDR 或 p_value
                fdr = row.get('fdr', row.get('FDR', row.get('p_value', row.get('Pvalue', 1.0))))
                is_significant = fdr < fdr_threshold
                
                if tf not in tf_significance:
                    tf_significance[tf] = {}
                    tf_info[tf] = {}
                
                tf_significance[tf][db_name] = is_significant
                tf_info[tf][db_name] = row.to_dict()
        
        # 计算一致性评分
        consistency_results = []
        
        for tf, sig_dict in tf_significance.items():
            # 统计在多少个库中显著
            n_significant = sum(sig_dict.values())
            n_total = len(sig_dict)
            
            # 一致性评分：显著库数 / 总库数
            consistency_score = n_significant / n_total if n_total > 0 else 0
            
            # 收集显著的数据库
            significant_dbs = [db for db, sig in sig_dict.items() if sig]
            non_significant_dbs = [db for db, sig in sig_dict.items() if not sig]
            
            # 获取最佳 p 值
            best_p = min(tf_info[tf][db].get('p_value', tf_info[tf][db].get('Pvalue', 1.0)) for db in sig_dict.keys())
            
            result_row = {
                'TF': tf,
                'consistency_score': consistency_score,
                'n_significant_databases': n_significant,
                'n_total_databases': n_total,
                'significant_in': ','.join(sorted(significant_dbs)) if significant_dbs else '',
                'not_significant_in': ','.join(sorted(non_significant_dbs)) if non_significant_dbs else '',
                'best_p_value': best_p,
            }
            
            consistency_results.append(result_row)
        
        if not consistency_results:
            return pd.DataFrame()
        
        result_df = pd.DataFrame(consistency_results)
        
        # 按一致性评分和显著库数排序
        result_df = result_df.sort_values(
            ['consistency_score', 'n_significant_databases', 'best_p_value'],
            ascending=[False, False, True]
        )
        
        return result_df
    
    def compare_modes(
        self,
        trrust_result: pd.DataFrame,
        other_result: pd.DataFrame,
        other_name: str = 'ChEA3'
    ) -> pd.DataFrame:
        """比较 TRRUST 和其他数据库的调控模式
        
        对于在两个库中都出现的 TF，比较其调控模式是否一致。
        TRRUST 提供 mode (activator/repressor)，其他库可能不提供。
        
        Args:
            trrust_result: TRRUST 分析结果
            other_result: 其他数据库（如 ChEA3、AnimalTFDB、hTFtarget）分析结果
            other_name: 其他数据库的名称
            
        Returns:
            比较结果 DataFrame
        """
        if trrust_result is None or trrust_result.empty or 'TF' not in trrust_result.columns:
            return pd.DataFrame()
        
        if other_result is None or other_result.empty or 'TF' not in other_result.columns:
            return pd.DataFrame()
        
        # 获取共同的 TF
        trrust_tfs = set(trrust_result['TF'].unique())
        other_tfs = set(other_result['TF'].unique())
        common_tfs = trrust_tfs & other_tfs
        
        if not common_tfs:
            return pd.DataFrame()
        
        # 为每个 TF 创建比较结果
        comparison_results = []
        
        for tf in common_tfs:
            trrust_row = trrust_result[trrust_result['TF'] == tf].iloc[0]
            other_row = other_result[other_result['TF'] == tf].iloc[0]
            
            # 获取 TRRUST 的调控模式
            trrust_mode = trrust_row.get('Mode', trrust_row.get('mode', 'unknown'))
            
            # 获取 p 值（兼容不同列名）
            trrust_p = trrust_row.get('p_value', trrust_row.get('Pvalue', 1.0))
            trrust_fdr = trrust_row.get('fdr', trrust_row.get('FDR', trrust_p))
            other_p = other_row.get('p_value', other_row.get('Pvalue', 1.0))
            other_fdr = other_row.get('fdr', other_row.get('FDR', other_p))
            
            result_row = {
                'TF': tf,
                'trrust_mode': trrust_mode,
                'trrust_p_value': trrust_p,
                'trrust_fdr': trrust_fdr,
                f'{other_name.lower()}_p_value': other_p,
                f'{other_name.lower()}_fdr': other_fdr,
                'mode_available': trrust_mode != 'unknown',
            }
            
            comparison_results.append(result_row)
        
        if not comparison_results:
            return pd.DataFrame()
        
        result_df = pd.DataFrame(comparison_results)
        
        # 按 TRRUST p 值排序
        result_df = result_df.sort_values('trrust_p_value', ascending=True)
        
        return result_df
    
    def get_high_confidence_tfs(
        self,
        results_dict: Dict[str, pd.DataFrame],
        min_databases: int = 2,
        fdr_threshold: float = 0.05
    ) -> pd.DataFrame:
        """获取高置信度 TF
        
        筛选在至少 min_databases 个库中显著的 TF。
        
        Args:
            results_dict: {'trrust': df1, 'chea3': df2, 'animaltfdb': df3, 'htftarget': df4}
            min_databases: 最小显著库数，默认为 2
            fdr_threshold: FDR 阈值，默认为 0.05
            
        Returns:
            高置信度 TF DataFrame
        """
        consistency_df = self.calculate_consistency_score(results_dict, fdr_threshold)
        
        if consistency_df.empty:
            return pd.DataFrame()
        
        # 筛选在至少 min_databases 个库中显著的 TF
        high_conf = consistency_df[
            consistency_df['n_significant_databases'] >= min_databases
        ].copy()
        
        return high_conf