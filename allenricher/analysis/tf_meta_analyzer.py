"""Cross-library consensus utilities for transcription-factor enrichment."""

import pandas as pd
from typing import Dict, List, Set, Optional, Tuple


class TFMetaAnalyzer:
    """Combine TF enrichment results from supported local libraries.
    
    TRRUST, ChEA3, AnimalTFDB, and hTFtarget results can be combined by union,
    intersection, MeanRank, or TopRank. P values from dependent libraries are not
    pooled as if they were independent experiments."""
    
    SUPPORTED_DATABASES = ['trrust', 'chea3', 'animaltfdb', 'htftarget']
    
    def __init__(self, analyzers: Dict[str, 'TFEnrichmentAnalyzer']):
        """Initialize the analyzer with preloaded TF-library analyzers."""
        self.analyzers = analyzers
    
    def get_available_databases(
        self,
        species: str,
        database_dir: str,
        version: Optional[str] = None
    ) -> List[str]:
        """Return TF databases that can be loaded for the selected species."""
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
                # Availability probing should not make the whole list operation fail.
                continue
        
        return available
    
    def analyze(
        self,
        gene_set: List[str],
        databases: Optional[List[str]] = None,
        method: str = 'ora',
        min_overlap: int = 3,
        background_size: int = 20000,
        combine_method: str = 'meanrank'
    ) -> pd.DataFrame:
        """Run local TF ORA and combine the resulting library tables."""
        from ..database.manager import DatabaseManager
        from .tf_enrichment import TFEnrichmentAnalyzer
        
        # Use all initialized analyzers when no subset is requested.
        if databases is None:
            if not self.analyzers:
                raise ValueError("No TF databases were requested and no analyzers are available")
            databases = list(self.analyzers.keys())
        
        results: Dict[str, pd.DataFrame] = {}
        
        # Run each initialized TF library independently.
        for db_name in databases:
            db_key = db_name.upper()
            if db_key in self.analyzers:
                analyzer = self.analyzers[db_key]
                if method == 'ora':
                    result = analyzer.ora(gene_set, min_overlap=min_overlap)
                    # Add normalized aliases required by consensus functions.
                    if 'Pvalue' in result.columns and 'p_value' not in result.columns:
                        result['p_value'] = result['Pvalue']
                    if 'FDR' in result.columns and 'fdr' not in result.columns:
                        result['fdr'] = result['FDR']
                    results[db_name.lower()] = result
                elif method == 'gsea':
                    raise NotImplementedError(
                        "TF GSEA requires a ranked gene table; call analyzer.gsea() directly"
                    )
        
        if not results:
            return pd.DataFrame()
        
        return self.combine_results(results, method=combine_method)
    
    def combine_results(
        self,
        results_dict: Dict[str, pd.DataFrame],
        method: str = "union"
    ) -> pd.DataFrame:
        """Combine TF result tables using the selected consensus strategy."""
        if method == "union":
            return self._combine_union(results_dict)
        elif method == "intersection":
            return self._combine_intersection(results_dict)
        elif method in {"meanrank", "toprank"}:
            return self.rank_consensus(results_dict, method=method)
        elif method == "meta":
            raise ValueError(
                "Stouffer meta-analysis is disabled because TF-library evidence "
                "overlaps. Use meanrank or toprank."
            )
        else:
            raise ValueError(
                f"Unknown method: {method}. Use union, intersection, meanrank, or toprank"
            )
    
    def _combine_union(
        self,
        results_dict: Dict[str, pd.DataFrame]
    ) -> pd.DataFrame:
        """Return the best row for every TF observed in any library."""
        combined_data = []
        
        for db_name, df in results_dict.items():
            if df is None or df.empty:
                continue
            df_copy = df.copy()
            df_copy['source_db'] = db_name
            combined_data.append(df_copy)
        
        if not combined_data:
            return pd.DataFrame()
        
        combined = pd.concat(combined_data, ignore_index=True)
        
        # Aggregate source libraries and retain the best row per TF.
        if 'TF' in combined.columns:
            tf_sources = combined.groupby('TF')['source_db'].apply(
                lambda x: ','.join(sorted(set(x)))
            ).reset_index()
            tf_sources.columns = ['TF', 'sources']
            
            best_results = combined.loc[combined.groupby('TF')['p_value'].idxmin()]
            best_results = best_results.drop(columns=['source_db'], errors='ignore')
            
            result = best_results.merge(tf_sources, on='TF', how='left')
            return result
        
        return combined
    
    def _combine_intersection(
        self,
        results_dict: Dict[str, pd.DataFrame]
    ) -> pd.DataFrame:
        """Return TFs observed in every non-empty input library."""
        # Identify TFs represented in every non-empty library.
        tf_sets = {}
        for db_name, df in results_dict.items():
            if df is not None and not df.empty and 'TF' in df.columns:
                tf_sets[db_name] = set(df['TF'].unique())
        
        if not tf_sets:
            return pd.DataFrame()
        
        common_tfs = set.intersection(*tf_sets.values())
        
        if not common_tfs:
            return pd.DataFrame()
        
        combined_data = []
        for db_name, df in results_dict.items():
            if df is not None and not df.empty:
                df_copy = df[df['TF'].isin(common_tfs)].copy()
                df_copy['source_db'] = db_name
                combined_data.append(df_copy)
        
        if not combined_data:
            return pd.DataFrame()
        
        combined = pd.concat(combined_data, ignore_index=True)
        
        tf_sources = combined.groupby('TF')['source_db'].apply(
            lambda x: ','.join(sorted(set(x)))
        ).reset_index()
        tf_sources.columns = ['TF', 'sources']
        
        best_results = combined.loc[combined.groupby('TF')['p_value'].idxmin()]
        best_results = best_results.drop(columns=['source_db'], errors='ignore')
        
        result = best_results.merge(tf_sources, on='TF', how='left')
        return result

    @staticmethod
    def rank_consensus(
        results_dict: Dict[str, pd.DataFrame],
        method: str = "meanrank",
    ) -> pd.DataFrame:
        """Combine within-library TF ranks using MeanRank or TopRank."""
        if method not in {"meanrank", "toprank"}:
            raise ValueError("method must be 'meanrank' or 'toprank'")

        ranked = []
        for source, frame in results_dict.items():
            if frame is None or frame.empty or 'TF' not in frame.columns:
                continue
            p_column = 'Pvalue' if 'Pvalue' in frame.columns else 'p_value'
            if p_column not in frame.columns:
                continue
            best = (
                frame[['TF', p_column]].dropna().sort_values([p_column, 'TF'])
                .drop_duplicates('TF', keep='first')
            )
            if best.empty:
                continue
            best = best.assign(
                Source=source,
                Scaled_Rank=best[p_column].rank(method='min') / len(best),
            )
            ranked.append(best[['TF', 'Source', 'Scaled_Rank']])

        if not ranked:
            return pd.DataFrame()
        long = pd.concat(ranked, ignore_index=True)
        grouped = long.groupby('TF', sort=True)
        score = grouped['Scaled_Rank'].mean() if method == 'meanrank' else grouped['Scaled_Rank'].min()
        result = score.rename('Consensus_Score').reset_index()
        result['Source_Count'] = result['TF'].map(grouped['Source'].nunique())
        result['Sources'] = result['TF'].map(
            grouped['Source'].apply(lambda values: ','.join(sorted(set(values))))
        )
        result['Consensus_Method'] = 'MeanRank' if method == 'meanrank' else 'TopRank'
        return result.sort_values(['Consensus_Score', 'TF']).reset_index(drop=True)
    
    def meta_analysis(
        self,
        results_dict: Dict[str, pd.DataFrame]
    ) -> pd.DataFrame:
        """Reject invalid P-value pooling across dependent TF libraries."""
        raise ValueError(
            "Stouffer meta-analysis is disabled because the TF libraries are "
            "not independent. Use rank_consensus(..., method='meanrank') or "
            "method='toprank'."
        )
    
    def calculate_consistency_score(
        self,
        results_dict: Dict[str, pd.DataFrame],
        fdr_threshold: float = 0.05
    ) -> pd.DataFrame:
        """Summarize how consistently each TF is significant across libraries."""
        # Record whether each TF passes the threshold in each library.
        tf_significance = {}  # {tf: {db: is_significant}}
        tf_info = {}          # {tf: {db: row_dict}}
        
        for db_name, df in results_dict.items():
            if df is None or df.empty or 'TF' not in df.columns:
                continue
            
            for _, row in df.iterrows():
                tf = row['TF']
                
                # Prefer FDR; fall back to a raw P value for legacy tables.
                fdr = row.get('fdr', row.get('FDR', row.get('p_value', row.get('Pvalue', 1.0))))
                is_significant = fdr < fdr_threshold
                
                if tf not in tf_significance:
                    tf_significance[tf] = {}
                    tf_info[tf] = {}
                
                tf_significance[tf][db_name] = is_significant
                tf_info[tf][db_name] = row.to_dict()
        
        consistency_results = []
        
        for tf, sig_dict in tf_significance.items():
            n_significant = sum(sig_dict.values())
            n_total = len(sig_dict)
            
            consistency_score = n_significant / n_total if n_total > 0 else 0
            
            significant_dbs = [db for db, sig in sig_dict.items() if sig]
            non_significant_dbs = [db for db, sig in sig_dict.items() if not sig]
            
            # Retain the strongest raw significance across available libraries.
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
        
        # Rank TFs by consistency, support count, and best P value.
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
        """Compare TRRUST regulatory modes with another TF result table."""
        if trrust_result is None or trrust_result.empty or 'TF' not in trrust_result.columns:
            return pd.DataFrame()
        
        if other_result is None or other_result.empty or 'TF' not in other_result.columns:
            return pd.DataFrame()
        
        trrust_tfs = set(trrust_result['TF'].unique())
        other_tfs = set(other_result['TF'].unique())
        common_tfs = trrust_tfs & other_tfs
        
        if not common_tfs:
            return pd.DataFrame()
        
        comparison_results = []
        
        for tf in common_tfs:
            trrust_row = trrust_result[trrust_result['TF'] == tf].iloc[0]
            other_row = other_result[other_result['TF'] == tf].iloc[0]
            
            # Regulatory mode is supplied by TRRUST.
            trrust_mode = trrust_row.get('Mode', trrust_row.get('mode', 'unknown'))
            
            # Resolve public and legacy significance column names.
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
        
        # Present the strongest TRRUST evidence first.
        result_df = result_df.sort_values('trrust_p_value', ascending=True)
        
        return result_df
    
    def get_high_confidence_tfs(
        self,
        results_dict: Dict[str, pd.DataFrame],
        min_databases: int = 2,
        fdr_threshold: float = 0.05
    ) -> pd.DataFrame:
        """Return TFs significant in at least the requested number of libraries."""
        consistency_df = self.calculate_consistency_score(results_dict, fdr_threshold)
        
        if consistency_df.empty:
            return pd.DataFrame()
        
        high_conf = consistency_df[
            consistency_df['n_significant_databases'] >= min_databases
        ].copy()
        
        return high_conf
