"""
同源映射引擎

通过 AnimalTFDB 的直系同源映射，将人类 TF-target 关系推断到目标物种。

映射流程：
1. 人类 TF → 靶基因 (hTFtarget)
2. 目标物种基因 → 人类基因 (AnimalTFDB ortholog_to_human, 反向)
3. 目标物种 TF → 人类 TF (AnimalTFDB ortholog_to_human)
4. 合并得到：目标物种 TF → 目标物种靶基因
"""

import gzip
from pathlib import Path
from typing import Any, Dict, List, Optional, Set, Tuple
from collections import defaultdict
import logging

import pandas as pd

logger = logging.getLogger(__name__)


class OrthologMapper:
    """同源映射引擎

    将人类 TF-target 关系通过直系同源映射推断到目标物种。
    """

    # 默认去重策略
    DEFAULT_DEDUP_STRATEGY = "none"

    def __init__(
        self,
        human_tf_to_targets: Dict[str, Set[str]],
        species_to_human: Dict[str, str],
        species_tf_set: Optional[Set[str]] = None,
        dedup_strategy: str = "none",
    ):
        """初始化同源映射引擎

        Args:
            human_tf_to_targets: 人类 TF→靶基因关系 {TF: {target1, target2, ...}}
            species_to_human: 目标物种基因→人类基因映射 {sp_gene: human_gene}
            species_tf_set: 目标物种的 TF 集合（可选，用于过滤）
            dedup_strategy: 多对一映射去重策略
                - "none": 保留所有映射（默认行为）
                - "first": 对每个人类基因只保留第一个（按字母排序）物种基因
                - "all": 保留所有映射但标记为重复
        """
        self.human_tf_to_targets = human_tf_to_targets
        self.species_to_human = species_to_human
        self.species_tf_set = species_tf_set

        # 验证并设置去重策略
        valid_strategies = {"none", "first", "all"}
        if dedup_strategy not in valid_strategies:
            logger.warning(
                f"无效的 dedup_strategy '{dedup_strategy}'，使用默认值 '{self.DEFAULT_DEDUP_STRATEGY}'"
            )
            dedup_strategy = self.DEFAULT_DEDUP_STRATEGY
        self.dedup_strategy = dedup_strategy

        # 构建反向映射：人类基因 → 目标物种基因
        self.human_to_species: Dict[str, Set[str]] = defaultdict(set)
        for sp_gene, hu_gene in species_to_human.items():
            self.human_to_species[hu_gene].add(sp_gene)

        # 构建目标物种 TF → 人类 TF 映射
        self.species_tf_to_human_tf: Dict[str, str] = {}
        if species_tf_set:
            for sp_tf in species_tf_set:
                if sp_tf in species_to_human:
                    self.species_tf_to_human_tf[sp_tf] = species_to_human[sp_tf]

    def get_duplicate_stats(self) -> Dict[str, Any]:
        """获取多对一映射统计

        Returns:
            包含多对一映射统计信息的字典:
            - total_human_genes: 人类基因总数
            - multi_mapping_count: 存在多对一映射的人类基因数量
            - multi_mapping_genes: {人类基因: 映射到的物种基因数量}
            - total_species_genes: 物种基因总数
        """
        # 统计有多少人类基因映射到多个物种基因
        multi_mappings = {
            hu: len(sp_genes)
            for hu, sp_genes in self.human_to_species.items()
            if len(sp_genes) > 1
        }
        return {
            "total_human_genes": len(self.human_to_species),
            "multi_mapping_count": len(multi_mappings),
            "multi_mapping_genes": multi_mappings,
            "total_species_genes": len(self.species_to_human),
        }

    def _get_deduplicated_targets(
        self, human_target: str, sp_targets: Set[str]
    ) -> Set[str]:
        """根据去重策略获取去重后的目标基因集合

        Args:
            human_target: 人类靶基因
            sp_targets: 映射到的物种基因集合

        Returns:
            去重后的物种基因集合
        """
        if self.dedup_strategy == "none" or len(sp_targets) <= 1:
            return sp_targets
        elif self.dedup_strategy == "first":
            # 按字母排序，只保留第一个
            return {sorted(sp_targets)[0]}
        else:  # "all" - 保留所有但标记（这里仍返回全部，调用方会记录）
            return sp_targets

    def map_tf_targets(self) -> Tuple[Dict[str, Set[str]], Dict[str, Set[str]]]:
        """执行同源映射，推断目标物种的 TF-target 关系

        Returns:
            (tf_to_targets, gene_to_tfs)
            - tf_to_targets: {目标物种TF: {目标物种靶基因1, ...}}
            - gene_to_tfs: {目标物种基因: {目标物种TF1, ...}}
        """
        tf_to_targets: Dict[str, Set[str]] = defaultdict(set)
        gene_to_tfs: Dict[str, Set[str]] = defaultdict(set)

        mapped_tf_count = 0
        unmapped_tf_count = 0
        dedup_removed_count = 0  # 因去重移除的映射数

        # 记录重复映射信息（仅当 dedup_strategy == "all" 时使用）
        duplicate_mappings: Dict[str, List[str]] = {}

        for sp_tf, human_tf in self.species_tf_to_human_tf.items():
            # 获取人类 TF 的靶基因
            human_targets = self.human_tf_to_targets.get(human_tf, set())
            if not human_targets:
                unmapped_tf_count += 1
                continue

            # 将人类靶基因映射回目标物种
            for human_target in human_targets:
                sp_targets = self.human_to_species.get(human_target, set())
                if not sp_targets:
                    continue

                # 应用去重策略
                original_count = len(sp_targets)
                deduped_targets = self._get_deduplicated_targets(human_target, sp_targets)
                dedup_removed_count += original_count - len(deduped_targets)

                # 记录重复映射（仅当策略为 "all" 且存在多对一时）
                if self.dedup_strategy == "all" and len(sp_targets) > 1:
                    duplicate_mappings[human_target] = sorted(sp_targets)

                for sp_target in deduped_targets:
                    # 排除自映射（TF不能调控自己）
                    if sp_target != sp_tf:
                        tf_to_targets[sp_tf].add(sp_target)
                        gene_to_tfs[sp_target].add(sp_tf)

            mapped_tf_count += 1

        # 记录去重统计
        if self.dedup_strategy != "none":
            dup_stats = self.get_duplicate_stats()
            logger.info(
                f"去重策略 '{self.dedup_strategy}': "
                f"共 {dup_stats['multi_mapping_count']} 个人类基因存在多对一映射, "
                f"移除了 {dedup_removed_count} 个重复映射"
            )
            if self.dedup_strategy == "all" and duplicate_mappings:
                logger.info(
                    f"标记为重复的映射示例（前5个）: "
                    f"{dict(list(duplicate_mappings.items())[:5])}"
                )

        logger.info(
            f"同源映射完成: {mapped_tf_count} 个TF成功映射, "
            f"{unmapped_tf_count} 个TF无人类对应关系, "
            f"共 {len(tf_to_targets)} 个TF有靶基因, "
            f"共 {len(gene_to_tfs)} 个靶基因"
        )

        return dict(tf_to_targets), dict(gene_to_tfs)

    @staticmethod
    def build_mapped_database(
        tf_to_targets: Dict[str, Set[str]],
        gene_to_tfs: Dict[str, Set[str]],
        species_tf_df: pd.DataFrame,
        output_dir: str,
        species: str,
    ) -> None:
        """构建同源映射后的数据库文件

        Args:
            tf_to_targets: 映射后的 TF→靶基因关系
            gene_to_tfs: 映射后的 基因→TF 关系
            species_tf_df: 目标物种 TF 信息 DataFrame
            output_dir: 输出目录
            species: 物种代码
        """
        outdir = Path(output_dir)
        outdir.mkdir(parents=True, exist_ok=True)

        all_genes = set(gene_to_tfs.keys())
        all_tfs = set(tf_to_targets.keys())
        gene_list = sorted(all_genes)
        tf_list = sorted(all_tfs)

        # 构建 Gene x TF 矩阵
        gene2tf_file = outdir / f"{species}.AnimalTFDB_2gene.tab.gz"
        logger.info(f"写入文件: {gene2tf_file}")

        with gzip.open(gene2tf_file, 'wt') as f:
            f.write('Gene\t' + '\t'.join(tf_list) + '\n')
            for gene in gene_list:
                regulating_tfs = gene_to_tfs.get(gene, set())
                values = ['1' if tf in regulating_tfs else '0' for tf in tf_list]
                f.write(gene + '\t' + '\t'.join(values) + '\n')

        # 构建 TF 描述文件
        disc_file = outdir / f"{species}.AnimalTFDB_mapped_2disc.gz"
        logger.info(f"写入文件: {disc_file}")

        # 构建 TF family 查找表
        tf_family_map = {}
        if species_tf_df is not None and 'Symbol' in species_tf_df.columns:
            for _, row in species_tf_df.iterrows():
                tf_family_map[row['Symbol']] = row.get('Family', 'Unknown')

        with gzip.open(disc_file, 'wt') as f:
            f.write("TF\ttarget_count\tFamily\tsource\n")
            for tf in tf_list:
                target_count = len(tf_to_targets[tf])
                family = tf_family_map.get(tf, 'Unknown')
                f.write(f"{tf}\t{target_count}\t{family}\tAnimalTFDB_ortholog\n")

        logger.info(f"OrthologMapper: 数据库构建完成")
