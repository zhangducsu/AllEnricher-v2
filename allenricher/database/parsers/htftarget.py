"""
hTFtarget 数据库解析器

解析 hTFtarget 的 TF-target 关系 TSV 文件。

输入文件：
- tf-target-information.txt: TF, target, tissue

输出格式：
- {species}.hTF_2gene.tab.gz: Gene x TF 0/1 矩阵
- {species}.hTF_2disc.gz: TF 描述信息
"""

import gzip
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple
from collections import defaultdict
import logging

import pandas as pd

logger = logging.getLogger(__name__)


class HTFtargetParser:
    """hTFtarget 数据库解析器

    解析 hTFtarget 的 TF-target 关系文件。
    """

    @staticmethod
    def parse_tsv(tsv_path: str) -> Tuple[Dict[str, Set[str]], Dict[str, Set[str]], Dict[str, Set[str]]]:
        """解析 hTFtarget TSV 文件

        Args:
            tsv_path: tf-target-information.txt 文件路径

        Returns:
            (tf_to_targets, gene_to_tfs, tf_to_tissues)
            - tf_to_targets: {TF: {target1, target2, ...}}
            - gene_to_tfs: {gene: {TF1, TF2, ...}}
            - tf_to_tissues: {TF: {tissue1, tissue2, ...}}
        """
        tf_to_targets: Dict[str, Set[str]] = defaultdict(set)
        gene_to_tfs: Dict[str, Set[str]] = defaultdict(set)
        tf_to_tissues: Dict[str, Set[str]] = defaultdict(set)

        count = 0
        with open(tsv_path, 'r') as f:
            header = f.readline()  # 跳过表头: TF\ttarget\ttissue
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) < 3:
                    continue

                tf = parts[0].strip()
                target = parts[1].strip()
                tissues = parts[2].strip()

                if not tf or not target:
                    continue

                tf_to_targets[tf].add(target)
                gene_to_tfs[target].add(tf)

                # tissue 可能是逗号分隔的多个值
                for tissue in tissues.split(','):
                    tissue = tissue.strip()
                    if tissue:
                        tf_to_tissues[tf].add(tissue)

                count += 1

        logger.info(f"hTFtarget: 解析 {count} 条 TF-target 关系, "
                     f"{len(tf_to_targets)} 个 TF, {len(gene_to_tfs)} 个靶基因")
        return dict(tf_to_targets), dict(gene_to_tfs), dict(tf_to_tissues)

    @staticmethod
    def build_database(tsv_path: str, output_dir: str, species: str,
                       valid_genes: Optional[Set[str]] = None) -> None:
        """构建 hTFtarget 数据库

        Args:
            tsv_path: tf-target-information.txt 路径
            output_dir: 输出目录
            species: 物种代码（如 hsa）
            valid_genes: 有效基因集合（可选）
        """
        outdir = Path(output_dir)
        outdir.mkdir(parents=True, exist_ok=True)

        logger.info(f"HTFtargetParser: 开始构建数据库 (species={species})")

        tf_to_targets, gene_to_tfs, tf_to_tissues = HTFtargetParser.parse_tsv(tsv_path)

        # 过滤有效基因
        if valid_genes:
            tf_to_targets = {
                tf: {g for g in targets if g in valid_genes}
                for tf, targets in tf_to_targets.items()
            }
            tf_to_targets = {tf: targets for tf, targets in tf_to_targets.items() if targets}
            logger.info(f"过滤后: {len(tf_to_targets)} 个 TF")

        # 获取所有基因
        all_genes = set()
        for targets in tf_to_targets.values():
            all_genes.update(targets)
        all_tfs = set(tf_to_targets.keys())

        # 构建 Gene x TF 矩阵
        gene_list = sorted(all_genes)
        tf_list = sorted(all_tfs)

        gene2tf_file = outdir / f"{species}.hTF_2gene.tab.gz"
        logger.info(f"写入文件: {gene2tf_file}")

        with gzip.open(gene2tf_file, 'wt') as f:
            f.write('Gene\t' + '\t'.join(tf_list) + '\n')
            for gene in gene_list:
                regulating_tfs = gene_to_tfs.get(gene, set())
                values = ['1' if tf in regulating_tfs else '0' for tf in tf_list]
                f.write(gene + '\t' + '\t'.join(values) + '\n')

        # 构建 TF 描述文件
        disc_file = outdir / f"{species}.hTF_2disc.gz"
        logger.info(f"写入文件: {disc_file}")

        with gzip.open(disc_file, 'wt') as f:
            f.write("TF\ttarget_count\ttissues\tsource\n")
            for tf in tf_list:
                target_count = len(tf_to_targets[tf])
                tissues = ','.join(sorted(tf_to_tissues.get(tf, set())))
                f.write(f"{tf}\t{target_count}\t{tissues}\thTFtarget\n")

        logger.info(f"HTFtargetParser: 数据库构建完成")
