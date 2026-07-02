"""
TRRUST v2 数据库解析器

解析 TRRUST TSV 文件（TF-target 转录调控关系），
生成 AllEnricher 标准的 TF2target.tab.gz、gene2TF.tab.gz 和 TF2disc.gz 文件。

输入文件格式 (TRRUST TSV):
    TF基因名\\ttarget基因名\\tmode_of_regulation(Activation/Repression/Unknown)

输出文件格式：
    - {species}.TF2target.tab.gz: TF\\ttarget1\\ttarget2\\t... (0/1 矩阵)
    - {species}.gene2TF.tab.gz: Gene\\tTF1\\tTF2\\t... (0/1 矩阵)
    - {species}.TF2disc.gz: TF\\tmode\\ttarget_count
"""

import gzip
import logging
from pathlib import Path
from typing import Dict, Optional, Set, Tuple

logger = logging.getLogger(__name__)


class TRRUSTParser:
    """TRRUST v2 数据库解析器

    解析 TRRUST TSV 文件，提取 TF-target 转录调控关系，
    生成 AllEnricher 标准格式的数据库文件。

    输入文件格式 (TRRUST TSV):
        TF基因名\\ttarget基因名\\tmode_of_regulation

        mode_of_regulation 取值: Activation / Repression / Unknown

    输出文件格式：
    - {species}.TF2target.tab.gz: TF\\ttarget1\\ttarget2\\t... (0/1 矩阵)
    - {species}.gene2TF.tab.gz: Gene\\tTF1\\tTF2\\t... (0/1 矩阵)
    - {species}.TF2disc.gz: TF\\tmode\\ttarget_count
    """

    @staticmethod
    def parse_tsv(tsv_path: str) -> Tuple[
        Dict[str, Set[str]],
        Dict[str, Set[str]],
        Dict[str, str]
    ]:
        """解析 TRRUST TSV 文件，提取 TF-target 关系

        读取 TRRUST TSV 文件，构建 TF->targets 映射、
        gene->TFs 反向映射，以及 TF 调控模式统计。

        Args:
            tsv_path: TRRUST TSV 文件路径（支持 .gz 压缩格式）

        Returns:
            三元组 (tf_to_targets, gene_to_tfs, tf_modes):
            - tf_to_targets: {TF: {target1, target2, ...}}
            - gene_to_tfs: {gene: {TF1, TF2, ...}}
            - tf_modes: {TF: 'activator'/'repressor'/'mixed'/'unknown'}

        Raises:
            FileNotFoundError: 输入文件不存在时抛出
            ValueError: 没有找到有效的 TF-target 关联时抛出
        """
        # 自动选择打开方式
        if tsv_path.endswith('.gz'):
            f_open = gzip.open(tsv_path, 'rt', encoding='utf-8')
        else:
            f_open = open(tsv_path, 'r', encoding='utf-8')

        tf_to_targets: Dict[str, Set[str]] = {}
        gene_to_tfs: Dict[str, Set[str]] = {}
        tf_mode_counts: Dict[str, Dict[str, int]] = {}  # {TF: {'activation': n, 'repression': m, 'unknown': k}}
        n = 0

        with f_open:
            for line in f_open:
                line = line.strip()
                if not line or line.startswith('#'):
                    continue
                parts = line.split('\t')
                if len(parts) < 3:
                    continue

                tf = parts[0].strip()
                target = parts[1].strip()
                mode = parts[2].strip()

                if not tf or not target:
                    continue

                # 构建 TF -> targets 映射
                if tf not in tf_to_targets:
                    tf_to_targets[tf] = set()
                tf_to_targets[tf].add(target)

                # 构建 gene -> TFs 反向映射
                if target not in gene_to_tfs:
                    gene_to_tfs[target] = set()
                gene_to_tfs[target].add(tf)

                # 统计 TF 的调控模式计数
                if tf not in tf_mode_counts:
                    tf_mode_counts[tf] = {'activation': 0, 'repression': 0, 'unknown': 0}
                mode_lower = mode.lower()
                if mode_lower == 'activation':
                    tf_mode_counts[tf]['activation'] += 1
                elif mode_lower == 'repression':
                    tf_mode_counts[tf]['repression'] += 1
                else:
                    tf_mode_counts[tf]['unknown'] += 1

                n += 1

        if n == 0:
            raise ValueError(
                "[错误] 在 TRRUST 文件中没有找到有效的 TF-target 关联！"
            )

        # 确定 TF 的主要调控模式
        tf_modes: Dict[str, str] = {}
        for tf, counts in tf_mode_counts.items():
            act = counts['activation']
            rep = counts['repression']
            unk = counts['unknown']

            if act > 0 and rep > 0:
                tf_modes[tf] = 'mixed'
            elif act > 0:
                tf_modes[tf] = 'activator'
            elif rep > 0:
                tf_modes[tf] = 'repressor'
            else:
                tf_modes[tf] = 'unknown'

        logger.info("TRRUSTParser: 共解析 %d 条 TF-target 关联", n)
        logger.info("TRRUSTParser: %d 个 TF, %d 个靶基因",
                    len(tf_to_targets), len(gene_to_tfs))

        return tf_to_targets, gene_to_tfs, tf_modes

    @staticmethod
    def build_database(tsv_path: str, output_dir: str, species: str,
                       valid_genes: Optional[Set[str]] = None) -> None:
        """构建 TRRUST 数据库文件

        解析 TRRUST TSV 文件，生成 TF2target.tab.gz、gene2TF.tab.gz
        和 TF2disc.gz 三个数据库文件。

        Args:
            tsv_path: TRRUST TSV 文件路径（支持 .gz 压缩格式）
            output_dir: 输出目录
            species: 物种缩写（如 hsa, mmu）
            valid_genes: 可选的有效基因集合，如果提供则过滤不在集合中的基因

        Raises:
            FileNotFoundError: 输入文件不存在时抛出
            ValueError: 没有找到有效的 TF-target 关联时抛出
        """
        outdir_path = Path(output_dir)
        outdir_path.mkdir(parents=True, exist_ok=True)

        logger.info("TRRUSTParser: 开始构建 TRRUST 数据库 (species=%s)", species)

        # 第一步：解析 TSV 文件
        tf_to_targets, gene_to_tfs, tf_modes = TRRUSTParser.parse_tsv(tsv_path)

        # 第二步：如果提供了 valid_genes，过滤不在集合中的基因
        if valid_genes is not None:
            logger.info("使用 valid_genes 过滤 (共 %d 个有效基因)", len(valid_genes))

            # 过滤 TF: TF 本身必须是有效基因
            filtered_tf_to_targets: Dict[str, Set[str]] = {}
            for tf, targets in tf_to_targets.items():
                if tf not in valid_genes:
                    continue
                filtered_targets = {t for t in targets if t in valid_genes}
                if filtered_targets:
                    filtered_tf_to_targets[tf] = filtered_targets

            # 重建 gene_to_tfs
            filtered_gene_to_tfs: Dict[str, Set[str]] = {}
            for tf, targets in filtered_tf_to_targets.items():
                for target in targets:
                    if target not in filtered_gene_to_tfs:
                        filtered_gene_to_tfs[target] = set()
                    filtered_gene_to_tfs[target].add(tf)

            # 过滤 tf_modes
            filtered_tf_modes: Dict[str, str] = {
                tf: mode for tf, mode in tf_modes.items()
                if tf in filtered_tf_to_targets
            }

            tf_to_targets = filtered_tf_to_targets
            gene_to_tfs = filtered_gene_to_tfs
            tf_modes = filtered_tf_modes

            logger.info("过滤后: %d 个 TF, %d 个靶基因",
                        len(tf_to_targets), len(gene_to_tfs))

        if not tf_to_targets:
            raise ValueError(
                "[错误] 过滤后没有有效的 TF-target 关联！"
            )

        # 第三步：收集所有基因（TF + target），作为 gene2TF 矩阵的行
        all_genes = set(gene_to_tfs.keys()) | set(tf_to_targets.keys())
        sorted_tfs = sorted(tf_to_targets.keys())
        sorted_genes = sorted(all_genes)

        # 第四步：写入 TF2target.tab.gz
        # 格式: TF\\ttarget1\\ttarget2\\t... (0/1 矩阵)
        sorted_targets = sorted(gene_to_tfs.keys())
        tab_file = outdir_path / f"{species}.TF2target.tab.gz"
        logger.info("写入文件: %s", tab_file)

        with gzip.open(tab_file, 'wt', encoding='utf-8') as f:
            header = ["TF"] + sorted_targets
            f.write('\t'.join(header) + '\n')

            for tf in sorted_tfs:
                row = [tf]
                for target in sorted_targets:
                    if target in tf_to_targets[tf]:
                        row.append('1')
                    else:
                        row.append('0')
                f.write('\t'.join(row) + '\n')

        # 第五步：写入 gene2TF.tab.gz
        # 格式: Gene\\tTF1\\tTF2\\t... (0/1 矩阵)
        gene2tf_file = outdir_path / f"{species}.gene2TF.tab.gz"
        logger.info("写入文件: %s", gene2tf_file)

        with gzip.open(gene2tf_file, 'wt', encoding='utf-8') as f:
            header = ["Gene"] + sorted_tfs
            f.write('\t'.join(header) + '\n')

            for gene in sorted_genes:
                row = [gene]
                for tf in sorted_tfs:
                    if tf in gene_to_tfs.get(gene, set()):
                        row.append('1')
                    else:
                        row.append('0')
                f.write('\t'.join(row) + '\n')

        # 第六步：写入 TF2disc.gz
        # 格式: TF\\tmode\\ttarget_count
        disc_file = outdir_path / f"{species}.TF2disc.gz"
        logger.info("写入文件: %s", disc_file)

        with gzip.open(disc_file, 'wt', encoding='utf-8') as f:
            for tf in sorted_tfs:
                mode = tf_modes.get(tf, 'unknown')
                target_count = len(tf_to_targets[tf])
                f.write(f"{tf}\t{mode}\t{target_count}\n")

        logger.info("TRRUSTParser: TRRUST 数据库构建完成")
