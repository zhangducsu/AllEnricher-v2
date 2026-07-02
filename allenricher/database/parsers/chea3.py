"""
ChEA3 数据库解析器

解析 ChEA3 GMT 格式文件和 API 返回结果，
生成 AllEnricher 标准的 ChEA3_2gene.tab.gz 和 ChEA3_2disc.gz 文件。

ChEA3 GMT 格式：TF_name\tdescription\ttarget1\ttarget2\t...
"""

import gzip
import logging
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple

logger = logging.getLogger(__name__)


class ChEA3Parser:
    """ChEA3 数据库解析器

    解析 ChEA3 GMT 格式文件和 API 返回结果，
    生成 AllEnricher 标准格式的 ChEA3 数据库文件。

    输入文件格式 (GMT):
        TF_name\tdescription\ttarget1\ttarget2\t...

    API 返回结果格式:
        {lib_name: [{TF, Rank, Pvalue, Overlap, TargetCount}, ...]}

    输出文件格式：
    - {species}.ChEA3_2gene.tab.gz: Gene\tTF1\tTF2\t... (0/1 矩阵)
    - {species}.ChEA3_2disc.gz: TF\tlib_count\ttarget_count
    """

    @staticmethod
    def parse_gmt(gmt_path: str, library_name: str = "unknown"
                  ) -> Tuple[Dict[str, Set[str]], Dict[str, str]]:
        """解析 ChEA3 GMT 格式文件

        Args:
            gmt_path: GMT 文件路径
            library_name: 库名称，用于日志标识

        Returns:
            (tf_to_targets, tf_descriptions)
            - tf_to_targets: {TF: {target1, target2, ...}}
            - tf_descriptions: {TF: description}

        Raises:
            FileNotFoundError: 输入文件不存在时抛出
        """
        gmt_file = Path(gmt_path)
        if not gmt_file.exists():
            raise FileNotFoundError(f"[错误] GMT 文件不存在: {gmt_path}")

        tf_to_targets: Dict[str, Set[str]] = {}
        tf_descriptions: Dict[str, str] = {}

        open_func = gzip.open if gmt_path.endswith('.gz') else open
        mode = 'rt' if gmt_path.endswith('.gz') else 'r'

        with open_func(gmt_path, mode, encoding='utf-8') as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue

                parts = line.split('\t')
                if len(parts) < 3:
                    continue

                tf_name = parts[0]
                description = parts[1]
                targets = set(parts[2:])

                tf_to_targets[tf_name] = targets
                tf_descriptions[tf_name] = description

        logger.info("ChEA3Parser: 解析 GMT [%s] 完成, 共 %d 个 TF",
                    library_name, len(tf_to_targets))

        return tf_to_targets, tf_descriptions

    @staticmethod
    def parse_api_result(api_result: Dict[str, List[Dict[str, str]]]
                        ) -> Dict[str, List[Dict[str, str]]]:
        """解析 ChEA3 API 返回结果

        将 API 返回的原始结果标准化为统一格式。

        Args:
            api_result: API 返回的原始结果
                {lib_name: [{TF, Rank, Pvalue, Overlap, TargetCount}, ...]}

        Returns:
            标准化后的结果，格式与输入相同：
            {lib_name: [{TF, Rank, Pvalue, Overlap, TargetCount}, ...]}
        """
        standardized: Dict[str, List[Dict[str, str]]] = {}

        for lib_name, entries in api_result.items():
            standardized[lib_name] = []
            for entry in entries:
                standardized[lib_name].append({
                    'TF': str(entry.get('TF', '')),
                    'Rank': str(entry.get('Rank', '')),
                    'Pvalue': str(entry.get('Pvalue', '')),
                    'Overlap': str(entry.get('Overlap', '')),
                    'TargetCount': str(entry.get('TargetCount', '')),
                })

        return standardized

    @staticmethod
    def merge_libraries(libraries: Dict[str, Dict[str, Set[str]]],
                        method: str = "union") -> Dict[str, Set[str]]:
        """合并多个库的 TF-target 数据

        Args:
            libraries: {lib_name: {TF: {targets}}}
            method: 合并方法
                - "union": 取所有库的并集
                - "intersection": 取交集（TF 需在 >=2 个库中存在）

        Returns:
            {TF: {targets}}

        Raises:
            ValueError: method 参数不合法时抛出
        """
        if method not in ("union", "intersection"):
            raise ValueError(
                f"[错误] 不支持的合并方法: {method}，"
                f"可选值: 'union', 'intersection'"
            )

        if not libraries:
            return {}

        # 统计每个 TF 出现在多少个库中
        tf_lib_count: Dict[str, int] = {}
        for lib_name, tf_data in libraries.items():
            for tf in tf_data:
                tf_lib_count[tf] = tf_lib_count.get(tf, 0) + 1

        merged: Dict[str, Set[str]] = {}

        if method == "union":
            # 并集：所有 TF 的 targets 合并
            for lib_name, tf_data in libraries.items():
                for tf, targets in tf_data.items():
                    if tf not in merged:
                        merged[tf] = set()
                    merged[tf].update(targets)

        elif method == "intersection":
            # 交集：TF 需在 >=2 个库中存在，targets 取所有库的并集
            for tf, count in tf_lib_count.items():
                if count >= 2:
                    merged[tf] = set()
                    for lib_name, tf_data in libraries.items():
                        if tf in tf_data:
                            merged[tf].update(tf_data[tf])

        logger.info("ChEA3Parser: 合并完成 (method=%s), 共 %d 个 TF",
                    method, len(merged))

        return merged

    @staticmethod
    def build_database(gmt_paths: List[str], output_dir: str, species: str,
                        merge_method: str = "union",
                        valid_genes: Optional[Set[str]] = None) -> None:
        """构建 ChEA3 数据库文件

        解析多个 GMT 文件，合并后生成 AllEnricher 标准格式文件。

        Args:
            gmt_paths: GMT 文件路径列表
            output_dir: 输出目录
            species: 物种标识（如 "hsa"）
            merge_method: 合并方法 ("union" 或 "intersection")
            valid_genes: 可选的有效基因集合，用于过滤 targets

        Raises:
            FileNotFoundError: 输入文件不存在时抛出
            ValueError: 没有有效的 TF 数据时抛出
        """
        outdir_path = Path(output_dir)
        outdir_path.mkdir(parents=True, exist_ok=True)

        logger.info("ChEA3Parser: 开始构建 ChEA3 数据库 (species=%s)", species)

        # 第一步：解析所有 GMT 文件
        libraries: Dict[str, Dict[str, Set[str]]] = {}
        for gmt_path in gmt_paths:
            lib_name = Path(gmt_path).stem
            if lib_name.endswith('.gmt'):
                lib_name = lib_name[:-4]
            tf_to_targets, _ = ChEA3Parser.parse_gmt(gmt_path, lib_name)
            libraries[lib_name] = tf_to_targets

        if not libraries:
            raise ValueError("[错误] 没有解析到任何有效的 TF 数据！")

        # 第二步：合并多个库
        merged = ChEA3Parser.merge_libraries(libraries, method=merge_method)

        if not merged:
            raise ValueError("[错误] 合并后没有有效的 TF 数据！")

        # 第三步：过滤有效基因（如果提供）
        if valid_genes is not None:
            filtered_merged: Dict[str, Set[str]] = {}
            for tf, targets in merged.items():
                filtered_targets = targets & valid_genes
                if filtered_targets:
                    filtered_merged[tf] = filtered_targets
            merged = filtered_merged
            logger.info("基因过滤后剩余 %d 个 TF", len(merged))

        # 第四步：收集所有基因和 TF
        all_genes: Set[str] = set()
        for targets in merged.values():
            all_genes.update(targets)

        sorted_genes = sorted(all_genes)
        sorted_tfs = sorted(merged.keys())

        # 统计每个 TF 出现在多少个库中
        tf_lib_count: Dict[str, int] = {}
        for lib_name, tf_data in libraries.items():
            for tf in tf_data:
                if tf in merged:
                    tf_lib_count[tf] = tf_lib_count.get(tf, 0) + 1

        # 第五步：写入 ChEA3_2gene.tab.gz
        # 格式：Gene\tTF1\tTF2\t... (0/1 矩阵)
        tab_file = outdir_path / f"{species}.ChEA3_2gene.tab.gz"
        logger.info("写入文件: %s", tab_file)

        with gzip.open(tab_file, 'wt', encoding='utf-8') as f:
            header = ["Gene"] + sorted_tfs
            f.write('\t'.join(header) + '\n')

            for gene in sorted_genes:
                row = [gene]
                for tf in sorted_tfs:
                    if gene in merged[tf]:
                        row.append('1')
                    else:
                        row.append('0')
                f.write('\t'.join(row) + '\n')

        # 第六步：写入 ChEA3_2disc.gz
        # 格式：TF\tlib_count\ttarget_count
        disc_file = outdir_path / f"{species}.ChEA3_2disc.gz"
        logger.info("写入文件: %s", disc_file)

        with gzip.open(disc_file, 'wt', encoding='utf-8') as f:
            for tf in sorted_tfs:
                lib_count = tf_lib_count.get(tf, 0)
                target_count = len(merged[tf])
                f.write(f"{tf}\t{lib_count}\t{target_count}\n")

        logger.info("共 %d 个 TF, %d 个基因", len(sorted_tfs), len(sorted_genes))
        logger.info("ChEA3Parser: ChEA3 数据库构建完成")
