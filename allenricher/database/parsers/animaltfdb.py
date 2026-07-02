"""
AnimalTFDB 4.0 数据库解析器

解析 AnimalTFDB 的 TF 列表和直系同源映射文件。

输入文件：
- {Species}_TF: TF 列表 (Species, Symbol, Ensembl, Family, Protein, Entrez_ID)
- {Species}_ortholog_to_human: 直系同源映射

输出格式：
- {species}.AnimalTFDB_2tf.tab.gz: TF 信息表
- {species}.AnimalTFDB_2disc.gz: TF 描述信息
- {species}.AnimalTFDB_ortholog.gz: 同源映射表
"""

import gzip
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple
from collections import defaultdict
import logging

import pandas as pd

logger = logging.getLogger(__name__)


class AnimalTFDBParser:
    """AnimalTFDB 4.0 数据库解析器

    解析 TF 列表和直系同源映射。
    """

    @staticmethod
    def parse_tf_list(tf_list_path: str) -> pd.DataFrame:
        """解析 TF 列表文件

        Args:
            tf_list_path: {Species}_TF 文件路径

        Returns:
            TF 信息 DataFrame，列: Species, Symbol, Ensembl, Family, Protein, Entrez_ID
        """
        df = pd.read_csv(tf_list_path, sep='\t', low_memory=False)

        # 标准化列名（去除可能的空格）
        df.columns = df.columns.str.strip()

        logger.info(f"AnimalTFDB TF列表: {len(df)} 个 TF")
        return df

    @staticmethod
    def parse_ortholog_to_human(ortholog_path: str) -> Dict[str, str]:
        """解析直系同源映射文件

        Args:
            ortholog_path: {Species}_ortholog_to_human 文件路径

        Returns:
            {物种基因Symbol: 人类基因Symbol} 映射字典
        """
        ortholog_map: Dict[str, str] = {}

        with open(ortholog_path, 'r') as f:
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) >= 2:
                    species_gene = parts[0].strip()
                    human_gene = parts[1].strip()
                    if species_gene and human_gene:
                        ortholog_map[species_gene] = human_gene

        logger.info(f"AnimalTFDB 同源映射: {len(ortholog_map)} 对")
        return ortholog_map

    @staticmethod
    def build_database(
        tf_list_path: str,
        ortholog_path: str,
        output_dir: str,
        species: str,
        valid_genes: Optional[Set[str]] = None
    ) -> Tuple[pd.DataFrame, Dict[str, str]]:
        """构建 AnimalTFDB 数据库

        Args:
            tf_list_path: TF 列表文件路径
            ortholog_path: 同源映射文件路径
            output_dir: 输出目录
            species: 物种代码（如 bta 代表牛）
            valid_genes: 有效基因集合（可选）

        Returns:
            (tf_df, ortholog_map) 用于后续同源映射
        """
        outdir = Path(output_dir)
        outdir.mkdir(parents=True, exist_ok=True)

        logger.info(f"AnimalTFDBParser: 开始构建数据库 (species={species})")

        # 解析 TF 列表
        tf_df = AnimalTFDBParser.parse_tf_list(tf_list_path)

        # 过滤有效基因
        if valid_genes:
            tf_df = tf_df[tf_df['Symbol'].isin(valid_genes)]
            logger.info(f"过滤后: {len(tf_df)} 个 TF")

        # 解析同源映射
        ortholog_map = AnimalTFDBParser.parse_ortholog_to_human(ortholog_path)

        # 保存 TF 信息表
        tf_file = outdir / f"{species}.AnimalTFDB_2tf.tab.gz"
        logger.info(f"写入文件: {tf_file}")

        with gzip.open(tf_file, 'wt') as f:
            f.write('\t'.join(tf_df.columns) + '\n')
            for _, row in tf_df.iterrows():
                f.write('\t'.join(str(v) for v in row.values) + '\n')

        # 保存 TF 描述文件
        disc_file = outdir / f"{species}.AnimalTFDB_2disc.gz"
        logger.info(f"写入文件: {disc_file}")

        with gzip.open(disc_file, 'wt') as f:
            f.write("TF\tFamily\tEntrez_ID\tEnsembl\tsource\n")
            for _, row in tf_df.iterrows():
                symbol = row.get('Symbol', '')
                family = row.get('Family', 'Unknown')
                entrez = row.get('Entrez_ID', 'NA')
                ensembl = row.get('Ensembl', 'NA')
                f.write(f"{symbol}\t{family}\t{entrez}\t{ensembl}\tAnimalTFDB\n")

        # 保存同源映射
        ortholog_file = outdir / f"{species}.AnimalTFDB_ortholog.gz"
        logger.info(f"写入文件: {ortholog_file}")

        with gzip.open(ortholog_file, 'wt') as f:
            f.write("Species_Gene\tHuman_Gene\n")
            for sp_gene, hu_gene in ortholog_map.items():
                f.write(f"{sp_gene}\t{hu_gene}\n")

        logger.info(f"AnimalTFDBParser: 数据库构建完成")
        return tf_df, ortholog_map
