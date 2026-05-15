"""
Disease Ontology (DO) 数据库解析器

解析 Jensen Lab 的 human_disease_*.tsv 文件，
生成 AllEnricher 标准的 DO2gene.tab.gz 和 DO2disc.gz 文件。

对应 v1 脚本：
- makeDB.do.v1.0.sh: 下载 DO 文件，提取 gene-DOID 关联，合并来源
- gene_filter.pl: 使用 gene_info 过滤有效基因
- gene2DO_extract.pl: 生成 DO2gene.tab.gz 和 DO2disc.gz
"""

import gzip
import re
from pathlib import Path
from typing import Dict, List, Optional, Set


class DOParser:
    """Disease Ontology 数据库解析器

    解析 Jensen Lab 的 human_disease_*.tsv 文件（textmining/knowledge/experiments），
    生成 AllEnricher 标准格式的 DO 数据库文件。

    输入文件格式 (human_disease_*.tsv):
        ..., gene_symbol, ..., DOID:xxxx, disease_name, ...

    输出文件格式：
    - hsa.DO2gene.tab.gz: Gene\\tDOID1\\tDOID2\\t... (0/1 矩阵)
    - hsa.DO2disc.gz: DOID\\tdisease_name (空格和连字符替换为下划线)
    """

    @staticmethod
    def _open_gz_or_text(filepath: str):
        """根据文件扩展名自动选择打开方式（gzip 或文本）

        Args:
            filepath: 文件路径

        Returns:
            文件对象
        """
        if filepath.endswith('.gz'):
            return gzip.open(filepath, 'rt', encoding='utf-8')
        else:
            return open(filepath, 'r', encoding='utf-8')

    @staticmethod
    def parse_disease_files(disease_files: List[str], gene_info_path: str,
                            taxid: int, outdir: str) -> None:
        """解析多个 human_disease_*.tsv 文件，生成 DO 数据库文件

        读取多个 Jensen Lab 的 human_disease_*.tsv 文件（textmining/knowledge/experiments），
        过滤 DOID 开头的行，合并所有来源，用 gene_info 过滤有效基因。

        对应 v1 的 makeDB.do.v1.0.sh + gene_filter.pl + gene2DO_extract.pl。

        v1 处理流程：
        1. 从 TSV 提取 gene_symbol, DOID, disease_name（过滤 DOID: 开头的行）
        2. disease_name 中的空格和连字符替换为下划线
        3. 合并 textmining/knowledge/experiments 三个来源，去重
        4. 使用 gene_info 过滤有效基因
        5. 生成 DO2gene.tab.gz 和 DO2disc.gz

        Args:
            disease_files: human_disease_*.tsv 文件路径列表
            gene_info_path: gene_info.gz 文件路径
            taxid: 物种分类学 ID（如 9606）
            outdir: 输出目录

        Raises:
            FileNotFoundError: 输入文件不存在时抛出
            ValueError: 没有找到有效的 gene-DO 关联时抛出
        """
        outdir_path = Path(outdir)
        outdir_path.mkdir(parents=True, exist_ok=True)

        print(f"|--- DOParser: 开始解析 Disease Ontology 文件 (taxid={taxid})")

        # 第一步：读取 gene_info.gz，建立有效基因集合
        # v1 逻辑：gene_filter.pl 读取 gene_info，建立 symbol 映射
        valid_genes: Set[str] = set()
        print(f"|--- 读取文件: {gene_info_path}")

        with DOParser._open_gz_or_text(gene_info_path) as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith('#'):
                    continue
                parts = line.split('\t')
                if len(parts) >= 3:
                    file_taxid = parts[0]
                    gene_id = parts[1]
                    symbol = parts[2]
                    if int(file_taxid) == taxid:
                        valid_genes.add(symbol.upper())

        print(f"|--- 找到 {len(valid_genes)} 个有效基因 (taxid={taxid})")

        # 第二步：读取所有 disease 文件，提取 gene-DOID 关联
        # v1 逻辑：
        #   列 [1]=gene_symbol, [2]=DOID, [3]=disease_name
        #   过滤 DOID: 开头的行
        #   disease_name 空格和连字符替换为下划线
        #   去除单引号
        all_doids: Set[str] = set()
        do_data: Dict[str, Dict[str, str]] = {}  # {symbol: {doid: disease_name}}
        all_symbols: Set[str] = set()
        n = 0

        for disease_file in disease_files:
            print(f"|--- 读取文件: {disease_file}")

            with DOParser._open_gz_or_text(disease_file) as f:
                for line in f:
                    line = line.strip()
                    if not line or line.startswith('#'):
                        continue

                    # 去除单引号（v1 逻辑：sed "s/'//g"）
                    line = line.replace("'", "")

                    parts = line.split('\t')
                    if len(parts) < 4:
                        continue

                    gene_symbol = parts[1]
                    doid = parts[2]
                    disease_name = parts[3]

                    # v1 逻辑：过滤 DOID: 开头的行
                    if not doid.startswith('DOID:'):
                        continue

                    # v1 逻辑：disease_name 空格和连字符替换为下划线
                    disease_name = disease_name.replace(' ', '_').replace('-', '_')

                    # 标准化 gene symbol（大写比较）
                    symbol_upper = gene_symbol.upper()

                    # v1 逻辑：使用 gene_info 过滤有效基因
                    if symbol_upper not in valid_genes:
                        continue

                    # 使用原始 symbol（保持大小写）
                    all_doids.add(doid)
                    all_symbols.add(gene_symbol)

                    if gene_symbol not in do_data:
                        do_data[gene_symbol] = {}
                    do_data[gene_symbol][doid] = disease_name
                    n += 1

        if n == 0:
            raise ValueError(
                "[错误] 在 Disease Ontology 文件中没有找到有效的 gene-DO 关联！"
            )

        print(f"|--- 共找到 {n} 条 gene-DO 注释")

        # 第三步：写入 DO2gene.tab.gz
        # v1 逻辑：表头 Gene\\tDOID1\\tDOID2\\t...
        sorted_doids = sorted(all_doids)
        tab_file = outdir_path / "hsa.DO2gene.tab.gz"
        print(f"|--- 写入文件: {tab_file}")

        # 收集实际存在的 DOID（有基因关联的）
        uniq_disc: Dict[str, str] = {}  # {doid: disease_name}

        with gzip.open(tab_file, 'wt', encoding='utf-8') as f:
            header = ["Gene"] + sorted_doids
            f.write('\t'.join(header) + '\n')

            for symbol in sorted(all_symbols):
                row = [symbol]
                for doid in sorted_doids:
                    if doid in do_data.get(symbol, {}):
                        row.append('1')
                        uniq_disc[doid] = do_data[symbol][doid]
                    else:
                        row.append('0')
                f.write('\t'.join(row) + '\n')

        # 第四步：写入 DO2disc.gz
        # v1 逻辑：DOID\\tdisease_name，仅包含有基因关联的 DOID
        disc_file = outdir_path / "hsa.DO2disc.gz"
        print(f"|--- 写入文件: {disc_file}")

        with gzip.open(disc_file, 'wt', encoding='utf-8') as f:
            for doid in sorted(uniq_disc.keys()):
                f.write(f"{doid}\t{uniq_disc[doid]}\n")

        print(f"|--- 共 {len(uniq_disc)} 个 DO Term")
        print(f"|--- DOParser: DO 数据库构建完成")
