"""
DisGeNET 数据库解析器

解析 DisGeNET all_gene_disease_associations.tsv.gz 文件，
生成 AllEnricher 标准的 CUI2gene.tab.gz 和 CUI2disc.gz 文件。

对应 v1 脚本：
- makeDB.DisGeNET.v1.0.sh: 下载 DisGeNET 文件，提取 gene-CUI 关联
- DisGeNET_gene_filter.pl: 使用 gene_info 过滤有效基因
- gene2CUI_extract.pl: 生成 CUI2gene.tab.gz 和 CUI2disc.gz
"""

import gzip
import re
from pathlib import Path
from typing import Dict, List, Optional, Set


class DisGeNETParser:
    """DisGeNET 数据库解析器

    解析 DisGeNET all_gene_disease_associations.tsv.gz 文件，
    生成 AllEnricher 标准格式的 DisGeNET 数据库文件。

    输入文件格式 (all_gene_disease_associations.tsv.gz):
        gene_symbol, ..., disease_id(CUI:xxx), disease_name, ...

    输出文件格式：
    - hsa.CUI2gene.tab.gz: Gene\\tCUI1\\tCUI2\\t... (0/1 矩阵)
    - hsa.CUI2disc.gz: CUI\\tdisease_name (空格和连字符替换为下划线)
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
    def parse_associations(assoc_path: str, gene_info_path: str,
                           taxid: int, outdir: str) -> None:
        """解析 DisGeNET 关联文件，生成 CUI 数据库文件

        读取 all_gene_disease_associations.tsv.gz，
        过滤 CUI 开头的行，用 gene_info 过滤有效基因。

        对应 v1 的 makeDB.DisGeNET.v1.0.sh + DisGeNET_gene_filter.pl + gene2CUI_extract.pl。

        v1 处理流程：
        1. 从 TSV 提取 gene_symbol, gene_id, CUI, disease_name
           列 [0]=gene_symbol, [1]=gene_id, [4]=CUI, [5]=disease_name
        2. 过滤 CUI: 开头的行（CUI 格式: Cxxxxxx）
        3. disease_name 空格和连字符替换为下划线
        4. 去除单引号
        5. 使用 gene_info 过滤有效基因
        6. 生成 CUI2gene.tab.gz 和 CUI2disc.gz

        Args:
            assoc_path: all_gene_disease_associations.tsv.gz 文件路径
            gene_info_path: gene_info.gz 文件路径
            taxid: 物种分类学 ID（如 9606）
            outdir: 输出目录

        Raises:
            FileNotFoundError: 输入文件不存在时抛出
            ValueError: 没有找到有效的 gene-disease 关联时抛出
        """
        outdir_path = Path(outdir)
        outdir_path.mkdir(parents=True, exist_ok=True)

        print(f"|--- DisGeNETParser: 开始解析 DisGeNET 文件 (taxid={taxid})")

        # 第一步：读取 gene_info.gz，建立有效基因集合
        # v1 逻辑：DisGeNET_gene_filter.pl 读取 gene_info，建立 symbol 映射
        valid_genes: Set[str] = set()
        print(f"|--- 读取文件: {gene_info_path}")

        with DisGeNETParser._open_gz_or_text(gene_info_path) as f:
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

        # 第二步：读取 DisGeNET 关联文件，提取 gene-CUI 关联
        # v1 逻辑：
        #   列 [0]=gene_symbol, [1]=gene_id, [4]=disease_id(CUI), [5]=disease_name
        #   过滤 CUI: 开头的行（CUI 格式: Cxxxxxx）
        #   disease_name 空格和连字符替换为下划线
        #   去除单引号
        all_cuis: Set[str] = set()
        cui_data: Dict[str, Dict[str, str]] = {}  # {symbol: {cui: disease_name}}
        all_symbols: Set[str] = set()
        n = 0

        print(f"|--- 读取文件: {assoc_path}")

        with DisGeNETParser._open_gz_or_text(assoc_path) as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith('#'):
                    continue

                # 去除单引号（v1 逻辑：sed "s/'//g"）
                line = line.replace("'", "")

                parts = line.split('\t')
                if len(parts) < 6:
                    continue

                gene_symbol = parts[0]
                gene_id = parts[1]
                disease_id = parts[4]
                disease_name = parts[5]

                # v1 逻辑：过滤 CUI: 开头的行（CUI 格式: Cxxxxxx）
                if not disease_id.startswith('CUI:'):
                    continue

                # 标准化 gene symbol（大写比较）
                symbol_upper = gene_symbol.upper()

                # v1 逻辑：使用 gene_info 过滤有效基因
                if symbol_upper not in valid_genes:
                    continue

                # v1 逻辑：disease_name 空格和连字符替换为下划线
                disease_name = disease_name.replace(' ', '_').replace('-', '_')

                all_cuis.add(disease_id)
                all_symbols.add(gene_symbol)

                if gene_symbol not in cui_data:
                    cui_data[gene_symbol] = {}
                cui_data[gene_symbol][disease_id] = disease_name
                n += 1

        if n == 0:
            raise ValueError(
                "[错误] 在 DisGeNET 文件中没有找到有效的 gene-disease 关联！"
            )

        print(f"|--- 共找到 {n} 条 gene-disease 注释")

        # 第三步：写入 CUI2gene.tab.gz
        # v1 逻辑：表头 Gene\\tCUI1\\tCUI2\\t...
        sorted_cuis = sorted(all_cuis)
        tab_file = outdir_path / "hsa.CUI2gene.tab.gz"
        print(f"|--- 写入文件: {tab_file}")

        # 收集实际存在的 CUI（有基因关联的）
        uniq_disc: Dict[str, str] = {}  # {cui: disease_name}

        with gzip.open(tab_file, 'wt', encoding='utf-8') as f:
            header = ["Gene"] + sorted_cuis
            f.write('\t'.join(header) + '\n')

            for symbol in sorted(all_symbols):
                row = [symbol]
                for cui in sorted_cuis:
                    if cui in cui_data.get(symbol, {}):
                        row.append('1')
                        uniq_disc[cui] = cui_data[symbol][cui]
                    else:
                        row.append('0')
                f.write('\t'.join(row) + '\n')

        # 第四步：写入 CUI2disc.gz
        # v1 逻辑：CUI\\tdisease_name，仅包含有基因关联的 CUI
        disc_file = outdir_path / "hsa.CUI2disc.gz"
        print(f"|--- 写入文件: {disc_file}")

        with gzip.open(disc_file, 'wt', encoding='utf-8') as f:
            for cui in sorted(uniq_disc.keys()):
                f.write(f"{cui}\t{uniq_disc[cui]}\n")

        print(f"|--- 共 {len(uniq_disc)} 个 CUI Term")
        print(f"|--- DisGeNETParser: DisGeNET 数据库构建完成")
