"""UniProt GOA proteomes 按需获取模块

通过 UniProt GOA proteomes FTP 获取物种 GO 注释数据（GAF 格式），
生成与 GOParser.parse_gene2go 兼容的 GO2gene.tab.gz 和 gene2go.txt 文件。

URL 模板：https://ftp.ebi.ac.uk/pub/databases/GO/goa/proteomes/{taxid}.goa

对应 v1 数据源：UniProt GOA (替代 NCBI gene2go.gz)
"""

from __future__ import annotations

import gzip
import logging
from pathlib import Path
from typing import Dict, Optional, Set, Tuple

import requests

logger = logging.getLogger(__name__)


class GOAFetcher:
    """UniProt GOA proteomes 按需获取器

    Usage::

        fetcher = GOAFetcher(cache_dir='./database/basic/goa')
        goa_file = fetcher.fetch_species_data(9606, 'Homo_sapiens')
        gene_to_go, all_genes = fetcher.parse_goa_file(goa_file, 9606)
    """

    BASE_URL = "https://ftp.ebi.ac.uk/pub/databases/GO/goa/proteomes"
    UA = "Mozilla/5.0 (Windows NT 10.0; Win64; x64)"
    TIMEOUT = 30

    def __init__(self, cache_dir: str, overwrite: bool = False):
        """
        Args:
            cache_dir: GOA 文件缓存目录
            overwrite: 是否覆盖已有缓存
        """
        self.cache_dir = Path(cache_dir)
        self.cache_dir.mkdir(parents=True, exist_ok=True)
        self.overwrite = overwrite

    def fetch_species_data(
        self,
        taxid: int,
        latin_name: str,
        goa_filename: Optional[str] = None,
    ) -> Path:
        """按 TaxID 下载物种的 GOA 文件

        Args:
            taxid: 物种分类学 ID（如 9606）
            latin_name: 物种拉丁名（如 Homo_sapiens）
            goa_filename: GOA 文件名，为 None 时自动构造

        Returns:
            本地缓存文件路径（.goa.gz）

        Raises:
            requests.RequestException: 下载失败
        """
        # EBI GOA 文件名格式为 {taxid}.goa（如 9606.goa）
        if goa_filename is None:
            goa_filename = f"{taxid}.goa"

        # 本地文件名与 URL 文件名保持一致
        local_file = self.cache_dir / f"{goa_filename}.gz"

        if local_file.exists() and not self.overwrite:
            logger.info("GOA 文件已缓存，跳过下载: %s", local_file)
            return local_file

        url = f"{self.BASE_URL}/{goa_filename}"
        logger.info("下载 GOA 文件: %s -> %s", url, local_file)

        resp = requests.get(
            url,
            headers={"User-Agent": self.UA},
            timeout=self.TIMEOUT,
        )
        resp.raise_for_status()

        local_file.write_bytes(resp.content)
        logger.info("GOA 文件下载完成: %s (%d bytes)", local_file, len(resp.content))

        return local_file

    def parse_goa_file(
        self,
        goa_file: Path,
        taxid: int,
    ) -> Tuple[Dict[str, Set[str]], Set[str]]:
        """解析 GOA (GAF 格式) 文件

        GAF 格式说明（制表符分隔，17 列）：
            第1列  DB             数据库名称（如 UniProtKB）
            第2列  DB_Object_ID   数据库对象 ID（如 P12345）
            第3列  DB_Object_Symbol Gene Symbol（如 GOT2）
            第5列  GO_ID          GO 标识符（如 GO:0006457）
            第7列  Evidence_Code  证据代码（如 IEA）
            第9列  Aspect         GO 类别：P(生物过程)/F(分子功能)/C(细胞组分)
            第13列 Taxon          分类学 ID（如 taxon:9606）

        Args:
            goa_file: GOA 文件路径（.goa 或 .goa.gz）
            taxid: 只保留该 taxid 的注释行

        Returns:
            (gene_to_go_terms, all_genes)
            - gene_to_go_terms: {gene_symbol: {GO_ID, ...}, ...}
            - all_genes: 所有基因 symbol 集合
        """
        if not goa_file.exists():
            raise FileNotFoundError(f"GOA 文件不存在: {goa_file}")

        opener: callable = gzip.open if str(goa_file).endswith(".gz") else open

        gene_to_go: Dict[str, Set[str]] = {}
        all_genes: Set[str] = set()
        skipped_symbol_eq_id = 0
        skipped_symbol_starts_digit = 0
        skipped_wrong_taxid = 0
        total_lines = 0

        with opener(goa_file, "rt", encoding="utf-8") as f:
            for line in f:
                line = line.rstrip("\n\r")
                # 跳过注释行
                if line.startswith("!"):
                    continue

                parts = line.split("\t")
                if len(parts) < 14:
                    continue

                total_lines += 1

                # 第13列: Taxon（可能包含 taxon:9606|taxon:xxxx 的格式）
                taxon_field = parts[13]
                expected_taxon = f"taxon:{taxid}"
                if taxon_field != expected_taxon:
                    skipped_wrong_taxid += 1
                    continue

                # 第2列: DB_Object_ID, 第3列: DB_Object_Symbol
                db_object_id = parts[1]
                symbol = parts[2]

                # 跳过 Symbol 等于 ID 的行
                if symbol == db_object_id:
                    skipped_symbol_eq_id += 1
                    continue

                # 跳过 Symbol 以数字开头的行
                if symbol and symbol[0].isdigit():
                    skipped_symbol_starts_digit += 1
                    continue

                # 第5列: GO_ID
                go_id = parts[4]
                if not go_id.startswith("GO:"):
                    continue

                gene_to_go.setdefault(symbol, set()).add(go_id)
                all_genes.add(symbol)

        logger.info(
            "GOA 解析完成: taxid=%d, 总行数=%d, 有效基因=%d, "
            "跳过(错误taxid)=%d, 跳过(Symbol=ID)=%d, 跳过(Symbol数字开头)=%d",
            taxid, total_lines, len(all_genes),
            skipped_wrong_taxid, skipped_symbol_eq_id, skipped_symbol_starts_digit,
        )

        return gene_to_go, all_genes

    @staticmethod
    def build_go2gene_matrix(
        gene_to_go: Dict[str, Set[str]],
        all_genes: Set[str],
        all_go_terms: Set[str],
        output_path: Path,
    ) -> None:
        """生成 GO2gene.tab.gz 文件

        输出格式与 GOParser.parse_gene2go 的 GO2gene.tab.gz 一致：
            第一行 header: Gene\\tGO:0000001\\tGO:0000002\\t...
            每行: gene_symbol\\t0/1\\t0/1\\t...

        Args:
            gene_to_go: {gene_symbol: {GO_ID, ...}, ...}
            all_genes: 所有基因 symbol 集合
            all_go_terms: 所有 GO ID 集合
            output_path: 输出文件路径（.tab.gz）
        """
        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)

        sorted_go = sorted(all_go_terms)
        sorted_genes = sorted(all_genes)

        with gzip.open(output_path, "wt", encoding="utf-8") as f:
            # 写入 header
            header = ["Gene"] + sorted_go
            f.write("\t".join(header) + "\n")

            # 写入数据行
            for symbol in sorted_genes:
                go_set = gene_to_go.get(symbol, set())
                row = [symbol] + ["1" if go_id in go_set else "0" for go_id in sorted_go]
                f.write("\t".join(row) + "\n")

        logger.info(
            "GO2gene.tab.gz 已生成: %s (%d genes x %d GO terms)",
            output_path, len(sorted_genes), len(sorted_go),
        )

    @staticmethod
    def build_gene2go_list(
        gene_to_go: Dict[str, Set[str]],
        go_names: Dict[str, str],
        output_path: Path,
    ) -> None:
        """生成 gene2go.txt 文件

        输出格式：
            每行: gene_symbol\\tGO_ID\\tEvidence\\tGO_Name

        Args:
            gene_to_go: {gene_symbol: {GO_ID, ...}, ...}
            go_names: {GO_ID: GO_name, ...}
            output_path: 输出文件路径（不压缩）
        """
        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)

        with open(output_path, "w", encoding="utf-8") as f:
            for symbol in sorted(gene_to_go):
                for go_id in sorted(gene_to_go[symbol]):
                    go_name = go_names.get(go_id, "")
                    f.write(f"{symbol}\t{go_id}\t\t{go_name}\n")

        total_annotations = sum(len(v) for v in gene_to_go.values())
        logger.info(
            "gene2go.txt 已生成: %s (%d genes, %d annotations)",
            output_path, len(gene_to_go), total_annotations,
        )
