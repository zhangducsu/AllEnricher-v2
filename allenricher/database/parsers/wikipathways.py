"""
WikiPathways 数据库解析器

解析 WikiPathways GMT 文件，生成 AllEnricher 标准的
WikiPathways2gene.tab.gz 和 WikiPathways2disc.gz 文件。

注意：WikiPathways GMT 格式使用 '/' 作为基因分隔符（不同于标准 GMT 的 tab 分隔）
"""

import gzip
from pathlib import Path
from typing import Dict, Optional, Set, Tuple


class WikiPathwaysParser:
    """WikiPathways 数据库解析器

    解析 WikiPathways GMT 文件，生成 AllEnricher 标准格式的 WikiPathways 数据库文件。

    输入文件格式 (GMT):
        WPID<TAB>Pathway Name<TAB>Gene1/Gene2/Gene3/...
        注意：基因之间使用 '/' 分隔，而不是标准 GMT 的 tab 分隔

    输出文件格式：
    - {species}.WikiPathways2gene.tab.gz: Gene\tpathway_id1\tpathway_id2\t... (0/1 矩阵)
    - {species}.WikiPathways2disc.gz: pathway_id\tpathway_name
    """

    @staticmethod
    def parse_gmt(gmt_path: str) -> Tuple[Dict[str, Set[str]], Dict[str, str]]:
        """解析 WikiPathways GMT 文件

        WikiPathways GMT 格式：WPID<TAB>Pathway Name<TAB>Gene1/Gene2/Gene3/...
        注意：基因之间使用 '/' 作为分隔符（与标准 GMT 不同）

        Args:
            gmt_path: GMT 文件路径

        Returns:
            Tuple[Dict[str, Set[str]], Dict[str, str]]:
                - gene_sets: {pathway_id: {gene1, gene2, ...}}
                - descriptions: {pathway_id: pathway_name}

        Raises:
            FileNotFoundError: 输入文件不存在时抛出
        """
        gmt_file = Path(gmt_path)
        if not gmt_file.exists():
            raise FileNotFoundError(f"GMT 文件不存在: {gmt_path}")

        gene_sets: Dict[str, Set[str]] = {}
        descriptions: Dict[str, str] = {}

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

                pathway_id = parts[0]
                pathway_name = parts[1]
                # WikiPathways 使用 '/' 分隔基因
                genes_str = parts[2]
                genes = set(genes_str.split('/')) if genes_str else set()

                # 过滤空基因
                genes = {g.strip() for g in genes if g.strip()}

                gene_sets[pathway_id] = genes
                descriptions[pathway_id] = pathway_name

        return gene_sets, descriptions

    @staticmethod
    def load_gene_id_mapping(gene_info_path: str, taxid: int) -> Dict[str, str]:
        """从 NCBI gene_info.gz 加载 NCBI Gene ID → Symbol 映射

        gene_info.gz 格式: tax_id\tGeneID\tSymbol\t...

        Args:
            gene_info_path: gene_info.gz 文件路径
            taxid: NCBI 物种分类学 ID，用于过滤特定物种

        Returns:
            Dict[str, str]: {ncbi_gene_id: gene_symbol} 映射字典

        Raises:
            FileNotFoundError: 输入文件不存在时抛出
        """
        gene_info_file = Path(gene_info_path)
        if not gene_info_file.exists():
            raise FileNotFoundError(f"gene_info 文件不存在: {gene_info_path}")

        id_mapping: Dict[str, str] = {}
        taxid_str = str(taxid)

        open_func = gzip.open if gene_info_path.endswith('.gz') else open
        mode = 'rt' if gene_info_path.endswith('.gz') else 'r'

        with open_func(gene_info_path, mode, encoding='utf-8') as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith('#'):
                    continue

                parts = line.split('\t')
                if len(parts) < 3:
                    continue

                # 检查 tax_id 是否匹配
                if parts[0] != taxid_str:
                    continue

                gene_id = parts[1]  # GeneID
                symbol = parts[2]   # Symbol

                if gene_id and symbol:
                    id_mapping[gene_id] = symbol

        return id_mapping

    @staticmethod
    def convert_ncbi_to_symbol(
        gene_sets: Dict[str, Set[str]],
        id_mapping: Dict[str, str]
    ) -> Dict[str, Set[str]]:
        """将 ncbigene:xxx 格式转换为 Symbol

        GMT 文件中的基因格式是 "ncbigene:1234/ncbigene:5678"（使用 / 分隔）
        需要提取数字 ID，通过 id_mapping 查找对应的 Symbol

        Args:
            gene_sets: {pathway_id: {gene1, gene2, ...}}，基因格式为 ncbigene:xxx
            id_mapping: {ncbi_gene_id: gene_symbol} 映射字典

        Returns:
            Dict[str, Set[str]]: {pathway_id: {symbol1, symbol2, ...}}，基因已转换为 Symbol
        """
        converted_sets: Dict[str, Set[str]] = {}

        for pathway_id, genes in gene_sets.items():
            converted_genes: Set[str] = set()

            for gene in genes:
                # 提取 ncbigene:xxx 中的数字 ID
                if gene.startswith('ncbigene:'):
                    ncbi_id = gene.split(':', 1)[1]
                else:
                    # 如果不是 ncbigene: 格式，保留原样
                    converted_genes.add(gene)
                    continue

                # 查找对应的 Symbol
                symbol = id_mapping.get(ncbi_id)
                if symbol:
                    converted_genes.add(symbol)
                else:
                    # 如果找不到映射，保留原始 ID 以便调试
                    converted_genes.add(gene)

            converted_sets[pathway_id] = converted_genes

        return converted_sets

    @staticmethod
    def build_database(
        gmt_path: str,
        output_dir: str,
        species: str,
        taxid: Optional[int] = None,
        gene_info_path: Optional[str] = None,
        valid_genes: Optional[Set[str]] = None
    ) -> None:
        """构建 WikiPathways 数据库

        从 GMT 文件读取基因-通路关联，生成 WikiPathways2gene.tab.gz 和 WikiPathways2disc.gz。

        Args:
            gmt_path: WikiPathways GMT 文件路径
            output_dir: 输出目录
            species: 物种缩写（如 hsa, mmu, rno）
            taxid: NCBI 物种分类学 ID（如 9606），用于 ID 转换
            gene_info_path: gene_info.gz 文件路径（可选），用于 ID 转换和获取有效基因列表
            valid_genes: 有效基因集合（可选），如果提供则直接用于过滤

        Raises:
            FileNotFoundError: 输入文件不存在时抛出
            ValueError: 没有找到有效的基因-通路关联时抛出
        """
        outdir_path = Path(output_dir)
        outdir_path.mkdir(parents=True, exist_ok=True)

        print(f"|--- WikiPathwaysParser: 开始构建 WikiPathways 数据库 (species={species})")

        # 第一步：解析 GMT 文件
        print(f"|--- 读取 GMT 文件: {gmt_path}")
        gene_sets, descriptions = WikiPathwaysParser.parse_gmt(gmt_path)

        if not gene_sets:
            raise ValueError("[错误] GMT 文件中没有找到有效的通路数据！")

        print(f"|--- 共找到 {len(gene_sets)} 个通路")

        # 第二步：如果提供了 gene_info_path 和 taxid，进行 NCBI Gene ID → Symbol 转换
        if gene_info_path and Path(gene_info_path).exists() and taxid is not None:
            print(f"|--- 加载 NCBI Gene ID 映射 (taxid={taxid})...")
            id_mapping = WikiPathwaysParser.load_gene_id_mapping(gene_info_path, taxid)
            print(f"|--- 找到 {len(id_mapping)} 个 Gene ID 映射")

            print(f"|--- 转换 Gene ID 到 Symbol...")
            gene_sets = WikiPathwaysParser.convert_ncbi_to_symbol(gene_sets, id_mapping)
            print(f"|--- ID 转换完成")

        # 第三步：获取有效基因集合（用于过滤）
        valid_gene_set: Set[str] = set()
        if valid_genes:
            valid_gene_set = valid_genes
            print(f"|--- 使用提供的有效基因集合: {len(valid_gene_set)} 个基因")
        elif gene_info_path and Path(gene_info_path).exists():
            print(f"|--- 从 gene_info 加载有效基因列表: {gene_info_path}")
            open_func = gzip.open if gene_info_path.endswith('.gz') else open
            mode = 'rt' if gene_info_path.endswith('.gz') else 'r'

            with open_func(gene_info_path, mode, encoding='utf-8') as f:
                for line in f:
                    line = line.strip()
                    if not line or line.startswith('#'):
                        continue
                    parts = line.split('\t')
                    if len(parts) >= 3:
                        valid_gene_set.add(parts[2])  # symbol 列

            print(f"|--- 找到 {len(valid_gene_set)} 个有效基因")

        # 第四步：构建基因-通路关联矩阵
        all_pathways: Set[str] = set(gene_sets.keys())
        all_genes: Set[str] = set()
        tab: Dict[str, Dict[str, int]] = {}  # {gene: {pathway_id: 1}}

        for pathway_id, genes in gene_sets.items():
            for gene in genes:
                # 如果提供了有效基因集合，则过滤
                if valid_gene_set and gene not in valid_gene_set:
                    continue

                all_genes.add(gene)
                if gene not in tab:
                    tab[gene] = {}
                tab[gene][pathway_id] = 1

        if not all_genes:
            raise ValueError("[错误] 没有找到有效的基因-通路关联！")

        print(f"|--- 共找到 {len(all_genes)} 个基因，{sum(len(v) for v in tab.values())} 条基因-通路关联")

        # 第五步：写入 WikiPathways2gene.tab.gz
        sorted_pathways = sorted(all_pathways)
        tab_file = outdir_path / f"{species}.WikiPathways2gene.tab.gz"
        print(f"|--- 写入文件: {tab_file}")

        with gzip.open(tab_file, 'wt', encoding='utf-8') as f:
            header = ["Gene"] + sorted_pathways
            f.write('\t'.join(header) + '\n')

            for gene in sorted(all_genes):
                row = [gene]
                for pid in sorted_pathways:
                    val = tab.get(gene, {}).get(pid, 0)
                    row.append(str(val))
                f.write('\t'.join(row) + '\n')

        # 第六步：写入 WikiPathways2disc.gz
        disc_file = outdir_path / f"{species}.WikiPathways2disc.gz"
        print(f"|--- 写入文件: {disc_file}")

        with gzip.open(disc_file, 'wt', encoding='utf-8') as f:
            for pid in sorted_pathways:
                pname = descriptions.get(pid, pid)
                # 将空格替换为下划线（与 Reactome 等保持一致）
                pname_underscore = pname.replace(' ', '_')
                f.write(f"{pid}\t{pname_underscore}\n")

        print(f"|--- WikiPathwaysParser: WikiPathways 数据库构建完成")
