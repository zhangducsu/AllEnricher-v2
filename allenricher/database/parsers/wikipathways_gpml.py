"""
WikiPathways GPML 解析器

解析 WikiPathways GPML ZIP 文件，从 GPML XML 中提取基因信息，
生成 AllEnricher 标准的 WikiPathways2gene.tab.gz 和 WikiPathways2disc.gz 文件。

GPML 文件结构:
- Pathway 元素包含 Name 属性（通路名称）
- DataNode 元素包含基因信息:
  - TextLabel 属性: 基因名称/标签
  - Xref 子元素: Database="Entrez Gene", ID="NCBI Gene ID"
"""

import gzip
import xml.etree.ElementTree as ET
import zipfile
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple


class WikiPathwaysGPMLParser:
    """WikiPathways GPML 解析器

    解析 WikiPathways GPML ZIP 文件，从 XML 中提取基因信息，
    生成 AllEnricher 标准格式的 WikiPathways 数据库文件。

    输入文件格式:
        GPML ZIP 文件，包含多个 .gpml 文件

    输出文件格式：
    - {species}.WikiPathways2gene.tab.gz: Gene\tpathway_id1\tpathway_id2\t... (0/1 矩阵)
    - {species}.WikiPathways2disc.gz: pathway_id\tpathway_name
    """

    GPML_NS = "http://pathvisio.org/GPML/2013a"

    def parse_gpml_zip(
        self,
        gpml_zip_path: str,
        gene_info_path: Optional[str] = None,
        taxid: Optional[str] = None
    ) -> Tuple[Dict[str, Set[str]], Dict[str, str]]:
        """解析 GPML ZIP 文件

        从 GPML ZIP 文件中提取所有通路的基因信息。

        Args:
            gpml_zip_path: GPML ZIP 文件路径
            gene_info_path: NCBI gene_info.gz 文件路径（可选），用于 Gene ID 到 Symbol 的转换
            taxid: 物种 Taxonomy ID（可选），用于过滤 gene_info

        Returns:
            Tuple[Dict[str, Set[str]], Dict[str, str]]:
                - gene_sets: {pathway_id: {gene1, gene2, ...}}
                - descriptions: {pathway_id: pathway_name}

        Raises:
            FileNotFoundError: 输入文件不存在时抛出
            zipfile.BadZipFile: ZIP 文件损坏时抛出
        """
        zip_path = Path(gpml_zip_path)
        if not zip_path.exists():
            raise FileNotFoundError(f"GPML ZIP 文件不存在: {gpml_zip_path}")

        # 加载 Gene ID 到 Symbol 的映射（如果提供了 gene_info）
        id_mapping: Dict[str, str] = {}
        if gene_info_path and taxid:
            id_mapping = self._load_gene_id_mapping(gene_info_path, taxid)

        gene_sets: Dict[str, Set[str]] = {}
        descriptions: Dict[str, str] = {}

        with zipfile.ZipFile(gpml_zip_path, 'r') as zf:
            for filename in zf.namelist():
                if not filename.endswith('.gpml'):
                    continue

                # 从文件名提取 pathway_id (例如: WP1234.gpml -> WP1234)
                pathway_id = Path(filename).stem

                try:
                    with zf.open(filename) as f:
                        tree = ET.parse(f)
                        root = tree.getroot()

                        # 获取通路名称
                        pathway_name = root.get('Name', pathway_id)
                        descriptions[pathway_id] = pathway_name

                        # 提取基因
                        genes = self._extract_genes_from_gpml(root, id_mapping)
                        gene_sets[pathway_id] = genes

                except ET.ParseError as e:
                    print(f"警告: 解析 {filename} 失败: {e}")
                    continue

        return gene_sets, descriptions

    def _extract_genes_from_gpml(
        self,
        root: ET.Element,
        id_mapping: Dict[str, str]
    ) -> Set[str]:
        """从 GPML XML 根元素中提取基因

        解析 DataNode 元素，提取基因信息。优先使用 NCBI Gene ID 映射到 Symbol，
        如果没有映射则使用 TextLabel。

        Args:
            root: GPML XML 根元素
            id_mapping: Gene ID 到 Symbol 的映射字典

        Returns:
            Set[str]: 基因 Symbol 集合
        """
        genes: Set[str] = set()
        ns = {'gpml': self.GPML_NS}

        # 查找所有 DataNode 元素
        for datanode in root.findall('.//gpml:DataNode', ns):
            # 获取 TextLabel
            text_label = datanode.get('TextLabel', '').strip()

            # 查找 Xref 子元素
            xref = datanode.find('gpml:Xref', ns)
            if xref is not None:
                database = xref.get('Database', '')
                gene_id = xref.get('ID', '')

                # 如果是 Entrez Gene 且提供了映射
                if database == 'Entrez Gene' and gene_id and gene_id in id_mapping:
                    genes.add(id_mapping[gene_id])
                elif text_label:
                    # 否则使用 TextLabel
                    genes.add(text_label)
            elif text_label:
                # 没有 Xref 时使用 TextLabel
                genes.add(text_label)

        return genes

    def _load_gene_id_mapping(
        self,
        gene_info_path: str,
        taxid: str
    ) -> Dict[str, str]:
        """加载 NCBI Gene ID 到 Symbol 的映射

        从 NCBI gene_info.gz 文件中加载指定物种的 Gene ID 到 Symbol 映射。

        Args:
            gene_info_path: gene_info.gz 文件路径
            taxid: 物种 Taxonomy ID

        Returns:
            Dict[str, str]: {gene_id: symbol}

        Raises:
            FileNotFoundError: 文件不存在时抛出
        """
        gene_info_file = Path(gene_info_path)
        if not gene_info_file.exists():
            raise FileNotFoundError(f"gene_info 文件不存在: {gene_info_path}")

        mapping: Dict[str, str] = {}

        open_func = gzip.open if gene_info_path.endswith('.gz') else open
        mode = 'rt' if gene_info_path.endswith('.gz') else 'r'

        with open_func(gene_info_path, mode, encoding='utf-8') as f:
            # 跳过标题行
            header = f.readline()

            for line in f:
                line = line.strip()
                if not line:
                    continue

                parts = line.split('\t')
                if len(parts) < 3:
                    continue

                # gene_info 格式: tax_id GeneID Symbol ...
                file_taxid = parts[0]
                gene_id = parts[1]
                symbol = parts[2]

                # 只加载指定物种
                if file_taxid == taxid:
                    mapping[gene_id] = symbol

        return mapping

    def build_database_from_gpml(
        self,
        gpml_zip_path: str,
        output_dir: str,
        species: str,
        taxid: Optional[str] = None,
        gene_info_path: Optional[str] = None
    ) -> Tuple[str, str]:
        """从 GPML ZIP 文件构建 WikiPathways 数据库

        解析 GPML ZIP 文件，生成 WikiPathways2gene.tab.gz 和 WikiPathways2disc.gz。

        Args:
            gpml_zip_path: GPML ZIP 文件路径
            output_dir: 输出目录
            species: 物种缩写（如 hsa, mmu, rno）
            taxid: 物种 Taxonomy ID（可选），用于 Gene ID 过滤
            gene_info_path: gene_info.gz 文件路径（可选），用于 Gene ID 到 Symbol 转换

        Returns:
            Tuple[str, str]: (tab_file_path, disc_file_path)

        Raises:
            FileNotFoundError: 输入文件不存在时抛出
            ValueError: 没有找到有效的基因-通路关联时抛出
        """
        outdir_path = Path(output_dir)
        outdir_path.mkdir(parents=True, exist_ok=True)

        print(f"|--- WikiPathwaysGPMLParser: 开始构建 WikiPathways 数据库 (species={species})")

        # 第一步：解析 GPML ZIP 文件
        print(f"|--- 读取 GPML ZIP 文件: {gpml_zip_path}")
        gene_sets, descriptions = self.parse_gpml_zip(
            gpml_zip_path,
            gene_info_path=gene_info_path,
            taxid=taxid
        )

        if not gene_sets:
            raise ValueError("[错误] GPML ZIP 文件中没有找到有效的通路数据！")

        print(f"|--- 共找到 {len(gene_sets)} 个通路")

        # 第二步：构建基因-通路关联矩阵
        all_pathways: Set[str] = set(gene_sets.keys())
        all_genes: Set[str] = set()
        tab: Dict[str, Dict[str, int]] = {}  # {gene: {pathway_id: 1}}

        for pathway_id, genes in gene_sets.items():
            for gene in genes:
                all_genes.add(gene)
                if gene not in tab:
                    tab[gene] = {}
                tab[gene][pathway_id] = 1

        if not all_genes:
            raise ValueError("[错误] 没有找到有效的基因-通路关联！")

        print(f"|--- 共找到 {len(all_genes)} 个基因，{sum(len(v) for v in tab.values())} 条基因-通路关联")

        # 第三步：写入 WikiPathways2gene.tab.gz
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

        # 第四步：写入 WikiPathways2disc.gz
        disc_file = outdir_path / f"{species}.WikiPathways2disc.gz"
        print(f"|--- 写入文件: {disc_file}")

        with gzip.open(disc_file, 'wt', encoding='utf-8') as f:
            for pid in sorted_pathways:
                pname = descriptions.get(pid, pid)
                # 将空格替换为下划线（与其他解析器保持一致）
                pname_underscore = pname.replace(' ', '_')
                f.write(f"{pid}\t{pname_underscore}\n")

        print(f"|--- WikiPathwaysGPMLParser: WikiPathways 数据库构建完成")

        return str(tab_file), str(disc_file)
