"""
KEGG 数据库解析器

解析已下载的 gene2pathway.txt 文件，生成 AllEnricher 标准的
kegg2gene.tab.gz 和 kegg2disc.gz 文件。

对应 v1 脚本：
- pathway2tab.pl: 从 pathway 基因列表生成 kegg2gene.tab
- makeDB.kegg.v1.1.sh: 生成 kegg2disc（Category|Subcategory|PathwayName 格式）

注意：v1 中 KEGG 数据需要实时抓取网页获取通路基因列表，
v2 中简化为接受已准备好的 gene2pathway.txt 文件。
"""

import gzip
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple


class KEGGParser:
    """KEGG 数据库解析器

    解析已准备好的 gene2pathway.txt 文件，生成 AllEnricher 标准格式的 KEGG 数据库文件。

    输入文件格式 (gene2pathway.txt):
        gene_symbol\\tentrez_id\\tpathway_id\\tpathway_name

    输入文件格式 (pathway_summary.txt, 可选):
        Category\\tSubcategory\\tpathway_id\\tpathway_name\\turl

    输出文件格式：
    - {species}.kegg2gene.tab.gz: Gene\\tpathway_id1\\tpathway_id2\\t... (0/1 矩阵)
    - {species}.kegg2disc.gz: pathway_id\\tCategory|Subcategory|PathwayName
    """

    @staticmethod
    def build_database(species: str, gene_info_path: str,
                       gene2pathway_path: str, outdir: str,
                       pathway_summary_path: Optional[str] = None) -> None:
        """构建 KEGG 数据库

        从 gene2pathway.txt 文件读取基因-通路关联，
        生成 kegg2gene.tab.gz 和 kegg2disc.gz。

        对应 v1 的 makeDB.kegg.v1.1.sh 中 pathway2tab.pl 和 kegg2disc 生成逻辑。

        Args:
            species: 物种缩写（如 hsa）
            gene_info_path: gene_info.gz 文件路径（用于验证基因）
            gene2pathway_path: gene2pathway.txt 文件路径
            outdir: 输出目录
            pathway_summary_path: pathway_summary.txt 文件路径（可选），
                格式为 Category\\tSubcategory\\tpathway_id\\tpathway_name\\turl
                如果不提供，kegg2disc 中的名称将使用 pathway_name

        Raises:
            FileNotFoundError: 输入文件不存在时抛出
            ValueError: 没有找到有效的基因-通路关联时抛出
        """
        outdir_path = Path(outdir)
        outdir_path.mkdir(parents=True, exist_ok=True)

        print(f"|--- KEGGParser: 开始构建 KEGG 数据库 (species={species})")

        # 第一步：读取 gene_info.gz，建立有效基因集合
        # v1 逻辑：过滤指定 taxid 的基因
        valid_genes: Set[str] = set()
        print(f"|--- 读取文件: {gene_info_path}")

        if gene_info_path.endswith('.gz'):
            f_open = gzip.open(gene_info_path, 'rt', encoding='utf-8')
        else:
            f_open = open(gene_info_path, 'r', encoding='utf-8')

        with f_open:
            for line in f_open:
                line = line.strip()
                if not line or line.startswith('#'):
                    continue
                parts = line.split('\t')
                if len(parts) >= 3:
                    valid_genes.add(parts[2])  # symbol 列

        print(f"|--- 找到 {len(valid_genes)} 个有效基因")

        # 第二步：读取 pathway_summary.txt（如果提供），建立通路分类信息
        # v1 逻辑：从 tmp03.{species}.kegg.xls 提取分类信息
        # 格式: Category\\tSubcategory\\tpathway_id\\tpathway_name\\turl
        pathway_categories: Dict[str, str] = {}  # {pathway_id: "Category|Subcategory|PathwayName"}

        if pathway_summary_path and Path(pathway_summary_path).exists():
            print(f"|--- 读取文件: {pathway_summary_path}")
            with open(pathway_summary_path, 'r', encoding='utf-8') as f:
                for line in f:
                    line = line.strip()
                    if not line or line.startswith('#'):
                        continue
                    parts = line.split('\t')
                    if len(parts) >= 4:
                        category = parts[0].replace(' ', '_')
                        subcategory = parts[1].replace(' ', '_')
                        pathway_id = parts[2]
                        pathway_name = parts[3].replace(' ', '_')
                        # 添加物种前缀
                        if not pathway_id.startswith(species):
                            pathway_id = f"{species}{pathway_id}"
                        pathway_categories[pathway_id] = (
                            f"{category}|{subcategory}|{pathway_name}"
                        )

        # 第三步：读取 gene2pathway.txt，构建基因-通路关联矩阵
        # v1 逻辑：pathway2tab.pl 从各个 pathway 的 glist.tab 文件读取基因列表
        # v2 简化：直接从 gene2pathway.txt 读取
        # 输入格式: gene_symbol\\tentrez_id\\tpathway_id\\tpathway_name
        all_pathways: Set[str] = set()
        all_genes: Set[str] = set()
        tab: Dict[str, Dict[str, int]] = {}  # {gene: {pathway_id: 1}}
        pathway_names: Dict[str, str] = {}   # {pathway_id: pathway_name}

        print(f"|--- 读取文件: {gene2pathway_path}")
        n = 0

        with open(gene2pathway_path, 'r', encoding='utf-8') as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith('#'):
                    continue
                parts = line.split('\t')
                if len(parts) < 4:
                    continue

                gene_symbol = parts[0]
                pathway_id = parts[2]
                pathway_name = parts[3]

                # 添加物种前缀（如 04110 -> hsa04110），与 v1 数据库格式一致
                if not pathway_id.startswith(species):
                    pathway_id = f"{species}{pathway_id}"

                # 如果 gene_info 中有该基因则使用，否则使用原始 symbol
                if gene_symbol in valid_genes:
                    symbol = gene_symbol
                else:
                    symbol = gene_symbol  # v1 中如果没有匹配则使用原始值

                all_genes.add(symbol)
                all_pathways.add(pathway_id)
                pathway_names[pathway_id] = pathway_name

                if symbol not in tab:
                    tab[symbol] = {}
                tab[symbol][pathway_id] = 1
                n += 1

        if n == 0:
            raise ValueError(
                "[错误] 在 gene2pathway.txt 文件中没有找到有效的基因-通路关联！"
            )

        print(f"|--- 共找到 {n} 条基因-通路关联")

        # 第四步：写入 kegg2gene.tab.gz
        # v1 逻辑：表头 Gene\\tpathway_id1\\tpathway_id2\\t...
        sorted_pathways = sorted(all_pathways)
        tab_file = outdir_path / f"{species}.kegg2gene.tab.gz"
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

        # 第五步：写入 kegg2disc.gz
        # v1 逻辑：格式为 pathway_id\\tCategory|Subcategory|PathwayName
        # 其中空格替换为下划线
        disc_file = outdir_path / f"{species}.kegg2disc.gz"
        print(f"|--- 写入文件: {disc_file}")

        with gzip.open(disc_file, 'wt', encoding='utf-8') as f:
            for pid in sorted_pathways:
                if pid in pathway_categories:
                    disc_name = pathway_categories[pid]
                else:
                    # 如果没有分类信息，使用 pathway_name（空格替换为下划线）
                    pname = pathway_names.get(pid, pid)
                    disc_name = pname.replace(' ', '_')
                f.write(f"{pid}\t{disc_name}\n")

        print(f"|--- KEGGParser: KEGG 数据库构建完成")
