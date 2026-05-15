"""
Reactome 数据库解析器

解析 Reactome NCBI2Reactome_All_Levels.txt.gz 文件，
生成 AllEnricher 标准的 Reactome2gene.tab.gz 和 Reactome2disc.gz 文件。

对应 v1 脚本：
- gene2ReactomePathway_extract.pl: 解析 NCBI2Reactome，生成 gene2pathway.txt 和 Reactome2gene.tab.gz
- makeDB.reactome.v1.0.sh: 生成 Reactome2disc（pathway_name 空格替换为下划线）
"""

import gzip
from pathlib import Path
from typing import Dict, Optional, Set


class ReactomeParser:
    """Reactome 数据库解析器

    解析 NCBI2Reactome_All_Levels.txt.gz 文件，
    生成 AllEnricher 标准格式的 Reactome 数据库文件。

    输入文件格式 (NCBI2Reactome_All_Levels.txt.gz):
        geneid\\tpathway_id\\tpathway_name\\turl\\t...

    输出文件格式：
    - {species}.gene2pathway.txt: gene_symbol\\tgene_id\\tpathway_id\\tpathway_name
    - {species}.Reactome2gene.tab.gz: Gene\\tpathway_id1\\tpathway_id2\\t... (0/1 矩阵)
    - {species}.Reactome2disc.gz: pathway_id\\tpathway_name (空格替换为下划线)
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
    def parse_ncbi2reactome(ncbi2reactome_path: str, gene_info_path: str,
                            taxid: int, species: str, outdir: str) -> None:
        """解析 NCBI2Reactome_All_Levels.txt.gz，生成 Reactome 数据库文件

        读取 NCBI2Reactome_All_Levels.txt.gz（tab分隔，列: geneid, pathway_id, pathway_name, url, ...）
        pathway_id 格式: "R-HSA-12345"，提取物种代码（大写）与 species 比较
        读取 gene_info.gz 获取 geneid -> symbol 映射
        输出 {species}.Reactome2gene.tab.gz 和 {species}.Reactome2disc.gz

        对应 v1 的 gene2ReactomePathway_extract.pl 和 makeDB.reactome.v1.0.sh。

        Args:
            ncbi2reactome_path: NCBI2Reactome_All_Levels.txt.gz 文件路径
            gene_info_path: gene_info.gz 文件路径
            taxid: 物种分类学 ID（如 9606）
            species: 物种缩写（如 hsa）
            outdir: 输出目录

        Raises:
            FileNotFoundError: 输入文件不存在时抛出
            ValueError: 没有找到指定物种的 Reactome 注释时抛出
        """
        outdir_path = Path(outdir)
        outdir_path.mkdir(parents=True, exist_ok=True)

        print(f"|--- ReactomeParser: 开始解析 NCBI2Reactome "
              f"(taxid={taxid}, species={species})")

        # 第一步：读取 gene_info.gz，建立 geneid -> symbol 映射
        # v1 逻辑：读取 gene_info，过滤指定 taxid，提取 (geneid, symbol)
        gene_id_to_symbol: Dict[str, str] = {}
        print(f"|--- 读取文件: {gene_info_path}")

        with ReactomeParser._open_gz_or_text(gene_info_path) as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith('#'):
                    continue
                parts = line.split('\t')
                file_taxid = parts[0]
                gene_id = parts[1]
                symbol = parts[2]
                if int(file_taxid) == taxid:
                    gene_id_to_symbol[gene_id] = symbol

        print(f"|--- 找到 {len(gene_id_to_symbol)} 个基因 (taxid={taxid})")

        # 第二步：读取 NCBI2Reactome 文件，过滤指定物种
        # v1 逻辑：
        #   列 [0]=geneid, [1]=pathway_id, [3]=pathway_name
        #   pathway_id 格式 "R-HSA-12345"，取 bb[1]（大写）与 uc(species) 比较
        species_upper = species.upper()
        gene2pathway_file = outdir_path / f"{species}.gene2pathway.txt"
        tab_file = outdir_path / f"{species}.Reactome2gene.tab.gz"

        all_pathways: Set[str] = set()
        all_symbols: Set[str] = set()
        tab: Dict[str, Dict[str, int]] = {}  # {symbol: {pathway_id: 1}}
        pathway_names: Dict[str, str] = {}   # {pathway_id: pathway_name}

        print(f"|--- 读取文件: {ncbi2reactome_path}")
        n = 0

        with ReactomeParser._open_gz_or_text(ncbi2reactome_path) as f_in, \
             open(gene2pathway_file, 'w', encoding='utf-8') as f_txt:
            for line in f_in:
                line = line.strip()
                if not line or line.startswith('#'):
                    continue
                parts = line.split('\t')
                gene_id = parts[0]
                pathway_id = parts[1]
                pathway_name = parts[3]

                # v1 逻辑：pathway_id 格式 "R-HSA-12345"，取 bb[1] 与 uc(name) 比较
                pathway_parts = pathway_id.split('-')
                if len(pathway_parts) < 2:
                    continue
                pathway_species_code = pathway_parts[1]
                if pathway_species_code != species_upper:
                    continue

                # 获取 gene symbol，如果没有则使用 gene_id
                symbol = gene_id_to_symbol.get(gene_id, gene_id)

                all_symbols.add(symbol)
                all_pathways.add(pathway_id)
                pathway_names[pathway_id] = pathway_name

                # 写入 gene2pathway.txt
                f_txt.write(f"{symbol}\t{gene_id}\t{pathway_id}\t{pathway_name}\n")

                # 构建 tab 矩阵数据
                if symbol not in tab:
                    tab[symbol] = {}
                tab[symbol][pathway_id] = 1
                n += 1

        if n == 0:
            raise ValueError(
                f"[错误] 在 NCBI2Reactome 文件中没有找到 species={species} 的通路注释！"
            )

        print(f"|--- 共找到 {n} 条基因-通路关联")

        # 第三步：写入 Reactome2gene.tab.gz
        # v1 逻辑：表头 Gene\\tpathway_id1\\tpathway_id2\\t...
        sorted_pathways = sorted(all_pathways)
        print(f"|--- 写入文件: {tab_file}")

        with gzip.open(tab_file, 'wt', encoding='utf-8') as f:
            header = ["Gene"] + sorted_pathways
            f.write('\t'.join(header) + '\n')

            for symbol in sorted(all_symbols):
                row = [symbol]
                for pid in sorted_pathways:
                    val = tab.get(symbol, {}).get(pid, 0)
                    row.append(str(val))
                f.write('\t'.join(row) + '\n')

        # 第四步：写入 Reactome2disc.gz
        # v1 逻辑：从 gene2pathway.txt 提取 pathway_id 和 pathway_name
        # pathway_name 空格替换为下划线，sort | uniq
        disc_file = outdir_path / f"{species}.Reactome2disc.gz"
        print(f"|--- 写入文件: {disc_file}")

        with gzip.open(disc_file, 'wt', encoding='utf-8') as f:
            for pid in sorted_pathways:
                pname = pathway_names.get(pid, pid)
                # v1 逻辑：空格替换为下划线
                pname_underscore = pname.replace(' ', '_')
                f.write(f"{pid}\t{pname_underscore}\n")

        print(f"|--- ReactomeParser: Reactome 数据库构建完成")
