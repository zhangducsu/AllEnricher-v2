"""
GO 数据库解析器

解析 NCBI gene2go.gz 和 Gene Ontology go-basic.obo 文件，
生成 AllEnricher 标准的 GO2gene.tab.gz 和 GO2disc.gz 文件。

对应 v1 脚本：
- gene2GO_extract.pl: 解析 gene2go.gz，生成 gene2go.txt 和 GO2gene.tab.gz
- obo2go.pl: 解析 go-basic.obo，生成 GO2disc.gz
"""

import gzip
import re
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple


class GOParser:
    """GO 数据库解析器

    解析 NCBI gene2go.gz 和 go-basic.obo 文件，
    生成 AllEnricher 标准格式的 GO 数据库文件。

    输出文件格式：
    - {species}.gene2go.txt: gene_symbol\\tgene_id\\tGO_ID\\tcategory\\tGO_name
    - {species}.GO2gene.tab.gz: Gene\\tGO_ID1\\tGO_ID2\\t... (0/1 矩阵)
    - {species}.GO2disc.gz: GO_ID\\tnamespace:name\\tfather1;father2;...
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
    def parse_gene2go(gene2go_path: str, gene_info_path: str,
                      taxid: int, species: str, outdir: str) -> None:
        """解析 gene2go.gz 和 gene_info.gz，生成 GO2gene.tab.gz

        读取 gene2go.gz（tab分隔，列: taxid, geneid, go_id, ..., go_name, ..., category）
        读取 gene_info.gz（tab分隔，列: taxid, geneid, symbol, ...）
        过滤指定 taxid，输出 {species}.GO2gene.tab.gz 和 {species}.gene2go.txt

        对应 v1 的 gene2GO_extract.pl 脚本。

        Args:
            gene2go_path: gene2go.gz 文件路径
            gene_info_path: gene_info.gz 文件路径
            taxid: 物种分类学 ID（如 9606）
            species: 物种缩写（如 hsa）
            outdir: 输出目录

        Raises:
            FileNotFoundError: 输入文件不存在时抛出
            ValueError: 没有找到指定 taxid 的 GO 注释时抛出
        """
        outdir_path = Path(outdir)
        outdir_path.mkdir(parents=True, exist_ok=True)

        print(f"|--- GOParser: 开始解析 gene2go (taxid={taxid}, species={species})")

        # 第一步：读取 gene_info.gz，建立 geneid -> symbol 映射
        # v1 逻辑：读取 gene_info，过滤指定 taxid，提取 (geneid, symbol)
        gene_id_to_symbol: Dict[str, str] = {}
        print(f"|--- 读取文件: {gene_info_path}")

        with GOParser._open_gz_or_text(gene_info_path) as f:
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

        # 第二步：读取 gene2go.gz，过滤指定 taxid
        # v1 逻辑：列 [0]=taxid, [1]=geneid, [2]=go_id, [5]=go_name, [7]=category
        # 2026 新格式新增 Qualifier 列: [0]=taxid, [1]=geneid, [2]=go_id,
        #   [3]=Evidence, [4]=Qualifier, [5]=GO_term, [6]=PubMed, [7]=Category
        # 兼容两种格式：通过列数自动判断
        # 输出 gene2go.txt: symbol\\tgeneid\\tgo_id\\tcategory\\tgo_name
        gene2go_file = outdir_path / f"{species}.gene2go.txt"
        tab_file = outdir_path / f"{species}.GO2gene.tab.gz"

        all_go: Set[str] = set()       # 所有 GO ID
        all_symbols: Set[str] = set()  # 所有基因 symbol
        tab: Dict[str, Dict[str, int]] = {}  # {symbol: {go_id: 1}}

        print(f"|--- 读取文件: {gene2go_path}")
        n = 0

        with GOParser._open_gz_or_text(gene2go_path) as f_in, \
             open(gene2go_file, 'w', encoding='utf-8') as f_txt:
            for line in f_in:
                line = line.strip()
                if not line or line.startswith('#'):
                    continue
                parts = line.split('\t')
                file_taxid = parts[0]
                gene_id = parts[1]
                go_id = parts[2]

                # 兼容新旧格式：
                # 新格式 (2026+): taxid GeneID GO_ID Evidence Qualifier GO_term PubMed Category
                #   parts[4] 为 Qualifier (enables/involved_in/located_in/part_of)
                # 旧格式: taxid GeneID GO_ID Evidence GO_term PubMed Category
                #   parts[4] 为 GO_term (长文本描述)
                _QUALIFIERS = {'enables', 'involved_in', 'located_in', 'part_of', 'acts_upstream_of',
                               'acts_upstream_of_negative_effect', 'acts_upstream_of_positive_effect',
                               'colocalizes_with', 'contributes_to'}
                if len(parts) >= 8 and parts[4] in _QUALIFIERS:
                    # 新格式
                    go_name = parts[5]
                    category = parts[7]
                else:
                    # 旧格式
                    go_name = parts[4] if len(parts) > 4 else ""
                    category = parts[6] if len(parts) > 6 else ""

                if int(file_taxid) != taxid:
                    continue

                # 获取 gene symbol，如果没有则使用 gene_id
                symbol = gene_id_to_symbol.get(gene_id, gene_id)

                all_symbols.add(symbol)
                all_go.add(go_id)

                # 写入 gene2go.txt
                f_txt.write(f"{symbol}\t{gene_id}\t{go_id}\t{category}\t{go_name}\n")

                # 构建 tab 矩阵数据
                if symbol not in tab:
                    tab[symbol] = {}
                tab[symbol][go_id] = 1
                n += 1

        if n == 0:
            raise ValueError(
                f"[错误] 在 NCBI gene2go.gz 文件中没有找到 taxid={taxid} 的 GO 注释信息！"
            )

        print(f"|--- 共找到 {n} 条 gene2go 注释")

        # 第三步：写入 GO2gene.tab.gz
        # v1 逻辑：表头 Gene\\tGO_ID1\\tGO_ID2\\t...
        # 数据行 symbol\\t0/1\\t0/1\\t...
        sorted_go = sorted(all_go)
        print(f"|--- 写入文件: {tab_file}")

        with gzip.open(tab_file, 'wt', encoding='utf-8') as f:
            # 写入表头
            header = ["Gene"] + sorted_go
            f.write('\t'.join(header) + '\n')

            # 写入数据行
            for symbol in sorted(all_symbols):
                row = [symbol]
                for go_id in sorted_go:
                    val = tab.get(symbol, {}).get(go_id, 0)
                    row.append(str(val))
                f.write('\t'.join(row) + '\n')

        print(f"|--- GOParser: gene2go 解析完成")

    @staticmethod
    def parse_obo(obo_path: str, outdir: str) -> None:
        """解析 go-basic.obo 文件，生成 GO2disc.gz

        解析 OBO 文件中的 [Term] 条目，提取 GO ID、name、namespace 和 is_a 关系。
        输出格式: GO_ID\\tnamespace:name\\tfather1;father2;...

        对应 v1 的 obo2go.pl 脚本。

        Args:
            obo_path: go-basic.obo 文件路径
            outdir: 输出目录

        Raises:
            FileNotFoundError: 输入文件不存在时抛出
        """
        outdir_path = Path(outdir)
        outdir_path.mkdir(parents=True, exist_ok=True)

        print(f"|--- GOParser: 开始解析 go-basic.obo")

        # v1 逻辑：按 [Term] 分割，提取 id, name, namespace, is_a
        # 输出格式: GO_ID\\tnamespace:name\\tfather1;father2;...
        disc_file = outdir_path / "GO2disc.gz"

        go_id_pattern = re.compile(r'^id:\s(GO:\d+)')
        name_pattern = re.compile(r'^name:\s(.*)')
        namespace_pattern = re.compile(r'^namespace:\s(.*)')
        is_a_pattern = re.compile(r'^is_a:\s(GO:\d+)')

        term_count = 0

        with open(obo_path, 'r', encoding='utf-8') as f_in, \
             gzip.open(disc_file, 'wt', encoding='utf-8') as f_out:

            content = f_in.read()
            # 按 [Term] 分割，跳过第一个空块
            # v1 使用 $/="[Term]" 来设置记录分隔符
            blocks = content.split('[Term]')

            for block in blocks[1:]:  # 跳过第一个空块
                # 移除 [Typedef] 及之后的内容
                block = re.sub(r'\[Typedef\].*$', '', block, flags=re.DOTALL)

                lines = block.strip().split('\n')
                go_id = None
                name = None
                namespace = None
                fathers = []

                for line in lines:
                    line = line.strip()
                    if not line:
                        continue

                    m = go_id_pattern.match(line)
                    if m:
                        go_id = m.group(1)
                        continue

                    m = name_pattern.match(line)
                    if m:
                        name = m.group(1)
                        continue

                    m = namespace_pattern.match(line)
                    if m:
                        namespace = m.group(1)
                        continue

                    m = is_a_pattern.match(line)
                    if m:
                        fathers.append(m.group(1))

                if go_id and name and namespace:
                    # 格式: GO_ID\\tnamespace:name\\tfather1;father2;...
                    father_str = ";".join(fathers) if fathers else ""
                    f_out.write(f"{go_id}\t{namespace}:{name}\t{father_str}\n")
                    term_count += 1

        print(f"|--- 共解析 {term_count} 个 GO Term")
        print(f"|--- 写入文件: {disc_file}")
        print(f"|--- GOParser: obo 解析完成")
