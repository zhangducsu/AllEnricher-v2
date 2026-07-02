"""
GMT 基因集文件生成器

从物种专属数据库产物（GO2gene.tab.gz / GO2disc.gz 等）中提取基因集信息，
生成 GMT 格式的基因集文件，供 GSEA / ssGSEA / GSVA 使用。

GMT 格式：
    pathway_name<TAB>description<TAB>gene1<TAB>gene2<TAB>gene3...

输出为 .gmt.gz 压缩文件。
"""

import gzip
from pathlib import Path
from typing import Dict, List, Optional, Tuple


class GMTGenerator:
    """GMT 基因集文件生成器

    从物种专属数据库目录中读取各数据库的基因-条目矩阵和描述文件，
    转换为标准 GMT 格式输出。

    Attributes:
        organism_dir: 物种数据库目录路径
    """

    def __init__(self, organism_dir: str):
        """初始化 GMT 生成器

        Args:
            organism_dir: 物种数据库目录路径，
                          通常为 database/organism/v{date}/{species}/
        """
        self.organism_dir = Path(organism_dir)

    # ============================
    # 内部工具方法
    # ============================
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

    def _read_tab_matrix(self, tab_path: str) -> Tuple[List[str], Dict[str, List[str]]]:
        """读取基因-条目 0/1 矩阵文件，提取每个条目关联的基因列表

        矩阵格式：
            表头: Gene<TAB>term1<TAB>term2<TAB>...
            数据: gene_symbol<TAB>0/1<TAB>0/1<TAB>...

        Args:
            tab_path: 矩阵文件路径（.tab.gz 或 .tab）

        Returns:
            (terms, term_to_genes): terms 为条目 ID 列表，
                term_to_genes 为 {term_id: [gene1, gene2, ...]}
        """
        terms = []
        term_to_genes: Dict[str, List[str]] = {}

        with self._open_gz_or_text(tab_path) as f:
            header_line = f.readline().strip()
            if not header_line:
                return terms, term_to_genes

            parts = header_line.split('\t')
            terms = parts[1:]  # 跳过第一列 "Gene"

            # 初始化
            for t in terms:
                term_to_genes[t] = []

            for line in f:
                line = line.strip()
                if not line:
                    continue
                cols = line.split('\t')
                gene = cols[0]
                for i, val in enumerate(cols[1:]):
                    if i < len(terms) and val == '1':
                        term_to_genes[terms[i]].append(gene)

        return terms, term_to_genes

    def _read_description(self, disc_path: str) -> Dict[str, str]:
        """读取描述文件，建立 term_id -> description 映射

        描述文件格式：
            term_id<TAB>description

        Args:
            disc_path: 描述文件路径（.gz 或文本）

        Returns:
            {term_id: description}
        """
        descriptions: Dict[str, str] = {}

        with self._open_gz_or_text(disc_path) as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                parts = line.split('\t', 1)
                if len(parts) >= 2:
                    descriptions[parts[0]] = parts[1]
                elif len(parts) == 1 and parts[0]:
                    descriptions[parts[0]] = ""

        return descriptions

    def _write_gmt(self, term_to_genes: Dict[str, List[str]],
                   descriptions: Dict[str, str],
                   output_path: str) -> str:
        """将基因集数据写入 GMT 格式压缩文件

        Args:
            term_to_genes: {term_id: [gene1, gene2, ...]}
            descriptions: {term_id: description}
            output_path: 输出文件路径（.gmt.gz）

        Returns:
            输出文件路径
        """
        out_path = Path(output_path)
        out_path.parent.mkdir(parents=True, exist_ok=True)

        count = 0
        with gzip.open(out_path, 'wt', encoding='utf-8') as f:
            for term_id, genes in sorted(term_to_genes.items()):
                if not genes:
                    continue
                desc = descriptions.get(term_id, "")
                line = f"{term_id}\t{desc}\t" + "\t".join(genes) + "\n"
                f.write(line)
                count += 1

        print(f"|--- GMTGenerator: 写入 {count} 个基因集 → {out_path}")
        return str(out_path)

    # ============================
    # 各数据库 GMT 生成
    # ============================
    def generate_go_gmt(self, species: str) -> str:
        """从 GO 数据库产物生成 GMT 文件

        读取 {species}.GO2gene.tab.gz 和 GO2disc.gz，
        生成 {species}.GO.gmt.gz。

        Args:
            species: 物种缩写（如 hsa）

        Returns:
            输出文件路径

        Raises:
            FileNotFoundError: 当所需数据库文件不存在时
        """
        tab_path = self.organism_dir / f"{species}.GO2gene.tab.gz"
        disc_path = self.organism_dir / "GO2disc.gz"

        if not tab_path.exists():
            raise FileNotFoundError(f"GO 基因矩阵文件不存在: {tab_path}")
        if not disc_path.exists():
            raise FileNotFoundError(f"GO 描述文件不存在: {disc_path}")

        print(f"|--- GMTGenerator: 生成 GO GMT (species={species})")

        terms, term_to_genes = self._read_tab_matrix(str(tab_path))
        descriptions = self._read_description(str(disc_path))

        output_path = str(self.organism_dir / f"{species}.GO.gmt.gz")
        return self._write_gmt(term_to_genes, descriptions, output_path)

    def generate_kegg_gmt(self, species: str) -> str:
        """从 KEGG 数据库产物生成 GMT 文件

        读取 {species}.kegg2gene.tab.gz 和 {species}.kegg2disc.gz，
        生成 {species}.KEGG.gmt.gz。

        Args:
            species: 物种缩写（如 hsa）

        Returns:
            输出文件路径

        Raises:
            FileNotFoundError: 当所需数据库文件不存在时
        """
        tab_path = self.organism_dir / f"{species}.kegg2gene.tab.gz"
        disc_path = self.organism_dir / f"{species}.kegg2disc.gz"

        if not tab_path.exists():
            raise FileNotFoundError(f"KEGG 基因矩阵文件不存在: {tab_path}")
        if not disc_path.exists():
            raise FileNotFoundError(f"KEGG 描述文件不存在: {disc_path}")

        print(f"|--- GMTGenerator: 生成 KEGG GMT (species={species})")

        terms, term_to_genes = self._read_tab_matrix(str(tab_path))
        descriptions = self._read_description(str(disc_path))

        output_path = str(self.organism_dir / f"{species}.KEGG.gmt.gz")
        return self._write_gmt(term_to_genes, descriptions, output_path)

    def generate_reactome_gmt(self, species: str) -> str:
        """从 Reactome 数据库产物生成 GMT 文件

        读取 {species}.Reactome2gene.tab.gz 和 {species}.Reactome2disc.gz，
        生成 {species}.Reactome.gmt.gz。

        Args:
            species: 物种缩写（如 hsa）

        Returns:
            输出文件路径

        Raises:
            FileNotFoundError: 当所需数据库文件不存在时
        """
        tab_path = self.organism_dir / f"{species}.Reactome2gene.tab.gz"
        disc_path = self.organism_dir / f"{species}.Reactome2disc.gz"

        if not tab_path.exists():
            raise FileNotFoundError(f"Reactome 基因矩阵文件不存在: {tab_path}")
        if not disc_path.exists():
            raise FileNotFoundError(f"Reactome 描述文件不存在: {disc_path}")

        print(f"|--- GMTGenerator: 生成 Reactome GMT (species={species})")

        terms, term_to_genes = self._read_tab_matrix(str(tab_path))
        descriptions = self._read_description(str(disc_path))

        output_path = str(self.organism_dir / f"{species}.Reactome.gmt.gz")
        return self._write_gmt(term_to_genes, descriptions, output_path)

    def generate_do_gmt(self, species: str = "hsa") -> str:
        """从 DO 数据库产物生成 GMT 文件

        读取 hsa.DO2gene.tab.gz 和 hsa.DO2disc.gz，
        生成 hsa.DO.gmt.gz。

        Args:
            species: 物种缩写（仅支持 hsa）

        Returns:
            输出文件路径

        Raises:
            FileNotFoundError: 当所需数据库文件不存在时
        """
        tab_path = self.organism_dir / "hsa.DO2gene.tab.gz"
        disc_path = self.organism_dir / "hsa.DO2disc.gz"

        if not tab_path.exists():
            raise FileNotFoundError(f"DO 基因矩阵文件不存在: {tab_path}")
        if not disc_path.exists():
            raise FileNotFoundError(f"DO 描述文件不存在: {disc_path}")

        print(f"|--- GMTGenerator: 生成 DO GMT (species={species})")

        terms, term_to_genes = self._read_tab_matrix(str(tab_path))
        descriptions = self._read_description(str(disc_path))

        output_path = str(self.organism_dir / f"{species}.DO.gmt.gz")
        return self._write_gmt(term_to_genes, descriptions, output_path)

    def generate_disgenet_gmt(self, species: str = "hsa") -> str:
        """从 DisGeNET 数据库产物生成 GMT 文件

        读取 hsa.CUI2gene.tab.gz 和 hsa.CUI2disc.gz，
        生成 hsa.DisGeNET.gmt.gz。

        Args:
            species: 物种缩写（仅支持 hsa）

        Returns:
            输出文件路径

        Raises:
            FileNotFoundError: 当所需数据库文件不存在时
        """
        tab_path = self.organism_dir / "hsa.CUI2gene.tab.gz"
        disc_path = self.organism_dir / "hsa.CUI2disc.gz"

        if not tab_path.exists():
            raise FileNotFoundError(f"DisGeNET 基因矩阵文件不存在: {tab_path}")
        if not disc_path.exists():
            raise FileNotFoundError(f"DisGeNET 描述文件不存在: {disc_path}")

        print(f"|--- GMTGenerator: 生成 DisGeNET GMT (species={species})")

        terms, term_to_genes = self._read_tab_matrix(str(tab_path))
        descriptions = self._read_description(str(disc_path))

        output_path = str(self.organism_dir / f"{species}.DisGeNET.gmt.gz")
        return self._write_gmt(term_to_genes, descriptions, output_path)

    def generate_wikipathways_gmt(self, species: str) -> str:
        """从 WikiPathways 数据库产物生成 GMT 文件

        读取 {species}.WikiPathways2gene.tab.gz 和 {species}.WikiPathways2disc.gz，
        生成 {species}.WikiPathways.gmt.gz。

        Args:
            species: 物种缩写（如 hsa）

        Returns:
            输出文件路径

        Raises:
            FileNotFoundError: 当所需数据库文件不存在时
        """
        tab_path = self.organism_dir / f"{species}.WikiPathways2gene.tab.gz"
        disc_path = self.organism_dir / f"{species}.WikiPathways2disc.gz"

        if not tab_path.exists():
            raise FileNotFoundError(f"WikiPathways 基因矩阵文件不存在: {tab_path}")
        if not disc_path.exists():
            raise FileNotFoundError(f"WikiPathways 描述文件不存在: {disc_path}")

        print(f"|--- GMTGenerator: 生成 WikiPathways GMT (species={species})")

        terms, term_to_genes = self._read_tab_matrix(str(tab_path))
        descriptions = self._read_description(str(disc_path))

        output_path = str(self.organism_dir / f"{species}.WikiPathways.gmt.gz")
        return self._write_gmt(term_to_genes, descriptions, output_path)

    def generate_all_gmt(self, species: str) -> Dict[str, str]:
        """生成所有可用的 GMT 文件

        依次尝试生成 GO、KEGG、Reactome、DO、DisGeNET 的 GMT 文件。
        如果某个数据库的产物文件不存在，则跳过该数据库。

        Args:
            species: 物种缩写（如 hsa）

        Returns:
            {数据库名称: 输出文件路径}，仅包含成功生成的条目
        """
        results: Dict[str, str] = {}

        generators = [
            ("GO", lambda: self.generate_go_gmt(species)),
            ("KEGG", lambda: self.generate_kegg_gmt(species)),
            ("Reactome", lambda: self.generate_reactome_gmt(species)),
            ("DO", lambda: self.generate_do_gmt(species)),
            ("DisGeNET", lambda: self.generate_disgenet_gmt(species)),
        ]

        for db_name, gen_func in generators:
            try:
                output_path = gen_func()
                results[db_name] = output_path
            except FileNotFoundError:
                print(f"|--- GMTGenerator: 跳过 {db_name}（数据文件不存在）")

        print(f"|--- GMTGenerator: 共生成 {len(results)} 个 GMT 文件")
        return results
