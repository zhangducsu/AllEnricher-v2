"""
自定义数据库构建器

从用户提供的条目注释文件构建 AllEnricher 自定义数据库，
自动生成基因矩阵（0/1）、描述文件和 GMT 基因集文件。

对应 v1 的自定义注释文件处理流程。

用法：
    builder = CustomDatabaseBuilder(root_dir="./database")
    outdir = builder.build_from_annotation(
        annotation_file="my_annotation.txt",
        species="hsa",
        taxid=9606,
        db_name="CustomDB"
    )
"""

import gzip
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional

import pandas as pd

from .parsers.annotation_parser import AnnotationParser


class CustomDatabaseBuilder:
    """自定义数据库构建器

    从用户提供的注释文件构建 AllEnricher 标准格式的数据库文件：
    - {species}.{db_name}2gene.tab.gz: 基因-条目 0/1 矩阵
    - {db_name}2disc.gz: 条目描述（含层级信息）
    - {species}.{db_name}.gmt.gz: GMT 基因集文件（自动生成）

    Attributes:
        root_dir: 数据库根目录
    """

    def __init__(self, root_dir: str = "./database"):
        """初始化构建器

        Args:
            root_dir: 数据库根目录
        """
        self.root_dir = Path(root_dir)

    def build_from_annotation(
        self,
        annotation_file: str,
        species: str,
        taxid: int,
        db_name: str,
        format_type: Optional[str] = None,
        hierarchy_separator: str = '|'
    ) -> str:
        """从条目注释文件构建数据库，自动生成GMT文件

        构建流程：
        1. 使用 AnnotationParser 解析注释文件
        2. 获取 term_genes 映射
        3. 创建输出目录 database/organism/v{YYYYMMDD}/{species}/
        4. 生成 {species}.{db_name}2gene.tab.gz（基因-条目 0/1 矩阵）
        5. 生成 {db_name}2disc.gz（条目描述，含层级信息）
        6. 自动生成 {species}.{db_name}.gmt.gz（GMT 基因集文件）

        Args:
            annotation_file: 注释文件路径
            species: 物种缩写（如 hsa, mmu）
            taxid: NCBI 物种分类学 ID
            db_name: 数据库名称（如 CustomDB, MyPathway）
            format_type: 文件格式类型（'2col', '3col', '4col'），
                         默认自动检测
            hierarchy_separator: 层级分隔符，默认 '|'

        Returns:
            str: 输出目录路径

        Raises:
            FileNotFoundError: 当注释文件不存在时
            ValueError: 当注释文件为空时
        """
        annotation_path = Path(annotation_file)
        if not annotation_path.exists():
            raise FileNotFoundError(f"注释文件不存在: {annotation_file}")

        # Step 1: 解析注释文件
        try:
            parser = AnnotationParser(
                filepath=str(annotation_path),
                format_type=format_type,
                hierarchy_separator=hierarchy_separator
            )
            parser.parse()
            term_genes = parser.get_term_genes()
            term_names = parser.get_term_names()
            term_hierarchies = parser.get_term_hierarchies()
        except (ValueError, FileNotFoundError):
            raise ValueError(
                f"注释文件中没有有效的基因-条目映射: {annotation_file}"
            )

        if not term_genes:
            raise ValueError(f"注释文件中没有有效的基因-条目映射: {annotation_file}")

        # Step 2: 创建输出目录
        date_str = datetime.now().strftime("%Y%m%d")
        outdir = self.root_dir / "organism" / f"v{date_str}" / species
        outdir.mkdir(parents=True, exist_ok=True)

        print(f"\n{'='*60}")
        print(f"构建自定义数据库: {species}.{db_name}")
        print(f"注释文件: {annotation_file}")
        print(f"输出目录: {outdir}")
        print(f"条目数: {len(term_genes)}")
        print(f"{'='*60}")

        # Step 3: 生成基因矩阵
        print("|--- Step 1/3: 生成基因-条目矩阵...")
        self._create_gene_matrix(term_genes, species, db_name, outdir)

        # Step 4: 生成描述文件
        print("|--- Step 2/3: 生成条目描述文件...")
        self._create_description_file(
            term_names, term_hierarchies, db_name, outdir
        )

        # Step 5: 自动生成 GMT 文件
        print("|--- Step 3/3: 自动生成 GMT 文件...")
        self._create_gmt_file(term_genes, term_names, species, db_name, outdir)

        # 验证输出文件
        expected_files = [
            f"{species}.{db_name}2gene.tab.gz",
            f"{db_name}2disc.gz",
            f"{species}.{db_name}.gmt.gz",
        ]
        for fname in expected_files:
            fpath = outdir / fname
            if fpath.exists():
                print(f"    [OK] {fname}")
            else:
                print(f"    [FAIL] {fname} - 未生成")

        print(f"\n自定义数据库构建完成 -> {outdir}")
        return str(outdir)

    def _create_gene_matrix(
        self,
        term_genes: Dict[str, List[str]],
        species: str,
        db_name: str,
        outdir: Path
    ) -> str:
        """生成基因-条目 0/1 矩阵

        构建 DataFrame，列: Gene, term1, term2, ...
        每行表示一个基因，值为 0 或 1 表示该基因是否属于对应条目。

        Args:
            term_genes: {term_id: [gene1, gene2, ...]}
            species: 物种缩写
            db_name: 数据库名称
            outdir: 输出目录

        Returns:
            str: 输出文件路径
        """
        # 收集所有基因（保持有序）
        all_genes: List[str] = []
        seen_genes = set()
        for genes in term_genes.values():
            for gene in genes:
                if gene not in seen_genes:
                    all_genes.append(gene)
                    seen_genes.add(gene)

        # 收集所有条目（排序）
        terms = sorted(term_genes.keys())

        # 构建 0/1 矩阵
        data = {"Gene": all_genes}
        for term in terms:
            gene_set = set(term_genes[term])
            data[term] = [1 if g in gene_set else 0 for g in all_genes]

        df = pd.DataFrame(data)

        # 保存为 gzip 压缩 TSV
        output_path = outdir / f"{species}.{db_name}2gene.tab.gz"
        df.to_csv(output_path, sep='\t', index=False, compression='gzip')

        print(f"    基因矩阵: {len(all_genes)} 基因 x {len(terms)} 条目 -> {output_path.name}")
        return str(output_path)

    def _create_description_file(
        self,
        term_names: Dict[str, str],
        term_hierarchies: Dict[str, str],
        db_name: str,
        outdir: Path
    ) -> str:
        """生成条目描述文件

        格式: term_id<TAB>term_name<TAB>hierarchy
        如果没有层级信息，使用 term_name 作为层级。

        Args:
            term_names: {term_id: term_name}
            term_hierarchies: {term_id: hierarchy_string}
            db_name: 数据库名称
            outdir: 输出目录

        Returns:
            str: 输出文件路径
        """
        output_path = outdir / f"{db_name}2disc.gz"

        with gzip.open(output_path, 'wt', encoding='utf-8') as f:
            for term_id in sorted(term_names.keys()):
                term_name = term_names[term_id]
                hierarchy = term_hierarchies.get(term_id, term_name)
                f.write(f"{term_id}\t{term_name}\t{hierarchy}\n")

        print(f"    描述文件: {len(term_names)} 条目 -> {output_path.name}")
        return str(output_path)

    def _create_gmt_file(
        self,
        term_genes: Dict[str, List[str]],
        term_names: Dict[str, str],
        species: str,
        db_name: str,
        outdir: Path
    ) -> str:
        """自动生成 GMT 基因集文件

        格式: term_id<TAB>term_name<TAB>gene1<TAB>gene2...
        从 term_genes 映射自动生成。

        Args:
            term_genes: {term_id: [gene1, gene2, ...]}
            term_names: {term_id: term_name}
            species: 物种缩写
            db_name: 数据库名称
            outdir: 输出目录

        Returns:
            str: 输出文件路径
        """
        output_path = outdir / f"{species}.{db_name}.gmt.gz"

        count = 0
        with gzip.open(output_path, 'wt', encoding='utf-8') as f:
            for term_id in sorted(term_genes.keys()):
                genes = term_genes[term_id]
                if not genes:
                    continue
                term_name = term_names.get(term_id, "")
                line = f"{term_id}\t{term_name}\t" + "\t".join(genes) + "\n"
                f.write(line)
                count += 1

        print(f"    GMT 文件: {count} 个基因集 -> {output_path.name}")
        return str(output_path)
