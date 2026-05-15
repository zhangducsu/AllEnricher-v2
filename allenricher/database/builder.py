"""
数据库构建器模块

对应 v1 脚本：make_speciesDB（总入口）+ makeDB.go.v1.0.sh / makeDB.reactome.v1.0.sh ...（各数据库构建）

从 database/basic/ 中的全体物种通用数据中，提取指定物种的数据，
格式化输出到 database/organism/v{date}/{species}/ 目录。

架构（与 v1 一致）：
  database/
    basic/                          ← 全体物种通用数据（download 阶段产出）
      go/GO{date}/
        gene2go.gz                  ← 全员基因-GO 映射
        gene_info.gz                ← 全员基因信息  
        go-basic.obo               ← GO 本体定义
      reactome/Reactome{date}/
        gene_info.gz
        NCBI2Reactome_All_Levels.txt.gz  ← 全员通路数据
    organism/                       ← 物种专属格式化数据（build 阶段产出）
      v{date}/{species}/
        {species}.GO2gene.tab.gz    ← 该物种的 GO 基因-条目矩阵
        GO2disc.gz                  ← GO 条目描述
        {species}.Reactome2gene.tab.gz
        {species}.Reactome2disc.gz
        ...

用法：
  builder = DatabaseBuilder(root_dir="./database")
  # 构建人类 GO 数据库
  org_dir = builder.build_go(species="hsa", taxid=9606)
  # 自动从 database/basic/go/ 最新版本读取，输出到 database/organism/v{date}/hsa/
"""

import os
from pathlib import Path
from datetime import datetime
from typing import Dict, List, Optional

from .parsers.go import GOParser
from .parsers.kegg import KEGGParser
from .parsers.reactome import ReactomeParser
from .parsers.do import DOParser
from .parsers.disgenet import DisGeNETParser
from .downloader import DataDownloader


class DatabaseBuilder:
    """物种专属数据库构建器

    从全体物种通用原始数据（database/basic/）中提取指定物种的数据，
    格式化后输出到 database/organism/v{date}/{species}/ 目录。

    对应 v1 的 make_speciesDB 脚本。

    Attributes:
        root_dir: 数据库根目录
        basic_dir: 基础数据目录（全体物种通用）
        organism_dir: 物种专属数据目录
    """

    def __init__(self, root_dir: str = "./database"):
        """初始化构建器

        Args:
            root_dir: 数据库根目录，必须已包含 basic/ 子目录
        """
        self.root_dir = Path(root_dir)
        self.basic_dir = self.root_dir / "basic"
        self.organism_dir = self.root_dir / "organism"

    # ============================
    # GO 数据库构建
    # ============================
    def build_go(self, species: str, taxid: int,
                 go_version: Optional[str] = None) -> str:
        """构建指定物种的 GO 数据库

        从 database/basic/go/GO{version}/ 中读取全体物种的：
        - gene2go.gz：过滤 taxid，生成物种专属矩阵
        - gene_info.gz：获取基因 ID→符号映射
        - go-basic.obo：生成 GO 条目描述

        输出到 database/organism/v{date}/{species}/ 目录。

        对应 v1 的 makeDB.go.v1.0.sh。

        Args:
            species: 物种缩写（如 hsa, mmu）
            taxid: NCBI 物种分类学 ID（如 9606）
            go_version: GO 基础数据版本号（如 "GO20250101"），
                        默认自动使用最新版本

        Returns:
            str: 输出目录路径

        Raises:
            FileNotFoundError: 当基础数据目录不存在时
        """
        # 确定 GO 版本
        if go_version is None:
            downloader = DataDownloader(root_dir=str(self.root_dir))
            go_version = downloader.get_latest_go_version()
            if go_version is None:
                raise FileNotFoundError(
                    "未找到 GO 基础数据。请先运行 download 下载数据：\n"
                    "  allenricher download go\n"
                    "  或\n"
                    "  python -m allenricher.cli download -d go"
                )

        # 输入目录
        go_basic = self.basic_dir / "go" / go_version
        if not go_basic.exists():
            raise FileNotFoundError(
                f"GO 基础数据目录不存在: {go_basic}\n"
                f"可用的 GO 版本: {DataDownloader(root_dir=str(self.root_dir)).list_go_versions()}"
            )

        # 输出目录
        date_str = datetime.now().strftime("%Y%m%d")
        outdir = self.organism_dir / f"v{date_str}" / species
        outdir.mkdir(parents=True, exist_ok=True)

        gene2go_path = go_basic / "gene2go.gz"
        gene_info_path = go_basic / "gene_info.gz"
        obo_path = go_basic / "go-basic.obo"

        for fpath in [gene2go_path, gene_info_path, obo_path]:
            if not fpath.exists():
                raise FileNotFoundError(f"缺失文件: {fpath}")

        print(f"\n{'='*60}")
        print(f"构建 GO 数据库: {species} (taxid={taxid})")
        print(f"数据源: {go_basic}")
        print(f"输出目录: {outdir}")
        print(f"{'='*60}")

        # Step 1: 解析 gene2go，生成 {species}.GO2gene.tab.gz 和 {species}.gene2go.txt
        print("|--- Step 1/2: 从 gene2go.gz 提取物种基因-GO 矩阵...")
        GOParser.parse_gene2go(
            gene2go_path=str(gene2go_path),
            gene_info_path=str(gene_info_path),
            taxid=taxid,
            species=species,
            outdir=str(outdir)
        )

        # Step 2: 解析 obo，生成 GO2disc.gz
        print("|--- Step 2/2: 从 go-basic.obo 生成 GO 条目描述...")
        GOParser.parse_obo(
            obo_path=str(obo_path),
            outdir=str(outdir)
        )

        # 验证输出文件
        expected_files = [
            f"{species}.GO2gene.tab.gz",
            f"{species}.gene2go.txt",
            "GO2disc.gz"
        ]
        for fname in expected_files:
            fpath = outdir / fname
            if fpath.exists():
                print(f"    ✅ {fname}")
            else:
                print(f"    ❌ {fname} - 未生成")

        print(f"\nGO 数据库构建完成 → {outdir}")
        return str(outdir)

    # ============================
    # Reactome 数据库构建
    # ============================
    def build_reactome(self, species: str, taxid: int,
                       reactome_version: Optional[str] = None) -> str:
        """构建指定物种的 Reactome 数据库

        从 database/basic/reactome/Reactome{version}/ 中读取全体物种的：
        - NCBI2Reactome_All_Levels.txt.gz：按物种代码过滤
        - gene_info.gz：获取基因 ID→符号映射

        输出到 database/organism/v{date}/{species}/ 目录。

        对应 v1 的 makeDB.reactome.v1.0.sh。

        支持 v1 的 16 种模式生物：
        bta, cel, cfa, dre, ddi, dme, gga, hsa, mmu, mtu,
        pfa, rno, sce, spo, ssc, xtr

        Args:
            species: 物种缩写（如 hsa, mmu）
            taxid: NCBI 物种分类学 ID
            reactome_version: Reactome 基础数据版本号，
                              默认自动使用最新版本

        Returns:
            str: 输出目录路径
        """
        if reactome_version is None:
            downloader = DataDownloader(root_dir=str(self.root_dir))
            reactome_version = downloader.get_latest_reactome_version()
            if reactome_version is None:
                raise FileNotFoundError(
                    "未找到 Reactome 基础数据。请先运行 download 下载数据：\n"
                    "  allenricher download reactome\n"
                    "  或\n"
                    "  python -m allenricher.cli download -d reactome"
                )

        re_basic = self.basic_dir / "reactome" / reactome_version
        if not re_basic.exists():
            raise FileNotFoundError(
                f"Reactome 基础数据目录不存在: {re_basic}"
            )

        date_str = datetime.now().strftime("%Y%m%d")
        outdir = self.organism_dir / f"v{date_str}" / species
        outdir.mkdir(parents=True, exist_ok=True)

        ncbi2reactome_path = re_basic / "NCBI2Reactome_All_Levels.txt.gz"
        gene_info_path = re_basic / "gene_info.gz"

        for fpath in [ncbi2reactome_path, gene_info_path]:
            if not fpath.exists():
                raise FileNotFoundError(f"缺失文件: {fpath}")

        print(f"\n{'='*60}")
        print(f"构建 Reactome 数据库: {species} (taxid={taxid})")
        print(f"数据源: {re_basic}")
        print(f"输出目录: {outdir}")
        print(f"{'='*60}")

        print("|--- 从 NCBI2Reactome 提取物种通路-基因矩阵...")
        ReactomeParser.parse_ncbi2reactome(
            ncbi2reactome_path=str(ncbi2reactome_path),
            gene_info_path=str(gene_info_path),
            taxid=taxid,
            species=species,
            outdir=str(outdir)
        )

        # 验证输出文件
        expected_files = [
            f"{species}.Reactome2gene.tab.gz",
            f"{species}.Reactome2disc.gz",
            f"{species}.gene2pathway.txt"
        ]
        for fname in expected_files:
            fpath = outdir / fname
            if fpath.exists():
                print(f"    ✅ {fname}")
            else:
                print(f"    ❌ {fname} - 未生成")

        print(f"\nReactome 数据库构建完成 → {outdir}")
        return str(outdir)

    # ============================
    # KEGG 数据库构建
    # ============================
    def build_kegg(self, species: str, taxid: int,
                   go_version: Optional[str] = None,
                   gene2pathway_path: Optional[str] = None,
                   pathway_summary_path: Optional[str] = None) -> str:
        """构建指定物种的 KEGG 数据库

        KEGG 需要先从 REST API 获取 gene2pathway.txt（基因-通路映射），
        这无法通过通用下载步骤完成，因为每个物种需要单独抓取。

        对应 v1 的 makeDB.kegg.v1.1.sh。

        Args:
            species: 物种缩写（如 hsa）
            taxid: NCBI 物种分类学 ID
            go_version: GO 基础数据版本号（用于获取 gene_info.gz）
            gene2pathway_path: 从 KEGG API 获取的基因-通路映射文件路径（外部）
            pathway_summary_path: 通路分类总结文件路径（可选）

        Returns:
            str: 输出目录路径
        """
        date_str = datetime.now().strftime("%Y%m%d")
        outdir = self.organism_dir / f"v{date_str}" / species
        outdir.mkdir(parents=True, exist_ok=True)

        # 从 GO basic 中获取 gene_info（如果没有专门提供的话）
        if go_version is None:
            downloader = DataDownloader(root_dir=str(self.root_dir))
            go_version = downloader.get_latest_go_version()
            if go_version is None:
                go_version = f"GO{date_str}"

        gene_info_path = self.basic_dir / "go" / go_version / "gene_info.gz"
        if not gene_info_path.exists():
            # 尝试从 reactome basic 获取
            re_ver = DataDownloader(root_dir=str(self.root_dir)).get_latest_reactome_version()
            if re_ver:
                alt_path = self.basic_dir / "reactome" / re_ver / "gene_info.gz"
                if alt_path.exists():
                    gene_info_path = alt_path

        print(f"\n{'='*60}")
        print(f"构建 KEGG 数据库: {species} (taxid={taxid})")
        print(f"输出目录: {outdir}")
        print(f"{'='*60}")

        if gene2pathway_path is None:
            # 自动通过 KEGG REST API 获取数据
            from .kegg_fetcher import KEGGFetcher
            fetcher = KEGGFetcher(
                cache_dir=str(self.basic_dir / "kegg"),
                overwrite=False,
            )
            try:
                gene2pathway_path, pathway_summary_path = fetcher.fetch_species_data(
                    species=species,
                    gene_info_path=str(gene_info_path),
                )
            except Exception as e:
                print(f"|--- [错误] KEGG REST API 获取失败: {e}")
                print("    请检查网络连接后重试")
                return str(outdir)

        KEGGParser.build_database(
            species=species,
            gene_info_path=str(gene_info_path),
            gene2pathway_path=gene2pathway_path,
            outdir=str(outdir),
            pathway_summary_path=pathway_summary_path
        )

        print(f"\nKEGG 数据库构建完成 → {outdir}")
        return str(outdir)

    # ============================
    # DO / DisGeNET（仅人类）
    #
    # 对应 v1 make_speciesDB L127-135:
    #   if [ $organism == "hsa" ]; then
    #     sh $bin/src/makeDB.do.v1.0.sh ...
    #     sh $bin/src/makeDB.DisGeNET.v1.0.sh ...
    #   fi
    #
    # DO/DisGeNET 数据源本身只含人类基因关联，无需按物种过滤。
    # 仅在 build_species_db() 中检查 species=="hsa" 后构建。
    # ============================
    def build_do(self, taxid: int,
                 go_version: Optional[str] = None) -> str:
        """构建人类 Disease Ontology 数据库（仅 hsa）

        从 database/basic/do/ 中读取 Jensen Lab 的 human_disease_*.tsv 文件，
        用 gene_info 过滤有效人类基因，生成 hsa.DO2gene.tab.gz 和 hsa.DO2disc.gz。

        对应 v1 的 makeDB.do.v1.0.sh。

        Args:
            taxid: 人类 taxid = 9606
            go_version: GO 基础数据版本号（用于获取 gene_info.gz）

        Returns:
            str: 输出目录路径
        """
        species = "hsa"
        date_str = datetime.now().strftime("%Y%m%d")
        outdir = self.organism_dir / f"v{date_str}" / species
        outdir.mkdir(parents=True, exist_ok=True)

        # 查找 disease 文件（支持 .tsv 和 .tsv.gz）
        do_dir = self.basic_dir / "do"
        disease_files = sorted(do_dir.glob("human_disease_*_filtered.tsv.gz"))
        if not disease_files:
            disease_files = sorted(do_dir.glob("human_disease_*_filtered.tsv"))
        if not disease_files:
            disease_files = sorted(do_dir.glob("human_disease_*.tsv.gz"))
        if not disease_files:
            disease_files = sorted(do_dir.glob("human_disease_*.tsv"))
        if not disease_files:
            raise FileNotFoundError(
                f"在 {do_dir} 中未找到 DO 数据文件。请先运行 download:\n"
                f"  allenricher download do"
            )

        # 获取 gene_info
        if go_version is None:
            downloader = DataDownloader(root_dir=str(self.root_dir))
            go_version = downloader.get_latest_go_version()
        gene_info_path = self.basic_dir / "go" / go_version / "gene_info.gz"

        print(f"\n{'='*60}")
        print(f"构建 DO 数据库 (taxid={taxid})")
        print(f"数据源: {do_dir}")
        print(f"输出目录: {outdir}")
        print(f"{'='*60}")

        DOParser.parse_disease_files(
            disease_files=[str(f) for f in disease_files],
            gene_info_path=str(gene_info_path),
            taxid=taxid,
            outdir=str(outdir)
        )

        print(f"\nDO 数据库构建完成 → {outdir}")
        return str(outdir)

    def build_disgenet(self, taxid: int,
                       go_version: Optional[str] = None) -> str:
        """构建人类 DisGeNET 数据库

        从 database/basic/disgenet/ 中读取 all_gene_disease_associations.tsv.gz，
        生成 hsa.CUI2gene.tab.gz 和 hsa.CUI2disc.gz。

        对应 v1 的 makeDB.DisGeNET.v1.0.sh。

        Args:
            taxid: 人类 taxid = 9606
            go_version: GO 基础数据版本号（用于获取 gene_info.gz）

        Returns:
            str: 输出目录路径
        """
        species = "hsa"
        date_str = datetime.now().strftime("%Y%m%d")
        outdir = self.organism_dir / f"v{date_str}" / species
        outdir.mkdir(parents=True, exist_ok=True)

        assoc_path = self.basic_dir / "disgenet" / "all_gene_disease_associations.tsv.gz"
        if not assoc_path.exists():
            raise FileNotFoundError(
                f"未找到 DisGeNET 数据: {assoc_path}。请先运行 download:\n"
                f"  allenricher download disgenet"
            )

        if go_version is None:
            downloader = DataDownloader(root_dir=str(self.root_dir))
            go_version = downloader.get_latest_go_version()
        gene_info_path = self.basic_dir / "go" / go_version / "gene_info.gz"

        print(f"\n{'='*60}")
        print(f"构建 DisGeNET 数据库 (taxid={taxid})")
        print(f"数据源: {assoc_path}")
        print(f"输出目录: {outdir}")
        print(f"{'='*60}")

        DisGeNETParser.parse_associations(
            assoc_path=str(assoc_path),
            gene_info_path=str(gene_info_path),
            taxid=taxid,
            outdir=str(outdir)
        )

        print(f"\nDisGeNET 数据库构建完成 → {outdir}")
        return str(outdir)

    # ============================
    # 一键构建（对应 make_speciesDB）
    # ============================
    def build_species_db(self, species: str, taxid: int,
                         databases: List[str] = None,
                         go_version: Optional[str] = None,
                         reactome_version: Optional[str] = None) -> str:
        """一键构建指定物种的所有数据库

        对应 v1 的 make_speciesDB 脚本。
        从 database/basic/ 中自动获取最新版本的通用数据，
        为指定物种构建 GO、KEGG、Reactome 等数据库。

        Args:
            species: 物种缩写（如 hsa）
            taxid: NCBI 物种分类学 ID（如 9606）
            databases: 要构建的数据库列表，默认 ['GO', 'Reactome']
            go_version: GO 基础数据版本号，默认自动最新
            reactome_version: Reactome 基础数据版本号，默认自动最新

        Returns:
            str: 输出目录路径
        """
        if databases is None:
            databases = ['GO', 'Reactome']

        date_str = datetime.now().strftime("%Y%m%d")
        outdir = self.organism_dir / f"v{date_str}" / species

        print(f"\n{'#'*60}")
        print(f"# AllEnricher v2 物种数据库构建 (make_speciesDB)")
        print(f"# 物种: {species} (taxid={taxid})")
        print(f"# 数据库: {', '.join(databases)}")
        print(f"# 输出目录: {outdir}")
        print(f"{'#'*60}")

        for db_name in databases:
            db_upper = db_name.upper().strip()
            try:
                if db_upper == 'GO':
                    self.build_go(species, taxid, go_version)
                elif db_upper == 'REACTOME':
                    self.build_reactome(species, taxid, reactome_version)
                elif db_upper == 'KEGG':
                    self.build_kegg(species, taxid, go_version)
                elif db_upper == 'DO':
                    if species.lower() == 'hsa':
                        self.build_do(taxid, go_version)
                    else:
                        print(f"|--- [跳过] DO 仅支持人类 (hsa)")
                elif db_upper == 'DISGENET':
                    if species.lower() == 'hsa':
                        self.build_disgenet(taxid, go_version)
                    else:
                        print(f"|--- [跳过] DisGeNET 仅支持人类 (hsa)")
                else:
                    print(f"|--- [警告] 未知数据库: {db_name}")
            except Exception as e:
                print(f"|--- [错误] 构建 {db_name} 失败: {e}")

        print(f"\n{'#'*60}")
        print(f"# 物种数据库构建完成 → {outdir}")
        print(f"{'#'*60}")
        return str(outdir)
