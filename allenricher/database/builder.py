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
import gzip
import logging
from pathlib import Path
from datetime import datetime
from typing import Dict, List, Optional

from .parsers.go import GOParser
from .parsers.kegg import KEGGParser
from .parsers.reactome import ReactomeParser
from .parsers.do import DOParser
from .parsers.disgenet import DisGeNETParser
from .parsers.wikipathways import WikiPathwaysParser
from .parsers.trrust import TRRUSTParser
from .parsers.chea3 import ChEA3Parser
from .downloader import DataDownloader
from .gmt_generator import GMTGenerator
from .species_registry import SpeciesRegistry, SpeciesEntry
from .goa_fetcher import GOAFetcher
from .wikipathways_fetcher import WikiPathwaysFetcher
from .trrust_fetcher import TRRUSTFetcher

logger = logging.getLogger(__name__)


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
        # 自动生成 GMT 文件
        self.generate_gmt_files(species, str(outdir))
        return str(outdir)

    # ============================
    # GO 数据库构建（GOA 来源）
    # ============================
    def _get_species_dir(self, taxid: int, latin_name: str) -> str:
        """生成物种目录名: taxid.物种拉丁名（空格替换为下划线）"""
        return f"{taxid}.{latin_name.replace(' ', '_')}"

    def _get_species_prefix(self, taxid: int, latin_name: str) -> str:
        """生成文件名前缀: taxid.物种拉丁名（空格替换为下划线）"""
        return f"{taxid}.{latin_name.replace(' ', '_')}"

    def build_go_from_goa(self, taxid: int, latin_name: str,
                          goa_filename: str, go_version: str = None) -> str:
        """从 UniProt GOA proteomes 构建 GO 数据库

        步骤:
        1. 确定版本号（使用 go_version 或自动检测）
        2. 获取物种目录名: f"{taxid}.{latin_name.replace(' ', '_')}"
        3. 创建输出目录: database/organism/v{date}/{species_dir}/
        4. 初始化 GOAFetcher，下载 GOA 文件
        5. 解析 GOA 文件，获取 gene_to_go 和 all_genes
        6. 从 go-basic.obo 获取 GO term 名称
        7. 生成 GO2gene.tab.gz（0/1矩阵）
        8. 生成 gene2go.txt（注释列表）
        9. 生成 GO2disc.gz（GO描述）
        10. 生成 GMT 文件

        参数:
            taxid: NCBI Taxonomy ID
            latin_name: 物种拉丁名
            goa_filename: GOA 文件名（如 "9606.Homo_sapiens.goa"）
            go_version: GO 版本号

        返回: 输出目录路径
        """
        # Step 1: 确定 GO 版本号
        if go_version is None:
            downloader = DataDownloader(root_dir=str(self.root_dir))
            go_version = downloader.get_latest_go_version()
            if go_version is None:
                date_str = datetime.now().strftime("%Y%m%d")
                go_version = f"GO{date_str}"
                logger.warning("未找到 GO 基础数据版本，使用默认版本: %s", go_version)

        # Step 2: 物种目录名和文件前缀
        species_dir = self._get_species_dir(taxid, latin_name)
        prefix = self._get_species_prefix(taxid, latin_name)

        # Step 3: 创建输出目录
        date_str = datetime.now().strftime("%Y%m%d")
        outdir = self.organism_dir / f"v{date_str}" / species_dir
        outdir.mkdir(parents=True, exist_ok=True)

        # 获取 go-basic.obo 路径（用于提取 GO term 名称和生成 GO2disc.gz）
        go_basic = self.basic_dir / "go" / go_version
        obo_path = go_basic / "go-basic.obo"

        print(f"\n{'='*60}")
        print(f"构建 GO 数据库 (GOA 来源): {species_dir} (taxid={taxid})")
        print(f"GOA 文件: {goa_filename}")
        print(f"GO 版本: {go_version}")
        print(f"输出目录: {outdir}")
        print(f"{'='*60}")

        # Step 4: 初始化 GOAFetcher，下载 GOA 文件
        goa_date = go_version.replace("GO", "") if go_version.startswith("GO") else date_str
        goa_cache_dir = self.basic_dir / "goa" / f"GOA{goa_date}"
        fetcher = GOAFetcher(cache_dir=str(goa_cache_dir), overwrite=False)

        latin_name_underscore = latin_name.replace(" ", "_")
        goa_file = fetcher.fetch_species_data(
            taxid=taxid,
            latin_name=latin_name_underscore,
            goa_filename=goa_filename,
        )
        logger.info("GOA 文件已就绪: %s", goa_file)

        # Step 5: 解析 GOA 文件
        print("|--- Step 1/4: 解析 GOA 文件...")
        gene_to_go, all_genes = fetcher.parse_goa_file(goa_file, taxid)

        if not all_genes:
            raise ValueError(
                f"GOA 文件中未找到 taxid={taxid} 的有效基因注释"
            )

        # 收集所有 GO term
        all_go_terms: set = set()
        for go_set in gene_to_go.values():
            all_go_terms.update(go_set)

        logger.info("共找到 %d 个基因, %d 个 GO term", len(all_genes), len(all_go_terms))

        # Step 6: 从 go-basic.obo 获取 GO term 名称
        go_names: Dict[str, str] = {}
        if obo_path.exists():
            print("|--- Step 2/4: 从 go-basic.obo 提取 GO term 名称...")
            go_names = self._extract_go_names_from_obo(str(obo_path))
            logger.info("从 obo 提取到 %d 个 GO term 名称", len(go_names))
        else:
            logger.warning("go-basic.obo 不存在: %s，GO term 名称将为空", obo_path)

        # Step 7: 生成 GO2gene.tab.gz
        print("|--- Step 3/4: 生成 GO2gene.tab.gz...")
        go2gene_path = outdir / f"{prefix}.GO2gene.tab.gz"
        GOAFetcher.build_go2gene_matrix(
            gene_to_go=gene_to_go,
            all_genes=all_genes,
            all_go_terms=all_go_terms,
            output_path=go2gene_path,
        )

        # Step 8: 生成 gene2go.txt
        print("|--- Step 3/4: 生成 gene2go.txt...")
        gene2go_path = outdir / f"{prefix}.gene2go.txt"
        GOAFetcher.build_gene2go_list(
            gene_to_go=gene_to_go,
            go_names=go_names,
            output_path=gene2go_path,
        )

        # Step 9: 生成 GO2disc.gz
        print("|--- Step 4/4: 生成 GO2disc.gz...")
        if obo_path.exists():
            GOParser.parse_obo(obo_path=str(obo_path), outdir=str(outdir))
        else:
            logger.warning("跳过 GO2disc.gz 生成（go-basic.obo 不存在）")

        # 验证输出文件
        expected_files = [
            f"{prefix}.GO2gene.tab.gz",
            f"{prefix}.gene2go.txt",
            "GO2disc.gz",
        ]
        for fname in expected_files:
            fpath = outdir / fname
            if fpath.exists():
                print(f"    [OK] {fname}")
            else:
                print(f"    [MISSING] {fname} - 未生成")

        print(f"\nGO 数据库构建完成 (GOA 来源) -> {outdir}")

        # Step 10: 生成 GMT 文件
        self._generate_gmt_for_species_dir(species_dir, str(outdir), prefix)

        return str(outdir)

    def _extract_go_names_from_obo(self, obo_path: str) -> Dict[str, str]:
        """从 go-basic.obo 提取 GO ID -> name 映射

        Args:
            obo_path: go-basic.obo 文件路径

        Returns:
            {GO_ID: GO_name}
        """
        import re

        go_names: Dict[str, str] = {}
        go_id_pattern = re.compile(r'^id:\s(GO:\d+)')
        name_pattern = re.compile(r'^name:\s(.*)')

        with open(obo_path, 'r', encoding='utf-8') as f:
            current_id = None
            for line in f:
                line = line.strip()
                m = go_id_pattern.match(line)
                if m:
                    current_id = m.group(1)
                    continue
                m = name_pattern.match(line)
                if m and current_id:
                    go_names[current_id] = m.group(1)
                    current_id = None

        return go_names

    def _generate_gmt_for_species_dir(self, species_dir: str, output_dir: str,
                                       prefix: str) -> Dict[str, str]:
        """为 species_dir 格式的数据库目录生成 GMT 文件

        与 generate_gmt_files 类似，但使用 prefix（如 9606.Homo_sapiens）
        而非 species 缩写（如 hsa）来定位矩阵文件。

        Args:
            species_dir: 物种目录名（如 "9606.Homo_sapiens"）
            output_dir: 数据库输出目录路径
            prefix: 文件名前缀（如 "9606.Homo_sapiens"）

        Returns:
            {数据库名称: GMT文件路径}
        """
        generator = GMTGenerator(organism_dir=output_dir)
        results: Dict[str, str] = {}

        # 尝试生成 GO GMT
        tab_path = Path(output_dir) / f"{prefix}.GO2gene.tab.gz"
        disc_path = Path(output_dir) / "GO2disc.gz"
        if tab_path.exists() and disc_path.exists():
            try:
                terms, term_to_genes = generator._read_tab_matrix(str(tab_path))
                descriptions = generator._read_description(str(disc_path))
                gmt_path = str(Path(output_dir) / f"{prefix}.GO.gmt.gz")
                generator._write_gmt(term_to_genes, descriptions, gmt_path)
                results["GO"] = gmt_path
            except Exception as e:
                logger.warning("GO GMT 生成失败: %s", e)

        return results

    # ============================
    # GO 数据库构建（带回退）
    # ============================
    def build_go_with_fallback(self, taxid: int, latin_name: str,
                               go_version: str = None) -> str:
        """带回退的 GO 构建：先尝试 gene2go，失败则使用 GOA

        步骤:
        1. 检查 gene2go.gz 是否包含该 taxid
        2. 如果包含，使用现有 build_go 方法
        3. 如果不包含，查找 GOA 注册表
        4. 使用 GOA 构建

        返回: 输出目录路径
        """
        logger.info("build_go_with_fallback: taxid=%d, latin_name=%s", taxid, latin_name)

        # Step 1: 检查 gene2go 是否包含该 taxid
        use_gene2go = False
        species_abbr = None

        # 方法 A: 读取 supported_species.tsv 检查 go_source
        registry_path = self.root_dir / "supported_species.tsv"
        if registry_path.exists():
            registry = SpeciesRegistry(registry_path=registry_path)
            registry.load()
            entry = registry.query_by_taxid(taxid)
            if entry and entry.has_go and entry.go_source:
                if entry.go_source.lower() == "gene2go":
                    use_gene2go = True
                    logger.info("物种 %d 在注册表中标记为 gene2go 来源", taxid)

        # 方法 B: 如果注册表没有信息，尝试直接检查 gene2go.gz
        if not use_gene2go:
            use_gene2go = self._check_gene2go_has_taxid(taxid, go_version)

        # Step 2: 如果 gene2go 有该物种，使用现有 build_go 方法
        if use_gene2go:
            # 尝试获取 species 缩写
            species_abbr = self._get_species_abbr(taxid, latin_name)
            if species_abbr:
                logger.info("使用 gene2go 构建 GO 数据库, species=%s", species_abbr)
                return self.build_go(species=species_abbr, taxid=taxid,
                                     go_version=go_version)
            else:
                logger.warning("无法确定 species 缩写，回退到 GOA")

        # Step 3: 如果 gene2go 没有，查找 GOA 注册表
        goa_filename = self._find_goa_filename(taxid, latin_name)
        if goa_filename is None:
            raise ValueError(
                f"无法为 taxid={taxid} ({latin_name}) 找到 GO 数据源。\n"
                f"gene2go.gz 中未找到该物种，GOA 注册表中也未找到对应条目。"
            )

        # Step 4: 使用 GOA 构建
        logger.info("使用 GOA 构建 GO 数据库, goa_filename=%s", goa_filename)
        return self.build_go_from_goa(taxid, latin_name, goa_filename, go_version)

    def _check_gene2go_has_taxid(self, taxid: int,
                                  go_version: str = None) -> bool:
        """检查 gene2go.gz 是否包含指定 taxid 的数据

        全文件扫描以确保找到所有物种（包括人类等数据量大的物种）。

        Args:
            taxid: NCBI Taxonomy ID
            go_version: GO 版本号

        Returns:
            True 如果找到匹配的 taxid
        """
        if go_version is None:
            downloader = DataDownloader(root_dir=str(self.root_dir))
            go_version = downloader.get_latest_go_version()
            if go_version is None:
                return False

        gene2go_path = self.basic_dir / "go" / go_version / "gene2go.gz"
        if not gene2go_path.exists():
            return False

        taxid_str = str(taxid)

        try:
            with gzip.open(gene2go_path, 'rt', encoding='utf-8') as f:
                for i, line in enumerate(f, 1):
                    line = line.strip()
                    if not line or line.startswith('#'):
                        continue
                    parts = line.split('\t')
                    if parts[0] == taxid_str:
                        logger.info("在 gene2go.gz 中找到 taxid=%s (第 %d 行)",
                                    taxid_str, i)
                        return True
        except Exception as e:
            logger.warning("检查 gene2go.gz 失败: %s", e)
            return False

        logger.info("在 gene2go.gz 中未找到 taxid=%s", taxid_str)
        return False

    def _get_species_abbr(self, taxid: int, latin_name: str) -> Optional[str]:
        """获取物种缩写代码

        优先从 supported_species.tsv 读取 kegg_code，
        如果没有则根据拉丁名自动生成。

        Args:
            taxid: NCBI Taxonomy ID
            latin_name: 物种拉丁名

        Returns:
            物种缩写（如 hsa），无法确定则返回 None
        """
        # 优先从注册表获取
        registry_path = self.root_dir / "supported_species.tsv"
        if registry_path.exists():
            registry = SpeciesRegistry(registry_path=registry_path)
            registry.load()
            entry = registry.query_by_taxid(taxid)
            if entry and entry.kegg_code:
                return entry.kegg_code

        # 回退: 自动生成
        return SpeciesRegistry.generate_kegg_abbreviation(latin_name)

    def _find_goa_filename(self, taxid: int, latin_name: str) -> Optional[str]:
        """查找物种对应的 GOA 文件名

        优先从 supported_species.tsv 读取 go_filename，
        如果没有则根据 taxid 和 latin_name 自动构造。

        Args:
            taxid: NCBI Taxonomy ID
            latin_name: 物种拉丁名

        Returns:
            GOA 文件名（如 "9606.goa"），未找到则返回 None
        """
        # 优先从注册表获取
        registry_path = self.root_dir / "supported_species.tsv"
        if registry_path.exists():
            registry = SpeciesRegistry(registry_path=registry_path)
            registry.load()
            entry = registry.query_by_taxid(taxid)
            if entry and entry.has_go and entry.go_filename:
                return entry.go_filename

        # 回退: 自动构造 - EBI GOA 文件名格式为 {taxid}.goa（如 9606.goa）
        return f"{taxid}.goa"

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
        # 自动生成 GMT 文件
        self.generate_gmt_files(species, str(outdir))
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
        # 自动生成 GMT 文件
        self.generate_gmt_files(species, str(outdir))
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
        # 自动生成 GMT 文件
        self.generate_gmt_files("hsa", str(outdir))
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
        # 自动生成 GMT 文件
        self.generate_gmt_files("hsa", str(outdir))
        return str(outdir)

    # ============================
    # 辅助方法：获取 GO 版本和 gene_info 路径
    # ============================
    def _get_go_version(self) -> Optional[str]:
        """获取最新的 GO 版本号

        Returns:
            str: GO 版本号（如 "GO20250101"），未找到则返回 None
        """
        downloader = DataDownloader(root_dir=str(self.root_dir))
        return downloader.get_latest_go_version()

    def _get_gene_info_path(self) -> Optional[Path]:
        """从 GO 目录查找 gene_info.gz 文件路径

        Returns:
            Path: gene_info.gz 文件路径，未找到则返回 None
        """
        go_version = self._get_go_version()
        if go_version is None:
            return None

        gene_info_path = self.basic_dir / "go" / go_version / "gene_info.gz"
        if gene_info_path.exists():
            return gene_info_path

        return None

    def _load_valid_genes(self, gene_info_path: str, taxid: int) -> Optional[set]:
        """从 gene_info.gz 加载指定 taxid 的有效基因符号集合

        Args:
            gene_info_path: gene_info.gz 文件路径
            taxid: NCBI Taxonomy ID

        Returns:
            有效基因符号集合，失败则返回 None
        """
        taxid_str = str(taxid)
        valid_genes = set()
        try:
            with gzip.open(gene_info_path, 'rt', encoding='utf-8') as f:
                for line in f:
                    line = line.strip()
                    if not line or line.startswith('#'):
                        continue
                    parts = line.split('\t')
                    if len(parts) < 3:
                        continue
                    if parts[0] == taxid_str:
                        gene_symbol = parts[2]
                        if gene_symbol and gene_symbol != '-':
                            valid_genes.add(gene_symbol)
            logger.info("从 gene_info.gz 加载了 %d 个有效基因 (taxid=%s)",
                        len(valid_genes), taxid_str)
            return valid_genes
        except Exception as e:
            logger.warning("加载 gene_info.gz 失败: %s", e)
            return None

    # ============================
    # WikiPathways 数据库构建
    # ============================
    def build_wikipathways(self, species: str, taxid: int,
                           wikipathways_version: Optional[str] = None) -> str:
        """构建指定物种的 WikiPathways 数据库

        从 database/basic/wikipathways/WP{version}/ 中读取下载的 GMT 文件，
        使用 WikiPathwaysParser 生成标准格式的数据库文件。

        对应 v1 的 wikipathways.R 脚本。

        Args:
            species: 物种缩写（如 hsa, mmu）
            taxid: NCBI 物种分类学 ID（如 9606）
            wikipathways_version: WikiPathways 版本号（YYYYMMDD 格式），
                                  默认自动使用最新版本

        Returns:
            str: 输出目录路径

        Raises:
            FileNotFoundError: 当基础数据目录或 GMT 文件不存在时
            ValueError: 当无法找到物种对应的拉丁名时
        """
        # 获取拉丁名
        latin_name = WikiPathwaysFetcher.get_latin_name(species)
        if latin_name is None:
            raise ValueError(
                f"无法找到物种 {species} 对应的拉丁名，"
                f"WikiPathways 不支持该物种"
            )

        # 确定版本号
        if wikipathways_version is None:
            fetcher = WikiPathwaysFetcher(basic_dir=str(self.basic_dir))
            versions = fetcher.list_cached_versions()
            if not versions:
                raise FileNotFoundError(
                    "未找到 WikiPathways 基础数据。请先运行 download 下载数据：\n"
                    "  allenricher download wikipathways\n"
                    "  或\n"
                    "  python -m allenricher.cli download -d wikipathways"
                )
            wikipathways_version = versions[-1]  # 使用最新版本

        # 输入目录
        wp_basic = self.basic_dir / "wikipathways" / f"WP{wikipathways_version}"
        if not wp_basic.exists():
            raise FileNotFoundError(
                f"WikiPathways 基础数据目录不存在: {wp_basic}\n"
                f"可用的版本: {WikiPathwaysFetcher(basic_dir=str(self.basic_dir)).list_cached_versions()}"
            )

        # 查找 GMT 文件
        species_filename = latin_name.replace(" ", "_")
        gmt_filename = f"wikipathways-{wikipathways_version}-gmt-{species_filename}.gmt"
        gmt_path = wp_basic / gmt_filename

        if not gmt_path.exists():
            raise FileNotFoundError(
                f"WikiPathways GMT 文件不存在: {gmt_path}\n"
                f"请确保已下载 {species} ({latin_name}) 的数据"
            )

        # 输出目录
        date_str = datetime.now().strftime("%Y%m%d")
        outdir = self.organism_dir / f"v{date_str}" / species
        outdir.mkdir(parents=True, exist_ok=True)

        # 获取 gene_info 路径用于 ID 转换和过滤有效基因
        gene_info_path = self._get_gene_info_path()

        print(f"\n{'='*60}")
        print(f"构建 WikiPathways 数据库: {species} (taxid={taxid})")
        print(f"拉丁名: {latin_name}")
        print(f"数据源: {gmt_path}")
        print(f"输出目录: {outdir}")
        print(f"{'='*60}")

        # 使用 WikiPathwaysParser 构建数据库
        print("|--- 从 WikiPathways GMT 提取通路-基因矩阵...")
        WikiPathwaysParser.build_database(
            gmt_path=str(gmt_path),
            output_dir=str(outdir),
            species=species,
            taxid=taxid,
            gene_info_path=str(gene_info_path) if gene_info_path else None
        )

        # 验证输出文件
        expected_files = [
            f"{species}.WikiPathways2gene.tab.gz",
            f"{species}.WikiPathways2disc.gz",
        ]
        for fname in expected_files:
            fpath = outdir / fname
            if fpath.exists():
                print(f"    ✅ {fname}")
            else:
                print(f"    ❌ {fname} - 未生成")

        print(f"\nWikiPathways 数据库构建完成 → {outdir}")
        # 自动生成 GMT 文件
        self.generate_gmt_files(species, str(outdir))
        return str(outdir)

    # ============================
    # TRRUST 数据库构建
    # ============================
    def build_trrust(self, species: str, taxid: int) -> str:
        """构建指定物种的 TRRUST 数据库

        从 database/basic/trrust/TRRUSTv2/ 中读取 TRRUST 原始 TSV 文件，
        使用 TRRUSTParser 生成标准格式的数据库文件。

        TRRUST 仅支持 Human (hsa) 和 Mouse (mmu)。

        Args:
            species: 物种缩写（如 hsa, mmu）
            taxid: NCBI 物种分类学 ID（如 9606）

        Returns:
            str: 输出目录路径

        Raises:
            ValueError: 物种不被 TRRUST 支持时抛出
            FileNotFoundError: 基础数据文件不存在时抛出
        """
        # 检查物种支持
        latin_name = TRRUSTFetcher.get_latin_name(species)
        if latin_name is None:
            raise ValueError(
                f"TRRUST 不支持物种 '{species}'，"
                f"支持的物种: {TRRUSTFetcher.get_supported_species()}"
            )

        # 查找 TRRUST 原始数据文件
        trrust_raw = self.basic_dir / "trrust" / "TRRUSTv2" / f"trrust_rawdata.{species}.tsv"
        if not trrust_raw.exists():
            raise FileNotFoundError(
                f"TRRUST 原始数据文件不存在: {trrust_raw}\n"
                f"请先运行 download 下载数据：\n"
                f"  allenricher download trrust"
            )

        # 输出目录
        date_str = datetime.now().strftime("%Y%m%d")
        outdir = self.organism_dir / f"v{date_str}" / species
        outdir.mkdir(parents=True, exist_ok=True)

        # 尝试加载有效基因集合用于过滤
        gene_info_path = self._get_gene_info_path()
        valid_genes = None
        if gene_info_path:
            valid_genes = self._load_valid_genes(str(gene_info_path), taxid)

        print(f"\n{'='*60}")
        print(f"构建 TRRUST 数据库: {species} (taxid={taxid})")
        print(f"拉丁名: {latin_name}")
        print(f"数据源: {trrust_raw}")
        print(f"输出目录: {outdir}")
        print(f"{'='*60}")

        # 使用 TRRUSTParser 构建数据库
        print("|--- 从 TRRUST TSV 提取 TF-target 转录调控关系...")
        TRRUSTParser.build_database(
            tsv_path=str(trrust_raw),
            output_dir=str(outdir),
            species=species,
            valid_genes=valid_genes
        )

        # 验证输出文件
        expected_files = [
            f"{species}.TF2target.tab.gz",
            f"{species}.gene2TF.tab.gz",
            f"{species}.TF2disc.gz",
        ]
        for fname in expected_files:
            fpath = outdir / fname
            if fpath.exists():
                print(f"    [OK] {fname}")
            else:
                print(f"    [MISSING] {fname} - 未生成")

        print(f"\nTRRUST 数据库构建完成 -> {outdir}")
        # 自动生成 GMT 文件
        self.generate_gmt_files(species, str(outdir))
        return str(outdir)

    # ============================
    # ChEA3 数据库构建
    # ============================
    def build_chea3(self, species: str, taxid: int,
                    merge_method: str = "union") -> str:
        """构建指定物种的 ChEA3 数据库

        从 database/basic/chea3/ChEA3v2024/ 中读取所有 *_tf.gmt 文件，
        使用 ChEA3Parser 合并后生成标准格式的数据库文件。

        ChEA3 支持 Human (hsa)、Mouse (mmu)、Rat (rno)。

        Args:
            species: 物种缩写（如 hsa, mmu, rno）
            taxid: NCBI 物种分类学 ID（如 9606）
            merge_method: 合并方法，"union"（并集）或 "intersection"（交集）

        Returns:
            str: 输出目录路径

        Raises:
            FileNotFoundError: 基础数据目录或 GMT 文件不存在时抛出
        """
        # 添加物种警告
        if species != "hsa":
            print(f"|--- [警告] ChEA3 数据主要针对人类 (hsa)，当前物种 {species} 的结果可能不准确")
            print(f"|--- 建议仅将 ChEA3 用于人类数据分析")

        # 查找 ChEA3 GMT 文件
        chea3_dir = self.basic_dir / "chea3" / "ChEA3v2024"
        if not chea3_dir.exists():
            raise FileNotFoundError(
                f"ChEA3 基础数据目录不存在: {chea3_dir}\n"
                f"请先运行 download 下载数据：\n"
                f"  allenricher download chea3"
            )

        gmt_files = sorted(chea3_dir.glob("*_tf.gmt"))
        if not gmt_files:
            raise FileNotFoundError(
                f"在 {chea3_dir} 中未找到 *_tf.gmt 文件。\n"
                f"请先运行 download 下载数据：\n"
                f"  allenricher download chea3"
            )

        # 输出目录
        date_str = datetime.now().strftime("%Y%m%d")
        outdir = self.organism_dir / f"v{date_str}" / species
        outdir.mkdir(parents=True, exist_ok=True)

        # 获取有效基因集合（从 gene_info）
        gene_info_path = self._get_gene_info_path()
        valid_genes = None
        if gene_info_path:
            valid_genes = self._load_valid_genes(str(gene_info_path), taxid)
            if valid_genes:
                print(f"|--- 有效基因过滤: {len(valid_genes)} 个基因")

        print(f"\n{'='*60}")
        print(f"构建 ChEA3 数据库: {species} (taxid={taxid})")
        print(f"数据源: {chea3_dir} ({len(gmt_files)} 个 GMT 库)")
        print(f"合并方法: {merge_method}")
        print(f"输出目录: {outdir}")
        print(f"{'='*60}")

        # 使用 ChEA3Parser 构建数据库
        print("|--- 从 ChEA3 GMT 文件提取 TF-target 关系...")
        ChEA3Parser.build_database(
            gmt_paths=[str(f) for f in gmt_files],
            output_dir=str(outdir),
            species=species,
            merge_method=merge_method,
            valid_genes=valid_genes  # 确保传递
        )

        # 验证输出文件
        expected_files = [
            f"{species}.ChEA3_2gene.tab.gz",
            f"{species}.ChEA3_2disc.gz",
        ]
        for fname in expected_files:
            fpath = outdir / fname
            if fpath.exists():
                print(f"    [OK] {fname}")
            else:
                print(f"    [MISSING] {fname} - 未生成")

        print(f"\nChEA3 数据库构建完成 -> {outdir}")
        # 自动生成 GMT 文件
        self.generate_gmt_files(species, str(outdir))
        return str(outdir)

    # ============================
    # GMT 基因集文件生成
    # ============================
    def generate_gmt_files(self, species: str,
                           output_dir: str = None) -> Dict[str, str]:
        """从物种数据库产物生成 GMT 基因集文件

        扫描物种数据库目录中的所有可用数据库产物，
        自动生成对应的 .gmt.gz 文件供 GSEA/ssGSEA/GSVA 使用。

        Args:
            species: 物种缩写（如 hsa）
            output_dir: 物种数据库目录路径，默认自动检测最新版本

        Returns:
            Dict[str, str]: {数据库名称: GMT文件路径}，仅包含成功生成的条目
        """
        if output_dir is None:
            # 自动查找最新的物种数据库目录
            if not self.organism_dir.exists():
                print("|--- [警告] 物种数据库目录不存在，跳过 GMT 生成")
                return {}
            species_dirs = sorted(
                self.organism_dir.glob(f"*/{species}"),
                key=lambda p: p.parent.name,
                reverse=True
            )
            if not species_dirs:
                print(f"|--- [警告] 未找到物种 {species} 的数据库目录，跳过 GMT 生成")
                return {}
            output_dir = str(species_dirs[0])

        generator = GMTGenerator(organism_dir=output_dir)
        return generator.generate_all_gmt(species)

    # ============================
    # hTFtarget 数据库构建（人类专用）
    # ============================
    def build_htftarget(self, species: str, taxid: int) -> str:
        """构建 hTFtarget 数据库（人类专用）

        Args:
            species: 物种代码（如 hsa）
            taxid: NCBI TaxID

        Returns:
            输出目录路径
        """
        from .parsers.htftarget import HTFtargetParser

        htftarget_file = self.basic_dir / "htftarget" / "tf-target-information.txt"

        if not htftarget_file.exists():
            print(f"|--- [跳过] 未找到 hTFtarget 文件: {htftarget_file}")
            print("|--- 请先运行: allenricher download --animaltfdb")
            return ""

        date_str = datetime.now().strftime("%Y%m%d")
        outdir = self.organism_dir / f"v{date_str}" / species
        outdir.mkdir(parents=True, exist_ok=True)

        print(f"\n{'='*60}")
        print(f"构建 hTFtarget 数据库: {species} (taxid={taxid})")
        print(f"数据源: {htftarget_file}")
        print(f"输出目录: {outdir}")
        print(f"{'='*60}")

        gene_info_path = self._get_gene_info_path()
        valid_genes = None
        if gene_info_path:
            valid_genes = self._load_valid_genes(str(gene_info_path), taxid)

        HTFtargetParser.build_database(
            tsv_path=str(htftarget_file),
            output_dir=str(outdir),
            species=species,
            valid_genes=valid_genes,
        )

        expected_files = [
            f"{species}.hTF_2gene.tab.gz",
            f"{species}.hTF_2disc.gz",
        ]
        for fname in expected_files:
            fpath = outdir / fname
            if fpath.exists():
                print(f"    ✅ {fname}")
            else:
                print(f"    ❌ {fname} - 未生成")

        # 对于人类，也更新注册表
        if species == 'hsa':
            try:
                # 从 hTFtarget 文件获取 TF 统计
                human_tf_to_targets, _, _ = HTFtargetParser.parse_tsv(str(htftarget_file))
                tf_count = len(human_tf_to_targets) if human_tf_to_targets else 0
                
                registry = SpeciesRegistry.load_default(str(self.root_dir))
                registry.update_animaltfdb_stats(
                    species_code='hsa',
                    tf_count=tf_count,
                    mapped_target_count=0,  # hTFtarget 直接数据，非映射
                    has_data=True
                )
                registry.save()
                print(f"✅ 物种注册表已更新: hsa -> hTFtarget")
            except Exception as e:
                print(f"⚠️  更新物种注册表失败: {e}")

        print(f"\nhTFtarget 数据库构建完成 → {outdir}")
        return str(outdir)

    # ============================
    # AnimalTFDB 数据库构建（同源映射）
    # ============================
    def build_animaltfdb(self, species: str, taxid: int,
                         species_latin: str = "") -> str:
        """构建 AnimalTFDB 数据库（通过同源映射）

        对于人类(hsa)：直接构建 hTFtarget 数据库
        对于其他物种：通过同源映射推断TF-target关系

        Args:
            species: 物种代码（如 bta 代表牛）
            taxid: NCBI TaxID
            species_latin: 物种拉丁名（下划线格式，如 Bos_taurus）

        Returns:
            输出目录路径
        """
        from .parsers.animaltfdb import AnimalTFDBParser
        from .parsers.htftarget import HTFtargetParser
        from .ortholog_mapper import OrthologMapper

        if species == 'hsa':
            return self.build_htftarget(species, taxid)

        if not species_latin:
            from .species_registry import SpeciesRegistry
            registry = SpeciesRegistry.load_default(str(self.root_dir))
            entry = registry.query_by_kegg_code(species)
            if entry and entry.latin_name:
                species_latin = entry.latin_name.replace(' ', '_')
            else:
                print(f"|--- [错误] 无法确定物种 {species} 的拉丁名，请使用 --latin-name 参数")
                return ""

        cache_dir = self.basic_dir / "animaltfdb" / "AnimalTFDBv4.0"
        tf_list_file = cache_dir / f"{species_latin}_TF"
        ortholog_file = cache_dir / f"{species_latin}_ortholog_to_human"
        htftarget_file = self.basic_dir / "htftarget" / "tf-target-information.txt"

        missing = []
        if not tf_list_file.exists():
            missing.append(f"TF列表: {tf_list_file}")
        if not ortholog_file.exists():
            missing.append(f"同源映射: {ortholog_file}")
        if not htftarget_file.exists():
            missing.append(f"hTFtarget: {htftarget_file}")

        if missing:
            print(f"|--- [跳过] 缺少必要文件:")
            for m in missing:
                print(f"|---   - {m}")
            print(f"|--- 请先运行: allenricher download --animaltfdb --species {species_latin}")
            return ""

        date_str = datetime.now().strftime("%Y%m%d")
        outdir = self.organism_dir / f"v{date_str}" / species
        outdir.mkdir(parents=True, exist_ok=True)

        print(f"\n{'='*60}")
        print(f"构建 AnimalTFDB 数据库: {species} ({species_latin}, taxid={taxid})")
        print(f"策略: 同源映射 (hTFtarget → {species_latin})")
        print(f"输出目录: {outdir}")
        print(f"{'='*60}")

        gene_info_path = self._get_gene_info_path()
        valid_genes = None
        if gene_info_path:
            valid_genes = self._load_valid_genes(str(gene_info_path), taxid)

        print("\n[1/3] 解析 AnimalTFDB 数据...")
        tf_df, ortholog_map = AnimalTFDBParser.build_database(
            tf_list_path=str(tf_list_file),
            ortholog_path=str(ortholog_file),
            output_dir=str(outdir),
            species=species,
            valid_genes=valid_genes,
        )

        print("[2/3] 解析 hTFtarget 人类TF-target关系...")
        human_tf_to_targets, _, _ = HTFtargetParser.parse_tsv(str(htftarget_file))

        print("[3/3] 执行同源映射...")
        species_tf_set = set(tf_df['Symbol'].values) if 'Symbol' in tf_df.columns else None

        mapper = OrthologMapper(
            human_tf_to_targets=human_tf_to_targets,
            species_to_human=ortholog_map,
            species_tf_set=species_tf_set,
        )

        tf_to_targets, gene_to_tfs = mapper.map_tf_targets()

        # 获取映射统计
        dedup_stats = mapper.get_duplicate_stats()
        mapping_stats = {
            'species': species,
            'species_latin': species_latin,
            'total_species_genes': len(ortholog_map),
            'total_human_tfs': len(human_tf_to_targets),
            'species_tfs_found': len(species_tf_set) if species_tf_set else 0,
            'mapped_tfs': len(tf_to_targets),
            'mapped_targets': len(gene_to_tfs),
            'avg_targets_per_tf': sum(len(t) for t in tf_to_targets.values()) / len(tf_to_targets) if tf_to_targets else 0,
            'multi_mapping_human_genes': dedup_stats.get('multi_mapping_count', 0),
            'coverage_ratio': len(tf_to_targets) / (len(species_tf_set) if species_tf_set else 1) * 100,
        }

        # 输出统计报告
        print("\n" + "="*60)
        print("映射质量统计报告")
        print("="*60)
        print(f"物种: {species_latin} ({species})")
        print(f"同源映射基因数: {mapping_stats['total_species_genes']}")
        print(f"人类TF总数: {mapping_stats['total_human_tfs']}")
        print(f"物种TF数: {mapping_stats['species_tfs_found']}")
        print(f"成功映射TF数: {mapping_stats['mapped_tfs']} ({mapping_stats['coverage_ratio']:.1f}% 覆盖率)")
        print(f"推断靶基因数: {mapping_stats['mapped_targets']}")
        print(f"平均每TF靶基因数: {mapping_stats['avg_targets_per_tf']:.1f}")
        print(f"多对一映射人类基因数: {mapping_stats['multi_mapping_human_genes']}")
        print("="*60)

        # 低覆盖率警告
        if mapping_stats['coverage_ratio'] < 50:
            print(f"\n警告: TF映射覆盖率仅 {mapping_stats['coverage_ratio']:.1f}%，")
            print(f"   可能原因: 该物种与人类亲缘关系较远，直系同源映射有限。")
            print(f"   建议: 考虑使用 motif 扫描方法补充预测。")

        # 保存统计到文件
        import json
        stats_file = outdir / f"{species}.AnimalTFDB_mapping_stats.json"
        with open(stats_file, 'w', encoding='utf-8') as f:
            json.dump(mapping_stats, f, indent=2, ensure_ascii=False)
        print(f"统计报告已保存: {stats_file}")

        if not tf_to_targets:
            print("|--- [警告] 同源映射结果为空，无法构建TF-target关系")
            return str(outdir)

        OrthologMapper.build_mapped_database(
            tf_to_targets=tf_to_targets,
            gene_to_tfs=gene_to_tfs,
            species_tf_df=tf_df,
            output_dir=str(outdir),
            species=species,
        )

        expected_files = [
            f"{species}.AnimalTFDB_2tf.tab.gz",
            f"{species}.AnimalTFDB_2disc.gz",
            f"{species}.AnimalTFDB_ortholog.gz",
            f"{species}.AnimalTFDB_2gene.tab.gz",
            f"{species}.AnimalTFDB_mapped_2disc.gz",
        ]
        for fname in expected_files:
            fpath = outdir / fname
            if fpath.exists():
                print(f"    ✅ {fname}")
            else:
                print(f"    ❌ {fname} - 未生成")

        # 自动更新物种注册表
        try:
            registry = SpeciesRegistry.load_default(str(self.root_dir))
            
            # 统计数据
            tf_count = len(tf_to_targets) if tf_to_targets else 0
            mapped_target_count = len(gene_to_tfs) if gene_to_tfs else 0
            
            registry.update_animaltfdb_stats(
                species_code=species,
                tf_count=tf_count,
                mapped_target_count=mapped_target_count,
                has_data=True
            )
            registry.save()
            print(f"✅ 物种注册表已更新: {species} -> AnimalTFDB")
        except Exception as e:
            print(f"⚠️  更新物种注册表失败: {e}")

        print(f"\nAnimalTFDB 数据库构建完成 → {outdir}")
        return str(outdir)

    # ============================
    # 一键构建（对应 make_speciesDB）
    # ============================
    def build_species_db(self, species: str, taxid: int,
                         databases: List[str] = None,
                         go_version: Optional[str] = None,
                         reactome_version: Optional[str] = None,
                         **kwargs) -> str:
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
                    # 获取 latin_name，优先从 SpeciesRegistry 查询
                    latin_name = species  # 默认回退
                    try:
                        from ..core.config import SPECIES_CONFIGS
                        if species in SPECIES_CONFIGS:
                            latin_name = SPECIES_CONFIGS[species].name
                    except:
                        pass

                    try:
                        from .species_registry import SpeciesRegistry
                        registry = SpeciesRegistry.load_default()
                        entry = registry.query_by_taxid(taxid)
                        if entry:
                            latin_name = entry.latin_name
                    except:
                        pass

                    # 使用回退构建
                    self.build_go_with_fallback(taxid, latin_name, go_version)
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
                elif db_upper == 'WIKIPATHWAYS':
                    self.build_wikipathways(species, taxid)
                elif db_upper == 'TRRUST':
                    self.build_trrust(species, taxid)
                elif db_upper == 'CHEA3':
                    self.build_chea3(species, taxid)
                elif db_upper == 'ANIMALTFDB':
                    latin_name = kwargs.get('latin_name', '')
                    self.build_animaltfdb(species, taxid, species_latin=latin_name)
                elif db_upper == 'HTFTARGET':
                    self.build_htftarget(species, taxid)
                else:
                    print(f"|--- [警告] 未知数据库: {db_name}")
            except Exception as e:
                print(f"|--- [错误] 构建 {db_name} 失败: {e}")

        print(f"\n{'#'*60}")
        print(f"# 物种数据库构建完成 → {outdir}")
        print(f"{'#'*60}")

        # 写入构建清单（血缘追踪）
        import json as _json
        from datetime import timezone as _tz

        try:
            _downloader = DataDownloader(root_dir=str(self.root_dir))

            _build_manifest = {
                "schema_version": "1.0",
                "built_at": datetime.now(_tz.utc).isoformat(),
                "allenricher_version": __import__("allenricher").__version__,
                "species": species,
                "taxid": taxid,
                "databases": sorted(databases),
                "dependencies": {},
                "source_versions": {},
            }

            # 从 DatabaseVersionManager 填充 source_versions
            try:
                from allenricher.database.version import DatabaseVersionManager
                _vm = DatabaseVersionManager(database_dir=str(self.root_dir))

                if "GO" in databases:
                    _go_obo_ver = _vm.get_local_version("go_obo")
                    if _go_obo_ver and _go_obo_ver.remote_version:
                        _build_manifest["source_versions"]["go_obo"] = _go_obo_ver.remote_version
                    _gene2go_ver = _vm.get_local_version("gene2go")
                    if _gene2go_ver and _gene2go_ver.remote_last_modified:
                        _build_manifest["source_versions"]["gene2go"] = _gene2go_ver.remote_last_modified

                if "Reactome" in databases:
                    _reactome_ver = _vm.get_local_version("reactome")
                    if _reactome_ver and _reactome_ver.remote_version:
                        _build_manifest["source_versions"]["reactome"] = _reactome_ver.remote_version

                if "KEGG" in databases:
                    _kegg_ver = _vm.get_local_version("kegg")
                    if _kegg_ver and _kegg_ver.remote_version:
                        _build_manifest["source_versions"]["kegg"] = _kegg_ver.remote_version
            except Exception as _sv_err:
                logger.warning("填充 source_versions 失败: %s", _sv_err)

            # GO 依赖
            if "GO" in databases:
                _go_ver = go_version if go_version else _downloader.get_latest_go_version()
                _build_manifest["dependencies"]["GO"] = {
                    "basic_dir": f"basic/go/{_go_ver}",
                    "files": ["gene2go.gz", "gene_info.gz", "go-basic.obo"],
                }

            # Reactome 依赖
            if "Reactome" in databases:
                _reactome_ver = reactome_version if reactome_version else _downloader.get_latest_reactome_version()
                _build_manifest["dependencies"]["Reactome"] = {
                    "basic_dir": f"basic/reactome/{_reactome_ver}",
                    "files": ["NCBI2Reactome_All_Levels.txt.gz", "gene_info.gz"],
                }

            # KEGG 依赖
            if "KEGG" in databases:
                _build_manifest["dependencies"]["KEGG"] = {
                    "source": "REST API (real-time)",
                    "gene_info_from": f"basic/go/{go_version or _downloader.get_latest_go_version()}/gene_info.gz",
                }

            # DO 依赖
            if "DO" in databases:
                _build_manifest["dependencies"]["DO"] = {
                    "basic_dir": "basic/do",
                    "files": ["human_disease_knowledge_filtered.tsv.gz",
                               "human_disease_experiments_filtered.tsv.gz",
                               "human_disease_textmining_filtered.tsv.gz"],
                }

            _manifest_path = outdir / "build_manifest.json"
            with open(_manifest_path, "w", encoding="utf-8") as _mf:
                _json.dump(_build_manifest, _mf, indent=2, ensure_ascii=False)
            logger.info("构建清单已保存: %s", _manifest_path)
        except Exception as _e:
            logger.warning("写入构建清单失败: %s", _e)

        return str(outdir)
