"""
数据下载器模块（重构版）

对应 v1 脚本：update_GOdb, update_ReactomeDB

从指定数据源下载全体物种的通用原始数据文件，存入 database/basic/ 目录。
使用 DownloadManager 实现多线程加速、镜像源自动切换、完整性校验。

架构（与 v1 一致）：
  database/basic/
    go/GO{date}/
      gene2go.gz          ← 全员基因-GO 映射
      gene_info.gz        ← 全员基因信息
      go-basic.obo        ← GO 本体定义
    reactome/Reactome{date}/
      gene_info.gz        ← 全员基因信息
      NCBI2Reactome_All_Levels.txt.gz  ← 全员通路数据

用法：
  downloader = DataDownloader(root_dir="./database")
  go_dir = downloader.download_go_basic()
  re_dir = downloader.download_reactome_basic()
"""

from __future__ import annotations

import csv
import gzip
import logging
import re
import shutil
import sys
import time
from pathlib import Path
from datetime import datetime
from typing import Dict, List, Optional

import requests

from .download_manager import DownloadManager
from .mirrors import get_mirrors, JENSEN_SOURCES
from .species_registry import SpeciesRegistry, SpeciesEntry

logger = logging.getLogger(__name__)


class DataDownloader:
    """全体物种通用数据下载器（重构版）

    使用 DownloadManager 提供多线程下载、镜像源切换、完整性校验。

    Attributes:
        root_dir: 数据库根目录
        basic_dir: 基础数据目录
        manager: DownloadManager 实例
    """

    def __init__(
        self,
        root_dir: str = "./database",
        overwrite: bool = False,
        max_workers: int = 4,
        use_multi_thread: bool = True,
        verify_integrity: bool = True,
    ):
        """初始化下载器

        Args:
            root_dir: 数据库根目录
            overwrite: 是否覆盖已存在文件
            max_workers: 多线程下载线程数
            use_multi_thread: 是否启用多线程下载大文件
            verify_integrity: 是否验证下载文件完整性
        """
        self.root_dir = Path(root_dir)
        self.basic_dir = self.root_dir / "basic"
        self.basic_dir.mkdir(parents=True, exist_ok=True)

        self.manager = DownloadManager(
            root_dir=root_dir,
            overwrite=overwrite,
            max_workers=max_workers,
            use_multi_thread=use_multi_thread,
            verify_integrity=verify_integrity,
            show_progress=True,
        )

    # ============================
    # GO 基础数据下载（全体物种）
    # ============================
    def download_go_basic(self, version: Optional[str] = None) -> str:
        """下载 GO 全体物种基础数据

        从 NCBI 镜像源下载 gene2go.gz、gene_info.gz，
        从 GO 镜像源下载 go-basic.obo。
        大文件自动使用多线程加速。

        Args:
            version: 版本号（如 "GO20250101"），默认当前日期

        Returns:
            GO 基础数据目录路径
        """
        if version is None:
            version = f"GO{datetime.now().strftime('%Y%m%d')}"

        go_dir = self.basic_dir / "go" / version
        go_dir.mkdir(parents=True, exist_ok=True)

        print(f"\n{'='*60}")
        print(f"下载 GO 基础数据 → {go_dir}")
        print(f"{'='*60}")

        ncbi_mirrors = get_mirrors('ncbi')

        self.manager.download_with_mirror_fallback(
            ncbi_mirrors, "gene2go.gz",
            go_dir / "gene2go.gz", desc="gene2go.gz"
        )
        self.manager.download_with_mirror_fallback(
            ncbi_mirrors, "gene_info.gz",
            go_dir / "gene_info.gz", desc="gene_info.gz"
        )

        go_mirrors = get_mirrors('go')
        self.manager.download_with_mirror_fallback(
            go_mirrors, "go-basic.obo",
            go_dir / "go-basic.obo", desc="go-basic.obo"
        )

        print(f"GO 基础数据下载完成 → {go_dir}")

        # 记录版本元数据到 versions.json
        try:
            from allenricher.database.version import DatabaseVersionManager, RemoteVersionChecker
            _vm = DatabaseVersionManager(database_dir=str(self.root_dir))
            _checker = RemoteVersionChecker()

            # 记录 gene2go
            _g2g_info = _checker.check_head("https://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2go.gz")
            if _g2g_info:
                _vm.record_download(
                    source="gene2go", local_version=version,
                    local_path=f"basic/go/{version}",
                    remote_last_modified=_g2g_info.get("last_modified"),
                )

            # 记录 go_obo
            _obo_info = _checker.check_go_obo_version()
            if _obo_info:
                _vm.record_download(
                    source="go_obo", local_version=version,
                    local_path=f"basic/go/{version}/go-basic.obo",
                    remote_version=_obo_info.get("remote_version"),
                    remote_last_modified=_obo_info.get("last_modified"),
                )

            # 记录 gene_info
            _gi_info = _checker.check_head("https://ftp.ncbi.nlm.nih.gov/gene/DATA/gene_info.gz")
            if _gi_info:
                _vm.record_download(
                    source="gene_info", local_version=version,
                    local_path=f"basic/go/{version}",
                    remote_last_modified=_gi_info.get("last_modified"),
                )
        except Exception as _e:
            logger.warning("记录 GO 版本元数据失败: %s", _e)

        return str(go_dir)

    # ============================
    # Reactome 基础数据下载（全体物种）
    # ============================
    def download_reactome_basic(self, version: Optional[str] = None,
                                  go_version: Optional[str] = None) -> str:
        """下载 Reactome 全体物种基础数据

        Args:
            version: 版本号，默认当前日期
            go_version: GO 版本号（用于复用 gene_info.gz）

        Returns:
            Reactome 基础数据目录路径
        """
        if version is None:
            version = f"Reactome{datetime.now().strftime('%Y%m%d')}"

        re_dir = self.basic_dir / "reactome" / version
        re_dir.mkdir(parents=True, exist_ok=True)

        print(f"\n{'='*60}")
        print(f"下载 Reactome 基础数据 → {re_dir}")
        print(f"{'='*60}")

        # gene_info.gz：优先复用 GO 已下载的文件
        re_gene_info = re_dir / "gene_info.gz"
        if re_gene_info.exists() and not self.manager.overwrite:
            if self.manager.verify_integrity:
                from .download_utils import verify_gzip_integrity
                valid, _ = verify_gzip_integrity(re_gene_info)
                if valid:
                    print(f"|--- 已存在且有效，跳过: gene_info.gz")
        else:
            # 尝试从 GO 目录复制
            if go_version is None:
                go_version = self.get_latest_go_version()
            if go_version:
                go_gene_info = self.basic_dir / "go" / go_version / "gene_info.gz"
                if go_gene_info.exists():
                    print(f"|--- 复用 GO 数据: {go_gene_info}")
                    shutil.copy2(go_gene_info, re_gene_info)
                    print(f"|--- 已复制: gene_info.gz ({re_gene_info.stat().st_size / 1024 / 1024:.1f} MB)")
                else:
                    # GO 目录也没有，从镜像下载
                    ncbi_mirrors = get_mirrors('ncbi')
                    self.manager.download_with_mirror_fallback(
                        ncbi_mirrors, "gene_info.gz",
                        re_gene_info, desc="gene_info.gz"
                    )
            else:
                # 没有 GO 版本，从镜像下载
                ncbi_mirrors = get_mirrors('ncbi')
                self.manager.download_with_mirror_fallback(
                    ncbi_mirrors, "gene_info.gz",
                    re_gene_info, desc="gene_info.gz"
                )

        # NCBI2Reactome（下载后 gzip 压缩）
        reactome_mirrors = get_mirrors('reactome')
        raw_file = re_dir / "NCBI2Reactome_All_Levels.txt"
        gz_file = re_dir / "NCBI2Reactome_All_Levels.txt.gz"

        if gz_file.exists() and not self.manager.overwrite:
            if self.manager.verify_integrity:
                from .download_utils import verify_gzip_integrity
                valid, _ = verify_gzip_integrity(gz_file)
                if valid:
                    print(f"|--- 已存在且有效，跳过: {gz_file.name}")
                    print(f"Reactome 基础数据下载完成 → {re_dir}")

                    self._record_reactome_version(version)

                    return str(re_dir)

        self.manager.download_with_mirror_fallback(
            reactome_mirrors, "NCBI2Reactome_All_Levels.txt",
            raw_file, desc="NCBI2Reactome"
        )

        if raw_file.exists():
            print(f"|--- 压缩: {raw_file.name} → {gz_file.name}")
            with open(raw_file, 'rb') as f_in:
                with gzip.open(gz_file, 'wb') as f_out:
                    shutil.copyfileobj(f_in, f_out)
            raw_file.unlink()

        print(f"Reactome 基础数据下载完成 → {re_dir}")

        self._record_reactome_version(version)

        return str(re_dir)

    def _record_reactome_version(self, version: str) -> None:
        """记录 Reactome 版本元数据到 versions.json"""
        try:
            from allenricher.database.version import DatabaseVersionManager, RemoteVersionChecker
            _vm = DatabaseVersionManager(database_dir=str(self.root_dir))
            _checker = RemoteVersionChecker()
            _re_info = _checker.check_reactome_version()
            if _re_info:
                _vm.record_download(
                    source="reactome", local_version=version,
                    local_path=f"basic/reactome/{version}",
                    remote_version=_re_info.get("remote_version"),
                    remote_last_modified=_re_info.get("last_modified"),
                )
        except Exception as _e:
            logger.warning("记录 Reactome 版本元数据失败: %s", _e)

    # ============================
    # DO / DisGeNET（仅人类）
    # ============================
    def download_do_files(self) -> Dict[str, str]:
        """下载 Jensen Lab Disease Ontology 文件（仅人类）

        3 个 TSV 文件并行下载，下载后自动 gzip 压缩存储。
        builder 侧通过 glob ``human_disease_*_filtered.tsv.gz`` 读取。

        Returns:
            {来源名: 文件路径}
        """
        do_dir = self.basic_dir / "do"
        do_dir.mkdir(parents=True, exist_ok=True)

        print(f"\n{'='*60}")
        print(f"下载 DO 数据 → {do_dir}")
        print(f"{'='*60}")

        files = {}
        for url in JENSEN_SOURCES:
            fname = url.split('/')[-1]
            raw_dest = do_dir / fname
            gz_dest = do_dir / f"{fname}.gz"

            # 已存在且有效的压缩文件 → 跳过
            if gz_dest.exists() and not self.manager.overwrite:
                if self.manager.verify_integrity:
                    from .download_utils import verify_gzip_integrity
                    valid, _ = verify_gzip_integrity(gz_dest)
                    if valid:
                        print(f"|--- 已存在且有效，跳过: {gz_dest.name}")
                        files[fname] = str(gz_dest)
                        continue

            # 下载原始 TSV
            self.manager.download_file(url, raw_dest, desc=fname)

            # gzip 压缩
            if raw_dest.exists():
                print(f"|--- 压缩: {raw_dest.name} → {gz_dest.name}")
                with open(raw_dest, 'rb') as f_in:
                    with gzip.open(gz_dest, 'wb') as f_out:
                        shutil.copyfileobj(f_in, f_out)
                raw_dest.unlink()
                files[fname] = str(gz_dest)
            else:
                files[fname] = str(raw_dest)

        # 记录 DO 版本元数据到 versions.json
        try:
            from allenricher.database.version import DatabaseVersionManager
            _vm = DatabaseVersionManager(database_dir=str(self.root_dir))
            _vm.record_download(
                source="do",
                local_version="cached",
                local_path="basic/do",
            )
        except Exception as _e:
            logger.warning("记录 DO 版本元数据失败: %s", _e)

        return files

    def download_disgenet(self) -> str:
        """下载 DisGeNET 数据（仅人类）

        .. warning::
            DisGeNET 已迁移至商业平台 (disgenet.com)，
            旧 URL 已失效。此方法保留接口兼容性。

        Returns:
            文件路径
        """
        dg_dir = self.basic_dir / "disgenet"
        dg_dir.mkdir(parents=True, exist_ok=True)

        print(f"\n{'='*60}")
        print(f"下载 DisGeNET 数据 → {dg_dir}")
        print(f"{'='*60}")

        url = (
            "http://www.disgenet.org/static/disgenet_ap1/files/"
            "downloads/all_gene_disease_associations.tsv.gz"
        )
        dest = dg_dir / "all_gene_disease_associations.tsv.gz"

        print("|--- [警告] DisGeNET 旧 URL 已失效，跳过下载")
        print("|--- 替代方案: 使用 CTD (Comparative Toxicogenomics Database)")
        return str(dest)

    # ============================
    # 批量下载
    # ============================
    def download_all(self, db_types: List[str] = None) -> Dict[str, str]:
        """下载所有基础数据

        Args:
            db_types: 要下载的类型列表，默认 ['go', 'reactome']

        Returns:
            {数据库类型: 目录路径}
        """
        if db_types is None:
            db_types = ['go', 'reactome']

        result = {}
        go_version = None  # 记录 GO 版本，供 Reactome 复用

        for db_type in db_types:
            db_type = db_type.lower().strip()
            try:
                if db_type in ('go',):
                    result['go'] = self.download_go_basic()
                    go_version = self.get_latest_go_version()  # 记录版本
                elif db_type in ('reactome',):
                    result['reactome'] = self.download_reactome_basic(go_version=go_version)
                elif db_type in ('do',):
                    result['do'] = str(self.basic_dir / "do")
                    self.download_do_files()
                elif db_type in ('kegg',):
                    # KEGG 通过 REST API 获取，需要 gene_info.gz
                    result['kegg'] = str(self.basic_dir / "kegg")
                    print("|--- KEGG 数据将在 build 阶段通过 REST API 自动获取")
                elif db_type in ('disgenet',):
                    result['disgenet'] = self.download_disgenet()
                else:
                    print(f"|--- [警告] 未知数据库类型: {db_type}")
            except Exception as e:
                print(f"|--- [错误] 下载 {db_type} 失败: {e}")

        # ============================
        # 注册表构建流水线
        # ============================
        print(f"\n{'='*60}")
        print("构建物种注册表")
        print(f"{'='*60}")

        go_registry: Optional[Path] = None
        goa_index: Optional[Path] = None
        kegg_registry: Optional[Path] = None
        reactome_registry: Optional[Path] = None
        do_registry: Optional[Path] = None

        # GO 注册表构建（如果下载了 GO）
        if 'go' in result:
            go_dir = Path(result['go'])
            gene2go_path = go_dir / "gene2go.gz"

            if gene2go_path.exists():
                try:
                    go_registry = self._build_go_registry(gene2go_path, go_dir)
                    logger.info(f"GO registry built: {go_registry}")
                except Exception as e:
                    logger.warning(f"Failed to build GO registry: {e}")

                try:
                    goa_index = self._download_goa_index(go_dir)
                    logger.info(f"GOA index downloaded: {goa_index}")
                except Exception as e:
                    logger.warning(f"Failed to download GOA index: {e}")

                if go_registry and goa_index:
                    try:
                        go_registry = self._merge_go_registries(go_registry, goa_index, go_dir)
                        logger.info(f"GO registries merged: {go_registry}")
                    except Exception as e:
                        logger.warning(f"Failed to merge GO registries: {e}")

        # KEGG 注册表构建（如果下载了 KEGG）
        if 'kegg' in result:
            kegg_dir = Path(result['kegg'])
            try:
                kegg_registry = self._build_kegg_registry(kegg_dir)
                logger.info(f"KEGG registry built: {kegg_registry}")
            except Exception as e:
                logger.warning(f"Failed to build KEGG registry: {e}")

        # Reactome 注册表构建（如果下载了 Reactome）
        if 'reactome' in result:
            reactome_dir = Path(result['reactome'])
            ncbi2reactome_path = reactome_dir / "NCBI2Reactome_All_Levels.txt.gz"

            if ncbi2reactome_path.exists():
                try:
                    reactome_registry = self._build_reactome_registry(ncbi2reactome_path, reactome_dir)
                    logger.info(f"Reactome registry built: {reactome_registry}")
                except Exception as e:
                    logger.warning(f"Failed to build Reactome registry: {e}")

        # DO 注册表构建（如果下载了 DO）
        if 'do' in result:
            do_dir = Path(result['do'])
            try:
                do_registry = self._build_do_registry(do_dir)
                logger.info(f"DO registry built: {do_registry}")
            except Exception as e:
                logger.warning(f"Failed to build DO registry: {e}")

        # 合并所有注册表
        supported_species_path = self.basic_dir / "supported_species.tsv"
        try:
            self._merge_all_registries(
                go_registry=go_registry or Path("/dev/null"),
                kegg_registry=kegg_registry or Path("/dev/null"),
                reactome_registry=reactome_registry or Path("/dev/null"),
                do_registry=do_registry or Path("/dev/null"),
                output_path=supported_species_path,
            )
            logger.info(f"All registries merged: {supported_species_path}")
            result['supported_species'] = str(supported_species_path)
        except Exception as e:
            logger.warning(f"Failed to merge all registries: {e}")

        # 打印下载统计摘要
        try:
            self._report_download_summary(supported_species_path)
        except Exception as e:
            logger.warning(f"Failed to report download summary: {e}")

        return result

    # ============================
    # 版本管理
    # ============================
    def list_go_versions(self) -> list:
        """列出已下载的 GO 基础数据版本"""
        go_basic = self.basic_dir / "go"
        if not go_basic.exists():
            return []
        return sorted([d.name for d in go_basic.iterdir() if d.is_dir()])

    def list_reactome_versions(self) -> list:
        """列出已下载的 Reactome 基础数据版本"""
        re_basic = self.basic_dir / "reactome"
        if not re_basic.exists():
            return []
        return sorted([d.name for d in re_basic.iterdir() if d.is_dir()])

    def get_latest_go_version(self) -> Optional[str]:
        """获取最新的 GO 基础数据版本"""
        versions = self.list_go_versions()
        return versions[-1] if versions else None

    def get_latest_reactome_version(self) -> Optional[str]:
        """获取最新的 Reactome 基础数据版本"""
        versions = self.list_reactome_versions()
        return versions[-1] if versions else None

    def _download_taxonomy_names(self, output_dir: Path) -> Optional[Path]:
        """下载 NCBI Taxonomy 物种名称文件

        从 NCBI FTP 下载 taxdump.tar.gz，提取 names.dmp 文件，
        生成 taxid → scientific_name 映射。

        Args:
            output_dir: 输出目录

        Returns:
            生成的 taxid_to_name.tsv 文件路径
        """
        logger.info("下载 NCBI Taxonomy 物种名称...")
        taxonomy_dir = self.basic_dir / "taxonomy"
        taxonomy_dir.mkdir(parents=True, exist_ok=True)

        # 检查是否已有 names.dmp
        names_dmp = taxonomy_dir / "names.dmp"
        if names_dmp.exists() and not self.manager.overwrite:
            logger.info("NCBI Taxonomy 已存在，跳过下载")
        else:
            # 下载 taxdump.tar.gz
            taxdump_url = "https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz"
            taxdump_tgz = taxonomy_dir / "taxdump.tar.gz"

            logger.info("从 NCBI 下载 taxonomy 数据: %s", taxdump_url)
            try:
                response = requests.get(taxdump_url, stream=True, timeout=300)
                response.raise_for_status()

                with open(taxdump_tgz, 'wb') as f:
                    for chunk in response.iter_content(chunk_size=8192):
                        f.write(chunk)
                logger.info("下载完成: %s", taxdump_tgz)

                # 解压 names.dmp
                import tarfile
                with tarfile.open(taxdump_tgz, 'r:gz') as tar:
                    names_member = tar.getmember('names.dmp')
                    tar.extract(names_member, taxonomy_dir)
                    logger.info("解压 names.dmp 成功")
            except Exception as e:
                logger.warning("下载 NCBI Taxonomy 失败: %s", e)
                return None

        if not names_dmp.exists():
            return None

        # 解析 names.dmp 生成 taxid → latin_name 映射
        # names.dmp 格式: tax_id \t name \t unique name \t name class
        taxid_to_name: Dict[int, str] = {}
        seen_taxids: set = set()

        with open(names_dmp, 'r', encoding='utf-8', errors='ignore') as f:
            for line in f:
                line = line.rstrip('\t|\n')
                if not line:
                    continue
                parts = line.split('\t|\t')
                if len(parts) < 4:
                    continue
                try:
                    taxid = int(parts[0])
                    name_class = parts[3].rstrip('\t|')

                    # 只取 scientific name
                    if name_class == 'scientific name':
                        taxid_to_name[taxid] = parts[1]
                except (ValueError, IndexError):
                    continue

        # 保存为 TSV
        output_path = output_dir / "taxid_to_name.tsv"
        with open(output_path, 'w', encoding='utf-8', newline='') as f:
            writer = csv.writer(f, delimiter='\t')
            writer.writerow(['taxid', 'latin_name'])
            for taxid in sorted(taxid_to_name.keys()):
                writer.writerow([taxid, taxid_to_name[taxid]])

        logger.info("NCBI Taxonomy 物种名称已生成: %s (%d 物种)", output_path, len(taxid_to_name))

        # 记录版本元数据到 versions.json
        try:
            from allenricher.database.version import DatabaseVersionManager, RemoteVersionChecker
            _vm = DatabaseVersionManager(database_dir=str(self.root_dir))
            _checker = RemoteVersionChecker()
            _tax_info = _checker.check_head("https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz")
            if _tax_info:
                _vm.record_download(
                    source="taxonomy", local_version="cached",
                    local_path="basic/taxonomy",
                    remote_last_modified=_tax_info.get("last_modified"),
                )
        except Exception as _e:
            logger.warning("记录 Taxonomy 版本元数据失败: %s", _e)

        return output_path

    def _load_taxid_to_name_map(self, taxonomy_dir: Path) -> Dict[int, str]:
        """从 taxid_to_name.tsv 加载 taxid → latin_name 映射

        Args:
            taxonomy_dir: taxonomy 目录路径

        Returns:
            taxid → latin_name 映射字典
        """
        taxid_map_path = taxonomy_dir / "taxid_to_name.tsv"
        if not taxid_map_path.exists():
            return {}

        result: Dict[int, str] = {}
        with open(taxid_map_path, 'r', encoding='utf-8') as f:
            reader = csv.DictReader(f, delimiter='\t')
            for row in reader:
                try:
                    result[int(row['taxid'])] = row['latin_name']
                except (ValueError, KeyError):
                    continue
        return result

    # ============================
    # 物种注册表构建
    # ============================

    def _build_go_registry(
        self, gene2go_path: Path, output_dir: Path
    ) -> Path:
        """从 gene2go.gz 提取物种列表，生成 GO 物种注册表

        读取 gene2go.gz 中的所有唯一 taxid，再从 gene_info.gz 提取
        taxid → latin_name 映射，输出 go_species_registry.tsv。

        Args:
            gene2go_path: gene2go.gz 文件路径
            output_dir: 输出目录

        Returns:
            生成的 go_species_registry.tsv 文件路径
        """
        logger.info("构建 GO 物种注册表: gene2go_path=%s", gene2go_path)

        # ---- 1. 从 gene2go.gz 提取所有唯一 taxid 及统计 ----
        taxid_stats: Dict[int, Dict[str, int]] = {}  # taxid -> {gene_count, term_count}
        with gzip.open(gene2go_path, "rt", encoding="utf-8") as fh:
            for line in fh:
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                parts = line.split("\t")
                if len(parts) < 6:
                    continue
                try:
                    taxid = int(parts[0])
                    gene_id = parts[1]
                    go_term = parts[2]
                except (ValueError, IndexError):
                    continue
                if taxid not in taxid_stats:
                    taxid_stats[taxid] = {"gene_count": set(), "term_count": set()}
                taxid_stats[taxid]["gene_count"].add(gene_id)
                taxid_stats[taxid]["term_count"].add(go_term)

        # ---- 2. 从 NCBI Taxonomy 获取 taxid → latin_name 映射 ----
        # 优先使用 taxonomy 数据库中的 scientific name
        taxonomy_dir = self.basic_dir / "taxonomy"
        taxid_to_name = self._load_taxid_to_name_map(taxonomy_dir)

        # 如果 taxonomy 映射为空，先下载
        if not taxid_to_name:
            logger.info("NCBI Taxonomy 映射为空，尝试下载...")
            taxonomy_tsv = self._download_taxonomy_names(taxonomy_dir)
            if taxonomy_tsv:
                taxid_to_name = self._load_taxid_to_name_map(taxonomy_dir)

        # 如果仍然为空，回退到 gene_info（仅用于无法获取 taxonomy 的情况）
        # 注意：gene_info 不包含物种拉丁名，仅用于记录 taxid 存在性
        if not taxid_to_name:
            logger.warning("无法获取 NCBI Taxonomy，gene_info 不包含物种拉丁名信息")
            logger.warning("将使用空 latin_name 生成注册表，建议手动检查 taxonomy 数据")

        # ---- 3. 写入 go_species_registry.tsv ----
        output_dir.mkdir(parents=True, exist_ok=True)
        output_path = output_dir / "go_species_registry.tsv"

        with open(output_path, "w", encoding="utf-8", newline="") as fh:
            writer = csv.writer(fh, delimiter="\t")
            writer.writerow(["taxid", "latin_name", "source", "gene_count", "term_count"])
            for taxid in sorted(taxid_stats):
                latin_name = taxid_to_name.get(taxid, "")
                stats = taxid_stats[taxid]
                writer.writerow([
                    taxid,
                    latin_name,
                    "ncbi_gene2go",
                    len(stats["gene_count"]),
                    len(stats["term_count"]),
                ])

        logger.info(
            "GO 物种注册表已生成: %s (%d 物种)", output_path, len(taxid_stats)
        )
        return output_path

    def _download_goa_index(self, output_dir: Path) -> Path:
        """从 UniProt GOA FTP 获取物种索引

        请求 EBI GOA proteomes 目录页面，解析 HTML 提取所有 .goa 文件链接，
        生成 goa_species_index.tsv。

        Args:
            output_dir: 输出目录

        Returns:
            生成的 goa_species_index.tsv 文件路径
        """
        logger.info("下载 GOA 物种索引...")

        url = "https://ftp.ebi.ac.uk/pub/databases/GO/goa/proteomes/"
        headers = {
            "User-Agent": (
                "AllEnricher/2.0 (https://github.com/allenricher; "
                "data download pipeline)"
            ),
        }

        resp = requests.get(url, headers=headers, timeout=30)
        resp.raise_for_status()

        # 解析 HTML，提取 .goa 文件链接
        goa_entries: List[Dict[str, str]] = []
        taxid_list: List[int] = []  # 收集所有 taxid 用于后续查询
        # 匹配 href="xxx.goa" 或 href="xxx.goa.gz"
        for match in re.finditer(r'href="([^"]+\.goa(?:\.gz)?)"', resp.text):
            filename = match.group(1)
            # 去掉 .gz 后缀以统一处理
            base_name = filename.removesuffix(".gz")
            # 文件名格式: {taxid}.{species_name}.goa
            # 例如: 9606.Homo_sapiens.goa
            if not base_name.endswith(".goa"):
                continue
            stem = base_name[:-4]  # 去掉 .goa
            dot_pos = stem.find(".")
            if dot_pos < 0:
                continue
            taxid_str = stem[:dot_pos]
            species_part = stem[dot_pos + 1:]
            # taxid 必须是纯数字
            if not taxid_str.isdigit():
                continue
            taxid = int(taxid_str)
            taxid_list.append(taxid)
            # 暂时使用文件名中的名称，后续会从 NCBI Taxonomy 获取
            latin_name = species_part.replace("_", " ")
            goa_entries.append({
                "taxid": taxid,
                "latin_name": latin_name,
                "filename": filename,
            })

        # 从 NCBI Taxonomy 获取 latin_name
        taxonomy_dir = self.basic_dir / "taxonomy"
        taxid_to_name = self._load_taxid_to_name_map(taxonomy_dir)

        # 更新 latin_name：优先使用 NCBI Taxonomy 中的名称
        for entry in goa_entries:
            taxid = entry["taxid"]
            if taxid in taxid_to_name:
                entry["latin_name"] = taxid_to_name[taxid]

        # 写入 goa_species_index.tsv
        output_dir.mkdir(parents=True, exist_ok=True)
        output_path = output_dir / "goa_species_index.tsv"

        with open(output_path, "w", encoding="utf-8", newline="") as fh:
            writer = csv.writer(fh, delimiter="\t")
            writer.writerow(["taxid", "latin_name", "filename"])
            for entry in sorted(goa_entries, key=lambda e: e["taxid"]):
                writer.writerow([
                    entry["taxid"],
                    entry["latin_name"],
                    entry["filename"],
                ])

        logger.info(
            "GOA 物种索引已生成: %s (%d 物种)", output_path, len(goa_entries)
        )
        return output_path

    def _merge_go_registries(
        self,
        gene2go_registry: Path,
        goa_index: Path,
        output_dir: Path,
    ) -> Path:
        """合并 gene2go 和 GOA 两个 GO 注册表

        gene2go 优先级高于 GOA，source 标记为:
        - "ncbi_gene2go": 仅 gene2go 有
        - "uniprot_goa": 仅 GOA 有
        - "both": 两者都有

        Args:
            gene2go_registry: go_species_registry.tsv 路径
            goa_index: goa_species_index.tsv 路径
            output_dir: 输出目录

        Returns:
            生成的 go_species_registry.tsv 文件路径
        """
        logger.info("合并 GO 注册表...")

        # 读取 gene2go 注册表
        g2g_data: Dict[int, Dict[str, str]] = {}
        if gene2go_registry.exists():
            with open(gene2go_registry, "r", encoding="utf-8") as fh:
                reader = csv.DictReader(fh, delimiter="\t")
                for row in reader:
                    taxid = int(row["taxid"])
                    g2g_data[taxid] = {
                        "latin_name": row.get("latin_name", ""),
                        "source": row.get("source", "ncbi_gene2go"),
                        "gene_count": row.get("gene_count", ""),
                        "term_count": row.get("term_count", ""),
                    }

        # 读取 GOA 索引
        goa_data: Dict[int, Dict[str, str]] = {}
        if goa_index.exists():
            with open(goa_index, "r", encoding="utf-8") as fh:
                reader = csv.DictReader(fh, delimiter="\t")
                for row in reader:
                    taxid = int(row["taxid"])
                    goa_data[taxid] = {
                        "latin_name": row.get("latin_name", ""),
                        "filename": row.get("filename", ""),
                    }

        # 合并 - 优先使用更可靠的物种名来源
        all_taxids = sorted(set(g2g_data) | set(goa_data))
        merged: List[Dict[str, str]] = []
        for taxid in all_taxids:
            in_g2g = taxid in g2g_data
            in_goa = taxid in goa_data

            if in_g2g and in_goa:
                source = "both"
                # 优先使用 gene2go 的 latin_name（来自 NCBI Taxonomy）
                # 只有当 gene2go 名称为空时，才使用 GOA 的名称
                g2g_name = g2g_data[taxid]["latin_name"]
                goa_name = goa_data[taxid]["latin_name"]
                if g2g_name:
                    latin_name = g2g_name
                else:
                    latin_name = goa_name
                gene_count = g2g_data[taxid]["gene_count"]
                term_count = g2g_data[taxid]["term_count"]
            elif in_g2g:
                source = "ncbi_gene2go"
                latin_name = g2g_data[taxid]["latin_name"]
                gene_count = g2g_data[taxid]["gene_count"]
                term_count = g2g_data[taxid]["term_count"]
            else:
                source = "uniprot_goa"
                latin_name = goa_data[taxid]["latin_name"]
                gene_count = ""
                term_count = ""

            merged.append({
                "taxid": str(taxid),
                "latin_name": latin_name,
                "source": source,
                "gene_count": gene_count,
                "term_count": term_count,
            })

        # 写入合并后的注册表
        output_dir.mkdir(parents=True, exist_ok=True)
        output_path = output_dir / "go_species_registry.tsv"

        with open(output_path, "w", encoding="utf-8", newline="") as fh:
            writer = csv.DictWriter(
                fh, fieldnames=["taxid", "latin_name", "source", "gene_count", "term_count"],
                delimiter="\t",
            )
            writer.writeheader()
            writer.writerows(merged)

        logger.info(
            "GO 注册表合并完成: %s (%d 物种, gene2go=%d, goa=%d, both=%d)",
            output_path,
            len(merged),
            sum(1 for m in merged if m["source"] in ("ncbi_gene2go", "both")),
            sum(1 for m in merged if m["source"] in ("uniprot_goa", "both")),
            sum(1 for m in merged if m["source"] == "both"),
        )
        return output_path

    def _build_kegg_registry(self, output_dir: Path) -> Path:
        """调用 KEGG API 构建物种注册表

        请求 KEGG list/organism 接口，解析返回的 TSV 数据，
        生成 kegg_species_registry.tsv。

        Args:
            output_dir: 输出目录

        Returns:
            生成的 kegg_species_registry.tsv 文件路径
        """
        logger.info("构建 KEGG 物种注册表...")

        url = "https://rest.kegg.jp/list/organism"
        headers = {
            "User-Agent": (
                "AllEnricher/2.0 (https://github.com/allenricher; "
                "data download pipeline)"
            ),
        }

        text = self._api_get_with_retry(
            url, headers=headers, timeout=30, max_retries=3
        )

        # 解析 TSV: kegg_code\tlatin_name\ttaxid\tdefinition\tgene_count\t...
        entries: List[Dict[str, str]] = []
        for line in text.strip().split("\n"):
            if not line.strip():
                continue
            parts = line.split("\t")
            if len(parts) < 5:
                continue
            kegg_code = parts[0].strip()
            latin_name = parts[1].strip()
            taxid_str = parts[2].strip()
            gene_count_str = parts[4].strip()

            if not taxid_str.isdigit():
                continue

            entries.append({
                "taxid": taxid_str,
                "latin_name": latin_name,
                "kegg_code": kegg_code,
                "kegg_code_source": "kegg",
                "gene_count": gene_count_str,
            })

        # 写入 kegg_species_registry.tsv
        output_dir.mkdir(parents=True, exist_ok=True)
        output_path = output_dir / "kegg_species_registry.tsv"

        with open(output_path, "w", encoding="utf-8", newline="") as fh:
            writer = csv.DictWriter(
                fh,
                fieldnames=["taxid", "latin_name", "kegg_code", "kegg_code_source", "gene_count"],
                delimiter="\t",
            )
            writer.writeheader()
            writer.writerows(entries)

        logger.info(
            "KEGG 物种注册表已生成: %s (%d 物种)", output_path, len(entries)
        )
        return output_path

    def _build_reactome_registry(
        self, ncbi2reactome_path: Path, output_dir: Path
    ) -> Path:
        """从 NCBI2Reactome 文件提取物种列表，生成 Reactome 注册表

        解析 NCBI2Reactome_All_Levels.txt.gz，从 pathway_id 中提取
        Reactome 物种代码，通过内置映射表关联到 taxid。

        Args:
            ncbi2reactome_path: NCBI2Reactome_All_Levels.txt.gz 路径
            output_dir: 输出目录

        Returns:
            生成的 reactome_species_registry.tsv 文件路径
        """
        logger.info("构建 Reactome 物种注册表...")

        # Reactome 物种代码 → taxid 内置映射
        REACTOME_CODE_TO_TAXID: Dict[str, int] = {
            "HSA": 9606, "MMU": 10090, "RNO": 10116, "CEL": 6239,
            "DME": 7227, "SCE": 4932, "ATH": 3702, "DDI": 44689,
            "GGA": 9031, "SSC": 9823, "BTA": 9913, "XTR": 8364,
            "CFA": 9615, "DRE": 7955, "PFA": 5833, "SPO": 4896,
            "MTU": 1772,
        }

        # 反向映射: taxid → reactome_code
        taxid_to_code: Dict[int, str] = {}

        # 从文件中提取 pathway_id 中的物种代码
        pathway_pattern = re.compile(r"R-([A-Z]{3})-\d+")
        with gzip.open(ncbi2reactome_path, "rt", encoding="utf-8") as fh:
            for line in fh:
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                parts = line.split("\t")
                if len(parts) < 3:
                    continue
                pathway_id = parts[2]
                m = pathway_pattern.match(pathway_id)
                if m:
                    code = m.group(1)
                    if code in REACTOME_CODE_TO_TAXID:
                        taxid = REACTOME_CODE_TO_TAXID[code]
                        taxid_to_code[taxid] = code

        # 写入 reactome_species_registry.tsv
        output_dir.mkdir(parents=True, exist_ok=True)
        output_path = output_dir / "reactome_species_registry.tsv"

        with open(output_path, "w", encoding="utf-8", newline="") as fh:
            writer = csv.writer(fh, delimiter="\t")
            writer.writerow(["taxid", "reactome_code"])
            for taxid in sorted(taxid_to_code):
                writer.writerow([taxid, taxid_to_code[taxid]])

        logger.info(
            "Reactome 物种注册表已生成: %s (%d 物种)",
            output_path, len(taxid_to_code),
        )
        return output_path

    def _build_do_registry(self, output_dir: Path) -> Path:
        """生成 DO（Disease Ontology）物种注册表

        DO 目前仅支持人类 (taxid=9606)。

        Args:
            output_dir: 输出目录

        Returns:
            生成的 do_species_registry.tsv 文件路径
        """
        logger.info("构建 DO 物种注册表...")

        output_dir.mkdir(parents=True, exist_ok=True)
        output_path = output_dir / "do_species_registry.tsv"

        with open(output_path, "w", encoding="utf-8", newline="") as fh:
            writer = csv.writer(fh, delimiter="\t")
            writer.writerow(["taxid", "latin_name"])
            writer.writerow([9606, "Homo sapiens"])

        logger.info("DO 物种注册表已生成: %s", output_path)
        return output_path

    def _merge_all_registries(
        self,
        go_registry: Path,
        kegg_registry: Path,
        reactome_registry: Path,
        do_registry: Path,
        output_path: Path,
    ) -> Path:
        """合并所有专用注册表为统一 supported_species.tsv

        读取 GO、KEGG、Reactome、DO 四个注册表，按 taxid 合并。
        为没有 kegg_code 的物种自动生成缩写。

        Args:
            go_registry: go_species_registry.tsv 路径
            kegg_registry: kegg_species_registry.tsv 路径
            reactome_registry: reactome_species_registry.tsv 路径
            do_registry: do_species_registry.tsv 路径
            output_path: 输出文件路径（supported_species.tsv）

        Returns:
            生成的 supported_species.tsv 文件路径
        """
        logger.info("合并所有物种注册表...")

        # ---- 读取 GO 注册表 ----
        go_data: Dict[int, Dict[str, str]] = {}
        if go_registry.exists():
            with open(go_registry, "r", encoding="utf-8") as fh:
                reader = csv.DictReader(fh, delimiter="\t")
                for row in reader:
                    taxid = int(row["taxid"])
                    go_data[taxid] = {
                        "latin_name": row.get("latin_name", ""),
                        "source": row.get("source", ""),
                        "gene_count": row.get("gene_count", ""),
                        "term_count": row.get("term_count", ""),
                    }

        # ---- 读取 KEGG 注册表 ----
        kegg_data: Dict[int, Dict[str, str]] = {}
        if kegg_registry and kegg_registry.exists():
            with open(kegg_registry, "r", encoding="utf-8") as fh:
                reader = csv.DictReader(fh, delimiter="\t")
                for row in reader:
                    taxid = int(row["taxid"])
                    kegg_data[taxid] = {
                        "latin_name": row.get("latin_name", ""),
                        "kegg_code": row.get("kegg_code", ""),
                        "kegg_code_source": row.get("kegg_code_source", "kegg"),
                        "gene_count": row.get("gene_count", ""),
                    }

        # ---- 读取 Reactome 注册表 ----
        reactome_data: Dict[int, Dict[str, str]] = {}
        if reactome_registry and reactome_registry.exists():
            with open(reactome_registry, "r", encoding="utf-8") as fh:
                reader = csv.DictReader(fh, delimiter="\t")
                for row in reader:
                    taxid = int(row["taxid"])
                    reactome_data[taxid] = {
                        "reactome_code": row.get("reactome_code", ""),
                    }

        # ---- 读取 DO 注册表 ----
        do_taxids: set = set()
        if do_registry and do_registry.exists():
            with open(do_registry, "r", encoding="utf-8") as fh:
                reader = csv.DictReader(fh, delimiter="\t")
                for row in reader:
                    try:
                        do_taxids.add(int(row["taxid"]))
                    except (ValueError, KeyError):
                        continue

        # ---- 合并 ----
        all_taxids = sorted(
            set(go_data) | set(kegg_data) | set(reactome_data) | do_taxids
        )

        entries: List[SpeciesEntry] = []
        for taxid in all_taxids:
            # 确定 latin_name: GO > KEGG > 空
            latin_name = ""
            if taxid in go_data and go_data[taxid]["latin_name"]:
                latin_name = go_data[taxid]["latin_name"]
            elif taxid in kegg_data and kegg_data[taxid]["latin_name"]:
                latin_name = kegg_data[taxid]["latin_name"]

            entry = SpeciesEntry(taxid=taxid, latin_name=latin_name)

            # GO
            if taxid in go_data:
                entry.has_go = True
                entry.go_source = go_data[taxid].get("source")
                gc = go_data[taxid].get("gene_count", "")
                tc = go_data[taxid].get("term_count", "")
                entry.go_gene_count = int(gc) if gc and gc.isdigit() else None
                entry.go_term_count = int(tc) if tc and tc.isdigit() else None

            # KEGG
            if taxid in kegg_data:
                entry.has_kegg = True
                entry.kegg_code = kegg_data[taxid].get("kegg_code", "")
                entry.kegg_code_source = kegg_data[taxid].get(
                    "kegg_code_source", "kegg"
                )
                gc = kegg_data[taxid].get("gene_count", "")
                entry.kegg_gene_count = int(gc) if gc and gc.isdigit() else None
            elif latin_name:
                # 自动生成 KEGG 缩写
                entry.kegg_code = SpeciesRegistry.generate_kegg_abbreviation(
                    latin_name
                )
                entry.kegg_code_source = "auto"
                entry.has_kegg = True  # 自动生成的 KEGG 代码也算作有 KEGG 支持

            # Reactome
            if taxid in reactome_data:
                entry.has_reactome = True
                entry.reactome_code = reactome_data[taxid].get("reactome_code", "")

            # DO
            if taxid in do_taxids:
                entry.has_do = True

            entries.append(entry)

        # ---- 写入 supported_species.tsv ----
        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)

        registry = SpeciesRegistry(registry_path=output_path)
        for entry in entries:
            registry.add_entry(entry)
        registry.save()

        logger.info(
            "统一物种注册表已生成: %s (%d 物种)", output_path, len(entries)
        )
        return output_path

    def _report_download_summary(self, registry_path: Path) -> None:
        """打印格式化的下载统计摘要到 stderr

        加载 supported_species.tsv，统计各数据库覆盖情况并输出。

        Args:
            registry_path: supported_species.tsv 文件路径
        """
        registry = SpeciesRegistry(registry_path=registry_path)
        registry.load()

        if not registry.entries:
            print("注册表为空，无统计数据。", file=sys.stderr)
            return

        summary = registry.get_summary()
        total = summary["total_species"]

        go_info = summary["go"]
        kegg_info = summary["kegg"]
        reactome_info = summary["reactome"]
        do_info = summary["do"]

        sep = "=" * 60
        print(f"\n{sep}", file=sys.stderr)
        print("  数据下载统计摘要", file=sys.stderr)
        print(sep, file=sys.stderr)
        print(f"  总物种数:            {total}", file=sys.stderr)
        print(f"{'─' * 60}", file=sys.stderr)
        print(f"  GO 支持物种数:       {go_info['count']}", file=sys.stderr)
        print(f"    - 有基因数统计:    {go_info['with_gene_count']}", file=sys.stderr)
        print(f"    - 有术语数统计:    {go_info['with_term_count']}", file=sys.stderr)
        print(f"{'─' * 60}", file=sys.stderr)
        print(f"  KEGG 支持物种数:     {kegg_info['count']}", file=sys.stderr)
        print(f"    - 有基因数统计:    {kegg_info['with_gene_count']}", file=sys.stderr)
        print(f"    - 有通路数统计:    {kegg_info['with_pathway_count']}", file=sys.stderr)
        print(f"{'─' * 60}", file=sys.stderr)
        print(f"  Reactome 支持物种数: {reactome_info['count']}", file=sys.stderr)
        print(f"    - 有基因数统计:    {reactome_info['with_gene_count']}", file=sys.stderr)
        print(f"    - 有通路数统计:    {reactome_info['with_pathway_count']}", file=sys.stderr)
        print(f"{'─' * 60}", file=sys.stderr)
        print(f"  DO 支持物种数:       {do_info['count']}", file=sys.stderr)
        print(f"    - 有基因数统计:    {do_info['with_gene_count']}", file=sys.stderr)
        print(f"    - 有术语数统计:    {do_info['with_term_count']}", file=sys.stderr)
        print(sep, file=sys.stderr)

    # ============================
    # 内部辅助方法
    # ============================

    def _api_get_with_retry(
        self,
        url: str,
        headers: Optional[Dict[str, str]] = None,
        timeout: int = 30,
        max_retries: int = 3,
    ) -> str:
        """带重试的 HTTP GET 请求

        参考 kegg_fetcher 的模式，使用 requests 库发起请求，
        失败时按指数退避重试。

        Args:
            url: 请求 URL
            headers: 请求头字典
            timeout: 超时秒数
            max_retries: 最大重试次数

        Returns:
            响应文本内容

        Raises:
            RuntimeError: 所有重试均失败
        """
        if headers is None:
            headers = {
                "User-Agent": (
                    "AllEnricher/2.0 (https://github.com/allenricher; "
                    "data download pipeline)"
                ),
            }

        last_error: Optional[Exception] = None
        for attempt in range(1, max_retries + 1):
            try:
                resp = requests.get(url, headers=headers, timeout=timeout)
                resp.raise_for_status()
                return resp.text
            except requests.RequestException as exc:
                last_error = exc
                logger.warning(
                    "请求失败 (第 %d/%d 次): %s - %s",
                    attempt, max_retries, url, exc,
                )
                if attempt < max_retries:
                    wait = 2 ** attempt
                    logger.info("等待 %d 秒后重试...", wait)
                    time.sleep(wait)

        raise RuntimeError(
            f"API 请求失败（已重试 {max_retries} 次）: {url}"
        ) from last_error
