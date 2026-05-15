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

import gzip
import shutil
from pathlib import Path
from datetime import datetime
from typing import Dict, List, Optional

from .download_manager import DownloadManager
from .mirrors import get_mirrors, JENSEN_SOURCES


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
        return str(re_dir)

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
