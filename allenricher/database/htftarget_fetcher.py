"""
hTFtarget 数据下载器

从 guolab.wchscu.cn 下载人类转录因子-靶基因关系数据。

hTFtarget 特性：
- 基于 ENCODE/SRA 的 ChIP-Seq 数据
- 659 个人类 TF 的 ~134 万条 TF-target 关系
- 30% 阈值过滤（至少 3 个数据集）
- 包含组织来源信息

数据源：
- hTFtarget: https://guolab.wchscu.cn/hTFtarget/
"""

from pathlib import Path
from typing import Dict, List, Optional
import requests
import logging

logger = logging.getLogger(__name__)

# hTFtarget 下载链接
HTFTARGET_DOWNLOAD_URL = (
    "https://guolab.wchscu.cn/static/hTFtarget/file_download/tf-target-infomation.txt"
)


class HTFtargetFetcher:
    """hTFtarget 数据下载器

    下载人类转录因子-靶基因关系数据。

    Usage::

        fetcher = HTFtargetFetcher(basic_dir='./database/basic')
        fetcher.download()
    """

    REQUEST_TIMEOUT = 300  # 56MB 文件需要较长时间

    def __init__(self, basic_dir: str):
        """
        Args:
            basic_dir: 基础缓存目录
        """
        self.basic_dir = Path(basic_dir)
        self.basic_dir.mkdir(parents=True, exist_ok=True)

    def _get_cache_dir(self) -> Path:
        """获取缓存目录"""
        cache_dir = self.basic_dir / "htftarget"
        cache_dir.mkdir(parents=True, exist_ok=True)
        return cache_dir

    def download(self, overwrite: bool = False) -> Path:
        """下载 hTFtarget TF-target 关系文件

        Args:
            overwrite: 是否覆盖已存在的文件

        Returns:
            下载的 TSV 文件路径
        """
        cache_dir = self._get_cache_dir()
        local_path = cache_dir / "tf-target-information.txt"

        if local_path.exists() and not overwrite:
            logger.info(f"hTFtarget 已缓存，跳过: {local_path}")
            return local_path

        logger.info(f"下载 hTFtarget: {HTFTARGET_DOWNLOAD_URL}")

        try:
            resp = requests.get(
                HTFTARGET_DOWNLOAD_URL,
                timeout=self.REQUEST_TIMEOUT,
                stream=True,
                headers={"User-Agent": "Mozilla/5.0"},
            )
            resp.raise_for_status()
        except requests.RequestException as e:
            raise RuntimeError(f"hTFtarget 下载失败: {e}") from e

        # 流式写入
        with open(local_path, "wb") as f:
            for chunk in resp.iter_content(chunk_size=8192):
                if chunk:
                    f.write(chunk)

        size_mb = local_path.stat().st_size / 1024 / 1024
        logger.info(f"hTFtarget 已保存: {local_path} ({size_mb:.1f} MB)")
        return local_path

    @staticmethod
    def get_info() -> Dict[str, str]:
        """获取数据库基本信息"""
        return {
            "name": "hTFtarget",
            "version": "2020",
            "url": "https://guolab.wchscu.cn/hTFtarget/",
            "species": "Homo sapiens",
            "description": "人类TF-target关系，基于ChIP-Seq",
        }
