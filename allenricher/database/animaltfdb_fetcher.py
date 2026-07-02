"""
AnimalTFDB 4.0 数据下载器

从 guolab.wchscu.cn/AnimalTFDB4/ 下载动物转录因子注释数据。

AnimalTFDB 4.0 特性：
- 183 个动物物种
- 274,633 个 TF，73 个家族
- TF 列表、直系同源映射、蛋白序列

数据源：
- AnimalTFDB 4.0: https://guolab.wchscu.cn/AnimalTFDB4/
"""

from pathlib import Path
from typing import Dict, List, Optional, Set
import requests
import logging

logger = logging.getLogger(__name__)

# AnimalTFDB 下载基础 URL
ANIMALTFDB_BASE_URL = (
    "https://guolab.wchscu.cn/AnimalTFDB4_static/download"
)

# 主要模式生物和经济动物（优先支持）
ANIMALTFDB_PRIORITY_SPECIES = [
    "Homo_sapiens", "Mus_musculus", "Rattus_norvegicus",
    "Danio_rerio", "Drosophila_melanogaster", "Caenorhabditis_elegans",
    "Bos_taurus", "Sus_scrofa", "Ovis_aries", "Capra_hircus",
    "Gallus_gallus", "Canis_lupus_familiaris", "Equus_caballus",
    "Felis_catus", "Macaca_mulatta", "Gorilla_gorilla_gorilla",
    "Pan_troglodytes", "Oryctolagus_cuniculus", "Xenopus_tropicalis",
]


class AnimalTFDBFetcher:
    """AnimalTFDB 4.0 数据下载器

    下载动物转录因子注释和直系同源映射数据。

    Usage::

        fetcher = AnimalTFDBFetcher(basic_dir='./database/basic')
        # 下载指定物种的TF列表和同源映射
        fetcher.download_species_data("Bos_taurus")
    """

    REQUEST_TIMEOUT = 120

    def __init__(self, basic_dir: str):
        """
        Args:
            basic_dir: 基础缓存目录
        """
        self.basic_dir = Path(basic_dir)
        self.basic_dir.mkdir(parents=True, exist_ok=True)

    def _get_cache_dir(self) -> Path:
        """获取缓存目录"""
        cache_dir = self.basic_dir / "animaltfdb" / "AnimalTFDBv4.0"
        cache_dir.mkdir(parents=True, exist_ok=True)
        return cache_dir

    def _download_url(self, url: str, local_path: Path, overwrite: bool = False) -> Path:
        """通用下载方法"""
        if local_path.exists() and not overwrite:
            logger.info(f"已缓存，跳过: {local_path.name}")
            return local_path

        logger.info(f"下载: {url}")

        try:
            resp = requests.get(
                url, timeout=self.REQUEST_TIMEOUT,
                stream=True, headers={"User-Agent": "Mozilla/5.0"},
            )
            resp.raise_for_status()
        except requests.RequestException as e:
            raise RuntimeError(f"下载失败: {e}") from e

        with open(local_path, "wb") as f:
            for chunk in resp.iter_content(chunk_size=8192):
                if chunk:
                    f.write(chunk)

        logger.info(f"已保存: {local_path} ({local_path.stat().st_size / 1024:.1f} KB)")
        return local_path

    def download_tf_list(self, species: str, overwrite: bool = False) -> Path:
        """下载指定物种的 TF 列表

        Args:
            species: 物种拉丁名（下划线格式，如 Bos_taurus）
            overwrite: 是否覆盖

        Returns:
            TF 列表文件路径
        """
        cache_dir = self._get_cache_dir()
        url = f"{ANIMALTFDB_BASE_URL}/TF_list_final/{species}_TF"
        local_path = cache_dir / f"{species}_TF"
        return self._download_url(url, local_path, overwrite)

    def download_ortholog_to_human(self, species: str, overwrite: bool = False) -> Path:
        """下载指定物种到人类的直系同源映射

        Args:
            species: 物种拉丁名（下划线格式）
            overwrite: 是否覆盖

        Returns:
            同源映射文件路径
        """
        cache_dir = self._get_cache_dir()
        url = f"{ANIMALTFDB_BASE_URL}/ortholog_to_human_download/{species}_ortholog_to_human"
        local_path = cache_dir / f"{species}_ortholog_to_human"
        return self._download_url(url, local_path, overwrite)

    def download_species_data(self, species: str, overwrite: bool = False) -> Dict[str, Path]:
        """下载指定物种的全部数据（TF列表 + 同源映射）

        Args:
            species: 物种拉丁名（下划线格式）
            overwrite: 是否覆盖

        Returns:
            {'tf_list': Path, 'ortholog': Path}
        """
        results = {}
        try:
            results['tf_list'] = self.download_tf_list(species, overwrite)
        except Exception as e:
            logger.error(f"下载 {species} TF列表失败: {e}")

        try:
            results['ortholog'] = self.download_ortholog_to_human(species, overwrite)
        except Exception as e:
            logger.error(f"下载 {species} 同源映射失败: {e}")

        return results

    def download_htftarget(self, overwrite: bool = False) -> Path:
        """下载 hTFtarget 人类TF-target关系数据

        Args:
            overwrite: 是否覆盖

        Returns:
            hTFtarget 文件路径
        """
        from .htftarget_fetcher import HTFtargetFetcher
        fetcher = HTFtargetFetcher(str(self.basic_dir))
        return fetcher.download(overwrite)

    @staticmethod
    def get_priority_species() -> List[str]:
        """获取优先支持的物种列表"""
        return ANIMALTFDB_PRIORITY_SPECIES.copy()

    @staticmethod
    def get_info() -> Dict[str, str]:
        """获取数据库基本信息"""
        return {
            "name": "AnimalTFDB",
            "version": "4.0",
            "url": "https://guolab.wchscu.cn/AnimalTFDB4/",
            "species_count": 183,
            "tf_count": 274633,
            "description": "183个动物物种的TF分类注释和直系同源映射",
        }
