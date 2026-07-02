"""TRRUST v2 转录调控关系数据下载器

通过 GRNpedia 获取 TRRUST v2 数据库中转录因子-靶基因调控关系。

数据源：
  - Human: https://www.grnpedia.org/trrust/data/trrust_rawdata.human.tsv
  - Mouse: https://www.grnpedia.org/trrust/data/trrust_rawdata.mouse.tsv

TSV 格式：TF基因名\\ttarget基因名\\tmode_of_regulation(Activation/Repression)

TRRUST 仅支持 Human 和 Mouse 两种物种。
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Dict, List, Optional

import requests

logger = logging.getLogger(__name__)

# 物种拉丁名 → TRRUST 物种代码
TRRUST_SPECIES_MAP: Dict[str, str] = {
    "Homo sapiens": "hsa",
    "Mus musculus": "mmu",
}

# 物种代码 → 下载 URL
TRRUST_DOWNLOAD_URLS: Dict[str, str] = {
    "hsa": "https://www.grnpedia.org/trrust/data/trrust_rawdata.human.tsv",
    "mmu": "https://www.grnpedia.org/trrust/data/trrust_rawdata.mouse.tsv",
}

# 请求超时（秒）
_TIMEOUT = 60

# User-Agent
_UA = "Mozilla/5.0 (Windows NT 10.0; Win64; x64)"


class TRRUSTFetcher:
    """TRRUST v2 数据下载器

    Usage::

        fetcher = TRRUSTFetcher(basic_dir='./database/basic')
        path = fetcher.download_species('Homo sapiens')
        all_paths = fetcher.download_all()
    """

    def __init__(self, basic_dir: str):
        """
        Args:
            basic_dir: 基础数据目录（缓存根目录）
        """
        self._basic_dir = Path(basic_dir)

    def _get_cache_dir(self) -> Path:
        """获取 TRRUST v2 缓存目录

        Returns:
            basic_dir/trrust/TRRUSTv2/
        """
        cache_dir = self._basic_dir / "trrust" / "TRRUSTv2"
        cache_dir.mkdir(parents=True, exist_ok=True)
        return cache_dir

    def download_species(self, species: str, overwrite: bool = False) -> Path:
        """下载指定物种的 TRRUST 数据

        Args:
            species: 物种拉丁名（如 'Homo sapiens'）
            overwrite: 是否覆盖已有缓存

        Returns:
            下载文件的本地路径

        Raises:
            ValueError: 不支持的物种
            requests.RequestException: 下载失败
        """
        species_code = self.get_species_code(species)
        if species_code is None:
            raise ValueError(
                f"TRRUST 不支持的物种: '{species}'，"
                f"支持的物种: {self.get_supported_species()}"
            )

        url = TRRUST_DOWNLOAD_URLS[species_code]
        cache_dir = self._get_cache_dir()
        local_file = cache_dir / f"trrust_rawdata.{species_code}.tsv"

        if local_file.exists() and not overwrite:
            logger.info("TRRUST 数据已缓存，跳过下载: %s", local_file)
            return local_file

        logger.info("下载 TRRUST 数据: %s -> %s", url, local_file)

        resp = requests.get(
            url,
            headers={"User-Agent": _UA},
            timeout=_TIMEOUT,
        )
        resp.raise_for_status()

        local_file.write_bytes(resp.content)
        logger.info(
            "TRRUST 数据下载完成: %s (%d bytes)",
            local_file,
            len(resp.content),
        )

        return local_file

    def download_all(self, overwrite: bool = False) -> Dict[str, Path]:
        """下载所有支持物种的 TRRUST 数据

        Args:
            overwrite: 是否覆盖已有缓存

        Returns:
            {物种拉丁名: 文件路径, ...}
        """
        results: Dict[str, Path] = {}
        for species in self.get_supported_species():
            results[species] = self.download_species(species, overwrite=overwrite)
        return results

    @staticmethod
    def get_species_code(latin_name: str) -> Optional[str]:
        """根据物种拉丁名获取 TRRUST 物种代码

        Args:
            latin_name: 物种拉丁名（如 'Homo sapiens'）

        Returns:
            物种代码（如 'hsa'），不支持则返回 None
        """
        return TRRUST_SPECIES_MAP.get(latin_name)

    @staticmethod
    def get_latin_name(species_code: str) -> Optional[str]:
        """根据 TRRUST 物种代码获取物种拉丁名

        Args:
            species_code: 物种代码（如 'hsa'）

        Returns:
            物种拉丁名（如 'Homo sapiens'），不支持则返回 None
        """
        code_to_name = {v: k for k, v in TRRUST_SPECIES_MAP.items()}
        return code_to_name.get(species_code)

    @staticmethod
    def get_supported_species() -> List[str]:
        """获取 TRRUST 支持的物种列表

        Returns:
            物种拉丁名列表
        """
        return list(TRRUST_SPECIES_MAP.keys())
