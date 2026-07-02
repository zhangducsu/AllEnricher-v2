"""ChEA3 数据下载器和 API 客户端

ChEA3 (ChIP-X Enrichment Analysis 3) 提供转录因子富集分析功能，
包含 6 个 GMT 文件库和在线富集 API。

GMT 库：
  - ENCODE: ChIP-seq 数据 (ENCODE 项目)
  - ReMap: ChIP-seq 数据 (ReMap 项目)
  - LiteratureChIP: 文献挖掘的 ChIP-seq 数据
  - GTExCoexpression: GTEx 共表达数据
  - ARCHS4Coexpression: ARCHS4 共表达数据
  - EnrichrQueries: Enrichr 查询数据

API: POST https://maayanlab.cloud/chea3/api/enrich/
  请求体: {"query_name": "...", "gene_set": ["GENE1", "GENE2", ...]}
"""

from __future__ import annotations

import logging
import time
from pathlib import Path
from typing import Any, Dict, List

import requests

logger = logging.getLogger(__name__)


class ChEA3Fetcher:
    """ChEA3 数据下载器和 API 客户端

    Usage::

        fetcher = ChEA3Fetcher(basic_dir='./database/basic')
        paths = fetcher.download_all_gmt_libraries()
        result = fetcher.enrich_api(['TP53', 'BRCA1'], query_name='test')
    """

    # 物种名称 → 代码映射
    CHEA3_SPECIES_MAP: Dict[str, str] = {
        "Homo sapiens": "hsa",
        "Mus musculus": "mmu",
        "Rattus norvegicus": "rno",
    }

    # 6 个 GMT 库的名称 → URL 映射
    CHEA3_GMT_LIBS: Dict[str, str] = {
        "ENCODE": "https://maayanlab.cloud/chea3/assets/tflibs/ENCODE_tf.gmt",
        "ReMap": "https://maayanlab.cloud/chea3/assets/tflibs/ReMap_tf.gmt",
        "LiteratureChIP": "https://maayanlab.cloud/chea3/assets/tflibs/LiteratureChIP_tf.gmt",
        "GTExCoexpression": "https://maayanlab.cloud/chea3/assets/tflibs/GTEx_tf.gmt",
        "ARCHS4Coexpression": "https://maayanlab.cloud/chea3/assets/tflibs/ARCHS4_tf.gmt",
        "EnrichrQueries": "https://maayanlab.cloud/chea3/assets/tflibs/EnrichrQueries_tf.gmt",
    }

    # API 端点
    API_URL = "https://maayanlab.cloud/chea3/api/enrich/"

    # 请求超时（秒）
    TIMEOUT = 120

    # User-Agent
    UA = "Mozilla/5.0 (Windows NT 10.0; Win64; x64)"

    def __init__(self, basic_dir: str):
        """
        Args:
            basic_dir: 数据库基础目录（如 ./database/basic）
        """
        self.basic_dir = Path(basic_dir)

    def _get_cache_dir(self) -> Path:
        """获取 ChEA3 缓存目录

        Returns:
            basic_dir/chea3/ChEA3v2024/
        """
        cache_dir = self.basic_dir / "chea3" / "ChEA3v2024"
        cache_dir.mkdir(parents=True, exist_ok=True)
        return cache_dir

    def download_gmt_library(self, lib_name: str, overwrite: bool = False) -> Path:
        """下载单个 GMT 库文件

        Args:
            lib_name: 库名称（如 'ENCODE'），必须是 CHEA3_GMT_LIBS 中的键
            overwrite: 是否覆盖已有缓存

        Returns:
            本地缓存文件路径

        Raises:
            ValueError: lib_name 不在支持的库列表中
            requests.RequestException: 下载失败
        """
        if lib_name not in self.CHEA3_GMT_LIBS:
            raise ValueError(
                f"未知库 '{lib_name}'，支持的库: {list(self.CHEA3_GMT_LIBS.keys())}"
            )

        url = self.CHEA3_GMT_LIBS[lib_name]
        cache_dir = self._get_cache_dir()
        local_file = cache_dir / f"{lib_name}_tf.gmt"

        # 缓存检查
        if local_file.exists() and not overwrite:
            logger.info("ChEA3 GMT 库已缓存，跳过下载: %s", local_file)
            return local_file

        logger.info("下载 ChEA3 GMT 库: %s -> %s", url, local_file)

        resp = requests.get(
            url,
            headers={"User-Agent": self.UA},
            timeout=self.TIMEOUT,
        )
        resp.raise_for_status()

        local_file.write_bytes(resp.content)
        logger.info(
            "ChEA3 GMT 库下载完成: %s (%d bytes)",
            local_file, len(resp.content),
        )

        return local_file

    def download_all_gmt_libraries(
        self, overwrite: bool = False
    ) -> Dict[str, Path]:
        """下载所有 6 个 GMT 库文件

        Args:
            overwrite: 是否覆盖已有缓存

        Returns:
            {库名称: 本地文件路径} 字典
        """
        results: Dict[str, Path] = {}
        for lib_name in self.CHEA3_GMT_LIBS:
            results[lib_name] = self.download_gmt_library(lib_name, overwrite=overwrite)
        return results

    def enrich_api(
        self,
        gene_set: List[str],
        query_name: str | None = None,
        max_retries: int = 3,
        backoff_factor: float = 1.0,
    ) -> Dict[str, Any]:
        """调用 ChEA3 富集分析 API

        Args:
            gene_set: 基因符号列表（如 ['TP53', 'BRCA1', 'MYC']）
            query_name: 查询名称（可选，默认使用前 5 个基因拼接）
            max_retries: 最大重试次数（默认 3）
            backoff_factor: 指数退避因子（默认 1.0）

        Returns:
            API 返回的 JSON 响应字典

        Raises:
            requests.RequestException: API 请求失败
            ValueError: gene_set 为空或响应格式无效
        """
        if not gene_set:
            raise ValueError("gene_set 不能为空")

        if query_name is None:
            query_name = ", ".join(gene_set[:5])

        payload = {
            "query_name": query_name,
            "gene_set": gene_set,
        }

        logger.info("调用 ChEA3 API: %d 个基因", len(gene_set))

        # 指数退避重试机制
        last_exception = None
        for attempt in range(max_retries):
            try:
                resp = requests.post(
                    self.API_URL,
                    json=payload,
                    headers={"User-Agent": self.UA},
                    timeout=self.TIMEOUT,
                )
                resp.raise_for_status()

                result = resp.json()

                # 响应格式验证
                self._validate_api_response(result)

                logger.info("ChEA3 API 返回结果，包含 %d 个库", len(result))
                return result

            except (requests.RequestException, ValueError) as e:
                last_exception = e
                logger.warning(
                    "ChEA3 API 请求失败 (尝试 %d/%d): %s",
                    attempt + 1, max_retries, str(e)
                )
                if attempt < max_retries - 1:
                    sleep_time = backoff_factor * (2 ** attempt)
                    logger.info("等待 %.1f 秒后重试...", sleep_time)
                    time.sleep(sleep_time)

        # 所有重试都失败
        logger.error("ChEA3 API 请求在 %d 次尝试后均失败", max_retries)
        raise last_exception

    def _validate_api_response(self, result: Any) -> None:
        """验证 ChEA3 API 响应格式

        Args:
            result: API 返回的 JSON 数据

        Raises:
            ValueError: 响应格式无效
        """
        # 验证是否为字典
        if not isinstance(result, dict):
            raise ValueError(
                f"API 响应格式错误: 期望 dict，实际为 {type(result).__name__}"
            )

        # 验证是否包含预期的库键
        expected_libs = set(self.CHEA3_GMT_LIBS.keys())
        actual_libs = set(result.keys())

        # 允许部分库返回，但至少需要一个库
        if not actual_libs:
            raise ValueError("API 响应为空，未返回任何库数据")

        # 验证每个库的数据格式
        for lib_name, lib_data in result.items():
            if not isinstance(lib_data, list):
                raise ValueError(
                    f"库 '{lib_name}' 的数据格式错误: 期望 list，"
                    f"实际为 {type(lib_data).__name__}"
                )
            # 验证列表中的每个条目是否为字典
            for i, entry in enumerate(lib_data):
                if not isinstance(entry, dict):
                    raise ValueError(
                        f"库 '{lib_name}' 的第 {i} 个条目格式错误: "
                        f"期望 dict，实际为 {type(entry).__name__}"
                    )

    @staticmethod
    def get_supported_species() -> List[str]:
        """获取支持的物种列表

        Returns:
            物种拉丁名列表（如 ['Homo sapiens', 'Mus musculus', ...]）
        """
        return list(ChEA3Fetcher.CHEA3_SPECIES_MAP.keys())

    @staticmethod
    def get_library_names() -> List[str]:
        """获取支持的 GMT 库名称列表

        Returns:
            库名称列表
        """
        return list(ChEA3Fetcher.CHEA3_GMT_LIBS.keys())
