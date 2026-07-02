"""WikiPathways 数据获取器

从 data.wikipathways.org 下载 GMT 格式的通路数据文件。

GMT 文件格式：
    每行代表一个通路，格式为：
    pathway_name\tpathway_id\tgene1\tgene2\t...\tgenen

对应 v1 脚本：wikipathways.R
"""

import re
from pathlib import Path
from typing import Dict, List, Optional, Set
import requests


# WikiPathways 拉丁名到 AllEnricher 物种代码的映射
# 包含 38 个物种：18 个有 GMT 文件 + 20 个只有 GPML
SPECIES_NAME_MAP: Dict[str, str] = {
    # === GMT 物种 (18) ===
    "Anopheles gambiae": "aga",
    "Arabidopsis thaliana": "ath",
    "Bos taurus": "bta",
    "Caenorhabditis elegans": "cel",
    "Canis familiaris": "cfa",
    "Danio rerio": "dre",
    "Drosophila melanogaster": "dme",
    "Equus caballus": "eca",
    "Gallus gallus": "gga",
    "Homo sapiens": "hsa",
    "Mus musculus": "mmu",
    "Pan troglodytes": "ptr",
    "Populus trichocarpa": "ptc",
    "Rattus norvegicus": "rno",
    "Saccharomyces cerevisiae": "sce",
    "Solanum lycopersicum": "sly",
    "Sus scrofa": "ssc",
    "Zea mays": "zma",
    # === GPML Only 物种 (20) ===
    "Acetobacterium woodii": "awo",
    "Bacillus subtilis": "bsu",
    "Beta vulgaris": "bvu",
    "Brassica napus": "bna",
    "Caulobacter vibrioides": "cvi",
    "Citrus sinensis": "csi",
    "Coffea arabica": "car",
    "Daphnia magna": "dma",
    "Escherichia coli": "eco",
    "Gibberella zeae": "gze",
    "Hordeum vulgare": "hvu",
    "Ilex paraguariensis": "ipa",
    "Mycobacterium tuberculosis": "mtu",
    "Oryza sativa": "osa",
    "Paullinia cupana": "pcu",
    "Perilla frutescens": "pfr",
    "Plasmodium falciparum": "pfa",
    "Theobroma cacao": "tcc",
    "Triticum aestivum": "tae",
    "Vitis vinifera": "vvi",
}

# 有 GMT 文件的 18 个物种（拉丁名集合）
SPECIES_WITH_GMT: Set[str] = {
    "Anopheles gambiae",
    "Arabidopsis thaliana",
    "Bos taurus",
    "Caenorhabditis elegans",
    "Canis familiaris",
    "Danio rerio",
    "Drosophila melanogaster",
    "Equus caballus",
    "Gallus gallus",
    "Homo sapiens",
    "Mus musculus",
    "Pan troglodytes",
    "Populus trichocarpa",
    "Rattus norvegicus",
    "Saccharomyces cerevisiae",
    "Solanum lycopersicum",
    "Sus scrofa",
    "Zea mays",
}

# 只有 GPML 文件的 20 个物种（拉丁名集合）
SPECIES_GPML_ONLY: Set[str] = {
    "Acetobacterium woodii",
    "Bacillus subtilis",
    "Beta vulgaris",
    "Brassica napus",
    "Caulobacter vibrioides",
    "Citrus sinensis",
    "Coffea arabica",
    "Daphnia magna",
    "Escherichia coli",
    "Gibberella zeae",
    "Hordeum vulgare",
    "Ilex paraguariensis",
    "Mycobacterium tuberculosis",
    "Oryza sativa",
    "Paullinia cupana",
    "Perilla frutescens",
    "Plasmodium falciparum",
    "Theobroma cacao",
    "Triticum aestivum",
    "Vitis vinifera",
}


class WikiPathwaysFetcher:
    """WikiPathways 数据获取器

    从 data.wikipathways.org 下载 GMT 格式的通路数据。

    Usage::

        fetcher = WikiPathwaysFetcher(basic_dir='./database/basic')
        species_list = fetcher.get_available_species('20240510')
        fetcher.download_all_gmt('20240510')
    """

    BASE_URL = "https://data.wikipathways.org"
    REQUEST_TIMEOUT = 60
    UA = "Mozilla/5.0 (Windows NT 10.0; Win64; x64)"

    def __init__(self, basic_dir: str):
        """
        Args:
            basic_dir: 基础缓存目录（将在此目录下创建 wikipathways 子目录）
        """
        self.basic_dir = Path(basic_dir)
        self.basic_dir.mkdir(parents=True, exist_ok=True)

    def _get_cache_dir(self, version: str) -> Path:
        """获取指定版本的缓存目录

        Args:
            version: 版本号（YYYYMMDD 格式）

        Returns:
            缓存目录路径
        """
        cache_dir = self.basic_dir / "wikipathways" / f"WP{version}"
        cache_dir.mkdir(parents=True, exist_ok=True)
        return cache_dir

    def _detect_latest_version(self) -> str:
        """检测最新版本号

        从 data.wikipathways.org 主页解析最新的日期文件夹。

        Returns:
            最新版本号（YYYYMMDD 格式）

        Raises:
            RuntimeError: 无法获取或解析版本号
        """
        url = self.BASE_URL
        headers = {"User-Agent": self.UA}

        try:
            resp = requests.get(url, headers=headers, timeout=self.REQUEST_TIMEOUT)
            resp.raise_for_status()
        except requests.RequestException as e:
            raise RuntimeError(f"无法访问 WikiPathways 数据站点: {e}") from e

        # 解析 HTML 中的日期文件夹（格式：YYYYMMDD 或 YYYYMMDD/）
        # 匹配类似：<a href="20240510">20240510</a> 或 <a href="20240510/">20240510/</a>
        pattern = r'href=["\'](\d{8})/?["\']'
        matches = re.findall(pattern, resp.text)

        if not matches:
            raise RuntimeError(f"无法从 {url} 解析版本号")

        # 返回最新的版本号（按字符串排序即可，因为都是 YYYYMMDD 格式）
        latest = sorted(matches)[-1]
        return latest

    def get_available_species(self, version: Optional[str] = None) -> List[str]:
        """获取指定版本可用的物种列表

        Args:
            version: 版本号（YYYYMMDD 格式），如果为 None 则自动检测最新版本

        Returns:
            可用的物种拉丁名列表
        """
        if version is None:
            version = self._detect_latest_version()

        url = f"{self.BASE_URL}/{version}/gmt/"
        headers = {"User-Agent": self.UA}

        try:
            resp = requests.get(url, headers=headers, timeout=self.REQUEST_TIMEOUT)
            resp.raise_for_status()
        except requests.RequestException as e:
            raise RuntimeError(f"无法获取物种列表: {e}") from e

        # 解析 GMT 文件链接
        # 文件名格式：wikipathways-{version}-gmt-{Species_Latin_Name}.gmt
        pattern = rf'wikipathways-{version}-gmt-([\w\s]+)\.gmt'
        matches = re.findall(pattern, resp.text)

        # 清理物种名（去除 URL 编码等）
        species_list: List[str] = []
        for match in matches:
            # 将下划线替换为空格
            species_name = match.replace("_", " ")
            species_list.append(species_name)

        return sorted(species_list)

    def download_gmt(
        self,
        species_latin_name: str,
        version: Optional[str] = None,
        overwrite: bool = False
    ) -> Path:
        """下载单个物种的 GMT 文件

        Args:
            species_latin_name: 物种拉丁名（如 "Homo sapiens"）
            version: 版本号（YYYYMMDD 格式），如果为 None 则自动检测最新版本
            overwrite: 是否覆盖已存在的文件

        Returns:
            下载的 GMT 文件路径
        """
        if version is None:
            version = self._detect_latest_version()

        # 构建文件名（空格替换为下划线）
        species_filename = species_latin_name.replace(" ", "_")
        filename = f"wikipathways-{version}-gmt-{species_filename}.gmt"

        url = f"{self.BASE_URL}/{version}/gmt/{filename}"
        cache_dir = self._get_cache_dir(version)
        local_path = cache_dir / filename

        # 检查是否已存在
        if local_path.exists() and not overwrite:
            print(f"|--- 已缓存，跳过: {species_latin_name}")
            return local_path

        # 下载文件
        print(f"|--- 下载: {species_latin_name} ({filename})")
        headers = {"User-Agent": self.UA}

        try:
            resp = requests.get(url, headers=headers, timeout=self.REQUEST_TIMEOUT)
            resp.raise_for_status()
        except requests.RequestException as e:
            raise RuntimeError(f"下载失败 {url}: {e}") from e

        # 保存文件
        with open(local_path, "wb") as f:
            f.write(resp.content)

        print(f"|--- 已保存: {local_path}")
        return local_path

    def download_all_gmt(
        self,
        version: Optional[str] = None,
        overwrite: bool = False
    ) -> Dict[str, Path]:
        """下载所有物种的 GMT 文件

        Args:
            version: 版本号（YYYYMMDD 格式），如果为 None 则自动检测最新版本
            overwrite: 是否覆盖已存在的文件

        Returns:
            物种拉丁名到文件路径的映射
        """
        if version is None:
            version = self._detect_latest_version()

        print(f"\n{'='*60}")
        print(f"WikiPathways GMT 下载 (version={version})")
        print(f"{'='*60}")

        # 获取可用物种列表
        available_species = self.get_available_species(version)
        print(f"|--- 可用物种: {len(available_species)} 个")

        # 下载每个物种的 GMT 文件
        results: Dict[str, Path] = {}
        for species in available_species:
            try:
                path = self.download_gmt(species, version, overwrite)
                results[species] = path
            except RuntimeError as e:
                print(f"|--- 错误: {e}")

        print(f"|--- 下载完成: {len(results)}/{len(available_species)} 个物种")
        return results

    @staticmethod
    def get_species_code(latin_name: str) -> Optional[str]:
        """将物种拉丁名转换为 AllEnricher 物种代码

        Args:
            latin_name: 物种拉丁名（如 "Homo sapiens"）

        Returns:
            AllEnricher 物种代码（如 "hsa"），如果未找到则返回 None
        """
        return SPECIES_NAME_MAP.get(latin_name)

    @staticmethod
    def get_latin_name(species_code: str) -> Optional[str]:
        """将 AllEnricher 物种代码转换为拉丁名

        Args:
            species_code: AllEnricher 物种代码（如 "hsa"）

        Returns:
            物种拉丁名（如 "Homo sapiens"），如果未找到则返回 None
        """
        for latin, code in SPECIES_NAME_MAP.items():
            if code == species_code:
                return latin
        return None

    @staticmethod
    def get_species_data_type(latin_name: str) -> str:
        """获取物种的数据类型

        Args:
            latin_name: 物种拉丁名（如 "Homo sapiens"）

        Returns:
            数据类型：'gmt'（有 GMT 文件）、'gpml'（只有 GPML）、'none'（不支持）
        """
        if latin_name in SPECIES_WITH_GMT:
            return "gmt"
        elif latin_name in SPECIES_GPML_ONLY:
            return "gpml"
        else:
            return "none"

    @staticmethod
    def get_supported_species() -> Dict[str, List[str]]:
        """获取分类后的支持物种列表

        Returns:
            包含分类物种的字典：
            {
                'gmt': [...],      # 18 个有 GMT 的物种拉丁名列表
                'gpml': [...],     # 20 个只有 GPML 的物种拉丁名列表
                'all': [...]       # 全部 38 个物种拉丁名列表
            }
        """
        return {
            "gmt": sorted(list(SPECIES_WITH_GMT)),
            "gpml": sorted(list(SPECIES_GPML_ONLY)),
            "all": sorted(list(SPECIES_NAME_MAP.keys())),
        }

    def get_local_gmt_path(self, species_latin_name: str, version: str) -> Optional[Path]:
        """获取本地 GMT 文件路径（不下载）

        Args:
            species_latin_name: 物种拉丁名
            version: 版本号

        Returns:
            本地文件路径，如果不存在则返回 None
        """
        species_filename = species_latin_name.replace(" ", "_")
        filename = f"wikipathways-{version}-gmt-{species_filename}.gmt"
        cache_dir = self._get_cache_dir(version)
        local_path = cache_dir / filename

        if local_path.exists():
            return local_path
        return None

    def list_cached_versions(self) -> List[str]:
        """列出本地缓存的所有版本

        Returns:
            版本号列表（YYYYMMDD 格式）
        """
        wp_dir = self.basic_dir / "wikipathways"
        if not wp_dir.exists():
            return []

        versions: List[str] = []
        for item in wp_dir.iterdir():
            if item.is_dir() and item.name.startswith("WP"):
                version = item.name[2:]  # 去掉 "WP" 前缀
                if re.match(r"^\d{8}$", version):
                    versions.append(version)

        return sorted(versions)
