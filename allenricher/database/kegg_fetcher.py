"""KEGG REST API 数据获取器

通过 KEGG REST API 批量接口获取物种通路数据，替代 v1 的 HTML 网页爬取。

API 接口：
  - list/pathway/{org}          → 物种所有通路 ID + 名称
  - link/{org}/pathway          → 所有基因-通路关联
  - conv/{org}/ncbi-geneid      → KEGG ID ↔ NCBI Gene ID 映射

性能：仅 3 次 API 请求即可获取全部数据（vs v1 的 ~340 次 HTML 请求）。

对应 v1 脚本：keggMapGrab.R + pathway2tab.pl + makeDB.kegg.v1.1.sh
"""

import gzip
import time
import urllib.request
import urllib.error
from pathlib import Path
from typing import Dict, List, Optional, Tuple


class KEGGFetcher:
    """KEGG REST API 数据获取器

    Usage::

        fetcher = KEGGFetcher(cache_dir='./database/basic/kegg')
        gene2pathway, pathway_summary = fetcher.fetch_species_data('hsa', 'gene_info.gz')
    """

    BASE_URL = "https://rest.kegg.jp"
    # KEGG 要求每秒不超过 10 次请求
    REQUEST_INTERVAL = 0.15  # 秒
    UA = "Mozilla/5.0 (Windows NT 10.0; Win64; x64)"

    def __init__(self, cache_dir: str, overwrite: bool = False):
        """
        Args:
            cache_dir: 缓存目录（存放下载的原始数据和生成的文件）
            overwrite: 是否覆盖已缓存的数据
        """
        self.cache_dir = Path(cache_dir)
        self.cache_dir.mkdir(parents=True, exist_ok=True)
        self.overwrite = overwrite

    # ================================================================
    # 公共接口
    # ================================================================

    def fetch_species_data(
        self,
        species: str,
        gene_info_path: str,
    ) -> Tuple[str, str]:
        """获取物种 KEGG 数据并生成标准格式文件

        通过 3 次 KEGG REST API 请求获取全部数据，结合 gene_info.gz
        中的基因 Symbol 映射，生成与 KEGGParser 兼容的输入文件。

        Args:
            species: KEGG 物种代码（如 'hsa', 'mmu'）
            gene_info_path: gene_info.gz 文件路径

        Returns:
            (gene2pathway_path, pathway_summary_path)
            可直接传给 KEGGParser.build_database()
        """
        print(f"\n{'='*60}")
        print(f"KEGG REST API 数据获取 (species={species})")
        print(f"{'='*60}")

        # Step 1: 获取通路列表
        pathways = self._list_pathways(species)
        print(f"|--- 通路列表: {len(pathways)} 个通路")

        # Step 2: 获取基因-通路关联
        gene_pathway_links = self._get_gene_pathway_links(species)
        print(f"|--- 基因-通路关联: {len(gene_pathway_links)} 个基因")

        # Step 3: 获取 KEGG ID → NCBI Gene ID 映射
        kegg_to_ncbi = self._get_kegg_ncbi_mapping(species)
        print(f"|--- ID 映射: {len(kegg_to_ncbi)} 个 KEGG 基因")

        # Step 4: 构建 NCBI Gene ID → Symbol 映射
        ncbi_to_symbol = self._ncbi_id_to_symbol(gene_info_path)
        print(f"|--- 基因 Symbol: {len(ncbi_to_symbol)} 个")

        # Step 5: 生成 gene2pathway.txt
        gene2pathway_path = self._build_gene2pathway(
            species, gene_pathway_links, kegg_to_ncbi, ncbi_to_symbol, pathways
        )

        # Step 6: 生成 pathway_summary.txt
        pathway_summary_path = self._build_pathway_summary(species, pathways)

        print(f"|--- KEGG 数据获取完成")
        return str(gene2pathway_path), str(pathway_summary_path)

    # ================================================================
    # API 调用
    # ================================================================

    def _api_get(self, endpoint: str) -> str:
        """调用 KEGG REST API

        Args:
            endpoint: API 端点（如 'list/pathway/hsa'）

        Returns:
            响应文本

        Raises:
            RuntimeError: API 请求失败
        """
        url = f"{self.BASE_URL}/{endpoint}"

        # 优先使用 requests（如果有）
        try:
            import requests
            resp = requests.get(url, headers={"User-Agent": self.UA}, timeout=60)
            resp.raise_for_status()
            return resp.text
        except ImportError:
            pass

        # 回退到 urllib
        req = urllib.request.Request(url)
        req.add_header("User-Agent", self.UA)

        try:
            with urllib.request.urlopen(req, timeout=60) as resp:
                return resp.read().decode("utf-8")
        except urllib.error.HTTPError as e:
            raise RuntimeError(f"KEGG API HTTP {e.code}: {url}") from e
        except urllib.error.URLError as e:
            raise RuntimeError(f"KEGG API 网络错误: {e.reason}: {url}") from e

    def _list_pathways(self, species: str) -> List[Tuple[str, str]]:
        """获取物种所有通路列表

        API: list/pathway/{species}
        格式: pathway_id\\tname (with hierarchy prefix)

        Returns:
            [(pathway_id, pathway_name), ...]
            其中 pathway_id 已去除物种前缀（如 '00010'）
        """
        cache_file = self.cache_dir / f"{species}_pathways.txt"

        if cache_file.exists() and not self.overwrite:
            print(f"|--- 已缓存，跳过: list/pathway/{species}")
            pathways = []
            with open(cache_file, "r") as f:
                for line in f:
                    parts = line.strip().split("\t", 1)
                    if len(parts) == 2:
                        pathways.append((parts[0], parts[1]))
            return pathways

        print(f"|--- API: list/pathway/{species}")
        data = self._api_get(f"list/pathway/{species}")
        time.sleep(self.REQUEST_INTERVAL)

        # 解析通路列表，提取分类层级
        pathways = []
        pathway_names: Dict[str, str] = {}  # pathway_id → full name

        for line in data.strip().split("\n"):
            parts = line.split("\t", 1)
            if len(parts) != 2:
                continue
            full_id = parts[0].strip()
            full_name = parts[1].strip()

            # 提取纯通路编号（去除物种前缀）
            pathway_id = full_id.replace(species, "", 1)
            pathway_names[pathway_id] = full_name

        # 保存缓存
        with open(cache_file, "w") as f:
            for pid, pname in pathway_names.items():
                f.write(f"{species}{pid}\t{pname}\n")

        return [(pid, pname) for pid, pname in pathway_names.items()]

    def _get_gene_pathway_links(self, species: str) -> Dict[str, List[str]]:
        """获取所有基因-通路关联（批量接口）

        API: link/{species}/pathway
        格式: path:{species}00010\\t{species}:10327

        Returns:
            {kegg_gene_id: [pathway_id, ...], ...}
        """
        cache_file = self.cache_dir / f"{species}_gene_pathway_links.txt"

        if cache_file.exists() and not self.overwrite:
            print(f"|--- 已缓存，跳过: link/{species}/pathway")
            links = {}
            with open(cache_file, "r") as f:
                for line in f:
                    parts = line.strip().split("\t")
                    if len(parts) == 2:
                        gene, pw = parts
                        links.setdefault(gene, []).append(pw)
            return links

        print(f"|--- API: link/{species}/pathway")
        data = self._api_get(f"link/{species}/pathway")
        time.sleep(self.REQUEST_INTERVAL)

        links: Dict[str, List[str]] = {}
        for line in data.strip().split("\n"):
            parts = line.split("\t")
            if len(parts) != 2:
                continue
            # path:hsa00010 → 00010
            pw = parts[0].replace(f"path:{species}", "")
            # hsa:10327 → 10327
            gene = parts[1].replace(f"{species}:", "")
            links.setdefault(gene, []).append(pw)

        # 保存缓存
        with open(cache_file, "w") as f:
            for gene, pws in links.items():
                for pw in pws:
                    f.write(f"{gene}\t{pw}\n")

        return links

    def _get_kegg_ncbi_mapping(self, species: str) -> Dict[str, str]:
        """获取 KEGG ID → NCBI Gene ID 映射

        API: conv/{species}/ncbi-geneid
        格式: ncbi-geneid:7157\\thsa:10327

        Returns:
            {kegg_gene_id: ncbi_gene_id, ...}
        """
        cache_file = self.cache_dir / f"{species}_kegg_ncbi_map.txt"

        if cache_file.exists() and not self.overwrite:
            print(f"|--- 已缓存，跳过: conv/{species}/ncbi-geneid")
            mapping = {}
            with open(cache_file, "r") as f:
                for line in f:
                    parts = line.strip().split("\t")
                    if len(parts) == 2:
                        ncbi, kegg = parts
                        mapping[kegg] = ncbi
            return mapping

        print(f"|--- API: conv/{species}/ncbi-geneid")
        data = self._api_get(f"conv/{species}/ncbi-geneid")
        time.sleep(self.REQUEST_INTERVAL)

        mapping: Dict[str, str] = {}
        for line in data.strip().split("\n"):
            parts = line.split("\t")
            if len(parts) != 2:
                continue
            ncbi = parts[0].replace("ncbi-geneid:", "")
            kegg = parts[1].replace(f"{species}:", "")
            mapping[kegg] = ncbi

        # 保存缓存
        with open(cache_file, "w") as f:
            for kegg, ncbi in mapping.items():
                f.write(f"{ncbi}\t{kegg}\n")

        return mapping

    def _ncbi_id_to_symbol(self, gene_info_path: str) -> Dict[str, str]:
        """从 gene_info.gz 构建 NCBI Gene ID → Symbol 映射

        Args:
            gene_info_path: gene_info.gz 文件路径

        Returns:
            {ncbi_gene_id: gene_symbol, ...}
        """
        mapping: Dict[str, str] = {}
        opener = gzip.open if gene_info_path.endswith(".gz") else open

        with opener(gene_info_path, "rt", encoding="utf-8") as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                parts = line.split("\t")
                if len(parts) >= 3:
                    mapping[parts[1]] = parts[2]  # gene_id → symbol

        return mapping

    # ================================================================
    # 文件生成
    # ================================================================

    def _build_gene2pathway(
        self,
        species: str,
        gene_pathway_links: Dict[str, List[str]],
        kegg_to_ncbi: Dict[str, str],
        ncbi_to_symbol: Dict[str, str],
        pathways: List[Tuple[str, str]],
    ) -> Path:
        """生成 gene2pathway.txt

        格式: gene_symbol\\tentrez_id\\tpathway_id\\tpathway_name
        与 KEGGParser.build_database() 兼容。
        """
        pathway_names = {pid: pname for pid, pname in pathways}
        out_file = self.cache_dir / f"{species}_gene2pathway.txt"

        n = 0
        with open(out_file, "w", encoding="utf-8") as f:
            for kegg_gene_id, pw_ids in sorted(gene_pathway_links.items()):
                # KEGG ID → NCBI ID → Symbol
                ncbi_id = kegg_to_ncbi.get(kegg_gene_id)
                if not ncbi_id:
                    continue
                symbol = ncbi_to_symbol.get(ncbi_id)
                if not symbol:
                    continue

                for pw_id in pw_ids:
                    pw_name = pathway_names.get(pw_id, pw_id)
                    # 去除名称中的 " - Homo sapiens (human)" 后缀
                    pw_name = self._clean_pathway_name(pw_name)
                    f.write(f"{symbol}\t{ncbi_id}\t{pw_id}\t{pw_name}\n")
                    n += 1

        print(f"|--- gene2pathway.txt: {n} 条关联")
        return out_file

    def _build_pathway_summary(
        self,
        species: str,
        pathways: List[Tuple[str, str]],
    ) -> Path:
        """生成 pathway_summary.txt

        格式: Category\\tSubcategory\\tpathway_id\\tpathway_name\\turl
        与 KEGGParser.build_database() 兼容。
        """
        out_file = self.cache_dir / f"{species}_pathway_summary.txt"

        # 解析通路层级（从 list/pathway 的缩进推断分类）
        # KEGG 的 list 输出中，分类行以 pathway ID 开头但无数字编号
        # 实际上 list/pathway/{org} 只返回通路条目，不包含分类层级
        # 分类信息需要从通路名称推断或使用全局分类

        # 从 KEGG API 获取真实分类信息
        brite_categories = self._get_brite_categories(species, pathways)

        n = 0
        with open(out_file, "w", encoding="utf-8") as f:
            for pw_id, pw_name in pathways:
                pw_name_clean = self._clean_pathway_name(pw_name)
                url = f"https://www.kegg.jp/entry/{species}{pw_id}"

                # 获取通路分类（从 KEGG API CLASS 字段）
                category, subcategory = brite_categories.get(pw_id, ("Uncategorized", "Uncategorized"))

                f.write(f"{category}\t{subcategory}\t{pw_id}\t{pw_name_clean}\t{url}\n")
                n += 1

        print(f"|--- pathway_summary.txt: {n} 个通路")
        return out_file

    @staticmethod
    def _clean_pathway_name(name: str) -> str:
        """清理通路名称

        去除 " - Homo sapiens (human)" 等物种后缀。
        """
        import re
        # 匹配 " - Species name (common name)" 后缀
        # 例如: "Glycolysis / Gluconeogenesis - Homo sapiens (human)"
        cleaned = re.sub(r"\s*-\s*[\w\s]+\(\w+\)\s*$", "", name)
        return cleaned.strip()

    def _get_brite_categories(self, species: str, pathways: List[Tuple[str, str]]) -> Dict[str, Tuple[str, str]]:
        """获取 KEGG 通路分类映射

        从 KEGG API 获取每个通路的 CLASS 信息，包含真实的层级分类。

        Args:
            species: 物种代码
            pathways: 通路列表 [(pathway_id, pathway_name), ...]

        Returns:
            {pathway_id: (category, subcategory), ...}
        """
        cache_file = self.cache_dir / f"{species}_pathway_classes.txt"

        # 如果已有缓存，直接读取
        if cache_file.exists() and not self.overwrite:
            print(f"|--- 从缓存加载通路分类")
            categories: Dict[str, Tuple[str, str]] = {}
            with open(cache_file, 'r', encoding='utf-8') as f:
                for line in f:
                    parts = line.strip().split('\t')
                    if len(parts) == 3:
                        pw_id, cat, subcat = parts
                        categories[pw_id] = (cat, subcat)
            return categories

        print(f"|--- 从 KEGG API 获取通路分类信息")
        categories: Dict[str, Tuple[str, str]] = {}

        try:
            # KEGG API 支持批量查询，每次最多 10 个通路
            # 格式: get/pathway/hsa04110+hsa00010
            pw_ids = [pw_id for pw_id, _ in pathways]
            batch_size = 10
            fetched = 0

            for i in range(0, len(pw_ids), batch_size):
                batch = pw_ids[i:i + batch_size]
                ids_str = '+'.join(batch)
                data = self._api_get(f"get/{ids_str}")

                # 解析 CLASS 字段
                # 批量返回格式：多个通路数据拼接，每个通路以 ENTRY 开始，以 /// 结束
                # ENTRY -> NAME -> CLASS -> ... -> ///
                lines = data.split('\n')

                for pw_id in batch:
                    found = False
                    in_entry = False
                    class_line = None

                    for line in lines:
                        # 检测 ENTRY 行，确认是当前通路
                        if line.startswith("ENTRY") and pw_id in line:
                            in_entry = True
                            continue

                        # 检测通路结束标记
                        if line.startswith("///"):
                            in_entry = False
                            continue

                        # 在当前通路中查找 CLASS
                        if in_entry and line.startswith("CLASS"):
                            class_line = line.replace("CLASS", "").strip()
                            break

                    if class_line:
                        # 格式: "Cellular Processes; Cell growth and death"
                        parts = class_line.split("; ")
                        category = parts[0].replace(" ", "_") if parts else "Uncategorized"
                        subcategory = parts[1].replace(" ", "_") if len(parts) > 1 else "Uncategorized"
                        categories[pw_id] = (category, subcategory)
                    else:
                        # 未找到 CLASS，使用默认分类
                        categories[pw_id] = ("Uncategorized", "Uncategorized")

                    fetched += 1
                    if fetched % 50 == 0:
                        print(f"|--- 已获取 {fetched}/{len(pw_ids)} 个通路分类")

                time.sleep(self.REQUEST_INTERVAL)

        except Exception as e:
            print(f"|--- API 获取失败，使用备用分类映射: {e}")
            # 使用硬编码的分类映射
            categories = self._get_hardcoded_categories()
            # 只保留有分类的
            for pw_id, _ in pathways:
                if pw_id not in categories:
                    categories[pw_id] = ("Uncategorized", "Uncategorized")

        # 保存缓存
        with open(cache_file, 'w', encoding='utf-8') as f:
            for pw_id, (cat, subcat) in categories.items():
                f.write(f"{pw_id}\t{cat}\t{subcat}\n")

        print(f"|--- 获取了 {len(categories)} 个通路分类")
        return categories

    @staticmethod
    def _get_hardcoded_categories() -> Dict[str, Tuple[str, str]]:
        """获取硬编码的 KEGG 通路分类映射

        作为 KEGG API 不可用时的备用方案。
        包含人类 KEGG 通路的常见分类。
        注意: 键使用 `hsa` 前缀（如 hsa00010）

        Returns:
            {pathway_id: (category, subcategory), ...}
        """
        # 硬编码的分类映射（使用 hsa 前缀）
        categories = {
            # 代谢通路
            "hsa00010": ("Metabolism", "Carbohydrate_Metabolism"),
            "hsa00020": ("Metabolism", "Carbohydrate_Metabolism"),
            "hsa00030": ("Metabolism", "Carbohydrate_Metabolism"),
            "hsa00040": ("Metabolism", "Carbohydrate_Metabolism"),
            "hsa00051": ("Metabolism", "Carbohydrate_Metabolism"),
            "hsa00052": ("Metabolism", "Carbohydrate_Metabolism"),
            "hsa00053": ("Metabolism", "Carbohydrate_Metabolism"),
            "hsa00061": ("Metabolism", "Lipid_Metabolism"),
            "hsa00062": ("Metabolism", "Lipid_Metabolism"),
            "hsa00071": ("Metabolism", "Lipid_Metabolism"),
            "hsa00100": ("Metabolism", "Lipid_Metabolism"),
            "hsa00120": ("Metabolism", "Lipid_Metabolism"),
            "hsa00130": ("Metabolism", "Lipid_Metabolism"),
            "hsa00140": ("Metabolism", "Lipid_Metabolism"),
            "hsa00190": ("Metabolism", "Energy_Metabolism"),
            "hsa00220": ("Metabolism", "Amino_Acid_Metabolism"),
            "hsa00230": ("Metabolism", "Nucleotide_Metabolism"),
            "hsa00232": ("Metabolism", "Metabolism_of_Cofactors"),
            "hsa00240": ("Metabolism", "Nucleotide_Metabolism"),
            "hsa00250": ("Metabolism", "Amino_Acid_Metabolism"),
            "hsa00260": ("Metabolism", "Amino_Acid_Metabolism"),
            "hsa00270": ("Metabolism", "Amino_Acid_Metabolism"),
            "hsa00280": ("Metabolism", "Amino_Acid_Metabolism"),
            "hsa00290": ("Metabolism", "Amino_Acid_Metabolism"),
            "hsa00300": ("Metabolism", "Amino_Acid_Metabolism"),
            "hsa00310": ("Metabolism", "Amino_Acid_Metabolism"),
            "hsa00330": ("Metabolism", "Amino_Acid_Metabolism"),
            "hsa00340": ("Metabolism", "Amino_Acid_Metabolism"),
            "hsa00350": ("Metabolism", "Amino_Acid_Metabolism"),
            "hsa00360": ("Metabolism", "Amino_Acid_Metabolism"),
            "hsa00380": ("Metabolism", "Amino_Acid_Metabolism"),
            "hsa00400": ("Metabolism", "Amino_Acid_Metabolism"),
            "hsa00410": ("Metabolism", "Metabolism_of_Other_Amino_Acids"),
            "hsa00430": ("Metabolism", "Metabolism_of_Other_Amino_Acids"),
            "hsa00440": ("Metabolism", "Metabolism_of_Other_Amino_Acids"),
            "hsa00450": ("Metabolism", "Metabolism_of_Other_Amino_Acids"),
            "hsa00460": ("Metabolism", "Metabolism_of_Other_Amino_Acids"),
            "hsa00470": ("Metabolism", "Metabolism_of_Other_Amino_Acids"),
            "hsa00471": ("Metabolism", "Metabolism_of_Other_Amino_Acids"),
            "hsa00472": ("Metabolism", "Metabolism_of_Other_Amino_Acids"),
            "hsa00480": ("Metabolism", "Metabolism_of_Other_Amino_Acids"),
            "hsa00500": ("Metabolism", "Glycan_Biosynthesis"),
            "hsa00510": ("Metabolism", "Glycan_Biosynthesis"),
            "hsa00511": ("Metabolism", "Glycan_Biosynthesis"),
            "hsa00512": ("Metabolism", "Glycan_Biosynthesis"),
            "hsa00513": ("Metabolism", "Glycan_Biosynthesis"),
            "hsa00514": ("Metabolism", "Glycan_Biosynthesis"),
            "hsa00520": ("Metabolism", "Glycan_Biosynthesis"),
            "hsa00531": ("Metabolism", "Glycan_Biosynthesis"),
            "hsa00532": ("Metabolism", "Glycan_Biosynthesis"),
            "hsa00533": ("Metabolism", "Glycan_Biosynthesis"),
            "hsa00534": ("Metabolism", "Glycan_Biosynthesis"),
            "hsa00540": ("Metabolism", "Glycan_Biosynthesis"),
            "hsa00550": ("Metabolism", "Glycan_Biosynthesis"),
            "hsa00561": ("Metabolism", "Lipid_Metabolism"),
            "hsa00562": ("Metabolism", "Lipid_Metabolism"),
            "hsa00563": ("Metabolism", "Lipid_Metabolism"),
            "hsa00564": ("Metabolism", "Lipid_Metabolism"),
            "hsa00565": ("Metabolism", "Lipid_Metabolism"),
            "hsa00590": ("Metabolism", "Lipid_Metabolism"),
            "hsa00591": ("Metabolism", "Lipid_Metabolism"),
            "hsa00592": ("Metabolism", "Lipid_Metabolism"),
            "hsa00600": ("Metabolism", "Carbohydrate_Metabolism"),
            "hsa00601": ("Metabolism", "Carbohydrate_Metabolism"),
            "hsa00620": ("Metabolism", "Carbohydrate_Metabolism"),
            "hsa00625": ("Metabolism", "Carbohydrate_Metabolism"),
            "hsa00626": ("Metabolism", "Carbohydrate_Metabolism"),
            "hsa00630": ("Metabolism", "One_Carbon_Metabolism"),
            "hsa00640": ("Metabolism", "Carbohydrate_Metabolism"),
            "hsa00650": ("Metabolism", "Carbohydrate_Metabolism"),
            "hsa00670": ("Metabolism", "One_Carbon_Metabolism"),
            "hsa00790": ("Metabolism", "Metabolism_of_Cofactors"),
            "hsa00830": ("Metabolism", "Metabolism_of_Cofactors"),
            "hsa00860": ("Metabolism", "Metabolism_of_Cofactors"),
            "hsa00900": ("Metabolism", "Terpenoid_Backbone_Biosynthesis"),
            "hsa00910": ("Metabolism", "Nitrogen_Metabolism"),
            "hsa00920": ("Metabolism", "Sulfur_Metabolism"),
            "hsa00970": ("Metabolism", "Translation"),
            "hsa00980": ("Metabolism", "Metabolism_of_Cofactors"),
            "hsa00982": ("Metabolism", "Xenobiotics_Biodegradation"),
            "hsa00983": ("Metabolism", "Xenobiotics_Biodegradation"),
            # 细胞过程
            "hsa04110": ("Cellular_Processes", "Cell_Cycle"),
            "hsa04111": ("Cellular_Processes", "Cell_Cycle"),
            "hsa04112": ("Cellular_Processes", "Cell_Cycle"),
            "hsa04113": ("Cellular_Processes", "Cell_Cycle"),
            "hsa04114": ("Cellular_Processes", "Cell_Cycle"),
            "hsa04115": ("Cellular_Processes", "Cell_Cycle"),
            "hsa04120": ("Cellular_Processes", "Cellular_Senescence"),
            "hsa04122": ("Cellular_Processes", "Transport"),
            "hsa04130": ("Cellular_Processes", "Folding_and_Degradation"),
            "hsa04140": ("Cellular_Processes", "Transport"),
            "hsa04141": ("Cellular_Processes", "Transport"),
            "hsa04142": ("Cellular_Processes", "Transport"),
            "hsa04144": ("Cellular_Processes", "Endocytosis"),
            "hsa04145": ("Cellular_Processes", "Phagocytosis"),
            "hsa04146": ("Cellular_Processes", "Autophagy"),
            "hsa04150": ("Cellular_Processes", "Signal_Transduction"),
            "hsa04151": ("Cellular_Processes", "Signal_Transduction"),
            # 信号传导
            "hsa04010": ("Environmental_Information_Processing", "Signal_Transduction"),
            "hsa04012": ("Environmental_Information_Processing", "Signal_Transduction"),
            "hsa04014": ("Environmental_Information_Processing", "Signal_Transduction"),
            "hsa04015": ("Environmental_Information_Processing", "Signal_Transduction"),
            "hsa04020": ("Environmental_Information_Processing", "Signal_Transduction"),
            "hsa04022": ("Environmental_Information_Processing", "Signal_Transduction"),
            "hsa04024": ("Environmental_Information_Processing", "Signal_Transduction"),
            "hsa04066": ("Environmental_Information_Processing", "Signal_Transduction"),
            "hsa04068": ("Environmental_Information_Processing", "Signal_Transduction"),
            "hsa04070": ("Environmental_Information_Processing", "Signal_Transduction"),
            "hsa04071": ("Environmental_Information_Processing", "Signal_Transduction"),
            "hsa04072": ("Environmental_Information_Processing", "Signal_Transduction"),
            "hsa04010": ("Environmental_Information_Processing", "Signal_Transduction"),
            "hsa04210": ("Environmental_Information_Processing", "Signal_Transduction"),
            "hsa04211": ("Environmental_Information_Processing", "Signal_Transduction"),
            "hsa04213": ("Environmental_Information_Processing", "Signal_Transduction"),
            "hsa04215": ("Environmental_Information_Processing", "Signal_Transduction"),
            "hsa04217": ("Environmental_Information_Processing", "Signal_Transduction"),
            "hsa04218": ("Environmental_Information_Processing", "Signal_Transduction"),
            # DNA 复制和修复
            "hsa03030": ("Genetic_Information_Processing", "Replication_and_Repair"),
            "hsa03040": ("Genetic_Information_Processing", "Replication_and_Repair"),
            "hsa03410": ("Genetic_Information_Processing", "Replication_and_Repair"),
            "hsa03420": ("Genetic_Information_Processing", "Replication_and_Repair"),
            "hsa03430": ("Genetic_Information_Processing", "Replication_and_Repair"),
            "hsa03440": ("Genetic_Information_Processing", "Replication_and_Repair"),
            "hsa03450": ("Genetic_Information_Processing", "Replication_and_Repair"),
            "hsa03460": ("Genetic_Information_Processing", "Replication_and_Repair"),
            # 转录
            "hsa03020": ("Genetic_Information_Processing", "Transcription"),
            # 翻译
            "hsa03010": ("Genetic_Information_Processing", "Translation"),
            "hsa03013": ("Genetic_Information_Processing", "Translation"),
            "hsa03015": ("Genetic_Information_Processing", "Translation"),
            # 折叠和降解
            "hsa04130": ("Genetic_Information_Processing", "Folding_and_Degradation"),
            # 免疫系统
            "hsa04612": ("Organismal_Systems", "Immune_System"),
            "hsa04620": ("Organismal_Systems", "Immune_System"),
            "hsa04621": ("Organismal_Systems", "Immune_System"),
            "hsa04622": ("Organismal_Systems", "Immune_System"),
            "hsa04623": ("Organismal_Systems", "Immune_System"),
            "hsa04625": ("Organismal_Systems", "Immune_System"),
            "hsa04630": ("Organismal_Systems", "Immune_System"),
            "hsa04640": ("Organismal_Systems", "Immune_System"),
            "hsa04650": ("Organismal_Systems", "Immune_System"),
            "hsa04611": ("Organismal_Systems", "Immune_System"),
            "hsa04657": ("Human_Diseases", "Infectious_Disease"),
            # 内分泌系统
            "hsa04910": ("Organismal_Systems", "Endocrine_System"),
            "hsa04911": ("Organismal_Systems", "Endocrine_System"),
            "hsa04912": ("Organismal_Systems", "Endocrine_System"),
            "hsa04913": ("Organismal_Systems", "Endocrine_System"),
            "hsa04914": ("Organismal_Systems", "Endocrine_System"),
            "hsa04915": ("Organismal_Systems", "Endocrine_System"),
            "hsa04916": ("Organismal_Systems", "Endocrine_System"),
            "hsa04917": ("Organismal_Systems", "Endocrine_System"),
            "hsa04918": ("Organismal_Systems", "Endocrine_System"),
            "hsa04919": ("Organismal_Systems", "Endocrine_System"),
            "hsa04920": ("Organismal_Systems", "Endocrine_System"),
            "hsa04921": ("Organismal_Systems", "Endocrine_System"),
            "hsa04922": ("Organismal_Systems", "Endocrine_System"),
            "hsa04923": ("Organismal_Systems", "Endocrine_System"),
            "hsa04924": ("Organismal_Systems", "Endocrine_System"),
            "hsa04925": ("Organismal_Systems", "Endocrine_System"),
            "hsa04926": ("Organismal_Systems", "Endocrine_System"),
            "hsa04927": ("Organismal_Systems", "Endocrine_System"),
            "hsa04928": ("Organismal_Systems", "Endocrine_System"),
            "hsa04929": ("Organismal_Systems", "Endocrine_System"),
            "hsa04930": ("Organismal_Systems", "Endocrine_System"),
            "hsa04931": ("Organismal_Systems", "Endocrine_System"),
            "hsa04932": ("Organismal_Systems", "Endocrine_System"),
            "hsa04933": ("Organismal_Systems", "Endocrine_System"),
            "hsa04934": ("Organismal_Systems", "Endocrine_System"),
            "hsa04935": ("Organismal_Systems", "Endocrine_System"),
            # 消化系统
            "hsa04970": ("Organismal_Systems", "Digestive_System"),
            "hsa04971": ("Organismal_Systems", "Digestive_System"),
            "hsa04972": ("Organismal_Systems", "Digestive_System"),
            "hsa04973": ("Organismal_Systems", "Digestive_System"),
            "hsa04974": ("Organismal_Systems", "Digestive_System"),
            "hsa04975": ("Organismal_Systems", "Digestive_System"),
            "hsa04976": ("Organismal_Systems", "Digestive_System"),
            "hsa04977": ("Organismal_Systems", "Digestive_System"),
            "hsa04978": ("Organismal_Systems", "Digestive_System"),
            "hsa04979": ("Organismal_Systems", "Digestive_System"),
            # 神经系统
            "hsa04710": ("Organismal_Systems", "Nervous_System"),
            "hsa04711": ("Organismal_Systems", "Nervous_System"),
            "hsa04720": ("Organismal_Systems", "Nervous_System"),
            "hsa04721": ("Organismal_Systems", "Nervous_System"),
            "hsa04722": ("Organismal_Systems", "Nervous_System"),
            "hsa04723": ("Organismal_Systems", "Nervous_System"),
            "hsa04724": ("Organismal_Systems", "Nervous_System"),
            "hsa04725": ("Organismal_Systems", "Nervous_System"),
            "hsa04726": ("Organismal_Systems", "Nervous_System"),
            "hsa04727": ("Organismal_Systems", "Nervous_System"),
            "hsa04728": ("Organismal_Systems", "Nervous_System"),
            "hsa04730": ("Organismal_Systems", "Nervous_System"),
            "hsa04740": ("Organismal_Systems", "Sensory_System"),
            "hsa04742": ("Organismal_Systems", "Sensory_System"),
            "hsa04744": ("Organismal_Systems", "Sensory_System"),
            "hsa04750": ("Organismal_Systems", "Sensory_System"),
            "hsa04710": ("Organismal_Systems", "Nervous_System"),
            # 传染病
            "hsa05160": ("Human_Diseases", "Infectious_Disease"),
            "hsa05161": ("Human_Diseases", "Infectious_Disease"),
            "hsa05162": ("Human_Diseases", "Infectious_Disease"),
            "hsa05164": ("Human_Diseases", "Infectious_Disease"),
            "hsa05165": ("Human_Diseases", "Infectious_Disease"),
            "hsa05166": ("Human_Diseases", "Infectious_Disease"),
            "hsa05167": ("Human_Diseases", "Infectious_Disease"),
            "hsa05168": ("Human_Diseases", "Infectious_Disease"),
            "hsa05169": ("Human_Diseases", "Infectious_Disease"),
            "hsa05170": ("Human_Diseases", "Infectious_Disease"),
            "hsa05171": ("Human_Diseases", "Infectious_Disease"),
            # 癌症
            "hsa05200": ("Human_Diseases", "Cancer"),
            "hsa05210": ("Human_Diseases", "Cancer"),
            "hsa05211": ("Human_Diseases", "Cancer"),
            "hsa05212": ("Human_Diseases", "Cancer"),
            "hsa05213": ("Human_Diseases", "Cancer"),
            "hsa05214": ("Human_Diseases", "Cancer"),
            "hsa05215": ("Human_Diseases", "Cancer"),
            "hsa05216": ("Human_Diseases", "Cancer"),
            "hsa05217": ("Human_Diseases", "Cancer"),
            "hsa05218": ("Human_Diseases", "Cancer"),
            "hsa05219": ("Human_Diseases", "Cancer"),
            "hsa05220": ("Human_Diseases", "Cancer"),
            "hsa05221": ("Human_Diseases", "Cancer"),
            "hsa05222": ("Human_Diseases", "Cancer"),
            "hsa05223": ("Human_Diseases", "Cancer"),
            "hsa05224": ("Human_Diseases", "Cancer"),
            "hsa05225": ("Human_Diseases", "Cancer"),
            "hsa05226": ("Human_Diseases", "Cancer"),
            "hsa05230": ("Human_Diseases", "Cancer"),
            "hsa05231": ("Human_Diseases", "Cancer"),
            "hsa05232": ("Human_Diseases", "Cancer"),
            "hsa05235": ("Human_Diseases", "Cancer"),
        }
        return categories
