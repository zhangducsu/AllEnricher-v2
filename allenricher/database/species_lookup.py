"""
Species lookup module for AllEnricher v2.0

物种检索模块，提供物种信息的查询功能。
支持通过 KEGG 代码、拉丁名、TaxID 等方式检索物种配置信息。
"""

from dataclasses import dataclass
from typing import Optional, Dict, List

from allenricher.core.config import SPECIES_CONFIGS, SpeciesConfig


@dataclass
class SpeciesInfo:
    """物种信息数据类
    
    存储单个物种的完整信息，包括 KEGG 代码、拉丁名、TaxID 等。
    
    Attributes:
        kegg_code: KEGG 物种代码（如 'hsa'）
        latin_name: 物种拉丁学名（如 'Homo sapiens'）
        taxonomy_id: NCBI 分类学 ID（如 9606）
        display_name: 显示名称（如 'Human'）
    """
    kegg_code: str
    latin_name: str
    taxonomy_id: int
    display_name: str = ""


# 内置 TaxID 映射（用于反向查找）
BUILTIN_TAXID_MAP: Dict[int, str] = {
    9606: "hsa",    # Human
    10090: "mmu",   # Mouse
    10116: "rno",   # Rat
    7955: "dre",    # Zebrafish
    7227: "dme",    # Fruit fly
    6239: "cel",    # C. elegans
    9823: "ssc",    # Pig
    9913: "bta",    # Cow
    9031: "gga",    # Chicken
    8364: "xtr",    # Xenopus
    9615: "cfa",    # Dog
    44689: "ddi",   # Dictyostelium
    1772: "mtu",    # M. tuberculosis
    5833: "pfa",    # P. falciparum
    4932: "sce",    # S. cerevisiae
    4896: "spo",    # S. pombe
}


class SpeciesLookup:
    """物种检索类
    
    提供物种信息的多种查询方式，支持离线模式（仅使用内置数据）
    和在线模式（从 KEGG API 获取更多物种信息）。
    
    Attributes:
        auto_load: 是否自动从网络加载物种数据
        loaded: 是否已完成数据加载
        species_data: 物种数据字典，键为 KEGG 代码
    """
    
    def __init__(self, auto_load: bool = True):
        """初始化物种检索器
        
        Args:
            auto_load: 是否自动从网络加载物种数据（默认 True）
        """
        self.auto_load = auto_load
        self.loaded = False
        self.species_data: Dict[str, SpeciesInfo] = {}
        
        if auto_load:
            self._load_builtin_species()
            self.loaded = True
    
    def _load_builtin_species(self) -> None:
        """加载内置物种数据
        
        从 SPECIES_CONFIGS 加载预定义的物种配置信息。
        """
        for kegg_code, config in SPECIES_CONFIGS.items():
            self.species_data[kegg_code] = SpeciesInfo(
                kegg_code=config.kegg_code,
                latin_name=config.name,
                taxonomy_id=config.taxonomy_id,
                display_name=config.display_name
            )
    
    def lookup_by_kegg_code(self, kegg_code: str) -> Optional[SpeciesInfo]:
        """通过 KEGG 代码检索物种
        
        Args:
            kegg_code: KEGG 物种代码（如 'hsa'）
        
        Returns:
            SpeciesInfo: 物种信息，如果未找到则返回 None
        """
        return self.species_data.get(kegg_code)
    
    def lookup_by_latin_name(self, latin_name: str) -> Optional[SpeciesInfo]:
        """通过拉丁学名检索物种
        
        Args:
            latin_name: 物种拉丁学名（如 'Homo sapiens'）
        
        Returns:
            SpeciesInfo: 物种信息，如果未找到则返回 None
        """
        for info in self.species_data.values():
            if info.latin_name.lower() == latin_name.lower():
                return info
        return None
    
    def lookup_by_taxid(self, taxonomy_id: int) -> Optional[SpeciesInfo]:
        """通过 NCBI TaxID 检索物种
        
        Args:
            taxonomy_id: NCBI 分类学 ID（如 9606）
        
        Returns:
            SpeciesInfo: 物种信息，如果未找到则返回 None
        """
        kegg_code = BUILTIN_TAXID_MAP.get(taxonomy_id)
        if kegg_code:
            return self.species_data.get(kegg_code)
        return None
    
    def get_all_species(self) -> List[SpeciesInfo]:
        """获取所有已加载的物种信息
        
        Returns:
            List[SpeciesInfo]: 所有物种信息列表
        """
        return list(self.species_data.values())
    
    def get_species_count(self) -> int:
        """获取已加载的物种数量
        
        Returns:
            int: 物种数量
        """
        return len(self.species_data)