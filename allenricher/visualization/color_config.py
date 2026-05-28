"""
AllEnricher 颜色配置系统

提供统一的颜色管理，支持多种预设配色方案。
所有图表颜色必须从此模块获取，禁止硬编码。
"""

from typing import Dict, List, Optional


# =============================================================================
# Paul Tol 色板 - 色盲友好、高对比度
# =============================================================================

TOL_BRIGHT = [
    "#4477AA", "#66CCEE", "#228833", "#CCBB44", "#EE6677",
    "#AA3377", "#BBBBBB"
]

TOL_HIGH_CONTRAST = [
    "#004488", "#DDAA33", "#BB5566"
]

TOL_VIBRANT = [
    "#0077BB", "#33BBEE", "#009988", "#EE7733", "#CC3311",
    "#EE3377", "#BBBBBB"
]

TOL_MUTED = [
    "#332288", "#88CCEE", "#44AA99", "#117733", "#999933",
    "#DDCC77", "#CC6677", "#882255", "#AA4499", "#DDDDDD"
]

TOL_MEDIUM_CONTRAST = [
    "#6699CC", "#004488", "#EECC66", "#994455", "#997700",
    "#EE99AA"
]

TOL_LIGHT = [
    "#77AADD", "#99DDFF", "#44BB99", "#BBCC77", "#AAAA00",
    "#EEDD88", "#EE8866", "#FFAABB", "#DDDDDD"
]

TOL_SUNSET = [
    "#364B9A", "#4A7BB7", "#6EA6CD", "#98CAE1", "#C2E4EF",
    "#EAECCC", "#FEDA8B", "#FDB366", "#F67E4B", "#DD3D2D",
    "#A50026"
]

TOL_BURGA = [
    "#F7F4F9", "#E7E1EF", "#D4B9DA", "#C994C7", "#DF65B0",
    "#E7298A", "#CE1256", "#980043", "#67001F"
]

TOL_PRGn = [
    "#762A83", "#9970AB", "#C2A5CF", "#E7D4E8", "#F7F7F7",
    "#D9F0D3", "#ACD39E", "#5AAE61", "#1B7837"
]


# =============================================================================
# Okabe-Ito 色板 - 色盲友好标准
# =============================================================================

OKABE_ITO = [
    "#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442",
    "#0072B2", "#D55E00", "#CC79A7"
]


# =============================================================================
# GO/KEGG 类别色 - 生物信息学专用
# =============================================================================

GO_BP_COLORS = [
    "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",
    "#FFFF33", "#A65628", "#F781BF", "#999999"
]

GO_CC_COLORS = [
    "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854",
    "#FFD92F", "#E5C494", "#B3B3B3"
]

GO_MF_COLORS = [
    "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3",
    "#FDB462", "#B3DE69", "#FCCDE5", "#D9D9D9"
]

KEGG_PATHWAY_COLORS = [
    "#1F77B4", "#FF7F0E", "#2CA02C", "#D62728", "#9467BD",
    "#8C564B", "#E377C2", "#7F7F7F", "#BCBD22", "#17BECF"
]


# =============================================================================
# 火山图专用色
# =============================================================================

VOLCANO_COLORS = {
    "up": "#DC143C",      # 上调 - 深红
    "down": "#4169E1",    # 下调 - 皇家蓝
    "ns": "#808080",      # 不显著 - 灰色
}


# =============================================================================
# 科研期刊风格色板
# =============================================================================

NATURE_COLORS = [
    "#0C5DA5", "#FF9500", "#00B945", "#FF2C00", "#845B97",
    "#474747", "#9E9E9E"
]

SCIENCE_COLORS = [
    "#1F77B4", "#FF7F0E", "#2CA02C", "#D62728", "#9467BD",
    "#8C564B", "#E377C2", "#7F7F7F", "#BCBD22", "#17BECF"
]

CELL_COLORS = [
    "#0072B2", "#D55E00", "#CC79A7", "#F0E442", "#009E73",
    "#56B4E9", "#E69F00", "#000000"
]

LANCET_COLORS = [
    "#00468B", "#ED0000", "#42B540", "#0099B4", "#925E9F",
    "#FDAF91", "#AD002A", "#ADB6B6"
]

NEJM_COLORS = [
    "#BC3C29", "#0072B5", "#E18727", "#20854E", "#7876B1",
    "#6F99AD", "#FFDC91", "#EE4C97"
]

JAMA_COLORS = [
    "#374E55", "#DF8F44", "#00A1D5", "#B24745", "#79AF97",
    "#6A6599", "#80796B"
]


# =============================================================================
# 生物信息学工具风格色板
# =============================================================================

GSEA_COLORS = [
    "#58ACFA", "#BC8F8F", "#FF6347", "#4682B4", "#9ACD32",
    "#DDA0DD", "#F0E68C", "#FF69B4"
]

Cytoscape_COLORS = [
    "#FF9900", "#66CC00", "#0099FF", "#FF0066", "#9900CC",
    "#00CC99", "#FFCC00", "#CC3300"
]

IGV_COLORS = [
    "#0000FF", "#00FF00", "#FF0000", "#00FFFF", "#FF00FF",
    "#FFFF00", "#FFA500", "#800080"
]

TBTOOLS_COLORS = [
    "#FF6B6B", "#4ECDC4", "#45B7D1", "#96CEB4", "#FFEAA7",
    "#DDA0DD", "#98D8C8", "#F7DC6F"
]

OMICSHARE_COLORS = [
    "#FF6B9D", "#C44569", "#F8B500", "#4ECDC4", "#556270",
    "#36D1DC", "#5AB9EA", "#8860D0"
]


# =============================================================================
# 中国风格色板
# =============================================================================

CHINA_STYLE_COLORS = [
    "#C23531", "#2F4554", "#61A0A8", "#D48265", "#91C7AE",
    "#749F83", "#CA8622", "#BDA29A", "#6E7074", "#546570"
]


# =============================================================================
# PALETTES 注册表 - 19组配色方案
# =============================================================================

PALETTES: Dict[str, List[str]] = {
    # 默认 (1组)
    "default": TOL_BRIGHT,
    
    # Paul Tol 系列 (8组)
    "tol_bright": TOL_BRIGHT,
    "tol_high_contrast": TOL_HIGH_CONTRAST,
    "tol_vibrant": TOL_VIBRANT,
    "tol_muted": TOL_MUTED,
    "tol_medium_contrast": TOL_MEDIUM_CONTRAST,
    "tol_light": TOL_LIGHT,
    "tol_sunset": TOL_SUNSET,
    "tol_burga": TOL_BURGA,
    
    # Okabe-Ito - 色盲友好标准 (1组)
    "okabe_ito": OKABE_ITO,
    
    # 科研期刊 (6组)
    "nature": NATURE_COLORS,
    "science": SCIENCE_COLORS,
    "cell": CELL_COLORS,
    "lancet": LANCET_COLORS,
    "nejm": NEJM_COLORS,
    "jama": JAMA_COLORS,
    
    # 生物信息学工具 (2组)
    "gsea": GSEA_COLORS,
    "omicshare": OMICSHARE_COLORS,
    
    # 中国风格 (1组)
    "china_style": CHINA_STYLE_COLORS,
}


class ColorConfig:
    """
    颜色配置类
    
    提供统一的颜色管理接口，所有图表颜色必须从此类获取。
    
    Usage:
        >>> config = ColorConfig()
        >>> colors = config.get_colors('nature', n=5)
        >>> go_colors = config.get_categorical_colors('go')
        >>> volcano_colors = config.get_volcano_colors()
    """
    
    def __init__(self):
        """初始化颜色配置"""
        self._palettes = PALETTES.copy()
    
    def get_available_palettes(self) -> List[str]:
        """
        获取所有可用的配色方案名称
        
        Returns:
            配色方案名称列表
        """
        return list(self._palettes.keys())
    
    def get_colors(self, palette_name: str = 'default', n: int = 8) -> List[str]:
        """
        获取指定配色方案的颜色列表
        
        Args:
            palette_name: 配色方案名称，默认为'default'
            n: 需要的颜色数量
            
        Returns:
            颜色列表（十六进制格式）
            
        Raises:
            ValueError: 如果配色方案不存在
        """
        if palette_name not in self._palettes:
            raise ValueError(f"未知的配色方案: {palette_name}。可用的配色方案: {list(self._palettes.keys())}")
        
        palette = self._palettes[palette_name]
        
        # 如果需要的颜色数量超过色板大小，循环使用
        colors = []
        for i in range(n):
            colors.append(palette[i % len(palette)])
        
        return colors
    
    def get_categorical_colors(self, category_type: str, palette: Optional[str] = None) -> Dict[str, str]:
        """
        获取分类颜色映射
        
        Args:
            category_type: 分类类型，支持 'go', 'kegg'
            palette: 色板名称，None则使用默认色板
            
        Returns:
            分类到颜色的映射字典
            
        Raises:
            ValueError: 如果分类类型不支持
        """
        if category_type.lower() == 'go':
            return self.get_go_category_colors(palette)
        elif category_type.lower() == 'kegg':
            return self.get_kegg_category_colors(palette)
        else:
            raise ValueError(f"不支持的分类类型: {category_type}。支持: 'go', 'kegg'")
    
    def get_go_category_colors(self, palette: Optional[str] = None) -> Dict[str, str]:
        """
        获取GO分类颜色映射 - 从指定配色方案动态生成
        
        GO三大分类：
        - biological_process: 生物过程
        - cellular_component: 细胞组分
        - molecular_function: 分子功能
        
        Args:
            palette: 色板名称，None则使用默认色板
            
        Returns:
            GO三大分类的颜色映射字典
        """
        colors = self.get_colors(palette or 'default', n=3)
        return {
            "biological_process": colors[0],
            "cellular_component": colors[1],
            "molecular_function": colors[2],
        }
    
    def get_kegg_category_colors(self, palette: Optional[str] = None) -> Dict[str, str]:
        """
        获取KEGG分类颜色映射 - 从指定配色方案动态生成
        
        KEGG六大分类：
        - Genetic Information Processing: 遗传信息处理
        - Human Diseases: 人类疾病
        - Metabolism: 代谢
        - Cellular Processes: 细胞过程
        - Organismal Systems: 生物体系统
        - Environmental Information Processing: 环境信息处理
        
        Args:
            palette: 色板名称，None则使用默认色板
            
        Returns:
            KEGG六大分类的颜色映射字典
        """
        colors = self.get_colors(palette or 'default', n=6)
        return {
            "Genetic_Information_Processing": colors[0],
            "Human_Diseases": colors[1],
            "Metabolism": colors[2],
            "Cellular_Processes": colors[3],
            "Organismal_Systems": colors[4],
            "Environmental_Information_Processing": colors[5],
        }
    
    def get_volcano_colors(self) -> Dict[str, str]:
        """
        获取火山图颜色配置
        
        Returns:
            火山图颜色映射字典，包含 'up', 'down', 'ns' 三个键
        """
        return VOLCANO_COLORS.copy()
    
    def get_palette_colors(self, palette_name: str) -> List[str]:
        """
        获取指定配色方案的完整颜色列表
        
        Args:
            palette_name: 配色方案名称
            
        Returns:
            完整颜色列表
        """
        if palette_name not in self._palettes:
            raise ValueError(f"未知的配色方案: {palette_name}")
        return self._palettes[palette_name].copy()


# 全局颜色配置实例（单例模式）
_color_config: Optional[ColorConfig] = None


def get_color_config() -> ColorConfig:
    """
    获取全局颜色配置实例
    
    Returns:
        ColorConfig 实例
    """
    global _color_config
    if _color_config is None:
        _color_config = ColorConfig()
    return _color_config
