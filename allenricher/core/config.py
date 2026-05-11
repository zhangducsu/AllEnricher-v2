"""
Configuration management for AllEnricher v2.0

AllEnricher v2.0 配置管理模块
==============================

本模块是 AllEnricher v2.0 的核心配置管理模块，负责定义和管理整个富集分析流程中
所需的所有配置参数。主要功能包括：

1. 定义支持的富集分析方法（Fisher精确检验、超几何检验、GSEA、ssGSEA）
2. 定义支持的多重检验校正方法（BH、BY、Bonferroni、Holm等）
3. 定义支持的数据库类型（GO、KEGG、Reactome、WikiPathways等）
4. 提供物种配置信息（物种名称、KEGG代码、分类学ID等）
5. 提供主配置类 Config，支持从YAML/JSON文件加载和保存配置
6. 提供配置参数校验功能，确保参数合法性

主要类说明：
- EnrichmentMethod: 富集分析方法枚举
- CorrectionMethod: 多重检验校正方法枚举
- DatabaseType: 数据库类型枚举
- SpeciesConfig: 物种配置数据类
- Config: 主配置类，包含所有可配置参数
"""

import os
import json
import yaml
from pathlib import Path
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Any
from enum import Enum


class EnrichmentMethod(Enum):
    """Supported enrichment methods

    支持的富集分析方法枚举类。

    定义了 AllEnricher 支持的所有富集分析算法，包括基于过表示分析（ORA）
    的方法和基于基因集富集分析（GSEA）的方法。

    Attributes:
        FISHER: Fisher精确检验，适用于基因列表的过表示分析
        HYPERGEOMETRIC: 超几何检验，另一种过表示分析方法
        GSEA: 基因集富集分析（Gene Set Enrichment Analysis），适用于排序基因列表
        SSGSEA: 单样本GSEA（Single Sample GSEA），适用于单样本的基因集富集分析
    """
    FISHER = "fisher"                    # Fisher精确检验
    HYPERGEOMETRIC = "hypergeometric"    # 超几何检验
    GSEA = "gsea"                        # 基因集富集分析（Gene Set Enrichment Analysis）
    SSGSEA = "ssgsea"                    # 单样本GSEA（Single Sample GSEA）


class CorrectionMethod(Enum):
    """Multiple testing correction methods

    多重检验校正方法枚举类。

    在富集分析中，由于同时检验大量基因集，需要进行多重检验校正
    以控制假阳性率。本枚举定义了支持的校正方法。

    Attributes:
        BH: Benjamini-Hochberg方法，控制错误发现率（FDR），最常用
        BY: Benjamini-Yekutieli方法，比BH更保守，适用于依赖性检验
        BONFERRONI: Bonferroni校正，最保守的方法，控制族错误率（FWER）
        HOLM: Holm逐步校正方法，比Bonferroni更有力
        NONE: 不进行多重检验校正
    """
    BH = "BH"                # Benjamini-Hochberg校正，控制错误发现率（FDR）
    BY = "BY"                # Benjamini-Yekutieli校正，适用于依赖性检验
    BONFERRONI = "bonferroni"  # Bonferroni校正，最保守的FWER控制方法
    HOLM = "holm"            # Holm逐步校正方法
    NONE = "none"            # 不进行多重检验校正


class DatabaseType(Enum):
    """Supported databases

    支持的数据库类型枚举类。

    定义了 AllEnricher 支持的所有功能注释数据库，这些数据库用于
    提供基因集/通路信息以进行富集分析。

    Attributes:
        GO: Gene Ontology基因本体数据库，包含生物过程（BP）、分子功能（MF）、细胞组分（CC）
        KEGG: KEGG通路数据库，提供代谢和信号通路信息
        REACTOME: Reactome通路数据库，专注于人类生物学通路
        WIKIPATHWAYS: WikiPathways通路数据库，社区维护的开源通路数据库
        MSIGDB: MSigDB分子特征数据库，包含多种基因集集合
        DO: Disease Ontology疾病本体数据库，用于疾病相关富集分析
        DISGENET: DisGeNET数据库，基因-疾病关联数据库
    """
    GO = "GO"                      # Gene Ontology基因本体数据库
    KEGG = "KEGG"                  # KEGG通路数据库
    REACTOME = "Reactome"          # Reactome通路数据库
    WIKIPATHWAYS = "WikiPathways"  # WikiPathways开源通路数据库
    MSIGDB = "MSigDB"              # MSigDB分子特征数据库
    DO = "DO"                      # Disease Ontology疾病本体数据库
    DISGENET = "DisGeNET"          # DisGeNET基因-疾病关联数据库


@dataclass
class SpeciesConfig:
    """Species configuration

    物种配置数据类。

    存储特定物种的配置信息，包括物种名称、KEGG物种代码、
    NCBI分类学ID和显示名称。用于在富集分析中正确映射基因和通路。

    Attributes:
        name: 物种拉丁学名（如 "Homo sapiens"）
        kegg_code: KEGG数据库中使用的物种代码（如 "hsa"）
        taxonomy_id: NCBI分类学数据库中的物种ID（如人类为9606）
        display_name: 物种的常用显示名称（如 "Human"），默认与name相同
    """
    name: str              # 物种拉丁学名
    kegg_code: str         # KEGG物种代码
    taxonomy_id: int       # NCBI分类学ID
    display_name: str = ""  # 显示名称，默认为空字符串
    
    def __post_init__(self):
        """初始化后处理：如果未指定显示名称，则使用物种拉丁学名作为显示名称"""
        if not self.display_name:
            self.display_name = self.name  # 默认使用拉丁学名作为显示名称


# Common species configurations
# 常用物种配置字典：键为KEGG物种代码，值为SpeciesConfig对象
# 涵盖了生物医学研究中最常用的模式生物和实验动物
SPECIES_CONFIGS: Dict[str, SpeciesConfig] = {
    "hsa": SpeciesConfig("Homo sapiens", "hsa", 9606, "Human"),          # 人类
    "mmu": SpeciesConfig("Mus musculus", "mmu", 10090, "Mouse"),          # 小鼠
    "rno": SpeciesConfig("Rattus norvegicus", "rno", 10116, "Rat"),       # 大鼠
    "dre": SpeciesConfig("Danio rerio", "dre", 7955, "Zebrafish"),        # 斑马鱼
    "dme": SpeciesConfig("Drosophila melanogaster", "dme", 7227, "Fruit fly"),  # 果蝇
    "cel": SpeciesConfig("Caenorhabditis elegans", "cel", 6239, "C. elegans"),  # 秀丽隐杆线虫
    "ssc": SpeciesConfig("Sus scrofa", "ssc", 9823, "Pig"),               # 猪
    "bta": SpeciesConfig("Bos taurus", "bta", 9913, "Cow"),               # 牛
    "gga": SpeciesConfig("Gallus gallus", "gga", 9031, "Chicken"),         # 鸡
    "xtr": SpeciesConfig("Xenopus tropicalis", "xtr", 8364, "Xenopus"),   # 爪蟾
    "cfa": SpeciesConfig("Canis familiaris", "cfa", 9615, "Dog"),           # 狗
    "ddi": SpeciesConfig("Dictyostelium discoideum", "ddi", 44689, "Dictyostelium"),  # 盘基网柄菌
    "mtu": SpeciesConfig("Mycobacterium tuberculosis", "mtu", 1772, "M. tuberculosis"),  # 结核分枝杆菌
    "pfa": SpeciesConfig("Plasmodium falciparum", "pfa", 5833, "P. falciparum"),  # 恶性疟原虫
    "sce": SpeciesConfig("Saccharomyces cerevisiae", "sce", 4932, "S. cerevisiae"),  # 酿酒酵母
    "spo": SpeciesConfig("Schizosaccharomyces pombe", "spo", 4896, "S. pombe"),  # 裂殖酵母
}


@dataclass
class Config:
    """Main configuration class for AllEnricher

    AllEnricher 主配置类。

    本类是 AllEnricher v2.0 的核心配置类，使用 dataclass 实现，
    包含了富集分析流程中所有可配置的参数。支持从YAML/JSON文件加载配置，
    也可以直接通过构造函数或属性赋值来设置参数。

    配置参数分为以下几大类：
    - 输入/输出设置：指定输入基因列表文件、输出目录、背景基因列表等
    - 物种设置：指定分析的目标物种
    - 分析设置：指定数据库、分析方法、校正方法、阈值等
    - GSEA专用设置：GSEA方法的特定参数（排列次数、基因集大小范围等）
    - 可视化设置：图表格式、分辨率、尺寸等
    - 性能设置：并行任务数、缓存配置等
    - 报告设置：报告生成、AI解读等
    - API设置：Web服务的相关配置
    - 数据库设置：本地数据库存储和更新配置
    """
    
    # Input/Output settings
    # 输入/输出设置
    input_file: Optional[str] = None       # 输入基因列表文件路径，支持TXT/CSV格式
    output_dir: str = "./results"          # 结果输出目录，默认为当前目录下的results文件夹
    background_file: Optional[str] = None  # 背景基因列表文件路径（可选），用于ORA分析
    
    # Species settings
    # 物种设置
    species: str = "hsa"                   # 目标物种的KEGG代码，默认为人类（hsa）
    
    # Analysis settings
    # 分析设置
    databases: List[str] = field(default_factory=lambda: ["GO", "KEGG"])  # 要使用的数据库列表，默认使用GO和KEGG
    method: str = "fisher"                 # 富集分析方法，默认为Fisher精确检验
    correction: str = "BH"                 # 多重检验校正方法，默认为BH（Benjamini-Hochberg）
    pvalue_cutoff: float = 0.05            # p值显著性阈值，默认为0.05
    qvalue_cutoff: float = 0.05            # 校正后q值（FDR）显著性阈值，默认为0.05
    min_genes: int = 2                     # 基因集最小基因数，少于该值的基因集将被过滤
    max_genes: float = float('inf')           # 基因集最大基因数，默认无限制（设为无穷大）；设为具体数值可过滤过于宽泛的条目
    
    # GSEA specific settings
    # GSEA专用设置（仅在使用GSEA或ssGSEA方法时生效）
    gsea_permutations: int = 1000          # GSEA排列检验次数，次数越多结果越精确但越慢
    gsea_min_size: int = 10                # GSEA基因集最小大小，参考clusterProfiler默认值（minGSSize=10）
    gsea_max_size: int = 500               # GSEA基因集最大大小，参考clusterProfiler默认值（maxGSSize=500）
    
    # Visualization settings
    # 可视化设置
    plot_formats: List[str] = field(default_factory=lambda: ["pdf", "png"])  # 图表输出格式列表，支持pdf和png
    plot_dpi: int = 300                    # 图表分辨率（DPI），默认为300
    top_terms: int = 20                    # 可视化展示的top富集条目数量
    plot_width: int = 10                   # 图表宽度（英寸）
    plot_height: int = 8                   # 图表高度（英寸）
    
    # Performance settings
    # 性能设置
    n_jobs: int = 1                        # 并行任务数，1表示不使用并行，-1表示使用所有CPU核心
    
    # Report settings
    # 报告设置
    generate_report: bool = True           # 是否自动生成分析报告
    report_format: str = "html"            # 报告格式，支持html
    ai_interpretation: bool = False        # 是否启用AI解读功能，需要配置AI API密钥
    ai_api_key: Optional[str] = None       # AI API密钥（如OpenAI API Key）
    ai_model: str = "gpt-4"                # AI模型名称，默认为gpt-4
    
    # API settings
    # API服务设置（用于启动Web API服务）
    api_host: str = "0.0.0.0"              # API服务监听地址，0.0.0.0表示监听所有网络接口
    api_port: int = 8000                   # API服务监听端口
    api_debug: bool = False                # 是否启用API调试模式
    
    # Database settings
    # 数据库设置
    database_dir: str = "./database"       # 本地数据库文件存储目录
    auto_update: bool = False              # 是否自动更新本地数据库
    
    @classmethod
    def from_file(cls, config_file: str) -> "Config":
        """Load configuration from YAML or JSON file

        从YAML或JSON配置文件加载配置参数，并创建Config实例。

        Args:
            config_file: 配置文件的路径，支持 .yaml/.yml（YAML格式）和 .json（JSON格式）

        Returns:
            Config: 根据配置文件内容创建的配置对象

        Raises:
            FileNotFoundError: 当配置文件不存在时抛出
            ValueError: 当配置文件格式不支持时抛出
        """
        path = Path(config_file)  # 将文件路径转换为Path对象
        
        if not path.exists():
            raise FileNotFoundError(f"Configuration file not found: {config_file}")
        
        with open(path, 'r', encoding='utf-8') as f:
            if path.suffix in ['.yaml', '.yml']:
                data = yaml.safe_load(f)      # 使用PyYAML安全加载YAML文件
            elif path.suffix == '.json':
                data = json.load(f)            # 使用json模块加载JSON文件
            else:
                raise ValueError(f"Unsupported config file format: {path.suffix}")
        
        return cls(**data)  # 将解析后的字典作为关键字参数创建Config实例
    
    def to_file(self, config_file: str) -> None:
        """Save configuration to file

        将当前配置保存到YAML或JSON文件中。

        Args:
            config_file: 输出配置文件的路径，支持 .yaml/.yml（YAML格式）和 .json（JSON格式）

        Raises:
            ValueError: 当配置文件格式不支持时抛出
        """
        path = Path(config_file)  # 将文件路径转换为Path对象
        
        # Convert dataclass to dict
        # 将dataclass转换为字典，枚举值转换为其实际值
        data = {
            k: v.value if isinstance(v, Enum) else v  # 如果值是枚举类型，取其.value属性
            for k, v in self.__dict__.items()          # 遍历实例的所有属性
        }
        
        with open(path, 'w', encoding='utf-8') as f:
            if path.suffix in ['.yaml', '.yml']:
                yaml.dump(data, f, default_flow_style=False)  # 以易读格式写入YAML文件
            elif path.suffix == '.json':
                json.dump(data, f, indent=2)                   # 以缩进格式写入JSON文件
            else:
                raise ValueError(f"Unsupported config file format: {path.suffix}")
    
    def get_species_config(self) -> SpeciesConfig:
        """Get species configuration

        根据当前配置中的物种代码获取对应的物种配置信息。

        Returns:
            SpeciesConfig: 与当前物种代码对应的物种配置对象

        Raises:
            ValueError: 当物种代码不在预定义的物种配置中时抛出，
                       错误信息中会列出所有可用的物种代码
        """
        if self.species not in SPECIES_CONFIGS:
            raise ValueError(f"Unknown species: {self.species}. "
                           f"Available: {list(SPECIES_CONFIGS.keys())}")
        return SPECIES_CONFIGS[self.species]  # 返回对应的物种配置对象
    
    def validate(self) -> List[str]:
        """Validate configuration and return list of errors

        验证当前配置参数的合法性，检查所有关键参数是否在有效范围内。

        验证内容包括：
        1. 物种代码是否在预定义的物种列表中
        2. 富集分析方法是否为支持的枚举值
        3. 多重检验校正方法是否为支持的枚举值
        4. 数据库类型是否为支持的枚举值
        5. p值和q值阈值是否在 (0, 1] 范围内

        Returns:
            List[str]: 错误信息列表，如果配置完全合法则返回空列表。
                      每个错误信息字符串描述了一个具体的配置问题。
        """
        errors = []  # 初始化错误列表
        
        # Validate input file
        # 验证输入文件是否存在（如果已指定）
        if self.input_file and not Path(self.input_file).exists():
            errors.append(f"Input file not found: {self.input_file}")
        
        # Validate background file
        # 验证背景基因文件是否存在（如果已指定）
        if self.background_file and not Path(self.background_file).exists():
            errors.append(f"Background file not found: {self.background_file}")
        
        # Validate species
        # 验证物种代码是否有效
        if self.species not in SPECIES_CONFIGS:
            errors.append(f"Unknown species: {self.species}")
        
        # Validate method
        # 验证富集分析方法是否有效
        valid_methods = [m.value for m in EnrichmentMethod]  # 获取所有有效的方法值
        if self.method not in valid_methods:
            errors.append(f"Invalid method: {self.method}. Valid: {valid_methods}")
        
        # Validate correction
        # 验证多重检验校正方法是否有效
        valid_corrections = [c.value for c in CorrectionMethod]  # 获取所有有效的校正方法值
        if self.correction not in valid_corrections:
            errors.append(f"Invalid correction: {self.correction}. Valid: {valid_corrections}")
        
        # Validate databases
        # 验证数据库列表中的每个数据库是否有效
        valid_databases = [d.value for d in DatabaseType]  # 获取所有有效的数据库值
        for db in self.databases:
            if db not in valid_databases:
                errors.append(f"Invalid database: {db}. Valid: {valid_databases}")
        
        # Validate cutoffs
        # 验证p值和q值阈值是否在有效范围 (0, 1] 内
        if not 0 < self.pvalue_cutoff <= 1:
            errors.append(f"pvalue_cutoff must be in (0, 1], got: {self.pvalue_cutoff}")
        if not 0 < self.qvalue_cutoff <= 1:
            errors.append(f"qvalue_cutoff must be in (0, 1], got: {self.qvalue_cutoff}")
        
        return errors  # 返回所有错误信息


# Default configuration file content
# 默认配置文件内容模板（YAML格式）
# 用户可以基于此模板创建自定义配置文件，修改需要的参数值
DEFAULT_CONFIG_YAML = """# AllEnricher v2.0 Configuration File
# =================================
# AllEnricher v2.0 配置文件模板

# Input/Output settings
# 输入/输出设置
input_file: null  # Path to gene list file / 基因列表文件路径
output_dir: "./results"  # Output directory / 结果输出目录
background_file: null  # Path to background gene list (optional) / 背景基因列表文件路径（可选）

# Species settings
# 物种设置
species: "hsa"  # KEGG species code / KEGG物种代码

# Analysis settings
# 分析设置
databases:
  - "GO"        # Gene Ontology / 基因本体数据库
  - "KEGG"      # KEGG pathways / KEGG通路数据库
method: "fisher"  # fisher, hypergeometric, gsea, ssgsea / 富集分析方法
correction: "BH"  # BH, BY, bonferroni, holm, none / 多重检验校正方法
pvalue_cutoff: 0.05  # P-value cutoff / p值显著性阈值
qvalue_cutoff: 0.05  # Q-value (FDR) cutoff / q值（FDR）显著性阈值
min_genes: 2   # Minimum genes in gene set / 基因集最小基因数
max_genes: .inf  # Maximum genes in gene set (.inf = unlimited) / 基因集最大基因数（.inf表示无限制）

# GSEA specific settings
# GSEA专用设置
gsea_permutations: 1000  # Number of permutations / 排列检验次数
gsea_min_size: 10        # Minimum gene set size / 基因集最小大小（参考clusterProfiler默认值）
gsea_max_size: 500       # Maximum gene set size / 基因集最大大小

# Visualization settings
# 可视化设置
plot_formats:
  - "pdf"  # PDF format / PDF格式
  - "png"  # PNG format / PNG格式
plot_dpi: 300    # DPI for raster images / 栅格图像分辨率
top_terms: 20    # Number of top terms to plot / 展示的top富集条目数
plot_width: 10   # Plot width in inches / 图表宽度（英寸）
plot_height: 8   # Plot height in inches / 图表高度（英寸）

# Performance settings
# 性能设置
n_jobs: 1  # Number of parallel jobs / 并行任务数

# Report settings
# 报告设置
generate_report: true   # Generate HTML report / 是否生成报告
report_format: "html"   # Report format / 报告格式
ai_interpretation: false  # Enable AI interpretation / 是否启用AI解读
ai_api_key: null        # AI API key / AI API密钥
ai_model: "gpt-4"       # AI model name / AI模型名称

# Database settings
# 数据库设置
database_dir: "./database"  # Database directory / 数据库存储目录
auto_update: false          # Auto update databases / 是否自动更新数据库
"""
