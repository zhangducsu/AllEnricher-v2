"""
Configuration management for AllEnricher v2.0

AllEnricher v2.0 配置管理模块
==============================

本模块是 AllEnricher v2.0 的核心配置管理模块，负责定义和管理整个富集分析流程中
所需的所有配置参数。主要功能包括：

1. 定义支持的富集分析方法（Fisher精确检验、超几何检验、GSEA、ssGSEA、GSVA）
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
- AIBackendConfig: AI后端配置（API密钥、模型等）
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
        HYPERGEOMETRIC: 超几何检验（ORA默认方法）
        GSEA: 基因集富集分析（Gene Set Enrichment Analysis），适用于排序基因列表
        SSGSEA: 单样本GSEA（Single Sample GSEA），适用于单样本的基因集富集分析
        GSVA: 基因集变异分析（Gene Set Variation Analysis），适用于样本级别的基因集活性评估
    """
    HYPERGEOMETRIC = "hypergeometric"    # 超几何检验（ORA默认方法）
    GSEA = "gsea"                        # 基因集富集分析（Gene Set Enrichment Analysis）
    SSGSEA = "ssgsea"                    # 单样本GSEA
    GSVA = "gsva"                        # 基因集变异分析（Gene Set Variation Analysis）


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
class AIBackendConfig:
    """AI后端配置

    存储单个AI后端的配置信息，包括API密钥、模型名称等。
    支持OpenAI、Claude、DeepSeek、GLM、MiniMax、Ollama等多种后端。

    Attributes:
        api_key: API密钥，从YAML配置文件读取或环境变量获取
        model: 模型名称，如gpt-4、claude-3-opus、deepseek-chat等
        base_url: 自定义API基础URL（可选，用于兼容OpenAI格式的第三方服务）
        group_id: MiniMax等后端需要的额外参数
        enabled: 是否启用该后端
    """
    api_key: Optional[str] = None
    model: Optional[str] = None
    base_url: Optional[str] = None
    group_id: Optional[str] = None
    enabled: bool = True


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
    - AI后端设置：各AI服务的API密钥和模型配置（从YAML读取）
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
    method: str = "hypergeometric"         # 富集分析方法，默认超几何检验（ORA）
    correction: str = "BH"                 # 多重检验校正方法，默认为BH（Benjamini-Hochberg）
    pvalue_cutoff: float = 0.05            # p值显著性阈值，默认为0.05
    qvalue_cutoff: float = 0.05            # 校正后q值（FDR）显著性阈值，默认为0.05
    min_genes: int = 2                     # 基因集最小基因数，少于该值的基因集将被过滤
    max_genes: float = float('inf')        # 基因集最大基因数，默认无限制（设为无穷大）；设为具体数值可过滤过于宽泛的条目
    output_all: bool = True                # 是否输出全部条目（默认True，与v1一致）；设为False则仅输出满足p/q阈值的显著条目

    # GSEA specific settings
    # GSEA专用设置（仅在使用GSEA或ssGSEA方法时生效）
    gsea_permutations: int = 1000          # GSEA排列检验次数，次数越多结果越精确但越慢
    gsea_min_size: int = 10                # GSEA基因集最小大小，参考clusterProfiler默认值（minGSSize=10）
    gsea_max_size: int = 500               # GSEA基因集最大大小，参考clusterProfiler默认值（maxGSSize=500）

    # GSVA specific settings
    # GSVA专用设置（仅在使用GSVA方法时生效）
    gsva_method: str = "gsva"              # GSVA 方法变体: gsva（默认）/ plage / zscore
    gsva_kcdf: str = "Gaussian"            # 核密度估计核函数: Gaussian（默认）/ Poisson
    gsva_tau: float = 1.0                  # 核密度带宽参数，默认1.0（仅Gaussian核时生效）

    # Visualization settings
    # 可视化设置
    plot_formats: List[str] = field(default_factory=lambda: ["pdf", "png"])  # 图表输出格式列表，支持pdf和png
    plot_dpi: int = 300                    # 图表分辨率（DPI），默认为300
    top_terms: int = 20                    # 可视化展示的top富集条目数量
    plot_width: int = 10                   # 图表宽度（英寸）
    plot_height: int = 8                   # 图表高度（英寸）
    plot_style: str = 'nature'             # 图表风格: nature, science, colorblind, presentation, omicshare
    plot_palette: Optional[str] = None     # 自定义配色方案名称（可选）

    # Performance settings
    # 性能设置
    n_jobs: int = 1                        # 并行任务数，1表示不使用并行，-1表示使用所有CPU核心

    # Report settings
    # 报告设置
    generate_report: bool = True           # 是否自动生成分析报告
    report_format: str = "html"            # 报告格式，支持html
    ai_interpretation: bool = False        # 是否启用AI解读功能，需要配置AI API密钥
    ai_backend: str = "openai"             # 默认使用的AI后端，默认为openai
    ai_api_key: Optional[str] = None       # 兼容旧配置：全局AI API密钥（建议使用ai_backends）
    ai_model: Optional[str] = None         # 兼容旧配置：AI模型名称（建议使用ai_backends）

    # AI Backends configuration
    # AI后端多密钥配置，从YAML读取
    ai_backends: Dict[str, AIBackendConfig] = field(default_factory=dict)
    # 格式: {
    #   "openai": {"api_key": "sk-xxx", "model": "gpt-4"},
    #   "claude": {"api_key": "sk-ant-xxx", "model": "claude-3-opus"},
    #   "deepseek": {"api_key": "ds-xxx"},
    #   "glm": {"api_key": "glm-xxx"},
    #   "minimax": {"api_key": "mm-xxx", "group_id": "xxx"},
    #   "ollama": {"model": "llama2", "base_url": "http://localhost:11434"}
    # }

    # API settings
    # API服务设置（用于启动Web API服务）
    api_host: str = "0.0.0.0"              # API服务监听地址，0.0.0.0表示监听所有网络接口
    api_port: int = 8000                   # API服务监听端口
    api_debug: bool = False                # 是否启用API调试模式

    # Database settings
    # 数据库设置
    database_dir: str = "./database"       # 本地数据库文件存储目录
    auto_update: bool = False              # 是否自动更新本地数据库
    use_version: Optional[str] = None      # 指定使用的数据库版本（如 v20260515），默认使用最新版本

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
        path = Path(config_file)

        if not path.exists():
            raise FileNotFoundError(f"Configuration file not found: {config_file}")

        with open(path, 'r', encoding='utf-8') as f:
            if path.suffix in ['.yaml', '.yml']:
                data = yaml.safe_load(f)      # 使用PyYAML安全加载YAML文件
            elif path.suffix == '.json':
                data = json.load(f)            # 使用json模块加载JSON文件
            else:
                raise ValueError(f"Unsupported config file format: {path.suffix}")

        # 如果配置中有 ai_backends，需要将字典转换为 AIBackendConfig 对象
        if 'ai_backends' in data and isinstance(data['ai_backends'], dict):
            ai_backends_config = {}
            for backend_name, backend_data in data['ai_backends'].items():
                if isinstance(backend_data, dict):
                    ai_backends_config[backend_name] = AIBackendConfig(**backend_data)
                else:
                    ai_backends_config[backend_name] = AIBackendConfig(api_key=backend_data)
            data['ai_backends'] = ai_backends_config

        return cls(**data)

    def to_file(self, config_file: str) -> None:
        """Save configuration to file

        将当前配置保存到YAML或JSON文件中。

        Args:
            config_file: 输出配置文件的路径，支持 .yaml/.yml（YAML格式）和 .json（JSON格式）

        Raises:
            ValueError: 当配置文件格式不支持时抛出
        """
        path = Path(config_file)

        # Convert dataclass to dict
        # 将dataclass转换为字典，枚举值转换为其实际值
        data = {}
        for k, v in self.__dict__.items():
            if isinstance(v, Enum):
                data[k] = v.value
            elif k == 'ai_backends' and isinstance(v, dict):
                # 转换 AIBackendConfig 为字典
                ai_backends_dict = {}
                for name, backend in v.items():
                    ai_backends_dict[name] = {
                        'api_key': backend.api_key,
                        'model': backend.model,
                        'base_url': backend.base_url,
                        'group_id': backend.group_id,
                        'enabled': backend.enabled,
                    }
                data[k] = ai_backends_dict
            else:
                data[k] = v

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

    def get_ai_backend_config(self, backend_name: str) -> Optional[AIBackendConfig]:
        """Get AI backend configuration by name

        获取指定AI后端的配置对象。

        Args:
            backend_name: AI后端名称（openai, claude, deepseek, glm, minimax, ollama, mock）

        Returns:
            AIBackendConfig: 配置对象，如果未配置则返回None
        """
        if backend_name in self.ai_backends:
            return self.ai_backends[backend_name]
        return None

    def get_ai_api_key(self, backend_name: str) -> Optional[str]:
        """Get API key for specified AI backend

        获取指定AI后端的API密钥，优先级：
        1. YAML配置中的 ai_backends[backend_name].api_key
        2. 全局 ai_api_key（兼容旧配置）
        3. 环境变量（{BACKEND_NAME}_API_KEY）

        Args:
            backend_name: AI后端名称

        Returns:
            str: API密钥，如果都未配置则返回None
        """
        # 优先级 1: YAML ai_backends 中的配置
        backend_config = self.get_ai_backend_config(backend_name)
        if backend_config and backend_config.api_key:
            return backend_config.api_key

        # 优先级 2: 全局 ai_api_key（兼容旧配置）
        if self.ai_api_key:
            return self.ai_api_key

        # 优先级 3: 环境变量
        env_var = f"{backend_name.upper()}_API_KEY"
        return os.getenv(env_var)

    def get_ai_model(self, backend_name: str, default: str) -> str:
        """Get model name for specified AI backend

        获取指定AI后端的模型名称，优先级：
        1. YAML配置中的 ai_backends[backend_name].model
        2. 全局 ai_model（兼容旧配置）
        3. 默认值

        Args:
            backend_name: AI后端名称
            default: 默认模型名称（各后端硬编码默认值）

        Returns:
            str: 模型名称
        """
        backend_config = self.get_ai_backend_config(backend_name)
        if backend_config and backend_config.model:
            return backend_config.model

        if self.ai_model:
            return self.ai_model

        return default

    def get_ai_base_url(self, backend_name: str) -> Optional[str]:
        """Get custom base URL for specified AI backend

        获取指定AI后端的自定义API基础URL。

        Args:
            backend_name: AI后端名称

        Returns:
            str: 自定义base URL，如果未配置则返回None
        """
        backend_config = self.get_ai_backend_config(backend_name)
        if backend_config:
            return backend_config.base_url
        return None

    def get_ai_group_id(self, backend_name: str) -> Optional[str]:
        """Get group ID for specified AI backend (for MiniMax)

        获取指定AI后端的group ID（MiniMax需要）。

        Args:
            backend_name: AI后端名称

        Returns:
            str: group ID，如果未配置则返回None
        """
        backend_config = self.get_ai_backend_config(backend_name)
        if backend_config:
            return backend_config.group_id
        return None

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
        species_valid = False
        
        # 1. 先检查 SPECIES_CONFIGS（保持向后兼容）
        if self.species in SPECIES_CONFIGS:
            species_valid = True
        else:
            # 2. 尝试从注册表查询（延迟导入避免循环依赖）
            try:
                from ..database.species_registry import SpeciesRegistry
                registry = SpeciesRegistry.load_default()
                
                # 尝试按 kegg_code 查询
                entry = registry.query_by_kegg_code(self.species)
                if entry:
                    species_valid = True
                else:
                    # 尝试按 taxid 查询（如果 species 是数字）
                    if self.species.isdigit():
                        entry = registry.query_by_taxid(int(self.species))
                        if entry:
                            species_valid = True
            except Exception:
                pass  # 注册表查询失败，静默回退到原有逻辑
        
        if not species_valid:
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
method: "hypergeometric"  # hypergeometric, gsea, ssgsea, gsva / 富集分析方法
correction: "BH"  # BH, BY, bonferroni, holm, none / 多重检验校正方法
pvalue_cutoff: 0.05  # P-value cutoff / p值显著性阈值
qvalue_cutoff: 0.05  # Q-value (FDR) cutoff / q值（FDR）显著性阈值
min_genes: 2   # Minimum genes in gene set / 基因集最小基因数
max_genes: .inf  # Maximum genes in gene set (.inf = unlimited) / 基因集最大基因数（.inf表示无限制）
output_all: true  # Output all terms (true=match v1, false=only significant) / 输出全部条目（true=与v1一致，false=仅显著）

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
ai_backend: "openai"     # Default AI backend to use / 默认使用的AI后端

# AI Backends configuration
# AI后端配置（从YAML读取多密钥配置）
# 每个后端可以配置 api_key、model、base_url、group_id 等参数
ai_backends:
  openai:
    api_key: null        # OpenAI API key / 获取方式: https://platform.openai.com/
    model: "gpt-4"       # Model name / 模型名称: gpt-4, gpt-3.5-turbo 等
    enabled: true

  claude:
    api_key: null        # Anthropic API key / 获取方式: https://console.anthropic.com/
    model: "claude-3-opus-20240229"  # Model name
    enabled: true

  deepseek:
    api_key: null        # DeepSeek API key / 获取方式: https://platform.deepseek.com/
    model: "deepseek-chat"  # Model name
    enabled: true

  glm:
    api_key: null        # 智谱AI GLM API key / 获取方式: https://open.bigmodel.cn/
    model: "glm-4"       # Model name
    enabled: true

  minimax:
    api_key: null        # MiniMax API key / 获取方式: https://www.minimaxi.com/
    group_id: null       # MiniMax Group ID / 需要从控制台获取
    model: "abab6.5s-chat"  # Model name
    enabled: true

  ollama:
    model: "llama2"      # Model name / 需要先执行 ollama pull <model>
    base_url: "http://localhost:11434"  # Ollama service address
    enabled: true

# Database settings
# 数据库设置
database_dir: "./database"  # Database directory / 数据库存储目录
auto_update: false          # Auto update databases / 是否自动更新数据库
"""
