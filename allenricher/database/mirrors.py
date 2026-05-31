"""镜像源配置模块

定义各数据库的镜像源列表，支持按优先级排序和自动切换。
主源失败时自动切换到备用源，提高下载可靠性。
"""
from dataclasses import dataclass, field
from typing import List


@dataclass
class MirrorSource:
    """镜像源定义

    Attributes:
        name: 镜像名称，如 "ncbi-official"
        base_url: 基础 URL（末尾带 /）
        priority: 优先级，数字越小优先级越高
        region: 地区标识，如 "US", "EU", "CN"
        enabled: 是否启用
    """
    name: str
    base_url: str
    priority: int
    region: str
    enabled: bool = True


# ============================
# NCBI gene_info / gene2go 镜像源
# ============================
NCBI_MIRRORS: List[MirrorSource] = [
    MirrorSource(
        "ncbi-official",
        "https://ftp.ncbi.nlm.nih.gov/gene/DATA/",
        1, "US"
    ),
    MirrorSource(
        "ebi-ftp",
        "https://ftp.ebi.ac.uk/pub/databases/ncbi/gene/DATA/",
        2, "EU"
    ),
    # 国内镜像待验证可用性后添加
]

# ============================
# GO OBO 文件镜像源
# ============================
GO_MIRRORS: List[MirrorSource] = [
    MirrorSource(
        "purl-obo",
        "http://purl.obolibrary.org/obo/go/",
        1, "US"
    ),
    MirrorSource(
        "geneontology-releases",
        "http://release.geneontology.org/",
        2, "US"
    ),
]

# ============================
# Reactome 镜像源
# ============================
REACTOME_MIRRORS: List[MirrorSource] = [
    MirrorSource(
        "reactome-official",
        "https://reactome.org/download/current/",
        1, "US"
    ),
    MirrorSource(
        "reactome-cshl",
        "https://download.reactome.org/",
        2, "US"
    ),
    MirrorSource(
        "reactome-ebi",
        "https://ftp.ebi.ac.uk/pub/databases/reactome/",
        3, "EU"
    ),
]

# ============================
# Jensen Lab DO 数据源（无公共镜像）
# ============================
JENSEN_SOURCES: List[str] = [
    "http://download.jensenlab.org/human_disease_textmining_filtered.tsv",
    "http://download.jensenlab.org/human_disease_knowledge_filtered.tsv",
    "http://download.jensenlab.org/human_disease_experiments_filtered.tsv",
]

# ============================
# WikiPathways 镜像源
# ============================
WIKIPATHWAYS_MIRRORS: List[MirrorSource] = [
    MirrorSource(
        name="wikipathways-official",
        base_url="https://data.wikipathways.org/",
        priority=1,
        region="US",
    ),
]


def get_mirrors(db_type: str) -> List[MirrorSource]:
    """获取指定数据库类型的镜像源列表（按优先级排序）

    Args:
        db_type: 数据库类型标识，如 'ncbi', 'go', 'reactome'

    Returns:
        按优先级排序的已启用镜像源列表
    """
    mirrors_map = {
        'ncbi': NCBI_MIRRORS,
        'go': GO_MIRRORS,
        'reactome': REACTOME_MIRRORS,
        'wikipathways': WIKIPATHWAYS_MIRRORS,
    }
    mirrors = mirrors_map.get(db_type, [])
    return sorted(
        [m for m in mirrors if m.enabled],
        key=lambda x: x.priority
    )
