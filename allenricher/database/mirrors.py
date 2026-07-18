"""Define ordered mirrors for supported upstream data sources."""
from dataclasses import dataclass, field
from typing import List


@dataclass
class MirrorSource:
    """Describe one downloadable mirror and its priority."""
    name: str
    base_url: str
    priority: int
    region: str
    enabled: bool = True


# ============================
# NCBI gene_info / gene2go mirror source
# ============================
NCBI_MIRRORS: List[MirrorSource] = [
    MirrorSource(
        "ncbi-official",
        "https://ftp.ncbi.nlm.nih.gov/gene/DATA/",
        1, "US"
    ),
    # Add domestic mirror to verify usability
]

# ============================
# GO OBO File Mirror Source
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
# Reactome Mirror Source
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
# Jensen Lab DO data source (no public mirror)
# ============================
JENSEN_SOURCES: List[str] = [
    "http://download.jensenlab.org/human_disease_textmining_filtered.tsv",
    "http://download.jensenlab.org/human_disease_knowledge_filtered.tsv",
    "http://download.jensenlab.org/human_disease_experiments_filtered.tsv",
]

# ============================
# WikiPathways Mirror Source
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
    """Return mirrors for a data source in priority order."""
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
