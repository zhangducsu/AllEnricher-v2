"""Resolve species identifiers and names from the bundled lookup data."""

import warnings
warnings.warn(
    "species_lookup is deprecated, use species_registry instead",
    DeprecationWarning,
    stacklevel=2
)

from dataclasses import dataclass
from typing import Optional, Dict, List

from allenricher.core.config import SPECIES_CONFIGS, SpeciesConfig


@dataclass
class SpeciesInfo:
    """Store normalized identifiers and names for one species."""
    kegg_code: str
    latin_name: str
    taxonomy_id: int
    display_name: str = ""


# Inline TaxID Map (for inverted search)
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
    """Resolve species by NCBI TaxID, KEGG code, or scientific name."""
    
    def __init__(self, auto_load: bool = True):
        """Initialization species searcher
        
        Args:
auto_load: Whether to automatically load species data from the network (Default True)
        """
        self.auto_load = auto_load
        self.loaded = False
        self.species_data: Dict[str, SpeciesInfo] = {}
        
        if auto_load:
            self._load_builtin_species()
            self.loaded = True
    
    def _load_builtin_species(self) -> None:
        """Load the bundled species lookup records."""
        for kegg_code, config in SPECIES_CONFIGS.items():
            self.species_data[kegg_code] = SpeciesInfo(
                kegg_code=config.kegg_code,
                latin_name=config.name,
                taxonomy_id=config.taxonomy_id,
                display_name=config.display_name
            )
    
    def lookup_by_kegg_code(self, kegg_code: str) -> Optional[SpeciesInfo]:
        """Resolve one species by KEGG organism code."""
        return self.species_data.get(kegg_code)
    
    def lookup_by_latin_name(self, latin_name: str) -> Optional[SpeciesInfo]:
        """Resolve one species by scientific name."""
        for info in self.species_data.values():
            if info.latin_name.lower() == latin_name.lower():
                return info
        return None
    
    def lookup_by_taxid(self, taxonomy_id: int) -> Optional[SpeciesInfo]:
        """Resolve one species by NCBI TaxID."""
        kegg_code = BUILTIN_TAXID_MAP.get(taxonomy_id)
        if kegg_code:
            return self.species_data.get(kegg_code)
        return None
    
    def get_all_species(self) -> List[SpeciesInfo]:
        """Return all loaded species records."""
        return list(self.species_data.values())
    
    def get_species_count(self) -> int:
        """Return the number of loaded species records."""
        return len(self.species_data)