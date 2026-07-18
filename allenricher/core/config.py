"""Configuration models and validation for AllEnricher.

This module defines supported methods, correction procedures, database
metadata, common species, runtime settings, and YAML/JSON serialization."""

import os
import json
import yaml
from pathlib import Path
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Any
from enum import Enum


class EnrichmentMethod(Enum):
    """Enrichment methods exposed by the public analysis workflow."""
    HYPERGEOMETRIC = "hypergeometric"
    GSEA = "gsea"
    SSGSEA = "ssgsea"
    GSVA = "gsva"


class CorrectionMethod(Enum):
    """Supported multiple-testing correction procedures."""
    BH = "BH"
    BY = "BY"
    BONFERRONI = "bonferroni"
    HOLM = "holm"
    NONE = "none"


class DatabaseType(Enum):
    """Databases supported by local analysis workflows."""
    GO = "GO"
    KEGG = "KEGG"
    REACTOME = "Reactome"
    WIKIPATHWAYS = "WikiPathways"
    DO = "DO"
    DISGENET = "DisGeNET"
    TRRUST = "TRRUST"
    CHEA3 = "ChEA3"
    ANIMALTFDB = "AnimalTFDB"
    HTFTARGET = "hTFtarget"
    CUSTOM = "CUSTOM"


# Shared public database catalog used by the CLI, API, reports, and Web UI.
# Runtime availability is still verified against local database files.
DISGENET_SNAPSHOT_VERSION = "v20190612"

DATABASE_CATALOG = (
    {"name": "GO", "description": "Gene Ontology", "species": "all"},
    {"name": "KEGG", "description": "KEGG Pathways", "species": "all"},
    {
        "name": "Reactome",
        "description": "Reactome Pathways",
        "species": "model_organisms",
    },
    {"name": "WikiPathways", "description": "WikiPathways", "species": "all"},
    {"name": "DO", "description": "Disease Ontology", "species": "hsa"},
    {
        "name": "DisGeNET",
        "display_name": f"DisGeNET ({DISGENET_SNAPSHOT_VERSION})",
        "description": "Disease-gene associations (AllEnricher-v1 free snapshot)",
        "species": "hsa",
        "source_version": DISGENET_SNAPSHOT_VERSION,
        "theoretical_species_count": 1,
    },
    {
        "name": "TRRUST",
        "description": "TF-target regulation",
        "species": "hsa,mmu",
        "theoretical_species_count": 2,
    },
    {
        "name": "ChEA3",
        "description": "Human TF-target libraries",
        "species": "hsa",
        "theoretical_species_count": 1,
    },
    {
        "name": "AnimalTFDB",
        "description": "Animal transcription factors with ortholog-mapped targets",
        "species": "non_human_animals",
        "theoretical_species_count": 183,
    },
    {
        "name": "hTFtarget",
        "description": "Human tissue-specific TF targets",
        "species": "hsa",
        "theoretical_species_count": 1,
    },
    {"name": "CUSTOM", "description": "User-supplied gene-set database", "species": "all"},
)


def database_catalog_entry(name: str) -> Dict[str, Any]:
    """Return shared display/source metadata for a database name."""
    target = str(name).casefold()
    return next(
        (dict(item) for item in DATABASE_CATALOG if item["name"].casefold() == target),
        {},
    )


def database_display_name(name: str) -> str:
    """Return a user-facing database label without changing its internal key."""
    item = database_catalog_entry(name)
    return str(item.get("display_name") or item.get("name") or name)


@dataclass
class SpeciesConfig:
    """Stable identifiers and display metadata for a species."""
    name: str
    kegg_code: str
    taxonomy_id: int
    display_name: str = ""

    def __post_init__(self):
        """Use the scientific name when no display name is supplied."""
        if not self.display_name:
            self.display_name = self.name


# Built-in shortcuts for commonly analyzed species. Taxonomy ID remains the
# authoritative identifier; the KEGG code is a convenient CLI alias.
SPECIES_CONFIGS: Dict[str, SpeciesConfig] = {
    "hsa": SpeciesConfig("Homo sapiens", "hsa", 9606, "Human"),
    "mmu": SpeciesConfig("Mus musculus", "mmu", 10090, "Mouse"),
    "rno": SpeciesConfig("Rattus norvegicus", "rno", 10116, "Rat"),
    "dre": SpeciesConfig("Danio rerio", "dre", 7955, "Zebrafish"),
    "dme": SpeciesConfig("Drosophila melanogaster", "dme", 7227, "Fruit fly"),
    "cel": SpeciesConfig("Caenorhabditis elegans", "cel", 6239, "C. elegans"),
    "ssc": SpeciesConfig("Sus scrofa", "ssc", 9823, "Pig"),
    "bta": SpeciesConfig("Bos taurus", "bta", 9913, "Cow"),
    "gga": SpeciesConfig("Gallus gallus", "gga", 9031, "Chicken"),
    "xtr": SpeciesConfig("Xenopus tropicalis", "xtr", 8364, "Xenopus"),
    "cfa": SpeciesConfig("Canis familiaris", "cfa", 9615, "Dog"),
    "ddi": SpeciesConfig("Dictyostelium discoideum", "ddi", 44689, "Dictyostelium"),
    "mtu": SpeciesConfig("Mycobacterium tuberculosis", "mtu", 1772, "M. tuberculosis"),
    "pfa": SpeciesConfig("Plasmodium falciparum", "pfa", 5833, "P. falciparum"),
    "sce": SpeciesConfig("Saccharomyces cerevisiae", "sce", 4932, "S. cerevisiae"),
    "spo": SpeciesConfig("Schizosaccharomyces pombe", "spo", 4896, "S. pombe"),
}


@dataclass
class AIBackendConfig:
    """Connection settings for one optional AI interpretation backend."""
    api_key: Optional[str] = None
    model: Optional[str] = None
    base_url: Optional[str] = None
    group_id: Optional[str] = None
    enabled: bool = True


@dataclass
class Config:
    """Validated configuration for an AllEnricher analysis run."""

    # Input and output
    input_file: Optional[str] = None
    output_dir: str = "./results"
    background_file: Optional[str] = None

    # Species
    species: str = "hsa"

    # Analysis
    databases: List[str] = field(default_factory=lambda: ["GO", "KEGG"])
    method: str = "hypergeometric"
    correction: str = "BH"
    pvalue_cutoff: float = 0.05
    qvalue_cutoff: float = 0.05
    min_genes: int = 3
    max_genes: float = float('inf')
    output_all: bool = True

    # GSEA and activity-method size limits
    gsea_permutations: int = 1000
    gsea_min_size: Optional[int] = None
    gsea_max_size: Optional[int] = None

    # GSVA
    gsva_method: str = "gsva"
    gsva_kcdf: str = "Gaussian"
    gsva_tau: float = 1.0

    # Visualization
    plot_formats: List[str] = field(default_factory=lambda: ["pdf", "png"])
    plot_dpi: int = 300
    top_terms: int = 20
    activity_heatmap_top_n: int = 40
    plot_width: Optional[float] = None
    plot_height: Optional[float] = None
    plot_style: str = 'nature'
    plot_palette: Optional[str] = None
    categorical_palette: Optional[str] = None
    sequential_palette: Optional[str] = None
    diverging_palette: Optional[str] = None

    # Performance
    n_jobs: int = 1

    # Reporting and optional AI interpretation
    generate_report: bool = True
    report_format: str = "html"
    ai_interpretation: bool = False
    ai_backend: str = "openai"
    ai_api_key: Optional[str] = None
    ai_model: Optional[str] = None

    # Per-backend AI settings loaded from YAML.
    ai_backends: Dict[str, AIBackendConfig] = field(default_factory=dict)
    # Format: {
    #   "openai": {"api_key": "sk-xxx", "model": "gpt-4"},
    #   "claude": {"api_key": "sk-ant-xxx", "model": "claude-3-opus"},
    #   "deepseek": {"api_key": "ds-xxx"},
    #   "glm": {"api_key": "glm-xxx"},
    #   "minimax": {"api_key": "mm-xxx", "group_id": "xxx"},
    #   "ollama": {"model": "llama2", "base_url": "http://localhost:11434"}
    # }

    # Local API service
    api_host: str = "127.0.0.1"
    api_port: int = 8000
    api_debug: bool = False

    # Local database storage
    database_dir: str = "./database"
    auto_update: bool = False
    use_version: Optional[str] = None

    @classmethod
    def from_file(cls, config_file: str) -> "Config":
        """Load configuration from a YAML or JSON file.
        
        Args:
            config_file: Path to a ``.yaml``, ``.yml``, or ``.json`` file.
        
        Returns:
            A populated configuration object.
        
        Raises:
            FileNotFoundError: If the file does not exist.
            ValueError: If the file type or contents are invalid."""
        path = Path(config_file)

        if not path.exists():
            raise FileNotFoundError(f"Configuration file not found: {config_file}")

        with open(path, 'r', encoding='utf-8') as f:
            if path.suffix in ['.yaml', '.yml']:
                data = yaml.safe_load(f)      # Load YAML files safely using PyYAML
            elif path.suffix == '.json':
                data = json.load(f)            # Load JSON files using json modules
            else:
                raise ValueError(f"Unsupported config file format: {path.suffix}")

        data = data or {}
        # Compatible with old version of the single format field, and is harmonized to plot_formats.
        if 'plot_format' in data and 'plot_formats' not in data:
            data['plot_formats'] = [data.pop('plot_format')]
        if isinstance(data.get('plot_formats'), str):
            data['plot_formats'] = [data['plot_formats']]

        # Normalize nested mappings into typed backend settings.
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
        """Write the configuration as YAML or JSON."""
        path = Path(config_file)

        # Convert dataclass to dict
        # Convert datacals to dictionary with enumerator values converted to actual values
        data = {}
        for k, v in self.__dict__.items():
            if isinstance(v, Enum):
                data[k] = v.value
            elif k == 'ai_backends' and isinstance(v, dict):
                # Serialize typed backend settings as plain mappings.
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
                yaml.dump(data, f, default_flow_style=False)  # Write YAML files in readable format
            elif path.suffix == '.json':
                json.dump(data, f, indent=2)                   # Write JSON files in indentation
            else:
                raise ValueError(f"Unsupported config file format: {path.suffix}")

    def get_species_config(self) -> SpeciesConfig:
        """Return metadata for the configured common species.
        
        Raises:
            ValueError: If the species code is not in the built-in common-species map."""
        if self.species not in SPECIES_CONFIGS:
            raise ValueError(f"Unknown species: {self.species}. "
                           f"Available: {list(SPECIES_CONFIGS.keys())}")
        return SPECIES_CONFIGS[self.species]  # Return the corresponding species configuration object

    def get_ai_backend_config(self, backend_name: str) -> Optional[AIBackendConfig]:
        """Return normalized settings for an AI backend."""
        if backend_name in self.ai_backends:
            return self.ai_backends[backend_name]
        return None

    def get_ai_api_key(self, backend_name: str) -> Optional[str]:
        """Resolve an AI API key from explicit settings or environment variables."""
        # Prefer the per-backend YAML setting.
        backend_config = self.get_ai_backend_config(backend_name)
        if backend_config and backend_config.api_key:
            return backend_config.api_key

        # Priority 2: Global
        if self.ai_api_key:
            return self.ai_api_key

        # Priority 3: Environmental variables
        env_var = f"{backend_name.upper()}_API_KEY"
        return os.getenv(env_var)

    def get_ai_model(self, backend_name: str, default: str) -> str:
        """Resolve the configured model name for an AI backend."""
        backend_config = self.get_ai_backend_config(backend_name)
        if backend_config and backend_config.model:
            return backend_config.model

        if self.ai_model:
            return self.ai_model

        return default

    def get_ai_base_url(self, backend_name: str) -> Optional[str]:
        """Resolve the optional base URL for an AI backend."""
        backend_config = self.get_ai_backend_config(backend_name)
        if backend_config:
            return backend_config.base_url
        return None

    def get_ai_group_id(self, backend_name: str) -> Optional[str]:
        """Resolve the optional group identifier for an AI backend."""
        backend_config = self.get_ai_backend_config(backend_name)
        if backend_config:
            return backend_config.group_id
        return None

    def validate(self) -> List[str]:
        """Validate method, thresholds, paths, plotting, and service settings."""
        errors = []  # Initializing error list

        # Validate input file
        # Validate optional input paths.
        if self.input_file and not Path(self.input_file).exists():
            errors.append(f"Input file not found: {self.input_file}")

        # Validate background file
        if self.background_file and not Path(self.background_file).exists():
            errors.append(f"Background file not found: {self.background_file}")

        # Validate species
        # Validate species code.
        species_valid = False
        
        # 1. Inspection of SPECIES_CONFIGS (separate backward compatibility)
        if self.species in SPECIES_CONFIGS:
            species_valid = True
        else:
            # 2. Query the species registry. Import lazily to avoid a circular dependency.
            try:
                from ..database.species_registry import SpeciesRegistry
                registry = SpeciesRegistry.load_default()
                
                # Try to press kegg_code query
                entry = registry.query_by_kegg_code(self.species)
                if entry:
                    species_valid = True
                else:
                    # Interpret a numeric species value as an NCBI TaxID.
                    if self.species.isdigit():
                        entry = registry.query_by_taxid(int(self.species))
                        if entry:
                            species_valid = True
            except Exception:
                pass  # Registry lookup is optional during configuration validation.
        
        if not species_valid:
            # Custom annotation databases may use species codes that are not in
            # the built-in registry. DatabaseManager performs the definitive
            # existence check, so validation only needs to reject unsafe paths.
            species_valid = (
                isinstance(self.species, str)
                and bool(self.species)
                and self.species not in {'.', '..'}
                and '..' not in self.species
                and all(char.isalnum() or char in '_.-' for char in self.species)
            )
        if not species_valid:
            errors.append(
                f"Invalid species code: {self.species}. Use a registered species or a "
                "custom code containing only letters, numbers, '_', '-' or '.'."
            )

        # Validate method
        # Validate the Eusenity Analysis
        valid_methods = [m.value for m in EnrichmentMethod]  # Get all valid method values
        if self.method not in valid_methods:
            errors.append(f"Invalid method: {self.method}. Valid: {valid_methods}")

        # Validate the multiple-testing correction method.
        valid_corrections = [c.value for c in CorrectionMethod]
        if self.correction not in valid_corrections:
            errors.append(f"Invalid correction: {self.correction}. Valid: {valid_corrections}")

        # Validate databases. Built-in names are enumerated; custom databases may
        # use any path-safe name because their existence is checked by DatabaseManager.
        valid_databases = [d.value for d in DatabaseType]  # Get all valid database values
        for db in self.databases:
            if db in valid_databases:
                continue
            if (
                not db or db in {'.', '..'} or '..' in db
                or not all(char.isalnum() or char in "_.-" for char in db)
            ):
                errors.append(
                    f"Invalid database name: {db}. Use a built-in database or a "
                    "custom name containing only letters, numbers, '_', '-' or '.'."
                )

        # Validate cutoffs
        # Verify whether the p and q thresholds are in the range (0, 1]
        if not 0 < self.pvalue_cutoff <= 1:
            errors.append(f"pvalue_cutoff must be in (0, 1], got: {self.pvalue_cutoff}")
        if not 0 < self.qvalue_cutoff <= 1:
            errors.append(f"qvalue_cutoff must be in (0, 1], got: {self.qvalue_cutoff}")

        valid_plot_formats = {"png", "pdf", "svg"}
        if not self.plot_formats:
            errors.append("plot_formats must contain at least one format")
        else:
            invalid_formats = [fmt for fmt in self.plot_formats if str(fmt).lower() not in valid_plot_formats]
            if invalid_formats:
                errors.append(f"Invalid plot_formats: {invalid_formats}. Valid: {sorted(valid_plot_formats)}")
        if self.plot_dpi <= 0:
            errors.append(f"plot_dpi must be positive, got: {self.plot_dpi}")
        if self.top_terms <= 0:
            errors.append(f"top_terms must be positive, got: {self.top_terms}")
        if self.activity_heatmap_top_n <= 0:
            errors.append(
                f"activity_heatmap_top_n must be positive, got: {self.activity_heatmap_top_n}"
            )
        if self.plot_width is not None and self.plot_width <= 0:
            errors.append(f"plot_width must be positive, got: {self.plot_width}")
        if self.plot_height is not None and self.plot_height <= 0:
            errors.append(f"plot_height must be positive, got: {self.plot_height}")

        valid_plot_styles = {"nature", "science", "presentation", "cell", "omicshare"}
        if not isinstance(self.plot_style, str) or self.plot_style not in valid_plot_styles:
            errors.append(
                f"Invalid plot_style: {self.plot_style}. Valid: {sorted(valid_plot_styles)}"
            )
        if self.plot_palette is not None:
            from allenricher.visualization.color_config import PALETTES, PUBLIC_PALETTES
            if not isinstance(self.plot_palette, str) or self.plot_palette not in PALETTES:
                errors.append(
                    f"Invalid plot_palette: {self.plot_palette}. Valid: {sorted(PUBLIC_PALETTES)}"
                )
        palette_fields = {
            "categorical_palette": (
                self.categorical_palette,
                "CATEGORICAL_PALETTES",
            ),
            "sequential_palette": (
                self.sequential_palette,
                "SEQUENTIAL_PALETTES",
            ),
            "diverging_palette": (
                self.diverging_palette,
                "DIVERGING_PALETTES",
            ),
        }
        if any(value is not None for value, _ in palette_fields.values()):
            from allenricher.visualization import color_config
            for field_name, (value, registry_name) in palette_fields.items():
                registry = getattr(color_config, registry_name)
                if value is not None and (
                    not isinstance(value, str) or value not in registry
                ):
                    errors.append(
                        f"Invalid {field_name}: {value}. Valid: {sorted(registry)}"
                    )

        return errors  # Returns all errors


# Template written by `allenricher config`.
DEFAULT_CONFIG_YAML = """# AllEnricher v2 configuration

# Input and output
input_file: null
output_dir: "./results"
background_file: null

# Analysis scope
species: "hsa"
databases:
  - "GO"
  - "KEGG"
method: "hypergeometric"  # hypergeometric, gsea, ssgsea, or gsva
correction: "BH"          # BH, BY, bonferroni, holm, or none
pvalue_cutoff: 0.05
qvalue_cutoff: 0.05
min_genes: 3
max_genes: .inf
output_all: true

# GSEA, ssGSEA, and GSVA
gsea_permutations: 1000
gsea_min_size: null        # method default: GSEA 15; ssGSEA/GSVA 1
gsea_max_size: null        # method default: GSEA 500; ssGSEA/GSVA unbounded

# Figures
plot_formats:
  - "pdf"
  - "png"
plot_dpi: 300
top_terms: 20
activity_heatmap_top_n: 40
plot_width: null           # adaptive when null
plot_height: null          # adaptive when null
plot_style: "nature"
plot_palette: null         # legacy semantic-role override
categorical_palette: null
sequential_palette: null
diverging_palette: null

# Performance and reports
n_jobs: 1
generate_report: true
report_format: "html"
ai_interpretation: false
ai_backend: "openai"

# Optional AI backends
ai_backends:
  openai:
    api_key: null
    model: "gpt-4"
    enabled: true
  claude:
    api_key: null
    model: "claude-3-opus-20240229"
    enabled: true
  deepseek:
    api_key: null
    model: "deepseek-chat"
    enabled: true
  glm:
    api_key: null
    model: "glm-4"
    enabled: true
  minimax:
    api_key: null
    group_id: null
    model: "abab6.5s-chat"
    enabled: true
  ollama:
    model: "llama2"
    base_url: "http://localhost:11434"
    enabled: true

# Database storage
database_dir: "./database"
auto_update: false
"""
