"""Define semantic categorical, sequential, and diverging color palettes."""

import logging
from dataclasses import dataclass
from typing import Dict, List, Optional, Union

from matplotlib.colors import to_rgb


logger = logging.getLogger(__name__)


MIN_GRADIENT_DISTANCE_FROM_WHITE = 0.18
VISIBLE_GRADIENT_NEUTRAL = "#D9D9D9"
DIVERGING_WHITE_MIDPOINT = "#FFFFFF"


# =============================================================================
# Paul Tol palettes: color-vision-deficiency friendly and high contrast.
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

COLORBREWER_PURD = [
    "#F7F4F9", "#E7E1EF", "#D4B9DA", "#C994C7", "#DF65B0",
    "#E7298A", "#CE1256", "#980043", "#67001F"
]

COLORBREWER_BLUES = [
    "#F7FBFF", "#DEEBF7", "#C6DBEF", "#9ECAE1", "#6BAED6",
    "#4292C6", "#2171B5", "#08519C", "#08306B"
]

VIRIDIS = [
    "#440154", "#46327E", "#365C8D", "#277F8E", "#1FA187",
    "#4AC16D", "#A0DA39", "#FDE725"
]

CIVIDIS = [
    "#00204C", "#2E3F6D", "#575D6D", "#7C7B78", "#A59C74",
    "#D2C060", "#FFEA46"
]

COLORBREWER_RDBU = [
    "#053061", "#2166AC", "#4393C3", "#92C5DE", "#D1E5F0",
    "#F7F7F7", "#FDDBC7", "#F4A582", "#D6604D", "#B2182B", "#67001F"
]

COLORBREWER_PRGN = [
    "#762A83", "#9970AB", "#C2A5CF", "#E7D4E8", "#F7F7F7",
    "#D9F0D3", "#ACD39E", "#5AAE61", "#1B7837"
]

COLORBREWER_BRBG = [
    "#543005", "#8C510A", "#BF812D", "#DFC27D", "#F6E8C3",
    "#F5F5F5", "#C7EAE5", "#80CDC1", "#35978F", "#01665E", "#003C30"
]


# =============================================================================
# Okabe-Ito palette: established color-vision-deficiency-friendly categories.
# =============================================================================

OKABE_ITO = [
    "#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442",
    "#0072B2", "#D55E00", "#CC79A7"
]


# =============================================================================
# Categorical palettes inspired by common journal figure conventions.
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
# Palettes retained for compatibility with common bioinformatics figures.
# =============================================================================

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
# ECharts 4 Default Palette
# =============================================================================

ECHARTS_V4_COLORS = [
    "#C23531", "#2F4554", "#61A0A8", "#D48265", "#91C7AE",
    "#749F83", "#CA8622", "#BDA29A", "#6E7074", "#546570"
]


# Internal fallback used only when a categorical palette contains too few colors.
HIGH_CARDINALITY_CATEGORICAL = [
    "#1F77B4", "#FF7F0E", "#2CA02C", "#D62728", "#9467BD",
    "#8C564B", "#E377C2", "#7F7F7F", "#BCBD22", "#17BECF",
    "#393B79", "#637939", "#8C6D31", "#843C39", "#7B4173",
    "#3182BD", "#31A354", "#756BB1", "#E6550D", "#969696",
]


# =============================================================================
# Palette registry organized by semantic data role.
# =============================================================================

CATEGORICAL_PALETTES: Dict[str, List[str]] = {
    "tol_bright": TOL_BRIGHT,
    "tol_high_contrast": TOL_HIGH_CONTRAST,
    "tol_vibrant": TOL_VIBRANT,
    "tol_muted": TOL_MUTED,
    "tol_medium_contrast": TOL_MEDIUM_CONTRAST,
    "tol_light": TOL_LIGHT,
    "okabe_ito": OKABE_ITO,
    "nature": NATURE_COLORS,
    "science": SCIENCE_COLORS,
    "cell": CELL_COLORS,
    "lancet": LANCET_COLORS,
    "nejm": NEJM_COLORS,
    "jama": JAMA_COLORS,
    "omicshare": OMICSHARE_COLORS,
    "echarts_v4": ECHARTS_V4_COLORS,
}

SEQUENTIAL_PALETTES: Dict[str, List[str]] = {
    "colorbrewer_blues": COLORBREWER_BLUES,
    "colorbrewer_purd": COLORBREWER_PURD,
    "viridis": VIRIDIS,
    "cividis": CIVIDIS,
}

# Sequential scales use two related anchors selected from each registered palette.
SEQUENTIAL_GRADIENT_ANCHORS: Dict[str, List[str]] = {
    "colorbrewer_blues": ["#9ECAE1", "#08519C"],
    "colorbrewer_purd": ["#DF65B0", "#980043"],
    "viridis": ["#365C8D", "#1FA187"],
    "cividis": ["#00204C", "#575D6D"],
}

DIVERGING_PALETTES: Dict[str, List[str]] = {
    "colorbrewer_rdbu": COLORBREWER_RDBU,
    "tol_sunset": TOL_SUNSET,
    "colorbrewer_prgn": COLORBREWER_PRGN,
    "colorbrewer_brbg": COLORBREWER_BRBG,
}

DEFAULT_CATEGORICAL_PALETTE = "tol_bright"
DEFAULT_SEQUENTIAL_PALETTE = "colorbrewer_blues"
DEFAULT_DIVERGING_PALETTE = "colorbrewer_rdbu"

# PALETTES is a compatibility index; new code should request a semantic role.
PALETTES: Dict[str, List[str]] = {
    "default": CATEGORICAL_PALETTES[DEFAULT_CATEGORICAL_PALETTE],
    **CATEGORICAL_PALETTES,
    **SEQUENTIAL_PALETTES,
    **DIVERGING_PALETTES,
}
PUBLIC_CATEGORICAL_PALETTES = tuple(CATEGORICAL_PALETTES)
PUBLIC_SEQUENTIAL_PALETTES = tuple(SEQUENTIAL_PALETTES)
PUBLIC_DIVERGING_PALETTES = tuple(DIVERGING_PALETTES)
PUBLIC_PALETTES = (
    PUBLIC_CATEGORICAL_PALETTES
    + PUBLIC_SEQUENTIAL_PALETTES
    + PUBLIC_DIVERGING_PALETTES
)
PALETTE_ROLES = {
    **{name: "categorical" for name in CATEGORICAL_PALETTES},
    **{name: "sequential" for name in SEQUENTIAL_PALETTES},
    **{name: "diverging" for name in DIVERGING_PALETTES},
}


@dataclass(frozen=True)
class PaletteSelection:
    """Store independent categorical, sequential, and diverging palette selections."""

    categorical: str = DEFAULT_CATEGORICAL_PALETTE
    sequential: str = DEFAULT_SEQUENTIAL_PALETTE
    diverging: str = DEFAULT_DIVERGING_PALETTE

    def for_role(self, role: str) -> str:
        if role not in {"categorical", "sequential", "diverging"}:
            raise ValueError(f"Unknown palette role: {role}")
        return getattr(self, role)


PaletteLike = Optional[Union[str, PaletteSelection]]


def get_palette_role(name: str) -> str:
    """Return the semantic role assigned to a palette name."""
    if name == "default":
        return "categorical"
    if name not in PALETTE_ROLES:
        raise ValueError(f"Unknown palette '{name}'. Available: {', '.join(PUBLIC_PALETTES)}")
    return PALETTE_ROLES[name]


def resolve_palette_selection(
    legacy_palette: Optional[str] = None,
    categorical_palette: Optional[str] = None,
    sequential_palette: Optional[str] = None,
    diverging_palette: Optional[str] = None,
) -> PaletteSelection:
    """Resolve legacy and role-specific palette arguments into one selection."""
    selected = {
        "categorical": DEFAULT_CATEGORICAL_PALETTE,
        "sequential": DEFAULT_SEQUENTIAL_PALETTE,
        "diverging": DEFAULT_DIVERGING_PALETTE,
    }

    if legacy_palette and legacy_palette != "default":
        selected[get_palette_role(legacy_palette)] = legacy_palette

    explicit = {
        "categorical": categorical_palette,
        "sequential": sequential_palette,
        "diverging": diverging_palette,
    }
    registries = {
        "categorical": CATEGORICAL_PALETTES,
        "sequential": SEQUENTIAL_PALETTES,
        "diverging": DIVERGING_PALETTES,
    }
    for role, name in explicit.items():
        if name is None:
            continue
        if name not in registries[role]:
            available = ", ".join(registries[role])
            raise ValueError(f"Palette '{name}' is not {role}. Available: {available}")
        selected[role] = name

    return PaletteSelection(**selected)


def coerce_palette_selection(palette: PaletteLike) -> PaletteSelection:
    if isinstance(palette, PaletteSelection):
        return palette
    return resolve_palette_selection(legacy_palette=palette)


def palette_name_for_role(palette: PaletteLike, role: str) -> str:
    return coerce_palette_selection(palette).for_role(role)


def visible_gradient_colors(
    colors: List[str], role: str, palette_name: Optional[str] = None
) -> List[str]:
    """Return perceptually distinguishable anchors for a sequential or diverging scale."""
    if role not in {"sequential", "diverging"}:
        raise ValueError(f"Gradient visibility requires a continuous role, got: {role}")

    protected = list(colors)

    def is_near_white(color: str) -> bool:
        rgb = to_rgb(color)
        distance = sum((1.0 - channel) ** 2 for channel in rgb) ** 0.5
        return distance < MIN_GRADIENT_DISTANCE_FROM_WHITE

    if role == "sequential":
        if palette_name in SEQUENTIAL_GRADIENT_ANCHORS:
            return list(SEQUENTIAL_GRADIENT_ANCHORS[palette_name])
        visible = [color for color in protected if not is_near_white(color)]
        if len(visible) >= 2:
            return visible

    if role == "diverging":
        return [DIVERGING_WHITE_MIDPOINT if is_near_white(color) else color for color in protected]

    return [VISIBLE_GRADIENT_NEUTRAL if is_near_white(color) else color for color in protected]


def categorical_colors(palette: PaletteLike, n: int) -> List[str]:
    """Return distinct categorical colors, extending with a fallback palette when necessary."""
    if n < 0:
        raise ValueError("n must be non-negative")
    name = palette_name_for_role(palette, "categorical")
    colors = CATEGORICAL_PALETTES[name]
    if n <= len(colors):
        return colors[:n]
    if n > len(HIGH_CARDINALITY_CATEGORICAL):
        raise ValueError(
            f"Categorical plot requires {n} distinct colors; maximum supported is "
            f"{len(HIGH_CARDINALITY_CATEGORICAL)}"
        )
    logger.warning(
        "Categorical palette '%s' has %d colors; using the %d-color fallback instead",
        name,
        len(colors),
        n,
    )
    return HIGH_CARDINALITY_CATEGORICAL[:n]


class ColorConfig:
    """Resolve semantic palettes into colors and Matplotlib colormaps."""
    
    def __init__(self):
        """Initialize semantic palettes from defaults or a configuration mapping."""
        self._palettes = PALETTES.copy()
    
    def get_available_palettes(self) -> List[str]:
        """Return all registered semantic palette names."""
        return list(PUBLIC_PALETTES)
    
    def get_colors(self, palette_name: str = 'default', n: int = 8) -> List[str]:
        """Return a requested number of colors from the active palette."""
        if palette_name not in self._palettes:
            raise ValueError(f"Unknown palette: {palette_name}. Available palettes: {list(PUBLIC_PALETTES)}")
        
        if palette_name == "default" or get_palette_role(palette_name) == "categorical":
            return categorical_colors(palette_name, n)

        palette = self._palettes[palette_name]
        if n > len(palette):
            raise ValueError(
                f"Continuous palette '{palette_name}' provides color anchors, not "
                "repeatable categories; use a sequential/diverging colormap"
            )
        return palette[:n]
    
    def get_categorical_colors(self, category_type: str, palette: PaletteLike = None) -> Dict[str, str]:
        """Return distinct colors for unordered categories."""
        if category_type.lower() == 'go':
            return self.get_go_category_colors(palette)
        elif category_type.lower() == 'kegg':
            return self.get_kegg_category_colors(palette)
        else:
            raise ValueError(f"Unsupported category type: {category_type}. Expected 'go' or 'kegg'.")
    
    def get_go_category_colors(self, palette: PaletteLike = None) -> Dict[str, str]:
        """Return stable colors for displayed GO namespaces."""
        colors = categorical_colors(palette, n=3)
        return {
            "biological_process": colors[0],
            "cellular_component": colors[1],
            "molecular_function": colors[2],
        }
    
    def get_kegg_category_colors(self, palette: PaletteLike = None) -> Dict[str, str]:
        """Return stable colors for displayed KEGG categories."""
        colors = categorical_colors(palette, n=6)
        return {
            "Genetic_Information_Processing": colors[0],
            "Human_Diseases": colors[1],
            "Metabolism": colors[2],
            "Cellular_Processes": colors[3],
            "Organismal_Systems": colors[4],
            "Environmental_Information_Processing": colors[5],
        }
    
    def get_palette_colors(self, palette_name: str) -> List[str]:
        """Return the discrete color sequence registered for a palette."""
        if palette_name not in self._palettes:
            raise ValueError(f"Unknown palette: {palette_name}")
        return self._palettes[palette_name].copy()


# Process-wide color configuration.
_color_config: Optional[ColorConfig] = None


def get_color_config() -> ColorConfig:
    """Return the process-wide color configuration."""
    global _color_config
    if _color_config is None:
        _color_config = ColorConfig()
    return _color_config
