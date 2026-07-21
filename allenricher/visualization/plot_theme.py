"""Define and resolve reusable publication figure styles."""

from __future__ import annotations

import contextlib
import os
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple, Union

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.collections import Collection
from matplotlib.colors import LinearSegmentedColormap, ListedColormap, to_hex
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from matplotlib.text import Text

from .color_config import (
    CATEGORICAL_PALETTES,
    DEFAULT_CATEGORICAL_PALETTE,
    DEFAULT_DIVERGING_PALETTE,
    DEFAULT_SEQUENTIAL_PALETTE,
    DIVERGING_PALETTES,
    PALETTES,
    PUBLIC_PALETTES,
    SEQUENTIAL_PALETTES,
    PaletteLike,
    categorical_colors,
    palette_name_for_role,
    visible_gradient_colors,
)


# =============================================================================
# Publication figure-style configuration.
# =============================================================================

@dataclass
class StylePreset:
    """Store typography, spacing, and line settings for a publication figure style."""

    # Identity.
    name: str                           # Style name
    # Typography.
    font_family: str = "sans-serif"     # Font family
    font_weight: str = "normal"         # Base font weight
    font_size: float = 10               # Base font size
    title_size: float = 12              # Figure-title size
    title_weight: str = "bold"          # Figure-title weight
    label_size: float = 10              # Axis-label size
    tick_label_size: float = 9          # Tick-label size
    legend_size: float = 9              # Legend-label size
    legend_title_size: float = 10       # Legend-title size

    # Lines.
    line_width: float = 1.0             # Width of data lines
    line_style: str = "solid"           # Line Style (solid, dashed, dotted)
    marker_size: float = 6.0            # Point Size
    marker_edge_width: float = 0.5      # Margin width of the point mark

    # Axes borders. =
    axes_line_width: float = 0.8        # Axis Width
    spine_top: bool = True              # Whether to display the upper border
    spine_right: bool = True            # Whether to show the right box
    spine_color: str = "#000000"        # Border Colour

    # Ticks.
    tick_direction: str = "out"         # Queue (in, out, inout)
    tick_major_size: float = 5.0        # Length of main scale
    tick_major_width: float = 0.8       # Width of main scale
    tick_minor_size: float = 3.0        # Sub-ticity Length
    tick_minor_visible: bool = False    # Whether or not to display the lower scale

    # Grid.
    grid: bool = False                  # Whether to show the grid
    grid_axis: str = "both"             # Grid axes (both, x, y)
    grid_alpha: float = 0.3             # Grid Transparency
    grid_linewidth: float = 0.5         # Wide Grid Line
    grid_color: str = "#CCCCCC"         # Grid Colour
    grid_linestyle: str = "dashed"      # Grid Style (solid, dashed, dotted)

    # Legend and figure output.
    figure_dpi: int = 300               # Graphic DPI
    figure_format: str = "png"          # Default Output Format
    savefig_dpi: int = 300              # Save DPI
    savefig_format: str = "png"         # Save Format
    savefig_bbox: str = "tight"         # Save Border Box (tight, standard)
    savefig_pad: float = 0.1            # Save Margin

    # Background.
    facecolor: str = "white"            # Graphic Background Color
    axes_facecolor: str = "white"       # Background color for the plot area

    # Layout.
    context: str = "paper"              # Seeaborn context (paper, notebook, talk, poster)

    # Optional Matplotlib overrides.
    rc_overrides: Dict[str, Any] = field(default_factory=dict)


# =============================================================================
# Registered style presets.
# =============================================================================

PRESETS: Dict[str, StylePreset] = {
    "nature": StylePreset(
        name="nature",
        # Typography.
        font_family="sans-serif",
        font_weight="normal",
        font_size=8,
        title_size=10,
        title_weight="bold",
        label_size=8,
        tick_label_size=7,
        legend_size=7,
        legend_title_size=8,
        # Lines.
        line_width=1.0,
        line_style="solid",
        marker_size=4.0,
        marker_edge_width=0.5,
        # Axes borders.
        axes_line_width=0.5,
        spine_top=False,
        spine_right=False,
        spine_color="#333333",
        # Ticks.
        tick_direction="in",
        tick_major_size=4.0,
        tick_major_width=0.5,
        tick_minor_size=2.0,
        tick_minor_visible=False,
        # Grid.
        grid=False,
        grid_axis="both",
        grid_alpha=0.3,
        grid_linewidth=0.5,
        grid_color="#CCCCCC",
        grid_linestyle="dashed",
        # Legend and figure output.
        figure_dpi=300,
        savefig_dpi=300,
        savefig_bbox="tight",
        savefig_pad=0.1,
        # Background.
        facecolor="white",
        axes_facecolor="white",
        # Layout.
        context="paper",
    ),

    "science": StylePreset(
        name="science",
        # Typography.
        font_family="serif",
        font_weight="normal",
        font_size=8.4,
        title_size=10.5,
        title_weight="bold",
        label_size=8.4,
        tick_label_size=7.4,
        legend_size=7.4,
        legend_title_size=8.4,
        # Lines.
        line_width=1.15,
        line_style="solid",
        marker_size=4.6,
        marker_edge_width=0.6,
        # Axes borders.
        axes_line_width=0.8,
        spine_top=True,
        spine_right=True,
        spine_color="#000000",
        # Ticks.
        tick_direction="out",
        tick_major_size=5.0,
        tick_major_width=0.7,
        tick_minor_size=3.0,
        tick_minor_visible=True,
        # Grid.
        grid=False,
        grid_axis="both",
        grid_alpha=0.3,
        grid_linewidth=0.5,
        grid_color="#CCCCCC",
        grid_linestyle="dashed",
        # Legend and figure output.
        figure_dpi=300,
        savefig_dpi=300,
        savefig_bbox="tight",
        savefig_pad=0.1,
        # Background.
        facecolor="white",
        axes_facecolor="white",
        # Layout.
        context="paper",
    ),

    "presentation": StylePreset(
        name="presentation",
        # Typography.
        font_family="sans-serif",
        font_weight="normal",
        font_size=10.4,
        title_size=13.0,
        title_weight="bold",
        label_size=10.4,
        tick_label_size=9.1,
        legend_size=9.1,
        legend_title_size=10.4,
        # Lines.
        line_width=1.4,
        line_style="solid",
        marker_size=7.0,
        marker_edge_width=0.9,
        # Axes borders.
        axes_line_width=1.0,
        spine_top=True,
        spine_right=True,
        spine_color="#000000",
        # Ticks.
        tick_direction="out",
        tick_major_size=6.0,
        tick_major_width=0.9,
        tick_minor_size=4.0,
        tick_minor_visible=True,
        # Grid.
        grid=True,
        grid_axis="both",
        grid_alpha=0.42,
        grid_linewidth=0.7,
        grid_color="#EEEEEE",
        grid_linestyle="solid",
        # Legend and figure output.
        figure_dpi=300,
        savefig_dpi=300,
        savefig_bbox="tight",
        savefig_pad=0.2,
        # Background.
        facecolor="white",
        axes_facecolor="white",
        # Layout.
        context="talk",
    ),

}

STYLE_ALIASES = {"cell": "nature", "omicshare": "science"}


def resolve_style(style: Optional[str]) -> str:
    """Resolve a style name while honoring documented compatibility aliases."""
    name = STYLE_ALIASES.get((style or "nature").lower(), (style or "nature").lower())
    if name not in PRESETS:
        available = ", ".join(PRESETS)
        raise ValueError(f"Unknown style '{style}'. Available: {available}")
    return name


def _resolved_font_family(preset: StylePreset) -> str:
    fallback = "DejaVu Serif" if preset.font_family.lower() == "serif" else "DejaVu Sans"
    for family in (preset.font_family, fallback):
        try:
            font_manager.findfont(family, fallback_to_default=False)
            return family
        except ValueError:
            continue
    return fallback


def apply_figure_style(
    fig,
    style: Optional[str],
    axes=None,
    grid_axis: Optional[str] = None,
    border: Optional[str] = "style",
    max_text_scale: float = 1.35,
) -> str:
    """Apply shared typography and axis styling to an existing Matplotlib figure."""
    canonical = resolve_style(style)
    preset = PRESETS[canonical]
    base = PRESETS["nature"]
    text_scale = min(max_text_scale, preset.font_size / base.font_size)
    line_scale = preset.line_width / base.line_width
    family = _resolved_font_family(preset)
    preset_text_sizes = {
        preset.font_size,
        preset.title_size,
        preset.label_size,
        preset.tick_label_size,
        preset.legend_size,
        preset.legend_title_size,
    }

    for text in fig.findobj(match=Text):
        text.set_fontfamily(family)
        if canonical != "nature" and not any(
            abs(text.get_fontsize() - size) < 1e-9 for size in preset_text_sizes
        ):
            text.set_fontsize(text.get_fontsize() * text_scale)
    for line in fig.findobj(match=Line2D):
        if line.get_linewidth() > 0 and abs(line.get_linewidth() - preset.line_width) > 1e-9:
            line.set_linewidth(line.get_linewidth() * line_scale)
    for patch in fig.findobj(match=Patch):
        if patch.get_linewidth() > 0:
            patch.set_linewidth(patch.get_linewidth() * line_scale)
    for collection in fig.findobj(match=Collection):
        widths = collection.get_linewidths()
        if len(widths):
            collection.set_linewidths(widths * line_scale)

    target_axes = list(fig.axes if axes is None else axes)
    for axis in target_axes:
        if not axis.axison:
            continue
        if border in {"style", "full"}:
            visibility = {
                "left": True,
                "bottom": True,
                "top": border == "full" or preset.spine_top,
                "right": border == "full" or preset.spine_right,
            }
            for side, visible in visibility.items():
                axis.spines[side].set_visible(visible)
                axis.spines[side].set_color(preset.spine_color)
                axis.spines[side].set_linewidth(preset.axes_line_width)
        if grid_axis is not None:
            if preset.grid:
                axis.grid(
                    True,
                    axis=grid_axis,
                    color=preset.grid_color,
                    alpha=preset.grid_alpha,
                    linewidth=preset.grid_linewidth,
                    linestyle=preset.grid_linestyle,
                )
            else:
                axis.grid(False, axis=grid_axis)
        axis.tick_params(
            direction=preset.tick_direction,
            length=preset.tick_major_size,
            width=preset.tick_major_width,
        )
    return canonical


# =============================================================================
# Figure output helpers.
# =============================================================================

def save_figure(
    fig,
    output_path: str,
    dpi: int = 300,
    bbox_inches: str = "tight",
    facecolor: Optional[str] = None,
    **kwargs
) -> str:
    """Save one figure as PNG, PDF, or SVG according to its filename."""
    path = Path(output_path)
    extension = path.suffix.lower().lstrip(".")
    if extension not in {"png", "pdf", "svg"}:
        raise ValueError("Figure output format must be png, pdf, or svg")

    path.parent.mkdir(parents=True, exist_ok=True)
    if facecolor is None:
        facecolor = fig.get_facecolor()

    save_kwargs = {
        "format": extension,
        "bbox_inches": bbox_inches,
        "facecolor": facecolor,
        **kwargs,
    }
    if extension == "png":
        save_kwargs["dpi"] = dpi
    fig.savefig(str(path), **save_kwargs)
    return str(path)

def save_figure_dual(
    fig,
    output_path: str,
    dpi: int = 300,
    bbox_inches: str = "tight",
    facecolor: Optional[str] = None,
    **kwargs
) -> Tuple[str, str]:
    """Save PNG and PDF copies, plus a separately requested third format."""
    requested_path = Path(output_path)
    requested_extension = requested_path.suffix.lower()
    # Remove the extension before creating both output filenames.
    base_path = output_path
    for ext in ['.png', '.pdf', '.jpg', '.jpeg', '.svg', '.eps']:
        if base_path.lower().endswith(ext):
            base_path = base_path[:-len(ext)]
            break

    png_path = base_path + '.png'
    pdf_path = base_path + '.pdf'

    # Preserve the requested background in both formats.
    if facecolor is None:
        facecolor = fig.get_facecolor()

    # Raster output.
    fig.savefig(
        png_path,
        format='png',
        dpi=dpi,
        bbox_inches=bbox_inches,
        facecolor=facecolor,
        **kwargs
    )

    # Vector output.
    fig.savefig(
        pdf_path,
        format='pdf',
        bbox_inches=bbox_inches,
        facecolor=facecolor,
        **kwargs
    )

    if requested_extension and requested_extension not in {".png", ".pdf"}:
        requested_path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(
            str(requested_path),
            format=requested_extension.lstrip("."),
            bbox_inches=bbox_inches,
            facecolor=facecolor,
            **kwargs,
        )

    return png_path, pdf_path


# =============================================================================
# PlotTheme facade.
# =============================================================================

class PlotTheme:
    """Resolve typography, spacing, and color defaults for one figure style."""

    _active_preset: Optional[str] = None
    _active_config: Optional[StylePreset] = None
    _original_rc: Optional[Dict[str, Any]] = None

    @classmethod
    def available_styles(cls) -> List[str]:
        """Return all registered figure-style names."""
        return list(PRESETS.keys())

    @classmethod
    def available_palettes(cls) -> List[str]:
        """Return all registered palette names."""
        return list(PUBLIC_PALETTES)

    @classmethod
    def apply(cls, style: str) -> None:
        """Apply the selected style to Matplotlib defaults."""
        style = resolve_style(style)
        preset = PRESETS[style]
        cls._active_preset = style
        cls._active_config = preset

        # Capture Matplotlib defaults once so reset() is deterministic.
        if cls._original_rc is None:
            cls._original_rc = plt.rcParams.copy()

        resolved_font = _resolved_font_family(preset)

        # Map the style preset to Matplotlib runtime parameters.
        rc_params = {
            # Typography.
            "font.family": resolved_font,
            "font.sans-serif": ["Helvetica", "Arial", "DejaVu Sans", "Liberation Sans",
                                 "Tahoma", "Verdana", "sans-serif"],
            "font.weight": preset.font_weight,
            "font.size": preset.font_size,
            "axes.titlesize": preset.title_size,
            "axes.titleweight": preset.title_weight,
            "axes.labelsize": preset.label_size,
            "xtick.labelsize": preset.tick_label_size,
            "ytick.labelsize": preset.tick_label_size,
            "legend.fontsize": preset.legend_size,
            "legend.title_fontsize": preset.legend_title_size,

            # Lines.
            "lines.linewidth": preset.line_width,
            "lines.linestyle": preset.line_style,
            "lines.markersize": preset.marker_size,
            "lines.markeredgewidth": preset.marker_edge_width,

            # Axes borders.
            "axes.linewidth": preset.axes_line_width,
            "axes.spines.top": preset.spine_top,
            "axes.spines.right": preset.spine_right,
            "axes.edgecolor": preset.spine_color,

            # Ticks.
            "xtick.direction": preset.tick_direction,
            "ytick.direction": preset.tick_direction,
            "xtick.major.size": preset.tick_major_size,
            "ytick.major.size": preset.tick_major_size,
            "xtick.major.width": preset.tick_major_width,
            "ytick.major.width": preset.tick_major_width,
            "xtick.minor.size": preset.tick_minor_size,
            "ytick.minor.size": preset.tick_minor_size,
            "xtick.minor.visible": preset.tick_minor_visible,
            "ytick.minor.visible": preset.tick_minor_visible,

            # Layout.
            "axes.grid": preset.grid,
            "axes.grid.axis": preset.grid_axis,
            "grid.alpha": preset.grid_alpha,
            "grid.linewidth": preset.grid_linewidth,
            "grid.color": preset.grid_color,
            "grid.linestyle": preset.grid_linestyle,

            # Legend and figure output.
            "figure.dpi": preset.figure_dpi,
            "savefig.dpi": preset.savefig_dpi,
            "savefig.bbox": preset.savefig_bbox,
            "savefig.pad_inches": preset.savefig_pad,

            # Background.
            "figure.facecolor": preset.facecolor,
            "axes.facecolor": preset.axes_facecolor,

            # Style controls layout; the default color cycle comes from the palette.
            "axes.prop_cycle": plt.cycler(
                color=CATEGORICAL_PALETTES[DEFAULT_CATEGORICAL_PALETTE]
            ),
        }

        # Explicit preset overrides take final precedence.
        rc_params.update(preset.rc_overrides)

        # Apply the resolved parameters process-wide.
        plt.rcParams.update(rc_params)

    @classmethod
    @contextlib.contextmanager
    def context(cls, style: str, palette: PaletteLike = None):
        """Temporarily apply a figure style inside a context manager."""
        # Record the active style for helper functions.
        previous_preset = cls._active_preset
        previous_rc = plt.rcParams.copy()

        try:
            cls.apply(style)
            # Matplotlib's property cycle is categorical by definition.
            plt.rcParams["axes.prop_cycle"] = plt.cycler(
                color=categorical_colors(palette, n=len(
                    CATEGORICAL_PALETTES[
                        palette_name_for_role(palette, "categorical")
                    ]
                ))
            )
            yield cls
        finally:
            # Restore the style that was active before entering the context.
            cls._active_preset = previous_preset
            if previous_preset:
                cls._active_config = PRESETS.get(previous_preset)
            else:
                cls._active_config = None
            plt.rcParams.update(previous_rc)

    @classmethod
    def get_active(cls) -> Optional[StylePreset]:
        """Return the active figure style."""
        return cls._active_config

    @classmethod
    def get_palette(cls, name: Optional[str] = None, n: Optional[int] = None, palette: Optional[str] = None) -> List[str]:
        """Return the named categorical palette."""
        # An explicit palette argument takes precedence over the style default.
        if palette is not None:
            name = palette

        if name is None:
            name = DEFAULT_CATEGORICAL_PALETTE

        if name not in PALETTES:
            available = ", ".join(cls.available_palettes())
            raise ValueError(f"Unknown palette '{name}'. Available: {available}")

        colors = PALETTES[name]
        if name in SEQUENTIAL_PALETTES:
            colors = visible_gradient_colors(colors, "sequential", name)
        elif name in DIVERGING_PALETTES:
            colors = visible_gradient_colors(colors, "diverging")

        if n is None:
            return colors

        if name == "default" or name in CATEGORICAL_PALETTES:
            return categorical_colors(name, n)
        if n <= len(colors):
            return colors[:n]
        cmap = LinearSegmentedColormap.from_list(name, colors)
        if n == 1:
            return [to_hex(cmap(0.5))]
        return [to_hex(cmap(index / (n - 1))) for index in range(n)]

    @classmethod
    def get_plot_colors(
        cls,
        style: Optional[str] = "nature",
        palette: PaletteLike = None,
        default: Optional[List[str]] = None,
        n: Optional[int] = None,
        divergent: bool = False,
        sequential: bool = False,
        role: Optional[str] = None,
    ) -> List[str]:
        """Resolve colors from semantic data roles independently of figure style."""
        style = resolve_style(style)

        if role is None:
            role = "diverging" if divergent else "sequential" if sequential else "categorical"
        if role not in {"categorical", "sequential", "diverging"}:
            raise ValueError(f"Unknown palette role: {role}")

        # The requested count is retained for API compatibility; semantic role selects the scale.
        _ = default
        name = palette_name_for_role(palette, role)
        if role == "categorical":
            count = len(CATEGORICAL_PALETTES[name]) if n is None else n
            return categorical_colors(palette, count)

        registry = SEQUENTIAL_PALETTES if role == "sequential" else DIVERGING_PALETTES
        colors = visible_gradient_colors(list(registry[name]), role, name)
        if n is None:
            return colors
        cmap = LinearSegmentedColormap.from_list(name, colors)
        if n == 1:
            return [to_hex(cmap(0.5))]
        return [to_hex(cmap(index / (n - 1))) for index in range(n)]

    @classmethod
    def get_sequential_cmap(
        cls,
        name: Optional[str] = None,
        colors: Optional[List[str]] = None
    ) -> LinearSegmentedColormap:
        """Return a sequential colormap for ordered magnitudes."""
        registered_name = None
        if colors is None:
            name = name or DEFAULT_SEQUENTIAL_PALETTE
            if name not in SEQUENTIAL_PALETTES:
                available = ", ".join(SEQUENTIAL_PALETTES)
                raise ValueError(f"Unknown sequential palette '{name}'. Available: {available}")
            colors = SEQUENTIAL_PALETTES[name]
            registered_name = name

        colors = visible_gradient_colors(list(colors), "sequential", registered_name)

        return LinearSegmentedColormap.from_list("sequential", colors)

    @staticmethod
    def get_continuous_cmap(
        palette: str = DEFAULT_SEQUENTIAL_PALETTE,
    ) -> LinearSegmentedColormap:
        """Return a continuous colormap for the requested semantic role."""
        if palette not in SEQUENTIAL_PALETTES:
            available = ", ".join(SEQUENTIAL_PALETTES)
            raise ValueError(f"Unknown sequential palette '{palette}'. Available: {available}")
        colors = visible_gradient_colors(
            list(SEQUENTIAL_PALETTES[palette]), "sequential", palette
        )
        return LinearSegmentedColormap.from_list(f"{palette}_continuous", colors)

    @classmethod
    def get_diverging_cmap(
        cls,
        name: Optional[str] = None,
        colors: Optional[List[str]] = None
    ) -> LinearSegmentedColormap:
        """Return a diverging colormap centered on a meaningful midpoint."""
        if colors is None:
            name = name or DEFAULT_DIVERGING_PALETTE
            if name not in DIVERGING_PALETTES:
                available = ", ".join(DIVERGING_PALETTES)
                raise ValueError(f"Unknown diverging palette '{name}'. Available: {available}")
            colors = DIVERGING_PALETTES[name]

        colors = visible_gradient_colors(list(colors), "diverging")

        return LinearSegmentedColormap.from_list("diverging", colors)

    @classmethod
    def get_category_colors(
        cls,
        categories: List[str],
        palette: PaletteLike = None
    ) -> Dict[str, str]:
        """Return stable colors for a sequence of category labels."""
        colors = categorical_colors(palette, n=len(categories))
        return dict(zip(categories, colors))

    @classmethod
    def reset(cls) -> None:
        """Restore Matplotlib default settings."""
        cls._active_preset = None
        cls._active_config = None
        if cls._original_rc is not None:
            plt.rcParams.update(cls._original_rc)
            cls._original_rc = None
        else:
            plt.rcdefaults()


# =============================================================================
# Convenience helpers.
# =============================================================================

def set_theme(style: str) -> None:
    """Apply a figure style to the process-wide theme."""
    PlotTheme.apply(style)


def get_colors(palette: Optional[str] = None, n: Optional[int] = None) -> List[str]:
    """Return a requested number of colors from the active palette."""
    return PlotTheme.get_palette(palette, n)


# =============================================================================
# Module initialization.
# =============================================================================

# Register the default style without mutating unrelated Matplotlib settings.
plt.rcParams["axes.prop_cycle"] = plt.cycler(color=PALETTES["default"])
