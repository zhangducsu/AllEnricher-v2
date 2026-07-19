"""Validate the R environment and invoke maintained R figure scripts."""
import logging
import os
import subprocess
import shutil
from pathlib import Path
from typing import Dict, List, Optional

from allenricher.core.bioconductor import build_wsl_r_command, should_retry_with_wsl

from .color_config import PaletteLike, coerce_palette_selection
from .plot_theme import resolve_style

logger = logging.getLogger(__name__)

R_SCRIPTS_DIR = Path(__file__).parent / "r_scripts"


def check_r_environment() -> bool:
    """Return whether Rscript and all required plotting packages are available."""
    return shutil.which("Rscript") is not None or (
        os.name == "nt" and (shutil.which("wsl.exe") is not None or shutil.which("wsl") is not None)
    )


def run_r_script(
    script_name: str,
    args: Dict[str, str],
    output_file: str,
    timeout: int = 300,
) -> bool:
    """Invoke one maintained R script with validated arguments and captured logs."""
    script_path = R_SCRIPTS_DIR / script_name
    if not script_path.exists():
        logger.error(f"R script not found: {script_path}")
        return False

    cmd = ["Rscript", str(script_path)]
    # Resolve only arguments that represent files. Values such as gene_set_ids
    # are comma-separated identifiers and must remain unchanged.
    _PATH_ARGS = {"tsv", "expr", "running_es", "scores", "metadata", "statistics"}
    for key, value in args.items():
        if key in _PATH_ARGS:
            cmd.extend([f"--{key}", str(Path(value).resolve())])
        else:
            cmd.extend([f"--{key}", str(value)])
    resolved_output = Path(output_file).resolve()
    cmd.extend(["--output", str(resolved_output)])

    logger.info(f"Running: {' '.join(cmd)}")
    env = os.environ.copy()
    if os.name == "nt":
        for key in ("LANG", "LC_ALL", "LC_CTYPE"):
            if env.get(key) == "C.UTF-8":
                env.pop(key)
    try:
        run_options = dict(
            capture_output=True, text=True, encoding="utf-8", errors="replace",
            timeout=timeout, cwd=str(R_SCRIPTS_DIR.parent), env=env,
        )
        try:
            result = subprocess.run(cmd, **run_options)
        except OSError as error:
            if os.name != "nt":
                raise
            fallback = build_wsl_r_command(script_path, cmd[2:])
            logger.info("Native Windows R could not start (%s); retrying with Rscript in WSL", error)
            result = subprocess.run(fallback, **run_options)
        if should_retry_with_wsl(result):
            fallback = build_wsl_r_command(script_path, cmd[2:])
            logger.info("Native Windows R is unavailable; retrying with Rscript in WSL")
            result = subprocess.run(fallback, **run_options)
        if result.returncode != 0:
            logger.error(
                "R script failed: %s\nstdout:\n%s\nstderr:\n%s",
                script_name,
                result.stdout,
                result.stderr,
            )
            return False
        if result.stdout:
            logger.debug("R script stdout (%s):\n%s", script_name, result.stdout)
        if result.stderr:
            logger.debug("R script stderr (%s):\n%s", script_name, result.stderr)
        if not resolved_output.is_file() or resolved_output.stat().st_size == 0:
            logger.error("R script exited successfully but did not create a non-empty output: %s", resolved_output)
            return False
        logger.info(f"R plot saved: {output_file}")
        return True
    except subprocess.TimeoutExpired:
        logger.error(f"R script timed out after {timeout}s")
        return False
    except Exception as e:
        logger.error(f"R script error: {e}")
        return False


# Shared display options passed to each maintained R plot.
def _style_args(style: str, palette: PaletteLike) -> Dict[str, str]:
    selection = coerce_palette_selection(palette)
    return {
        "style": resolve_style(style),
        "categorical_palette": selection.categorical,
        "sequential_palette": selection.sequential,
        "diverging_palette": selection.diverging,
    }


def plot_gsea_barplot_r(
    tsv_path: str, output_file: str, top_n: int = 20, dpi: int = 300,
    style: str = "nature", palette: PaletteLike = None,
) -> bool:
    args = {"tsv": tsv_path, "top_n": str(top_n), "dpi": str(dpi)}
    args.update(_style_args(style, palette))
    return run_r_script(
        "gsea_barplot.R", args, output_file
    )

def plot_gsea_lollipop_r(
    tsv_path: str, output_file: str, top_n: int = 20, dpi: int = 300,
    style: str = "nature", palette: PaletteLike = None,
) -> bool:
    args = {"tsv": tsv_path, "top_n": str(top_n), "dpi": str(dpi)}
    args.update(_style_args(style, palette))
    return run_r_script(
        "gsea_lollipop.R", args, output_file
    )

def plot_gsea_ridgeplot_r(
    tsv_path: str,
    output_file: str,
    top_n: int = 10,
    running_es_path: str = "",
    dpi: int = 300,
    style: str = "nature",
    palette: PaletteLike = None,
) -> bool:
    args = {"tsv": tsv_path, "top_n": str(top_n), "dpi": str(dpi)}
    args.update(_style_args(style, palette))
    if running_es_path:
        args["running_es"] = running_es_path
    return run_r_script("gsea_ridgeplot.R", args, output_file)

def plot_gsea_emapplot_r(
    tsv_path: str,
    output_file: str,
    top_n: int = 30,
    qvalue: float = 0.05,
    min_count: int = 3,
    dpi: int = 300,
    style: str = "nature",
    palette: PaletteLike = None,
) -> bool:
    args = {
        "tsv": tsv_path,
        "top_n": str(top_n),
        "qvalue": str(qvalue),
        "min_count": str(min_count),
        "dpi": str(dpi),
    }
    args.update(_style_args(style, palette))
    return run_r_script(
        "gsea_emapplot.R",
        args,
        output_file,
    )

def plot_gsea_enrichment_r(
    tsv_path: str, gene_set_id: str, output_file: str,
    running_es_path: str = "", dpi: int = 300,
    style: str = "nature", palette: PaletteLike = None,
) -> bool:
    args = {"tsv": tsv_path, "gene_set_id": gene_set_id, "dpi": str(dpi)}
    args.update(_style_args(style, palette))
    if running_es_path:
        args["running_es"] = running_es_path
    return run_r_script("gsea_enrichment_plot.R", args, output_file)

def plot_gsea_enrichment2_r(
    tsv_path: str, gene_set_ids: List[str], output_file: str,
    running_es_path: str = "", dpi: int = 300,
    style: str = "nature", palette: PaletteLike = None,
) -> bool:
    args = {"tsv": tsv_path, "gene_set_ids": ",".join(gene_set_ids), "dpi": str(dpi)}
    args.update(_style_args(style, palette))
    if running_es_path:
        args["running_es"] = running_es_path
    return run_r_script("gsea_enrichment_plot2.R", args, output_file)


def plot_activity_heatmap_r(
    scores_path: str,
    metadata_path: str,
    output_file: str,
    analysis_method: str = "",
    scale: str = "row",
    top_n: int = 40,
    dpi: int = 300,
    style: str = "nature",
    palette: PaletteLike = None,
) -> bool:
    args = {
        "scores": scores_path, "metadata": metadata_path,
        "analysis_method": analysis_method, "scale": scale,
        "top_n": str(top_n), "dpi": str(dpi),
    }
    args.update(_style_args(style, palette))
    return run_r_script(
        "activity_heatmap.R",
        args,
        output_file,
    )


def plot_sample_correlation_r(
    scores_path: str,
    metadata_path: str,
    output_file: str,
    method: str = "pearson",
    dpi: int = 300,
    style: str = "nature",
    palette: PaletteLike = None,
) -> bool:
    args = {"scores": scores_path, "metadata": metadata_path, "method": method, "dpi": str(dpi)}
    args.update(_style_args(style, palette))
    return run_r_script(
        "sample_correlation.R",
        args,
        output_file,
    )


def plot_group_comparison_r(
    scores_path: str,
    metadata_path: str,
    output_file: str,
    statistics_file: str,
    top_n: int = 6,
    dpi: int = 300,
    style: str = "nature",
    palette: PaletteLike = None,
) -> bool:
    args = {
        "scores": scores_path,
        "metadata": metadata_path,
        "statistics": statistics_file,
        "top_n": str(top_n),
        "dpi": str(dpi),
    }
    args.update(_style_args(style, palette))
    generated = run_r_script(
        "group_comparison.R",
        args,
        output_file,
    )
    if generated and not Path(statistics_file).is_file():
        logger.error("R group comparison did not create statistics table: %s", statistics_file)
        return False
    return generated

# R plot types supported by the Python-to-R bridge.
R_PLOT_TYPES = [
    "barplot", "lollipop", "ridgeplot", "emapplot", "enrichment", "enrichment2",
]

R_PLOT_FUNC_MAP = {
    "barplot": plot_gsea_barplot_r,
    "lollipop": plot_gsea_lollipop_r,
    "ridgeplot": plot_gsea_ridgeplot_r,
    "emapplot": plot_gsea_emapplot_r,
    "enrichment": plot_gsea_enrichment_r,
    "enrichment2": plot_gsea_enrichment2_r,
}
