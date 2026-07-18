"""Minimal bridge to the official Bioconductor fgsea and GSVA implementations."""

from __future__ import annotations

import logging
import ntpath
import os
import shutil
import subprocess
import tempfile
from pathlib import Path
from typing import Dict, Iterable, Mapping, Optional, Sequence, Tuple

import numpy as np
import pandas as pd


logger = logging.getLogger(__name__)

R_SCRIPTS_DIR = Path(__file__).parent / "r_scripts"
FGSEA_COLUMNS = [
    "pathway",
    "pval",
    "padj",
    "log2err",
    "ES",
    "NES",
    "size",
    "leadingEdge",
]

_WINDOWS_R_RUNTIME_FAILURE_MARKERS = (
    "package 'grdevices' in options(\"defaultpackages\") was not found",
    "unable to load shared object",
    "loadlibrary failure",
    "application control policy has blocked this file",
    "Application control policy has prevented this file",
)


def find_rscript() -> str:
    """Return the configured Rscript executable or fail with an actionable error."""
    configured = os.environ.get("ALLENRICHER_RSCRIPT")
    executable = configured or shutil.which("Rscript")
    if not executable:
        raise RuntimeError(
            "Rscript was not found. Install R and add Rscript to PATH, "
            "or set ALLENRICHER_RSCRIPT to the Rscript executable."
        )
    return executable


def windows_to_wsl_path(value: object) -> str:
    """Translate an absolute Windows path for WSL interop; leave other values intact."""
    text = str(value)
    if text.startswith("\\\\?\\"):
        text = text[4:]
    drive, tail = ntpath.splitdrive(text)
    if len(drive) == 2 and drive[1] == ":":
        normalized_tail = tail.lstrip("\\/").replace("\\", "/")
        return f"/mnt/{drive[0].lower()}/{normalized_tail}"
    return text


def build_wsl_r_command(script_path: object, arguments: Sequence[object]) -> list[str]:
    """Build a WSL Rscript command for Windows hosts where native R cannot launch."""
    wsl = shutil.which("wsl.exe") or shutil.which("wsl")
    if not wsl:
        raise RuntimeError("Native Windows R could not be started, and no WSL R runtime is available")
    return [
        wsl,
        "Rscript",
        windows_to_wsl_path(script_path),
        *(windows_to_wsl_path(value) for value in arguments),
    ]


def _r_environment() -> Dict[str, str]:
    env = os.environ.copy()
    if os.name == "nt":
        for key in ("LANG", "LC_ALL", "LC_CTYPE"):
            if env.get(key) == "C.UTF-8":
                env.pop(key)
    return env


def should_retry_with_wsl(result: subprocess.CompletedProcess) -> bool:
    """Detect a blocked/broken native Windows R runtime, not script-level errors."""
    if os.name != "nt" or result.returncode == 0:
        return False
    details = f"{result.stdout}\n{result.stderr}".lower()
    return any(marker in details for marker in _WINDOWS_R_RUNTIME_FAILURE_MARKERS)


def _run_r_script(script_name: str, arguments: Sequence[object], timeout: int = 900) -> None:
    script_path = R_SCRIPTS_DIR / script_name
    command = [find_rscript(), str(script_path), *(str(value) for value in arguments)]
    run_options = dict(
        capture_output=True,
        text=True,
        encoding="utf-8",
        errors="replace",
        timeout=timeout,
        env=_r_environment(),
    )
    try:
        result = subprocess.run(command, **run_options)
    except OSError as error:
        if os.name != "nt":
            raise
        fallback = build_wsl_r_command(script_path, arguments)
        logger.info("Native Windows R could not start (%s); retrying with Rscript in WSL", error)
        result = subprocess.run(fallback, **run_options)
    if should_retry_with_wsl(result):
        fallback = build_wsl_r_command(script_path, arguments)
        logger.info("Native Windows R is unavailable; retrying with Rscript in WSL")
        result = subprocess.run(fallback, **run_options)
    if result.returncode != 0:
        details = "\n".join(part.strip() for part in (result.stdout, result.stderr) if part.strip())
        raise RuntimeError(f"Bioconductor analysis failed ({script_name}): \n{details}")
    if result.stdout.strip():
        logger.debug("R stdout (%s):\n%s", script_name, result.stdout.strip())
    if result.stderr.strip():
        logger.debug("R stderr (%s):\n%s", script_name, result.stderr.strip())


def _write_gmt(path: Path, gene_sets: Mapping[str, Iterable[str]]) -> None:
    with path.open("w", encoding="utf-8", newline="\n") as handle:
        for pathway, genes in gene_sets.items():
            pathway = str(pathway)
            if any(char in pathway for char in "\t\r\n"):
                raise ValueError(f"Gene-set name contains tabs or line breaks: {pathway!r}")
            clean_genes = sorted({str(gene).strip() for gene in genes if str(gene).strip()})
            if clean_genes:
                handle.write("\t".join([pathway, pathway, *clean_genes]) + "\n")


def run_fgsea(
    ranked_genes: Sequence[Tuple[str, float]],
    gene_sets: Mapping[str, Iterable[str]],
    min_size: int = 15,
    max_size: int = 500,
    seed: int = 42,
) -> pd.DataFrame:
    """Run fgseaMultilevel and return the unmodified official result columns."""
    if not ranked_genes:
        raise ValueError("GSEA requires non-empty sequenced list of genes")

    genes = [str(gene).strip() for gene, _ in ranked_genes]
    weights = np.asarray([weight for _, weight in ranked_genes], dtype=float)
    if any(not gene for gene in genes):
        raise ValueError("The GSEA ranked gene list contains an empty gene ID")
    if len(genes) != len(set(genes)):
        raise ValueError("The GSEA ranked gene list contains duplicate gene IDs")
    if not np.isfinite(weights).all():
        raise ValueError("GSEA Sort weight contains NaN or infinite value")
    if np.unique(weights).size < 2:
        raise ValueError("GSEA ranking weights must contain at least two distinct values")
    if not gene_sets:
        return pd.DataFrame(columns=FGSEA_COLUMNS)

    with tempfile.TemporaryDirectory(prefix="allenricher_fgsea_") as temp_dir:
        temp_path = Path(temp_dir)
        ranking_path = temp_path / "ranking.tsv"
        gmt_path = temp_path / "gene_sets.gmt"
        output_path = temp_path / "fgsea.tsv"
        pd.DataFrame({"gene": genes, "weight": weights}).to_csv(
            ranking_path, sep="\t", index=False, lineterminator="\n"
        )
        _write_gmt(gmt_path, gene_sets)
        _run_r_script(
            "fgsea_analysis.R",
            [ranking_path, gmt_path, output_path, min_size, max_size, seed],
        )
        result = pd.read_csv(output_path, sep="\t")

    missing = [column for column in FGSEA_COLUMNS if column not in result.columns]
    if missing:
        raise RuntimeError(f"fgsea output is missing: {', '.join(missing)}")
    return result.loc[:, FGSEA_COLUMNS]


def run_gsva(
    expression_matrix: pd.DataFrame,
    gene_sets: Mapping[str, Iterable[str]],
    method: str,
    kcdf: str = "Gaussian",
    tau: float = 1.0,
    min_size: int = 1,
    max_size: Optional[int] = None,
) -> pd.DataFrame:
    """Run an official GSVA-family method and return its pathway-by-sample matrix."""
    if expression_matrix.empty:
        return pd.DataFrame()
    if expression_matrix.index.has_duplicates:
        raise ValueError("The expression matrix contains duplicated gene IDs")
    if expression_matrix.columns.has_duplicates:
        raise ValueError("The expression matrix contains duplicate sample names")
    numeric = expression_matrix.apply(pd.to_numeric, errors="raise")
    if not np.isfinite(numeric.to_numpy(dtype=float)).all():
        raise ValueError("Expression matrix contains NaN or infinite value")
    if method not in {"gsva", "ssgsea", "plage", "zscore"}:
        raise ValueError(f"Unsupported GSVA method: {method}")
    if not gene_sets:
        return pd.DataFrame(columns=numeric.columns)

    available_genes = set(map(str, numeric.index))
    effective_max_size = max_size if max_size is not None else max(len(available_genes), 1)
    filtered_gene_sets = {
        str(pathway): {str(gene) for gene in genes} & available_genes
        for pathway, genes in gene_sets.items()
    }
    filtered_gene_sets = {
        pathway: genes
        for pathway, genes in filtered_gene_sets.items()
        if min_size <= len(genes) <= effective_max_size
    }
    if not filtered_gene_sets:
        empty = pd.DataFrame(columns=numeric.columns, dtype=float)
        empty.index.name = "Pathway"
        return empty

    with tempfile.TemporaryDirectory(prefix="allenricher_gsva_") as temp_dir:
        temp_path = Path(temp_dir)
        expression_path = temp_path / "expression.tsv"
        gmt_path = temp_path / "gene_sets.gmt"
        output_path = temp_path / "activity.tsv"
        numeric.to_csv(expression_path, sep="\t", index=True, index_label="gene", lineterminator="\n")
        _write_gmt(gmt_path, filtered_gene_sets)
        _run_r_script(
            "gsva_analysis.R",
            [
                expression_path,
                gmt_path,
                output_path,
                method,
                kcdf,
                tau,
                min_size,
                effective_max_size,
            ],
        )
        result = pd.read_csv(output_path, sep="\t", index_col=0)

    result.index.name = "Pathway"
    return result
