"""R 脚本调用层 — 通过 subprocess 调用 R 生成发表级图表"""
import logging
import subprocess
import shutil
from pathlib import Path
from typing import Dict, List

logger = logging.getLogger(__name__)

R_SCRIPTS_DIR = Path(__file__).parent / "r_scripts"


def check_r_environment() -> bool:
    """检测 R 环境是否可用"""
    return shutil.which("Rscript") is not None


def run_r_script(
    script_name: str,
    args: Dict[str, str],
    output_file: str,
    timeout: int = 300,
) -> bool:
    """运行 R 脚本生成图表

    Args:
        script_name: R 脚本文件名
        args: 传递给 R 脚本的参数
        output_file: 输出图表路径
        timeout: 超时时间（秒）
    Returns:
        bool: 是否成功
    """
    script_path = R_SCRIPTS_DIR / script_name
    if not script_path.exists():
        logger.error(f"R script not found: {script_path}")
        return False

    cmd = ["Rscript", str(script_path)]
    # 只解析真实文件路径参数；gene_set_ids 是逗号分隔的 ID 字符串，不能当路径处理。
    _PATH_ARGS = {"tsv", "expr", "running_es"}
    for key, value in args.items():
        if key in _PATH_ARGS:
            cmd.extend([f"--{key}", str(Path(value).resolve())])
        else:
            cmd.extend([f"--{key}", str(value)])
    cmd.extend(["--output", str(Path(output_file).resolve())])

    logger.info(f"Running: {' '.join(cmd)}")
    try:
        result = subprocess.run(
            cmd, capture_output=True, text=True, timeout=timeout,
            cwd=str(R_SCRIPTS_DIR.parent)
        )
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
        logger.info(f"R plot saved: {output_file}")
        return True
    except subprocess.TimeoutExpired:
        logger.error(f"R script timed out after {timeout}s")
        return False
    except Exception as e:
        logger.error(f"R script error: {e}")
        return False


# 便捷函数
def plot_gsea_dotplot_r(tsv_path: str, output_file: str, top_n: int = 20) -> bool:
    return run_r_script("gsea_dotplot.R", {"tsv": tsv_path, "top_n": str(top_n)}, output_file)

def plot_gsea_barplot_r(tsv_path: str, output_file: str, top_n: int = 20) -> bool:
    return run_r_script("gsea_barplot.R", {"tsv": tsv_path, "top_n": str(top_n)}, output_file)

def plot_gsea_nes_plot_r(tsv_path: str, output_file: str, top_n: int = 30) -> bool:
    return run_r_script("gsea_nes_plot.R", {"tsv": tsv_path, "top_n": str(top_n)}, output_file)

def plot_gsea_ridgeplot_r(
    tsv_path: str,
    output_file: str,
    top_n: int = 15,
    running_es_path: str = "",
) -> bool:
    args = {"tsv": tsv_path, "top_n": str(top_n)}
    if running_es_path:
        args["running_es"] = running_es_path
    return run_r_script("gsea_ridgeplot.R", args, output_file)

def plot_gsea_heatmap_r(
    expr_path: str,
    output_file: str,
    tsv_path: str = "",
    top_n: int = 12,
) -> bool:
    args = {"expr": expr_path, "top_n": str(top_n)}
    if tsv_path:
        args["tsv"] = tsv_path
    return run_r_script("gsea_heatmap.R", args, output_file)

def plot_gsea_emapplot_r(tsv_path: str, output_file: str, top_n: int = 30) -> bool:
    return run_r_script("gsea_emapplot.R", {"tsv": tsv_path, "top_n": str(top_n)}, output_file)

def plot_gsea_cnetplot_r(tsv_path: str, output_file: str, top_n: int = 10) -> bool:
    return run_r_script("gsea_cnetplot.R", {"tsv": tsv_path, "top_n": str(top_n)}, output_file)

def plot_gsea_circos_r(tsv_path: str, output_file: str, top_n: int = 30) -> bool:
    return run_r_script("gsea_circos.R", {"tsv": tsv_path, "top_n": str(top_n)}, output_file)

def plot_gsea_enrichment_r(tsv_path: str, gene_set_id: str, output_file: str, running_es_path: str = "") -> bool:
    args = {"tsv": tsv_path, "gene_set_id": gene_set_id}
    if running_es_path:
        args["running_es"] = running_es_path
    return run_r_script("gsea_enrichment_plot.R", args, output_file)

def plot_gsea_enrichment2_r(tsv_path: str, gene_set_ids: List[str], output_file: str, running_es_path: str = "") -> bool:
    args = {"tsv": tsv_path, "gene_set_ids": ",".join(gene_set_ids)}
    if running_es_path:
        args["running_es"] = running_es_path
    return run_r_script("gsea_enrichment_plot2.R", args, output_file)

# 全部 R 图表类型
R_PLOT_TYPES = [
    "dotplot", "barplot", "nes_plot", "ridgeplot", "heatmap",
    "emapplot", "cnetplot", "circos", "enrichment", "enrichment2",
]

R_PLOT_FUNC_MAP = {
    "dotplot": plot_gsea_dotplot_r,
    "barplot": plot_gsea_barplot_r,
    "nes_plot": plot_gsea_nes_plot_r,
    "ridgeplot": plot_gsea_ridgeplot_r,
    "heatmap": plot_gsea_heatmap_r,
    "emapplot": plot_gsea_emapplot_r,
    "cnetplot": plot_gsea_cnetplot_r,
    "circos": plot_gsea_circos_r,
    "enrichment": plot_gsea_enrichment_r,
    "enrichment2": plot_gsea_enrichment2_r,
}
