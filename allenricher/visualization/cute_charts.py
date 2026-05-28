"""cutecharts 手绘风格适配层

将 ORA 富集分析的 barplot/bubble 映射到 cutecharts 手绘风格，
输出为交互式 HTML 文件。cutecharts 为可选依赖。
"""
from __future__ import annotations

from pathlib import Path
from typing import Dict, List, Optional


def _check_cutecharts() -> bool:
    """检查 cutecharts 是否已安装"""
    try:
        import cutecharts  # noqa: F401
        return True
    except ImportError:
        return False


def plot_cute_barplot(
    data: List[Dict],
    output_file: str,
    db_name: str = 'GO',
    top_n: int = 20,
) -> Path:
    """cutecharts 手绘风格水平柱状图

    Args:
        data: 富集结果列表
        output_file: 输出 HTML 文件路径
        db_name: 数据库名称
        top_n: 最大条目数

    Returns:
        输出文件路径
    """
    if not _check_cutecharts():
        raise ImportError(
            "cutecharts 未安装。请运行: pip install cutecharts"
        )

    from cutecharts.charts import Bar

    if not data:
        return Path(output_file)

    df = data[:top_n]

    # 提取数据（取最后 top_n 条，因为 cutecharts bar 是垂直的，从下到上）
    names = [d.get('Term_Name', d.get('term_name', ''))[-30:] for d in reversed(df)]
    neg_log_q = [
        -round(max(0, __import__('math').log10(
            max(d.get('Adjusted_P_Value', d.get('adjusted_p_value', 1)), 1e-300)
        )), 2)
        for d in reversed(df)
    ]

    chart = Bar(f"{db_name} Enrichment (Top {len(names)})")
    chart.set_options(
        labels=names,
        x_label="-log10(Q-value)",
        y_label="",
        colors=[
            "#EE6677", "#4477AA", "#228833", "#CCBB44",
            "#66CCEE", "#AA3377", "#BBBBBB", "#EE7733",
        ],
    )
    chart.add_series("-log10(Q-value)", neg_log_q)

    # 确保输出为 .html
    out = Path(output_file)
    if out.suffix != '.html':
        out = out.with_suffix('.html')

    chart.render(str(out))
    return out


def plot_cute_bubble(
    data: List[Dict],
    output_file: str,
    db_name: str = 'GO',
    top_n: int = 20,
) -> Path:
    """cutecharts 手绘风格散点图（近似气泡图）

    注意：cutecharts 没有原生气泡图（大小映射），此处用 Scatter 近似。
    点大小固定，颜色区分 -log10(Q-value) 区间。

    Args:
        data: 富集结果列表
        output_file: 输出 HTML 文件路径
        db_name: 数据库名称
        top_n: 最大条目数

    Returns:
        输出文件路径
    """
    if not _check_cutecharts():
        raise ImportError(
            "cutecharts 未安装。请运行: pip install cutecharts"
        )

    from cutecharts.charts import Scatter
    import math

    if not data:
        return Path(output_file)

    df = data[:top_n]

    names = [d.get('Term_Name', d.get('term_name', ''))[-25:] for d in df]
    gene_counts = [d.get('Gene_Count', d.get('gene_count', 1)) for d in df]
    bg_counts = [d.get('Background_Count', d.get('background_count', 1)) for d in df]
    pvalues = [d.get('Adjusted_P_Value', d.get('adjusted_p_value', 1)) for d in df]

    rich_factors = []
    for gc, bc in zip(gene_counts, bg_counts):
        if bc and bc > 0:
            rich_factors.append(round(gc / bc, 4))
        else:
            rich_factors.append(0)

    neg_log_q = [
        -round(max(0, math.log10(max(p, 1e-300))), 2)
        for p in pvalues
    ]

    chart = Scatter(f"{db_name} Enrichment (Top {len(names)})")
    chart.set_options(
        x_label="Rich Factor",
        y_label="-log10(Q-value)",
        dot_size=3,
    )
    chart.add_series(
        "Enrichment",
        [(rf, nlq) for rf, nlq in zip(rich_factors, neg_log_q)],
    )

    out = Path(output_file)
    if out.suffix != '.html':
        out = out.with_suffix('.html')

    chart.render(str(out))
    return out
