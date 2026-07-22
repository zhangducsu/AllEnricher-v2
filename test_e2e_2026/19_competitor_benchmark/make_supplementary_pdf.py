#!/usr/bin/env python3
"""Build the single-file AllEnricher Supplementary Data PDF."""

from __future__ import annotations

import argparse
import csv
from html import escape
from pathlib import Path
from typing import Iterable

from PIL import Image as PILImage
from reportlab.lib import colors
from reportlab.lib.enums import TA_CENTER
from reportlab.lib.pagesizes import A4, landscape
from reportlab.lib.styles import ParagraphStyle, getSampleStyleSheet
from reportlab.lib.units import mm
from reportlab.platypus import (
    Image, LongTable, PageBreak, Paragraph, SimpleDocTemplate, Spacer, Table, TableStyle,
)

PAGE = landscape(A4)
BLUE = colors.HexColor("#0072B2")
PURPLE = colors.HexColor("#7B61A8")
LIGHT_BLUE = colors.HexColor("#E8F3F8")
LIGHT_GREY = colors.HexColor("#F3F3F3")


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def paragraph(value: object, style: ParagraphStyle) -> Paragraph:
    text = "" if value is None else str(value)
    return Paragraph(escape(text).replace("\n", "<br/>"), style)


def styled_table(rows: list[list[object]], widths: list[float], repeat_rows: int = 1) -> LongTable:
    table = LongTable(rows, colWidths=widths, repeatRows=repeat_rows, hAlign="LEFT")
    table.setStyle(TableStyle([
        ("BACKGROUND", (0, 0), (-1, 0), LIGHT_BLUE),
        ("TEXTCOLOR", (0, 0), (-1, 0), colors.black),
        ("FONTNAME", (0, 0), (-1, 0), "Helvetica-Bold"),
        ("FONTSIZE", (0, 0), (-1, -1), 7.5),
        ("LEADING", (0, 0), (-1, -1), 9),
        ("VALIGN", (0, 0), (-1, -1), "TOP"),
        ("GRID", (0, 0), (-1, -1), 0.35, colors.HexColor("#A8A8A8")),
        ("ROWBACKGROUNDS", (0, 1), (-1, -1), [colors.white, LIGHT_GREY]),
        ("LEFTPADDING", (0, 0), (-1, -1), 3),
        ("RIGHTPADDING", (0, 0), (-1, -1), 3),
        ("TOPPADDING", (0, 0), (-1, -1), 3),
        ("BOTTOMPADDING", (0, 0), (-1, -1), 3),
    ]))
    return table


def figure(path: Path, max_width: float, max_height: float) -> Image:
    with PILImage.open(path) as image:
        width, height = image.size
    scale = min(max_width / width, max_height / height)
    return Image(str(path), width=width * scale, height=height * scale)


def page_header_footer(canvas, doc) -> None:
    canvas.saveState()
    canvas.setFont("Helvetica", 7.5)
    canvas.setFillColor(colors.HexColor("#555555"))
    canvas.drawString(15 * mm, PAGE[1] - 10 * mm, "AllEnricher v2 - Supplementary Data")
    canvas.drawRightString(PAGE[0] - 15 * mm, 9 * mm, f"Page {doc.page}")
    canvas.restoreState()


def table_s1(story: list[object], rows: list[dict[str, str]], styles: dict[str, ParagraphStyle]) -> None:
    story += [Paragraph("Table S1. Datasets, inputs, annotations and licences", styles["Heading2"])]
    core = [
        ["Dataset", "Accession", "Species", "DB", "Query", "Background", "Ranking", "Query mapping", "Background mapping"]
    ]
    for row in rows:
        core.append([
            row["dataset"], row["accession"], row["species"], row["database"], row["query_genes"],
            row["background_genes"], row["ranked_genes"], f"{float(row['query_mapping_rate']):.4f}",
            f"{float(row['background_mapping_rate']):.4f}",
        ])
    story += [styled_table(core, [42 * mm, 25 * mm, 17 * mm, 14 * mm, 17 * mm, 22 * mm, 19 * mm, 24 * mm, 28 * mm]), Spacer(1, 4 * mm)]
    sources = [["Dataset / DB", "Archive-relative source GMT", "SHA-256", "Expression Atlas licence", "Annotation distribution", "Snapshot"]]
    for row in rows:
        sources.append([
            paragraph(f"{row['dataset']} / {row['database']}", styles["Cell"]),
            paragraph(row["source_gmt"], styles["Cell"]),
            paragraph(row["source_gmt_sha256"], styles["Cell"]),
            row["expression_atlas_license"], row["annotation_distribution"], row["annotation_snapshot"],
        ])
    story += [styled_table(sources, [32 * mm, 59 * mm, 66 * mm, 30 * mm, 45 * mm, 27 * mm]), Spacer(1, 2 * mm)]
    story.append(Paragraph("Paths are archive-relative provenance identifiers. GO annotation files may be redistributed under CC BY 4.0. KEGG GMT files and gene-set membership are not included in the public archive; their frozen-source paths and hashes are retained for audit and reconstruction. Numerical input statistics are unchanged from the frozen run.", styles["Note"]))


def table_s2(story: list[object], rows: list[dict[str, str]], styles: dict[str, ParagraphStyle]) -> None:
    story += [PageBreak(), Paragraph("Table S2. Evidence-linked capability matrix", styles["Heading2"])]
    tools = list(dict.fromkeys(row["tool"] for row in rows))
    groups = list(dict.fromkeys(row["group"] for row in rows))
    lookup = {(row["feature"], row["tool"]): row["value"] for row in rows}
    for group in groups:
        features = list(dict.fromkeys(row["feature"] for row in rows if row["group"] == group))
        matrix = [["Feature", *tools]]
        matrix.extend([[feature, *(lookup[(feature, tool)] for tool in tools)] for feature in features])
        story += [Paragraph(group, styles["Heading3"]), styled_table(matrix, [64 * mm, *([39 * mm] * len(tools))]), Spacer(1, 3 * mm)]
    definitions = [["Group", "Feature", "Operational definition"]]
    for row in rows:
        if row["tool"] == tools[0]:
            definitions.append([
                row["group"], row["feature"], paragraph(row["definition"], styles["Cell"]),
            ])
    story += [PageBreak(), Paragraph("Capability definitions", styles["Heading3"])]
    story += [styled_table(definitions, [52 * mm, 55 * mm, 147 * mm]), Spacer(1, 4 * mm)]
    grouped_evidence: dict[tuple[str, str], list[str]] = {}
    for row in rows:
        grouped_evidence.setdefault((row["tool"], row["evidence"]), []).append(row["feature"])
    evidence = [["Tool", "Assessed features", "Evidence"]]
    for (tool, source), features in grouped_evidence.items():
        evidence.append([
            tool, paragraph(", ".join(features), styles["Cell"]), paragraph(source, styles["Cell"]),
        ])
    story += [Paragraph("Evidence sources", styles["Heading3"])]
    story += [styled_table(evidence, [31 * mm, 103 * mm, 120 * mm]), Spacer(1, 2 * mm)]
    story.append(Paragraph("Evidence entries document the basis for Yes, Partial and No assessments; listing a feature here does not imply support. Cells are restricted to Yes, Partial, No or N/A. No composite score is calculated.", styles["Note"]))


def table_s3(story: list[object], rows: list[dict[str, str]], styles: dict[str, ParagraphStyle]) -> None:
    story += [PageBreak(), Paragraph("Table S3. Tool versions and reproducibility evidence", styles["Heading2"])]
    summary = [["Tool", "Version", "Role", "Environment", "Access date", "Commit SHA", "Source URL"]]
    for row in rows:
        summary.append([
            row["tool"], paragraph(row["version"], styles["Cell"]), paragraph(row["role"], styles["Cell"]),
            paragraph(row["execution_environment"], styles["Cell"]), row["access_date"],
            paragraph(row["commit_sha"], styles["Cell"]), paragraph(row["source_url"], styles["Cell"]),
        ])
    story += [styled_table(summary, [24 * mm, 30 * mm, 47 * mm, 48 * mm, 25 * mm, 48 * mm, 45 * mm]), Spacer(1, 4 * mm)]
    detail = [["Tool", "Portable command template", "Evidence archive path"]]
    for row in rows:
        detail.append([
            row["tool"], paragraph(row["command_template"], styles["Cell"]),
            paragraph(row["evidence_archive_path"], styles["Cell"]),
        ])
    story += [styled_table(detail, [28 * mm, 155 * mm, 84 * mm]), Spacer(1, 2 * mm)]
    story.append(Paragraph("The 56 case-specific commands are provided as source_data/Data_S1_full_case_commands.tsv. {REPO_ROOT} and {ARCHIVE_ROOT} denote the released repository and evidence-archive roots.", styles["Note"]))


def compact_metrics(row: dict[str, str]) -> str:
    labels = [
        ("reference_terms", "reference terms"), ("comparator_terms", "comparator terms"),
        ("common_terms", "common terms"), ("term_jaccard", "term Jaccard"),
        ("spearman", "NES/q Spearman"), ("p_spearman", "P-value Spearman"),
        ("max_abs_p_diff", "max |delta P|"), ("median_abs_p_diff", "median |delta P|"),
        ("max_abs_q_diff", "max |delta q|"), ("median_abs_q_diff", "median |delta q|"),
        ("valid_p_pairs", "valid P pairs"), ("valid_q_pairs", "valid q pairs"),
        ("valid_nes_pairs", "valid NES pairs"), ("max_abs_nes_diff", "max |delta NES|"),
        ("median_abs_nes_diff", "median |delta NES|"), ("sign_concordance", "direction"),
        ("significant_jaccard", "significant Jaccard"), ("positive_significant_jaccard", "positive Jaccard"),
        ("negative_significant_jaccard", "negative Jaccard"), ("top20_jaccard", "top-20 Jaccard"),
        ("leading_edge_jaccard", "leading-edge Jaccard"),
    ]
    values = []
    for key, label in labels:
        value = row.get(key, "")
        if value not in ("", None, "nan"):
            values.append(f"<b>{escape(label)}:</b> {escape(str(value))}")
    return "; ".join(values) if values else "No numeric comparison metrics"


def table_s4(story: list[object], rows: list[dict[str, str]], styles: dict[str, ParagraphStyle]) -> None:
    story += [PageBreak(), Paragraph("Table S4. All case-level metrics and failure records", styles["Heading2"])]
    data = [["Case", "Predefined metrics", "Status and retained reason"]]
    for row in rows:
        comparator = row.get("comparator") or row.get("tool") or "N/A"
        case = "<b>{}</b><br/>{} / {} / {}<br/>{}".format(
            escape(comparator), escape(row.get("dataset", "")), escape(row.get("database", "")),
            escape(row.get("method", "")), escape(row.get("record_type", "")),
        )
        status = f"<b>{escape(row.get('status', ''))}</b>"
        reason = row.get("reason", "")
        if reason:
            status += f"<br/>{escape(reason)}"
        data.append([
            Paragraph(case, styles["Cell"]), Paragraph(compact_metrics(row), styles["Cell"]),
            Paragraph(status, styles["Cell"]),
        ])
    story += [styled_table(data, [52 * mm, 139 * mm, 76 * mm]), Spacer(1, 2 * mm)]
    story.append(Paragraph("Blank fields are not applicable to that record type. GSEA P-value comparisons are descriptive across executed engines. Maximum adjusted-P differences can be amplified by BH ranks over large testing families and should be interpreted with median differences, valid-pair counts, NES ranks, directions and set-overlap metrics. Expression Atlas rankings were retained without artificial tie-breaking. Failures, missing outputs and semantic incompatibilities are retained rather than removed or imputed.", styles["Note"]))


def build(paper_dir: Path, output: Path) -> None:
    supplementary = paper_dir / "supplementary"
    styles0 = getSampleStyleSheet()
    styles = {
        "Title": ParagraphStyle("Title", parent=styles0["Title"], fontName="Helvetica-Bold", fontSize=20, leading=24, textColor=BLUE, alignment=TA_CENTER, spaceAfter=8 * mm),
        "Heading2": ParagraphStyle("Heading2", parent=styles0["Heading2"], fontName="Helvetica-Bold", fontSize=12, leading=15, textColor=BLUE, spaceAfter=4 * mm),
        "Heading3": ParagraphStyle("Heading3", parent=styles0["Heading3"], fontName="Helvetica-Bold", fontSize=9, leading=11, textColor=PURPLE, spaceBefore=2 * mm, spaceAfter=2 * mm),
        "Body": ParagraphStyle("Body", parent=styles0["BodyText"], fontName="Helvetica", fontSize=8, leading=10, spaceAfter=3 * mm),
        "Cell": ParagraphStyle("Cell", parent=styles0["BodyText"], fontName="Helvetica", fontSize=7.5, leading=9),
        "Note": ParagraphStyle("Note", parent=styles0["BodyText"], fontName="Helvetica-Oblique", fontSize=7.5, leading=9, textColor=colors.HexColor("#444444")),
        "Caption": ParagraphStyle("Caption", parent=styles0["BodyText"], fontName="Helvetica", fontSize=8, leading=10, spaceBefore=2 * mm, spaceAfter=2 * mm),
    }
    output.parent.mkdir(parents=True, exist_ok=True)
    doc = SimpleDocTemplate(
        str(output), pagesize=PAGE, rightMargin=14 * mm, leftMargin=14 * mm,
        topMargin=16 * mm, bottomMargin=14 * mm, title="AllEnricher v2 Supplementary Data",
        author="Jieling Xiao; Qianwen Zhou; Chunyu Wang; Du Zhang", subject="Figures S1-S2 and Tables S1-S4",
    )
    story: list[object] = [Spacer(1, 12 * mm), Paragraph("AllEnricher v2", styles["Title"]), Paragraph("Supplementary Data", styles["Title"])]
    story += [
        Paragraph("This document contains Figures S1-S2 and Tables S1-S4. Machine-readable TSV files and the complete 56-case command ledger must accompany the stable Figshare evidence archive before submission (FIGSHARE_DOI_PENDING). Paths in publication tables are relative to that archive root.", styles["Body"]),
        Paragraph("Internal quality-control files, local run directories, database logos and internal design references are not part of the public evidence archive.", styles["Body"]),
        PageBreak(),
    ]
    figures = [
        ("Figure S1. Term-level ORA agreement across four species", "Figure_S1_ORA_term_agreement.png", "Dashed lines indicate equality. AllEnricher and clusterProfiler use the same positive-overlap BH family; WebGestaltR and g:Profiler retain their documented custom-annotation universe semantics.", "Multiple small scatter plots compare adjusted ORA values for human, cattle, fly and yeast GO and KEGG cases."),
        ("Figure S2. Term-level GSEA agreement and case metrics", "Figure_S2_GSEA_term_agreement.png", "Scatter plots compare NES values; the heat map summarizes predefined case-level metrics. Nominal P-value and BH-adjusted-P comparisons, including valid-pair denominators, are reported in Table S4.", "Scatter plots across four species are followed by a heat map of NES rank, direction and significant-set agreement."),
    ]
    for index, (title, filename, caption, alt) in enumerate(figures):
        story += [Paragraph(title, styles["Heading2"]), figure(supplementary / filename, 260 * mm, 155 * mm), Paragraph(caption, styles["Caption"]), Paragraph(f"Alt text: {alt}", styles["Note"])]
        if index < len(figures) - 1:
            story.append(PageBreak())
    table_s1(story, read_tsv(supplementary / "Table_S1_datasets_inputs_databases.tsv"), styles)
    table_s2(story, read_tsv(supplementary / "Table_S2_capability_evidence.tsv"), styles)
    table_s3(story, read_tsv(supplementary / "Table_S3_versions_commands_access.tsv"), styles)
    table_s4(story, read_tsv(supplementary / "Table_S4_case_metrics_failures.tsv"), styles)
    doc.build(story, onFirstPage=page_header_footer, onLaterPages=page_header_footer)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--paper-dir", type=Path, required=True)
    parser.add_argument("--output", type=Path)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    paper_dir = args.paper_dir.resolve()
    output = args.output.resolve() if args.output else paper_dir / "oup_submission" / "Supplementary_Data.pdf"
    build(paper_dir, output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())