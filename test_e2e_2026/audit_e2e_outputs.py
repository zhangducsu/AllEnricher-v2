#!/usr/bin/env python3
"""A light, re-recoverable machine audit of E2E charts and logs."""

from __future__ import annotations

import argparse
import json
import re
import xml.etree.ElementTree as ET
from pathlib import Path

from PIL import Image, ImageChops


def audit_png(path: Path) -> tuple[dict, list[dict]]:
    issues = []
    metrics = {"path": str(path), "size": path.stat().st_size}
    try:
        with Image.open(path) as image:
            image.load()
            rgb = image.convert("RGB")
            metrics["dimensions"] = list(rgb.size)
            aspect_ratio = max(rgb.width / rgb.height, rgb.height / rgb.width)
            metrics["aspect_ratio"] = round(aspect_ratio, 3)
            if rgb.width < 300 or rgb.height < 200:
                issues.append({"type": "small_image", "path": str(path), "detail": f"{rgb.width}x{rgb.height}"})
            if aspect_ratio > 4:
                issues.append({"type": "extreme_aspect_ratio", "path": str(path), "detail": round(aspect_ratio, 3)})

            sample = rgb.copy()
            sample.thumbnail((256, 256))
            background = Image.new("RGB", sample.size, sample.getpixel((0, 0)))
            bbox = ImageChops.difference(sample, background).getbbox()
            metrics["content_bbox"] = list(bbox) if bbox else None
            if not bbox:
                issues.append({"type": "blank_image", "path": str(path)})
            else:
                content_ratio = ((bbox[2] - bbox[0]) * (bbox[3] - bbox[1])) / (sample.width * sample.height)
                metrics["content_bbox_ratio"] = round(content_ratio, 4)
                if content_ratio < 0.1:
                    issues.append({"type": "excessive_whitespace", "path": str(path), "detail": content_ratio})
    except Exception as exc:
        issues.append({"type": "unreadable_png", "path": str(path), "detail": str(exc)})
    return metrics, issues


def audit_vector(path: Path) -> list[dict]:
    try:
        if path.suffix.lower() == ".pdf":
            if not path.read_bytes().startswith(b"%PDF-"):
                raise ValueError("missing PDF signature")
        else:
            root = ET.parse(path).getroot()
            if not root.tag.lower().endswith("svg"):
                raise ValueError(f"unexpected root element: {root.tag}")
    except Exception as exc:
        return [{"type": f"unreadable_{path.suffix.lower()[1:]}", "path": str(path), "detail": str(exc)}]
    return []


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("root", nargs="?", default=str(Path(__file__).resolve().parent))
    args = parser.parse_args()
    root = Path(args.root).resolve()
    images = [p for p in root.rglob("*") if p.is_file() and p.suffix.lower() in {".png", ".pdf", ".svg"}]
    logs = [p for p in root.rglob("*") if p.is_file() and p.suffix.lower() == ".log"]
    issues = []
    png_metrics = []
    for p in images:
        if p.stat().st_size == 0:
            issues.append({"type": "empty_image", "path": str(p)})
        if re.search(r"[<>:\"/\\|?*]|U\\+F0", p.name):
            issues.append({"type": "unsafe_filename", "path": str(p)})
        if p.suffix.lower() == ".png" and p.stat().st_size:
            metrics, png_issues = audit_png(p)
            png_metrics.append(metrics)
            issues.extend(png_issues)
        elif p.suffix.lower() in {".pdf", ".svg"} and p.stat().st_size:
            issues.extend(audit_vector(p))
    for p in logs:
        text = p.read_text(encoding="utf-8", errors="replace")
    for pattern in ("Traceback", "findfont", "empty distance matrix", "\uFFFD"):
            if pattern in text:
                issues.append({"type": "log_pattern", "pattern": pattern, "path": str(p)})
    report = {
        "root": str(root),
        "images": len(images),
        "nonempty_images": sum(p.stat().st_size > 0 for p in images),
        "logs": len(logs),
        "png_metrics": png_metrics,
        "issues": issues,
    }
    (root / "E2E_VISUAL_AUDIT.json").write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8")
    md = ["# E2E Visual Audit", "", f"- images: {report['images']}", f"- nonempty: {report['nonempty_images']}", f"- issues: {len(issues)}", ""]
    md += ["| Type | Path | Detail |", "|---|---|---|"]
    md += [f"| {x['type']} | {x['path']} | {x.get('pattern', x.get('detail', ''))} |" for x in issues]
    (root / "E2E_VISUAL_AUDIT.md").write_text("\n".join(md) + "\n", encoding="utf-8")
    print(json.dumps(report, ensure_ascii=True, indent=2))
    return 1 if issues else 0


if __name__ == "__main__":
    raise SystemExit(main())
