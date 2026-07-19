"""Checks that shipped application text is English-only."""

from __future__ import annotations

import re
import subprocess
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
TEXT_SUFFIXES = {
    ".cfg",
    ".html",
    ".ini",
    ".json",
    ".md",
    ".py",
    ".pyi",
    ".r",
    ".rst",
    ".sh",
    ".toml",
    ".txt",
    ".yaml",
    ".yml",
}
EXCLUDED_PARTS = {
    ".git",
    ".pytest_cache",
    "__pycache__",
    "build",
    "deepseek_test_output",
    "dist",
    "test_data",
    "test_output",
}
EXCLUDED_PREFIXES = (
    "docs/",
    "examples/",
    "tests/",
    "test_e2e_2026/",
    "test_e2e_2026/00_input_data/",
    "test_e2e_2026/99_runs/",
)
INCLUDED_PREFIXES = ("allenricher/",)
HAN_RE = re.compile(r"[\u3400-\u4dbf\u4e00-\u9fff\uf900-\ufaff]")
PLACEHOLDER_RE = re.compile(r"(?:PLACEHOLDER|SPANPLACEHOLDER)\d+")
MACHINE_TRANSLATION_ARTIFACTS = (
    "allencer",
    "animalltfdb",
    "build blood",
    "co-source map",
    "dise ontoroy",
    "diseese",
    "disese",
    "end-to-endpoint",
    "escort analysis",
    "gene collection",
    "gene-access connection",
    "gene set collection analysis analysis",
    "genetic collection",
    "genetic condominium",
    "genetic =",
    "genetic fusion",
    "genetically covered",
    "genome collection",
    "homogeneity map",
    "litin_name",
    "no, i'm not",
    "other organiser",
    "passes x",
    "please.2",
    "please transmit",
    "registration form",
    "retransmission factor",
    "richness analysis",
    "roadway",
    "sorting gene",
    "sorting number",
    "specialies",
    "supergeometry",
    "the blog is a blog",
    "the road",
    "thermal map",
    "thermogram genetic",
    "time-consuming:",
    "transposition factor",
    "texid",
    "ttrutt",
    "we found it together",
    "visible access",
    "♪",
    "\u9225",
    "\u922b",
    "\u923e",
    "\u923a",
    "\u9241",
    "\u9983",
    "\u8133",
)


def _repository_text_files() -> list[Path]:
    result = subprocess.run(
        ["git", "ls-files", "--cached", "--others", "--exclude-standard"],
        cwd=ROOT,
        check=True,
        capture_output=True,
        text=True,
        encoding="utf-8",
    )
    files: list[Path] = []
    for relative_name in result.stdout.splitlines():
        normalized = relative_name.replace("\\", "/")
        if not normalized.startswith(INCLUDED_PREFIXES):
            continue
        path = ROOT / relative_name
        if not path.is_file() or path.suffix.lower() not in TEXT_SUFFIXES:
            continue
        if any(part.lower().endswith(".egg-info") for part in path.parts):
            continue
        if any(part in EXCLUDED_PARTS for part in path.parts):
            continue
        if normalized.startswith(EXCLUDED_PREFIXES):
            continue
        files.append(path)
    return files


def test_repository_maintained_text_is_english_only() -> None:
    violations: list[str] = []
    for path in _repository_text_files():
        text = path.read_text(encoding="utf-8")
        for line_number, line in enumerate(text.splitlines(), start=1):
            if HAN_RE.search(line):
                violations.append(f"{path.relative_to(ROOT)}:{line_number}: contains Han text")
            if "\ufffd" in line:
                violations.append(f"{path.relative_to(ROOT)}:{line_number}: contains U+FFFD")
            if PLACEHOLDER_RE.search(line):
                violations.append(f"{path.relative_to(ROOT)}:{line_number}: contains translation placeholder")
    assert not violations, "\n" + "\n".join(violations[:100])


def test_web_and_report_defaults_are_english() -> None:
    web = (ROOT / "allenricher/api/static/index.html").read_text(encoding="utf-8")
    report = (ROOT / "allenricher/report/generator.py").read_text(encoding="utf-8")

    assert '<html lang="en">' in web
    assert "zh-CN" not in web
    assert '<html lang="zh-CN">' not in report


def test_repository_avoids_known_machine_translation_artifacts() -> None:
    violations: list[str] = []
    for path in _repository_text_files():
        if path == Path(__file__).resolve():
            continue
        for line_number, line in enumerate(
            path.read_text(encoding="utf-8").splitlines(), start=1
        ):
            normalized = line.casefold().replace("_", " ")
            for artifact in MACHINE_TRANSLATION_ARTIFACTS:
                if artifact in normalized:
                    violations.append(
                        f"{path.relative_to(ROOT)}:{line_number}: contains '{artifact}'"
                    )
    assert not violations, "\n" + "\n".join(violations[:100])
