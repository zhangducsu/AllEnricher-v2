"""Shared helpers for plot labels, filenames, and compact sizing."""

from __future__ import annotations

import re
import unicodedata
from typing import Tuple


_SMALL_WORDS = {"a", "an", "and", "as", "at", "by", "for", "from", "in", "of", "on", "or", "the", "to", "via", "with"}
_GO_CATEGORY_PREFIXES = {
    "biological_process",
    "biological process",
    "cellular_component",
    "cellular component",
    "molecular_function",
    "molecular function",
}


def _format_label_word(word: str, first: bool) -> str:
    if not word:
        return word
    if len(word) > 1 and word.upper() == word and any(ch.isalpha() for ch in word):
        return word
    lowered = word.lower()
    if not first:
        return lowered
    return lowered[:1].upper() + lowered[1:]


def normalize_pathway_case(label: str) -> str:
    """Return readable sentence-style pathway labels without touching acronyms."""
    text = re.sub(r"\s+", " ", str(label or "")).strip()
    if not text:
        return text
    words = text.split(" ")
    starts_upper = sum(1 for word in words if word[:1].isupper())
    looks_title_case = len(words) > 1 and starts_upper >= max(2, len(words) // 2)
    if text[:1].islower() or looks_title_case:
        return " ".join(_format_label_word(word, idx == 0) for idx, word in enumerate(words))
    return text


def clean_pathway_label(value: object, normalize_case: bool = True) -> str:
    """Clean database/category prefixes from display labels while preserving IDs elsewhere."""
    text = re.sub(r"\s+", " ", str(value or "")).strip()
    if not text:
        return ""

    if "|" in text:
        parts = [part.strip() for part in text.split("|") if part.strip()]
        if parts:
            text = parts[-1]
    elif ":" in text:
        prefix, suffix = text.split(":", 1)
        if prefix.strip().lower().replace(" ", "_") in _GO_CATEGORY_PREFIXES:
            text = suffix.strip()

    text = text.replace("_", " ")
    text = re.sub(r"\s+", " ", text).strip()
    return normalize_pathway_case(text) if normalize_case else text


def safe_plot_stem(value: object, fallback: str = "term") -> str:
    """Create an ASCII-only filesystem-safe stem for plot IDs."""
    text = re.sub(r"[\ue000-\uf8ff]+", "_", str(value or ""))
    text = unicodedata.normalize("NFKD", text)
    text = text.encode("ascii", "ignore").decode("ascii")
    text = re.sub(r"[^A-Za-z0-9._-]+", "_", text)
    text = re.sub(r"_+", "_", text).strip("._-")
    return text or fallback


def term_figure_size(
    n_terms: int,
    width: float = 8.0,
    min_height: float = 2.8,
    row_height: float = 0.42,
    max_height: float = 12.0,
) -> Tuple[float, float]:
    """Return compact term-list figure dimensions."""
    n = max(int(n_terms), 1)
    height = min(max_height, max(min_height, 1.6 + n * row_height))
    return width, height
