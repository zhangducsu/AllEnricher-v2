#!/usr/bin/env python3
"""Audit and load only the R packages directly used by retained scripts."""

from __future__ import annotations

import argparse
import json
import os
import re
import shutil
import subprocess
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[2]
R_ROOTS = [
    PROJECT_ROOT / "allenricher" / "core" / "r_scripts",
    PROJECT_ROOT / "allenricher" / "visualization" / "r_scripts",
]
BASE_PACKAGES = {"base", "grDevices", "grid", "stats", "tools", "utils"}
INSTALL_ONLY_PACKAGES = {"remotes"}


def discover_direct_dependencies() -> set[str]:
    dependencies = set()
    patterns = [
        re.compile(r"\blibrary\(\s*[\"']?([A-Za-z][A-Za-z0-9.]*)"),
        re.compile(r"\brequireNamespace\(\s*[\"']([A-Za-z][A-Za-z0-9.]*)"),
        re.compile(r"\b([A-Za-z][A-Za-z0-9.]*)::"),
    ]
    for root in R_ROOTS:
        for path in root.glob("*.R"):
            text = path.read_text(encoding="utf-8")
            for pattern in patterns:
                dependencies.update(pattern.findall(text))
    return dependencies - BASE_PACKAGES - INSTALL_ONLY_PACKAGES


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    rscript = shutil.which("Rscript")
    if not rscript:
        raise RuntimeError("Rscript not found")
    packages = sorted(discover_direct_dependencies())
    expression = (
        "pkgs <- strsplit(tail(commandArgs(TRUE), 1), ',', fixed=TRUE)[[1]]; "
        "ok <- vapply(pkgs, requireNamespace, logical(1), quietly=TRUE); "
        "cat(paste(pkgs, ok, sep='='), sep='\\n'); "
        "if (any(!ok)) quit(status=1)"
    )
    env = {
        key: value
        for key, value in os.environ.items()
        if key not in {"LANG", "LC_ALL", "LC_CTYPE"}
    }
    result = subprocess.run(
        [rscript, "-e", expression, "--args", ",".join(packages)],
        capture_output=True,
        text=True,
        encoding="utf-8",
        errors="replace",
        env=env,
    )
    report = {
        "rscript": rscript,
        "direct_packages": packages,
        "base_packages_excluded": sorted(BASE_PACKAGES),
        "install_only_packages_excluded": sorted(INSTALL_ONLY_PACKAGES),
        "returncode": result.returncode,
        "stdout": result.stdout,
        "stderr": result.stderr,
    }
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(result.stdout, end="")
    if result.stderr:
        print(result.stderr, end="", file=__import__("sys").stderr)
    return result.returncode


if __name__ == "__main__":
    raise SystemExit(main())
