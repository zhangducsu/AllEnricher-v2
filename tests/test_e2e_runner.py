"""Local E2E runner regression tests."""

from __future__ import annotations

import importlib.util
from pathlib import Path


ROOT = Path(__file__).parents[1]
SPEC = importlib.util.spec_from_file_location(
    "run_all_e2e", ROOT / "test_e2e_2026" / "run_all_e2e.py"
)
assert SPEC and SPEC.loader
RUNNER = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(RUNNER)


def test_fixed_human_version_is_not_forced_onto_nonhuman_analysis(tmp_path):
    hsa = RUNNER.redirect_output(
        "python -m allenricher analyze -i genes.txt -s hsa",
        tmp_path / "hsa-output",
        "T_HSA",
        tmp_path,
    )
    mmu = RUNNER.redirect_output(
        "python -m allenricher analyze -i genes.txt -s mmu",
        tmp_path / "mmu-output",
        "T_MMU",
        tmp_path,
    )

    assert "--use-version v20260515" in hsa
    assert "--use-version" not in mmu
