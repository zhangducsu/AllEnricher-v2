import json
from pathlib import Path

from test_e2e_2026.tools.audit_e2e_outputs import audit_run


def _write_run(tmp_path: Path, stderr: str):
    case_dir = tmp_path / "01_cmd_analyze" / "CASE001"
    case_dir.mkdir(parents=True)
    (case_dir / "stdout.log").write_text("", encoding="utf-8")
    (case_dir / "stderr.log").write_text(stderr, encoding="utf-8")
    summary = {
        "status_counts": {"PASS": 1},
        "results": [{
            "id": "CASE001",
            "status": "PASS",
            "run_dir": str(case_dir),
        }],
    }
    (tmp_path / "E2E_SUMMARY.json").write_text(json.dumps(summary), encoding="utf-8")


def test_audit_flags_findfont_in_pass_case(tmp_path):
    _write_run(tmp_path, "WARNING - findfont: Font family 'Helvetica' not found\n")

    audit = audit_run(tmp_path)

    messages = [issue["message"] for issue in audit["issues"]]
    assert any("findfont" in message for message in messages)
    assert audit["issue_count"] >= 1


def test_audit_allows_declared_skip_warning(tmp_path):
    _write_run(tmp_path, "WARNING - 样本相关性图跳过：pathway 数不足 2，无法计算可靠样本相关性\n")

    audit = audit_run(tmp_path)

    assert audit["issues"] == []
