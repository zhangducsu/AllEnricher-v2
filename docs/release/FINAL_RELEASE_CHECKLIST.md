# AllEnricher-v2 Final Release Checklist

This checklist is the release gate for the current v2 closing phase. From this point, avoid new features unless a check below exposes a blocking defect.

## Release Rule

- Scope is frozen: fix blockers only.
- Download-heavy and live-network tests are optional unless the release explicitly changes downloader or online-service code.
- AI tests must use mock backends unless a temporary real key is intentionally provided for a separate evidence run.
- Any failure must be recorded as `PASS`, `FAIL`, `SKIP`, or `KNOWN_LIMITATION`; do not silently ignore it.

## Required Gate

- `git status` is clean and local `main` matches `origin/main`.
- Python files compile.
- Focused pytest passes for CLI, API, AI, report, database registry, and plotting integration.
- CLI smoke commands return successfully.
- API/Web smoke can start the local service and expose the expected endpoints.
- HTML report smoke includes summary, result table, plots, AI interpretation state, and methods-writing reference when available.
- No user-facing Chinese remains in package code, web UI, generated reports, or CLI output. Tests and internal documentation are exempt.

## Manual Review

- Web workbench shows unsupported databases as disabled instead of hiding tool capability.
- Species selection uses taxid as the stable identifier.
- Database support counts include TF databases and special human-only or human/mouse databases.
- DisGeNET is labeled `v20190612` wherever it appears.
- ORA result tables include term/pathway/disease/TF IDs and readable names or descriptions.
- Hierarchical databases include hierarchy text when available.
- GSEA tables keep fgsea-compatible columns and show readable pathway names.
- ssGSEA/GSVA tables keep official activity-matrix structure.
- AI interpretation cites valid evidence IDs and avoids result discussion beyond the provided enrichment evidence.
- AI failure does not fail the analysis job; it is reported as an AI-specific warning with the error artifact linked.
- Figures are readable, compact, and use the current palette/style controls consistently.

## One-command Validation

Run the required local release gate:

```powershell
powershell -ExecutionPolicy Bypass -File docs/release/run_final_release_checks.ps1
```

Run the extended gate before tagging a public release:

```powershell
powershell -ExecutionPolicy Bypass -File docs/release/run_final_release_checks.ps1 -Full
```

Outputs are saved under:

```text
test_e2e_2026/99_runs/final_release_<timestamp>/
```

## Release Decision

- `GO`: required gate passes, manual review finds no release blocker.
- `GO_WITH_NOTES`: required gate passes, only documented limitations remain.
- `NO_GO`: any required gate fails, generated reports are misleading, or analysis correctness is uncertain.
