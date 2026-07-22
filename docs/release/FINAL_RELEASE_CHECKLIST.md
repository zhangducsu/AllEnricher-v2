# AllEnricher-v2 Final Release Checklist

This checklist is the release gate for AllEnricher v2.1.1 maintenance and patch releases. Use it with `docs/CURRENT_IMPLEMENTATION.md`, which is the current implementation matrix.

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
- The release container builds from the pinned base digest and passes the CLI and R-version smoke checks.

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
- Current default figures match `docs/CURRENT_IMPLEMENTATION.md`: ORA `barplot,lollipop`; GSEA `enrichment,enrichment2,barplot,lollipop,ridgeplot` plus R-only `emapplot`; ssGSEA/GSVA `heatmap,group_comparison,correlation`.
- Removed or non-current figure names are not advertised as active defaults: `bubble`, `dotplot`, `network`, `upset`, `volcano`, `method_comparison`, `cnetplot`, and `circos`.

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

Build and verify the release container before tagging:

```powershell
docker build --pull=false -t allenricher:2.1.1 .
docker run --rm allenricher:2.1.1 --version
docker run --rm --entrypoint Rscript allenricher:2.1.1 -e "stopifnot(as.character(packageVersion('fgsea')) == '1.38.0', as.character(packageVersion('GSVA')) == '2.6.2')"
```

Record the published image digest with the release and use the digest, not a
mutable tag, in archived workflow commands.

## Release Decision

- `GO`: required gate passes, manual review finds no release blocker.
- `GO_WITH_NOTES`: required gate passes, only documented limitations remain.
- `NO_GO`: any required gate fails, generated reports are misleading, or analysis correctness is uncertain.
