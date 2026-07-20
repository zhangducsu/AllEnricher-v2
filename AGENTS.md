# AllEnricher v2 Agent Guide

This file defines repository-specific expectations for contributors and coding
agents. Read the relevant implementation and tests before changing behavior.

## Documentation Authority

Use `README.md`, `docs/CURRENT_IMPLEMENTATION.md`, and `docs/release/FINAL_RELEASE_CHECKLIST.md` as the current implementation contract. Files under `docs/plans`, `docs/reports`, `docs/specs`, and `docs/superpowers` are historical archive records unless they explicitly say otherwise.

## Project Scope

AllEnricher v2 is a Python implementation of multi-species gene set enrichment
analysis. It provides a CLI, a local Web workbench, a REST API, database
construction utilities, Python/R visualization, and HTML reports.

- Repository: <https://github.com/zhangducsu/AllEnricher-v2>
- Package version source: `pyproject.toml`
- Python requirement: 3.8 or newer
- License: MIT
- v1 compatibility: v2 can read the established `.tab.gz` and `.2desc.gz`
  species-database layout, but v1 numerical behavior is not an oracle when it
  conflicts with the documented statistical model.

This repository performs dry-data computation only. It does not define or run
wet-laboratory workflows.

## Language Policy

- All committed source code, comments, docstrings, documentation, examples,
  configuration descriptions, test descriptions, CLI output, API messages,
  Web text, and generated report templates must be written in English.
- User-facing scientific terms must follow common bioinformatics usage. Use
  `gene set`, `enrichment analysis`, `pathway`, `heatmap`, and `transcription
  factor`; do not invent literal alternatives.
- Coding agents may communicate with the repository owner in Chinese, but no
  Chinese text may be committed to the repository.
- Run `python -m pytest tests/test_english_only.py -q` before completing a
  content or interface change.

## Architecture

```text
User input
  -> CLI (`allenricher/cli.py`)
  -> configuration (`allenricher/core/config.py`)
  -> database loading (`allenricher/database/manager.py`)
  -> analysis (`allenricher/core/enrichment.py`, Bioconductor bridge)
  -> figures (`allenricher/visualization/`)
  -> HTML report (`allenricher/report/generator.py`)
```

The Web and REST API in `allenricher/api/server.py` prepare validated inputs and
invoke the same CLI path. Do not implement a second statistical workflow in the
API layer.

## Analysis Contracts

### ORA

- Method: one-sided hypergeometric test.
- Default gene-set size: 3 to unbounded.
- The tested universe must be explicit and internally consistent.
- `annotated` uses the database annotation universe.
- `custom` uses the supplied background after intersection with the database
  annotation universe.
- The query, term, and background sets must use the same identifier space.
- Result tables must include term ID, readable term name, and hierarchy when
  the source database provides it.

### GSEA

- Runtime: Bioconductor `fgseaMultilevel` through
  `allenricher/core/bioconductor.py`.
- Input: unique gene IDs with signed numeric ranking weights.
- Default gene-set size: 15 to 500 after intersection with ranked genes.
- TF-target GSEA uses a default maximum size of 5,000.
- Output must retain the official-compatible `fgsea` columns.
- Enrichment figures must be derived from deterministic running-ES data. Never
  generate simulated curves for a real analysis.

### ssGSEA and GSVA

- Runtime: Bioconductor `GSVA`.
- Input: a numeric gene-by-sample expression matrix.
- Default gene-set size: 1 to unbounded after identifier matching.
- Output must remain a gene-set-by-sample activity matrix.
- Group-comparison figures require at least two groups with complete sample
  assignment. Correlation figures require at least two samples.

## Database Contracts

Supported public database keys are defined once in
`allenricher/core/config.py::DATABASE_CATALOG`. CLI listings, API responses,
Web controls, reports, and the species registry must use this catalog instead of
maintaining independent product lists.

- GO, KEGG, Reactome, WikiPathways: functional/pathway gene sets.
- DO and DisGeNET: disease gene sets.
- TRRUST, ChEA3, AnimalTFDB, hTFtarget: TF-target gene sets.
- CUSTOM: user-built gene sets.

DisGeNET uses the frozen AllEnricher v1 free snapshot `v20190612`. The database
is still maintained upstream, but later releases are not freely available
through the previous distribution route. Every user-facing DisGeNET label must
show the frozen snapshot version.

TRRUST supports human and mouse. ChEA3 and hTFtarget are human-only. AnimalTFDB
species support is prepared by its downloader and incorporated into the shared
species registry. Small fixed species lists may be registered explicitly, but
the query path must remain the same as for dynamically generated registry data.

Use TaxID as the unique species identity. KEGG codes and names are lookup aliases,
not primary identifiers.

## Visualization Contracts

Current default figure types:

| Method | Figure types |
| --- | --- |
| ORA | `barplot`, `lollipop` |
| GSEA | `enrichment`, `enrichment2`, `barplot`, `lollipop`, `ridgeplot`, `emapplot` |
| ssGSEA/GSVA | `heatmap`, `group_comparison`, `correlation` |

Visualization requirements:

- Use readable term names on figures; retain stable IDs in tables and file
  metadata.
- Apply categorical, sequential, and diverging palettes only to compatible data
  roles.
- Keep lollipop background bands independent of the data palette.
- Continuous legends must remain distinguishable at both ends; do not use an
  effectively invisible endpoint.
- Heatmap dimensions must adapt to the displayed rows and columns. The default
  display is limited to the most variable or most relevant rows while the full
  matrix remains in the result file.
- Legends on heatmaps belong on the right, grouped by annotation role, compact
  within a group, and separated between groups.
- Figures must support the configured PNG, PDF, or SVG output without changing
  statistical content.
- R and Python implementations must honor the same style and semantic palette
  configuration where both backends exist.

## Web and Report Contracts

- The Web interface is a progressive workflow: method, method-specific input,
  species/database scope, and optional advanced settings.
- Unsupported databases remain visible but disabled for the selected species.
- The Web interface must not expose backend runtime choices that ordinary users
  cannot evaluate safely.
- Input controls must describe the accepted file structure in concise English.
- The Web application and HTML report must display the package version, the
  GitHub repository URL, and the verified AllEnricher citation.
- Reports must include all generated figures, readable term names, accurate tab
  labels, output files, and an English Materials and Methods writing reference.
- Materials and Methods prose may use only metadata recorded by the current run.
  Missing versions or references must be shown as `To be added`.

## Editing Rules

1. Read the implementation, call path, and relevant tests before editing.
2. Choose the smallest implementation that satisfies the business contract.
3. Keep changes within the responsible module; avoid unrelated refactors.
4. Preserve user changes in a dirty worktree.
5. Do not hide failures, skipped tests, placeholders, or unverified behavior.
6. Use structured parsers for structured data.
7. Keep deterministic logic in code, not in an AI model.
8. Add comments only where the behavior is not self-explanatory.
9. Do not change statistical output contracts merely to make a test pass.
10. Never replace missing biological data with random or fabricated values.

## Verification

Choose tests based on the affected surface. The minimum repository checks are:

```bash
python -m compileall -q allenricher tests test_e2e_2026
python -m pytest tests/test_english_only.py -q
python -m pytest -q
git diff --check
```

For Web changes, also extract the inline JavaScript and run `node --check`, then
start a local server and inspect the desktop and narrow layouts in a browser.

For R changes, parse every maintained `.R` file with `Rscript` and run the
relevant smoke tests in PNG, PDF, and SVG formats. The current R package audit is:

```bash
python test_e2e_2026/18_real_world_sci/verify_r_dependencies.py
```

For full deterministic E2E validation:

```bash
python test_e2e_2026/run_all_e2e.py --mode local --keep-going
```

Download tests remain disabled unless the repository owner explicitly restores
them. Do not overwrite historical E2E evidence or production database files.

## Completion Criteria

A change is complete only when:

- the requested behavior is implemented end to end;
- business outputs, not only exit codes, are validated;
- all failures and skipped checks are reported;
- committed user-facing text is English and scientifically accurate;
- no unrelated worktree changes were reverted; and
- generated evidence is stored in the documented E2E location when requested.
