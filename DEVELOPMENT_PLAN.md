# AllEnricher v2 Development Status

> Current implementation source of truth: use `docs/CURRENT_IMPLEMENTATION.md`, `README.md`, `AGENTS.md`, and `docs/release/FINAL_RELEASE_CHECKLIST.md` for v2.1.0 behavior. This file is a planning aid, not the release contract.


This document summarizes the implemented product surface and the remaining
release work. Detailed historical plans and evidence remain under `docs/` and
`test_e2e_2026/`.

## Implemented Product Surface


### Analysis


- ORA with a one-sided hypergeometric test and configurable background.
- GSEA through Bioconductor `fgseaMultilevel`.
- ssGSEA and GSVA through Bioconductor `GSVA`.
- Method-specific gene-set size defaults with CLI/configuration overrides.
- Official-compatible output contracts for GSEA and activity matrices.

### Databases


- GO, KEGG, Reactome, WikiPathways, Disease Ontology, and the frozen DisGeNET
  `v20190612` snapshot.
- TRRUST, ChEA3, AnimalTFDB, and hTFtarget TF-target gene sets.
- Custom database construction from flat or hierarchical annotations.
- Shared TaxID-based species registry used by CLI, API, and Web queries.
- Version metadata and local database discovery.

### Visualization and Reports


- ORA bar and lollipop plots.
- GSEA single-pathway, multi-pathway, bar, lollipop, ridge, and pathway-network
  plots.
- ssGSEA/GSVA activity heatmaps, group comparisons, and sample-correlation
  plots.
- Python and R plotting backends with semantic categorical, sequential, and
  diverging palette roles.
- PNG, PDF, and SVG output.
- HTML reports containing result tables, figures, metadata, and an English
  Materials and Methods writing reference.

### Interfaces


- CLI subcommands for analysis, database preparation and inspection, species
  queries, service startup, configuration generation, and TF enrichment.
- FastAPI REST endpoints and a local Web workbench using the same CLI workflow.
- Persistent API job metadata and downloadable analysis artifacts.

### Validation Infrastructure


- Unit and integration tests under `tests/`.
- Deterministic command/API/R-plot E2E infrastructure under `test_e2e_2026/`.
- Real-world analysis matrix and visual-review contact sheets.
- Repository-wide English-only source-text gate.

## Release Priorities


1. Complete the English semantic rewrite of maintained source comments,
   documentation, templates, messages, and tests.
2. Pass Python compilation, JavaScript syntax validation, R parsing, targeted
   interface tests, and the complete pytest suite.
3. Run the deterministic local E2E suite without restoring disabled download
   cases.
4. Review generated result tables, HTML reports, and representative figures.
5. Build and inspect the release wheel to confirm that static assets, R scripts,
   database registry data, and templates are included.
6. Record any external-service limitations separately from deterministic release
   readiness.

## Release Acceptance Criteria


- No unexpected test failures.
- No Chinese text or translation placeholders in maintained Git content.
- No syntax errors in Python, inline Web JavaScript, or maintained R scripts.
- CLI, REST API, Web, and report behavior remain aligned.
- Result tables retain stable IDs, readable names, hierarchy metadata where
  available, and method-specific official output contracts.
- Figures are non-empty, readable, and use compatible semantic palette roles.
- Download tests remain explicitly skipped until the repository owner restores
  them.
- All skipped, expected-failure, and external-service cases are visible in the
  final test report.
