# Real-World Scientific E2E Suite

This suite validates AllEnricher against four public EMBL-EBI Expression Atlas studies. The
registered matrix contains 108 independent `dataset x database x method` cases across ORA, GSEA,
ssGSEA, and GSVA. Accessions, contrasts, sample filters, thresholds, ranked statistics, and
expression matrices are fixed before any enrichment result is inspected.

The studies cover human, *Drosophila melanogaster*, cattle, and *Saccharomyces cerevisiae*. The
cattle activity matrix aggregates technical runs from each animal and compares distinct groups at
the same time point. Dataset and database-specific decisions are documented in
`PROTOCOL_AMENDMENTS.md`.

## Registered Inputs

- **ORA:** the differential-gene list distributed by Expression Atlas for the registered contrast.
- **GSEA:** the complete Expression Atlas analysis table ranked by its recorded statistic.
- **ssGSEA and GSVA:** the public count matrix with groups derived from the study design.
- **GSVA kernel:** `kcdf=Poisson` for count data.
- **Databases:** the isolated snapshot under `00_input_data/real_world_sci/database_snapshot/`.
- **DisGeNET:** human only, using the archived AllEnricher-v1 free snapshot. DisGeNET remains an
  active resource, but newer releases are not freely downloadable under the same terms.

Preparing public fixtures is separate from testing the `allenricher download` command. Download
CLI cases remain outside this suite.

## Execution

```powershell
# 1. Verify the R packages used by the current implementation.
python test_e2e_2026/18_real_world_sci/verify_r_dependencies.py

# 2. Download and convert Expression Atlas inputs, or validate the existing cache.
python test_e2e_2026/18_real_world_sci/prepare_expression_atlas.py

# 3. Build the isolated E2E database snapshot.
python test_e2e_2026/18_real_world_sci/prepare_database_snapshot.py

# 4. Reproduce one case while debugging.
python test_e2e_2026/18_real_world_sci/run_real_world_sci.py `
  --case human_airway_dex__GO__hypergeometric

# 5. Run the registered 108-case primary matrix.
python test_e2e_2026/18_real_world_sci/run_real_world_sci.py --mode primary

# 6. Validate cached inputs and reproduce the matrix without network access.
python test_e2e_2026/18_real_world_sci/prepare_expression_atlas.py --offline
python test_e2e_2026/18_real_world_sci/prepare_database_snapshot.py --offline
python test_e2e_2026/18_real_world_sci/run_real_world_sci.py `
  --mode offline --compare-to <PRIMARY_RUN_DIR>
```

## Acceptance Criteria

- `MATRIX_COVERAGE.tsv` contains exactly 108 data rows and every case is `PASS`.
- ORA P values are independently reproduced with the SciPy hypergeometric distribution and
  statsmodels Benjamini-Hochberg adjustment.
- GSEA is independently reproduced with Bioconductor fgsea using identical ranks, gene sets, and
  random seed.
- ssGSEA and GSVA are independently reproduced with Bioconductor GSVA to a maximum absolute error
  of `1e-8`.
- Official result tables contain analysis values only; runtime provenance is stored in
  `analysis_metadata.json`.
- Every case archives its command, logs, exit code, manifests, result tables, HTML report,
  applicable plots, and numerical oracle output.
- Non-allowlisted errors, warnings, encoding defects, font warnings, or invalid clustering
  messages fail the case.
- `E2E_VISUAL_AUDIT.json` reports zero issues. Contact sheets remain mandatory for manual review.

## Failure Handling

Failed runs are immutable evidence. Reproduce a case with `--case` in a new run directory; do not
overwrite, delete, or manually edit the failed output. A real failure must not be relabeled as
`SKIP` or `EXPECTED_FAIL`. After a fix, rerun the complete primary and offline matrices in addition
to the focused case.
