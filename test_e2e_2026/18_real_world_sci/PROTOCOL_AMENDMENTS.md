# Real-World E2E Protocol Decisions

## Registered Data Matrix

The primary matrix uses four public, non-infectious Expression Atlas studies. Each accession,
contrast, sample filter, database list, and analysis method is fixed before results are inspected.
The 108 cases are defined in `case_matrix.yaml` and cover ORA, GSEA, ssGSEA, and GSVA.

- **Human:** E-GEOD-52778, dexamethasone versus untreated airway smooth-muscle cells.
- **Model organism:** E-MTAB-5069, *Drosophila melanogaster* `fas P218` versus wild type.
- **Non-model organism:** E-MTAB-5838, high-risk versus low-risk cattle at one week.
- **Microorganism:** E-MTAB-9355, *Saccharomyces cerevisiae* `Tsa1` deletion versus wild type
  under the registered protein-folding stress condition.

For E-MTAB-5838, technical runs from the same animal are aggregated before sample-level
activity analysis. They are not treated as independent biological replicates. The registered
contrast compares distinct high-risk and low-risk animals at the same time point; the earlier
longitudinal comparison is not part of the matrix.

Two previously considered studies are excluded from the primary matrix. E-MTAB-5879 does not
provide enough native WikiPathways gene sets to satisfy the preregistered per-database coverage
criterion, and E-GEOD-73681 yields too few query genes at the study's fixed published filter.
Neither study is substituted or relaxed after results are observed.

## Transcription-Factor Database Boundary

AnimalTFDB 4.0 supplies species-specific TF and transcription cofactor classifications, but it
does not publish complete TF-target regulatory networks for every supported animal. Consequently:

- Human analyses use TRRUST, ChEA3, and hTFtarget as direct TF-target sources.
- Drosophila and cattle AnimalTFDB cases use the explicitly recorded ortholog-mapping workflow
  defined by the E2E database snapshot.
- AnimalTFDB-derived collections remain distinct from hTFtarget and are never presented as a
  native AnimalTFDB TF-target network.

All source provenance and mapping steps are retained with the isolated database snapshot. Inputs
must not be replaced, thresholds must not be relaxed, and database coverage must not be changed in
response to enrichment significance.
