# AllEnricher v2

AllEnricher is a multi-species gene set enrichment toolkit for command-line,
Python, REST API, and local Web workflows. Version 2.1.0 preserves compatibility
with the AllEnricher v1 database layout while providing deterministic analysis,
publication-oriented figures, auditable HTML reports, and method-aware AI
interpretation.

[![Python](https://img.shields.io/badge/Python-3.8%2B-blue.svg)](https://www.python.org/)
[![License](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)

## Current Documentation

The maintained implementation matrix is in [`docs/CURRENT_IMPLEMENTATION.md`](docs/CURRENT_IMPLEMENTATION.md). Older planning and audit files under `docs/` are retained for traceability and are marked as historical archive records.

## Features

- Over-representation analysis (ORA) using the hypergeometric test.
- Gene set enrichment analysis (GSEA) through Bioconductor `fgsea`.
- Single-sample GSEA (ssGSEA) and GSVA through Bioconductor `GSVA`.
- GO, KEGG, Reactome, WikiPathways, Disease Ontology, DisGeNET, and
  transcription factor gene-set databases.
- Custom database construction from user annotations.
- Multi-species registry and TaxID-based species resolution.
- R-backed publication figures in PNG, PDF, or SVG format, with Python kept as a minimal fallback when R is unavailable.
- Interactive HTML reports containing results, figures, run metadata, and a
  Materials and Methods writing reference.
- A local Web workbench and REST API backed by the same CLI analysis path.

## Supported Databases

| Database | Content | Species scope |
| --- | --- | --- |
| GO | Gene Ontology annotations | Registry-defined species |
| KEGG | KEGG pathways | KEGG organisms |
| Reactome | Reactome pathways | Reactome-supported species |
| WikiPathways | Community-curated pathways | WikiPathways-supported species |
| DO | Disease Ontology associations | Human |
| DisGeNET (`v20190612`) | Disease-gene associations from the frozen AllEnricher v1 free snapshot | Human |
| TRRUST v2 | Curated TF-target regulatory interactions | Human and mouse |
| ChEA3 | TF-target gene-set libraries | Human |
| AnimalTFDB | Animal TF annotations with ortholog-mapped target sets | AnimalTFDB species |
| hTFtarget | Tissue-specific TF-target associations | Human |
| CUSTOM | User-built gene sets | User-defined |

Database availability depends on the files installed in the selected database
directory. Use `allenricher list-species`, `allenricher query-species`, or the
Web workbench to inspect local support.

## Analysis Methods

| Method | Required input | Default gene-set size |
| --- | --- | --- |
| ORA (`hypergeometric`) | Query gene list | 3 to unbounded |
| GSEA (`gsea`) | Ranked genes with signed numeric weights | 15 to 500 |
| TF GSEA | Ranked genes with signed numeric weights | 15 to 5,000 |
| ssGSEA (`ssgsea`) | Gene-by-sample expression matrix | 1 to unbounded |
| GSVA (`gsva`) | Gene-by-sample expression matrix | 1 to unbounded |

Gene-set size filters are applied after intersecting the database gene sets with
the genes available to the selected method. Users can override the defaults
through the existing CLI or configuration options.

## Installation

```bash
git clone https://github.com/zhangducsu/AllEnricher-v2.git
cd AllEnricher-v2
python -m pip install -e ".[visualization,api]"
```

Optional extras are available for AI clients and development tools:

```bash
python -m pip install -e ".[ai]"
python -m pip install -e ".[dev]"
```

### R Analysis Dependencies

GSEA, ssGSEA, and GSVA require a working `Rscript` installation with the
Bioconductor packages `fgsea` and `GSVA`.

```r
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
BiocManager::install(c("fgsea", "GSVA"), ask = FALSE, update = FALSE)
```

GSEA, ssGSEA, and GSVA publication figures are rendered by R by default. The Web workbench uses this default automatically and does not expose an R/Python plotting switch. Python plotting code is retained only for workflow helpers, visual-audit utilities, report fallback figures, and explicit CLI `--python-plots` runs when R is unavailable. ORA uses the maintained Python bar/lollipop renderer. R figures use `ggplot2`, `dplyr`, `tidyr`, and `scales`. The R pathway-network figure additionally requires `aPEAR` and its dependencies. Run the repository
dependency audit to check the packages required by the current R scripts:

```bash
python test_e2e_2026/18_real_world_sci/verify_r_dependencies.py
```

## Input Formats

### ORA Gene List

A plain-text file with one gene ID per line and no header:

```text
TP53
BRCA1
EGFR
```

### GSEA Ranked Gene Table

A TSV or CSV file containing `gene` and numeric `weight` columns. Weights should
be directional and should normally include both positive and negative values.

```text
gene	weight
STAT1	4.82
IRF7	3.91
MYC	-2.74
```

The Web workbench can also convert an existing differential-results table into
this two-column format after the user selects the gene and ranking-statistic
columns. AllEnricher does not perform upstream differential-expression analysis.

### Expression Matrix

A TSV or CSV matrix with gene IDs in the first column, sample names in the
header, and numeric values in all remaining cells:

```text
gene	Control_1	Control_2	Treatment_1	Treatment_2
TP53	8.2	7.9	10.4	10.1
BRCA1	5.1	5.3	4.8	4.9
```

Sample groups can be supplied as
`Control:Control_1,Control_2;Treatment:Treatment_1,Treatment_2`. The Web
workbench reads sample names from the matrix and provides an interactive group
editor.

### GMT Gene Sets

Each tab-delimited row contains a gene-set ID, a description, and one or more
genes:

```text
SET_001	Cell-cycle genes	CDK1	CCNB1	CDC20
```

## Command-Line Examples

### ORA

```bash
allenricher analyze \
  --input genes.txt \
  --species hsa \
  --databases GO,KEGG \
  --method hypergeometric \
  --output results/ora
```

Use a custom background only when it represents the genes that could have been
selected in the upstream experiment:

```bash
allenricher analyze \
  --input genes.txt \
  --species hsa \
  --databases GO \
  --method hypergeometric \
  --background-mode custom \
  --background measured_genes.txt \
  --output results/ora_custom_background
```

### GSEA

```bash
allenricher analyze \
  --ranked-genes ranked_genes.tsv \
  --species hsa \
  --databases GO,KEGG \
  --method gsea \
  --plot-types enrichment,enrichment2,barplot,lollipop,ridgeplot,emapplot \
  --output results/gsea
```

### ssGSEA and GSVA

```bash
allenricher analyze \
  --expression-matrix expression.tsv \
  --groups "Control:Control_1,Control_2;Treatment:Treatment_1,Treatment_2" \
  --species hsa \
  --databases GO \
  --method ssgsea \
  --plot-types heatmap,group_comparison,correlation \
  --output results/ssgsea
```

Replace `ssgsea` with `gsva` to run GSVA. The activity matrix retains the
official Bioconductor-compatible matrix structure.

### Custom Database

```bash
allenricher build \
  --database custom \
  --species custom_species \
  --custom-annot annotations.tsv \
  --custom-db-name CUSTOM \
  --annot-format auto \
  --database-dir ./database
```

Run `allenricher build --help` for annotation-column and hierarchy options.

## Local Web Workbench

Install the API extra and start the local service:

```bash
python -m pip install -e ".[api,visualization]"
allenricher serve --host 127.0.0.1 --port 8000
```

Open [http://127.0.0.1:8000/](http://127.0.0.1:8000/) in a browser. OpenAPI
documentation is available at [http://127.0.0.1:8000/docs](http://127.0.0.1:8000/docs).

The Web interface exposes method-specific inputs, checks installed database
support for the selected species, and stores each task under the configured API
job directory. Runtime-internal options remain controlled by the backend.

## Output

A typical analysis directory contains:

```text
results/
|-- GO_enrichment.tsv
|-- KEGG_enrichment.tsv
|-- analysis_metadata.json
|-- report.html
`-- plots/
    |-- GO_barplot.png
    `-- GO_lollipop.png
```

Result tables include stable term identifiers and readable term names. Where a
database supplies hierarchy metadata, the hierarchy is retained in a separate
column. GSEA output follows the `fgsea` result contract; ssGSEA and GSVA output
is an activity matrix with gene sets as rows and samples as columns.

The HTML report includes all generated result tables and figures, recorded run
metadata, and an English Materials and Methods writing reference based only on
the values stored for that run. Missing versions or references are shown as
`To be added` rather than inferred.

## Configuration

Generate a documented YAML configuration file:

```bash
allenricher config --output config.yaml
```

Command-line arguments take precedence over values loaded from a configuration
file. See [`config.example.yaml`](config.example.yaml) for the maintained
example.

## Testing

Run the unit and integration suite:

```bash
python -m pytest -q
```

Run the deterministic local E2E suite:

```bash
python test_e2e_2026/run_all_e2e.py --mode local --keep-going
```

Download cases remain excluded from the default E2E gate unless explicitly
enabled. Generated E2E evidence is stored under `test_e2e_2026/99_runs/` and is
not part of the maintained source-text internationalization gate.

## Citation

Zhang D, Hu Q, Liu X, Zou K, Sarkodie EK, Liu X, et al. AllEnricher: a
comprehensive gene set function enrichment tool for both model and non-model
species. *BMC Bioinformatics*. 2020;21:106.
[https://doi.org/10.1186/s12859-020-3408-y](https://doi.org/10.1186/s12859-020-3408-y)

## License

AllEnricher is distributed under the MIT License. See [`LICENSE`](LICENSE).
