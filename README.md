# AllEnricher v2

AllEnricher is a multi-species gene set enrichment toolkit for command-line,
Python, REST API, and local Web workflows. Version 2.1.0 keeps compatibility
with the AllEnricher v1 database layout while adding deterministic GSEA,
ssGSEA, GSVA, publication-oriented figures, auditable HTML reports, species
registry queries, and method-aware AI interpretation.

[![Python](https://img.shields.io/badge/Python-3.8%2B-blue.svg)](https://www.python.org/)
[![License](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)

## Documentation Status

This README is the user-facing entry point for the current v2 implementation.
The maintained implementation matrix is in
[`docs/CURRENT_IMPLEMENTATION.md`](docs/CURRENT_IMPLEMENTATION.md). Older files
under `docs/plans`, `docs/reports`, `docs/specs`, and `docs/superpowers` are
historical planning or audit records unless they explicitly point back to the
current implementation document.

## Main Features

- ORA with a one-sided hypergeometric test.
- GSEA through Bioconductor `fgsea`.
- ssGSEA and GSVA through Bioconductor `GSVA`.
- GO, KEGG, Reactome, WikiPathways, Disease Ontology, DisGeNET, TF, and custom
  gene-set databases.
- Database download, version inspection, cleanup, species-database build, and
  species-registry query commands.
- TaxID-centered species registry with KEGG code and Latin-name lookup.
- R-backed publication figures for GSEA, ssGSEA, and GSVA, with a minimal
  Python fallback for CLI runs. ORA uses the maintained Python bar/lollipop
  renderer.
- HTML reports with result tables, figures, metadata, AI interpretation state,
  and a Materials and Methods writing reference.
- Local Web workbench and REST API backed by the same CLI analysis path.
- Optional AI interpretation with evidence IDs linking narrative claims back to
  result-table rows.

## Supported Databases

Database availability depends on the files installed in the selected database
directory. Use `allenricher list-species`, `allenricher query-species`,
`allenricher list-versions`, or the Web workbench to inspect local support.

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
| AnimalTFDB | Animal TF annotations with ortholog-mapped target sets | AnimalTFDB species with prepared local files |
| hTFtarget | Tissue-specific TF-target associations | Human |
| CUSTOM | User-built gene sets | User-defined |

DisGeNET is not downloaded from current DisGeNET releases because later data are
not freely downloadable through the same public route. AllEnricher v2 can reuse
the frozen v1 free snapshot and labels it as `v20190612` in user-facing output.

## Analysis Methods

| Method | Required input | Default gene-set size after gene intersection |
| --- | --- | --- |
| ORA (`hypergeometric`) | Query gene list | 3 to unbounded |
| GSEA (`gsea`) | Ranked genes with signed numeric weights | 15 to 500 |
| TF GSEA | Ranked genes with signed numeric weights | 15 to 5,000 |
| ssGSEA (`ssgsea`) | Gene-by-sample expression matrix | 1 to unbounded |
| GSVA (`gsva`) | Gene-by-sample expression matrix | 1 to unbounded |

Gene-set size filters are applied after intersecting database gene sets with the
genes available to the selected method. Users can override these defaults
through CLI options or a YAML/JSON configuration file.

## Installation

```bash
git clone https://github.com/zhangducsu/AllEnricher-v2.git
cd AllEnricher-v2
python -m pip install -e ".[visualization,api]"
```

Optional extras:

```bash
python -m pip install -e ".[ai]"   # AI client libraries
python -m pip install -e ".[dev]"  # pytest and development checks
```

## R Dependencies

GSEA, ssGSEA, and GSVA require a working `Rscript` installation with
Bioconductor packages `fgsea` and `GSVA`.

```r
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
BiocManager::install(c("fgsea", "GSVA"), ask = FALSE, update = FALSE)
```

R publication figures use `ggplot2`, `dplyr`, `tidyr`, and `scales`. The R
pathway-network figure additionally requires `aPEAR` and its dependencies. Check
the current R script requirements with:

```bash
python test_e2e_2026/18_real_world_sci/verify_r_dependencies.py
```

## Input Formats

### ORA Gene List

Plain text with one gene ID per line and no header:

```text
TP53
BRCA1
EGFR
```

### GSEA Ranked Gene Table

TSV or CSV with `gene` and numeric `weight` columns. Weights should be signed
and should normally include both positive and negative values.

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

TSV or CSV with gene IDs in the first column, sample names in the header, and
numeric values in all remaining cells:

```text
gene	Control_1	Control_2	Treatment_1	Treatment_2
TP53	8.2	7.9	10.4	10.1
BRCA1	5.1	5.3	4.8	4.9
```

Sample groups can be supplied as:

```text
Control:Control_1,Control_2;Treatment:Treatment_1,Treatment_2
```

The Web workbench reads sample names from the matrix and provides an interactive
group editor, so users do not need to prepare this string manually.

### GMT Gene Sets

Each tab-delimited row contains a gene-set ID, a description, and one or more
genes:

```text
SET_001	Cell-cycle genes	CDK1	CCNB1	CDC20
```

## Database Download and Version Management

Download shared source files into a database directory:

```bash
allenricher download \
  --databases GO,KEGG,Reactome,DO,WikiPathways \
  --species hsa \
  --database-dir ./database
```

Useful download options:

```bash
allenricher download -d GO,KEGG --database-dir ./database --workers 8
allenricher download -d GO --database-dir ./database --force
allenricher download -d GO --database-dir ./database --no-verify
allenricher download -d GO --database-dir ./database --no-multi-thread
```

Download TF source data:

```bash
allenricher download -d TRRUST --database-dir ./database
allenricher download -d ChEA3 --database-dir ./database
allenricher download -d AnimalTFDB --species Bos_taurus,Sus_scrofa --database-dir ./database
```

The compatibility flags remain available:

```bash
allenricher download -d TF --trrust --database-dir ./database
allenricher download -d TF --chea3 --database-dir ./database
allenricher download -d TF --animaltfdb --species Bos_taurus --database-dir ./database
```

AnimalTFDB species names use the official underscore Latin-name form, for
example `Bos_taurus` or `Sus_scrofa`. When no AnimalTFDB species is supplied,
the command refreshes TF registry coverage and downloads the human hTFtarget
library.

Inspect and maintain local versions:

```bash
allenricher check-update --database-dir ./database
allenricher check-update --database-dir ./database --json
allenricher list-versions --database-dir ./database
allenricher list-versions --database-dir ./database --json
allenricher list-versions --database-dir ./database --lineage
allenricher cleanup --database-dir ./database --dry-run --keep 2
```

`cleanup` can delete old database snapshots when `--dry-run` is omitted. Review
its preview before running a real cleanup.

## Building Species Databases

After downloading shared sources, build analysis-ready database artifacts for a
species. TaxID is the stable species identity; KEGG code is a convenient alias
when available.

```bash
allenricher build \
  --species hsa \
  --taxonomy 9606 \
  --databases GO,KEGG,Reactome,WikiPathways,DO,DisGeNET \
  --database-dir ./database \
  --gene-info ./database/basic/gene_info.gz \
  --latin-name Homo_sapiens
```

Build from user-provided annotations:

```bash
allenricher build \
  --species custom_species \
  --taxonomy 999999 \
  --databases custom \
  --custom-annot annotations.tsv \
  --custom-db-name CUSTOM \
  --annot-format auto \
  --hierarchy-sep "|" \
  --database-dir ./database
```

Custom GO or KEGG annotations can be supplied with `--go-annot` or
`--kegg-annot`. Annotation rows may include hierarchy text such as
`Metabolism|Amino acid metabolism|Arginine biosynthesis`; when hierarchy is
available, result tables retain it and ORA barplots can use it for category
coloring.

## Species Registry and Database Queries

List supported species from the unified registry:

```bash
allenricher list-species --summary
allenricher list-species --format table
allenricher list-species --format tsv
allenricher list-species --format json
```

Filter by database support:

```bash
allenricher list-species --go --kegg --summary
allenricher list-species --reactome --wikipathways
allenricher list-species --trrust
allenricher list-species --chea3
allenricher list-species --animaltfdb
allenricher list-species --htftarget
```

Query one species:

```bash
allenricher query-species --taxid 9606
allenricher query-species --name Homo_sapiens
allenricher query-species --kegg hsa
```

The Web workbench uses the same registry data. Unsupported databases remain
visible but disabled for the selected species, so users can see which functions
exist and which require an additional species build.

## Running Enrichment Analyses

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

Replace `ssgsea` with `gsva` to run GSVA. The activity matrix keeps the
Bioconductor-compatible matrix structure.

## Transcription Factor Analysis

TF databases can be used through the regular `analyze` command by passing TF
names in `--databases`:

```bash
allenricher analyze \
  --input genes.txt \
  --species hsa \
  --databases TRRUST,ChEA3,hTFtarget \
  --method hypergeometric \
  --tf-library ARCHS4,ENCODE \
  --tf-tissue liver,blood \
  --output results/tf_ora
```

For GSEA with TF gene sets:

```bash
allenricher analyze \
  --ranked-genes ranked_genes.tsv \
  --species hsa \
  --databases TRRUST,ChEA3 \
  --method gsea \
  --tf-max-size 5000 \
  --output results/tf_gsea
```

The compatibility entry point `tf-enrich` is retained for TF ORA and TF GSEA:

```bash
allenricher tf-enrich \
  --input genes.txt \
  --species hsa \
  --database trrust \
  --method ora \
  --report \
  --output results/tf_enrich
```

TF options include `--tf-library` for ChEA3, `--tf-tissue` for hTFtarget,
`--tf-regulation` for TRRUST, `--tf-min-size`, `--tf-max-size`, and optional
`--tf-combine` consensus ranking.

## Figures and Color Controls

Default plot types:

| Method | Figure types |
| --- | --- |
| ORA | `barplot`, `lollipop` |
| GSEA | `enrichment`, `enrichment2`, `barplot`, `lollipop`, `ridgeplot`; `emapplot` when requested |
| ssGSEA / GSVA | `heatmap`, `group_comparison`, `correlation` |

GSEA, ssGSEA, and GSVA publication figures use R by default. Use
`--python-plots` only when the minimal Python fallback is explicitly desired.
The Web workbench does not expose an R/Python switch and uses the backend
default.

```bash
allenricher analyze ... --plot-format png --plot-dpi 300
allenricher analyze ... --style nature
allenricher analyze ... --categorical-palette okabe_ito
allenricher analyze ... --sequential-palette blues
allenricher analyze ... --diverging-palette blue_red
```

Palette roles are separated: categorical palettes are used for groups and
categories, sequential palettes for one-direction continuous values such as
significance, and diverging palettes for centered values such as NES, activity,
and correlation.

## AI Interpretation

AI interpretation is optional. It is method-aware and uses structured evidence
prepared by code before any model call:

- ORA evidence includes term ID, term name, adjusted P value, EnrichFactor,
  gene count, and hit genes.
- GSEA evidence separates positive and negative NES terms and includes ES, NES,
  p value, adjusted P value, size, and leading-edge genes.
- ssGSEA/GSVA evidence summarizes pathway activity matrices, group means,
  group differences, and outlier samples when groups are available.
- TF evidence includes TF name, database source, significance, target-set size,
  matched targets, and ranking information when available.

Each selected evidence row receives a stable `evidence_id` such as `GO:R001` or
`GSEA_Reactome:R002`. HTML reports show these IDs and link them back to result
rows. AI failures are recorded separately and do not turn a completed enrichment
analysis into a failed analysis.

Available AI profiles:

| Mode | Purpose |
| --- | --- |
| `summary` | Researcher-oriented biological pattern summary |
| `reviewer` | Statistical and interpretation-risk review |
| `caption` | Concise figure-caption style text |

CLI examples:

```bash
allenricher analyze \
  --input genes.txt \
  --species hsa \
  --databases GO,KEGG \
  --ai deepseek \
  --ai-mode summary \
  --ai-top-n 15 \
  --output results/ora_ai
```

```bash
allenricher analyze \
  --ranked-genes ranked_genes.tsv \
  --species hsa \
  --databases KEGG \
  --method gsea \
  --ai openai \
  --ai-mode reviewer \
  --ai-top-n 10 \
  --output results/gsea_ai
```

Configure AI credentials through CLI flags, environment variables, or YAML:

```bash
export DEEPSEEK_API_KEY="your-key"
allenricher analyze ... --ai deepseek
```

```yaml
ai_interpretation: true
ai_backend: deepseek
ai_backends:
  deepseek:
    api_key: "your-key"
    model: "deepseek-chat"
    enabled: true
  mock:
    enabled: true
```

Supported backends are `openai`, `claude`, `deepseek`, `glm`, `minimax`,
`ollama`, and `mock`. `mock` is intended for validation and tests, not for
scientific interpretation.

## Configuration Files

Generate a documented YAML configuration file:

```bash
allenricher config --output allenricher.yaml
```

Run with a configuration file:

```bash
allenricher analyze --config allenricher.yaml --input genes.txt
```

Explicit CLI arguments take precedence over values loaded from a configuration
file. See [`config.example.yaml`](config.example.yaml) for the maintained
example.

## Local Web Workbench

Install the API extra and start the local service:

```bash
python -m pip install -e ".[api,visualization]"
allenricher serve --host 127.0.0.1 --port 8000
```

Open [http://127.0.0.1:8000/](http://127.0.0.1:8000/) in a browser. OpenAPI
documentation is available at
[http://127.0.0.1:8000/docs](http://127.0.0.1:8000/docs).

The Web workbench provides:

- Method-specific input panels for ORA, GSEA, ssGSEA, and GSVA.
- File upload and paste modes for ORA gene lists.
- Ranked-table upload or differential-table-to-rank conversion for GSEA.
- Expression-matrix upload and interactive sample grouping for ssGSEA/GSVA.
- Species-specific database support checks through the registry.
- Disabled-but-visible unsupported databases, so users know which features need
  an additional species build.
- Plot style and palette controls with palette previews.
- Optional AI interpretation controls using configured backends.
- Result, figure, report, artifact, AI, and Methods-reference views.

Server-side settings can be supplied with:

```bash
allenricher serve --port 8000 --config allenricher.yaml
```

Useful environment variables:

```bash
ALLENRICHER_DATABASE_DIR=./database
ALLENRICHER_API_JOB_DIR=./api_jobs
ALLENRICHER_CONFIG=./allenricher.yaml
DEEPSEEK_API_KEY=your-key
OPENAI_API_KEY=your-key
```

## REST API

The API submits jobs to the same CLI workflow used by the command line and Web
workbench.

Main endpoints:

| Endpoint | Purpose |
| --- | --- |
| `GET /` | Web workbench |
| `GET /api/species` | Installed species and local database support |
| `GET /api/species/summary` | Registry summary counts |
| `GET /api/species/{species}/databases` | Databases available for one species |
| `GET /api/databases` | Public database catalog |
| `GET /api/ai/backends` | Configured AI backend status |
| `POST /api/analyze` | JSON analysis request |
| `POST /api/upload` | Multipart upload analysis request |
| `GET /api/status/{job_id}` | Job status and AI error state |
| `GET /api/results/{job_id}?format=json` | Result tables as JSON |
| `GET /api/results/{job_id}?format=tsv` | Combined result table as TSV |
| `GET /api/results/{job_id}/plot` | One plot file |
| `GET /api/results/{job_id}/report` | HTML report |
| `GET /api/results/{job_id}/ai-interpretation` | Structured AI interpretation |
| `GET /api/results/{job_id}/methods-reference` | Materials and Methods writing reference |
| `GET /api/results/{job_id}/artifacts` | Output files, logs, and reports |
| `GET /api/results/{job_id}/files/{file_path}` | Download one managed artifact |
| `DELETE /api/jobs/{job_id}` | Delete one completed or failed job |

Minimal JSON request:

```bash
curl -X POST http://127.0.0.1:8000/api/analyze \
  -H "Content-Type: application/json" \
  -d '{
    "genes": ["TP53", "BRCA1", "EGFR"],
    "species": "hsa",
    "databases": ["GO", "KEGG"],
    "method": "hypergeometric"
  }'
```

Minimal upload request:

```bash
curl -X POST http://127.0.0.1:8000/api/upload \
  -F method=hypergeometric \
  -F species=hsa \
  -F databases=GO,KEGG \
  -F gene_file=@genes.txt
```

## Output Files

A typical analysis directory contains:

```text
results/
|-- GO_enrichment.tsv
|-- KEGG_enrichment.tsv
|-- analysis_metadata.json
|-- report.html
|-- ai_interpretation.json
`-- plots/
    |-- GO_barplot.png
    `-- GO_lollipop.png
```

Result tables include stable term identifiers and readable term names. Where a
database supplies hierarchy metadata, hierarchy is retained in a separate
column. GSEA output follows the `fgsea` result contract; ssGSEA and GSVA output
is an activity matrix with gene sets as rows and samples as columns.

The HTML report includes generated result tables and figures, recorded run
metadata, optional AI interpretation, and an English Materials and Methods
writing reference based only on values stored for that run. Missing versions or
references are shown as `To be added` rather than inferred.

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
