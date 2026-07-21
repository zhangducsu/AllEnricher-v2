# AllEnricher competitor benchmark

This directory contains the preregistered, case-level competitor benchmark for
the Bioinformatics Application Note. It does not modify the production API,
CLI, or statistical implementation.

## Run

```text
python run_benchmark.py --config benchmark_matrix.yaml --output <new-run-directory>
python run_benchmark.py --config benchmark_matrix.yaml --output <new-run-directory> --case <case-id>
python run_benchmark.py --config benchmark_matrix.yaml --from-raw-run <archived-run> --output <new-replay-directory>
python make_publication_figures.py --run-dir <completed-run-directory> --paper-dir <Paper>
```

The output directory must not already exist. Each run archives normalized GMT
files, input hashes, raw tool output, requests and responses, standard output,
standard error, tool sessions, case status, system information, and a manifest.
A failed external service or unavailable runtime remains an explicit failed or
unavailable record; the runner never substitutes another tool.

## Result contracts

`normalized_results.tsv` stores one term per row using the common schema in the
benchmark plan. `benchmark_metrics.tsv` stores one row per dataset-by-database
comparison, so term rows are not treated as biological replicates.
`benchmark_metrics_detail.tsv` contains supporting counts, input mapping, and
acceptance diagnostics. `getenrich_workflow_audit.tsv` is a workflow audit, not
an independent numerical oracle.

## Interpretation boundaries

The exact ORA acceptance threshold applies only when both tools use the same
tested universe and multiple-testing family. WebGestaltR custom ORA uses the
annotated reference universe, while a g:Profiler custom GMT defines the service
annotation universe; those cases are retained as informative but incomparable.
The exact GSEA threshold applies to AllEnricher and a direct fgseaMultilevel
reference. WebGestaltR GSEA is descriptive because it uses a different
permutation implementation. Current clusterProfiler behavior is recorded from
the installed release rather than assumed from older interfaces.
