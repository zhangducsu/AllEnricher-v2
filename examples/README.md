# AllEnricher v2 Examples

This directory contains the compact figure gallery used by the main README.
Most SVG figures are generated from deterministic example tables; selected PNG
previews are copied from E2E visual-review outputs for figure types that need
full workflow context, such as GSEA pathway networks and group-comparison plots.

Run from the repository root to regenerate the deterministic SVG examples:

```bash
python examples/run_examples.py
```

Generated and curated preview figures are written to `examples/output/figures/`.
