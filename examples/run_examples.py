"""Generate the small README figure gallery from fixed example tables."""

from pathlib import Path

import numpy as np
import pandas as pd

from allenricher.visualization.barplot import plot_barplot
from allenricher.visualization.gsea_plots import plot_gsea_enrichment, plot_gsea_lollipop
from allenricher.visualization.gsva_plots import plot_pathway_heatmap, plot_sample_correlation

ROOT = Path(__file__).resolve().parent
DATA = ROOT / "data"
FIGURES = ROOT / "output" / "figures"


def _ranked_example():
    genes = [f"G{i:03d}" for i in range(1, 121)]
    featured = ["CDK1", "CCNB1", "CDC20", "AURKB", "PLK1", "MCM2", "PCNA", "RAD51"]
    ranked = featured[:5] + genes[:55] + featured[5:] + genes[55:]
    weights = np.linspace(3.0, -3.0, len(ranked))
    return ranked, dict(zip(ranked, weights))


def main():
    FIGURES.mkdir(parents=True, exist_ok=True)

    ora = pd.read_csv(DATA / "ora_results.tsv", sep="\t")
    plot_barplot(
        ora,
        top_n=8,
        style="nature",
        palette="okabe_ito",
        title="Example ORA: KEGG enriched pathways",
        output_file=str(FIGURES / "ora_kegg_barplot.svg"),
        term_col="term",
        term_id_col="term_id",
        hierarchy_col="hierarchy",
        qvalue_col="qvalue",
        gene_count_col="gene_count",
        rich_factor_col="rich_factor",
    )
    plot_gsea_lollipop(
        ora.rename(columns={"term": "Term_Name", "qvalue": "Adjusted_P_Value", "gene_count": "Gene_Count", "rich_factor": "EnrichFactor"}),
        top_n=8,
        style="nature",
        palette="okabe_ito",
        title="Example ORA lollipop: KEGG enriched pathways",
        output_file=str(FIGURES / "ora_kegg_lollipop.svg"),
    )

    gsea = pd.read_csv(DATA / "gsea_results.tsv", sep="\t")
    plot_gsea_lollipop(
        gsea,
        top_n=6,
        style="nature",
        palette="colorbrewer_rdbu",
        title="Example GSEA: positive and negative NES pathways",
        output_file=str(FIGURES / "gsea_kegg_lollipop.svg"),
    )
    ranked, weights = _ranked_example()
    gene_set = {"CDK1", "CCNB1", "CDC20", "AURKB", "PLK1", "MCM2", "PCNA", "RAD51"}
    plot_gsea_enrichment(
        ranked,
        weights,
        gene_set,
        es=0.71,
        nes=2.24,
        pvalue=0.0008,
        padj=0.006,
        title="Cell cycle",
        style="nature",
        palette="colorbrewer_rdbu",
        output_file=str(FIGURES / "gsea_cell_cycle_enrichment.svg"),
    )

    activity = pd.read_csv(DATA / "activity_scores.tsv", sep="\t", index_col=0)
    groups = pd.DataFrame(
        {"Group": ["Control", "Control", "Control", "Treatment", "Treatment", "Treatment"]},
        index=activity.columns,
    )
    plot_pathway_heatmap(
        activity,
        annotation_col=groups,
        title="Example ssGSEA/GSVA pathway activity",
        style="nature",
        palette="colorbrewer_rdbu",
        output_file=str(FIGURES / "activity_heatmap.svg"),
        top_n=10,
    )
    plot_sample_correlation(
        activity,
        annotation_col=groups,
        title="Example sample correlation",
        style="nature",
        palette="colorbrewer_rdbu",
        output_file=str(FIGURES / "sample_correlation.svg"),
    )

    print(f"Wrote example figures to {FIGURES}")


if __name__ == "__main__":
    main()
