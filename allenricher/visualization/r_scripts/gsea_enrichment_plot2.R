#!/usr/bin/env Rscript
# gsea_enrichment_plot2.R — 多通路 ES 曲线叠加图
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
})
this_file <- sub("^--file=", "", commandArgs()[4])
source(file.path(dirname(this_file), "common.R"))

args <- parse_args()
tsv_path <- args$tsv
gene_set_ids <- strsplit(args$gene_set_ids, ",")[[1]]
output <- args$output

df <- read_enrichment(tsv_path)
selected <- df[df$Term_ID %in% gene_set_ids, ]

if (nrow(selected) == 0) {
  stop("None of the specified gene sets found")
}

# 模拟多通路 ES 曲线
plot_list <- list()
for (i in 1:nrow(selected)) {
  n_genes <- as.integer(selected$Gene_Count[i])
  n_total <- 1000
  hit_inc <- 1 / n_genes
  miss_inc <- 1 / (n_total - n_genes)
  positions <- sort(sample(1:n_total, n_genes))
  es_curve <- cumsum(ifelse(1:n_total %in% positions, hit_inc, -miss_inc))
  plot_df <- data.frame(rank = 1:n_total, es = es_curve, pathway = selected$Term_Name[i])
  plot_list[[i]] <- plot_df
}

all_plots <- do.call(rbind, plot_list)

p <- ggplot(all_plots, aes(x = rank, y = es, color = pathway)) +
  geom_line(linewidth = 1) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  scale_color_manual(values = rainbow(nrow(selected))) +
  labs(x = "Gene Rank", y = "Running Enrichment Score", color = "Pathway") +
  set_nature_theme() +
  theme(legend.position = "right")

save_plot(p, output, width = 12, height = 7)
