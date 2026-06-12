#!/usr/bin/env Rscript
# gsea_ridgeplot.R — 峰峦图/山脊图
suppressPackageStartupMessages({
  library(ggplot2)
  library(ggridges)
  library(dplyr)
})
this_file <- sub("^--file=", "", commandArgs()[4])
source(file.path(dirname(this_file), "common.R"))

args <- parse_args()
tsv_path <- args$tsv
output <- args$output
top_n <- as.integer(args$top_n %||% 15)

df <- read_enrichment(tsv_path) %>%
  arrange(desc(abs(NES))) %>%
  head(top_n)

# 使用 NES 模拟分布
p <- ggplot(df, aes(x = NES, y = reorder(Term_Name, NES), fill = after_stat(x))) +
  geom_density_ridges_gradient(scale = 3, rel_min_height = 0.01, from = -3, to = 3) +
  scale_fill_gradient2(low = "#377EBA", mid = "#FFFFCC", high = "#E41A1C", midpoint = 0) +
  labs(x = "Normalized Enrichment Score (NES)", y = "") +
  theme_ridges() +
  theme(legend.position = "none")

save_plot(p, output, width = 10, height = max(6, top_n * 0.4))
