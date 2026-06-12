#!/usr/bin/env Rscript
# gsea_nes_plot.R — NES 排序散点图
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
})
this_file <- sub("^--file=", "", commandArgs()[4])
source(file.path(dirname(this_file), "common.R"))

args <- parse_args()
tsv_path <- args$tsv
output <- args$output

df <- read_enrichment(tsv_path) %>%
  mutate(significance = ifelse(Adjusted_P_Value < 0.05, "Significant", "Not Significant"))

p <- ggplot(df, aes(x = reorder(Term_Name, NES), y = NES, color = significance)) +
  geom_point(size = 3) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  coord_flip() +
  scale_color_manual(values = c("Significant" = "#E41A1C", "Not Significant" = "gray70")) +
  labs(x = "", y = "Normalized Enrichment Score (NES)", color = "") +
  set_nature_theme()

save_plot(p, output, width = 12, height = min(max(8, nrow(df) * 0.15), 40))
