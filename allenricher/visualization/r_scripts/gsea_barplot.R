#!/usr/bin/env Rscript
# gsea_barplot.R — 双向柱状图 (-log10(FDR)*sign(NES))
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
})
this_file <- sub("^--file=", "", commandArgs()[4])
source(file.path(dirname(this_file), "common.R"))

args <- parse_args()
tsv_path <- args$tsv
output <- args$output
top_n <- as.integer(args$top_n %||% 20)

df <- read_enrichment(tsv_path) %>%
  mutate(bar_value = -log10(pmax(Adjusted_P_Value, 1e-10)) * sign(NES),
         direction = ifelse(NES > 0, "Up-regulated", "Down-regulated")) %>%
  arrange(desc(abs(NES))) %>%
  head(top_n)

p <- ggplot(df, aes(x = reorder(Term_Name, bar_value), y = bar_value, fill = direction)) +
  geom_col() +
  coord_flip() +
  scale_fill_manual(values = c("Up-regulated" = "#E41A1C", "Down-regulated" = "#377EBA")) +
  labs(x = "", y = expression(-log[10](FDR) %*% sign(NES)), fill = "Regulation") +
  set_nature_theme() +
  theme(legend.position = "top")

save_plot(p, output, width = 10, height = max(6, top_n * 0.3))
