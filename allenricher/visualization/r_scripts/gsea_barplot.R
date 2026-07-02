#!/usr/bin/env Rscript
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
  prepare_gsea_terms(top_n = top_n, sort_by = "q", label_width = 36, max_label_chars = 88) %>%
  arrange(NES, Adjusted_P_Value)

df$Term_Label <- factor(df$Term_Label, levels = df$Term_Label)

p <- ggplot(df, aes(x = Term_Label, y = NES, fill = Direction)) +
  geom_hline(yintercept = 0, color = "#777777", linewidth = 0.35) +
  geom_col(width = 0.72, alpha = 0.92) +
  coord_flip() +
  scale_fill_direction(name = NULL) +
  scale_y_continuous(expand = expansion(mult = c(0.08, 0.12))) +
  labs(
    x = NULL,
    y = "Normalized enrichment score (NES)",
    subtitle = "Direction and magnitude of enriched pathways"
  ) +
  set_nature_theme(base_size = 9) +
  theme(
    axis.text.y = element_text(size = 8.4, lineheight = 0.94),
    legend.position = "top",
    panel.grid.major.x = element_line(color = AE_COL_GRID, linewidth = 0.25)
  )

save_plot(p, output, width = 8.4, height = max(4.8, nrow(df) * 0.32 + 1.8))
