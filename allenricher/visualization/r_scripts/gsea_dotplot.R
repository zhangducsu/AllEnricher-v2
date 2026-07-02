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

p <- ggplot(df, aes(x = NES, y = Term_Label))
if (min(df$NES, na.rm = TRUE) < 0 && max(df$NES, na.rm = TRUE) > 0) {
  p <- p + geom_vline(xintercept = 0, linetype = "dashed", color = "#A0A4A8", linewidth = 0.35)
}
p <- p +
  geom_point(aes(size = Gene_Count, color = nlogFDR), alpha = 0.92) +
  scale_color_nlogfdr(name = expression(-log[10](FDR))) +
  scale_size_area(max_size = 6.4, name = "Gene count") +
  guides(
    color = guide_colorbar(order = 1, barheight = grid::unit(30, "pt")),
    size = guide_legend(order = 2)
  ) +
  scale_x_continuous(expand = expansion(mult = c(0.08, 0.12))) +
  labs(
    x = "Normalized enrichment score (NES)",
    y = NULL,
    subtitle = "Top pathways ranked by FDR, with point size showing gene count"
  ) +
  set_nature_theme(base_size = 9) +
  theme(
    axis.text.y = element_text(size = 8.4, lineheight = 0.94),
    panel.grid.major.x = element_line(color = AE_COL_GRID, linewidth = 0.25),
    legend.box = "vertical"
  )

save_plot(p, output, width = 8.6, height = max(4.8, nrow(df) * 0.32 + 1.8))
