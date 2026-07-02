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
top_n <- as.integer(args$top_n %||% 30)

df <- read_enrichment(tsv_path) %>%
  prepare_gsea_terms(top_n = top_n, sort_by = "abs_nes", label_width = 34, max_label_chars = 82) %>%
  arrange(NES, Adjusted_P_Value)

df$Term_Label <- factor(df$Term_Label, levels = df$Term_Label)

if (nrow(df) == 1) {
  label_text <- sprintf(
    "NES %.2f   FDR %.2g   Genes %s",
    df$NES[[1]], df$Adjusted_P_Value[[1]], df$Gene_Count[[1]]
  )
  p <- ggplot(df, aes(y = Term_Label)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "#A0A4A8", linewidth = 0.35) +
    geom_segment(aes(x = 0, xend = NES, yend = Term_Label, color = Direction), linewidth = 0.9, alpha = 0.86) +
    geom_point(aes(x = NES, fill = Direction), shape = 21, color = "white", size = 4.4, stroke = 0.45, alpha = 0.98) +
    geom_text(aes(x = NES, label = label_text), nudge_x = ifelse(df$NES[[1]] >= 0, 0.14, -0.14),
              hjust = ifelse(df$NES[[1]] >= 0, 0, 1), size = 3.0, color = AE_COL_TEXT) +
    scale_color_direction(guide = "none") +
    scale_fill_direction(guide = "none") +
    scale_x_continuous(expand = expansion(mult = c(0.18, 0.42))) +
    labs(x = "Normalized enrichment score (NES)", y = NULL, subtitle = "Single enriched pathway summary") +
    set_nature_theme(base_size = 9) +
    theme(axis.text.y = element_text(size = 8.2, lineheight = 0.94), legend.position = "none")
  save_plot(p, output, width = 7.4, height = 2.8)
} else {
  p <- ggplot(df, aes(y = Term_Label)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "#A0A4A8", linewidth = 0.35) +
    geom_segment(aes(x = 0, xend = NES, yend = Term_Label, color = Direction), linewidth = 0.65, alpha = 0.78) +
    geom_point(aes(x = NES, fill = Direction, size = nlogFDR), shape = 21, color = "white", stroke = 0.35, alpha = 0.96) +
    scale_color_direction(guide = "none") +
    scale_fill_direction(name = NULL) +
    scale_size_area(max_size = 5.8, name = expression(-log[10](FDR))) +
    guides(size = guide_legend(order = 2)) +
    scale_x_continuous(expand = expansion(mult = c(0.08, 0.12))) +
    labs(
      x = "Normalized enrichment score (NES)",
      y = NULL,
      subtitle = "NES ranking with FDR encoded by point size"
    ) +
    set_nature_theme(base_size = 9) +
    theme(
      axis.text.y = element_text(size = 7.8, lineheight = 0.94),
      legend.position = "right"
    )
  save_plot(p, output, width = 9.2, height = max(5.2, nrow(df) * 0.26 + 1.8))
}
