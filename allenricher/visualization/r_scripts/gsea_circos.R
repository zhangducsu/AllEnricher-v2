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
top_n <- min(as.integer(args$top_n %||% 20), 20)

df <- read_enrichment(tsv_path) %>%
  prepare_gsea_terms(top_n = top_n, sort_by = "q", label_width = 18, max_label_chars = 42) %>%
  arrange(NES) %>%
  mutate(
    id = row_number(),
    abs_nes = abs(NES),
    angle = 90 - 360 * (id - 0.5) / n(),
    hjust = ifelse(angle < -90, 1, 0),
    angle = ifelse(angle < -90, angle + 180, angle)
  )

df$Term_Label <- factor(df$Term_Label, levels = df$Term_Label)
inner_radius <- max(df$abs_nes, na.rm = TRUE) * 0.35

p <- ggplot(df, aes(x = Term_Label, y = abs_nes, fill = Direction)) +
  geom_col(width = 0.78, alpha = 0.92) +
  geom_hline(yintercept = 0, color = "white", linewidth = 0.4) +
  geom_text(
    aes(y = abs_nes + inner_radius * 0.16, label = Term_Label, angle = angle, hjust = hjust),
    size = 2.25,
    lineheight = 0.82,
    color = AE_COL_TEXT
  ) +
  scale_fill_direction(name = NULL) +
  scale_y_continuous(limits = c(-inner_radius, max(df$abs_nes, na.rm = TRUE) * 1.28)) +
  coord_polar(clip = "off") +
  labs(subtitle = "Circular summary of top pathways by absolute NES") +
  theme_void(base_size = 9) +
  theme(
    text = element_text(color = AE_COL_TEXT),
    legend.position = "bottom",
    plot.subtitle = element_text(size = 9, color = "#555555", hjust = 0),
    plot.margin = margin(18, 18, 18, 18)
  )

save_plot(p, output, width = 7.4, height = 7.4)
