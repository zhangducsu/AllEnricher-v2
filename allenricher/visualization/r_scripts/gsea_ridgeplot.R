#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(ggplot2)
  library(ggridges)
  library(dplyr)
})
this_file <- sub("^--file=", "", commandArgs()[4])
source(file.path(dirname(this_file), "common.R"))

args <- parse_args()
tsv_path <- args$tsv
running_es_path <- args$running_es %||% ""
output <- args$output
top_n <- as.integer(args$top_n %||% 15)

df <- read_enrichment(tsv_path) %>%
  prepare_gsea_terms(top_n = top_n, sort_by = "abs_nes", label_width = 32, max_label_chars = 76) %>%
  arrange(NES)

df$Term_Label <- factor(df$Term_Label, levels = df$Term_Label)

if (running_es_path != "" && file.exists(running_es_path)) {
  hit_data <- read_running_es(running_es_path) %>%
    filter(Term_ID %in% df$Term_ID, Hit) %>%
    left_join(df[, c("Term_ID", "Term_Label", "Direction")], by = "Term_ID")

  if (nrow(hit_data) > 0) {
    max_rank <- max(hit_data$Rank, na.rm = TRUE)
    hit_data$Rank_Percent <- 100 * hit_data$Rank / max_rank
    hit_data$Term_Label <- factor(hit_data$Term_Label, levels = levels(df$Term_Label))
    counts <- table(hit_data$Term_Label)

    if (all(counts >= 3)) {
      p <- ggplot(hit_data, aes(x = Rank_Percent, y = Term_Label, fill = Direction)) +
        geom_density_ridges(
          scale = 1.25,
          bandwidth = 5,
          alpha = 0.86,
          color = "white",
          linewidth = 0.25,
          rel_min_height = 0.01
        ) +
        scale_fill_direction(name = NULL) +
        labs(
          x = "Gene rank percentile",
          y = NULL,
          subtitle = "Distribution of hit genes across the ranked list"
        ) +
        set_nature_theme(base_size = 9) +
        theme(
          legend.position = "top",
          axis.text.y = element_text(size = 7.8, lineheight = 0.94),
          panel.grid.major.x = element_line(color = AE_COL_GRID, linewidth = 0.25)
        )
    } else {
      p <- ggplot(hit_data, aes(x = Rank_Percent, y = Term_Label, color = Direction)) +
        geom_point(size = 1.6, alpha = 0.72, position = position_jitter(height = 0.08, width = 0)) +
        scale_color_direction(name = NULL) +
        labs(
          x = "Gene rank percentile",
          y = NULL,
          subtitle = "Hit-gene positions; density is skipped for sparse pathways"
        ) +
        set_nature_theme(base_size = 9) +
        theme(
          legend.position = "top",
          axis.text.y = element_text(size = 7.8, lineheight = 0.94)
        )
    }

    save_plot(p, output, width = 8.6, height = max(5.2, nrow(df) * 0.36 + 1.8))
    quit(save = "no", status = 0)
  }
}

message("running_es is missing or has no hits; rendering deterministic NES summary instead")
p <- ggplot(df, aes(y = Term_Label)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "#A0A4A8", linewidth = 0.35) +
  geom_segment(aes(x = 0, xend = NES, yend = Term_Label, color = Direction), linewidth = 0.7) +
  geom_point(aes(x = NES, fill = Direction), shape = 21, color = "white", stroke = 0.35, size = 3.2) +
  scale_color_direction(guide = "none") +
  scale_fill_direction(name = NULL) +
  labs(
    x = "Normalized enrichment score (NES)",
    y = NULL,
    subtitle = "Fallback NES summary; pass running-ES data for hit-rank ridges"
  ) +
  set_nature_theme(base_size = 9) +
  theme(
    legend.position = "top",
    axis.text.y = element_text(size = 7.8, lineheight = 0.94)
  )

save_plot(p, output, width = 8.6, height = max(5.2, nrow(df) * 0.34 + 1.8))
