#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
})
this_file <- sub("^--file=", "", commandArgs()[4])
source(file.path(dirname(this_file), "common.R"))

args <- parse_args()
tsv_path <- args$tsv
gene_set_ids <- trimws(strsplit(args$gene_set_ids, ",")[[1]])
running_es_path <- args$running_es
output <- args$output

df <- read_enrichment(tsv_path)
selected <- df[df$Term_ID %in% gene_set_ids, ] %>%
  prepare_gsea_terms(top_n = length(gene_set_ids), sort_by = "abs_nes", label_width = 32, max_label_chars = 72)

if (nrow(selected) == 0) {
  stop("None of the specified gene sets found")
}

all_plots <- read_running_es(running_es_path)
all_plots <- all_plots[all_plots$Term_ID %in% selected$Term_ID, ]

if (nrow(all_plots) == 0) {
  stop("No running ES data found for specified gene sets")
}

term_names <- selected$Term_Label
names(term_names) <- selected$Term_ID
all_plots$Pathway <- term_names[all_plots$Term_ID]
all_plots$Pathway <- factor(all_plots$Pathway, levels = term_names)

palette <- pathway_palette(length(levels(all_plots$Pathway)))
names(palette) <- levels(all_plots$Pathway)
n_pathways <- length(levels(all_plots$Pathway))

p <- ggplot(all_plots, aes(x = Rank, y = Running_ES, color = Pathway)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "#A0A4A8", linewidth = 0.32) +
  geom_line(linewidth = 0.78, alpha = 0.92) +
  scale_color_manual(values = palette, name = NULL) +
  scale_x_continuous(expand = expansion(mult = c(0.02, 0.12))) +
  labs(
    x = "Gene rank",
    y = "Running enrichment score",
    subtitle = "Running-ES trajectories for the top pathways"
  ) +
  set_nature_theme(base_size = 9) +
  theme(
    panel.grid.major.x = element_line(color = AE_COL_GRID, linewidth = 0.25)
  )

if (n_pathways <= 7) {
  label_data <- all_plots %>%
    group_by(Pathway) %>%
    filter(Rank == max(Rank, na.rm = TRUE)) %>%
    slice_tail(n = 1) %>%
    ungroup()
  p <- p +
    geom_text(
      data = label_data,
      aes(label = Pathway),
      hjust = -0.02,
      size = 2.35,
      lineheight = 0.88,
      show.legend = FALSE
    ) +
    scale_x_continuous(expand = expansion(mult = c(0.02, 0.26))) +
    coord_cartesian(clip = "off") +
    theme(
      legend.position = "none",
      plot.margin = margin(8, 72, 8, 8)
    )
} else {
  p <- p +
    guides(color = guide_legend(nrow = 2, byrow = TRUE)) +
    theme(
      legend.position = "bottom",
      legend.text = element_text(size = 7.2, lineheight = 0.9),
      legend.key.width = grid::unit(1.1, "lines")
    )
}

save_plot(p, output, width = 8.4, height = 6.2)
