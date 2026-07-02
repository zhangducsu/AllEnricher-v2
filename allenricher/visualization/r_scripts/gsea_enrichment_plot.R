#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(gridExtra)
})
this_file <- sub("^--file=", "", commandArgs()[4])
source(file.path(dirname(this_file), "common.R"))

args <- parse_args()
tsv_path <- args$tsv
gene_set_id <- args$gene_set_id
running_es_path <- args$running_es
output <- args$output

df <- read_enrichment(tsv_path)
row_data <- df[df$Term_ID == gene_set_id, ]

if (nrow(row_data) == 0) {
  stop(paste("Gene set", gene_set_id, "not found in TSV"))
}

nes <- row_data$NES[1]
es <- row_data$ES[1]
fdr <- row_data$Adjusted_P_Value[1]
term_label <- wrap_label(shorten_label(clean_pathway_name(row_data$Term_Name[1]), 100), width = 58)
line_color <- ifelse(nes >= 0, AE_COL_UP, AE_COL_DOWN)

plot_data <- read_running_es(running_es_path)
plot_data <- plot_data[plot_data$Term_ID == gene_set_id, ]

if (nrow(plot_data) == 0) {
  stop(paste("Running ES data for", gene_set_id, "not found"))
}

hit_data <- plot_data[plot_data$Hit, ]
x_range <- range(plot_data$Rank, na.rm = TRUE)
if (nes >= 0) {
  peak_row <- plot_data[which.max(plot_data$Running_ES), ]
} else {
  peak_row <- plot_data[which.min(plot_data$Running_ES), ]
}
label_y <- ifelse(
  nes >= 0,
  max(plot_data$Running_ES, na.rm = TRUE),
  min(plot_data$Running_ES, na.rm = TRUE)
)
label_vjust <- ifelse(nes >= 0, 1, 0)

p1 <- ggplot(plot_data, aes(x = Rank, y = Running_ES)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "#A0A4A8", linewidth = 0.32) +
  geom_line(color = line_color, linewidth = 0.85) +
  geom_point(data = peak_row, aes(x = Rank, y = Running_ES), color = line_color, size = 2.2) +
  annotate(
    "label",
    x = x_range[1] + diff(x_range) * 0.04,
    y = label_y,
    hjust = 0,
    vjust = label_vjust,
    label = sprintf("NES %.2f | ES %.2f | FDR %.2g", nes, es, fdr),
    size = 3,
    linewidth = 0.18,
    fill = "white"
  ) +
  coord_cartesian(xlim = x_range, clip = "off") +
  labs(title = term_label, x = NULL, y = "Running ES") +
  set_nature_theme(base_size = 9) +
  theme(
    plot.title = element_text(size = 10.5, lineheight = 0.95),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )

p2 <- ggplot() +
  geom_segment(
    data = hit_data,
    aes(x = Rank, xend = Rank, y = 0, yend = 1),
    color = "#222222",
    linewidth = 0.28,
    alpha = 0.75
  ) +
  coord_cartesian(xlim = x_range, ylim = c(0, 1), clip = "off") +
  labs(x = NULL, y = "Hits") +
  set_nature_theme(base_size = 9) +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.line.y = element_blank(),
    panel.grid = element_blank()
  )

p3 <- ggplot(plot_data, aes(x = Rank, y = Weight)) +
  geom_hline(yintercept = 0, color = "#A0A4A8", linewidth = 0.25) +
  geom_col(aes(fill = Weight), width = 1, alpha = 0.86) +
  scale_fill_gradient2(low = AE_COL_DOWN, mid = "#F2F2F2", high = AE_COL_UP, midpoint = 0, guide = "none") +
  coord_cartesian(xlim = x_range, clip = "off") +
  labs(x = "Gene rank", y = "Rank metric") +
  set_nature_theme(base_size = 9) +
  theme(panel.grid.major.x = element_blank())

combined <- grid.arrange(p1, p2, p3, ncol = 1, heights = c(3.2, 0.75, 1.75))

save_plot(combined, output, width = 7.8, height = 6.6)
