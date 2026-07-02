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
top_n <- min(as.integer(args$top_n %||% 18), 18)

df <- read_enrichment(tsv_path) %>%
  prepare_gsea_terms(top_n = top_n, sort_by = "q", label_width = 20, max_label_chars = 52)

gene_lists <- parse_gene_list(df$Genes)
names(gene_lists) <- df$Term_ID

pair_rows <- list()
for (i in seq_len(nrow(df) - 1)) {
  for (j in seq((i + 1), nrow(df))) {
    shared <- length(intersect(gene_lists[[i]], gene_lists[[j]]))
    union_size <- length(union(gene_lists[[i]], gene_lists[[j]]))
    jaccard <- ifelse(union_size > 0, shared / union_size, 0)
    if (shared > 0) {
      pair_rows[[length(pair_rows) + 1]] <- data.frame(
        From = df$Term_Label[i],
        To = df$Term_Label[j],
        Shared = shared,
        Jaccard = jaccard
      )
    }
  }
}

levels_y <- rev(df$Term_Label)
if (length(pair_rows) == 0) {
  df$Term_Label <- factor(df$Term_Label, levels = rev(df$Term_Label))
  summary_plot <- ggplot(df, aes(y = Term_Label)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "#A0A4A8", linewidth = 0.35) +
    geom_segment(aes(x = 0, xend = NES, yend = Term_Label, color = Direction), linewidth = 0.65, alpha = 0.82) +
    geom_point(aes(x = NES, fill = Direction, size = Gene_Count), shape = 21, color = "white", stroke = 0.35, alpha = 0.95) +
    geom_text(aes(x = NES, label = paste0("genes=", Gene_Count)),
              nudge_x = ifelse(df$NES >= 0, 0.12, -0.12),
              hjust = ifelse(df$NES >= 0, 0, 1), size = 2.6, color = AE_COL_TEXT) +
    scale_color_direction(guide = "none") +
    scale_fill_direction(guide = "none") +
    scale_size_area(max_size = 5.2, guide = "none") +
    scale_x_continuous(expand = expansion(mult = c(0.16, 0.28))) +
    labs(
      x = "Normalized enrichment score (NES)",
      y = NULL,
      subtitle = "No leading-gene overlap among selected pathways; showing pathway summaries"
    ) +
    set_nature_theme(base_size = 8) +
    theme(axis.text.y = element_text(size = 7.2, lineheight = 0.9), legend.position = "none")
  save_plot(summary_plot, output, width = 8.2, height = max(3.2, nrow(df) * 0.34 + 1.7))
  quit(save = "no", status = 0)
}

pairs <- bind_rows(pair_rows)
pairs$From <- factor(pairs$From, levels = df$Term_Label)
pairs$To <- factor(pairs$To, levels = levels_y)

p <- ggplot(pairs, aes(x = From, y = To)) +
  geom_point(aes(size = Shared, fill = Jaccard), shape = 21, color = "white", stroke = 0.35, alpha = 0.95) +
  scale_fill_gradient(low = "#D8E8F4", high = "#B2182B", name = "Jaccard") +
  scale_size_area(max_size = 7.2, name = "Shared genes") +
  coord_equal(clip = "off") +
  labs(
    x = NULL,
    y = NULL,
    subtitle = "Pathway overlap map: bubble color shows Jaccard similarity"
  ) +
  set_nature_theme(base_size = 8) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 6.6, lineheight = 0.88),
    axis.text.y = element_text(size = 6.6, lineheight = 0.88),
    axis.ticks = element_blank(),
    axis.line = element_blank(),
    panel.grid.major.x = element_line(color = "#F0F1F3", linewidth = 0.25),
    panel.grid.major.y = element_line(color = "#F0F1F3", linewidth = 0.25),
    legend.position = "right",
    plot.margin = margin(8, 14, 28, 8)
  )

save_plot(p, output, width = 8.6, height = 7.8)
