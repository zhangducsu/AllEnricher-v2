#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
})
this_file <- sub("^--file=", "", commandArgs()[4])
source(file.path(dirname(this_file), "common.R"))

args <- parse_args()
tsv_path <- args$tsv
output <- args$output
top_n <- as.integer(args$top_n %||% 10)

df <- read_enrichment(tsv_path) %>%
  prepare_gsea_terms(top_n = top_n, sort_by = "q", label_width = 30, max_label_chars = 68) %>%
  arrange(desc(abs(NES)))

gene_sets <- parse_gene_list(df$Genes)
names(gene_sets) <- df$Term_ID
all_genes <- unique(unlist(gene_sets))
if (length(all_genes) == 0) {
  stop("No genes available for cnetplot")
}

gene_freq <- sort(table(unlist(gene_sets)), decreasing = TRUE)
selected_genes <- names(gene_freq[gene_freq >= 2])
if (nrow(df) == 1) {
  selected_genes <- names(gene_freq)[seq_len(min(18, length(gene_freq)))]
} else if (length(selected_genes) < 8) {
  selected_genes <- names(gene_freq)[seq_len(min(35, length(gene_freq)))]
} else {
  selected_genes <- selected_genes[seq_len(min(40, length(selected_genes)))]
}

edge_rows <- list()
for (i in seq_along(gene_sets)) {
  term_id <- names(gene_sets)[i]
  genes <- intersect(gene_sets[[i]], selected_genes)
  if (length(genes) == 0) {
    next
  }
  edge_rows[[length(edge_rows) + 1]] <- data.frame(Term_ID = term_id, Gene = genes)
}
edges <- bind_rows(edge_rows)
if (nrow(edges) == 0) {
  stop("No gene-pathway edges available for cnetplot")
}

path_nodes <- df %>%
  mutate(x = 0, y = rev(seq_len(n()))) %>%
  select(Term_ID, Term_Label, NES, Gene_Count, Direction, x, y)

gene_order <- names(sort(table(edges$Gene), decreasing = TRUE))
if (nrow(path_nodes) == 1) {
  path_nodes$x <- 0
  path_nodes$y <- 0
  grid_cols <- min(6, length(gene_order))
  grid_col <- (seq_along(gene_order) - 1) %% grid_cols
  grid_row <- floor((seq_along(gene_order) - 1) / grid_cols)
  gene_nodes <- data.frame(
    Gene = gene_order,
    x = 0.72 + grid_col * 0.23,
    y = 0.40 - grid_row * 0.28
  )
  x_limits <- c(-0.58, 2.30)
  plot_height <- max(3.6, 2.5 + max(grid_row) * 0.34)
} else {
  gene_nodes <- data.frame(
    Gene = gene_order,
    x = 1,
    y = seq(from = max(path_nodes$y), to = 1, length.out = length(gene_order))
  )
  x_limits <- c(-0.72, 1.42)
  plot_height <- max(6.2, max(nrow(path_nodes), nrow(gene_nodes)) * 0.18 + 2.2)
}

edges <- edges %>%
  left_join(path_nodes[, c("Term_ID", "x", "y")], by = "Term_ID")
gene_nodes_for_join <- gene_nodes
colnames(gene_nodes_for_join) <- c("Gene", "xend", "yend")
edges <- edges %>%
  left_join(gene_nodes_for_join, by = "Gene")

p <- ggplot() +
  geom_segment(
    data = edges,
    aes(x = x, y = y, xend = xend, yend = yend),
    color = "#B7BDC5",
    linewidth = 0.28,
    alpha = 0.52
  ) +
  geom_point(
    data = path_nodes,
    aes(x = x, y = y, size = Gene_Count, fill = NES),
    shape = 21,
    color = "white",
    stroke = 0.45
  ) +
  geom_point(data = gene_nodes, aes(x = x, y = y), size = 1.55, color = "#30363D") +
  geom_text(
    data = path_nodes,
    aes(x = x - 0.035, y = y, label = Term_Label),
    hjust = 1,
    size = 2.55,
    lineheight = 0.88,
    color = AE_COL_TEXT
  ) +
  geom_text(
    data = gene_nodes,
    aes(x = x + 0.035, y = y, label = Gene),
    hjust = 0,
    size = 2.25,
    color = "#30363D"
  ) +
  scale_fill_nes(name = "NES") +
  scale_size_area(max_size = 7, name = "Gene count") +
  coord_cartesian(xlim = x_limits, clip = "off") +
  labs(subtitle = "Concept-gene network; shared genes are prioritized for readability") +
  theme_void(base_size = 9) +
  theme(
    text = element_text(color = AE_COL_TEXT),
    plot.subtitle = element_text(size = 9, color = "#555555", hjust = 0),
    legend.position = "right",
    plot.margin = margin(12, 36, 12, 36)
  )

save_plot(p, output, width = ifelse(nrow(path_nodes) == 1, 8.6, 10.2), height = plot_height)
