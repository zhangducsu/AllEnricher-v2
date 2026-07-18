#!/usr/bin/env Rscript
this_file <- sub("^--file=", "", commandArgs()[4])
source(file.path(dirname(this_file), "common.R"))

args <- parse_args()
configure_plot_style(args)
tsv_path <- args$tsv
gene_set_ids <- trimws(strsplit(args$gene_set_ids, ",")[[1]])
running_es_path <- args$running_es
output <- args$output
plot_dpi <- as.integer(args$dpi %||% 300)

heat_colors <- ae_diverging_colors(256)

format_p <- function(x) {
  if (length(x) == 0 || is.na(x) || !is.finite(x)) return("NA")
  if (x < 0.001) return("<0.001")
  if (x < 0.01) return(formatC(x, format = "e", digits = 2))
  formatC(x, format = "f", digits = 3)
}

df <- read_enrichment(tsv_path)
selected <- df[match(gene_set_ids, df$Term_ID), , drop = FALSE]
selected <- selected[!is.na(selected$Term_ID), , drop = FALSE]
if (nrow(selected) == 0) stop("None of the specified gene sets found")
if (nrow(selected) > 8) stop("Multi-pathway GSEA supports at most 8 pathways")

all_plots <- read_running_es(running_es_path)
all_plots <- all_plots[all_plots$Term_ID %in% selected$Term_ID, , drop = FALSE]
if (nrow(all_plots) == 0) stop("No running ES data found for specified gene sets")

nsets <- nrow(selected)
used_colors <- ae_palette(nsets)
descriptions <- vapply(selected$Term_Name, clean_pathway_name, character(1))
strip_labels <- wrap_label(descriptions, width = 24)
curve_list <- vector("list", nsets)
hit_list <- vector("list", nsets)
legend_labels <- character(nsets)

for (i in seq_len(nsets)) {
  term_id <- selected$Term_ID[i]
  term_data <- all_plots[all_plots$Term_ID == term_id, , drop = FALSE]
  term_data <- term_data[order(term_data$Rank), , drop = FALSE]
  if (nrow(term_data) == 0) stop("No running ES data found for ", term_id)
  curve_list[[i]] <- as.numeric(term_data$Running_ES)
  hit_list[[i]] <- which(as.logical(term_data$Hit))
  legend_labels[i] <- paste0(
    descriptions[i], " | NES ", sprintf("%.2f", as.numeric(selected$NES[i])),
    " | P ", format_p(as.numeric(selected$p_value[i])),
    " | FDR ", format_p(as.numeric(selected$Adjusted_P_Value[i]))
  )
}

n <- length(curve_list[[1]])
fig_width <- 7.6
fig_height <- 3.10 + 0.46 * nsets
panel_heights <- c(2.7, rep(0.50, max(0, nsets - 1)), 0.95)
left_margin <- 7.2

draw_multi_gsea <- function() {
  par(family = AE_FONT_FAMILY)
  layout(matrix(seq_len(nsets + 1), nrow = nsets + 1), heights = panel_heights)
  dense_text_scale <- min(AE_TEXT_SCALE, 1.15)

  par(mar = c(0.10, left_margin, 0.12, 1.00) * AE_SPACING_SCALE, mgp = c(2.65, 0.70, 0),
      xaxs = "i", yaxs = "r", cex = dense_text_scale)
  global_min <- min(vapply(curve_list, min, numeric(1)))
  global_max <- max(vapply(curve_list, max, numeric(1)))
  data_span <- max(global_max - global_min, 0.10)
  ypad <- max(0.04, data_span * 0.06)
  legend_pad <- max(ypad, data_span * (0.08 + 0.055 * nsets))
  plot(
    NA,
    xlim = c(1, n), ylim = c(global_min - ypad, global_max + legend_pad),
    xlab = "", ylab = "Running ES", xaxt = "n", yaxt = "n", bty = "n",
    cex.lab = 1.18
  )
  axis(2, las = 1, cex.axis = 0.95)
  ae_background_grid(h = pretty(c(global_min - ypad, global_max + legend_pad)))
  abline(h = 0, lty = 2, lwd = ae_lwd(0.8))
  for (i in seq_len(nsets)) {
    lines(seq_len(n), curve_list[[i]], col = used_colors[i], lwd = ae_lwd(2.15))
  }
  legend(
    if (mean(as.numeric(selected$NES), na.rm = TRUE) >= 0) "topright" else "topleft",
    legend = legend_labels, col = used_colors, lwd = ae_lwd(2.2),
    bty = "n", cex = 0.78, x.intersp = 0.65, y.intersp = 0.85
  )
  ae_panel_box(0.75)

  heat_palette <- grDevices::colorRampPalette(heat_colors)(256)
  heat_edges <- seq(1, n + 1, length.out = 257)
  for (i in seq_len(nsets)) {
    is_last <- i == nsets
    par(
      mar = c(if (is_last) 2.65 else 0.03, left_margin, 0.01, 1.00) * AE_SPACING_SCALE,
      mgp = c(2.10, 0.68, 0), xaxs = "i", yaxs = "i"
    )
    plot(
      NA, xlim = c(1, n), ylim = c(0, 1),
      xlab = if (is_last) "Rank in Ordered Dataset" else "", ylab = "",
      xaxt = if (is_last) "s" else "n", yaxt = "n", bty = "n", cex.lab = 1.18
    )
    for (j in seq_len(256)) {
      rect(heat_edges[j], 0, heat_edges[j + 1], 0.30, col = heat_palette[j], border = NA)
    }
    segments(hit_list[[i]], 0.30, hit_list[[i]], 1.00, col = used_colors[i], lwd = ae_lwd(0.65))
    usr <- par("usr")
    text(
      usr[1] - 0.012 * (usr[2] - usr[1]), 0.53,
      labels = strip_labels[i], adj = 1, cex = 0.78,
      font = 2, col = used_colors[i], xpd = TRUE
    )
    ae_panel_box(0.75)
  }
}

ext <- tolower(tools::file_ext(output))
if (ext == "pdf") {
  pdf(output, width = fig_width, height = fig_height, useDingbats = FALSE, bg = "white")
} else if (ext == "svg") {
  svg(output, width = fig_width, height = fig_height, bg = "white")
} else {
  png(output, width = fig_width, height = fig_height, units = "in", res = plot_dpi, bg = "white")
}
draw_multi_gsea()
dev.off()
