#!/usr/bin/env Rscript
this_file <- sub("^--file=", "", commandArgs()[4])
source(file.path(dirname(this_file), "common.R"))

args <- parse_args()
configure_plot_style(args)
tsv_path <- args$tsv
gene_set_id <- args$gene_set_id
running_es_path <- args$running_es
output <- args$output
png_dpi <- as.integer(args$dpi %||% 300)

heat_col <- ae_diverging_colors(256)
rank_col <- heat_col
figure_width <- 6.60
figure_height <- 4.40
panel_heights <- c(2.85, 0.60, 1.22)
max_rank_bars <- 3000
title_cex <- 1.48 * AE_TEXT_SCALE
axis_title_cex <- 1.28 * AE_TEXT_SCALE
axis_text_cex <- 1.04 * AE_TEXT_SCALE
stats_cex <- 1.00 * AE_TEXT_SCALE

format_p <- function(x) {
  if (length(x) == 0 || is.na(x) || !is.finite(x)) return("NA")
  if (x < 0.001) return("< 0.001")
  if (x < 0.01) return(formatC(x, format = "e", digits = 2))
  formatC(x, format = "f", digits = 3)
}

gradient_index <- function(values, limits, n_colors) {
  vmin <- limits[1]
  vmax <- limits[2]
  if (!is.finite(vmin) || !is.finite(vmax) || vmax <= vmin) {
    return(rep(ceiling(n_colors / 2), length(values)))
  }
  scaled <- (values - vmin) / (vmax - vmin)
  scaled <- pmax(0, pmin(1, scaled))
  round(scaled * (n_colors - 1)) + 1
}

draw_gsea_base <- function(plot_data, row_data) {
  par(family = AE_FONT_FAMILY)
  plot_data <- plot_data[order(plot_data$Rank), , drop = FALSE]
  n <- nrow(plot_data)
  x <- seq_len(n)
  running_es <- as.numeric(plot_data$Running_ES)
  hit_position <- which(as.logical(plot_data$Hit))
  gene_weight <- as.numeric(plot_data$Weight)

  nes <- as.numeric(row_data$NES[1])
  pvalue <- as.numeric(row_data$p_value[1])
  padj <- as.numeric(row_data$Adjusted_P_Value[1])
  title <- clean_pathway_name(row_data$Term_Name[1])

  layout(matrix(c(1, 2, 3), nrow = 3, ncol = 1), heights = panel_heights)

  par(mar = c(0.20, 4.65, 1.75, 0.25) * AE_SPACING_SCALE, mgp = c(2.70, 0.72, 0), tcl = -0.25, xaxs = "i", yaxs = "r")
  y_range <- range(running_es, finite = TRUE)
  y_pad <- max(0.045, diff(y_range) * 0.065)
  plot(NA, xlim = c(1, n), ylim = c(y_range[1] - y_pad, y_range[2] + y_pad),
       xlab = "", ylab = "Running Enrichment Score", xaxt = "n", yaxt = "n",
       bty = "n", cex.lab = axis_title_cex, main = title, cex.main = title_cex, font.main = 2)
  axis(2, las = 1, cex.axis = axis_text_cex, tck = -0.018)
  ae_background_grid(h = pretty(c(y_range[1] - y_pad, y_range[2] + y_pad)))
  abline(h = 0, lty = 2, lwd = ae_lwd(0.8), col = "black")
  curve_color <- if (is.finite(nes) && nes < 0) AE_COL_DOWN else AE_COL_UP
  lines(x, running_es, col = curve_color, lwd = ae_lwd(2.15))
  ae_panel_box(0.75)
  usr <- par("usr")
  text(
    x = usr[2] - 0.015 * (usr[2] - usr[1]),
    y = usr[4] - 0.030 * (usr[4] - usr[3]),
    labels = paste0("NES: ", sprintf("%.2f", nes), "\n",
                    "P value: ", format_p(pvalue), "\n",
                    "Adjusted P value: ", format_p(padj)),
    adj = c(1, 1), cex = stats_cex, font = 3, col = "grey20"
  )

  par(mar = c(0.08, 4.65, 0.02, 0.25) * AE_SPACING_SCALE, xaxs = "i", yaxs = "i")
  plot(NA, xlim = c(1, n), ylim = c(0, 1), xlab = "", ylab = "", xaxt = "n", yaxt = "n", bty = "n")
  heat_palette <- grDevices::colorRampPalette(heat_col)(256)
  bin_edges <- seq(1, n + 1, length.out = 257)
  for (i in seq_len(256)) {
    rect(bin_edges[i], 0, bin_edges[i + 1], 0.28, col = heat_palette[i], border = NA)
  }
  segments(hit_position, 0.28, hit_position, 1.00, col = "black", lwd = ae_lwd(0.60))
  ae_panel_box(0.75)

  par(mar = c(3.40, 4.65, 0.08, 0.25) * AE_SPACING_SCALE, mgp = c(2.45, 0.72, 0), tcl = -0.25, xaxs = "i", yaxs = "r")
  score_limit <- max(abs(gene_weight), na.rm = TRUE)
  if (!is.finite(score_limit) || score_limit <= 0) score_limit <- 1
  plot(NA, xlim = c(1, n), ylim = c(-score_limit, score_limit),
       xlab = "Rank in Ordered Dataset", ylab = "Ranked List",
       xaxt = "n", yaxt = "n", bty = "n", cex.lab = axis_title_cex)
  axis(1, cex.axis = axis_text_cex, tck = -0.018)
  axis(2, las = 1, cex.axis = axis_text_cex, tck = -0.018)
  ae_background_grid(h = pretty(c(-score_limit, score_limit)))
  abline(h = 0, lty = 2, lwd = ae_lwd(0.8), col = "black")
  rank_palette <- grDevices::colorRampPalette(rank_col)(256)
  rank_idx <- gradient_index(gene_weight, c(-score_limit, score_limit), length(rank_palette))
  step <- max(1, ceiling(n / max_rank_bars))
  selected <- seq(1, n, by = step)
  rect(pmax(1, selected - step / 2), pmin(0, gene_weight[selected]),
       pmin(n, selected + step / 2), pmax(0, gene_weight[selected]),
       col = rank_palette[rank_idx[selected]], border = NA)
  ae_panel_box(0.75)
}

save_one <- function(output, plot_data, row_data) {
  ext <- tolower(tools::file_ext(output))
  if (ext == "pdf") {
    pdf(output, width = figure_width, height = figure_height, useDingbats = FALSE, bg = "white")
  } else if (ext == "svg") {
    svg(output, width = figure_width, height = figure_height, bg = "white")
  } else {
    png(output, width = figure_width, height = figure_height, units = "in", res = png_dpi, bg = "white")
  }
  on.exit(dev.off(), add = TRUE)
  draw_gsea_base(plot_data, row_data)
}

df <- read_enrichment(tsv_path)
row_data <- df[df$Term_ID == gene_set_id, , drop = FALSE]
if (nrow(row_data) == 0) {
  stop(paste("Gene set", gene_set_id, "not found in TSV"))
}

plot_data <- read_running_es(running_es_path)
plot_data <- plot_data[plot_data$Term_ID == gene_set_id, , drop = FALSE]
if (nrow(plot_data) == 0) {
  stop(paste("Running ES data for", gene_set_id, "not found"))
}

save_one(output, plot_data, row_data)
