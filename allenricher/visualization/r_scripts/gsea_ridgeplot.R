#!/usr/bin/env Rscript

this_file <- sub("^--file=", "", commandArgs()[4])
source(file.path(dirname(this_file), "common.R"))

args <- parse_args()
configure_plot_style(args)
tsv_path <- args$tsv
running_es_path <- args$running_es %||% ""
output <- args$output
top_n <- as.integer(args$top_n %||% 10)
min_genes <- as.integer(args$min_genes %||% 5)
density_adjust <- as.numeric(args$density_adjust %||% 1.0)
plot_dpi <- as.integer(args$dpi %||% 300)

if (!is.finite(top_n) || top_n < 1) stop("top_n must be a positive integer")
if (!is.finite(min_genes) || min_genes < 1) stop("min_genes must be a positive integer")
if (!is.finite(density_adjust) || density_adjust <= 0) stop("density_adjust must be positive")

gsea <- read_enrichment(tsv_path)
running <- read_running_es(running_es_path)

sig_col <- if ("p_value" %in% colnames(gsea)) "p_value" else "Adjusted_P_Value"
gsea$significance <- as_num(gsea[[sig_col]], NA_real_)
gsea$Term_ID <- as.character(gsea$Term_ID)
gsea$Description <- clean_pathway_name(gsea$Term_Name)

running$Term_ID <- as.character(running$Term_ID)
running$Weight <- as_num(running$Weight, NA_real_)

matched_scores <- lapply(gsea$Term_ID, function(term_id) {
  values <- running$Weight[running$Term_ID == term_id & running$Hit]
  values[is.finite(values)]
})
gsea$matchedGeneCount <- lengths(matched_scores)
gsea$minusLog10Significance <- -log10(gsea$significance)

keep <- is.finite(gsea$significance) & gsea$significance > 0 &
  gsea$matchedGeneCount >= min_genes
gsea <- gsea[keep, , drop = FALSE]
matched_scores <- matched_scores[keep]

if (nrow(gsea) < 1) {
  stop("No pathways contain enough ranked genes for ridgeplot")
}

ordering <- order(gsea$significance, -abs(gsea$NES), na.last = TRUE)
gsea <- gsea[ordering, , drop = FALSE]
matched_scores <- matched_scores[ordering]
unique_terms <- !duplicated(gsea$Description)
gsea <- gsea[unique_terms, , drop = FALSE]
matched_scores <- matched_scores[unique_terms]

if (nrow(gsea) > top_n) {
  gsea <- head(gsea, top_n)
  matched_scores <- head(matched_scores, top_n)
}
rownames(gsea) <- NULL

all_scores <- running$Weight[is.finite(running$Weight)]
if (length(all_scores) < 2) stop("running-ES data contains fewer than two ranking statistics")
score_range <- range(all_scores)
score_padding <- diff(score_range) * 0.025
if (!is.finite(score_padding) || score_padding <= 0) score_padding <- 0.5
x_limits <- score_range + c(-score_padding, score_padding)

density_list <- lapply(matched_scores, function(values) {
  if (length(unique(values)) >= 2) {
    stats::density(
      values,
      from = x_limits[1],
      to = x_limits[2],
      n = 512,
      adjust = density_adjust
    )
  } else {
    bandwidth <- max(diff(x_limits) / 45, 0.05)
    x <- seq(x_limits[1], x_limits[2], length.out = 512)
    list(x = x, y = stats::dnorm(x, mean = values[1], sd = bandwidth))
  }
})

color_values <- gsea$minusLog10Significance
color_min <- floor(min(color_values, na.rm = TRUE))
color_max <- ceiling(max(color_values, na.rm = TRUE))
if (!is.finite(color_min)) color_min <- 0
if (!is.finite(color_max) || color_max <= color_min) color_max <- color_min + 1

ridge_palette <- ae_gradient_colors(256)

value_to_color <- function(value) {
  scaled <- (value - color_min) / (color_max - color_min)
  ridge_palette[1 + floor(pmax(0, pmin(1, scaled)) * 255)]
}
ridge_colors <- vapply(color_values, value_to_color, character(1))

draw_ridgeplot <- function() {
  n_terms <- nrow(gsea)
  layout(matrix(c(1, 2), nrow = 1), widths = c(6.3, 1.65))

  par(
    mar = c(3.9, 13.5, 0.55, 0.4) * AE_SPACING_SCALE, mgp = c(2.5, 0.7, 0),
    xaxs = "i", yaxs = "i", family = AE_FONT_FAMILY, cex = AE_TEXT_SCALE
  )
  plot(
    NA, xlim = x_limits, ylim = c(0.52, n_terms + 0.88),
    xlab = "", ylab = "", axes = FALSE, bty = "n"
  )
  x_ticks <- pretty(x_limits, n = 5)
  x_ticks <- x_ticks[x_ticks >= x_limits[1] & x_ticks <= x_limits[2]]
  axis(1, at = x_ticks, labels = format(x_ticks, trim = TRUE), cex.axis = 1.0,
       lwd = ae_lwd(0.7), tck = -0.016)
  ae_background_grid(v = x_ticks)
  mtext("Ranking statistic", side = 1, line = 2.35, cex = 1.12)

  for (i in seq_len(n_terms)) {
    baseline <- n_terms - i + 1
    density_object <- density_list[[i]]
    density_y <- density_object$y
    scaled_y <- if (length(density_y) && is.finite(max(density_y)) && max(density_y) > 0) {
      density_y / max(density_y) * 0.78
    } else {
      rep(0, length(density_y))
    }

    segments(x_limits[1], baseline, x_limits[2], baseline, col = "black", lwd = ae_lwd(0.65))
    polygon(
      c(density_object$x, rev(density_object$x)),
      c(rep(baseline, length(density_object$x)), rev(baseline + scaled_y)),
      col = ridge_colors[i], border = NA
    )
    lines(density_object$x, baseline + scaled_y, col = "black", lwd = ae_lwd(0.72))
    segments(
      matched_scores[[i]], baseline - 0.105,
      matched_scores[[i]], baseline + 0.025,
      col = ridge_palette[length(ridge_palette)], lwd = ae_lwd(0.72)
    )
    text(
      x_limits[1], baseline,
      labels = wrap_label(gsea$Description[i], width = 47),
      adj = c(1.03, 0.5), cex = 0.92, xpd = NA
    )
  }
  ae_panel_box(0.72)

  par(mar = c(3.9, 0.65, 3.5, 1.0) * AE_SPACING_SCALE, xaxs = "i", yaxs = "i",
      family = AE_FONT_FAMILY, cex = AE_TEXT_SCALE)
  plot(NA, xlim = c(0, 1.45), ylim = c(0, 1), axes = FALSE, xlab = "", ylab = "", bty = "n")
  legend_left <- 0.16
  legend_right <- 0.49
  legend_bottom <- 0.18
  legend_top <- 0.74
  text(0.05, legend_top + 0.085, labels = paste0("-log10(", sig_col, ")"), adj = 0, cex = 1.05)
  legend_edges <- seq(legend_bottom, legend_top, length.out = length(ridge_palette) + 1)
  for (i in seq_along(ridge_palette)) {
    rect(legend_left, legend_edges[i], legend_right, legend_edges[i + 1],
         col = ridge_palette[i], border = NA)
  }
  rect(legend_left, legend_bottom, legend_right, legend_top, border = "grey35", lwd = ae_lwd(0.55))
  legend_ticks <- pretty(c(color_min, color_max), n = 5)
  legend_ticks <- legend_ticks[legend_ticks >= color_min & legend_ticks <= color_max]
  for (tick in legend_ticks) {
    yy <- legend_bottom + (tick - color_min) / (color_max - color_min) *
      (legend_top - legend_bottom)
    segments(legend_right, yy, legend_right + 0.07, yy, lwd = ae_lwd(0.55))
    text(legend_right + 0.11, yy, labels = format(tick, trim = TRUE), adj = 0, cex = 0.88)
  }
}

dir.create(dirname(output), recursive = TRUE, showWarnings = FALSE)
figure_width <- 8.8
figure_height <- max(4.8, 1.4 + 0.48 * nrow(gsea))
ext <- tolower(tools::file_ext(output))
if (ext == "png") {
  png(output, width = figure_width, height = figure_height, units = "in", res = plot_dpi, bg = "white")
} else if (ext == "pdf") {
  pdf(output, width = figure_width, height = figure_height, useDingbats = FALSE, bg = "white")
} else if (ext == "svg") {
  svg(output, width = figure_width, height = figure_height, bg = "white")
} else {
  stop(paste("Unsupported output format:", ext))
}
draw_ridgeplot()
dev.off()
message("Saved: ", output)
