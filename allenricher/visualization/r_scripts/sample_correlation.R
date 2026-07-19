#!/usr/bin/env Rscript

this_file <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])
source(file.path(dirname(this_file), "common.R"))

raw_args <- commandArgs(trailingOnly = TRUE)
args <- list()
index <- 1L
while (index <= length(raw_args)) {
  args[[sub("^--", "", raw_args[index])]] <- if (index < length(raw_args)) raw_args[index + 1L] else ""
  index <- index + 2L
}
configure_plot_style(args)

score_df <- read.delim(args$scores, check.names = FALSE, stringsAsFactors = FALSE)
pathways <- as.character(score_df[[1]])
scores <- as.matrix(data.frame(
  lapply(score_df[-1], function(value) suppressWarnings(as.numeric(value))),
  check.names = FALSE
))
rownames(scores) <- make.unique(pathways)
method <- if (is.null(args$method)) "pearson" else args$method
plot_dpi <- if (is.null(args$dpi)) 300L else as.integer(args$dpi)
if (!method %in% c("pearson", "spearman")) stop("method must be pearson or spearman")

correlation <- cor(scores, method = method, use = "pairwise.complete.obs")
correlation[!is.finite(correlation)] <- 0
diag(correlation) <- 1
distance <- 1 - correlation
distance[distance < 0] <- 0
tree <- hclust(as.dist(distance), method = "average")
order <- tree$order
correlation <- correlation[order, order, drop = FALSE]
samples <- colnames(correlation)

metadata <- read.delim(args$metadata, check.names = FALSE, stringsAsFactors = FALSE)
sample_col <- colnames(metadata)[1]
if (anyDuplicated(metadata[[sample_col]])) stop("metadata sample IDs must be unique")
rownames(metadata) <- as.character(metadata[[sample_col]])
metadata <- metadata[samples, setdiff(colnames(metadata), sample_col), drop = FALSE]
annotation_columns <- colnames(metadata)

annotation_count <- sum(vapply(metadata, function(values) {
  length(unique(as.character(values)))
}, integer(1)))
annotation_palette <- ae_palette(max(2, annotation_count))
annotation_maps <- list()
offset <- 1L
for (column in annotation_columns) {
  categories <- unique(as.character(metadata[[column]]))
  colors <- annotation_palette[offset + seq_along(categories) - 1L]
  names(colors) <- categories
  annotation_maps[[column]] <- colors
  offset <- offset + length(categories)
}

ellipse_points <- function(value, center_x, center_y) {
  theta <- seq(0, 2 * pi, length.out = 180)
  points <- rbind(cos(theta), sin(theta))
  scales <- diag(c(sqrt(max(1 + value, 0)), sqrt(max(1 - value, 0))))
  angle <- pi / 4
  rotation <- matrix(c(cos(angle), -sin(angle), sin(angle), cos(angle)), 2, 2, byrow = TRUE)
  result <- rotation %*% scales %*% points * 0.40
  result[1, ] <- result[1, ] + center_x
  result[2, ] <- result[2, ] + center_y
  result
}

heat_palette <- ae_diverging_colors(201)
map_color <- function(value) heat_palette[max(1, min(201, round((value + 1) / 2 * 200) + 1))]
n <- ncol(correlation)
n_annotations <- ncol(metadata)
heat_units <- max(4.9, n * 0.295)
label_units <- max(0.78, min(1.60, max(nchar(samples)) * 0.070))
legend_units <- 1.25 * max(1, AE_TEXT_SCALE)
widths <- c(0.65, heat_units, label_units, legend_units)
heights <- c(0.72, max(0.06, 0.11 * n_annotations), heat_units)
figure_width <- max(9.2, 1.9 + heat_units + label_units + legend_units)
figure_height <- figure_width * sum(heights) / sum(widths)
title_x <- (widths[1] + widths[2] / 2) / sum(widths)

draw_plot <- function() {
  par(oma = c(0, 0, 0.9, 0) * AE_SPACING_SCALE, family = AE_FONT_FAMILY, cex = AE_TEXT_SCALE, lwd = AE_LINE_SCALE)
  layout(
    matrix(c(0, 1, 0, 3, 0, 2, 0, 3, 4, 5, 6, 3), nrow = 3, byrow = TRUE),
    widths = widths, heights = heights
  )
  par(cex = AE_TEXT_SCALE)

  par(mar = rep(0.05, 4))
  plot(as.dendrogram(tree), axes = FALSE, leaflab = "none", edgePar = list(lwd = ae_lwd(0.8)))

  par(mar = rep(0.05, 4), xaxs = "i", yaxs = "i")
  plot.new(); plot.window(c(-0.5, n - 0.5), c(max(n_annotations, 1) - 0.5, -0.5))
  if (n_annotations) for (row in seq_len(n_annotations)) {
    mapping <- annotation_maps[[row]]
    for (column in seq_len(n)) {
      rect(column - 1.5, row - 1.5, column - 0.5, row - 0.5,
           col = mapping[[as.character(metadata[column, row])]], border = NA)
    }
  }

  par(mar = c(0.05, 0.2, 0.05, 0.2), xpd = NA)
  plot.new(); plot.window(c(0, 1), c(0, 1))
  legend_rows <- sum(vapply(annotation_maps, function(mapping) 1L + length(mapping), integer(1)))
  legend_step <- min(0.045, 0.24 / figure_height, 0.54 / max(1, legend_rows))
  heading_gap <- min(0.05, 0.28 / figure_height)
  section_gap <- min(0.06, 0.34 / figure_height)
  legend_y <- 0.94
  for (column in annotation_columns) {
    text(0, legend_y, column, adj = c(0, 1), cex = 0.82, font = 2)
    mapping <- annotation_maps[[column]]
    categories <- names(mapping)
    first_y <- legend_y - heading_gap
    for (index in seq_along(categories)) {
      category <- categories[index]
      y_pos <- first_y - legend_step * (index - 1)
      rect(0.02, y_pos - 0.009, 0.095, y_pos + 0.009,
           col = mapping[[category]], border = NA)
      text(0.12, y_pos, category, adj = c(0, 0.5), cex = 0.72)
    }
    legend_y <- first_y - legend_step * max(0, length(categories) - 1) - section_gap
  }
  text(0, legend_y, "Correlation", adj = c(0, 1), cex = 0.82, font = 2)
  top <- legend_y - 0.28 / figure_height
  bottom <- max(0.04, top - min(0.18, 1.35 / figure_height))
  edges <- seq(bottom, top, length.out = 202)
  for (i in seq_len(201)) rect(0.02, edges[i], 0.095, edges[i + 1], col = heat_palette[i], border = NA)
  rect(0.02, bottom, 0.095, top, border = "grey30", lwd = ae_lwd(0.6))
  for (tick in c(-1, 0, 1)) {
    y <- bottom + (tick + 1) / 2 * (top - bottom)
    segments(0.095, y, 0.115, y)
    text(0.13, y, tick, adj = c(0, 0.5), cex = 0.70)
  }

  par(mar = c(2.2, 0.05, 0.05, 0.05) * AE_SPACING_SCALE)
  plot(as.dendrogram(tree), horiz = TRUE, axes = FALSE, leaflab = "none", edgePar = list(lwd = ae_lwd(0.8)))

  par(mar = c(2.2, 0.05, 0.05, 0.05) * AE_SPACING_SCALE, xpd = NA, xaxs = "i", yaxs = "i")
  plot.new(); plot.window(c(-0.5, n - 0.5), c(n - 0.5, -0.5), asp = 1)
  for (row in seq_len(n)) for (column in seq_len(n)) {
    if (row == column) {
      text(column - 1, row - 1, "1", cex = 0.72, font = 2, col = "#404040")
    } else if (row > column) {
      ellipse <- ellipse_points(correlation[row, column], column - 1, row - 1)
      polygon(ellipse[1, ], ellipse[2, ], col = map_color(correlation[row, column]), border = NA)
    } else {
      text(column - 1, row - 1, sprintf("%.2f", correlation[row, column]),
           cex = 0.66, font = 2, col = "#404040")
    }
  }
  ae_matrix_box(0.8)
  for (column in seq_len(n)) {
    text(column - 1, n - 0.4, samples[column], srt = 90, adj = c(1, 0.5), cex = 0.74)
  }

  par(mar = c(2.2, 0.05, 0.05, 0.05) * AE_SPACING_SCALE, xaxs = "i", yaxs = "i")
  plot.new(); plot.window(c(0, 1), c(n - 0.5, -0.5))
  for (row in seq_len(n)) text(0, row - 1, samples[row], adj = c(0, 0.5), cex = 0.78)
  mtext("Sample Correlation", side = 3, outer = TRUE, at = title_x,
        line = -0.25, cex = 1.05, font = 2)
}

output <- args$output
dir.create(dirname(output), recursive = TRUE, showWarnings = FALSE)
extension <- tolower(tools::file_ext(output))
if (extension == "png") {
  png(output, width = figure_width, height = figure_height, units = "in", res = plot_dpi, bg = "white")
} else if (extension == "pdf") {
  pdf(output, width = figure_width, height = figure_height, useDingbats = FALSE, bg = "white")
} else if (extension == "svg") {
  svg(output, width = figure_width, height = figure_height, bg = "white")
} else {
  stop("Unsupported output format: ", extension)
}
draw_plot()
dev.off()
message("Saved: ", output)
