#!/usr/bin/env Rscript

this_file <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])
source(file.path(dirname(this_file), "common.R"))

args_raw <- commandArgs(trailingOnly = TRUE)
args <- list()
index <- 1L
while (index <= length(args_raw)) {
  key <- sub("^--", "", args_raw[index])
  args[[key]] <- if (index < length(args_raw)) args_raw[index + 1L] else ""
  index <- index + 2L
}
configure_plot_style(args)

scores_path <- args$scores
metadata_path <- args$metadata
output <- args$output
scale_mode <- if (is.null(args$scale)) "row" else args$scale
top_n <- if (is.null(args$top_n)) 40L else as.integer(args$top_n)
analysis_method <- tolower(if (is.null(args$analysis_method)) "" else args$analysis_method)
plot_title <- if (analysis_method == "ssgsea") {
  "ssGSEA Pathway Activity"
} else if (analysis_method == "gsva") {
  "GSVA Pathway Activity"
} else {
  "Pathway Activity"
}
plot_dpi <- if (is.null(args$dpi)) 300L else as.integer(args$dpi)
if (!scale_mode %in% c("row", "column", "none")) stop("scale must be row, column, or none")
if (!is.finite(top_n) || top_n < 1L) stop("top_n must be a positive integer")

score_df <- read.delim(scores_path, check.names = FALSE, stringsAsFactors = FALSE)
metadata <- read.delim(metadata_path, check.names = FALSE, stringsAsFactors = FALSE)
pathways <- as.character(score_df[[1]])
samples <- colnames(score_df)[-1]
scores <- as.matrix(data.frame(lapply(score_df[-1], as.numeric), check.names = FALSE))
rownames(scores) <- make.unique(pathways)
colnames(scores) <- samples

sample_col <- colnames(metadata)[1]
rownames(metadata) <- as.character(metadata[[sample_col]])
metadata <- metadata[samples, setdiff(colnames(metadata), sample_col), drop = FALSE]

for (row in seq_len(nrow(scores))) {
  good <- is.finite(scores[row, ])
  scores[row, !good] <- if (any(good)) median(scores[row, good]) else 0
}
if (nrow(scores) > top_n) {
  row_variance <- apply(scores, 1, var, na.rm = TRUE)
  row_variance[!is.finite(row_variance)] <- 0
  display_score <- row_variance
  if (ncol(metadata)) {
    group_col <- if ("Group" %in% colnames(metadata)) "Group" else colnames(metadata)[1]
    group_values <- as.character(metadata[[group_col]])
    groups <- unique(group_values[!is.na(group_values)])
    if (length(groups) >= 2L) {
      group_means <- sapply(groups, function(group) {
        rowMeans(scores[, group_values == group, drop = FALSE], na.rm = TRUE)
      })
      display_score <- apply(group_means, 1, function(values) diff(range(values, na.rm = TRUE)))
      display_score[!is.finite(display_score)] <- 0
    }
  }
  selected <- head(order(-display_score, -row_variance, rownames(scores)), top_n)
  scores <- scores[selected, , drop = FALSE]
}
if (scale_mode == "row") {
  means <- rowMeans(scores)
  sds <- apply(scores, 1, sd)
  sds[!is.finite(sds) | sds == 0] <- 1
  scores <- sweep(sweep(scores, 1, means, "-"), 1, sds, "/")
} else if (scale_mode == "column") {
  means <- colMeans(scores)
  sds <- apply(scores, 2, sd)
  sds[!is.finite(sds) | sds == 0] <- 1
  scores <- sweep(sweep(scores, 2, means, "-"), 2, sds, "/")
}

corr_dist <- function(values) {
  correlation <- suppressWarnings(cor(t(values), use = "pairwise.complete.obs"))
  correlation[!is.finite(correlation)] <- 0
  diag(correlation) <- 1
  distance <- 1 - correlation
  distance[distance < 0] <- 0
  as.dist(distance)
}
row_tree <- if (nrow(scores) > 1) hclust(corr_dist(scores), method = "average") else NULL
col_tree <- if (ncol(scores) > 1) hclust(corr_dist(t(scores)), method = "average") else NULL
row_order <- if (is.null(row_tree)) seq_len(nrow(scores)) else row_tree$order
col_order <- if (is.null(col_tree)) seq_len(ncol(scores)) else col_tree$order
scores <- scores[row_order, col_order, drop = FALSE]
metadata <- metadata[col_order, , drop = FALSE]
display_pathways <- clean_pathway_name(rownames(scores))
max_label_chars <- max(nchar(display_pathways, type = "width"), na.rm = TRUE)
n_pathways <- nrow(scores)
n_samples <- ncol(scores)
n_annotations <- max(1, ncol(metadata))
label_panel_width <- max(1.8, min(4.2, 0.055 * max_label_chars))
legend_panel_width <- 1.25 * max(1, AE_TEXT_SCALE)
heatmap_panel_width <- max(3.2, min(9.2, 0.23 * n_samples))
heatmap_panel_height <- max(3.0, min(9.0, 0.28 * n_pathways))
sample_label_capacity <- max(1L, floor(heatmap_panel_width * 72 / 8.0))
pathway_label_capacity <- max(1L, floor(heatmap_panel_height * 72 / 8.0))
sample_label_stride <- max(1L, ceiling(n_samples / sample_label_capacity))
pathway_label_stride <- max(1L, ceiling(n_pathways / pathway_label_capacity))
sample_label_indices <- seq.int(1L, n_samples, by = sample_label_stride)
pathway_label_indices <- seq.int(1L, n_pathways, by = pathway_label_stride)
sample_label_cex <- max(0.38, min(0.72, heatmap_panel_width * 72 / n_samples / 10))
pathway_label_cex <- max(0.38, min(0.84, heatmap_panel_height * 72 / n_pathways / 13))

heat_palette <- ae_diverging_colors(256)
limit <- as.numeric(quantile(abs(scores), 0.98, na.rm = TRUE))
if (!is.finite(limit) || limit <= 0) limit <- 1
annotation_palette <- ae_palette(max(1, sum(vapply(metadata, function(values) {
  length(unique(as.character(values)))
}, integer(1)))))
annotation_offset <- 0L
annotation_maps <- lapply(metadata, function(values) {
  categories <- unique(as.character(values))
  indices <- annotation_offset + seq_along(categories)
  annotation_offset <<- annotation_offset + length(categories)
  setNames(annotation_palette[indices], categories)
})
color_index <- function(value) 1L + floor(pmax(0, pmin(1, (value + limit) / (2 * limit))) * 255)

rounded_tile <- function(x, y, half_x, half_y, radius_x, radius_y, fill) {
  angles <- c(
    seq(-pi / 2, 0, length.out = 5), seq(0, pi / 2, length.out = 5),
    seq(pi / 2, pi, length.out = 5), seq(pi, 3 * pi / 2, length.out = 5)
  )
  centers_x <- rep(c(x + half_x - radius_x, x + half_x - radius_x,
                     x - half_x + radius_x, x - half_x + radius_x), each = 5)
  centers_y <- rep(c(y - half_y + radius_y, y + half_y - radius_y,
                     y + half_y - radius_y, y - half_y + radius_y), each = 5)
  polygon(centers_x + radius_x * cos(angles), centers_y + radius_y * sin(angles),
          col = fill, border = NA)
}

draw_heatmap <- function() {
  par(oma = c(0, 0, 1.0, 0) * AE_SPACING_SCALE, family = AE_FONT_FAMILY, cex = AE_TEXT_SCALE, lwd = AE_LINE_SCALE)
  layout(
    matrix(c(0, 1, 0, 3, 0, 2, 0, 3, 4, 5, 6, 3), nrow = 3, byrow = TRUE),
    widths = c(0.70, heatmap_panel_width, label_panel_width, legend_panel_width),
    heights = c(0.55, 0.14 * n_annotations, heatmap_panel_height)
  )
  par(cex = AE_TEXT_SCALE)

  par(mar = rep(0, 4))
  if (is.null(col_tree)) plot.new() else plot(as.dendrogram(col_tree), axes = FALSE, leaflab = "none")

  par(mar = rep(0, 4), xaxs = "i", yaxs = "i")
  plot.new(); plot.window(c(0.5, n_samples + 0.5), c(0.5, n_annotations + 0.5))
  if (ncol(metadata)) for (annotation in seq_len(ncol(metadata))) {
    y <- n_annotations - annotation + 1
    mapping <- annotation_maps[[annotation]]
    for (sample in seq_len(n_samples)) {
      rect(sample - 0.5, y - 0.5, sample + 0.5, y + 0.5,
           col = mapping[[as.character(metadata[sample, annotation])]], border = NA)
    }
  }

  par(mar = c(0, 0.3, 0, 0) * AE_SPACING_SCALE)
  plot.new(); plot.window(c(0, 1), c(0, 1))
  legend_rows <- sum(vapply(annotation_maps, function(mapping) 1L + length(mapping), integer(1)))
  legend_step <- min(0.045, 0.24 / height, 0.54 / max(1, legend_rows))
  heading_gap <- min(0.05, 0.28 / height)
  section_gap <- min(0.06, 0.34 / height)
  legend_y <- 0.94
  if (ncol(metadata)) for (annotation in seq_len(ncol(metadata))) {
    text(0, legend_y, colnames(metadata)[annotation], adj = c(0, 1), font = 2, cex = 0.72)
    mapping <- annotation_maps[[annotation]]
    categories <- names(mapping)
    legend_y <- legend_y - heading_gap
    for (index in seq_along(categories)) {
      category <- categories[index]
      y_pos <- legend_y - legend_step * (index - 1)
      rect(0.02, y_pos - 0.009, 0.095, y_pos + 0.009,
           col = mapping[[category]], border = NA)
      text(0.12, y_pos, category, adj = c(0, 0.5), cex = 0.67)
    }
    legend_y <- legend_y - legend_step * max(0, length(categories) - 1) - section_gap
  }

  color_title <- if (scale_mode == "row") "Activity score\n(row Z-score)" else "Activity score"
  text(0, legend_y, color_title, adj = c(0, 1), cex = 0.68, font = 2)
  bar_top <- legend_y - if (scale_mode == "row") 0.42 / height else 0.24 / height
  bar_bottom <- max(0.04, bar_top - min(0.18, 1.35 / height))
  bar_left <- 0.02; bar_right <- 0.095
  edges <- seq(bar_bottom, bar_top, length.out = 257)
  for (i in seq_len(256)) rect(bar_left, edges[i], bar_right, edges[i + 1], col = heat_palette[i], border = NA)
  rect(bar_left, bar_bottom, bar_right, bar_top, border = "grey30", lwd = ae_lwd(0.6))
  ticks <- c(-limit, 0, limit)
  tick_y <- bar_bottom + (ticks + limit) / (2 * limit) * (bar_top - bar_bottom)
  segments(bar_right, tick_y, bar_right + 0.025, tick_y, lwd = ae_lwd(0.6))
  tick_labels <- formatC(ticks, format = "fg", digits = 2)
  for (index in seq_along(ticks)) {
    text(bar_right + 0.045, tick_y[index], tick_labels[index],
         adj = c(0, c(0, 0.5, 1)[index]), cex = 0.72)
  }

  par(mar = c(7.0, 0, 0.1, 0) * AE_SPACING_SCALE)
  if (is.null(row_tree)) plot.new() else plot(rev(as.dendrogram(row_tree)), horiz = TRUE, axes = FALSE, leaflab = "none")

  par(mar = c(7.0, 0.1, 0.1, 0.1) * AE_SPACING_SCALE, xaxs = "i", yaxs = "i")
  plot.new(); plot.window(c(0.5, n_samples + 0.5), c(0.5, n_pathways + 0.5))
  display <- scores[n_pathways:1, , drop = FALSE]
  abline(v = seq(0.5, n_samples + 0.5), h = seq(0.5, n_pathways + 0.5), col = AE_COL_GRID, lwd = ae_lwd(0.4))
  for (row in seq_len(n_pathways)) for (column in seq_len(n_samples)) {
    rounded_tile(column, row, 0.46, 0.44, 0.10, 0.10,
                 heat_palette[color_index(display[row, column])])
  }
  axis(
    1,
    sample_label_indices,
    colnames(scores)[sample_label_indices],
    las = 2,
    tick = FALSE,
    cex.axis = sample_label_cex
  )
  ae_matrix_box(0.8)

  par(mar = c(7.0, 0, 0.1, 0) * AE_SPACING_SCALE, xaxs = "i", yaxs = "i")
  plot.new(); plot.window(c(0, 1), c(0.5, n_pathways + 0.5))
  display_labels <- rev(display_pathways)
  text(
    0.02,
    pathway_label_indices,
    display_labels[pathway_label_indices],
    adj = c(0, 0.5),
    cex = pathway_label_cex
  )

  mtext(plot_title, side = 3, outer = TRUE,
        line = 0.05, cex = 1.05, font = 2)
}

dir.create(dirname(output), recursive = TRUE, showWarnings = FALSE)
width <- max(8.4, min(15.0, 1.2 + 0.70 + heatmap_panel_width + label_panel_width + legend_panel_width))
height <- max(5.2, min(12.0, 1.55 + 0.55 + 0.14 * n_annotations + heatmap_panel_height))
ext <- tolower(tools::file_ext(output))
if (ext == "png") {
  png(output, width = width, height = height, units = "in", res = plot_dpi, bg = "white")
} else if (ext == "pdf") {
  pdf(output, width = width, height = height, useDingbats = FALSE, bg = "white")
} else if (ext == "svg") {
  svg(output, width = width, height = height, bg = "white")
} else {
  stop("Unsupported output format: ", ext)
}
draw_heatmap()
dev.off()
message("Saved: ", output)
