#!/usr/bin/env Rscript

this_file <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])
source(file.path(dirname(this_file), "common.R"))

raw_args <- commandArgs(trailingOnly = TRUE)
args <- list()
i <- 1L
while (i <= length(raw_args)) {
  args[[sub("^--", "", raw_args[i])]] <- if (i < length(raw_args)) raw_args[i + 1L] else ""
  i <- i + 2L
}
configure_plot_style(args)

scores_file <- args$scores
metadata_file <- args$metadata
output_file <- args$output
statistics_file <- args$statistics
top_n <- if (is.null(args$top_n)) 6L else as.integer(args$top_n)
global_method <- if (is.null(args$global_test)) "kruskal" else args$global_test
pairwise_method <- if (is.null(args$pairwise_test)) "mannwhitney" else args$pairwise_test
p_adjust_method <- if (is.null(args$p_adjust)) "BH" else args$p_adjust
global_label <- if (global_method == "anova") "ANOVA" else "Kruskal-Wallis"
pairwise_label <- if (pairwise_method == "ttest") "Welch t-test" else "Wilcoxon"
comparison_mode <- if (is.null(args$comparison_mode)) "all" else args$comparison_mode
reference_group <- if (is.null(args$reference_group)) "" else args$reference_group
ncols <- if (is.null(args$ncols)) 3L else as.integer(args$ncols)
plot_dpi <- if (is.null(args$dpi)) 300L else as.integer(args$dpi)

if (is.null(scores_file) || is.null(metadata_file) || is.null(output_file) || is.null(statistics_file)) {
  stop("scores, metadata, output and statistics are required")
}

scores_df <- read.delim(scores_file, check.names = FALSE, stringsAsFactors = FALSE)
pathways <- as.character(scores_df[[1]])
scores <- as.matrix(data.frame(lapply(scores_df[-1], as.numeric), check.names = FALSE))
rownames(scores) <- make.unique(pathways)

metadata <- read.delim(metadata_file, check.names = FALSE, stringsAsFactors = FALSE)
sample_col <- colnames(metadata)[1]
rownames(metadata) <- as.character(metadata[[sample_col]])
if (!"Group" %in% colnames(metadata)) stop("metadata must contain Group")
missing_samples <- setdiff(colnames(scores), rownames(metadata))
if (length(missing_samples)) stop(paste("metadata missing samples:", paste(missing_samples, collapse = ", ")))
metadata <- metadata[colnames(scores), , drop = FALSE]
group_order <- unique(as.character(metadata$Group))
if (length(group_order) < 2L) stop("group comparison requires at least two groups")
two_groups <- length(group_order) == 2L

colors <- setNames(ae_palette(length(group_order)), group_order)

format_p <- function(p) {
  if (!is.finite(p)) return("NA")
  if (p < 0.001) formatC(p, format = "e", digits = 1) else formatC(p, format = "f", digits = 3)
}

global_p <- function(groups) {
  groups <- lapply(groups, function(x) x[is.finite(x)])
  groups <- groups[lengths(groups) > 0L]
  if (length(groups) < 2L) return(NA_real_)
  values <- unlist(groups)
  factors <- factor(rep(names(groups), lengths(groups)), levels = names(groups))
  result <- try(
    if (global_method == "anova") summary(aov(values ~ factors))[[1]][["Pr(>F)"]][1]
    else kruskal.test(values, factors)$p.value,
    silent = TRUE
  )
  if (inherits(result, "try-error")) NA_real_ else as.numeric(result)
}

pair_p <- function(left, right) {
  left <- left[is.finite(left)]
  right <- right[is.finite(right)]
  if (!length(left) || !length(right)) return(NA_real_)
  result <- try(
    if (pairwise_method == "ttest") t.test(left, right, var.equal = FALSE)$p.value
    else wilcox.test(left, right, exact = FALSE)$p.value,
    silent = TRUE
  )
  if (inherits(result, "try-error")) NA_real_ else as.numeric(result)
}

pathway_groups <- function(pathway) {
  row <- scores[pathway, ]
  groups <- lapply(group_order, function(group) row[as.character(metadata$Group) == group])
  names(groups) <- group_order
  groups
}

pvalues <- vapply(rownames(scores), function(pathway) {
  groups <- pathway_groups(pathway)
  if (two_groups) pair_p(groups[[1]], groups[[2]]) else global_p(groups)
}, numeric(1))
selected <- rownames(scores)[head(order(pvalues, na.last = TRUE), max(1L, min(top_n, nrow(scores))))]

if (comparison_mode == "none") {
  comparisons <- list()
} else if (comparison_mode == "reference") {
  reference <- if (nzchar(reference_group)) reference_group else group_order[1]
  if (!reference %in% group_order) stop("reference group is absent")
  comparisons <- lapply(setdiff(group_order, reference), function(group) c(reference, group))
} else {
  comparisons <- combn(group_order, 2L, simplify = FALSE)
}

half_violin <- function(values, x, color, ygrid, width = 0.28) {
  values <- values[is.finite(values)]
  if (!length(values)) return()
  density_result <- if (length(unique(values)) >= 2L) {
    density(values, from = min(ygrid), to = max(ygrid), n = length(ygrid))
  } else {
    list(x = ygrid, y = dnorm(ygrid, values[1], max(diff(range(ygrid)) / 30, 0.05)))
  }
  keep <- density_result$y >= max(density_result$y) * 0.02
  if (sum(keep) < 2L) return()
  density_result$x <- density_result$x[keep]
  density_result$y <- density_result$y[keep]
  width_values <- density_result$y / max(density_result$y) * width
  polygon(c(x - width_values, rep(x, length(width_values))),
          c(density_result$x, rev(density_result$x)),
          col = adjustcolor(color, alpha.f = 0.18), border = NA)
  lines(x - width_values, density_result$x, col = AE_COL_NEUTRAL, lwd = ae_lwd(0.85))
}

draw_bracket <- function(x1, x2, y, height, label) {
  segments(x1, y, x1, y + height, lwd = ae_lwd(0.85))
  segments(x1, y + height, x2, y + height, lwd = ae_lwd(0.85))
  segments(x2, y + height, x2, y, lwd = ae_lwd(0.85))
  text((x1 + x2) / 2, y + height, label, adj = c(0.5, -0.12), cex = 0.80)
}

wrap_panel_label <- function(label, width = 28L, max_lines = 2L) {
  lines <- strwrap(clean_pathway_name(label), width = width, simplify = FALSE)[[1]]
  if (!length(lines)) return("")
  if (length(lines) > max_lines) {
    remainder <- paste(lines[max_lines:length(lines)], collapse = " ")
    remainder <- if (nchar(remainder) > width) {
      paste0(substr(remainder, 1L, max(1L, width - 3L)), "...")
    } else {
      remainder
    }
    lines <- c(lines[seq_len(max_lines - 1L)], remainder)
  }
  paste(lines, collapse = "\n")
}

ncols <- max(1L, min(ncols, length(selected)))
nrows <- ceiling(length(selected) / ncols)
figure_width <- max(7.8, 3.05 * ncols)
figure_height <- max(4.7, 3.60 * nrows)
stats_rows <- list()
stats_index <- 1L

draw_plot <- function() {
  par(mfrow = c(nrows, ncols), mar = c(2.8, 3.15, 2.90, 0.35) * AE_SPACING_SCALE,
      mgp = c(1.85, 0.52, 0), family = AE_FONT_FAMILY, cex = AE_TEXT_SCALE, cex.axis = 0.94,
      cex.lab = 1.03, tcl = -0.22)
  for (pathway in selected) {
    groups <- pathway_groups(pathway)
    all_values <- unlist(lapply(groups, function(x) x[is.finite(x)]))
    ymin <- min(all_values)
    ymax <- max(all_values)
    value_range <- max(ymax - ymin, 1)
    pathway_global_p <- if (two_groups) NA_real_ else global_p(groups)
    raw_pair <- vapply(comparisons, function(pair) pair_p(groups[[pair[1]]], groups[[pair[2]]]), numeric(1))
    adjusted_pair <- if (two_groups || !length(raw_pair) || tolower(p_adjust_method) == "none") raw_pair
                     else p.adjust(raw_pair, method = p_adjust_method)
    step <- value_range * 0.13
    bracket_height <- value_range * 0.022
    start <- ymax + value_range * 0.08
    global_padding <- if (two_groups) 0.10 else 0.18
    upper <- start + max(0, length(comparisons) - 1) * step + bracket_height + value_range * global_padding
    lower <- ymin - value_range * 0.06
    plot(NA, xlim = c(0.45, length(group_order) + 0.55), ylim = c(lower, upper),
         xaxt = "n", xlab = "", ylab = "Activity score", bty = "n")
    ae_background_grid(h = pretty(c(lower, ymax)), lwd = 0.7)
    ygrid <- seq(ymin - value_range * 0.04, ymax + value_range * 0.04, length.out = 256)

    for (group_index in seq_along(group_order)) {
      group <- group_order[group_index]
      finite <- groups[[group]][is.finite(groups[[group]])]
      half_violin(finite, group_index - 0.05, colors[group], ygrid)
      stripchart(finite, at = group_index - 0.03, method = "jitter", jitter = 0.08,
                 vertical = TRUE, pch = 21, bg = adjustcolor(colors[group], alpha.f = 0.72),
                 col = "white", cex = 0.80, add = TRUE)
      boxplot(finite, at = group_index + 0.10, boxwex = 0.22, add = TRUE,
              axes = FALSE, outline = FALSE, col = adjustcolor(colors[group], alpha.f = 0.82),
              border = "#555555", medcol = "white", medlwd = ae_lwd(1.4),
              boxlwd = ae_lwd(0.9), whisklwd = ae_lwd(0.9), staplelwd = ae_lwd(0.9),
              whiskcol = "#555555", staplecol = "#555555")
    }
    axis(
      1,
      at = seq_along(group_order),
      labels = group_order,
      tick = FALSE,
      cex.axis = 0.90,
      gap.axis = -1
    )

    for (pair_index in seq_along(comparisons)) {
      pair <- comparisons[[pair_index]]
      label <- if (two_groups || tolower(p_adjust_method) == "none") {
        paste0(pairwise_label, ", P=", format_p(adjusted_pair[pair_index]))
      } else {
        paste0(pairwise_label, ", ", p_adjust_method, " P=", format_p(adjusted_pair[pair_index]))
      }
      draw_bracket(match(pair[1], group_order), match(pair[2], group_order),
                   start + (pair_index - 1) * step, bracket_height, label)
      stats_rows[[stats_index]] <<- data.frame(
        Pathway = pathway, test_type = "pairwise", group1 = pair[1], group2 = pair[2],
        raw_pvalue = raw_pair[pair_index],
        adjusted_pvalue = if (two_groups) NA_real_ else adjusted_pair[pair_index],
        p_adjust_method = if (two_groups) "" else p_adjust_method
      )
      stats_index <<- stats_index + 1L
    }
    if (!two_groups) {
      stats_rows[[stats_index]] <<- data.frame(
        Pathway = pathway, test_type = "global", group1 = "", group2 = "",
        raw_pvalue = pathway_global_p, adjusted_pvalue = NA_real_, p_adjust_method = ""
      )
      stats_index <<- stats_index + 1L
    }

    limits <- par("usr")
    strip_label <- wrap_panel_label(pathway)
    available_width <- 0.94 * diff(limits[1:2])
    label_width <- max(strwidth(strsplit(strip_label, "\n", fixed = TRUE)[[1]], cex = 0.90, units = "user"))
    strip_cex <- if (is.finite(label_width) && label_width > 0) {
      max(0.72, min(0.90, 0.90 * available_width / label_width))
    } else {
      0.90
    }
    label_height <- strheight(strip_label, cex = strip_cex, units = "user")
    strip_height <- max(0.070 * diff(limits[3:4]), 1.22 * label_height)
    rect(limits[1], limits[4], limits[2], limits[4] + strip_height,
         col = "#E8E8E8", border = "#707070", xpd = NA)
    text(mean(limits[1:2]), limits[4] + strip_height / 2,
         strip_label, cex = strip_cex, font = 2, xpd = NA)
    if (!two_groups) {
      text(limits[1] + 0.02 * diff(limits[1:2]), limits[4] - 0.025 * diff(limits[3:4]),
           paste0(global_label, ", P=", format_p(pathway_global_p)),
           adj = c(0, 1), cex = 0.84)
    }
    ae_matrix_box(0.85, col = "#404040")
  }
  if (length(selected) < nrows * ncols) {
    for (unused in seq_len(nrows * ncols - length(selected))) plot.new()
  }
}

extension <- tolower(tools::file_ext(output_file))
if (extension == "png") {
  png(output_file, width = figure_width, height = figure_height, units = "in", res = plot_dpi, bg = "white")
} else if (extension == "pdf") {
  pdf(output_file, width = figure_width, height = figure_height, useDingbats = FALSE)
} else if (extension == "svg") {
  svg(output_file, width = figure_width, height = figure_height)
} else {
  stop("output extension must be png, pdf, or svg")
}
set.seed(123)
draw_plot()
dev.off()
write.table(do.call(rbind, stats_rows), statistics_file, sep = "\t", quote = FALSE, row.names = FALSE)
cat("Saved:", output_file, "and", statistics_file, "\n")
