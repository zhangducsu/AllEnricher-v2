#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(scales)
  library(grid)
})

this_file <- sub("^--file=", "", commandArgs()[4])
source(file.path(dirname(this_file), "common.R"))

args <- parse_args()
configure_plot_style(args)
tsv_path <- args$tsv
output <- args$output
top_n <- as.integer(args$top_n %||% 20)
plot_dpi <- as.integer(args$dpi %||% 300)

format_fdr_label <- function(score) {
  vapply(score, function(value) {
    parts <- strsplit(formatC(10^(-value), format = "e", digits = 1), "e", fixed = TRUE)[[1]]
    paste0(parts[1], "e", as.integer(parts[2]))
  }, character(1))
}

build_fdr_breaks <- function(scores) {
  rng <- range(scores, na.rm = TRUE)
  if (!all(is.finite(rng))) {
    return(c(1, 2, 3))
  }
  if (rng[1] == rng[2]) {
    rng <- rng + c(-0.5, 0.5)
  }
  seq(rng[1], rng[2], length.out = 5)
}

draw_string_style_lollipop <- function(df, metric_col, metric_label, title, is_nes = FALSE) {
  df <- df %>%
    mutate(
      metric = .data[[metric_col]],
      fdr = Adjusted_P_Value,
      fdr_score = -log10(pmax(fdr, .Machine$double.xmin)),
      term_wrap = wrap_label(Term_Short, width = 34)
    ) %>%
    filter(is.finite(metric), is.finite(fdr_score), is.finite(Gene_Count)) %>%
    arrange(metric)

  if (nrow(df) == 0) {
    stop("No valid rows available for lollipop plot")
  }

  df$y <- seq_len(nrow(df))

  is_bidirectional <- is_nes && any(df$metric > 0, na.rm = TRUE) && any(df$metric < 0, na.rm = TRUE)

  if (is_bidirectional) {
    xmax <- max(abs(df$metric), na.rm = TRUE)
    if (!is.finite(xmax) || xmax <= 0) {
      xmax <- 1
    }
    x0 <- 0
    xmin_main <- -xmax * 1.08
    xmax_main <- xmax * 1.08
    x_limits <- c(-xmax * 1.18, xmax * 1.18)
    x_breaks <- pretty(x_limits, n = 7)
  } else {
    xmin <- min(df$metric, na.rm = TRUE)
    xmax <- max(df$metric, na.rm = TRUE)
    span <- max(xmax - xmin, 0.1)
    if (is_nes && xmax <= 0) {
      x0 <- xmax + span * 0.06
      xmin_main <- xmin - span * 0.02
      xmax_main <- x0
      x_limits <- c(xmin - span * 0.12, x0 + span * 0.02)
    } else {
      x0 <- xmin - span * 0.06
      xmin_main <- x0
      xmax_main <- xmax + span * 0.02
      x_limits <- c(x0 - span * 0.02, xmax + span * 0.12)
    }
    x_breaks <- pretty(x_limits, n = 7)
  }

  fdr_cols <- ae_gradient_colors(6)
  fdr_breaks <- build_fdr_breaks(df$fdr_score)
  fdr_limits <- range(fdr_breaks)
  gene_breaks <- pretty(df$Gene_Count, n = 3)
  gene_breaks <- gene_breaks[gene_breaks > 0]

  p <- ggplot(df) +
    geom_rect(
      aes(
        xmin = xmin_main,
        xmax = xmax_main,
        ymin = y - 0.45,
        ymax = y + 0.45,
        fill = fdr_score
      ),
      alpha = 0.12,
      colour = NA
    ) +
    geom_segment(
      aes(
        x = x0,
        xend = metric,
        y = y,
        yend = y,
        colour = fdr_score
      ),
      linewidth = ae_lwd(1.2),
      lineend = "butt"
    ) +
    geom_segment(
      aes(
        x = x0,
        xend = metric,
        y = y,
        yend = y
      ),
      colour = "white",
      linewidth = ae_lwd(1.2),
      alpha = 0.18,
      lineend = "butt"
    ) +
    geom_segment(
      aes(
        x = metric,
        xend = metric,
        y = y - 0.26,
        yend = y + 0.26
      ),
      colour = "grey55",
      linewidth = ae_lwd(0.25),
      linetype = "22",
      alpha = 0.8
    ) +
    geom_point(
      aes(
        x = metric,
        y = y,
        size = Gene_Count,
        fill = fdr_score
      ),
      shape = 21,
      colour = "grey55",
      stroke = ae_lwd(1.1),
      alpha = 1
    ) +
    scale_fill_gradientn(
      name = "FDR",
      colours = fdr_cols,
      values = scales::rescale(fdr_breaks),
      limits = fdr_limits,
      breaks = fdr_breaks,
      labels = format_fdr_label(fdr_breaks),
      guide = guide_colourbar(
        order = 1,
        barheight = unit(48, "pt"),
        barwidth = unit(8, "pt"),
        title.position = "top"
      )
    ) +
    scale_colour_gradientn(
      colours = fdr_cols,
      values = scales::rescale(fdr_breaks),
      limits = fdr_limits,
      guide = "none"
    ) +
    scale_size_area(
      name = "Gene count",
      max_size = 8.8,
      breaks = gene_breaks,
      guide = guide_legend(
        order = 2,
        override.aes = list(
          fill = tail(fdr_cols, 1),
          colour = AE_COL_NEUTRAL,
          alpha = 1
        )
      )
    ) +
    scale_x_continuous(
      name = metric_label,
      limits = x_limits,
      breaks = x_breaks,
      expand = expansion(mult = 0)
    ) +
    scale_y_continuous(
      name = NULL,
      breaks = df$y,
      labels = df$term_wrap,
      limits = c(0.5, nrow(df) + 0.5),
      expand = expansion(mult = 0)
    ) +
    coord_cartesian(clip = "off") +
    labs(title = title) +
    set_nature_theme(base_size = 9) +
    theme(
      panel.grid.major.x = if (AE_SHOW_GRID) element_line(
        colour = AE_COL_GRID,
        linewidth = ae_lwd(0.35),
        linetype = "dashed"
      ) else element_blank(),
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(size = 8.5 * AE_TEXT_SCALE, colour = "black"),
      axis.text.y = element_text(
        size = 8.5 * AE_TEXT_SCALE,
        colour = "black",
        hjust = 1,
        margin = margin(r = 4)
      ),
      axis.title.x = element_text(
        size = 9 * AE_TEXT_SCALE,
        colour = "black",
        margin = margin(t = 8)
      ),
      plot.title = element_text(
        size = 11 * AE_TEXT_SCALE,
        hjust = 0.5,
        face = "plain",
        margin = margin(b = 6)
      ),
      plot.margin = ae_margin(8, 70, 8, 6),
      legend.position = "right",
      legend.box = "vertical",
      legend.title = element_text(size = 9 * AE_TEXT_SCALE, colour = "black"),
      legend.text = element_text(size = 8.5 * AE_TEXT_SCALE, colour = "black"),
      legend.background = element_rect(fill = "white", colour = NA),
      legend.key = element_rect(fill = "white", colour = NA)
    )

  if (is_bidirectional) {
    p <- p + geom_vline(xintercept = 0, linewidth = ae_lwd(0.35), colour = "grey55")
  }

  save_plot(p, output, width = 9.2, height = max(5.1, 1.5 + 0.30 * nrow(df)), dpi = plot_dpi)
}

raw_df <- read_enrichment(tsv_path)
database_label <- if ("Database" %in% colnames(raw_df) && nrow(raw_df) > 0) {
  trimws(as.character(raw_df$Database[1]))
} else {
  ""
}
df <- raw_df %>%
  prepare_gsea_terms(top_n = top_n, sort_by = "q", label_width = 34)

has_nes <- "NES" %in% colnames(df) && any(is.finite(df$NES) & df$NES != 0)
has_enrich_factor <- "EnrichFactor" %in% colnames(df) && any(is.finite(df$EnrichFactor) & df$EnrichFactor != 0)

if (has_nes) {
  plot_df <- raw_df %>%
    prepare_gsea_diverging_terms(top_n = top_n, label_width = 34)
  draw_string_style_lollipop(
    plot_df,
    metric_col = "NES",
    metric_label = "NES",
    title = trimws(paste(database_label, "GSEA Lollipop Plot")),
    is_nes = TRUE
  )
} else if (has_enrich_factor) {
  draw_string_style_lollipop(
    df,
    metric_col = "EnrichFactor",
    metric_label = "EnrichFactor",
    title = trimws(paste(database_label, "Enrichment Lollipop Plot")),
    is_nes = FALSE
  )
} else {
  draw_string_style_lollipop(
    df,
    metric_col = "NES",
    metric_label = "NES",
    title = trimws(paste(database_label, "GSEA Lollipop Plot")),
    is_nes = TRUE
  )
}
