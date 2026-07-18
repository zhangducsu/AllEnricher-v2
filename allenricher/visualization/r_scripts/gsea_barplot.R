#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
})

this_file <- sub("^--file=", "", commandArgs()[4])
source(file.path(dirname(this_file), "common.R"))

args <- parse_args()
configure_plot_style(args)
tsv_path <- args$tsv
output <- args$output
top_n <- as.integer(args$top_n %||% 20)
plot_dpi <- as.integer(args$dpi %||% 300)

safe_neg_log10 <- function(x) {
  -log10(pmax(as.numeric(x), .Machine$double.xmin))
}

format_gene_rich <- function(gene_n, rich_factor) {
  paste0(
    ifelse(gene_n == round(gene_n), as.character(round(gene_n)), as.character(gene_n)),
    "/",
    sprintf("%.2f", as.numeric(rich_factor))
  )
}

split_hierarchy <- function(value) {
  if (is.na(value) || !nzchar(trimws(as.character(value)))) return(character(0))
  parts <- trimws(strsplit(as.character(value), "\\|", fixed = FALSE)[[1]])
  parts[nzchar(parts)]
}

select_hierarchy_level <- function(values, max_categories = 6) {
  paths <- lapply(values, split_hierarchy)
  paths <- paths[vapply(paths, length, integer(1)) >= 2]
  if (length(paths) == 0) return(list(level = NA_integer_, counts = integer(0)))
  max_level <- max(vapply(paths, length, integer(1)) - 1)
  counts <- vapply(seq_len(max_level), function(level) {
    length(unique(vapply(paths[vapply(paths, length, integer(1)) > level], `[[`, character(1), level)))
  }, integer(1))
  selected <- which(counts >= 2 & counts <= max_categories)[1]
  list(level = ifelse(length(selected), selected, NA_integer_), counts = counts)
}

category_at_level <- function(value, level) {
  parts <- split_hierarchy(value)
  if (is.na(level) || length(parts) <= level) return(NA_character_)
  parts[level]
}

prepare_ora_bar_data <- function(df, top_n = 20) {
  database <- if ("Database" %in% colnames(df) && nrow(df) > 0) df$Database[1] else ""
  hierarchy_col <- intersect(c("Hierarchy", "hierarchy"), colnames(df))[1]
  hierarchy_values <- if (!is.na(hierarchy_col)) {
    as.character(df[[hierarchy_col]])
  } else {
    ifelse(grepl("|", df$Term_Name, fixed = TRUE), as.character(df$Term_Name), "")
  }
  hierarchy <- select_hierarchy_level(hierarchy_values)
  if (length(hierarchy$counts)) {
    message("ORA hierarchy category counts: ", paste0("level ", seq_along(hierarchy$counts), "=", hierarchy$counts, collapse = ", "))
  }
  categories <- if (!is.na(hierarchy$level)) {
    vapply(hierarchy_values, category_at_level, character(1), level = hierarchy$level)
  } else {
    rep(NA_character_, length(hierarchy_values))
  }
  labels <- vapply(df$Term_Name, function(term_name) {
    parts <- split_hierarchy(term_name)
    if (length(parts) >= 2) tail(parts, 1) else unname(clean_pathway_name(term_name))
  }, character(1))

  rich_factor <- if ("EnrichFactor" %in% colnames(df)) df$EnrichFactor else df$Rich_Factor
  plotdf <- data.frame(
    database = database,
    category = categories,
    label = labels,
    hierarchy_level = hierarchy$level,
    has_category = !is.na(hierarchy$level) && all(!is.na(categories)),
    gene_n = as.numeric(df$Gene_Count),
    rich_factor = as.numeric(rich_factor),
    qvalue = as.numeric(df$Adjusted_P_Value),
    stringsAsFactors = FALSE
  )

  plotdf <- plotdf %>%
    filter(is.finite(gene_n), is.finite(rich_factor), is.finite(qvalue), qvalue > 0) %>%
    mutate(
      score = safe_neg_log10(qvalue),
      gene_rich_label = format_gene_rich(gene_n, rich_factor)
    ) %>%
    arrange(qvalue, desc(score))

  if (nrow(plotdf) == 0) {
    stop("No valid ORA rows available for barplot")
  }
  head(plotdf, top_n)
}

draw_compact_barplot <- function(plotdf, output_file) {
  database <- toupper(as.character(plotdf$database[1]))
  use_category <- any(plotdf$has_category) && all(!is.na(plotdf$category))

  row_gap <- 0.62
  plotdf <- plotdf %>% arrange(score)
  plotdf$y <- seq_len(nrow(plotdf)) * row_gap
  plotdf$label_wrap <- wrap_label(plotdf$label, width = ifelse(database == "KEGG", 44, 46))
  xmax <- max(plotdf$score, na.rm = TRUE)
  if (!is.finite(xmax) || xmax <= 0) {
    xmax <- 1
  }
  bar_h <- row_gap * 0.62
  bar_layer <- if (use_category) {
    geom_rect(aes(xmin = 0, xmax = score, ymin = y - bar_h / 2, ymax = y + bar_h / 2, fill = category), color = "grey30", linewidth = ae_lwd(0.35))
  } else {
    geom_rect(aes(xmin = 0, xmax = score, ymin = y - bar_h / 2, ymax = y + bar_h / 2, fill = score), color = "grey30", linewidth = ae_lwd(0.35))
  }

  p <- ggplot(plotdf, aes(y = y)) +
    geom_vline(xintercept = 0, color = "grey20", linewidth = ae_lwd(1.05)) +
    bar_layer +
    geom_text(aes(x = score + xmax * 0.020, label = gene_rich_label), hjust = 1, size = 2.95) +
    geom_text(aes(x = -xmax * 0.020, label = label_wrap), hjust = 0, size = 3.05, fontface = "italic", lineheight = 0.84) +
    scale_x_reverse(
      name = "-log10(Q-value)",
      limits = c(xmax * 1.18, -xmax * 0.78),
      breaks = pretty(c(0, xmax), n = 5),
      expand = expansion(mult = 0)
    ) +
    scale_y_continuous(NULL, breaks = NULL, limits = c(min(plotdf$y) - 0.40, max(plotdf$y) + 0.50), expand = expansion(mult = 0)) +
    coord_cartesian(clip = "off") +
    labs(title = ifelse(database == "KEGG", "KEGG Pathways (Gene# / Rich Factor)", paste0(database, " (Gene# / Rich Factor)"))) +
    set_nature_theme(base_size = 9) +
    theme(
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.line.y = element_blank(),
      axis.title.x = element_text(size = 10.4 * AE_TEXT_SCALE),
      axis.text.x = element_text(size = 9.3 * AE_TEXT_SCALE),
      plot.title = element_text(size = 11.4 * AE_TEXT_SCALE, hjust = 0.5, face = "plain", margin = margin(b = 4)),
      plot.margin = ae_margin(4, ifelse(database == "KEGG", 76, 68), 4, 6),
      legend.position = "bottom",
      legend.key.size = grid::unit(4.2, "mm"),
      legend.text = element_text(size = ifelse(database == "KEGG", 7.4, 8.0) * AE_TEXT_SCALE)
    )

  if (use_category) {
    categories <- sort(unique(plotdf$category))
    category_colors <- setNames(ae_palette(length(categories)), categories)
    p <- p +
      scale_fill_manual(values = category_colors, breaks = categories, name = paste("Hierarchy level", plotdf$hierarchy_level[1])) +
      guides(fill = guide_legend(ncol = ifelse(length(categories) > 3, 2, 1)))
  } else {
    p <- p +
      scale_fill_gradientn(
        colors = ae_sequential_colors(6),
        name = "-log10(Q-value)",
        breaks = pretty(range(plotdf$score, na.rm = TRUE), n = 3),
        guide = guide_colorbar(
          title.position = "top", title.hjust = 0.5,
          barwidth = grid::unit(44, "mm"), barheight = grid::unit(4, "mm")
        )
      )
  }

  save_plot(p, output_file, width = 8.0, height = max(4.0, 0.25 * nrow(plotdf) + 1.15), dpi = plot_dpi)
}

df <- read_enrichment(tsv_path)
has_nes <- "NES" %in% colnames(df) && any(is.finite(df$NES) & df$NES != 0)

if (has_nes) {
  database_label <- if ("Database" %in% colnames(df) && nrow(df) > 0) {
    trimws(as.character(df$Database[1]))
  } else {
    ""
  }
  chart_title <- trimws(paste(database_label, "GSEA NES Ranking"))
  plot_gsea_diverging_bar(df, output, top_n = top_n, title = chart_title, dpi = plot_dpi)
} else {
  plotdf <- prepare_ora_bar_data(df, top_n = top_n)
  draw_compact_barplot(plotdf, output)
}
