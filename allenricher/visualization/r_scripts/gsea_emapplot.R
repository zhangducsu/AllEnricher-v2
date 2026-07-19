#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(ggplot2))

if (!requireNamespace("aPEAR", quietly = TRUE)) {
  stop("R package 'aPEAR' is required for emapplot. Install with: remotes::install_github('ievaKer/aPEAR')")
}
if (!requireNamespace("scales", quietly = TRUE)) {
  stop("R package 'scales' is required for emapplot")
}

this_file <- sub("^--file=", "", commandArgs()[4])
source(file.path(dirname(this_file), "common.R"))

args <- parse_args()
configure_plot_style(args)
tsv_path <- args$tsv
output <- args$output
top_n <- as.integer(args$top_n %||% 30)
qvalue_cutoff <- as.numeric(args$qvalue %||% 0.05)
min_count <- as.integer(args$min_count %||% 3)
plot_dpi <- as.integer(args$dpi %||% 300)

if (is.na(top_n) || top_n < 1) stop("top_n must be a positive integer")
if (is.na(qvalue_cutoff) || qvalue_cutoff < 0 || qvalue_cutoff > 1) {
  stop("qvalue must be between 0 and 1")
}
if (is.na(min_count) || min_count < 0) stop("min_count must be a non-negative integer")

save_empty <- function(message) {
  warning(message, call. = FALSE)
  p <- ggplot() +
    labs(title = "Pathway network unavailable", subtitle = message) +
    theme_void(base_size = 11) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
      plot.subtitle = element_text(hjust = 0.5, size = 11),
      plot.margin = margin(24, 16, 24, 16)
    )
  save_plot(p, output, width = 5.4, height = 2.0, dpi = plot_dpi)
  quit(save = "no", status = 0)
}

balanced_top_terms <- function(df, limit) {
  order_terms <- function(x) {
    x[order(x$Adjusted_P_Value, -abs(x$NES)), , drop = FALSE]
  }
  if (nrow(df) <= limit) return(order_terms(df))

  up <- order_terms(df[df$NES > 0, , drop = FALSE])
  down <- order_terms(df[df$NES < 0, , drop = FALSE])
  selected <- rbind(head(up, ceiling(limit / 2)), head(down, floor(limit / 2)))

  if (nrow(selected) < limit) {
    remainder <- order_terms(df[!df$Term_ID %in% selected$Term_ID, , drop = FALSE])
    selected <- rbind(selected, head(remainder, limit - nrow(selected)))
  }
  selected
}

df <- read_enrichment(tsv_path)
df$NES <- suppressWarnings(as.numeric(df$NES))
df$Adjusted_P_Value <- suppressWarnings(as.numeric(df$Adjusted_P_Value))
df$Gene_Count <- suppressWarnings(as.numeric(df$Gene_Count))
df <- df[
  is.finite(df$NES) &
    is.finite(df$Adjusted_P_Value) &
    df$Adjusted_P_Value > 0 &
    df$Adjusted_P_Value < qvalue_cutoff &
    is.finite(df$Gene_Count) &
    df$Gene_Count >= min_count,
  ,
  drop = FALSE
]

if (nrow(df) < 2) {
  save_empty(sprintf(
    "Fewer than 2 pathways after filtering (FDR < %s, setSize >= %s)",
    qvalue_cutoff,
    min_count
  ))
}

df <- balanced_top_terms(df, top_n)
gene_lists <- parse_gene_list(df$Genes)
has_overlap <- any(vapply(seq_len(nrow(df) - 1), function(i) {
  any(vapply(seq.int(i + 1, nrow(df)), function(j) {
    length(intersect(gene_lists[[i]], gene_lists[[j]])) > 0
  }, logical(1)))
}, logical(1)))
if (!has_overlap) save_empty("No leading-gene overlap among selected pathways")

labels <- clean_pathway_name(df$Term_Name)
network_input <- data.frame(
  ID = as.character(df$Term_ID),
  Description = make_unique_labels(labels, df$Term_ID),
  pathwayGenes = vapply(gene_lists, function(x) paste(x, collapse = "/"), character(1)),
  NES = df$NES,
  setSize = df$Gene_Count,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

set.seed(123)
p <- tryCatch(
  withCallingHandlers(
    aPEAR::enrichmentNetwork(
      enrichment = network_input,
      simMethod = "jaccard",
      clustMethod = "markov",
      clustNameMethod = "pagerank",
      colorBy = "NES",
      nodeSize = "setSize",
      colorType = "nes",
      innerCutoff = 0.10,
      outerCutoff = 0.65,
      minClusterSize = 2,
      fontSize = 4.0 * AE_TEXT_SCALE,
      repelLabels = TRUE,
      drawEllipses = TRUE,
      verbose = FALSE
    ),
    warning = function(w) {
      if (grepl("Ignoring unknown aesthetics", conditionMessage(w), fixed = TRUE)) {
        invokeRestart("muffleWarning")
      }
    }
  ),
  error = function(e) save_empty(paste("aPEAR could not build a network:", conditionMessage(e)))
)

nes_limit <- max(abs(network_input$NES), na.rm = TRUE)
if (!is.finite(nes_limit) || nes_limit <= 0) nes_limit <- 1
size_breaks <- pretty(range(network_input$setSize), n = 4)
size_breaks <- size_breaks[
  size_breaks >= min(network_input$setSize) & size_breaks <= max(network_input$setSize)
]

# APEAR with default/sizeScale; remove first and avoid repeating scalewarning.
p$scales$scales <- Filter(function(scale) {
  !any(scale$aesthetics %in% c("colour", "color", "size"))
}, p$scales$scales)
p <- p +
  scale_color_gradientn(
    name = "NES",
    colours = ae_diverging_colors(6),
    values = scales::rescale(seq(-nes_limit, nes_limit, length.out = 6)),
    limits = c(-nes_limit, nes_limit),
    oob = scales::squish
  ) +
  scale_size_continuous(
    name = "Pathway size",
    range = c(3.2, 10.5) * sqrt(AE_LINE_SCALE),
    breaks = size_breaks
  ) +
  theme_void(base_size = 13 * AE_TEXT_SCALE, base_family = AE_FONT_FAMILY) +
  theme(
    text = element_text(family = AE_FONT_FAMILY, colour = AE_COL_TEXT),
    legend.position = "right",
    legend.box = "vertical",
    legend.title = element_text(size = 13 * AE_TEXT_SCALE),
    legend.text = element_text(size = 11.5 * AE_TEXT_SCALE),
    plot.margin = ae_margin(4, 6, 4, 4)
  )

save_plot(p, output, width = 10.2, height = 7.4, dpi = plot_dpi)
