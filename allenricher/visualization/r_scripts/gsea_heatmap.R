#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(ComplexHeatmap)
  library(circlize)
})
this_file <- sub("^--file=", "", commandArgs()[4])
source(file.path(dirname(this_file), "common.R"))
ht_opt$message <- FALSE

args <- parse_args()
expr_path <- args$expr
tsv_path <- args$tsv %||% ""
output <- args$output
top_n <- as.integer(args$top_n %||% 12)

expr <- read.delim(expr_path, row.names = 1, sep = "\t", check.names = FALSE)
expr_all <- as.matrix(expr)
mode(expr_all) <- "numeric"

select_variable_genes <- function(mat, max_genes = 80) {
  row_var <- apply(mat, 1, var, na.rm = TRUE)
  mat <- mat[is.finite(row_var) & row_var > 0, , drop = FALSE]
  if (nrow(mat) > max_genes) {
    keep <- order(row_var[rownames(mat)], decreasing = TRUE)[seq_len(max_genes)]
    mat <- mat[keep, , drop = FALSE]
  }
  mat
}

selected_genes <- character(0)
heatmap_title <- "Most variable genes"

if (!is.null(tsv_path) && tsv_path != "" && file.exists(tsv_path)) {
  pathway_df <- read_enrichment(tsv_path) %>%
    prepare_gsea_terms(top_n = top_n, sort_by = "abs_nes", label_width = 34)
  gene_lists <- parse_gene_list(pathway_df$Genes)
  selected_genes <- unique(unlist(gene_lists))
  selected_genes <- intersect(selected_genes, rownames(expr_all))
  if (length(selected_genes) > 0) {
    heatmap_title <- paste0("Leading-edge genes from top ", nrow(pathway_df), " GSEA pathways")
  }
}

if (length(selected_genes) > 0) {
  expr_mat <- expr_all[selected_genes, , drop = FALSE]
  expr_mat <- select_variable_genes(expr_mat, max_genes = 80)
  if (nrow(expr_mat) == 0) {
    warning("Selected pathway genes have no usable variance; falling back to most variable genes")
    expr_mat <- select_variable_genes(expr_all, max_genes = 80)
    heatmap_title <- "Most variable genes"
  }
} else {
  expr_mat <- select_variable_genes(expr_all, max_genes = 80)
}

if (nrow(expr_mat) == 0) {
  stop("No variable genes available for heatmap")
}

expr_z <- t(scale(t(expr_mat)))
expr_z[!is.finite(expr_z)] <- 0
expr_z <- pmax(pmin(expr_z, 2.5), -2.5)

ht <- Heatmap(
  expr_z,
  name = "Z-score",
  col = colorRamp2(c(-2.5, 0, 2.5), c(AE_COL_DOWN, "white", AE_COL_UP)),
  show_row_names = nrow(expr_z) <= 60,
  row_names_gp = grid::gpar(fontsize = 5.8),
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  column_names_gp = grid::gpar(fontsize = 8),
  column_title = heatmap_title,
  column_title_gp = grid::gpar(fontsize = 10, fontface = "bold"),
  heatmap_legend_param = list(title_gp = grid::gpar(fontsize = 8, fontface = "bold"), labels_gp = grid::gpar(fontsize = 7)),
  use_raster = nrow(expr_z) > 50
)

ext <- tolower(tools::file_ext(output))
if (ext == "png") {
  png(output, width = 7.2, height = 6.4, units = "in", res = 300, bg = "white")
} else if (ext == "pdf") {
  pdf(output, width = 7.2, height = 6.4, bg = "white")
} else if (ext == "svg") {
  svg(output, width = 7.2, height = 6.4, bg = "white")
} else {
  stop(paste("Unsupported output format:", ext))
}
draw(ht, heatmap_legend_side = "right")
dev.off()
message(paste("Saved:", output))
