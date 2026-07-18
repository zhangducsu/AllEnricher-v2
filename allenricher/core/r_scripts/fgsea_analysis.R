args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 6) {
  stop("Usage: fgsea_analysis.R ranking.tsv gene_sets.gmt output.tsv min_size max_size seed")
}
if (!requireNamespace("fgsea", quietly = TRUE)) {
  stop("The Bioconductor package 'fgsea' is required but not installed")
}

read_gmt <- function(path) {
  records <- strsplit(readLines(path, warn = FALSE, encoding = "UTF-8"), "\t", fixed = TRUE)
  pathways <- lapply(records, function(record) unique(record[-c(1, 2)]))
  names(pathways) <- vapply(records, function(record) record[[1]], character(1))
  pathways
}

ranking <- read.delim(args[[1]], check.names = FALSE, stringsAsFactors = FALSE)
if (!all(c("gene", "weight") %in% colnames(ranking))) {
  stop("The ranking table must contain 'gene' and 'weight' columns")
}
if (anyDuplicated(ranking$gene)) {
  stop("The ranked gene list contains duplicate gene identifiers")
}

stats <- as.numeric(ranking$weight)
names(stats) <- as.character(ranking$gene)
if (any(!is.finite(stats))) {
  stop("Ranking weights must be finite numeric values")
}
stats <- sort(stats, decreasing = TRUE)
pathways <- read_gmt(args[[2]])
set.seed(as.integer(args[[6]]))

fgsea_result <- fgsea::fgseaMultilevel(
  pathways = pathways,
  stats = stats,
  minSize = as.integer(args[[4]]),
  maxSize = as.integer(args[[5]]),
  eps = 0,
  nproc = 1
)

result <- data.frame(
  pathway = as.character(fgsea_result$pathway),
  pval = as.numeric(fgsea_result$pval),
  padj = as.numeric(fgsea_result$padj),
  log2err = as.numeric(fgsea_result$log2err),
  ES = as.numeric(fgsea_result$ES),
  NES = as.numeric(fgsea_result$NES),
  size = as.integer(fgsea_result$size),
  leadingEdge = vapply(fgsea_result$leadingEdge, paste, collapse = ";", character(1)),
  check.names = FALSE,
  stringsAsFactors = FALSE
)
write.table(result, args[[3]], sep = "\t", row.names = FALSE, quote = FALSE, na = "NA", fileEncoding = "UTF-8")
