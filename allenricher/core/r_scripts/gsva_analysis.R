args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 8) {
  stop("Usage: gsva_analysis.R expression.tsv gene_sets.gmt output.tsv method kcdf tau min_size max_size")
}
if (!requireNamespace("GSVA", quietly = TRUE)) {
  stop("No Bioconductor package GSVA available. Please install GSVA")
}

read_gmt <- function(path) {
  records <- strsplit(readLines(path, warn = FALSE, encoding = "UTF-8"), "\t", fixed = TRUE)
  pathways <- lapply(records, function(record) unique(record[-c(1, 2)]))
  names(pathways) <- vapply(records, function(record) record[[1]], character(1))
  pathways
}

expression <- read.delim(args[[1]], check.names = FALSE, stringsAsFactors = FALSE)
if (ncol(expression) < 2 || colnames(expression)[[1]] != "gene") {
  stop("Expression.tsv first column must be gene, and subsequent column must be sampled")
}
if (anyDuplicated(expression$gene)) {
  stop("The expression matrix contains duplicated gene IDs")
}

genes <- as.character(expression$gene)
expr <- as.matrix(expression[, -1, drop = FALSE])
storage.mode(expr) <- "double"
rownames(expr) <- genes
if (any(!is.finite(expr))) {
  stop("Expression matrix contains NaN or infinite value")
}

pathways <- read_gmt(args[[2]])
method <- args[[4]]
min_size <- as.integer(args[[7]])
max_size <- as.integer(args[[8]])

if (method == "ssgsea") {
  parameter <- GSVA::ssgseaParam(
    exprData = expr,
    geneSets = pathways,
    alpha = as.numeric(args[[6]]),
    normalize = TRUE,
    minSize = min_size,
    maxSize = max_size,
    verbose = FALSE
  )
} else if (method == "gsva") {
  parameter <- GSVA::gsvaParam(
    exprData = expr,
    geneSets = pathways,
    kcdf = args[[5]],
    tau = as.numeric(args[[6]]),
    minSize = min_size,
    maxSize = max_size,
    filterRows = ncol(expr) > 1,
    verbose = FALSE
  )
} else if (method == "plage") {
  parameter <- GSVA::plageParam(
    exprData = expr,
    geneSets = pathways,
    minSize = min_size,
    maxSize = max_size,
    verbose = FALSE
  )
} else if (method == "zscore") {
  parameter <- GSVA::zscoreParam(
    exprData = expr,
    geneSets = pathways,
    minSize = min_size,
    maxSize = max_size,
    verbose = FALSE
  )
} else {
  stop(paste("Unsupported GSVA method:", method))
}

scores <- GSVA::gsva(parameter, verbose = FALSE)
result <- data.frame(Pathway = rownames(scores), scores, check.names = FALSE, stringsAsFactors = FALSE)
write.table(result, args[[3]], sep = "\t", row.names = FALSE, quote = FALSE, na = "NA", fileEncoding = "UTF-8")
