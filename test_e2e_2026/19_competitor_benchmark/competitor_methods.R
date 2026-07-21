#!/usr/bin/env Rscript

# Minimal adapters for the two local R competitor implementations.
suppressPackageStartupMessages(library(data.table))

parse_args <- function(values) {
  result <- list()
  i <- 1
  while (i <= length(values)) {
    key <- sub("^--", "", values[[i]])
    if (i == length(values) || startsWith(values[[i + 1]], "--")) {
      result[[key]] <- TRUE
      i <- i + 1
    } else {
      result[[key]] <- values[[i + 1]]
      i <- i + 2
    }
  }
  result
}

read_gmt <- function(path) {
  lines <- readLines(path, warn = FALSE)
  fields <- strsplit(lines[nzchar(lines)], "\t", fixed = FALSE)
  fields <- fields[lengths(fields) >= 3]
  term2gene <- rbindlist(lapply(fields, function(row) {
    data.table(term_id = row[[1]], gene = unique(row[3:length(row)]))
  }))
  term2name <- unique(rbindlist(lapply(fields, function(row) {
    data.table(term_id = row[[1]], term_name = row[[2]])
  })))
  list(term2gene = term2gene, term2name = term2name)
}

read_genes <- function(path) {
  unique(trimws(readLines(path, warn = FALSE)))
}

read_ranked <- function(path) {
  ranking <- fread(path)
  if (ncol(ranking) < 2) stop("ranked input requires two columns")
  genes <- as.character(ranking[[1]])
  scores <- as.numeric(ranking[[2]])
  keep <- !duplicated(genes) & nzchar(genes) & is.finite(scores)
  scores <- scores[keep]
  names(scores) <- genes[keep]
  sort(scores, decreasing = TRUE)
}

empty_result <- function() {
  data.table(
    term_id = character(), term_name = character(), p_value = numeric(),
    adjusted_p_value = numeric(), es = numeric(), nes = numeric(),
    leading_edge = character(), overlap_count = integer(), term_size = integer()
  )
}

run_clusterprofiler <- function(method, opts, gmt) {
  suppressPackageStartupMessages(library(clusterProfiler))
  if (method == "ORA") {
    options(enrichment_force_universe = TRUE)
    result <- enricher(
      gene = read_genes(opts$query), universe = read_genes(opts$background),
      TERM2GENE = gmt$term2gene, TERM2NAME = gmt$term2name,
      pvalueCutoff = 1, qvalueCutoff = 1, pAdjustMethod = "BH",
      minGSSize = as.integer(opts$`min-size`), maxGSSize = Inf
    )
    frame <- as.data.table(as.data.frame(result))
    if (!nrow(frame)) return(empty_result())
    return(data.table(
      term_id = frame$ID, term_name = frame$Description,
      p_value = frame$pvalue, adjusted_p_value = frame$p.adjust,
      es = NA_real_, nes = NA_real_, leading_edge = "",
      overlap_count = frame$Count,
      term_size = as.integer(sub("/.*$", "", frame$BgRatio))
    ))
  }
  result <- GSEA(
    geneList = read_ranked(opts$ranked), TERM2GENE = gmt$term2gene,
    TERM2NAME = gmt$term2name, exponent = 1, minGSSize = as.integer(opts$`min-size`),
    maxGSSize = as.integer(opts$`max-size`), eps = 0, pvalueCutoff = 1,
    pAdjustMethod = "BH", by = "fgsea", seed = TRUE, verbose = FALSE
  )
  frame <- as.data.table(as.data.frame(result))
  if (!nrow(frame)) return(empty_result())
  data.table(
    term_id = frame$ID, term_name = frame$Description,
    p_value = frame$pvalue, adjusted_p_value = frame$p.adjust,
    es = frame$enrichmentScore, nes = frame$NES,
    leading_edge = gsub("/", ";", frame$core_enrichment, fixed = TRUE),
    overlap_count = NA_integer_, term_size = frame$setSize
  )
}

run_webgestalt <- function(method, opts) {
  suppressPackageStartupMessages(library(WebGestaltR))
  description_file <- paste0(opts$gmt, ".des")
  fwrite(gmt$term2name, description_file, sep = "\t", col.names = FALSE)
  common <- list(
    organism = "others", enrichDatabaseFile = opts$gmt,
    enrichDatabaseDescriptionFile = description_file,
    isOutput = FALSE, sigMethod = "top", topThr = 1000000,
    minNum = as.numeric(opts$`min-size`), fdrMethod = "BH",
    nThreads = 1
  )
  if (method == "ORA") {
    raw <- do.call(WebGestaltR, c(common, list(
      enrichMethod = "ORA", interestGene = read_genes(opts$query),
      referenceGene = read_genes(opts$background),
      maxNum = length(read_genes(opts$background))
    )))
  } else {
    ranking <- read_ranked(opts$ranked)
    raw <- do.call(WebGestaltR, c(common, list(
      enrichMethod = "GSEA",
      interestGene = data.frame(gene = names(ranking), score = unname(ranking)),
      maxNum = as.numeric(opts$`max-size`), perNum = as.numeric(opts$permutations),
      gseaP = 1, saveRawGseaResult = FALSE
    )))
  }
  frame <- as.data.table(raw)
  if (!nrow(frame)) return(empty_result())
  pick <- function(candidates, default = NA) {
    name <- candidates[candidates %in% names(frame)][1]
    if (is.na(name)) rep(default, nrow(frame)) else frame[[name]]
  }
  data.table(
    term_id = as.character(pick(c("geneSet", "geneSetID", "ID"), "")),
    term_name = as.character(pick(c("description", "geneSet", "name"), "")),
    p_value = as.numeric(pick(c("pValue", "pvalue", "p.value"), NA_real_)),
    adjusted_p_value = as.numeric(pick(c("FDR", "fdr", "p.adjust"), NA_real_)),
    es = as.numeric(pick(c("enrichmentScore", "ES"), NA_real_)),
    nes = as.numeric(pick(c("normalizedEnrichmentScore", "NES"), NA_real_)),
    leading_edge = vapply(pick(c("leadingEdgeId", "leadingEdge", "leading_edge", "userId"), ""), function(x) paste(x, collapse = ";"), character(1)),
    overlap_count = as.integer(pick(c("overlap", "overlap_count", "size"), NA_integer_)),
    term_size = as.integer(pick(c("geneSetSize", "setSize", "size"), NA_integer_))
  )
}

opts <- parse_args(commandArgs(trailingOnly = TRUE))
required <- c("tool", "method", "gmt", "output", "session-info", "min-size")
missing <- required[!required %in% names(opts)]
if (length(missing)) stop(paste("missing arguments:", paste(missing, collapse = ", ")))
set.seed(as.integer(ifelse(is.null(opts$seed), 42, opts$seed)))
gmt <- read_gmt(opts$gmt)
result <- if (opts$tool == "clusterProfiler") {
  run_clusterprofiler(opts$method, opts, gmt)
} else if (opts$tool == "WebGestaltR") {
  run_webgestalt(opts$method, opts)
} else {
  stop(paste("unsupported tool:", opts$tool))
}
result[, tool_version := as.character(packageVersion(opts$tool))]
fwrite(result, opts$output, sep = "\t", na = "")
writeLines(capture.output(sessionInfo()), opts$`session-info`)
