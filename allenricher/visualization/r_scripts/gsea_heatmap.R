#!/usr/bin/env Rscript
# gsea_heatmap.R — 通路基因热图
suppressPackageStartupMessages({
  library(ComplexHeatmap)
  library(circlize)
  library(dplyr)
})
this_file <- sub("^--file=", "", commandArgs()[4])
source(file.path(dirname(this_file), "common.R"))

args <- parse_args()
expr_path <- args$expr
output <- args$output

expr <- read.delim(expr_path, row.names = 1, sep = "\t", check.names = FALSE)

# Z-score 标准化
expr_z <- t(scale(t(as.matrix(expr))))

# 过滤方差为0的行
expr_z <- expr_z[apply(expr_z, 1, var) > 0, , drop = FALSE]

ht <- Heatmap(expr_z,
              name = "Z-score",
              col = colorRamp2(c(-2, 0, 2), c("#377EBA", "white", "#E41A1C")),
              show_row_names = FALSE,
              cluster_rows = TRUE,
              cluster_columns = TRUE,
              column_names_gp = gpar(fontsize = 8),
              column_title = "Gene Expression Heatmap")

png(output, width = 10, height = 8, units = "in", res = 300)
draw(ht)
dev.off()
message(paste("Saved:", output))
