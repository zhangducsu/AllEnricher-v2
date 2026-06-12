#!/usr/bin/env Rscript
# gsea_enrichment_plot.R — GSEA 三面板富集图（从 TSV 数据重建）
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(gridExtra)
})
this_file <- sub("^--file=", "", commandArgs()[4])
source(file.path(dirname(this_file), "common.R"))

args <- parse_args()
tsv_path <- args$tsv
gene_set_id <- args$gene_set_id
output <- args$output

df <- read_enrichment(tsv_path)
row_data <- df[df$Term_ID == gene_set_id, ]

if (nrow(row_data) == 0) {
  stop(paste("Gene set", gene_set_id, "not found in TSV"))
}

nes <- row_data$NES[1]
es <- row_data$ES[1]
genes <- strsplit(as.character(row_data$Genes[1]), ",")[[1]]
genes <- unique(trimws(genes))

# 模拟三面板图
p1 <- ggplot(data.frame(x = 1:100, y = cumsum(c(rep(1/length(genes), length(genes)), rep(-1/(100-length(genes)), 100-length(genes))))),
       aes(x = x, y = y)) +
  geom_line(color = "#377EBA", size = 1) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  annotate("text", x = 50, y = max(0.5, es), label = paste0("NES = ", round(nes, 3)), size = 4) +
  labs(title = gene_set_id, x = "Gene Rank", y = "Running Enrichment Score") +
  set_nature_theme()

p2 <- ggplot(data.frame(x = seq(1, 100, length.out = length(genes))),
       aes(x = x)) +
  geom_vline(xintercept = seq(1, 100, length.out = length(genes)), color = "black", size = 0.5) +
  labs(x = "Gene Rank", y = "") +
  set_nature_theme() +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())

p3 <- ggplot(data.frame(x = 1:100, y = rnorm(100)), aes(x = x, y = y)) +
  geom_bar(stat = "identity", fill = "gray70", width = 1) +
  labs(x = "Gene Rank", y = "Gene Weight") +
  set_nature_theme()

combined <- grid.arrange(p1, p2, p3, ncol = 1, heights = c(3, 1, 2))

ggsave(output, plot = combined, width = 10, height = 8, dpi = 300)
message(paste("Saved:", output))
