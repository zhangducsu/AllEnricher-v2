#!/usr/bin/env Rscript
# gsea_cnetplot.R — 通路-基因关联网络图 (简化版: 基于 ggplot2 的基因-通路矩阵)
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
})
this_file <- sub("^--file=", "", commandArgs()[4])
source(file.path(dirname(this_file), "common.R"))

args <- parse_args()
tsv_path <- args$tsv
output <- args$output
top_n <- as.integer(args$top_n %||% 10)

df <- read_enrichment(tsv_path) %>%
  arrange(desc(abs(NES))) %>%
  head(top_n)

# 解析基因集 (用分号分隔)
gene_sets <- list()
for (i in 1:nrow(df)) {
  genes <- strsplit(as.character(df$Genes[i]), ";")[[1]]
  genes <- unique(trimws(genes))
  gene_sets[[df$Term_Name[i]]] <- genes
}

# 构建基因-通路关联矩阵
all_genes <- unique(unlist(gene_sets))
pathway_names <- names(gene_sets)

mat <- data.frame(matrix(0, nrow = length(all_genes), ncol = length(pathway_names),
                         dimnames = list(all_genes, pathway_names)), check.names = FALSE)
for (j in seq_along(pathway_names)) {
  mat[gene_sets[[j]], j] <- 1
}

# 统计每个基因出现的通路数
gene_freq <- rowSums(mat)
# 只保留出现在 >= 2 个通路的基因 (连接基因)
conn_genes <- names(gene_freq[gene_freq >= 2])
if (length(conn_genes) == 0) {
  # fallback: 取 top 基因
  conn_genes <- names(sort(gene_freq, decreasing = TRUE))[1:min(30, length(gene_freq))]
}

mat_sub <- mat[conn_genes, , drop = FALSE]
mat_sub$Gene <- rownames(mat_sub)

# 使用 tidyr::pivot_longer 代替 reshape2::melt
melted <- mat_sub %>%
  pivot_longer(cols = -Gene, names_to = "Pathway", values_to = "Present") %>%
  filter(Present == 1)

# 绘制基因-通路关联图
p <- ggplot(melted, aes(x = Pathway, y = Gene, fill = Pathway)) +
  geom_tile(color = "white", linewidth = 0.5) +
  scale_fill_brewer(palette = "Set3", guide = "none") +
  labs(x = "Pathway", y = "Gene",
       title = "Gene-Pathway Network") +
  set_nature_theme() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
        axis.text.y = element_text(size = 6),
        plot.title = element_text(hjust = 0.5, face = "bold"))

save_plot(p, output, width = 12, height = max(8, nrow(mat_sub) * 0.2))
