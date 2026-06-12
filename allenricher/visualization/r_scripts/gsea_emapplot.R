#!/usr/bin/env Rscript
# gsea_emapplot.R — 通路富集网络图 (简化版: 基于 ggplot2 的通路关联散点图)
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
})
this_file <- sub("^--file=", "", commandArgs()[4])
source(file.path(dirname(this_file), "common.R"))

args <- parse_args()
tsv_path <- args$tsv
output <- args$output
top_n <- as.integer(args$top_n %||% 30)

df <- read_enrichment(tsv_path)

# 兼容不同列名: FDR / Adjusted_P_Value
if (!"Adjusted_P_Value" %in% colnames(df) && "FDR" %in% colnames(df)) {
  df$Adjusted_P_Value <- df$FDR
}

df <- df %>%
  mutate(nlogFDR = -log10(pmax(Adjusted_P_Value, 1e-10))) %>%
  arrange(desc(abs(NES))) %>%
  head(top_n)

# 使用 ggplot2 绘制富集图谱 (气泡图变体)
p <- ggplot(df, aes(x = NES, y = reorder(Term_Name, NES))) +
  geom_point(aes(size = Gene_Count, color = nlogFDR), alpha = 0.8) +
  geom_segment(aes(xend = 0, yend = Term_Name), alpha = 0.3, color = "gray50") +
  scale_color_gradient2(low = "blue", mid = "yellow", high = "red", midpoint = 1.3) +
  scale_size_continuous(range = c(3, 12)) +
  labs(x = "Normalized Enrichment Score (NES)", y = "",
       size = "Gene Count", color = expression(-log[10](FDR)),
       title = "Enrichment Map") +
  set_nature_theme() +
  theme(axis.text.y = element_text(size = 8),
        plot.title = element_text(hjust = 0.5, face = "bold"))

save_plot(p, output, width = 12, height = max(8, top_n * 0.3))
