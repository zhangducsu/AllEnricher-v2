#!/usr/bin/env Rscript
# gsea_dotplot.R — GSEA 气泡图 (Y=通路, X=NES, 大小=基因数, 颜色=-log10(FDR))
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
})
this_file <- sub("^--file=", "", commandArgs()[4])
source(file.path(dirname(this_file), "common.R"))

args <- parse_args()
tsv_path <- args$tsv
output <- args$output
top_n <- as.integer(args$top_n %||% 20)

df <- read_enrichment(tsv_path)
df <- df %>%
  mutate(nlogFDR = -log10(pmax(Adjusted_P_Value, 1e-10))) %>%
  arrange(desc(abs(NES))) %>%
  head(top_n)

p <- ggplot(df, aes(x = NES, y = reorder(Term_Name, NES))) +
  geom_point(aes(size = Gene_Count, color = nlogFDR)) +
  scale_color_gradient2(low = "blue", mid = "yellow", high = "red", midpoint = 1.3) +
  scale_size_continuous(range = c(3, 10)) +
  labs(x = "Normalized Enrichment Score (NES)", y = "",
       size = "Gene Count", color = expression(-log[10](FDR))) +
  set_nature_theme() +
  theme(axis.text.y = element_text(size = 9))

save_plot(p, output, width = 10, height = max(6, top_n * 0.3))
