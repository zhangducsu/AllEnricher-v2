#!/usr/bin/env Rscript
# gsea_circos.R — 环形富集图 (基于 ggplot2 的极坐标柱状图)
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

df <- read_enrichment(tsv_path) %>%
  arrange(desc(abs(NES))) %>%
  head(top_n) %>%
  mutate(
    # 取 Term_Name 最后一段作为短标签
    short_name = sub(".*\\|", "", Term_Name),
    direction = ifelse(NES > 0, "Up", "Down"),
    abs_nes = abs(NES)
  )

# 极坐标柱状图
p <- ggplot(df, aes(x = reorder(short_name, abs_nes), y = abs_nes, fill = direction)) +
  geom_col(width = 0.8) +
  coord_polar() +
  scale_fill_manual(values = c("Up" = "#E41A1C", "Down" = "#377EBA")) +
  labs(x = "", y = "|NES|", fill = "Regulation",
       title = "Circos Plot - Pathway Enrichment") +
  set_nature_theme() +
  theme(
    axis.text.x = element_text(size = 6, angle = 0),
    legend.position = "top",
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

save_plot(p, output, width = 10, height = 10)
