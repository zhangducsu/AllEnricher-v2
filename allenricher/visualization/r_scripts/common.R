# common.R — GSEA R 可视化公共函数
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
})

# 解析命令行参数
parse_args <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  arg_list <- list()
  for (i in seq(1, length(args), 2)) {
    if (i < length(args) && grepl("^--", args[i])) {
      arg_list[[gsub("^--", "", args[i])]] <- args[i + 1]
    }
  }
  return(arg_list)
}

# 读取 AllEnricher 输出的 TSV (自动标准化列名)
read_enrichment <- function(tsv_path) {
  df <- read.delim(tsv_path, sep = "\t", comment.char = "#", stringsAsFactors = FALSE, check.names = FALSE)

  # 标准化列名: FDR -> Adjusted_P_Value
  if ("FDR" %in% colnames(df) && !"Adjusted_P_Value" %in% colnames(df)) {
    df$Adjusted_P_Value <- df$FDR
  }
  # p_value / pvalue / NOM p-val -> p_value
  if ("pvalue" %in% colnames(df) && !"p_value" %in% colnames(df)) {
    df$p_value <- df$pvalue
  }
  if ("NOM p-val" %in% colnames(df) && !"p_value" %in% colnames(df)) {
    df$p_value <- df$`NOM p-val`
  }
  # Gene_Count / setSize -> Gene_Count
  if ("setSize" %in% colnames(df) && !"Gene_Count" %in% colnames(df)) {
    df$Gene_Count <- df$setSize
  }
  # Genes / Lead_genes / matched_genes -> Genes
  if ("Genes" %in% colnames(df)) {
    # 已有标准列名
  } else if ("Lead_genes" %in% colnames(df)) {
    df$Genes <- df$Lead_genes
  } else if ("matched_genes" %in% colnames(df)) {
    df$Genes <- df$matched_genes
  }

  return(df)
}

# 保存图表
save_plot <- function(plot_obj, output_path, width = 8, height = 6, dpi = 300) {
  ext <- tools::file_ext(output_path)
  if (ext == "png") {
    ggsave(output_path, plot = plot_obj, width = width, height = height, units = "in", dpi = dpi, limitsize = FALSE)
  } else if (ext == "pdf") {
    ggsave(output_path, plot = plot_obj, width = width, height = height, device = "pdf")
  } else if (ext == "svg") {
    ggsave(output_path, plot = plot_obj, width = width, height = height, device = "svg")
  }
  message(paste("Saved:", output_path))
}

# Nature 风格主题
set_nature_theme <- function() {
  theme_set(theme_bw(base_size = 12) +
    theme(
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "grey80", fill = NA),
      axis.text = element_text(color = "black"),
      axis.title = element_text(size = 13),
      legend.title = element_text(size = 11),
      legend.text = element_text(size = 10)
    ))
}
