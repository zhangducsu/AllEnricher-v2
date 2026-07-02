# Shared helpers for AllEnricher GSEA R plots.
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
})

AE_COL_UP <- "#D55E00"
AE_COL_DOWN <- "#0072B2"
AE_COL_NEUTRAL <- "#4D4D4D"
AE_COL_GRID <- "#E6E8EB"
AE_COL_TEXT <- "#222222"
AE_COL_LOW <- "#D8E8F4"
AE_COL_MID <- "#F7F7F7"
AE_COL_HIGH <- "#B2182B"

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0 || all(is.na(x)) || identical(x, "")) y else x
}

parse_args <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  arg_list <- list()
  if (length(args) == 0) {
    return(arg_list)
  }

  i <- 1
  while (i <= length(args)) {
    if (grepl("^--", args[i])) {
      key <- gsub("^--", "", args[i])
      value <- if (i < length(args) && !grepl("^--", args[i + 1])) args[i + 1] else ""
      arg_list[[key]] <- value
      i <- i + 2
    } else {
      i <- i + 1
    }
  }
  arg_list
}

as_num <- function(x, default = NA_real_) {
  value <- suppressWarnings(as.numeric(x))
  value[is.na(value)] <- default
  value
}

parse_gene_list <- function(values) {
  lapply(as.character(values), function(value) {
    genes <- unlist(strsplit(value, "[,;/|]"))
    genes <- trimws(genes)
    unique(genes[genes != "" & !is.na(genes)])
  })
}

read_enrichment <- function(tsv_path) {
  df <- read.delim(tsv_path, sep = "\t", comment.char = "#", stringsAsFactors = FALSE, check.names = FALSE)

  if ("FDR" %in% colnames(df) && !"Adjusted_P_Value" %in% colnames(df)) {
    df$Adjusted_P_Value <- df$FDR
  }
  if ("FDR q-val" %in% colnames(df) && !"Adjusted_P_Value" %in% colnames(df)) {
    df$Adjusted_P_Value <- df$`FDR q-val`
  }
  if (!"Adjusted_P_Value" %in% colnames(df)) {
    df$Adjusted_P_Value <- 1
  }

  if ("pvalue" %in% colnames(df) && !"p_value" %in% colnames(df)) {
    df$p_value <- df$pvalue
  }
  if ("P_Value" %in% colnames(df) && !"p_value" %in% colnames(df)) {
    df$p_value <- df$P_Value
  }
  if ("NOM p-val" %in% colnames(df) && !"p_value" %in% colnames(df)) {
    df$p_value <- df$`NOM p-val`
  }

  if (!"Genes" %in% colnames(df)) {
    if ("Lead_genes" %in% colnames(df)) {
      df$Genes <- df$Lead_genes
    } else if ("matched_genes" %in% colnames(df)) {
      df$Genes <- df$matched_genes
    } else if ("core_enrichment" %in% colnames(df)) {
      df$Genes <- df$core_enrichment
    } else {
      df$Genes <- ""
    }
  }

  if ("setSize" %in% colnames(df) && !"Gene_Count" %in% colnames(df)) {
    df$Gene_Count <- df$setSize
  }
  if (!"Gene_Count" %in% colnames(df)) {
    df$Gene_Count <- lengths(parse_gene_list(df$Genes))
  }

  if (!"Term_ID" %in% colnames(df)) {
    id_col <- intersect(c("ID", "id", "term_id"), colnames(df))[1]
    df$Term_ID <- if (!is.na(id_col)) df[[id_col]] else seq_len(nrow(df))
  }
  if (!"Term_Name" %in% colnames(df)) {
    name_col <- intersect(c("Description", "pathway", "term_name"), colnames(df))[1]
    df$Term_Name <- if (!is.na(name_col)) df[[name_col]] else df$Term_ID
  }

  if (!"NES" %in% colnames(df)) {
    if ("nes" %in% colnames(df)) {
      df$NES <- df$nes
    } else if ("ES" %in% colnames(df)) {
      df$NES <- df$ES
    } else {
      df$NES <- 0
    }
  }

  if (!"ES" %in% colnames(df)) {
    df$ES <- df$NES
  }

  df$NES <- as_num(df$NES, 0)
  df$ES <- as_num(df$ES, df$NES)
  df$Adjusted_P_Value <- pmax(as_num(df$Adjusted_P_Value, 1), .Machine$double.xmin)
  df$Gene_Count <- pmax(as_num(df$Gene_Count, 0), 0)
  df
}

clean_pathway_name <- function(values) {
  vapply(as.character(values), function(value) {
    original <- trimws(value)
    label <- trimws(sub(".*\\|", "", original))
    if (label == "" || is.na(label)) {
      label <- original
    }
    label <- gsub("_", " ", label)
    label <- gsub("\\s+", " ", label)
    label <- normalize_pathway_case(label)
    label
  }, character(1))
}

normalize_pathway_case <- function(value) {
  label <- trimws(gsub("\\s+", " ", as.character(value)))
  if (label == "" || is.na(label)) {
    return(label)
  }
  words <- strsplit(label, "\\s+")[[1]]
  starts_upper <- sum(grepl("^[A-Z]", words))
  looks_title <- length(words) > 1 && starts_upper >= max(2, floor(length(words) / 2))
  if (!grepl("^[a-z]", label) && !looks_title) {
    return(label)
  }
  small_words <- c("a", "an", "and", "as", "at", "by", "for", "from", "in", "of", "on", "or", "the", "to", "via", "with")
  formatted <- vapply(seq_along(words), function(i) {
    word <- words[[i]]
    if (nchar(word) > 1 && toupper(word) == word && grepl("[A-Z]", word)) {
      return(word)
    }
    lower <- tolower(word)
    if (i > 1) {
      return(lower)
    }
    paste0(toupper(substr(lower, 1, 1)), substr(lower, 2, nchar(lower)))
  }, character(1))
  paste(formatted, collapse = " ")
}

shorten_label <- function(values, max_chars = 80) {
  if (is.null(max_chars) || is.na(max_chars) || max_chars <= 0) {
    return(as.character(values))
  }
  vapply(as.character(values), function(value) {
    value <- trimws(gsub("\\s+", " ", value))
    if (nchar(value) <= max_chars) {
      return(value)
    }
    paste0(substr(value, 1, max_chars - 3), "...")
  }, character(1))
}

wrap_label <- function(values, width = 34) {
  vapply(as.character(values), function(value) {
    paste(strwrap(value, width = width), collapse = "\n")
  }, character(1))
}

make_unique_labels <- function(labels, ids) {
  labels <- as.character(labels)
  ids <- as.character(ids)
  duplicated_labels <- labels %in% labels[duplicated(labels)]
  labels[duplicated_labels] <- paste0(labels[duplicated_labels], " (", ids[duplicated_labels], ")")
  labels
}

prepare_gsea_terms <- function(df, top_n = 20, sort_by = "abs_nes", label_width = 34, max_label_chars = 84) {
  df <- df %>%
    mutate(
      nlogFDR = pmin(-log10(pmax(Adjusted_P_Value, .Machine$double.xmin)), 50),
      Direction = ifelse(NES >= 0, "Positive NES", "Negative NES"),
      Term_Short = shorten_label(clean_pathway_name(Term_Name), max_label_chars)
    )

  if (sort_by == "q") {
    df <- df %>% arrange(Adjusted_P_Value, desc(abs(NES)))
  } else {
    df <- df %>% arrange(desc(abs(NES)), Adjusted_P_Value)
  }

  df <- head(df, top_n)
  df$Term_Label <- make_unique_labels(wrap_label(df$Term_Short, label_width), df$Term_ID)
  df
}

is_hit_value <- function(values) {
  as.character(values) %in% c("TRUE", "True", "true", "1", "T") | values == TRUE
}

read_running_es <- function(running_es_path) {
  if (is.null(running_es_path) || running_es_path == "" || !file.exists(running_es_path)) {
    stop("Real running-ES file is required for this plot")
  }
  df <- read.delim(running_es_path, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
  df$Rank <- as_num(df$Rank, 0)
  df$Running_ES <- as_num(df$Running_ES, 0)
  df$Weight <- as_num(df$Weight, 0)
  df$Hit <- is_hit_value(df$Hit)
  df
}

pathway_palette <- function(n) {
  base <- c(
    "#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00",
    "#56B4E9", "#F0E442", "#000000", "#8DA0CB", "#A6761D"
  )
  rep(base, length.out = n)
}

scale_fill_direction <- function(...) {
  scale_fill_manual(
    values = c("Positive NES" = AE_COL_UP, "Negative NES" = AE_COL_DOWN),
    ...
  )
}

scale_color_direction <- function(...) {
  scale_color_manual(
    values = c("Positive NES" = AE_COL_UP, "Negative NES" = AE_COL_DOWN),
    ...
  )
}

scale_color_nlogfdr <- function(...) {
  scale_color_gradientn(
    colors = c("#D8E8F4", "#74A9CF", "#F1A340", "#B2182B"),
    ...
  )
}

scale_fill_nes <- function(...) {
  scale_fill_gradient2(
    low = AE_COL_DOWN,
    mid = AE_COL_MID,
    high = AE_COL_UP,
    midpoint = 0,
    ...
  )
}

set_nature_theme <- function(base_size = 9) {
  theme_classic(base_size = base_size, base_family = "sans") +
    theme(
      text = element_text(color = AE_COL_TEXT),
      axis.text = element_text(color = AE_COL_TEXT),
      axis.title = element_text(color = AE_COL_TEXT),
      axis.line = element_line(color = "#5F6368", linewidth = 0.35),
      axis.ticks = element_line(color = "#5F6368", linewidth = 0.3),
      panel.grid.major.x = element_line(color = AE_COL_GRID, linewidth = 0.28),
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "right",
      legend.title = element_text(size = base_size, face = "bold"),
      legend.text = element_text(size = base_size - 1),
      plot.title = element_text(face = "bold", size = base_size + 2, hjust = 0),
      plot.subtitle = element_text(size = base_size, color = "#555555", hjust = 0),
      plot.caption = element_text(size = base_size - 1, color = "#777777", hjust = 1),
      plot.margin = margin(8, 12, 8, 8)
    )
}

save_plot <- function(plot_obj, output_path, width = 8, height = 6, dpi = 300) {
  dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
  ext <- tolower(tools::file_ext(output_path))
  if (ext == "png") {
    ggsave(output_path, plot = plot_obj, width = width, height = height, units = "in", dpi = dpi, limitsize = FALSE, bg = "white")
  } else if (ext == "pdf") {
    ggsave(output_path, plot = plot_obj, width = width, height = height, device = "pdf", limitsize = FALSE, bg = "white")
  } else if (ext == "svg") {
    grDevices::svg(output_path, width = width, height = height, bg = "white")
    on.exit(grDevices::dev.off(), add = TRUE)
    print(plot_obj)
  } else {
    stop(paste("Unsupported output format:", ext))
  }
  message(paste("Saved:", output_path))
}
