# Shared helpers for AllEnricher GSEA R plots.
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
})

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0 || all(is.na(x)) || identical(x, "")) y else x
}

AE_CATEGORICAL_PALETTES <- list(
  tol_bright = c("#4477AA", "#66CCEE", "#228833", "#CCBB44", "#EE6677", "#AA3377", "#BBBBBB"),
  tol_high_contrast = c("#004488", "#DDAA33", "#BB5566"),
  tol_vibrant = c("#0077BB", "#33BBEE", "#009988", "#EE7733", "#CC3311", "#EE3377", "#BBBBBB"),
  tol_muted = c("#332288", "#88CCEE", "#44AA99", "#117733", "#999933", "#DDCC77", "#CC6677", "#882255", "#AA4499", "#DDDDDD"),
  tol_medium_contrast = c("#6699CC", "#004488", "#EECC66", "#994455", "#997700", "#EE99AA"),
  tol_light = c("#77AADD", "#99DDFF", "#44BB99", "#BBCC77", "#AAAA00", "#EEDD88", "#EE8866", "#FFAABB", "#DDDDDD"),
  okabe_ito = c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7"),
  nature = c("#0C5DA5", "#FF9500", "#00B945", "#FF2C00", "#845B97", "#474747", "#9E9E9E"),
  science = c("#1F77B4", "#FF7F0E", "#2CA02C", "#D62728", "#9467BD", "#8C564B", "#E377C2", "#7F7F7F", "#BCBD22", "#17BECF"),
  cell = c("#0072B2", "#D55E00", "#CC79A7", "#F0E442", "#009E73", "#56B4E9", "#E69F00", "#000000"),
  lancet = c("#00468B", "#ED0000", "#42B540", "#0099B4", "#925E9F", "#FDAF91", "#AD002A", "#ADB6B6"),
  nejm = c("#BC3C29", "#0072B5", "#E18727", "#20854E", "#7876B1", "#6F99AD", "#FFDC91", "#EE4C97"),
  jama = c("#374E55", "#DF8F44", "#00A1D5", "#B24745", "#79AF97", "#6A6599", "#80796B"),
  omicshare = c("#FF6B9D", "#C44569", "#F8B500", "#4ECDC4", "#556270", "#36D1DC", "#5AB9EA", "#8860D0"),
  echarts_v4 = c("#C23531", "#2F4554", "#61A0A8", "#D48265", "#91C7AE", "#749F83", "#CA8622", "#BDA29A", "#6E7074", "#546570")
)

AE_SEQUENTIAL_PALETTES <- list(
  colorbrewer_blues = c("#F7FBFF", "#DEEBF7", "#C6DBEF", "#9ECAE1", "#6BAED6", "#4292C6", "#2171B5", "#08519C", "#08306B"),
  colorbrewer_purd = c("#F7F4F9", "#E7E1EF", "#D4B9DA", "#C994C7", "#DF65B0", "#E7298A", "#CE1256", "#980043", "#67001F"),
  viridis = c("#440154", "#46327E", "#365C8D", "#277F8E", "#1FA187", "#4AC16D", "#A0DA39", "#FDE725"),
  cividis = c("#00204C", "#2E3F6D", "#575D6D", "#7C7B78", "#A59C74", "#D2C060", "#FFEA46")
)

AE_SEQUENTIAL_GRADIENT_ANCHORS <- list(
  colorbrewer_blues = c("#9ECAE1", "#08519C"),
  colorbrewer_purd = c("#DF65B0", "#980043"),
  viridis = c("#365C8D", "#1FA187"),
  cividis = c("#00204C", "#575D6D")
)

AE_DIVERGING_PALETTES <- list(
  colorbrewer_rdbu = c("#053061", "#2166AC", "#4393C3", "#92C5DE", "#D1E5F0", "#F7F7F7", "#FDDBC7", "#F4A582", "#D6604D", "#B2182B", "#67001F"),
  tol_sunset = c("#364B9A", "#4A7BB7", "#6EA6CD", "#98CAE1", "#C2E4EF", "#EAECCC", "#FEDA8B", "#FDB366", "#F67E4B", "#DD3D2D", "#A50026"),
  colorbrewer_prgn = c("#762A83", "#9970AB", "#C2A5CF", "#E7D4E8", "#F7F7F7", "#D9F0D3", "#ACD39E", "#5AAE61", "#1B7837"),
  colorbrewer_brbg = c("#543005", "#8C510A", "#BF812D", "#DFC27D", "#F6E8C3", "#F5F5F5", "#C7EAE5", "#80CDC1", "#35978F", "#01665E", "#003C30")
)

AE_PALETTES <- c(
  list(default = AE_CATEGORICAL_PALETTES$tol_bright),
  AE_CATEGORICAL_PALETTES,
  AE_SEQUENTIAL_PALETTES,
  AE_DIVERGING_PALETTES
)

AE_HIGH_CARDINALITY_CATEGORICAL <- c(
  "#1F77B4", "#FF7F0E", "#2CA02C", "#D62728", "#9467BD",
  "#8C564B", "#E377C2", "#7F7F7F", "#BCBD22", "#17BECF",
  "#393B79", "#637939", "#8C6D31", "#843C39", "#7B4173",
  "#3182BD", "#31A354", "#756BB1", "#E6550D", "#969696"
)

AE_STYLE_PRESETS <- list(
  nature = list(family = "sans", scale = 1.00, line = 1.00, spacing = 1.00, grid = FALSE, border = FALSE),
  science = list(family = "serif", scale = 1.05, line = 1.15, spacing = 1.00, grid = FALSE, border = TRUE),
  presentation = list(family = "sans", scale = 1.30, line = 1.40, spacing = 1.08, grid = TRUE, border = TRUE)
)
AE_STYLE_ALIASES <- c(cell = "nature", omicshare = "science")

palette_role <- function(name) {
  if (name %in% names(AE_CATEGORICAL_PALETTES)) return("categorical")
  if (name %in% names(AE_SEQUENTIAL_PALETTES)) return("sequential")
  if (name %in% names(AE_DIVERGING_PALETTES)) return("diverging")
  stop(paste("Unknown plot palette:", name))
}

configure_plot_style <- function(args = list()) {
  style <- tolower(args$style %||% "nature")
  if (style %in% names(AE_STYLE_ALIASES)) style <- unname(AE_STYLE_ALIASES[[style]])
  if (!style %in% names(AE_STYLE_PRESETS)) {
    stop(paste("Unknown plot style:", style))
  }
  preset <- AE_STYLE_PRESETS[[style]]

  palette_names <- list(
    categorical = "tol_bright",
    sequential = "colorbrewer_blues",
    diverging = "colorbrewer_rdbu"
  )
  legacy_name <- tolower(args$palette %||% "")
  if (legacy_name != "" && legacy_name != "default") {
    palette_names[[palette_role(legacy_name)]] <- legacy_name
  }
  explicit_names <- list(
    categorical = tolower(args$categorical_palette %||% ""),
    sequential = tolower(args$sequential_palette %||% ""),
    diverging = tolower(args$diverging_palette %||% "")
  )
  registries <- list(
    categorical = AE_CATEGORICAL_PALETTES,
    sequential = AE_SEQUENTIAL_PALETTES,
    diverging = AE_DIVERGING_PALETTES
  )
  for (role in names(explicit_names)) {
    name <- explicit_names[[role]]
    if (name == "") next
    if (!name %in% names(registries[[role]])) {
      stop(paste("Palette", name, "is not", role))
    }
    palette_names[[role]] <- name
  }

  AE_STYLE <<- style
  AE_CATEGORICAL_PALETTE_NAME <<- palette_names$categorical
  AE_SEQUENTIAL_PALETTE_NAME <<- palette_names$sequential
  AE_DIVERGING_PALETTE_NAME <<- palette_names$diverging
  AE_CATEGORICAL_PALETTE <<- AE_CATEGORICAL_PALETTES[[palette_names$categorical]]
  AE_SEQUENTIAL_PALETTE <<- AE_SEQUENTIAL_PALETTES[[palette_names$sequential]]
  AE_DIVERGING_PALETTE <<- AE_DIVERGING_PALETTES[[palette_names$diverging]]
  AE_PALETTE_NAME <<- AE_CATEGORICAL_PALETTE_NAME
  AE_PALETTE <<- AE_CATEGORICAL_PALETTE
  AE_FONT_FAMILY <<- preset$family
  AE_TEXT_SCALE <<- preset$scale
  AE_LINE_SCALE <<- preset$line
  AE_SPACING_SCALE <<- preset$spacing
  AE_SHOW_GRID <<- preset$grid
  AE_FULL_BORDER <<- preset$border
  visible_diverging <- ae_visible_gradient_colors(AE_DIVERGING_PALETTE, "diverging")
  AE_COL_DOWN <<- visible_diverging[1]
  AE_COL_UP <<- visible_diverging[length(visible_diverging)]
  AE_COL_NEUTRAL <<- "#4D4D4D"
  AE_COL_GRID <<- if (preset$grid) "#E1E4E8" else "#ECEEF0"
  AE_COL_TEXT <<- "#222222"
  AE_COL_LOW <<- AE_COL_DOWN
  AE_COL_MID <<- visible_diverging[ceiling(length(visible_diverging) / 2)]
  AE_COL_HIGH <<- AE_COL_UP
  invisible(NULL)
}

ae_palette <- function(n) {
  if (n <= length(AE_CATEGORICAL_PALETTE)) {
    return(AE_CATEGORICAL_PALETTE[seq_len(n)])
  }
  if (n > length(AE_HIGH_CARDINALITY_CATEGORICAL)) {
    stop(paste("Categorical plot requires", n, "distinct colors; maximum supported is", length(AE_HIGH_CARDINALITY_CATEGORICAL)))
  }
  warning(paste("Categorical palette", AE_CATEGORICAL_PALETTE_NAME, "has too few colors; using high-cardinality fallback"))
  AE_HIGH_CARDINALITY_CATEGORICAL[seq_len(n)]
}

ae_visible_gradient_colors <- function(colors, role, palette_name = NULL) {
  if (!role %in% c("sequential", "diverging")) {
    stop(paste("Gradient visibility requires a continuous role, got:", role))
  }
  distance_from_white <- vapply(colors, function(color) {
    rgb <- grDevices::col2rgb(color)[, 1] / 255
    sqrt(sum((1 - rgb)^2))
  }, numeric(1))
  near_white <- distance_from_white < 0.18
  if (role == "sequential") {
    if (!is.null(palette_name) && palette_name %in% names(AE_SEQUENTIAL_GRADIENT_ANCHORS)) {
      return(AE_SEQUENTIAL_GRADIENT_ANCHORS[[palette_name]])
    }
    visible <- colors[!near_white]
    if (length(visible) >= 2) return(visible)
  }
  if (role == "diverging") {
    colors[near_white] <- "#FFFFFF"
    return(colors)
  }
  colors[near_white] <- "#D9D9D9"
  colors
}

ae_sequential_colors <- function(n = 256) {
  anchors <- ae_visible_gradient_colors(AE_SEQUENTIAL_PALETTE, "sequential", AE_SEQUENTIAL_PALETTE_NAME)
  grDevices::colorRampPalette(anchors)(n)
}

ae_gradient_colors <- ae_sequential_colors

ae_diverging_colors <- function(n = 256) {
  anchors <- ae_visible_gradient_colors(AE_DIVERGING_PALETTE, "diverging")
  grDevices::colorRampPalette(anchors)(n)
}

configure_plot_style()

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
  if ("padj" %in% colnames(df) && !"Adjusted_P_Value" %in% colnames(df)) {
    df$Adjusted_P_Value <- df$padj
  }
  if (!"Adjusted_P_Value" %in% colnames(df)) {
    df$Adjusted_P_Value <- 1
  }

  if ("pvalue" %in% colnames(df) && !"p_value" %in% colnames(df)) {
    df$p_value <- df$pvalue
  }
  if ("pval" %in% colnames(df) && !"p_value" %in% colnames(df)) {
    df$p_value <- df$pval
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
    } else if ("leadingEdge" %in% colnames(df)) {
      df$Genes <- df$leadingEdge
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
  if ("size" %in% colnames(df) && !"Gene_Count" %in% colnames(df)) {
    df$Gene_Count <- df$size
  }
  if (!"Gene_Count" %in% colnames(df)) {
    df$Gene_Count <- lengths(parse_gene_list(df$Genes))
  }

  if (!"Term_ID" %in% colnames(df)) {
    id_col <- intersect(c("pathway", "ID", "id", "term_id"), colnames(df))[1]
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

  if (!"EnrichFactor" %in% colnames(df)) {
    enrich_col <- intersect(c("enrichfactor", "RichFactor", "rich_factor", "Rich_Factor"), colnames(df))[1]
    if (!is.na(enrich_col)) {
      df$EnrichFactor <- df[[enrich_col]]
    }
  }

  df$NES <- as_num(df$NES, 0)
  df$ES <- as_num(df$ES, df$NES)
  df$Adjusted_P_Value <- pmax(as_num(df$Adjusted_P_Value, 1), .Machine$double.xmin)
  df$Gene_Count <- pmax(as_num(df$Gene_Count, 0), 0)
  if ("EnrichFactor" %in% colnames(df)) {
    df$EnrichFactor <- as_num(df$EnrichFactor, NA_real_)
  }
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
    label
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

prepare_gsea_terms <- function(df, top_n = 20, sort_by = "abs_nes", label_width = 34) {
  df <- df %>%
    mutate(
      nlogFDR = pmin(-log10(pmax(Adjusted_P_Value, .Machine$double.xmin)), 50),
      Direction = ifelse(NES >= 0, "Positive NES", "Negative NES"),
      Term_Short = clean_pathway_name(Term_Name)
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

format_sig_short <- function(x, digits = 1) {
  vapply(x, function(value) {
    if (is.na(value) || !is.finite(value)) {
      return("NA")
    }
    text <- formatC(value, format = "e", digits = digits)
    text <- sub("e([+-])0+", "e\\1", text)
    text <- sub("e\\+", "e", text)
    sub("\\.0e", "e", text)
  }, character(1))
}

prepare_gsea_diverging_terms <- function(df, top_n = 20, label_width = 30, sig_cutoff = 0.05) {
  df <- df %>%
    mutate(
      NES = as_num(NES, NA_real_),
      Adjusted_P_Value = pmax(as_num(Adjusted_P_Value, NA_real_), .Machine$double.xmin),
      Gene_Count = as_num(Gene_Count, NA_real_),
      nlogFDR = pmin(-log10(Adjusted_P_Value), 50),
      Term_Short = clean_pathway_name(Term_Name)
    ) %>%
    filter(is.finite(NES), is.finite(Adjusted_P_Value), Adjusted_P_Value > 0)

  if (nrow(df) == 0) {
    stop("No valid GSEA rows available for diverging barplot")
  }

  plot_df <- df %>% filter(Adjusted_P_Value < sig_cutoff)
  if (nrow(plot_df) == 0) {
    plot_df <- df %>% arrange(Adjusted_P_Value) %>% head(top_n)
  }

  top_each <- max(1, floor(top_n / 2))
  up <- plot_df %>% filter(NES > 0) %>% arrange(desc(NES)) %>% head(top_each)
  down <- plot_df %>% filter(NES < 0) %>% arrange(NES) %>% head(top_each) %>% arrange(desc(NES))
  plot_df <- bind_rows(up, down)

  if (nrow(plot_df) == 0) {
    stop("No positive or negative NES terms to plot")
  }

  plot_df$Term_Label <- make_unique_labels(wrap_label(plot_df$Term_Short, label_width), plot_df$Term_ID)
  plot_df$Value_Label <- paste0(
    ifelse(is.na(plot_df$Gene_Count), "NA", as.integer(plot_df$Gene_Count)),
    " / ",
    format_sig_short(plot_df$Adjusted_P_Value)
  )
  plot_df
}

plot_gsea_diverging_bar <- function(df, output, top_n = 20, title = "GSEA NES Ranking", dpi = 300) {
  plot_df <- prepare_gsea_diverging_terms(df, top_n = top_n)
  max_abs <- max(abs(plot_df$NES), na.rm = TRUE)
  if (!is.finite(max_abs) || max_abs <= 0) {
    max_abs <- 1
  }
  offset <- max_abs * 0.03
  plot_df$Term_Label <- factor(plot_df$Term_Label, levels = rev(plot_df$Term_Label))
  plot_df$Term_X <- ifelse(plot_df$NES >= 0, -offset, offset)
  plot_df$Term_Hjust <- ifelse(plot_df$NES >= 0, 1, 0)
  plot_df$Value_X <- ifelse(plot_df$NES >= 0, plot_df$NES + offset, plot_df$NES - offset)
  plot_df$Value_Hjust <- ifelse(plot_df$NES >= 0, 0, 1)
  fill_limits <- range(plot_df$nlogFDR, na.rm = TRUE)
  if (!all(is.finite(fill_limits))) {
    fill_limits <- c(0, 1)
  } else if (diff(fill_limits) < .Machine$double.eps^0.5) {
    fill_limits <- fill_limits + c(-0.5, 0.5)
  }

  p <- ggplot(plot_df, aes(y = Term_Label)) +
    geom_vline(xintercept = 0, color = "grey20", linewidth = ae_lwd(0.55)) +
    geom_col(aes(x = NES, fill = pmin(nlogFDR, 10)), width = 0.62, color = "grey25", linewidth = ae_lwd(0.25)) +
    geom_text(aes(x = Term_X, label = Term_Label, hjust = Term_Hjust), size = 2.55 * AE_TEXT_SCALE, lineheight = 0.92) +
    geom_text(aes(x = Value_X, label = Value_Label, hjust = Value_Hjust), size = 2.25 * AE_TEXT_SCALE) +
    annotate("text", x = max_abs * 0.62, y = nrow(plot_df) + 0.75, label = "UP", color = AE_COL_UP, fontface = "bold", size = 3.0 * AE_TEXT_SCALE) +
    annotate("text", x = -max_abs * 0.62, y = nrow(plot_df) + 0.75, label = "DOWN", color = AE_COL_DOWN, fontface = "bold", size = 3.0 * AE_TEXT_SCALE) +
    scale_fill_gradientn(colors = ae_gradient_colors(6), limits = fill_limits, name = expression(-log[10](FDR))) +
    guides(fill = guide_colorbar(barheight = grid::unit(18, "mm"), barwidth = grid::unit(3, "mm"), title.position = "top")) +
    scale_x_continuous(limits = c(-max_abs * 1.48, max_abs * 1.56), expand = expansion(mult = c(0, 0))) +
    coord_cartesian(clip = "off") +
    labs(x = "Normalized enrichment score (NES)", y = NULL, title = paste0(title, " (core# / FDR)")) +
    set_nature_theme(base_size = 9) +
    theme(
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.line.y = element_blank(),
      legend.position = "right",
      panel.grid.major.x = if (AE_SHOW_GRID) element_line(
        color = AE_COL_GRID, linewidth = ae_lwd(0.35)
      ) else element_blank(),
      plot.margin = ae_margin(14, 18, 8, 16)
    )

  save_plot(p, output, width = 7.45, height = max(4.25, nrow(plot_df) * 0.205 + 1.35), dpi = dpi)
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
  ae_palette(n)
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
    colors = ae_gradient_colors(6),
    ...
  )
}

scale_fill_nlogfdr <- function(...) {
  scale_fill_gradientn(
    colors = ae_gradient_colors(6),
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

ae_lwd <- function(value) value * AE_LINE_SCALE

ae_margin <- function(top, right, bottom, left) {
  margin(
    top * AE_SPACING_SCALE, right * AE_SPACING_SCALE,
    bottom * AE_SPACING_SCALE, left * AE_SPACING_SCALE
  )
}

ae_panel_box <- function(lwd = 0.75, col = AE_COL_NEUTRAL) {
  box(bty = if (AE_FULL_BORDER) "o" else "l", lwd = ae_lwd(lwd), col = col)
}

ae_matrix_box <- function(lwd = 0.75, col = AE_COL_NEUTRAL) {
  box(bty = "o", lwd = ae_lwd(lwd), col = col)
}

ae_background_grid <- function(h = NULL, v = NULL, lwd = 0.7, col = AE_COL_GRID) {
  if (!AE_SHOW_GRID) return(invisible(NULL))
  if (!is.null(h)) abline(h = h, col = col, lwd = ae_lwd(lwd))
  if (!is.null(v)) abline(v = v, col = col, lwd = ae_lwd(lwd))
  invisible(NULL)
}

set_nature_theme <- function(base_size = 9) {
  scaled_size <- base_size * AE_TEXT_SCALE
  base_theme <- if (AE_SHOW_GRID) theme_minimal else theme_classic
  base_theme(base_size = scaled_size, base_family = AE_FONT_FAMILY) +
    theme(
      text = element_text(color = AE_COL_TEXT),
      axis.text = element_text(color = AE_COL_TEXT),
      axis.title = element_text(color = AE_COL_TEXT),
      axis.line = if (AE_FULL_BORDER) element_blank() else element_line(color = "#5F6368", linewidth = ae_lwd(0.35)),
      axis.ticks = element_line(color = "#5F6368", linewidth = ae_lwd(0.3)),
      panel.border = if (AE_FULL_BORDER) element_rect(color = "#5F6368", fill = NA, linewidth = ae_lwd(0.4)) else element_blank(),
      panel.grid.major.x = if (AE_SHOW_GRID) element_line(color = AE_COL_GRID, linewidth = ae_lwd(0.35)) else element_blank(),
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "right",
      legend.title = element_text(size = scaled_size, face = "bold"),
      legend.text = element_text(size = max(6, scaled_size - 1)),
      plot.title = element_text(face = "bold", size = scaled_size + 2, hjust = 0),
      plot.subtitle = element_text(size = scaled_size, color = "#555555", hjust = 0),
      plot.caption = element_text(size = max(6, scaled_size - 1), color = "#777777", hjust = 1),
      plot.margin = ae_margin(8, 12, 8, 8)
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
    svg(output_path, width = width, height = height, bg = "white")
    print(plot_obj)
    dev.off()
  } else {
    stop(paste("Unsupported output format:", ext))
  }
  message(paste("Saved:", output_path))
}
