#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(ggupset)
  library(ggpattern)
  library(ggnewscale)
  library(readr)
  library(dplyr)
  library(tidyr)
  library(grid)
  library(stringr)
  library(patchwork)
  library(purrr)
  library(svglite)
})

# -------------------------
# Helpers
# -------------------------

# Classify overlap category between methods
classify_overlap <- function(methods) {
  if (setequal(methods, c("Coinfinder"))) return("Coinfinder only")
  if (setequal(methods, c("Goldfinder"))) return("Goldfinder only")
  if (setequal(methods, c("PanForest")))  return("PanForest only")
  if (setequal(methods, c("Coinfinder", "Goldfinder")))             return("Coinfinder + Goldfinder")
  if (setequal(methods, c("Coinfinder", "PanForest")))              return("Coinfinder + PanForest")
  if (setequal(methods, c("Goldfinder", "PanForest")))              return("Goldfinder + PanForest")
  if (setequal(methods, c("Coinfinder", "Goldfinder", "PanForest")))return("All three methods")
  stop("Unexpected MethodsPresent combination: ", paste(methods, collapse = ", "))
}

# Safe loader for D-value CSVs
safe_read <- function(path, method_name) {
  if (!is.null(path) && file.exists(path)) {
    readr::read_csv(path, show_col_types = FALSE) %>%
      dplyr::filter(Level %in% c("Gene", "Pair")) %>%
      dplyr::transmute(
        Level,
        Gene,
        Gene_1,
        Gene_2,
        Pair = dplyr::if_else(Level == "Pair",
                              paste(Gene_1, Gene_2, sep = "_"),
                              NA_character_),
        D_value,
        Method = method_name
      )
  } else {
    message("Warning: file not found for ", method_name, " (", path, ") - continuing without values")
    tibble::tibble(
      Level  = character(),
      Gene   = character(),
      Gene_1 = character(),
      Gene_2 = character(),
      Pair   = character(),
      D_value= numeric(),
      Method = character()
    )
  }
}

# Build pairs dataframe with PairID and PairType
prepare_pairs_df <- function(df_all) {
  df_all %>%
    dplyr::filter(Level == "Pair") %>%
    dplyr::mutate(
      Gene_1_clean = stringr::str_remove(Gene_1, "_dup$"),
      Gene_2_clean = stringr::str_remove(Gene_2, "_dup$"),
      is_dup_1     = stringr::str_detect(Gene_1, "_dup$"),
      is_dup_2     = stringr::str_detect(Gene_2, "_dup$"),
      PairType = dplyr::case_when(
        Gene_1_clean == Gene_2_clean & (is_dup_1 | is_dup_2) ~ "Correct_dup_pair",
        is_dup_1 | is_dup_2                                   ~ "Incorrect_dup_pair",
        TRUE                                                  ~ "Non_dup_pair"
      ),
      PairID = paste(pmin(Gene_1, Gene_2), pmax(Gene_1, Gene_2), sep = "-")
    )
}

# Generalized builder for UpSet input from a dataframe and an ID column
make_upset_df <- function(df, id_col, filter_expr = TRUE, method_order = c("Coinfinder","Goldfinder","PanForest")) {
  df %>%
    dplyr::filter({{ filter_expr }}) %>%
    dplyr::distinct(Method, {{ id_col }}) %>%
    dplyr::group_by({{ id_col }}) %>%
    dplyr::summarise(sets = list(sort(unique(Method))), .groups = "drop") %>%
    dplyr::count(sets, name = "value") %>%
    dplyr::mutate(
      sets_ord    = lapply(sets, function(s) intersect(method_order, s)),
      n_sets      = lengths(sets_ord),
      base_tool   = vapply(sets_ord, function(s) s[1], character(1)),
      overlay_tool= dplyr::case_when(
        n_sets == 2 ~ vapply(sets_ord, function(s) s[2], character(1)),
        n_sets == 3 ~ "TRIPLE_WHITE",
        TRUE        ~ NA_character_
      ),
      pattern     = dplyr::case_when(
        n_sets == 1 ~ "none",
        n_sets == 2 ~ "stripe",
        n_sets == 3 ~ "crosshatch"
      )
    )
}

# UpSet plot generator
make_upset_plot <- function(df2, title, ylab, outfile, method_colors, out_dir, show_pattern_legend = TRUE) {
  pattern_fill_colors <- c(
    Coinfinder   = method_colors[["Coinfinder"]],
    Goldfinder   = method_colors[["Goldfinder"]],
    PanForest    = method_colors[["PanForest"]],
    TRIPLE_WHITE = "#FFFFFF"
  )

  df2 <- df2 %>%
    dplyr::mutate(
      base_tool    = factor(base_tool, levels = names(method_colors)),
      overlay_tool = factor(overlay_tool, levels = c(names(method_colors), "TRIPLE_WHITE")),
      pattern      = factor(pattern, levels = c("none","stripe","crosshatch"))
    )

  p <- ggplot(df2, aes(x = sets, y = value)) +
    geom_col(aes(fill = base_tool), colour = "black", width = 0.7, show.legend = c(fill = TRUE)) +
    geom_col(data = dplyr::filter(df2, n_sets == 3), aes(y = value),
             fill = "#999999", colour = "black", width = 0.7, show.legend = FALSE) +
    geom_col_pattern(
      aes(pattern = pattern, pattern_fill = overlay_tool),
      fill = NA, colour = "black", pattern_colour = "black",
      width = 0.7,
      show.legend = c(pattern = show_pattern_legend, pattern_fill = FALSE, fill = FALSE),
      key_glyph = ggpattern::draw_key_polygon_pattern,
      pattern_angle = 45, pattern_density = 0.5, pattern_spacing = 0.05
    ) +
    geom_label(aes(y = value, label = value),
               vjust = -0.3, size = 3, label.size = 0, fill = "white") +
    scale_x_upset(order_by = "degree") +
    scale_fill_manual(name = "Method", values = method_colors, limits = names(method_colors),
                      breaks = names(method_colors), drop = FALSE) +
    scale_pattern_manual(values = c(none = "none", stripe = "stripe", crosshatch = "crosshatch"), guide = "none") +
    scale_pattern_fill_manual(values = pattern_fill_colors, na.value = NA, guide = "none") +
    labs(x = "Method intersections", y = ylab, title = title) +
    theme_minimal() +
    theme(legend.key = element_rect(fill = "white", colour = NA))

  ggsave(file.path(out_dir, outfile), p, width = 8, height = 5, dpi = 600, bg = "white")
}

# Relative counts panel builder for proportions
make_relative_counts <- function(df, id_col, label, method_colors) {
  all_methods <- tibble::tibble(Method = names(method_colors))

  totals <- df %>%
    dplyr::distinct(Method, {{ id_col }}) %>%
    dplyr::count(Method, name = "total")

  dups <- df %>%
    dplyr::filter(grepl("_dup$", {{ id_col }})) %>%
    dplyr::distinct(Method, {{ id_col }}) %>%
    dplyr::count(Method, name = "dup")

  totals %>%
    dplyr::left_join(dups, by = "Method") %>%
    dplyr::mutate(
      total = tidyr::replace_na(total, 0),
      dup   = tidyr::replace_na(dup, 0),
      value = dplyr::if_else(total > 0, dup / total, 0),
      panel = label
    ) %>%
    dplyr::right_join(all_methods, by = "Method") %>%   # ensure all methods included
    dplyr::mutate(
      value  = tidyr::replace_na(value, 0),
      panel  = label
    )
}

# Relative plot (single or faceted) saver
make_relative_plot <- function(df_both, title, ylab, outfile, method_colors, out_dir, facet_cols = NULL) {
  df_both <- df_both %>%
    dplyr::mutate(Method = factor(Method, levels = names(method_colors)))

  p <- ggplot(df_both, aes(x = Method, y = value, fill = Method)) +
    geom_col(colour = "black", width = 0.7) +
    geom_text(aes(label = scales::percent(value, accuracy = 0.1)), vjust = -0.3, size = 3) +
    scale_fill_manual(values = method_colors, drop = FALSE) +
    labs(x = NULL, y = ylab, fill = "Method", title = title) +
    theme_minimal() +
    theme(legend.position = "bottom")

  if (!is.null(facet_cols)) {
    p <- p + facet_wrap(facet_cols, ncol = 2)
  }

  ggsave(file.path(out_dir, outfile), p, width = 10, height = 4, dpi = 600, bg = "white")
}

# Metrics plot (Precision/Recall/F1)
make_metrics_plot <- function(dup_summary_run, method_colors, unique_id, out_dir) {
  all_methods <- names(method_colors)

  df_metrics <- dup_summary_run %>%
    dplyr::select(tool, precision, recall, f1) %>%
    tidyr::pivot_longer(cols = c(precision, recall, f1),
                        names_to = "Metric", values_to = "Value") %>%
    dplyr::mutate(
      Method = tool,
      Metric = recode(Metric, precision = "Precision", recall = "Recall", f1 = "F1")
    ) %>%
    tidyr::complete(Method = all_methods, Metric, fill = list(Value = 0)) %>%
    dplyr::mutate(
      Method = factor(Method, levels = all_methods),
      Metric = factor(Metric, levels = c("Precision", "Recall", "F1"))
    )

  p <- ggplot(df_metrics, aes(x = Method, y = Value, fill = Method)) +
    geom_col(position = position_dodge(width = 0.8), colour = "black", width = 0.7) +
    geom_text(aes(label = scales::percent(Value, accuracy = 0.1)),
              position = position_dodge(width = 0.8),
              vjust = -0.3, size = 3) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
    scale_fill_manual(values = method_colors, drop = FALSE) +
    labs(x = NULL, y = "Performance", fill = "Method",
         title = paste("Precision, Recall, and F1 -", unique_id)) +
    facet_wrap(~Metric, nrow = 1, scales = "free_y") +
    theme_minimal() +
    theme(legend.position = "bottom")

  ggsave(file.path(out_dir, paste0(unique_id, "_precision_recall_f1.png")),
         p, width = 10, height = 4, dpi = 600, bg = "white")
}

# Build binned counts (for histograms)
make_bin_counts <- function(df, id_col, level, binwidth, d_cutoff, classify_fun) {
  bin_left <- function(x, bw = binwidth, a = d_cutoff) a + bw * floor((x - a) / bw)

  df %>%
    dplyr::filter(Level == level) %>%
    dplyr::mutate(
      BinLeft   = round(bin_left(D_value), 10),
      BinCenter = BinLeft + binwidth/2
    ) %>%
    dplyr::distinct(Level, BinLeft, BinCenter, Method, {{ id_col }}) %>%
    dplyr::group_by(Level, BinLeft, BinCenter, {{ id_col }}) %>%
    dplyr::summarise(MethodsPresent = list(sort(unique(Method))), .groups = "drop") %>%
    dplyr::mutate(OverlapCategory = vapply(MethodsPresent, classify_fun, character(1))) %>%
    dplyr::group_by(Level, BinLeft, BinCenter, OverlapCategory) %>%
    dplyr::summarise(Count = n(), .groups = "drop")
}

# Combined histogram + boxplot saver for a given dataset
save_hist_box_combo <- function(df_all_like,
                                overlap_levels,
                                overlap_palette,
                                method_palette,
                                d_cutoff,
                                binwidth,
                                run_id,
                                out_dir,
                                suffix = "unfiltered") {
  # Prepare IDs for binning
  df_ids_gene <- df_all_like %>%
    dplyr::filter(Level == "Gene") %>%
    dplyr::mutate(GeneID = toupper(trimws(Gene)))

  df_ids_pair <- df_all_like %>%
    dplyr::filter(Level == "Pair") %>%
    dplyr::mutate(
      Gene_1 = toupper(trimws(Gene_1)),
      Gene_2 = toupper(trimws(Gene_2)),
      PairID = paste(pmin(Gene_1, Gene_2), pmax(Gene_1, Gene_2), sep = "-")
    )

  # Binned counts
  bin_counts_gene <- make_bin_counts(df_ids_gene, GeneID, "Gene", binwidth, d_cutoff, classify_overlap)
  bin_counts_pair <- make_bin_counts(df_ids_pair, PairID, "Pair", binwidth, d_cutoff, classify_overlap)

  # Histograms
  p_hist_genes <- ggplot(
    bin_counts_gene %>% dplyr::mutate(OverlapCategory = factor(OverlapCategory, levels = overlap_levels)),
    aes(x = BinCenter, y = Count, fill = OverlapCategory)
  ) +
    geom_col(width = binwidth, position = "stack", color = NA, alpha = 0.85) +
    geom_vline(aes(xintercept = d_cutoff, linetype = "D-cutoff"), color = "black") +
    facet_wrap(~Level, scales = "free_y") +
    scale_fill_manual(name = "Overlap category", values = overlap_palette, drop = FALSE) +
    scale_linetype_manual(name = "", values = c("D-cutoff" = "dashed")) +
    labs(title = paste("Genes per D-value bin -", run_id),
         subtitle = sprintf("D-value cutoff: %.2f", d_cutoff),
         x = "D-value", y = "Genes") +
    theme_minimal() +
    guides(fill = "none")

  p_hist_pairs <- ggplot(
    bin_counts_pair %>% dplyr::mutate(OverlapCategory = factor(OverlapCategory, levels = overlap_levels)),
    aes(x = BinCenter, y = Count, fill = OverlapCategory)
  ) +
    geom_col(width = binwidth, position = "stack", color = NA, alpha = 0.85) +
    geom_vline(aes(xintercept = d_cutoff, linetype = "D-cutoff"), color = "black") +
    facet_wrap(~Level, scales = "free_y") +
    scale_fill_manual(
      name   = "Overlap category",
      values = overlap_palette,
      drop   = FALSE,
      guide  = guide_legend(
        order = 1,
        override.aes = list(color = "black", size = 0.8, alpha = 1),
        keywidth  = 1.2,
        keyheight = 0.8
      )
    ) +
    scale_linetype_manual(
      name   = "",
      values = c("D-cutoff" = "dashed"),
      guide  = guide_legend(order = 2)
    ) +
    labs(title = paste("Pairs per D-value bin -", run_id),
         subtitle = sprintf("D-value cutoff: %.2f", d_cutoff),
         x = "D-value", y = "Gene pairs") +
    theme_minimal() +
    theme(legend.position = "bottom")

  # Boxplots
  p_box_genes <- ggplot(df_all_like %>% dplyr::filter(!is.na(Gene), Level == "Gene"),
                        aes(x = Method, y = D_value, fill = Method)) +
    geom_boxplot(outlier.alpha = 0.3) +
    geom_hline(aes(yintercept = d_cutoff, linetype = "D-cutoff"), color = "black") +
    scale_linetype_manual(name = "", values = c("D-cutoff" = "dashed")) +
    scale_fill_manual(values = method_palette) +
    labs(title = "D-value distributions (genes)", x = "Method", y = "D-value") +
    theme_minimal() +
    theme(legend.position = "none")

  p_box_pairs <- ggplot(df_all_like %>% dplyr::filter(!is.na(Gene_1) & !is.na(Gene_2), Level == "Pair"),
                        aes(x = Method, y = D_value, fill = Method)) +
    geom_boxplot(outlier.alpha = 0.3) +
    geom_hline(aes(yintercept = d_cutoff, linetype = "D-cutoff"), color = "black") +
    scale_linetype_manual(name = "", values = c("D-cutoff" = "dashed")) +
    scale_fill_manual(values = method_palette) +
    labs(title = "D-value distributions (gene pairs)", x = "Method", y = "D-value") +
    theme_minimal() +
    theme(legend.position = "none")

  # Combine plots
  top_row <- (p_hist_genes | p_hist_pairs) +
    plot_layout(guides = "collect") &
    theme(legend.position = "bottom")

  bottom_row <- (p_box_genes | p_box_pairs) +
    plot_layout(guides = "collect") &
    theme(legend.position = "bottom")

  combined_plot <- top_row / bottom_row + plot_layout(heights = c(2, 1))

  ggsave(file.path(out_dir, paste0("combined_hist_box_", suffix, "_", run_id, ".png")),
         plot = combined_plot, width = 12, height = 8, dpi = 300, bg = "white")
}

# D-value binning function for F1/precision/recall over D-value plots
compute_bin_metrics <- function(pairs_df, binwidth = 0.25, d_cutoff = 0, total_dups = NULL) {
  bin_left <- function(x, bw = binwidth, a = d_cutoff) a + bw * floor((x - a) / bw)

  # Bin pairs
  df_bins <- pairs_df %>%
    mutate(
      BinLeft   = round(bin_left(D_value), 10),
      BinCenter = BinLeft + binwidth/2,
      is_true   = PairType == "Correct_dup_pair"
    )

  # Count per bin per method
  metrics <- df_bins %>%
    group_by(Method, BinCenter) %>%
    summarise(
      total_pairs = n(),
      found       = sum(is_true),
      .groups = "drop"
    ) %>%
    rowwise() %>%
    mutate(
      TP = found,
      FN = if (!is.null(total_dups)) max(total_dups - found, 0) else 0,
      FP = max(total_pairs - found, 0),
      precision = ifelse(TP + FP > 0, TP / (TP + FP), 0),
      recall    = ifelse(TP + FN > 0, TP / (TP + FN), 0),
      f1        = ifelse(precision + recall > 0, 2 * precision * recall / (precision + recall), 0)
    ) %>%
    ungroup()

  metrics
}

# F1/precision/recall over D-value plotting function 
plot_bin_metrics <- function(metrics, method_colors, unique_id, out_dir) {
  df_long <- metrics %>%
    select(Method, BinCenter, precision, recall, f1) %>%
    pivot_longer(cols = c(precision, recall, f1),
                 names_to = "Metric", values_to = "Value") %>%
    mutate(
      Metric = recode(Metric,
        precision = "Precision",
        recall    = "Recall",
        f1        = "F1"
      )
    )

  p <- ggplot(df_long, aes(x = BinCenter, y = Value, color = Method)) +
    geom_line(linewidth = 1) +
    geom_point(size = 1.5) +
    scale_color_manual(values = method_colors) +
    labs(
      title = paste("Precision, Recall, and F1 across D-value -", unique_id),
      x = "D-value bin center",
      y = "Score",
      color = "Method"
    ) +
    facet_wrap(~Metric, nrow = 1, scales = "free_y") +
    theme_minimal() +
    theme(legend.position = "bottom")

  ggsave(file.path(out_dir, paste0(unique_id, "_metrics_over_bins.png")),
         p, width = 10, height = 4, dpi = 600, bg = "white")
}

# -------------------------
# Parse arguments
# -------------------------
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 4) {
  stop("Usage: Rscript compare_dvalues.R
        /path/to/{run_id}/coinfinder_dvalues.csv
        /path/to/{run_id}/goldfinder_dvalues.csv
        /path/to/{run_id}/panforest_dvalues.csv
        /path/to/coinfinder/{run_id}_d_cutoff.txt
        [/path/to/simulated/dup_match_summary.tsv]
        /path/to/out_dir")
}

coin_path       <- args[1]
gold_path       <- args[2]
pan_path        <- args[3]
dcutoff_file    <- args[4]
dcutoff_value   <- as.numeric(readLines(dcutoff_file, n = 1))
dup_summary_path<- if (length(args) > 5) args[5] else NULL
out_dir         <- args[length(args)]
unique_id       <- stringr::str_extract(basename(coin_path), "\\d+(?=\\.csv$)")

# -------------------------
# Data loading
# -------------------------
coin <- safe_read(coin_path, "Coinfinder")
gold <- safe_read(gold_path, "Goldfinder")
pan  <- safe_read(pan_path,  "PanForest")

df_all <- dplyr::bind_rows(coin, gold, pan)

# Goldfinder-filtered master table
df_all_filtered <- df_all %>%
  dplyr::filter(Method != "Goldfinder" | D_value >= dcutoff_value)

# Pairs table enriched
pairs_df <- prepare_pairs_df(df_all)

# -------------------------
# Style constants
# -------------------------
method_colors <- c(
  Coinfinder = "#D95F5F",
  Goldfinder = "#FFD966",
  PanForest  = "#5F9ED1"
)

overlap_levels <- c(
  "Coinfinder only",
  "Goldfinder only",
  "PanForest only",
  "Coinfinder + Goldfinder",
  "Coinfinder + PanForest",
  "Goldfinder + PanForest",
  "All three methods"
)

overlap_palette <- c(
  "Coinfinder only"         = "#D95F5F",
  "Goldfinder only"         = "#FFD966",
  "PanForest only"          = "#5F9ED1",
  "Coinfinder + Goldfinder" = "#F2A65A",
  "Coinfinder + PanForest"  = "#B07CC6",
  "Goldfinder + PanForest"  = "#7FB77E",
  "All three methods"       = "#8C564B"
)

method_palette <- c(
  Coinfinder = unname(overlap_palette["Coinfinder only"]),
  Goldfinder = unname(overlap_palette["Goldfinder only"]),
  PanForest  = unname(overlap_palette["PanForest only"])
)

set_order <- c("Coinfinder", "Goldfinder", "PanForest")
binwidth  <- 0.25
d_cutoff  <- dcutoff_value
run_id    <- unique_id

# -------------------------
# UpSet plots
# -------------------------

# Gene-level (all genes)
df_gene <- make_upset_df(df_all %>% dplyr::filter(Level == "Gene"),
                         id_col = Gene, filter_expr = TRUE, method_order = set_order)
make_upset_plot(df_gene,
                paste("Gene counts across the methods -", unique_id),
                "Gene count",
                paste0(unique_id, "_upset.png"),
                method_colors, out_dir, show_pattern_legend = TRUE)

# Gene-level (dup genes only), skip if none
df_dup_gene <- make_upset_df(df_all %>% dplyr::filter(Level == "Gene"),
                             id_col = Gene, filter_expr = grepl("_dup$", Gene), method_order = set_order)
if (nrow(df_dup_gene) > 0) {
  make_upset_plot(df_dup_gene,
                  paste("Gene counts across the methods (duplicated genes only) -", unique_id),
                  "Gene count",
                  paste0(unique_id, "_upset_dup.png"),
                  method_colors, out_dir, show_pattern_legend = FALSE)
} else {
  message("No _dup genes found across any methods - skipping _dup UpSet plot.")
}

# Pair-level (all pairs)
df_pairs_upset <- make_upset_df(pairs_df, id_col = PairID, filter_expr = TRUE, method_order = set_order)
make_upset_plot(df_pairs_upset,
                paste("Pair counts across the methods -", unique_id),
                "Pair count",
                paste0(unique_id, "_upset_pairs.png"),
                method_colors, out_dir, show_pattern_legend = TRUE)

# Pair-level (dup pairs only), skip if none
df_dup_pairs_upset <- make_upset_df(pairs_df, id_col = PairID,
                                    filter_expr = PairType == "Correct_dup_pair",
                                    method_order = set_order)
if (nrow(df_dup_pairs_upset) > 0) {
  make_upset_plot(df_dup_pairs_upset,
                  paste("Pair counts across the methods (duplicated genes only) -", unique_id),
                  "Pair count",
                  paste0(unique_id, "_upset_dup_pairs.png"),
                  method_colors, out_dir, show_pattern_legend = FALSE)
} else {
  message("No _dup pairs found across any methods - skipping _dup UpSet plot.")
}

# -------------------------
# Relative counts panels (genes and pairs)
# -------------------------
if (!is.null(dup_summary_path) && file.exists(dup_summary_path)) {
  dup_summary <- readr::read_tsv(dup_summary_path, show_col_types = FALSE)

  dup_summary_run <- dup_summary %>%
    dplyr::filter(pangenome_id == unique_id) %>%
    dplyr::mutate(
      tool = dplyr::recode(tool,
        "coinfinder" = "Coinfinder",
        "goldfinder" = "Goldfinder",
        "panforest"  = "PanForest"
      )
    )

  if (nrow(dup_summary_run) > 0) {
    # Ensure all three tools present for summary-driven panels
    all_tools_tbl <- tibble::tibble(tool = names(method_colors))
    dup_summary_run <- all_tools_tbl %>%
      dplyr::left_join(dup_summary_run, by = "tool") %>%
      dplyr::mutate(
        duplicate_found_pct = tidyr::replace_na(duplicate_found_pct, 0),
        tool = factor(tool, levels = names(method_colors))
      )

    # Panel 1: from summary file (duplicate_found_pct)
    df_panel1 <- dup_summary_run %>%
      dplyr::transmute(
        Method = tool,
        value  = duplicate_found_pct / 100,
        panel  = "Relative duplicate gene counts vs all duplicated genes"
      )

    # Panel 2: computed from df_all (dup / total genes)
    df_genes_all <- df_all %>% dplyr::filter(Level == "Gene")

    df_panel2 <- make_relative_counts(df_genes_all, Gene,
                                      "Relative duplicate gene counts vs all genes",
                                      method_colors) %>%
      dplyr::mutate(Method = dplyr::recode(Method,
        "coinfinder" = "Coinfinder",
        "goldfinder" = "Goldfinder",
        "panforest"  = "PanForest"
      )) %>%
      dplyr::mutate(Method = factor(Method, levels = names(method_colors)))

    df_both <- dplyr::bind_rows(df_panel1, df_panel2) %>%
      dplyr::mutate(
        panel = factor(panel, levels = c(
          "Relative duplicate gene counts vs all duplicated genes",
          "Relative duplicate gene counts vs all genes"
        ))
      )

    make_relative_plot(df_both,
                       paste("Relative duplicate gene counts -", unique_id),
                       "Proportion of duplicate genes",
                       paste0(unique_id, "_dup_relative.png"),
                       method_colors, out_dir,
                       facet_cols = ~panel)

    # Pair-level relative panel from summary (dup_as_pct_of_pairs)
    df_panel_pairs <- dup_summary_run %>%
      dplyr::transmute(
        Method = tool,
        value  = dup_as_pct_of_pairs / 100,
        panel  = "Relative duplicate pair counts vs all pairs"
      ) %>%
      tidyr::complete(Method = names(method_colors),
                      fill = list(value = 0, panel = "Relative duplicate pair counts vs all pairs")) %>%
      dplyr::mutate(
        Method = factor(Method, levels = names(method_colors)),
        panel  = factor(panel, levels = c("Relative duplicate pair counts vs all pairs"))
      )

    make_relative_plot(df_panel_pairs,
                       paste("Relative duplicate pair counts -", unique_id),
                       "Proportion of duplicate pairs",
                       paste0(unique_id, "_dup_relative_pairs.png"),
                       method_colors, out_dir,
                       facet_cols = NULL)

    # Precision/Recall/F1
    make_metrics_plot(dup_summary_run, method_colors, unique_id, out_dir)

  } else {
    message("No matching pangenome_id found in summary file for ", unique_id)
  }
}

# Count all Correct_dup_pairs in pairs_df
total_dups <- sum(pairs_df$PairType == "Correct_dup_pair")

bin_metrics <- compute_bin_metrics(pairs_df,
                                   binwidth = 0.25,
                                   d_cutoff = dcutoff_value,
                                   total_dups = total_dups)

plot_bin_metrics(bin_metrics, method_colors, unique_id, out_dir)

# -------------------------
# Combined histograms + boxplots (unfiltered and filtered)
# -------------------------

# Unfiltered
save_hist_box_combo(
  df_all_like       = df_all,
  overlap_levels    = overlap_levels,
  overlap_palette   = overlap_palette,
  method_palette    = method_palette,
  d_cutoff          = d_cutoff,
  binwidth          = binwidth,
  run_id            = run_id,
  out_dir           = out_dir,
  suffix            = "unfiltered"
)

# Filtered (Goldfinder D-value cutoff applied)
save_hist_box_combo(
  df_all_like       = df_all_filtered,
  overlap_levels    = overlap_levels,
  overlap_palette   = overlap_palette,
  method_palette    = method_palette,
  d_cutoff          = d_cutoff,
  binwidth          = binwidth,
  run_id            = run_id,
  out_dir           = out_dir,
  suffix            = "filtered"
)