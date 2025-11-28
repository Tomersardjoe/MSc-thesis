#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(ggupset)
  library(ggpattern)
  library(ggpubr)
  library(patchwork)
  library(readr)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(lme4)
  library(lmerTest)
  library(emmeans)
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
message("Safe loading D-value CSVs - Do not worry if files are missing! This is most likely due to a tool not having produced significant results (always check). This is expected and does not break the script.")
safe_read <- function(path, method_name) {
  if (!is.null(path) && file.exists(path)) {
    # Check if file has more than just the header
    n_lines <- length(readr::read_lines(path))
    if (n_lines <= 1) {
      message("Skipping empty file (header only): ", path)
      return(tibble::tibble(
        Level  = character(),
        Gene   = character(),
        Gene_1 = character(),
        Gene_2 = character(),
        Pair   = character(),
        D_value= numeric(),
        Method = character()
      ))
    }
    
    readr::read_csv(
      path,
      show_col_types = FALSE,
      col_types = readr::cols(
        Level   = readr::col_character(),
        Gene    = readr::col_character(),
        Gene_1  = readr::col_character(),
        Gene_2  = readr::col_character(),
        D_value = readr::col_double()
      )
    ) %>%
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
make_upset_plot <- function(df2, title, ylab, outfile, method_colors, out_dir,
                            show_pattern_legend = FALSE, relative = FALSE) {
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
      pattern      = factor(pattern, levels = c("none","stripe","crosshatch")),
      sets         = lapply(sets, function(x) factor(x,
                            levels = c("Coinfinder","Goldfinder","PanForest")))
    )

  plot_df <- df2
  if (relative) {
    total <- sum(plot_df$value)
    plot_df <- plot_df %>%
      dplyr::mutate(value = value / total)
  }

  p <- ggplot(plot_df, aes(x = sets, y = value)) +
    geom_col(aes(fill = base_tool), colour = "black", linewidth = 0.2, width = 0.7,
             show.legend = c(fill = TRUE)) +
    geom_col(data = dplyr::filter(plot_df, n_sets == 3), aes(y = value),
             fill = "#999999", colour = NA, width = 0.7, show.legend = FALSE) +
    geom_col_pattern(
      aes(pattern = pattern, pattern_fill = overlay_tool),
      fill = NA, colour = NA, pattern_colour = "black",
      linewidth = 0.2, width = 0.7, pattern_size = 0.2,
      show.legend = c(pattern = show_pattern_legend, pattern_fill = FALSE, fill = FALSE),
      key_glyph = ggpattern::draw_key_polygon_pattern,
      pattern_angle = 45, pattern_density = 0.5, pattern_spacing = 0.05
    ) +
    geom_label(aes(y = value,
                   label = if (relative) scales::percent(value) else value),
               vjust = -0.3, size = 3, label.size = 0, fill = "white") +
    scale_y_continuous(
      expand = expansion(mult = c(0, 0.2)),
      labels = if (relative) scales::percent_format() else scales::label_number()
    ) +
    scale_x_upset(order_by = "degree", sets = c("Coinfinder","Goldfinder","PanForest")) +
    scale_fill_manual(name = "Method", values = method_colors, limits = names(method_colors),
                      breaks = names(method_colors), drop = FALSE) +
    scale_pattern_manual(values = c(none = "none", stripe = "stripe", crosshatch = "crosshatch"), guide = "none") +
    scale_pattern_fill_manual(values = pattern_fill_colors, na.value = NA, guide = "none") +
    labs(x = NULL, y = ylab, title = title) +
    theme_minimal() +
    theme_combmatrix(
      combmatrix.panel.point.size = 1.5,
      combmatrix.panel.line.size  = 0.3
    ) +
    theme(plot.title = element_text(size = 8))

  return(p)
}

make_metrics_data <- function(dup_summary_run) {
  dup_summary_run %>%
    group_by(Method, category) %>%
    summarise(
      precision = mean(precision, na.rm = TRUE),
      recall    = mean(recall, na.rm = TRUE),
      f1        = mean(f1, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    pivot_longer(
      cols = c(precision, recall, f1),
      names_to = "Metric", values_to = "Value"
    ) %>%
    mutate(
      Metric   = recode(Metric,
                        precision = "Precision",
                        recall    = "Recall",
                        f1        = "F1"),
      Metric   = factor(Metric, levels = c("Precision", "Recall", "F1")),
      category = factor(category, levels = c("Closed", "Moderate", "Open"))
    )
}

make_metrics_plot_grouped <- function(df_metrics, method_colors, sig_levels) {
  ymax <- max(df_metrics$Value, na.rm = TRUE)

  ggplot(df_metrics, aes(x = Metric, y = Value, fill = Method)) +
    geom_col(position = position_dodge(width = 0.8), colour = "black", linewidth = 0.2, width = 0.7) +
    geom_text(aes(label = sprintf("%.2f", Value)),
              position = position_dodge(width = 0.8),
              vjust = -0.3, size = 2.5) +
    geom_text(data = sig_levels,
              aes(x = Metric, y = Value, label = signif),
              inherit.aes = FALSE,
              vjust = -0.3, size = 2.5) +
    scale_y_continuous(limits = c(0, ymax + 0.15),
                       labels = scales::number_format(accuracy = 0.01),
                       breaks = scales::pretty_breaks(n = 4)) +
    scale_fill_manual(values = method_colors) +
    facet_wrap(~ category, ncol = 1) +
    labs(x = NULL, y = "Performance",
         title = "Precision, Recall, and F1 per tool across categories") +
    theme_minimal() +
    theme(legend.position = "none",
          plot.title = element_text(size = 8),
          strip.text = element_text(size = 8),
          axis.text.y = element_text(size = 5),
          axis.text.x = element_text(size = 8)
    )
}

p_to_star <- function(p) {
  if (p < 0.001) return("***")
  else if (p < 0.01) return("**")
  else if (p < 0.05) return("*")
  else return("ns")
}

# D-value binning function for F1/precision/recall over D-value plots
compute_bin_metrics_cat <- function(pairs_df, n_bins = 20) {
  
  # Define quantile cut points across observed D-values
  breaks <- quantile(pairs_df$D_value, probs = seq(0, 1, length.out = n_bins + 1), na.rm = TRUE)
  breaks <- unique(breaks)

  df_bins <- pairs_df %>%
    dplyr::mutate(
      Bin = cut(D_value, breaks = breaks, include.lowest = TRUE, labels = FALSE),
      BinCenter = (breaks[Bin] + breaks[Bin + 1]) / 2,
      is_true   = PairType == "Correct_dup_pair"
    )

  # Retrieve total duplicates per category
  total_dups_df <- pairs_df %>%
    dplyr::filter(PairType == "Correct_dup_pair") %>%
    dplyr::count(category, name = "total_dups")

  metrics <- df_bins %>%
    dplyr::group_by(category, Method, BinCenter) %>%
    dplyr::summarise(
      total_pairs = dplyr::n(),
      found       = sum(is_true),
      .groups = "drop"
    ) %>%
    dplyr::left_join(total_dups_df, by = "category") %>%
    dplyr::rowwise() %>%
    dplyr::mutate(
      TP = found,
      FN = max(total_dups - found, 0),
      FP = max(total_pairs - found, 0),
      precision = ifelse(TP + FP > 0, TP / (TP + FP), 0),
      recall    = ifelse(TP + FN > 0, TP / (TP + FN), 0),
      f1        = ifelse(precision + recall > 0,
                         2 * precision * recall / (precision + recall), 0)
    ) %>%
    dplyr::ungroup()

  metrics
}

# F1/precision/recall over D-value plotting function 
plot_bin_metrics_cat <- function(metrics, method_colors, metric_name) {
  df_long <- metrics %>%
    tidyr::pivot_longer(cols = c("precision", "recall", "f1"),
                        names_to = "Metric", values_to = "Value") %>%
    dplyr::mutate(
      Metric = dplyr::recode(Metric,
        precision = "Precision",
        recall    = "Recall",
        f1        = "F1"
      )
    ) %>%
    dplyr::filter(Metric == dplyr::recode(metric_name,
                                          precision = "Precision",
                                          recall    = "Recall",
                                          f1        = "F1"))

  p <- ggplot(df_long, aes(x = BinCenter, y = Value, color = Method)) +
    geom_line(linewidth = 1, alpha = 0.8) +
    geom_point(size = 0.5, alpha = 0.8) +
    scale_color_manual(values = method_colors) +
    labs(
      title = paste(unique(df_long$Metric), "across pangenome structure categories"),
      x = "D-value",
      y = unique(df_long$Metric),
      color = "Method"
    ) +
    facet_wrap(~category, nrow = 1, scales = "fixed") +
    theme_minimal() +
    theme(
      legend.position = "bottom",
      plot.title = element_text(size = 8)
    )

  return(p)
}

# -------------------------
# Style constants
# -------------------------
method_colors <- c(
  Coinfinder = "#D95F5F",
  Goldfinder = "#FFD966",
  PanForest  = "#5F9ED1"
)

prop_colors <- c(
  Closed = "#D44949",
  Moderate = "#EC9837",
  Open = "#58D168"
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

set_order <- c("Coinfinder", "Goldfinder", "PanForest")
category_order <- c("Closed", "Moderate", "Open")

# -------------------------
# Parse arguments
# -------------------------
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 5) {
  stop("Usage: Rscript aggregated_upset.R
        /path/to/{dataset}/coinfinder_runs(_all)
        /path/to/{dataset}/goldfinder_runs(_all)
        /path/to/{dataset}/panforest_runs(_all)
        /path/to/real_pangenomes/[species_categories or matched_all].csv
        /path/to/out_dir
        [optional: /path/to/summary_file]")
}

coin_dir        <- args[1]
gold_dir        <- args[2]
pan_dir         <- args[3]
categories_path <- args[4]
out_dir         <- args[5]

# Optional summary file (only for simulated datasets)
summary_file <- NULL
if (length(args) >= 6) {
  summary_file <- args[6]
  dup_summary_run <- readr::read_tsv(summary_file, show_col_types = FALSE)
}

# Detect simulation type if summary_file is provided
sim_type <- NULL
if (!is.null(summary_file)) {
  if (grepl("flip", summary_file)) {
    sim_type <- "flip"
  } else if (grepl("perfect", summary_file)) {
    sim_type <- "perfect"
  } else {
    sim_type <- "unknown"
  }
  message("Detected simulated dataset type: ", sim_type)
}

# -------------------------
# Data loading
# -------------------------

species_taxids <- list.dirs(coin_dir, recursive = FALSE, full.names = FALSE)

all_runs <- purrr::map_dfr(seq_along(species_taxids), function(i) {
  species_taxid <- species_taxids[i]

  tryCatch({
    coin_path <- file.path(coin_dir, species_taxid, "d_cutoff",
                           paste0("coinfinder_dvalues_", species_taxid, ".csv"))
    gold_path <- file.path(gold_dir, species_taxid, "d_distribution",
                           paste0("goldfinder_dvalues_", species_taxid, ".csv"))
    pan_path  <- file.path(pan_dir,  species_taxid, "imp_cutoff",
                           paste0("panforest_dvalues_", species_taxid, ".csv"))

    coin <- safe_read(coin_path, "Coinfinder")
    gold <- safe_read(gold_path, "Goldfinder")
    pan  <- safe_read(pan_path,  "PanForest")

    dplyr::bind_rows(coin, gold, pan) %>%
      dplyr::mutate(species_taxid = species_taxid)

  }, error = function(e) {
    message("Failed at index ", i, " (species_taxid = ", species_taxid, "): ", conditionMessage(e))
    tibble::tibble(
      Level  = character(),
      Gene   = character(),
      Gene_1 = character(),
      Gene_2 = character(),
      Pair   = character(),
      D_value= numeric(),
      Method = character(),
      species_taxid = character()
    )
  })
})

categories <- readr::read_csv(
  categories_path,
  col_names = c("category", "species_taxid"),
  col_select = 1:2,
  show_col_types = FALSE
)

all_runs <- all_runs %>%
  dplyr::left_join(categories, by = "species_taxid")
  
if ("category" %in% names(all_runs)) {
  all_runs <- all_runs %>% dplyr::filter(!is.na(category))
} else {
  message("No 'category' column found after join - skipping filter.")
}    

# -------------------------
# UpSet plots
# -------------------------

for (cat in unique(all_runs$category)) {
  df_cat <- all_runs %>% dplyr::filter(category == cat)

  # Gene-level
  df_gene <- make_upset_df(df_cat %>% dplyr::filter(Level == "Gene"),
                           id_col = Gene, filter_expr = TRUE, method_order = set_order)
  p_gene <- make_upset_plot(df_gene,
                            paste("Gene counts -", cat, "pangenomes"),
                            NULL,
                            paste0(cat, "_upset_genes.png"),
                            method_colors, out_dir)

  assign(paste0("p_gene_", cat), p_gene)

  # Pair-level
  pairs_df <- prepare_pairs_df(df_cat)
  df_pairs <- make_upset_df(pairs_df, id_col = PairID, filter_expr = TRUE, method_order = set_order)
  p_pair <- make_upset_plot(df_pairs,
                            paste("Pair counts -", cat, "pangenomes"),
                            NULL,
                            paste0(cat, "_upset_pairs.png"),
                            method_colors, out_dir)

  assign(paste0("p_pair_", cat), p_pair)
}

# Proportion of all three tool agreement to total gene pool of each category
prop_df <- all_runs %>%
  filter(Level == "Gene") %>%
  group_by(category, Gene) %>%
  summarise(methods = list(unique(Method)), .groups = "drop") %>%
  mutate(overlap = purrr::map_lgl(methods, ~ setequal(.x, c("Coinfinder","Goldfinder","PanForest")))) %>%
  group_by(category) %>%
  summarise(prop_all_three = mean(overlap))

# Statistics

# Build contingency table
gene_cont_table <- all_runs %>%
  filter(Level == "Gene") %>%
  group_by(category, Gene) %>%
  summarise(methods = list(unique(Method)), .groups = "drop") %>%
  mutate(overlap = purrr::map_lgl(methods,
                                  ~ setequal(.x, c("Coinfinder","Goldfinder","PanForest")))) %>%
  group_by(category) %>%
  summarise(
    n_all_three_genes = sum(overlap),
    n_total_genes     = n(),
    .groups = "drop"
  ) %>%
  mutate(n_not_all_three = n_total_genes - n_all_three_genes) %>%
  select(category, n_all_three_genes, n_not_all_three)

# Convert to matrix
gene_cont_table <- as.matrix(gene_cont_table[, -1])
rownames(gene_cont_table) <- all_runs %>%
  filter(Level == "Gene") %>%
  pull(category) %>%
  unique()
  
# Global chi-squared test
chisq_gene <-  suppressWarnings(chisq.test(gene_cont_table))

# Check expected counts
if (any(chisq_gene$expected < 5)) {
  message("Chi-squared approximation may be incorrect; using simulated p-values.")
  chisq_gene <- chisq.test(gene_cont_table, simulate.p.value = TRUE, B = 10000)
}

cat(sprintf(
  "The shared gene proportions%s differ significantly (p %s)\n",
  ifelse(chisq_gene$p.value < 0.05, "", "do not"),
  ifelse(chisq_gene$p.value < 0.0001,
        "< 0.0001",
        sprintf("- %.4f", chisq_gene$p.value))
))

gene_overall_mean <- sum(gene_cont_table[, "n_all_three_genes"]) /
                     sum(rowSums(gene_cont_table))

# Category prop test vs overall mean
gene_category_tests <- data.frame(
  category = rownames(gene_cont_table),
  n_all_three_genes = gene_cont_table[, "n_all_three_genes"],
  n_total_genes = rowSums(gene_cont_table)
) %>%
  rowwise() %>%
  mutate(
    p_value = {
      # Run chi-squared proportion test first
      chisq_res <- suppressWarnings(
        prop.test(
          x = n_all_three_genes,
          n = n_total_genes,
          p = gene_overall_mean
        )
      )
      
      # Check expected counts
      if (any(chisq_res$expected < 5)) {
        # Fall back to exact binomial test
        binom.test(
          x = n_all_three_genes,
          n = n_total_genes,
          p = gene_overall_mean
        )$p.value
      } else {
        chisq_res$p.value
      }
    }
  ) %>%
  ungroup() %>%
  mutate(
    p_adj = p.adjust(p_value, method = "BH"),
    signif = case_when(
      p_adj < 0.001 ~ "***",
      p_adj < 0.01  ~ "**",
      p_adj < 0.05  ~ "*",
      TRUE          ~ "ns"
    )
  )

# Pairwise Fisher's exact test
gene_categories <- category_order[category_order %in% rownames(gene_cont_table)]
gene_results <- data.frame(
  category1 = character(),
  category2 = character(),
  p_value   = numeric(),
  stringsAsFactors = FALSE
)

for (i in 1:(length(gene_categories)-1)) {
  for (j in (i+1):length(gene_categories)) {
    sub_table <- gene_cont_table[c(i,j), ]
    test <- fisher.test(sub_table)
    gene_results <- rbind(gene_results, data.frame(
      category1 = gene_categories[i],
      category2 = gene_categories[j],
      p_value   = test$p.value
    ))
  }
}

gene_results <- gene_results %>%
  mutate(
    p_adj = p.adjust(p_value, method = "BH"),
    signif = case_when(
      p_adj < 0.001 ~ "***",
      p_adj < 0.01  ~ "**",
      p_adj < 0.05  ~ "*",
      TRUE          ~ "ns"
    ),
    p_value = format(p_value, format = "g", digits = 3),
    p_adj   = formatC(p_adj,   format = "g", digits = 3)
  )

# Save as TSV
write.table(gene_results,
            file = file.path(out_dir, "gene_pairwise_fisher.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

# Gene plot
p_gene_prop <- ggplot(prop_df, aes(x = category, y = prop_all_three)) +
  geom_col(fill = prop_colors[prop_df$category],
           colour = "black", linewidth = 0.2, width = 0.7) +
  geom_label(aes(label = scales::percent(prop_all_three, accuracy = 0.1)),
             vjust = -0.3, size = 3, label.size = 0, fill = "white") +
  geom_text(data = gene_category_tests,
            aes(x = category, y = (n_all_three_genes / n_total_genes) + 0.05, label = signif),
            inherit.aes = FALSE, size = 3) +
  scale_x_discrete(labels = stringr::str_to_title) +
  scale_y_continuous(labels = scales::percent_format(),
                     expand = expansion(mult = c(0, 0.2))) +
  labs(title = "Shared gene proportion of its category gene pool", y = NULL, x = NULL) +
  theme_minimal() +
  theme(plot.title = element_text(size = 8))

# Combine into a facetted figure
combined <- (p_gene_Closed | p_gene_Moderate) /
            (p_gene_Open   | p_gene_prop)

final_plot <- combined + plot_layout(guides = "collect", heights = c(1, 1))
  
final_plot <- final_plot + plot_annotation(tag_levels = "A")

ggsave(file.path(out_dir, "gene_category.png"), final_plot,
       width = 8, height = 5, dpi = 600, bg = "white")

# Save as PDF
ggsave(file.path(out_dir, "gene_category.pdf"), final_plot,
       width = 8, height = 5, bg = "white")

# Proportion of all three tool agreement to total pair pool of each category
prop_pair_df <- all_runs %>%
  filter(Level == "Pair") %>%
  mutate(
    Pair_norm = purrr::map_chr(strsplit(Pair, "_"),
                               ~ paste(sort(.x), collapse = "_"))
  ) %>%
  group_by(category, Pair_norm) %>%
  summarise(methods = list(unique(Method)), .groups = "drop") %>%
  mutate(overlap = purrr::map_lgl(
    methods,
    ~ all(c("Coinfinder","Goldfinder","PanForest") %in% .x)
  )) %>%
  group_by(category) %>%
  summarise(
    prop_all_three_pairs = mean(overlap),
    n_all_three_pairs    = sum(overlap),
    n_total_pairs        = n()
  )

# Statistics

# Build contingency table
pair_cont_table <- prop_pair_df %>%
  mutate(n_not_all_three = n_total_pairs - n_all_three_pairs) %>%
  select(category, n_all_three_pairs, n_not_all_three)

pair_cont_table <- as.matrix(pair_cont_table[, -1])
rownames(pair_cont_table) <- prop_pair_df$category

# Global chi-squared test
chisq_pair <-  suppressWarnings(chisq.test(pair_cont_table))

# Check expected counts
if (any(chisq_pair$expected < 5)) {
  message("Chi-squared approximation may be incorrect; using simulated p-values.")
  chisq_pair <- chisq.test(pair_cont_table, simulate.p.value = TRUE, B = 10000)
}

cat(sprintf(
  "The shared pair proportions%s differ significantly (p %s)\n",
  ifelse(chisq_pair$p.value < 0.05, "", "do not"),
  ifelse(chisq_pair$p.value < 0.0001,
         "< 0.0001",
         sprintf("= %.4f", chisq_pair$p.value))
))

# Category prop test vs overall mean
overall_mean <- sum(prop_pair_df$n_all_three_pairs) / sum(prop_pair_df$n_total_pairs)

category_tests <- prop_pair_df %>%
  rowwise() %>%
  mutate(
    p_value = {
      # Run chi-squared proportion test
      chisq_res <- suppressWarnings(
        prop.test(
          x = n_all_three_pairs,
          n = n_total_pairs,
          p = overall_mean
        )
      )
      
      # Check expected counts
      if (any(chisq_res$expected < 5)) {
        # Fall back to exact binomial test
        binom.test(
          x = n_all_three_pairs,
          n = n_total_pairs,
          p = overall_mean
        )$p.value
      } else {
        chisq_res$p.value
      }
    }
  ) %>%
  ungroup() %>%
  mutate(
    p_adj = p.adjust(p_value, method = "BH"),
    signif = case_when(
      p_adj < 0.001 ~ "***",
      p_adj < 0.01  ~ "**",
      p_adj < 0.05  ~ "*",
      TRUE          ~ "ns"
    )
  )

# Pairwise Fisher's exact test
categories <- category_order[category_order %in% rownames(pair_cont_table)]
results <- data.frame(
  category1 = character(),
  category2 = character(),
  p_value   = numeric(),
  stringsAsFactors = FALSE
)

for (i in 1:(length(categories)-1)) {
  for (j in (i+1):length(categories)) {
    sub_table <- pair_cont_table[c(i,j), ]
    test <- fisher.test(sub_table)
    results <- rbind(results, data.frame(
      category1 = categories[i],
      category2 = categories[j],
      p_value   = test$p.value
    ))
  }
}

results <- results %>%
  mutate(
    p_adj = p.adjust(p_value, method = "BH"),
    signif = case_when(
      p_adj < 0.001 ~ "***",
      p_adj < 0.01  ~ "**",
      p_adj < 0.05  ~ "*",
      TRUE          ~ "ns"
    ),
    p_value = formatC(p_value, format = "g", digits = 3),
    p_adj   = formatC(p_adj,   format = "g", digits = 3)
  )

write.table(results,
            file = file.path(out_dir, "pair_pairwise_fisher.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

# Pair plot
p_pair_prop <- ggplot(prop_pair_df, aes(x = category, y = prop_all_three_pairs)) +
  geom_col(fill = prop_colors[prop_pair_df$category],
           colour = "black", linewidth = 0.2, width = 0.7) +
  geom_label(aes(label = scales::percent(prop_all_three_pairs, accuracy = 0.001)),
             vjust = -0.3, size = 3, label.size = 0, fill = "white") +
  geom_text(data = category_tests,
            aes(x = category, y = prop_all_three_pairs + 0.01, label = signif),
            inherit.aes = FALSE, size = 3) +
  scale_x_discrete(labels = stringr::str_to_title) +
  scale_y_continuous(labels = scales::percent_format(),
                     expand = expansion(mult = c(0, 0.2))) +
  labs(title = "Shared pair proportion of its category pair pool", y = NULL, x = NULL) +
  theme_minimal() +
  theme(plot.title = element_text(size = 8))

combined <- (p_pair_Closed | p_pair_Moderate) /
            (p_pair_Open   | p_pair_prop)

final_plot <- combined + plot_layout(guides = "collect", heights = c(1, 1)) +
              plot_annotation(tag_levels = "A")

# Save as PNG
ggsave(file.path(out_dir, "pair_category.png"), final_plot,
       width = 8, height = 5, dpi = 600, bg = "white")

# Save as PDF
ggsave(file.path(out_dir, "pair_category.pdf"), final_plot,
       width = 8, height = 5, bg = "white")

##### PSEUDO-SIMULATED DATASETS ONLY #####

# -------------------------
# UpSet plots (duplicated genes)
# -------------------------

if (!is.null(summary_file)) {

  for (cat in unique(all_runs$category)) {
    df_cat <- all_runs %>% dplyr::filter(category == cat)
  
    # Gene-level (dup genes only)
    df_dup_gene <- make_upset_df(
      df_cat %>% dplyr::filter(Level == "Gene"),
      id_col = Gene,
      filter_expr = grepl("_dup$", Gene),
      method_order = set_order
    )
    if (nrow(df_dup_gene) > 0) {
      p_gene_dup <- make_upset_plot(
        df_dup_gene,
        paste("Duplicated gene counts -", cat, "pangenomes"),
        "Gene count",
        paste0(cat, "_upset_genes_dup.png"),
        method_colors, out_dir, show_pattern_legend = FALSE
      )
      assign(paste0("p_gene_dup_", cat), p_gene_dup)
    } else {
      message("No _dup genes found in category ", cat, " - skipping _dup UpSet plot.")
    }
  
    # Pair-level (dup pairs only)
    pairs_df <- prepare_pairs_df(df_cat)
    
    df_dup_pairs_upset <- make_upset_df(
      pairs_df, id_col = PairID,
      filter_expr = PairType == "Correct_dup_pair",
      method_order = set_order
    )
    
    if (nrow(df_dup_pairs_upset) > 0) {
      p_pair_dup <- make_upset_plot(
        df_dup_pairs_upset,
        paste("Proportion of duplicated pairs found -", cat, "pangenomes"),
        "Pair proportion",
        paste0(cat, "_upset_pairs_dup.png"),
        method_colors, out_dir,
        show_pattern_legend = FALSE,
        relative = TRUE
      )
      assign(paste0("p_pair_dup_", cat), p_pair_dup)
    }
  }
  
  # Statistics
  
  df_metrics <- make_metrics_data(dup_summary_run)
  
  # Linear mixed models for each metric
  metrics <- c("Precision", "Recall", "F1")
  
  lmm_results <- lapply(metrics, function(m) {
    df_sub <- df_metrics %>% filter(Metric == m)
  
    model <- lmer(Value ~ category + (1|Method), data = df_sub)
    anova_res <- anova(model)
  
    list(metric = m, model = model, anova = anova_res)
  })

  # ANOVA p-values for plotting
  pvals_df <- do.call(rbind, lapply(lmm_results, function(res) {
    row <- which(rownames(res$anova) == "category")
    if (length(row) == 1) {
      pval <- res$anova[row, "Pr(>F)"]
      star <- p_to_star(pval)

      data.frame(
        Metric  = res$metric,
        pval    = pval,
        signif  = star,
        stringsAsFactors = FALSE
      )
    } else {
      data.frame(Metric = res$metric, pval = NA, signif = "", stringsAsFactors = FALSE)
    }
  }))

  sig_levels <- data.frame(
    Metric   = pvals_df$Metric,
    Value    = rep(max(df_metrics$Value) + 0.005, length(pvals_df)),
    signif   = pvals_df$signif,
    category = "Closed"
  )
  
  # Post-hoc contrasts only if ANOVA significant
  posthoc_df <- do.call(rbind, lapply(lmm_results, function(res) {
    row <- which(rownames(res$anova) == "category")
    if (length(row) == 1) {
      pval <- res$anova[row, "Pr(>F)"]
      star <- p_to_star(pval)
  
      # Run post-hoc only if significant
      if (!is.na(pval) && pval < 0.05) {
        contr <- summary(emmeans(res$model, pairwise ~ category)$contrasts)
        data.frame(
          Metric   = res$metric,
          ANOVA_p  = formatC(pval, format = "g", digits = 3),
          ANOVA_sig= star,
          Contrast = contr$contrast,
          Estimate = contr$estimate,
          SE       = contr$SE,
          df       = contr$df,
          t.ratio  = contr$t.ratio,
          p.value  = contr$p.value,
          signif   = sapply(contr$p.value, p_to_star),
          stringsAsFactors = FALSE
        )
      } else {
        NULL
      }
    } else {
      NULL
    }
  }))
  
  if (!is.null(posthoc_df)) {
    posthoc_df <- posthoc_df %>%
      mutate(
        p_adj = p.adjust(p.value, method = "BH"),
        signif = case_when(
          p_adj < 0.001 ~ "***",
          p_adj < 0.01  ~ "**",
          p_adj < 0.05  ~ "*",
          TRUE          ~ "ns"
        ),
        p.value = formatC(p.value, format = "g", digits = 3),
        p_adj   = formatC(p_adj,   format = "g", digits = 3)
      ) %>%
      select(Metric, ANOVA_p, ANOVA_sig, Contrast, Estimate, SE, df, t.ratio,
             p.value, p_adj, signif)
  
    write.table(posthoc_df,
                file = file.path(out_dir, "lmm_posthoc.tsv"),
                sep = "\t", row.names = FALSE, quote = FALSE)
  }
  
  # Plotting
  p_metrics <- make_metrics_plot_grouped(df_metrics, method_colors, sig_levels)
  
  final_plot <- (p_pair_dup_Closed | p_pair_dup_Moderate) /
                (p_pair_dup_Open   | p_metrics)

  final_plot <- final_plot +
              plot_layout(guides = "collect") +
              plot_annotation(tag_levels = "A")
  
  ggsave(file.path(out_dir, paste0("dup_", sim_type, "_pair.png")), final_plot,
         width = 10, height = 6, dpi = 600, bg = "white")
         
  # Save as PDF
  ggsave(file.path(out_dir, paste0("dup_", sim_type, "_pair.pdf")), final_plot,
         width = 10, height = 6, bg = "white")
         
  # Metrics over D-value plots
  pair_metrics_df <- prepare_pairs_df(all_runs)
  metrics <- compute_bin_metrics_cat(pair_metrics_df, n_bins = 20)
  
  # Build plots per metric
  p_precision <- plot_bin_metrics_cat(metrics, method_colors, "precision")
  p_recall    <- plot_bin_metrics_cat(metrics, method_colors, "recall")
  p_f1        <- plot_bin_metrics_cat(metrics, method_colors, "f1")

  p_metrics_cat <- p_precision / p_recall / p_f1 +
    plot_annotation(tag_levels = "A")
    plot_layout(guides = "collect") &
    theme(legend.position = "bottom")
  
  ggsave(file.path(out_dir, "metrics_over_bins.png"),
         p_metrics_cat, width = 9, height = 10, dpi = 600, bg = "white")
  
  # Save as PDF
  ggsave(file.path(out_dir, "metrics_over_bins.pdf"),
         p_metrics_cat, width = 9, height = 10, bg = "white")
         
  # Save just recall
  ggsave(file.path(out_dir, "recall_over_bins.pdf"),
         p_recall, width = 9, height = 4, bg = "white")
}