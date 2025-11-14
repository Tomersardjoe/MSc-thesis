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
    pattern      = factor(pattern, levels = c("none","stripe","crosshatch")),
    sets         = lapply(sets, function(x) factor(x,
                          levels = c("Coinfinder","Goldfinder","PanForest")))
  )

  p <- ggplot(df2, aes(x = sets, y = value)) +
    geom_col(aes(fill = base_tool), colour = "black", linewidth = 0.2, width = 0.7, show.legend = c(fill = TRUE)) +
    geom_col(data = dplyr::filter(df2, n_sets == 3), aes(y = value),
             fill = "#999999", colour = NA, width = 0.7, show.legend = FALSE) +
    geom_col_pattern(
      aes(pattern = pattern, pattern_fill = overlay_tool),
      fill = NA, colour = NA, pattern_colour = "black",
      linewidth = 0.2, width = 0.7, pattern_size = 0.2,
      show.legend = c(pattern = show_pattern_legend, pattern_fill = FALSE, fill = FALSE),
      key_glyph = ggpattern::draw_key_polygon_pattern,
      pattern_angle = 45, pattern_density = 0.5, pattern_spacing = 0.05
    ) +
    geom_label(aes(y = value, label = value),
               vjust = -0.3, size = 3, label.size = 0, fill = "white") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.2)), labels = scales::label_number()) +
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

# -------------------------
# Style constants
# -------------------------
method_colors <- c(
  Coinfinder = "#D95F5F",
  Goldfinder = "#FFD966",
  PanForest  = "#5F9ED1"
)

prop_colors <- c(
  closed = "#D44949",
  moderate = "#EC9837",
  open = "#58D168",
  wildcard = "#824992"
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
category_order <- c("closed", "moderate", "open", "wildcard")

# -------------------------
# Parse arguments
# -------------------------
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 5) {
  stop("Usage: Rscript aggregated_upset.R
        /path/to/{dataset}/coinfinder_runs(_all)
        /path/to/{dataset}/goldfinder_runs(_all)
        /path/to/{dataset}/panforest_runs(_all)
        /path/to/real_pangenomes/species_categories.csv
        /path/to/out_dir")
}

coin_dir        <- args[1]
gold_dir        <- args[2]
pan_dir         <- args[3]
categories_path <- args[4]
out_dir         <- args[5]

# -------------------------
# Data loading
# -------------------------

species_taxids <- list.dirs(coin_dir, recursive = FALSE, full.names = FALSE)

all_runs <- purrr::map_dfr(species_taxids, function(species_taxid) {
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
})

categories <- readr::read_csv(
  categories_path,
  col_names = c("category", "species_taxid"),
  col_select = 1:2,
  show_col_types = FALSE
)


all_runs <- all_runs %>%
  dplyr::left_join(categories, by = "species_taxid")

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
chisq_gene <- chisq.test(gene_cont_table)

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
    p_value = prop.test(
      x = n_all_three_genes,
      n = n_total_genes,
      p = gene_overall_mean
    )$p.value
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
combined <- (p_gene_closed | p_gene_moderate) /
            (p_gene_open   | p_gene_prop)

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
chisq_pair <- chisq.test(pair_cont_table)

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
    p_value = prop.test(
      x = n_all_three_pairs,
      n = n_total_pairs,
      p = overall_mean
    )$p.value
  ) %>%
  ungroup() %>%
  mutate(
    p_adj = p.adjust(p_value, method = "BH")
  )

category_tests <- category_tests %>%
  mutate(
    signif = case_when(
      p_adj < 0.001 ~ "***",
      p_adj < 0.01  ~ "**",
      p_adj < 0.05  ~ "*",
      TRUE            ~ "ns"
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

combined <- (p_pair_closed | p_pair_moderate) /
            (p_pair_open   | p_pair_prop)

final_plot <- combined + plot_layout(guides = "collect", heights = c(1, 1)) +
              plot_annotation(tag_levels = "A")

# Save as PNG
ggsave(file.path(out_dir, "pair_category.png"), final_plot,
       width = 8, height = 5, dpi = 600, bg = "white")

# Save as PDF
ggsave(file.path(out_dir, "pair_category.pdf"), final_plot,
       width = 8, height = 5, bg = "white")
