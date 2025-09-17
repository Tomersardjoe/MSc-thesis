#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(readr)
  library(tidyr)
  library(patchwork)
})

# Function to classify overlap (or no overlap) between the methods
classify_overlap <- function(methods) {
  if (setequal(methods, c("Coinfinder"))) return("Coinfinder only")
  if (setequal(methods, c("Goldfinder"))) return("Goldfinder only")
  if (setequal(methods, c("PanForest"))) return("PanForest only")
  if (setequal(methods, c("Coinfinder", "Goldfinder"))) return("Coinfinder + Goldfinder")
  if (setequal(methods, c("Coinfinder", "PanForest"))) return("Coinfinder + PanForest")
  if (setequal(methods, c("Goldfinder", "PanForest"))) return("Goldfinder + PanForest")
  if (setequal(methods, c("Coinfinder", "Goldfinder", "PanForest"))) return("All three methods")
  stop("Unexpected MethodsPresent combination: ", paste(methods, collapse = ", "))
}

# -------------------------
# Parse arguments
# -------------------------
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Usage: Rscript compare_dvalues.R coinfinder_dvalues.csv goldfinder_dvalues.csv panforest_dvalues.csv")
}

coin_path <- args[1]
gold_path <- args[2]
pan_path  <- args[3]

# HARD CODED FOR TESTING
# coin_path <- "1307/coinfinder_dvalues_1307.csv"
# gold_path <- "1307/goldfinder_dvalues_1307.csv"
# pan_path  <- "1307/panforest_dvalues_1307.csv"

# -------------------------
# Extract run_id from coin_path
# -------------------------
run_id <- sub("^.*_(\\d+)\\.csv$", "\\1", basename(coin_path))

# -------------------------
# Read data
# -------------------------
coin_df <- read_csv(coin_path, show_col_types = FALSE) %>% mutate(Method = "Coinfinder")
gold_df <- read_csv(gold_path, show_col_types = FALSE) %>% mutate(Method = "Goldfinder")
pan_df  <- read_csv(pan_path,  show_col_types = FALSE) %>% mutate(Method = "PanForest")

# -------------------------
# Combine into one data frame
# -------------------------
all_df <- bind_rows(coin_df, gold_df, pan_df)

#--------------------------
# Binning
#--------------------------

binwidth <- 0.25

# Flatten to gene level
singles <- all_df %>%
  filter(!is.na(Gene)) %>%
  transmute(Level, Method, GeneID = toupper(trimws(Gene)), D_value)

pairs_long <- all_df %>%
  filter(!is.na(Gene_1) & !is.na(Gene_2)) %>%
  pivot_longer(cols = c(Gene_1, Gene_2), names_to = "which", values_to = "GeneID") %>%
  mutate(GeneID = toupper(trimws(GeneID))) %>%
  select(Level, Method, GeneID, D_value)

gene_occ <- bind_rows(singles, pairs_long) %>%
  mutate(Bin = floor(D_value / binwidth) * binwidth) %>%
  distinct(Level, Bin, Method, GeneID)

# Flatten to pair level
pairs_df <- all_df %>%
  filter(!is.na(Gene_1) & !is.na(Gene_2)) %>%
  mutate(
    Gene_1 = toupper(trimws(Gene_1)),
    Gene_2 = toupper(trimws(Gene_2)),
    PairID = paste(pmin(Gene_1, Gene_2), pmax(Gene_1, Gene_2), sep = "-"),
    Bin    = floor(D_value / binwidth) * binwidth
  ) %>%
  distinct(Level, Bin, Method, PairID)

# Per-ID overlap tables
gene_method_overlap <- gene_occ %>%
  distinct(GeneID, Method, Bin) %>% # remove duplicates globally
  arrange(GeneID, Method, Bin) %>%  # sort once
  group_by(GeneID) %>%
  summarise(
    MethodBinMap   = list(tibble(Method = Method, Bin = Bin)),
    MethodsPresent = list(unique(Method)),
    BinsPresent    = list(unique(Bin)),
    .groups        = "drop"
  ) %>%
  mutate(
    OverlapCategory = vapply(MethodsPresent, classify_overlap, character(1)),
    SameBin         = lengths(BinsPresent) == 1
  )

pair_method_overlap <- pairs_df %>%
  distinct(PairID, Method, Bin) %>%
  arrange(PairID, Method, Bin) %>%
  group_by(PairID) %>%
  summarise(
    MethodBinMap   = list(tibble(Method = Method, Bin = Bin)),
    MethodsPresent = list(unique(Method)),
    BinsPresent    = list(unique(Bin)),
    .groups        = "drop"
  ) %>%
  mutate(
    OverlapCategory = vapply(MethodsPresent, classify_overlap, character(1)),
    SameBin         = lengths(BinsPresent) == 1
  )

# Per-bin counts for histograms
bin_counts <- gene_occ %>%
  group_by(Level, Bin, GeneID) %>%
  summarise(MethodsPresent = list(sort(unique(Method))), .groups = "drop") %>%
  mutate(OverlapCategory = vapply(MethodsPresent, classify_overlap, character(1))) %>%
  group_by(Level, Bin, OverlapCategory) %>%
  summarise(UniqueGenes = n(), .groups = "drop")

pair_bin_counts <- pairs_df %>%
  group_by(Level, Bin, PairID) %>%
  summarise(MethodsPresent = list(sort(unique(Method))), .groups = "drop") %>%
  mutate(OverlapCategory = vapply(MethodsPresent, classify_overlap, character(1))) %>%
  group_by(Level, Bin, OverlapCategory) %>%
  summarise(UniquePairs = n(), .groups = "drop")

# -------------------------
# Plotting
# -------------------------

# Fixed category order for consistent legends/colours
overlap_levels <- c(
  "Coinfinder only",
  "Goldfinder only",
  "PanForest only",
  "Coinfinder + Goldfinder",
  "Coinfinder + PanForest",
  "Goldfinder + PanForest",
  "All three methods"
)

# Fixed colour palette for those categories
overlap_palette <- c(
  "Coinfinder only"         = "#F28E2B",  # orange
  "Goldfinder only"         = "#59A14F",  # green
  "PanForest only"          = "#4E79A7",  # blue
  "Coinfinder + Goldfinder" = "#B6992D",  # olive
  "Coinfinder + PanForest"  = "#A76BAA",  # purple
  "Goldfinder + PanForest"  = "#2C9FBF",  # teal
  "All three methods"       = "#8C564B"   # brown
)

method_palette <- c(
  Coinfinder = unname(overlap_palette["Coinfinder only"]),
  Goldfinder = unname(overlap_palette["Goldfinder only"]),
  PanForest  = unname(overlap_palette["PanForest only"])
)

legend_df <- data.frame(
  OverlapCategory = factor(overlap_levels, levels = overlap_levels)
)

# Gene histogram
p_hist_genes <- ggplot(
  bin_counts %>% filter(Level == "Gene") %>%
    mutate(OverlapCategory = factor(OverlapCategory, levels = overlap_levels)),
  aes(x = Bin, y = UniqueGenes, fill = OverlapCategory)
) +
  geom_col(position = "stack", color = "black", alpha = 0.85) +
  facet_wrap(~Level, scales = "free_y") +
  scale_fill_manual(values = overlap_palette, drop = FALSE) +
  labs(
    title = paste("Unique genes per D-value bin (width =", binwidth, ") - Run", run_id),
    x = "D-value",
    y = "Number of unique genes",
    fill = "Overlap category"
  ) +
  theme_minimal() +
  guides(fill = "none")

# Pair histogram
p_hist_pairs <- ggplot(
  pair_bin_counts %>% filter(Level == "Pair") %>%
    mutate(OverlapCategory = factor(OverlapCategory, levels = overlap_levels)),
  aes(x = Bin, y = UniquePairs, fill = OverlapCategory)
) +
  geom_col(position = "stack", color = "black", alpha = 0.85) +
  
  # Force legend keys to be full-size rectangles
  guides(fill = guide_legend(
    override.aes = list(
      color = "black",
      size = 0.8,
      alpha = 1
    ),
    keywidth = 1.2,
    keyheight = 0.8
  )) +
  
  facet_wrap(~Level, scales = "free_y") +
  scale_fill_manual(values = overlap_palette, drop = FALSE) +
  labs(
    title = paste("Unique gene pairs per D-value bin (width =", binwidth, ") - Run", run_id),
    x = "D-value",
    y = "Number of unique gene pairs",
    fill = "Overlap category"
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")

# Gene boxplot
p_box_genes <- ggplot(all_df %>% filter(!is.na(Gene)), aes(x = Method, y = D_value, fill = Method)) +
  geom_boxplot(outlier.alpha = 0.3) +
  facet_wrap(~Level, scales = "free_y") +
  scale_fill_manual(values = method_palette) +
  labs(
    title = "D-value distributions (genes)",
    x = "Method",
    y = "D-value"
  ) +
  theme_minimal() +
  theme(legend.position = "none")

# Pair boxplot
p_box_pairs <- ggplot(all_df %>% filter(!is.na(Gene_1) & !is.na(Gene_2)), aes(x = Method, y = D_value, fill = Method)) +
  geom_boxplot(outlier.alpha = 0.3) +
  scale_fill_manual(values = method_palette) +
  facet_wrap(~Level, scales = "free_y") +
  labs(
    title = "D-value distributions (pairs)",
    x = "Method",
    y = "D-value"
  ) +
  theme_minimal() +
  theme(legend.position = "none")

# Combine plots

# Top row
top_row <- (p_hist_genes | p_hist_pairs) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

# Bottom row
bottom_row <- (p_box_genes | p_box_pairs) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

# Stack the two rows
combined_plot <- top_row / bottom_row +
  plot_layout(heights = c(2,1))

# -------------------------
# Summary statistics
# -------------------------

# Summary stats per category for genes
genes_id <- gene_method_overlap %>%
  group_by(OverlapCategory) %>%
  summarise(
    n         = n(),  # number of unique genes in that category
    mean_bin  = mean(unlist(BinsPresent), na.rm = TRUE),
    median_bin= median(unlist(BinsPresent), na.rm = TRUE),
    min_bin   = min(unlist(BinsPresent), na.rm = TRUE),
    max_bin   = max(unlist(BinsPresent), na.rm = TRUE),
    .groups   = "drop"
  )

genes_bin <- gene_occ %>%
  group_by(Level, Bin, GeneID) %>%
  summarise(MethodsPresent = list(sort(unique(Method))), .groups = "drop") %>%
  mutate(OverlapCategory = vapply(MethodsPresent, classify_overlap, character(1))) %>%
  group_by(OverlapCategory, Level) %>%
  summarise(
    n         = n(),  # number of unique genes in that category
    mean_d    = mean(Bin, na.rm = TRUE),
    median_d  = median(Bin, na.rm = TRUE),
    min_d     = min(Bin, na.rm = TRUE),
    max_d     = max(Bin, na.rm = TRUE),
    .groups   = "drop"
  )

# Summary stats per category for pairs
pairs_id <- pair_method_overlap %>%
  group_by(OverlapCategory) %>%
  summarise(
    n         = n(),  # number of unique pairs in that category
    mean_bin  = mean(unlist(BinsPresent), na.rm = TRUE),
    median_bin= median(unlist(BinsPresent), na.rm = TRUE),
    min_bin   = min(unlist(BinsPresent), na.rm = TRUE),
    max_bin   = max(unlist(BinsPresent), na.rm = TRUE),
    .groups   = "drop"
  )

pairs_bin <- pairs_df %>%
  group_by(Level, Bin, PairID) %>%
  summarise(MethodsPresent = list(sort(unique(Method))), .groups = "drop") %>%
  mutate(OverlapCategory = vapply(MethodsPresent, classify_overlap, character(1))) %>%
  group_by(OverlapCategory, Level) %>%
  summarise(
    n         = n(),  # number of unique pairs in that category
    mean_d    = mean(Bin, na.rm = TRUE),
    median_d  = median(Bin, na.rm = TRUE),
    min_d     = min(Bin, na.rm = TRUE),
    max_d     = max(Bin, na.rm = TRUE),
    .groups   = "drop"
  )

# Summary stats by method
summary_by_method <- all_df %>%
  group_by(Method, Level) %>%
  summarise(
    n = n(),
    mean_d = mean(D_value, na.rm = TRUE),
    median_d = median(D_value, na.rm = TRUE),
    min_d = min(D_value, na.rm = TRUE),
    max_d = max(D_value, na.rm = TRUE),
    .groups = "drop"
  )

# -------------------------
# Save outputs
# -------------------------
out_dir <- file.path("d_values_comparison", run_id)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

write_csv(genes_id, file.path(out_dir, paste0("dvalues_genes_id_", run_id, ".csv")))
write_csv(pairs_id, file.path(out_dir, paste0("dvalues_pairs_id_", run_id, ".csv")))
write_csv(genes_bin, file.path(out_dir, paste0("dvalues_genes_bin_", run_id, ".csv")))
write_csv(pairs_bin, file.path(out_dir, paste0("dvalues_pairs_bin_", run_id, ".csv")))
write_csv(summary_by_method, file.path(out_dir, paste0("dvalues_stats_methods_", run_id, ".csv")))

ggsave(file.path(out_dir, paste0("dvalues_comparison_plots_", run_id, ".png")),
       plot = combined_plot, width = 16, height = 10, dpi = 300, bg = "white")

message("Saved: summary stats, gene overlap details and plots to: ", out_dir)
