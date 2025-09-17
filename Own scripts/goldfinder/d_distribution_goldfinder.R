#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(readr)
  library(patchwork)
  library(stringdist)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript goldfinder_d_distribution.R /path/to/coinfinder/nodes_all.tsv /path/to/simultaneous_association_significant_pairs.csv")
}

nodes_path <- args[1]
pairs_path <- args[2]
unique_id  <- basename(dirname(nodes_path))

# Output directory
out_dir <- file.path(dirname(pairs_path), "d_distribution")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# -------------------------
# Read Coinfinder nodes_all.tsv (gene-level D-values)
# -------------------------
nodes <- read.table(nodes_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
nodes$Result <- as.numeric(nodes$Result)

if (all(is.na(nodes$Result))) {
  stop("No valid D-values in nodes file.")
}

# -------------------------
# Read Goldfinder significant pairs CSV
# -------------------------
pairs <- read_csv(pairs_path, show_col_types = FALSE) %>%
  mutate(
    Gene_1 = as.character(Gene_1),
    Gene_2 = as.character(Gene_2)
  )

# Merge D-values for Gene_1 and Gene_2
pairs <- pairs %>%
  left_join(nodes %>% select(ID, Result) %>% rename(Gene_1 = ID, d_gene1 = Result), by = "Gene_1") %>%
  left_join(nodes %>% select(ID, Result) %>% rename(Gene_2 = ID, d_gene2 = Result), by = "Gene_2")

# -------------------------
# Gene-level D-values for all genes in significant pairs
# -------------------------
genes_in_pairs <- unique(c(pairs$Gene_1, pairs$Gene_2))
nodes_in_pairs <- nodes %>% filter(ID %in% genes_in_pairs)

#missing_genes <- setdiff(genes_in_pairs, nodes$ID)

#closest_matches <- lapply(missing_genes, function(mg) {
#  dists <- stringdist(mg, nodes$ID, method = "osa")
#  best_idx <- which.min(dists)
#  data.frame(
#    missing_gene = mg,
#    closest_match = nodes$ID[best_idx],
#    distance = dists[best_idx]
#  )
#})

#closest_matches_df <- do.call(rbind, closest_matches)
#closest_matches_df

bin_w <- 0.25

p_hist_genes <- ggplot(nodes_in_pairs, aes(x = Result)) +
  geom_histogram(
    binwidth = bin_w,
    boundary = 0,  # align bins to 0 so 0.25 intervals are clean
    fill = "#6baed6", color = "#08519c", alpha = 0.6
  ) +
  scale_x_continuous(
    breaks = seq(
      floor(min(nodes_in_pairs$Result, na.rm = TRUE)),
      ceiling(max(nodes_in_pairs$Result, na.rm = TRUE)),
      by = bin_w
    )
  ) +
  labs(
    title = paste("Distribution of D-values - Genes in significant pairs (", unique_id, ")", sep = ""),
    x = "D-value", y = "Gene count"
  ) +
  theme_minimal()

# -------------------------
# Pair-level D-values (min of the two genes)
# -------------------------
pairs$d_pair <- pmin(pairs$d_gene1, pairs$d_gene2, na.rm = TRUE)

p_hist_pairs <- ggplot(pairs, aes(x = d_pair)) +
  geom_histogram(
    binwidth = bin_w,
    boundary = 0,
    fill = "#fdae6b", color = "#e6550d", alpha = 0.6
  ) +
  scale_x_continuous(
    breaks = seq(
      floor(min(pairs$d_pair, na.rm = TRUE)),
      ceiling(max(pairs$d_pair, na.rm = TRUE)),
      by = bin_w
    )
  ) +
  labs(
    title = "Distribution of D-values - Significant gene pairs",
    x = "D-value (pair-level, min of Gene_1 & Gene_2)",
    y = "Pair count"
  ) +
  theme_minimal()

# -------------------------
# Combine and save (matching Coinfinder style)
# -------------------------
final_plot <- p_hist_genes + p_hist_pairs

ggsave(file.path(out_dir, paste0("goldfinder_d_distributions_", unique_id, ".png")),
       plot = final_plot, width = 12, height = 6, dpi = 300, bg = "white")

message("Goldfinder D-value distribution plots saved to: ", out_dir)

# -------------------------
# Write output data for comparative plots
# -------------------------

# Gene-level data
df_genes <- nodes %>%
  transmute(
    Gene = ID,
    Gene_1 = NA_character_,
    Gene_2 = NA_character_,
    D_value = Result,
    Method = "Goldfinder",
    Level = "Gene"
  ) %>%
  filter(!is.na(D_value))

# Pair-level data
df_pairs <- pairs %>%
  transmute(
    Gene = NA_character_,
    Gene_1 = Gene_1,
    Gene_2 = Gene_2,
    D_value = d_pair,
    Method = "Goldfinder",
    Level = "Pair"
  ) %>%
  filter(!is.na(D_value))

# Combine
df_out <- bind_rows(df_genes, df_pairs)

# Write to CSV
write_csv(df_out, file.path(out_dir, paste0("goldfinder_dvalues_", unique_id, ".csv")))

message("D-value CSV saved to: ",
        file.path(out_dir, paste0("goldfinder_dvalues_", unique_id, ".csv")))
