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

bin_w <- 0.25

# Define common breaks for both plots
dval_breaks <- seq(
  floor(min(c(nodes_in_pairs$Result, pairs$d_gene1, pairs$d_gene2), na.rm = TRUE)),
  ceiling(max(c(nodes_in_pairs$Result, pairs$d_gene1, pairs$d_gene2), na.rm = TRUE)),
  by = bin_w
)

p_hist_genes <- ggplot(nodes_in_pairs, aes(x = Result)) +
  geom_histogram(
    binwidth = bin_w,
    boundary = 0,  # align bins to 0 so 0.25 intervals are clean
    fill = "#6baed6", color = "#08519c", alpha = 0.6
  ) +
  scale_x_continuous(
    breaks = dval_breaks,
    labels = ifelse(dval_breaks %% 0.5 == 0, dval_breaks, "")
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
    breaks = dval_breaks,
    labels = ifelse(dval_breaks %% 0.5 == 0, dval_breaks, "")
  ) +
  labs(
    title = "Distribution of D-values - Significant gene pairs",
    x = "D-value (min of source & target)",
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
# Histogram: number of partners per gene
# -------------------------

# Build undirected partner list (count unique partners per gene)
partners_long <- bind_rows(
  pairs %>% transmute(gene = Gene_1, partner = Gene_2),
  pairs %>% transmute(gene = Gene_2, partner = Gene_1)
) %>%
  filter(!is.na(gene), !is.na(partner), gene != partner)

partners_per_gene <- partners_long %>%
  distinct(gene, partner) %>%   # unique partners only
  count(gene, name = "n_partners")

# Plot histogram of partner counts
p_partners <- ggplot(partners_per_gene, aes(x = n_partners)) +
  geom_histogram(binwidth = 1, boundary = 0, closed = "right",
                 fill = "#6baed6", color = "#08519c", alpha = 0.6) +
  labs(
    title = paste("Distribution of gene partners -", unique_id),
    x = "Number of partners",
    y = "Gene count"
  ) +
  theme_minimal()

# Save
ggsave(file.path(out_dir, paste0("partners_hist_", unique_id, ".png")),
       plot = p_partners, width = 8, height = 6, dpi = 300, bg = "white")

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
