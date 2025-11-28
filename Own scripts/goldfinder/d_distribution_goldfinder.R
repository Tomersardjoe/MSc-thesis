#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(readr)
  library(patchwork)
  library(stringdist)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Usage: Rscript goldfinder_d_distribution.R /path/to/coinfinder/nodes_all.tsv /path/to/coinfinder/d_cutoff.txt /path/to/simultaneous_association_significant_pairs.csv")
}

nodes_path     <- args[1]
cutoff_value   <- as.numeric(readLines(args[2], warn = FALSE))
pairs_path     <- args[3]
unique_id      <- basename(dirname(nodes_path))

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

# If no pairs, exit with status 1
if (nrow(pairs) == 0) {
  message("Exiting run ", unique_id, " - no valid pairs in", pairs_path)
  quit(save = "no", status = 1)
}

# Merge D-values for Gene_1 and Gene_2
pairs <- pairs %>%
  left_join(nodes %>% select(ID, Result) %>% rename(Gene_1 = ID, d_source = Result), by = "Gene_1") %>%
  left_join(nodes %>% select(ID, Result) %>% rename(Gene_2 = ID, d_target = Result), by = "Gene_2")

# -------------------------
# Gene-level D-values for all genes in significant pairs
# -------------------------
gene_ids <- unique(c(pairs$Gene_1, pairs$Gene_2))
nodes <- nodes %>% filter(ID %in% gene_ids)

# calculate percentage of genes above cutoff
total_genes   <- sum(!is.na(nodes$Result))
retained_genes <- sum(nodes$Result >= cutoff_value, na.rm = TRUE)
pct_retained   <- 100 * retained_genes / total_genes

# -------------------------
# Histogram of gene-level D-values
# -------------------------

# Calculate outliers above and below D-value of 5 and -5
n_low  <- sum(nodes$Result < -5, na.rm = TRUE)
n_high <- sum(nodes$Result >  5, na.rm = TRUE)

p_hist_genes <- ggplot(nodes, aes(x = Result)) +
  geom_histogram(
    binwidth = 0.25, boundary = cutoff_value,
    fill = "#6baed6", color = "#08519c", alpha = 0.6
  ) +
  geom_vline(
    xintercept = cutoff_value,
    linetype = "dashed", color = "#e31a1c", linewidth = 1
  ) +
    annotate("rect",
           xmin = cutoff_value, xmax = Inf, ymin = 0, ymax = Inf,
           alpha = 0.1, fill = "#e31a1c") +
  scale_x_continuous(
    breaks = seq(-5, 5, by = 1.0)
  ) +
  labs(
    title    = paste("Distribution of gene D-values -", unique_id),
    subtitle = paste0("D-value cutoff = ", round(cutoff_value, 3),
      " | Genes retained = ", retained_genes, "/", total_genes,
      " (", round(pct_retained, 1), "%)"
    ),
    x = "D-value", y = "Gene count"
  ) +
  coord_cartesian(xlim = c(-5, 5)) +
  theme_minimal()

# If outliers, annotate in plot
if (n_low > 0) {
  p_hist_genes <- p_hist_genes +
    annotate(
      "text",
      x = -Inf, y = Inf,
      label = paste0(n_low, " genes < D: -5"),
      hjust = -0.1, vjust = 2,
      color = "blue"
    )
}

if (n_high > 0) {
  p_hist_genes <- p_hist_genes +
    annotate(
      "text",
      x = -Inf, y = Inf,
      label = paste0(n_high, " genes > D: 5"),
      hjust = -0.1, vjust = 3.5,
      color = "red"
    )
}

# -------------------------
# Histogram of pair-level D-values
# -------------------------

retained_pairs <- sum(pairs$d_source >= cutoff_value & pairs$d_target >= cutoff_value, na.rm = TRUE)
total_pairs    <- nrow(pairs)
pct_retained   <- 100 * retained_pairs / total_pairs

# Pair-level D-value = minimum of the two genes' D-values
pairs$d_pair <- pmin(pairs$d_source, pairs$d_target, na.rm = TRUE)

# Calculate outliers above and below D-value of 5 and -5 at the pair level
n_low_pairs  <- sum(pairs$d_pair < -5, na.rm = TRUE)
n_high_pairs <- sum(pairs$d_pair >  5, na.rm = TRUE)

p_hist_pairs <- ggplot(pairs, aes(x = d_pair)) +
  geom_histogram(
    binwidth = 0.25, boundary = cutoff_value,
    fill = "#fdae6b", color = "#e6550d", alpha = 0.6
  ) +
  geom_vline(
    xintercept = cutoff_value,
    linetype = "dashed", color = "#e31a1c", linewidth = 1
  ) +
    annotate("rect",
           xmin = cutoff_value, xmax = Inf, ymin = 0, ymax = Inf,
           alpha = 0.1, fill = "#e31a1c") +
  scale_x_continuous(
    breaks = seq(-5, 5, by = 1.0)
  ) +
  labs(
    title    = paste("Distribution of gene pairs D-values -", unique_id),
    subtitle = paste0("D-value cutoff = ", round(cutoff_value, 3),
      " | Gene pairs retained = ", retained_pairs, "/", total_pairs,
      " (", round(pct_retained, 1), "%)"
    ),
    x = "D-value (min of source & target)",
    y = "Pair count"
  ) +
  coord_cartesian(xlim = c(-5, 5)) +
  theme_minimal()
  
# If outliers, annotate in plot
if (n_low_pairs > 0) {
  p_hist_pairs <- p_hist_pairs +
    annotate(
      "text",
      x = -Inf, y = Inf,
      label = paste0(n_low_pairs, " pairs < D: -5"),
      hjust = -0.1, vjust = 2,
      color = "blue"
    )
}

if (n_high_pairs > 0) {
  p_hist_pairs <- p_hist_pairs +
    annotate(
      "text",
      x = -Inf, y = Inf,
      label = paste0(n_high_pairs, " pairs > D: 5"),
      hjust = -0.1, vjust = 3.5,
      color = "red"
    )
}

# -------------------------
# Combine and save plots
# -------------------------
final_plot <- p_hist_genes + p_hist_pairs
output_path <- file.path(out_dir, paste0("d_distribution_", unique_id, ".png"))
ggsave(output_path, plot = final_plot, width = 12, height = 6, dpi = 300, bg = "white")

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

mean_partners <- mean(partners_per_gene$n_partners, na.rm = TRUE)

# Plot histogram of partner counts
p_partners <- ggplot(partners_per_gene, aes(x = n_partners)) +
  geom_histogram(binwidth = 1, boundary = 0, closed = "right",
                 fill = "#6baed6", color = "#08519c", alpha = 0.6) +
  labs(
    title = paste("Distribution of gene partners -", unique_id),
    subtitle = paste0(
      "Mean gene partners = ", round(mean_partners, 2)
    ),
    x = "Number of partners",
    y = "Gene count"
  ) +
  theme_minimal()

# Save
ggsave(file.path(out_dir, paste0("partners_hist_", unique_id, ".png")),
       plot = p_partners, width = 8, height = 6, dpi = 300, bg = "white")
       
# -------------------------
# Histogram: number of partners per gene with cutoff
# -------------------------

# Keep only pairs where both genes meet the cutoff (i.e., d_pair >= cutoff_value)
pairs_kept <- pairs %>%
  filter(!is.na(d_pair), d_pair >= cutoff_value)

# Build undirected partner list (count unique partners per gene)
partners_long <- bind_rows(
  pairs_kept %>% transmute(gene = Gene_1, partner = Gene_2),
  pairs_kept %>% transmute(gene = Gene_2, partner = Gene_1)
) %>%
  filter(!is.na(gene), !is.na(partner), gene != partner)

partners_per_gene <- partners_long %>%
  distinct(gene, partner) %>%               # unique partners only
  count(gene, name = "n_partners")          # degree per gene above cutoff

mean_partners <- mean(partners_per_gene$n_partners, na.rm = TRUE)

# Plot histogram of partner counts with cutoff
p_partners <- ggplot(partners_per_gene, aes(x = n_partners)) +
  geom_histogram(binwidth = 1, boundary = 0, closed = "right",
                 fill = "#6baed6", color = "#08519c", alpha = 0.6) +
  labs(
    title = paste("Distribution of gene partners >= D-value cutoff -", unique_id),
    subtitle = paste0(
      "D-value cutoff = ", round(cutoff_value, 3),
      " | Mean gene partners = ", round(mean_partners, 2)
    ),
    x = "Number of partners >= cutoff",
    y = "Gene count"
  ) +
  theme_minimal()

# Save
ggsave(file.path(out_dir, paste0("partners_hist_cutoff_", unique_id, ".png")),
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
