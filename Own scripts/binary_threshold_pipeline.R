#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(tidyverse)
})

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 4) {
  stop("Usage: Rscript binary_threshold_pipeline.R <coinfinder.tsv> <gold_assoc.csv> <gold_diss.csv> <panforest.csv>")
}

coin_file  <- args[1]
gold_assoc_file <- args[2]
gold_diss_file  <- args[3]
pan_file   <- args[4]

coin <- read_tsv(
  coin_file,
  col_types = cols(
    Source = col_character(),
    Target = col_character(),
    .default = col_double()
  )
) %>% select(Source, Target)

gold_assoc <- read_csv(
  gold_assoc_file,
  col_types = cols(
    Gene_1 = col_character(),
    Gene_2 = col_character(),
    `p-value unadj` = col_double(),
    `p-value adj` = col_double(),
    Cluster = col_character()
  )
) %>% select(Gene_1, Gene_2)

gold_diss <- read_csv(
  gold_diss_file,
  col_types = cols(
    Gene_1 = col_character(),
    Gene_2 = col_character(),
    .default = col_double()
  )
) %>% select(Gene_1, Gene_2)

gold <- bind_rows(gold_assoc, gold_diss)

pan <- suppressMessages(
  read_csv(
    pan_file,
    col_types = cols(
      `...1` = col_character(),
      .default = col_double()
    )
  )
) %>%
  rename(Genes = `...1`) %>%
  column_to_rownames("Genes") %>%
  as.matrix()

genes_coin <- unique(c(coin$Source, coin$Target))
genes_gold <- unique(c(gold$Gene_1, gold$Gene_2))
genes_pan  <- rownames(pan)
genes <- sort(unique(c(genes_coin, genes_gold, genes_pan)))

pairs_to_matrix <- function(df, g1, g2, genes) {
  mat <- matrix(0, nrow = length(genes), ncol = length(genes),
                dimnames = list(genes, genes))
  for (i in seq_len(nrow(df))) {
    a <- df[[g1]][i]
    b <- df[[g2]][i]
    if (a %in% genes && b %in% genes) {
      mat[a, b] <- 1
      mat[b, a] <- 1
    }
  }
  diag(mat) <- 0
  return(mat)
}

M_coin <- pairs_to_matrix(coin, "Source", "Target", genes)
M_gold <- pairs_to_matrix(gold, "Gene_1", "Gene_2", genes)

threshold <- quantile(pan[upper.tri(pan)], 0.95, na.rm = TRUE)
M_pan_raw <- ifelse(pan >= threshold, 1, 0)
diag(M_pan_raw) <- 0
M_pan <- matrix(0, nrow = length(genes), ncol = length(genes),
                dimnames = list(genes, genes))
common_genes <- intersect(genes, genes_pan)
M_pan[common_genes, common_genes] <- M_pan_raw[common_genes, common_genes]

jaccard_index <- function(mat1, mat2) {
  m1 <- mat1[upper.tri(mat1)]
  m2 <- mat2[upper.tri(mat2)]
  intersection <- sum(m1 == 1 & m2 == 1)
  union <- sum(m1 == 1 | m2 == 1)
  if (union == 0) return(NA)
  return(intersection / union)
}

# Calculate Jaccard values
j_coin_gold <- jaccard_index(M_coin, M_gold)
j_coin_pan  <- jaccard_index(M_coin, M_pan)
j_gold_pan  <- jaccard_index(M_gold, M_pan)

# Convert adjacency matrix to a vector of "geneA_geneB" edge strings
pairs_from_matrix <- function(M) {
  g <- rownames(M)
  ut <- which(upper.tri(M) & M == 1, arr.ind = TRUE)
  paste(g[ut[,1]], g[ut[,2]], sep = "_")
}

coin_pairs <- pairs_from_matrix(M_coin)
gold_pairs <- pairs_from_matrix(M_gold)
pan_pairs  <- pairs_from_matrix(M_pan)

# Count overlaps and unions from edge sets
overlap_coin_gold <- length(intersect(coin_pairs, gold_pairs))
union_coin_gold   <- length(union(coin_pairs, gold_pairs))

overlap_coin_pan  <- length(intersect(coin_pairs, pan_pairs))
union_coin_pan    <- length(union(coin_pairs, pan_pairs))

overlap_gold_pan  <- length(intersect(gold_pairs, pan_pairs))
union_gold_pan    <- length(union(gold_pairs, pan_pairs))

# Print summary
cat("\n=== Pairwise Comparison Summary ===\n")
cat(sprintf("Coinfinder vs Goldfinder : %d shared / %d total  (Jaccard = %.5f)\n",
            overlap_coin_gold, union_coin_gold, overlap_coin_gold / union_coin_gold))
cat(sprintf("Coinfinder vs PanForest  : %d shared / %d total  (Jaccard = %.5f)\n",
            overlap_coin_pan, union_coin_pan, overlap_coin_pan / union_coin_pan))
cat(sprintf("Goldfinder vs PanForest  : %d shared / %d total  (Jaccard = %.5f)\n",
            overlap_gold_pan, union_gold_pan, overlap_gold_pan / union_gold_pan))
cat("====================================\n\n")

# Uncomment to write binary matrix files
#write.csv(M_coin, "binary_coinfinder.csv", row.names = TRUE)
#write.csv(M_coin, "binary_goldfinder.csv", row.names = TRUE)
#write.csv(M_coin, "binary_panforest.csv", row.names = TRUE)

# Read ground truth
truth_df <- read_tsv("simulation/old/sim_large/simulation/pairs.txt", # Hardcoded for testing
                     col_names = FALSE,
                     col_types = cols(.default = col_character()))

# Assign column names
colnames(truth_df) <- c("Gene1", "Gene2")

# Canonicalise so order doesn't matter
canonical_pair <- function(a, b) {
  ifelse(a < b, paste(a, b, sep = "_"), paste(b, a, sep = "_"))
}
truth_pairs <- canonical_pair(truth_df$Gene1, truth_df$Gene2)

# Extract edge sets from matrices
pairs_from_matrix <- function(M) {
  g <- rownames(M)
  ut <- which(upper.tri(M) & M == 1, arr.ind = TRUE)
  paste(g[ut[,1]], g[ut[,2]], sep = "_")
}

coin_pairs <- pairs_from_matrix(M_coin)
gold_pairs <- pairs_from_matrix(M_gold)
pan_pairs  <- pairs_from_matrix(M_pan)

# Compare a tool's edges to truth
compare_to_truth <- function(tool_pairs, truth_pairs) {
  tp <- length(intersect(tool_pairs, truth_pairs))
  fp <- length(setdiff(tool_pairs, truth_pairs))
  fn <- length(setdiff(truth_pairs, tool_pairs))
  precision <- if ((tp + fp) == 0) NA else tp / (tp + fp)
  recall    <- if ((tp + fn) == 0) NA else tp / (tp + fn)
  f1        <- if (is.na(precision) || is.na(recall) || (precision + recall) == 0) NA else 2 * precision * recall / (precision + recall)
  list(tp = tp, fp = fp, fn = fn, precision = precision, recall = recall, f1 = f1)
}

# Run comparisons
coin_vs_truth <- compare_to_truth(coin_pairs, truth_pairs)
gold_vs_truth <- compare_to_truth(gold_pairs, truth_pairs)
pan_vs_truth  <- compare_to_truth(pan_pairs, truth_pairs)

# Print results
cat("\n=== Overlap with Ground Truth ===\n")
print_tool <- function(name, res) {
  cat(sprintf("%-12s: TP=%d, FP=%d, FN=%d, Precision=%.3f, Recall=%.3f, F1=%.3f\n",
              name, res$tp, res$fp, res$fn, res$precision, res$recall, res$f1))
}
print_tool("Coinfinder", coin_vs_truth)
print_tool("Goldfinder", gold_vs_truth)
print_tool("PanForest",  pan_vs_truth)
cat("=================================\n\n")

