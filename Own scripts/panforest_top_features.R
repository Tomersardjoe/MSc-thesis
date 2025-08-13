#!/usr/bin/env Rscript

library(data.table)
suppressPackageStartupMessages(library(igraph))
library(Matrix)

# Normalize percentile function
normalize_percentile <- function(mat) {
  mat <- (mat + t(mat)) / 2  # enforce symmetry
  diag(mat) <- 0             # remove self-loop influence
  w <- as.numeric(mat[upper.tri(mat)])
  r <- rank(w, ties.method = "average") / length(w)  # percentile (0..1)
  norm_mat <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  norm_mat[upper.tri(norm_mat)] <- r
  norm_mat <- norm_mat + t(norm_mat)
  diag(norm_mat) <- 0
  norm_mat
}

# Threshold top 1% edges using percentile rank
threshold_by_percentile <- function(mat, q = 0.99) {
  thr <- quantile(mat[upper.tri(mat)], probs = q, na.rm = TRUE)
  keep <- mat >= thr
  keep[lower.tri(keep)] <- t(keep)[lower.tri(keep)]
  diag(keep) <- FALSE
  mat * keep
}

# Get file paths from arguments: presence file, importance matrix, graphml output
args <- commandArgs(trailingOnly = TRUE)
verbose <- FALSE

# Look for --verbose flag anywhere in the args
if ("--verbose" %in% args) {
  verbose <- TRUE
  args <- setdiff(args, "--verbose")
}

if (length(args) < 3) {
  cat("\nUsage: Rscript panforest_top_features.R <presence_csv> <input_csv> <output_graphml> [--verbose]\n",
      "       <presence_csv>    Path to gene_presence_absence.csv\n",
      "       <input_csv>       Path to input importance matrix CSV file\n",
      "       <output_graphml>  Path to output GraphML file\n",
      "       [--verbose]       Optional flag to print extra info to stderr\n\n")
  stop("Error: Missing required arguments.\n", call. = FALSE)
}

presence_file <- args[1]
input_file <- args[2]
output_file <- args[3]

# Ensure .graphml extension
if (!grepl("\\.graphml$", output_file, ignore.case = TRUE)) {
  output_file <- paste0(output_file, ".graphml")
}

if (verbose) {
cat("Presence file:", presence_file, "\n", file = stderr())
cat("Input file   :", input_file, "\n", file = stderr())
cat("Output file  :", output_file, "\n", file = stderr())
}

# Derive dataset name from directory containing the presence file
dataset_name <- basename(normalizePath(dirname(presence_file)))

# Total pangenome genes (rows minus header)
if (!file.exists(presence_file)) {
  stop("Presence file not found: ", presence_file)
}
total_pangenome_genes <- max(0, nrow(fread(presence_file))) # fread skips header

# Import importance matrix
imp_mat <- as.matrix(read.csv(input_file, row.names=1, check.names=FALSE))

# Normalize importance matrix to percentile ranks
imp_mat <- normalize_percentile(imp_mat)

# Apply 0.99 percentile edge thresholding
imp_mat_thr <- threshold_by_percentile(imp_mat, q = 0.99)

# Sort genes by their aggregate pairwise importance
imp_scores <- rowSums(imp_mat_thr)
imp_scores <- scale(imp_scores) # z-score normalization
genes <- sort(imp_scores[,1], decreasing=TRUE)

# Convert importance matrix to data.table
dt <- as.data.table(as.table(imp_mat_thr))
dt <- dt[V1 != V2 & N > 0] # remove self-interactions and zero-weight edges
setnames(dt, c("V1", "V2", "N"), c("Gene_1", "Gene_2", "Importance"))

# Sort gene pairs by their importance
pairs <- dt[order(-Importance)]
colnames(pairs) <- c("Gene_1", "Gene_2", "Importance")

# Remove parallel edges from gene pairs
dt[, gene_pair := paste(pmin(Gene_1, Gene_2), pmax(Gene_1, Gene_2), sep = "_")]
dt <- dt[!duplicated(gene_pair)]
dt[, gene_pair := NULL] # optional cleanup

# Create igraph object from filtered edge list
g <- graph_from_data_frame(pairs, directed = FALSE)
E(g)$weight <- pairs$Importance

# Export graph in GraphML format
V(g)$name <- V(g)$name
V(g)$zscore <- genes[V(g)$name]
write_graph(g, output_file, format = "graphml")

# Network statistics from ML network
network_nodes <- gorder(g)
possible_pairs <- if (network_nodes > 1) network_nodes * (network_nodes - 1) / 2 else 0
significant_associations <- gsize(g)
association_rate <- if (possible_pairs > 0) {
  round((significant_associations / possible_pairs) * 100, 2)
} else 0.00

# Cluster statistics (Louvain)
if (network_nodes > 0 && significant_associations > 0) {
  comm <- cluster_louvain(g)
  module_count <- length(unique(membership(comm)))
  avg_genes_per_module <- mean(sizes(comm))
  modularity_val <- modularity(comm)
} else {
  module_count <- 0
  avg_genes_per_module <- 0
  modularity_val <- NA_real_
}

# Average degree
avg_degree <- if (network_nodes > 0) mean(degree(g)) else 0

# Emit metrics line (no header)
metrics_line <- data.table(
  Dataset = dataset_name,
  Total_Pangenome_Genes = total_pangenome_genes,
  Total_Networked_Genes = network_nodes,
  Possible_Gene_Pairs = possible_pairs,
  Significant_Associations = significant_associations,
  Association_Rate = round(association_rate, 2),
  Module_Count = module_count,
  Avg_Genes_per_Module = round(avg_genes_per_module, 2),
  Avg_Degree = round(avg_degree, 2),
  Modularity = round(modularity_val, 2)
)

cat(paste(metrics_line, collapse = ","), "\n")
