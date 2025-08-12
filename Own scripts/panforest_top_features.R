library(data.table)
library(igraph)
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

# Import importance matrix
imp_mat <- as.matrix(read.csv("imp.csv", row.names=1, check.names=FALSE))

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
write_graph(g, "filtered_network.graphml", format = "graphml")

