# fluidity_calc.R

suppressPackageStartupMessages(library(micropan))

# Load arguments: first argument is the CSV file path
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  stop("Please provide the path to your presence–absence CSV file.\nUsage: Rscript fluidity_calc.R gene_presence_absence.csv")
}
pa_path <- args[1]

# Read matrix
pa <- read.csv(pa_path, header = TRUE, row.names = 1, check.names = FALSE)
pa <- as.matrix(pa)
stopifnot(all(pa %in% c(0,1)))
storage.mode(pa) <- "integer"

# Build genome sets
genome_sets <- lapply(seq_len(ncol(pa)), function(j) which(pa[,j] == 1))

# Pairwise fluidity
pair_fluidity <- function(a, b) {
  inter <- length(intersect(a,b))
  ua    <- length(a)
  ub    <- length(b)
  union <- ua + ub - inter
  unique_only <- (ua - inter) + (ub - inter)
  unique_only / union
}

# Calculate mean fluidity
phi_vals <- c()
for (i in 1:(length(genome_sets)-1)) {
  for (j in (i+1):length(genome_sets)) {
    phi_vals <- c(phi_vals, pair_fluidity(genome_sets[[i]], genome_sets[[j]]))
  }
}
genomic_fluidity <- mean(phi_vals)

# Openness calculation
pa_t <- t(pa)  # genomes in rows

# uncomment and adjust values for debugging
#slow_cols <- 1:320
#fast_cols <- 321:ncol(pa_t)

#slow <- pa_t[, slow_cols, drop = FALSE]
#fast <- pa_t[, fast_cols, drop = FALSE]

# Drop all-zero genes within each set
#slow <- slow[, colSums(slow) > 0, drop = FALSE]
#fast <- fast[, colSums(fast) > 0, drop = FALSE]
pa_t_nz <- pa_t[, colSums(pa_t) > 0, drop = FALSE]

gamma_all  <- heaps(pa_t_nz, n.perm=500)[["alpha"]]
#gamma_slow <- heaps(slow, n.perm=500)[["alpha"]]
#gamma_fast <- heaps(fast, n.perm=500)[["alpha"]]

cat("Genomic fluidity:", genomic_fluidity, "\n")
cat("Pangenome openness - all genes (alpha):", gamma_all, "\n")
#cat("Pangenome openness - slow partition (alpha):", gamma_slow, "\n")
#cat("Pangenome openness - fast partition (alpha):", gamma_fast, "\n")
