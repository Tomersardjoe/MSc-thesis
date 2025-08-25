#!/usr/bin/env Rscript

# Explicit, reproducible p-adjustment for Coinfinder <dataset>_pairs.tsv
# Usage: Rscript p_adj_coinfinder.R path to <dataset> (e.g. Coinfinder/sim_small)

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 1) {
    stop("Usage: p_adj_coinfinder.R <dataset_directory> [significance_threshold, default=0.05]")
}

dataset_dir <- args[1]

if (!dir.exists(dataset_dir)) {
    stop(paste("Directory does not exist:", dataset_dir))
}

# Find *_pairs.tsv in that directory
pairs_files <- list.files(dataset_dir, pattern = "_pairs\\.tsv$", full.names = TRUE)

if (length(pairs_files) == 0) {
    stop("No file ending in _pairs.tsv found in ", dataset_dir)
} else if (length(pairs_files) > 1) {
    stop("Multiple _pairs.tsv files found in ", dataset_dir, 
         ". Please ensure only one exists or adjust pattern matching.")
}

pairs_path <- pairs_files[1]

# Set default significance threshold if not provided
if (length(args) >= 2) {
    fdr_threshold <- as.numeric(args[2])
} else {
    fdr_threshold <- 0.05
}

if (!file.exists(pairs_path)) {
  stop("Input file not found: ", pairs_path)
}
if (is.na(fdr_threshold) || fdr_threshold < 0 || fdr_threshold > 1) {
  stop("FDR threshold must be a numeric value between 0 and 1")
}

# Read exactly as written in the TSV
pairs <- tryCatch(
  utils::read.table(pairs_path,
                    header = TRUE,
                    sep = "\t",
                    check.names = FALSE,
                    stringsAsFactors = FALSE),
  error = function(e) stop("Failed to read pairs TSV: ", e$message)
)

# Define the p-value column name from Coinfinder output
pcol <- "p"
if (!pcol %in% names(pairs)) {
  stop("Expected p-value column '", pcol, "' not found in: ", pairs_path)
}

# Check for existing FDR_BH column
if ("FDR_BH" %in% names(pairs)) {
  message("Column 'FDR_BH' already exists in: ", pairs_path,
       " — multiple testing correction already applied.")
  quit(status=0)
}

# Apply BH correction
fdr_values <- p.adjust(pairs[[pcol]], method = "BH")

# Insert FDR_BH directly after the p-value column
p_index <- which(names(pairs) == pcol)
pairs <- data.frame(
  pairs[1:p_index],
  FDR_BH = fdr_values,
  pairs[(p_index + 1):ncol(pairs)],
  check.names = FALSE
)

# Filter rows according to the passed threshold
pairs <- pairs[pairs$FDR_BH <= fdr_threshold, , drop = FALSE]

# Overwrite the original file
utils::write.table(pairs,
                   pairs_path,
                   sep = "\t",
                   quote = FALSE,
                   row.names = FALSE)

message("FDR_BH column inserted and filtered at threshold ", fdr_threshold,
        ". Updated file: ", pairs_path,
        " | Remaining rows: ", nrow(pairs))
