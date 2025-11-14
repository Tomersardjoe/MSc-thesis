#!/usr/bin/env Rscript

# Get command-line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Usage: Rscript prep_d_calc.R <collapsed_matrix.csv>")
}

input_file <- args[1]
input_dir  <- dirname(normalizePath(input_file))
run_id     <- basename(input_dir)

# Define output directory inside the run directory
out_dir <- file.path(input_dir, "imp_cutoff")
if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

nodes_in_file     <- file.path(out_dir, paste0(run_id, "_nodes_in.csv"))
clean_matrix_file <- file.path(out_dir, paste0(run_id, "_collapsed_matrix_clean.csv"))

# Read the collapsed presence/absence matrix
collapsed <- read.csv(input_file, header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)

# Remove any columns without headers
empty_cols <- which(colnames(collapsed) == "")
if (length(empty_cols) > 0) {
  message("Removing ", length(empty_cols), " column(s) with empty headers: ", paste(empty_cols, collapse = ", "))
  collapsed <- collapsed[, colnames(collapsed) != ""]
} else {
  message("No empty header columns found.")
}

# Save the cleaned matrix
write.csv(collapsed, clean_matrix_file, row.names = FALSE, quote = FALSE)
message("Wrote cleaned matrix to: ", clean_matrix_file)

# Grab only the first column values (excluding the header)
ids <- collapsed[[1]]

# Create a one-row data frame with IDs as column names
nodes_in <- as.data.frame(matrix(nrow = 1, ncol = length(ids)))
colnames(nodes_in) <- ids
write.csv(nodes_in, nodes_in_file, row.names = FALSE, quote = FALSE)

cat("Wrote", length(ids), "IDs to", nodes_in_file, "\n")
