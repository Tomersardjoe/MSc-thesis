#!/usr/bin/env Rscript

# Get command-line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Usage: Rscript prep_d_calc.R <collapsed_matrix.csv>")
}

input_file  <- args[1]

run_id <- basename(dirname(input_file))

nodes_in_file <- paste0(run_id, "_nodes_in.csv")
clean_matrix_file <- paste0(run_id, "_collapsed_matrix_clean.csv")

# Read the collapsed presence/absence matrix (header=TRUE skips the header row)
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
