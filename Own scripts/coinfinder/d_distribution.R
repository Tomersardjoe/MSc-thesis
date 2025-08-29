#!/usr/bin/env Rscript

library(ggplot2)

args <- commandArgs(trailingOnly = TRUE)

if (length(args) == 0) {
  stop("Please provide the path to a coincident_nodes.tsv file.
       \nUsage: Rscript plot_distribution.R path to coincident_nodes.tsv, e.g., coinfinder/uniqueID/coincident_nodes.tsv")
}

file_path <- args[1]

# Extract unique ID from parent directory
unique_id <- basename(dirname(file_path))

# Read the TSV file
data <- read.table(file_path, 
                   header = TRUE, 
                   sep = "\t", 
                   stringsAsFactors = FALSE)

# Ensure D-value is numeric
data$Result <- as.numeric(data$Result)

# Calculate summary stats
med_val <- median(data$Result, na.rm = TRUE)
q25 <- quantile(data$Result, 0.25, na.rm = TRUE)
q75 <- quantile(data$Result, 0.75, na.rm = TRUE)

# Create plot
p <- ggplot(data, aes(x = Result)) +
  geom_histogram(binwidth = 0.1, 
                 fill = "#69b3a2", 
                 color = "black", 
                 alpha = 0.7) +
  geom_vline(xintercept = med_val, color = "red", linetype = "solid", linewidth = 1) +
  geom_vline(xintercept = q25, color = "blue", linetype = "dashed", linewidth = 1) +
  geom_vline(xintercept = q75, color = "blue", linetype = "dashed", linewidth = 1) +
  annotate("text", x = med_val, y = Inf, label = "Median", vjust = -0.5, color = "red") +
  annotate("text", x = q25, y = Inf, label = "Q1", vjust = -0.5, color = "blue") +
  annotate("text", x = q75, y = Inf, label = "Q3", vjust = -0.5, color = "blue") +
  labs(title = paste("Distribution of Result Values —", unique_id),
       x = "Result",
       y = "Count") +
  theme_minimal()

# Build output path in the same directory as the TSV
output_path <- file.path(dirname(file_path), paste0("d_distribution_", unique_id, ".png"))

# Save plot
ggsave(output_path, plot = p, width = 8, height = 6)

message("Plot saved as: ", output_path)
