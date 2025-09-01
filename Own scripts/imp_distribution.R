#!/usr/bin/env Rscript

library(ggplot2)
library(dplyr)

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 1) {
  stop("Usage: Rscript imp_distribution.R /path/to/imp.csv")
}

imp_path <- args[1]
imp_dir <- dirname(imp_path)

# 1. Read the CSV
imp <- read.csv(imp_path, row.names = 1, check.names = FALSE)

# 2. Convert to matrix
imp_mat <- as.matrix(imp)

# 3. Remove the diagonal
diag(imp_mat) <- NA

# 4. Flatten and drop NAs
scores <- as.numeric(imp_mat)
scores <- scores[!is.na(scores)]

# Remove zeros
scores_nozero <- scores[scores > 0]

# Summary stats
med <- median(scores_nozero)
q1  <- quantile(scores_nozero, 0.25)
q3  <- quantile(scores_nozero, 0.75)

# Build the plot object
p <- ggplot(data.frame(score = scores_nozero), aes(x = score)) +
  geom_histogram(bins = 30, fill = "steelblue", color = "white") +
  scale_x_log10() +
  geom_vline(xintercept = med, color = "red", linetype = "solid", size = 1) +
  geom_vline(xintercept = q1, color = "orange", linetype = "dashed", size = 1) +
  geom_vline(xintercept = q3, color = "orange", linetype = "dashed", size = 1) +
  annotate("text", x = med, y = Inf, label = paste0("Median: ", signif(med, 3)),
           vjust = -0.5, hjust = -0.1, color = "red") +
  annotate("text", x = q1, y = Inf, label = paste0("Q1: ", signif(q1, 3)),
           vjust = -0.5, hjust = -0.1, color = "orange") +
  annotate("text", x = q3, y = Inf, label = paste0("Q3: ", signif(q3, 3)),
           vjust = -0.5, hjust = -0.1, color = "orange") +
  labs(
    title = "Log-Scaled Histogram with Median and Quantiles",
    x = "Importance Score (log10 scale)",
    y = "Frequency"
  ) +
  theme_minimal()

# Save the plot in the same directory as imp.csv
ggsave(
  filename = file.path(imp_dir, "importance_histogram.png"),
  plot = p,
  width = 8,
  height = 6,
  dpi = 300
)
