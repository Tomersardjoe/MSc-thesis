# Load required library
library(ggplot2)
library(dplyr)

# 1. Read the CSV, treating first column as row names if it's a labeled matrix
imp <- read.csv("imp.csv", row.names = 1, check.names = FALSE)

# 2. Convert to matrix (ensures numeric handling)
imp_mat <- as.matrix(imp)

# 3. Remove the diagonal (self-importance scores)
#    diag(imp_mat) <- NA replaces diagonal with NA
diag(imp_mat) <- NA

# 4. Flatten to a vector and drop NAs
scores <- as.numeric(imp_mat)
scores <- scores[!is.na(scores)]

# Remove zeros entirely
scores_nozero <- scores[scores > 0]

# Calculate summary stats
med <- median(scores_nozero)
q1  <- quantile(scores_nozero, 0.25)
q3  <- quantile(scores_nozero, 0.75)

# Plot
ggplot(data.frame(score = scores_nozero), aes(x = score)) +
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

# Display separately
print(p_high)
print(p_full_log)

      