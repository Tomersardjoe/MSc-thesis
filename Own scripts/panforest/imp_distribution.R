#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(ggrepel)
  library(dplyr)
  library(patchwork)
  library(readr)
})

safe_approx <- function(x, y, xout) {
  keep <- !is.na(x) & !is.na(y)
  x <- x[keep]; y <- y[keep]
  ord <- order(x)
  x <- x[ord]; y <- y[ord]
  if (any(duplicated(x))) {
    sel <- !duplicated(x)
    x <- x[sel]; y <- y[sel]
  }
  if (any(duplicated(xout))) {
    idx   <- ave(seq_along(xout), xout, FUN = seq_along)
    xout  <- xout + (idx - 1) * 1e-9
  }
  approx(x, y, xout = xout, rule = 2)$y
}

make_symmetric_mean <- function(mat) {
  mat <- as.matrix(mat)
  rn <- rownames(mat)
  diag(mat) <- NA
  sym_mat <- (mat + t(mat)) / 2
  na_idx  <- is.na(sym_mat)
  sym_mat[na_idx] <- mat[na_idx] + t(mat)[na_idx]
  rownames(sym_mat) <- rn
  colnames(sym_mat) <- rn
  sym_mat
}

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {stop("Usage: Rscript imp_distribution.R /path/to/panforest/{run_id}/imp_fixed.csv /path/to/panforest/{run_id}_nodes.tsv")} # NOT THE COINFINDER nodes_all.tsv but the output of PanForest's calculate_d.R

imp_path <- args[1]
dval_path <- args[2]
imp_dir  <- dirname(imp_path)

unique_id <- basename(dirname(imp_path))

imp            <- read.csv(imp_path, row.names = 1, check.names = FALSE)
imp_mat_sym    <- make_symmetric_mean(as.matrix(imp))
scores         <- as.numeric(imp_mat_sym)
scores         <- scores[!is.na(scores) & scores > 0]
sorted_scores  <- sort(scores, decreasing = TRUE)

n              <- length(sorted_scores)
all_points     <- cbind(1:n, sorted_scores)
first_point    <- matrix(all_points[1, ], nrow = 1)
last_point     <- matrix(all_points[n, ], nrow = 1)
line_vec       <- last_point - first_point
line_vec       <- line_vec / sqrt(sum(line_vec^2))
vec_from_first <- sweep(all_points, 2, first_point)
scalar_proj    <- vec_from_first %*% t(line_vec)
proj_point     <- scalar_proj %*% line_vec + matrix(rep(first_point, n), nrow = n, byrow = TRUE)
distances      <- sqrt(rowSums((all_points - proj_point)^2))
elbow_idx      <- which.max(distances)
cutoff_value   <- sorted_scores[elbow_idx]

# Save cutoff_value to a text file
cutoff_file <- file.path(imp_dir, paste0(unique_id, "_cutoff_value.txt"))
write(cutoff_value, file = cutoff_file)

summary_stats <- data.frame(
  Stat  = c("Min", "Q1", "Median", "Mean", "Q3", "Max", "N", "ElbowCutoff"),
  Value = c(
    min(sorted_scores),
    quantile(sorted_scores, 0.25),
    median(sorted_scores),
    mean(sorted_scores),
    quantile(sorted_scores, 0.75),
    max(sorted_scores),
    n,
    cutoff_value
  )
)

stats_file <- file.path(imp_dir, "hist_stats.csv")
write.csv(summary_stats, stats_file, row.names = FALSE)
stats_vals <- setNames(summary_stats$Value, summary_stats$Stat)

ref_lines <- data.frame(
  Stat  = c("Q1", "Median", "Q3", "Elbow"),
  Value = c(stats_vals["Q1"], stats_vals["Median"], stats_vals["Q3"], stats_vals["ElbowCutoff"])
)

stat_colors   <- c("Q1" = "#1f78b4", "Median" = "#33a02c", "Q3" = "#1f78b4", "Elbow" = "#e31a1c")
stat_linetypes <- c("Q1" = "dashed", "Median" = "solid", "Q3" = "dashed", "Elbow" = "dotdash")

df_sorted <- data.frame(score = log10(sorted_scores))
dens_obj  <- density(df_sorted$score)

ref_lines_density <- ref_lines %>%
  mutate(
    log_value = log10(Value),
    y         = safe_approx(dens_obj$x, dens_obj$y, log_value)
  )

pct_above <- mean(sorted_scores >= cutoff_value) * 100

p_density <- ggplot(df_sorted, aes(x = score)) +
  geom_density(fill = "#6baed6", alpha = 0.4, color = "#08519c", linewidth = 1) +
  geom_vline(data = ref_lines_density,
             aes(xintercept = log_value, color = Stat, linetype = Stat),
             linewidth = 1) +
  annotate("point",
           x = log10(cutoff_value),
           y = safe_approx(dens_obj$x, dens_obj$y, log10(cutoff_value)),
           color = "#e31a1c", size = 3) +
  annotate("rect",
           xmin = log10(cutoff_value), xmax = Inf, ymin = 0, ymax = Inf,
           alpha = 0.1, fill = "#e31a1c") +
  geom_label_repel(
    data = ref_lines_density,
    aes(
      x     = log_value,
      y     = y,
      label = paste0(Stat, "\n", format(Value, scientific = TRUE, digits = 3)),
      fill  = I("white"),
      color = Stat
    ),
    inherit.aes   = FALSE,
    nudge_y       = 0.05,
    nudge_x       = c(-0.25,0,0.25,0.25),
    size          = 3.5,
    max.overlaps  = Inf,
    segment.color = "grey50",
    box.padding   = 0.4,
    label.size    = 0
  ) +
  scale_color_manual(values = stat_colors) +
  scale_linetype_manual(values = stat_linetypes) +
  labs(
    title    = "Log-scaled density of symmetric gene-gene importance scores\nwith elbow & quartiles",
    subtitle = paste0("Elbow at ", round(cutoff_value, 3),
                      " (", round(pct_above, 1), "% of unique gene pairs above cutoff)"),
    x        = "Importance Score (log10 scale)",
    y        = "Density"
  ) +
  theme_minimal() +
  theme(legend.position = "none")

df_cdf  <- data.frame(score = sort(sorted_scores),
                      cum_prop = seq_along(sorted_scores) / length(sorted_scores))
elbow_y <- safe_approx(df_cdf$score, df_cdf$cum_prop, cutoff_value)

p_cdf <- ggplot(df_cdf, aes(x = score, y = cum_prop)) +
  geom_line(color = "#08519c", linewidth = 1) +
  geom_vline(data = subset(ref_lines, Stat == "Elbow"),
             aes(xintercept = Value, color = Stat, linetype = Stat),
             linewidth = 1) +
  annotate("point",
           x = cutoff_value, y = elbow_y,
           color = "#e31a1c", size = 3) +
  geom_label_repel(
    data = data.frame(Stat = "Elbow", Value = cutoff_value, y = elbow_y),
    aes(
      x     = Value,
      y     = y,
      label = paste0("Elbow\n", format(Value, scientific = TRUE, digits = 3)),
      fill  = I("white"),
      color = Stat
    ),
    inherit.aes   = FALSE,
    nudge_y       = -0.05,
    nudge_x       = 0.005,
    size          = 3.5,
    segment.color = "grey50",
    box.padding   = 0.4,
    label.size    = 0
  ) +
  scale_color_manual(values = stat_colors) +
  scale_linetype_manual(values = stat_linetypes) +
  labs(
    title    = "Cumulative distribution of symmetric gene-gene importance\nscores with elbow",
    subtitle = paste0(round(pct_above, 1), "% of unique gene pairs above cutoff"),
    x        = "Importance Score",
    y        = "Cumulative Proportion"
  ) +
  theme_minimal() +
  theme(legend.position = "none")

combined <- p_density + p_cdf

output_file_combined <- file.path(imp_dir, "importance_two_panel.png")
ggsave(output_file_combined, plot = combined, width = 12, height = 6, dpi = 300)

# Per-gene summaries
diag(imp_mat_sym) <- NA
strong_mat <- imp_mat_sym >= cutoff_value

gene_counts <- rowSums(strong_mat, na.rm = TRUE)
gene_mean_scores <- apply(imp_mat_sym, 1, function(x) {
  vals <- x[!is.na(x) & x >= cutoff_value]
  if (length(vals) == 0) return(NA)
  mean(vals)
})

# Data frames for plotting
df_counts <- data.frame(
  metric = "Strong partner count",
  value  = gene_counts
)

df_means <- data.frame(
  metric = "Mean strong score",
  value  = gene_mean_scores
)

# Histograms
p_counts <- ggplot(df_counts, aes(x = value)) +
  geom_histogram(binwidth = 1, fill = "#6baed6", color = "#08519c", alpha = 0.6) +
  labs(
    title = "Distribution of strong partner counts per gene",
    subtitle = paste0("Cutoff = ", round(cutoff_value, 3)),
    x = "Number of partners >= cutoff",
    y = "Gene count"
  ) +
  theme_minimal()

p_means <- ggplot(df_means, aes(x = value)) +
  geom_histogram(fill = "#fc9272", color = "#cb181d", alpha = 0.6, bins = 30) +
  labs(
    title = "Distribution of mean strong scores per gene",
    subtitle = paste0("Cutoff = ", round(cutoff_value, 3)),
    x = "Mean importance score (>= cutoff)",
    y = "Gene count"
  ) +
  theme_minimal()

# Combine into two-panel figure
combined_gene_plots <- p_counts + p_means

output_file_genes <- file.path(imp_dir, "per_gene_histograms.png")
ggsave(output_file_genes, plot = combined_gene_plots, width = 12, height = 6, dpi = 300)

# Read D-values (genes)
dvals <- read.delim(dval_path, header = TRUE, stringsAsFactors = FALSE)
if (!all(c("ID", "Result") %in% colnames(dvals))) {
  stop("D-value file must have columns: ID, Result")
}

# Map gene -> D-value
dval_map <- setNames(dvals$Result, dvals$ID)

# Keep only the union between the D-value file and the importance matrix genes
gene_names <- rownames(strong_mat)
common_genes <- intersect(gene_names, names(dval_map))

if (length(common_genes) < length(gene_names)) {
  warning("Filtered out ", length(gene_names) - length(common_genes),
          " genes not found in D-value file before plotting.")
}

upper_idx <- upper.tri(strong_mat, diag = FALSE)
gene_names <- rownames(strong_mat)

if (is.null(gene_names)) stop("No row names found in importance matrix - check your input CSV.")

missing_genes <- setdiff(gene_names, names(dval_map))
if (length(missing_genes) > 0) {
  warning("Some genes in importance matrix not found in D-value file: ",
          paste(head(missing_genes, 10), collapse = ", "),
          if (length(missing_genes) > 10) " ...")
}

# Histogram 1: all associating genes vs D-value
p_dval_genes <- ggplot(
  data.frame(Dvalue = dval_map[gene_names]),
  aes(x = Dvalue)
) +
  geom_histogram(fill = "#6baed6", color = "#08519c", alpha = 0.6, bins = 50) +
  labs(
    title = "Associating genes vs D-value",
    x = "D-value",
    y = "Number of genes"
  ) +
  theme_minimal()

# Histogram 2: all associating gene pairs vs D-value

pair_list <- which(strong_mat & upper_idx, arr.ind = TRUE)

if (nrow(pair_list) > 0) {
  g1 <- gene_names[pair_list[, 1]]
  g2 <- gene_names[pair_list[, 2]]
  df_pair_dvals <- data.frame(
    Dvalue = c(dval_map[g1], dval_map[g2]),
    stringsAsFactors = FALSE
  )
} else {
  warning("No strong gene pairs found above cutoff - pair D-value histogram will be empty.")
  df_pair_dvals <- data.frame(Dvalue = numeric(0))
}

p_dval_pairs <- ggplot(df_pair_dvals, aes(x = Dvalue)) +
  geom_histogram(fill = "#fc9272", color = "#cb181d", alpha = 0.6, bins = 50) +
  labs(
    title = "Associating gene pairs vs D-value",
    subtitle = "Each gene in a pair contributes its own D-value",
    x = "D-value",
    y = "Number of genes in pairs"
  ) +
  theme_minimal()

# Combine into two-panel figure
combined_dval_plots <- p_dval_genes + p_dval_pairs
output_file_dvals <- file.path(imp_dir, "dvalue_histograms.png")
ggsave(output_file_dvals, plot = combined_dval_plots, width = 12, height = 6, dpi = 300)

cat("Elbow cutoff figures saved to:", output_file_combined, "\n")
cat("Per-gene histograms saved to:", output_file_genes, "\n")
cat("D-value histograms saved to:", output_file_dvals, "\n")
cat("Summary stats saved to:", stats_file, "\n")

# -------------------------
# Write output data for comparative plots
# -------------------------

# Gene-level data
df_genes <- data.frame(
  Gene     = gene_names,
  Gene_1   = NA_character_,
  Gene_2   = NA_character_,
  D_value  = dval_map[gene_names],
  Method   = "PanForest",
  Level    = "Gene",
  stringsAsFactors = FALSE
) %>%
  filter(!is.na(D_value))

# Pair-level data
if (nrow(pair_list) > 0) {
  g1 <- gene_names[pair_list[, 1]]
  g2 <- gene_names[pair_list[, 2]]

  # For pairs, take the min D-value of the two genes
  d_pair <- pmin(dval_map[g1], dval_map[g2], na.rm = TRUE)

  df_pairs <- data.frame(
    Gene     = NA_character_,
    Gene_1   = g1,
    Gene_2   = g2,
    D_value  = d_pair,
    Method   = "PanForest",
    Level    = "Pair",
    stringsAsFactors = FALSE
  ) %>%
    filter(!is.na(D_value))
} else {
  df_pairs <- data.frame(
    Gene     = NA_character_,
    Gene_1   = character(0),
    Gene_2   = character(0),
    D_value  = numeric(0),
    Method   = character(0),
    Level    = character(0),
    stringsAsFactors = FALSE
  )
}

# Combine
df_out <- bind_rows(df_genes, df_pairs)

# Write to CSV
write_csv(df_out, file.path(imp_dir, paste0("panforest_dvalues_", unique_id, ".csv")))

message("D-value CSV saved to: ",
        file.path(imp_dir, paste0("panforest_dvalues_", unique_id, ".csv")))
