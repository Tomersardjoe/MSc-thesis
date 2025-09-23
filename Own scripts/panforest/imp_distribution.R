#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(ggrepel)
  library(dplyr)
  library(patchwork)
  library(readr)
})

# Helper functions

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

make_gene_summary <- function(values, metric, genes, filter_genes = NULL) {
  df <- data.frame(Gene = genes, value = values, metric = metric, stringsAsFactors = FALSE)
  if (!is.null(filter_genes)) df <- df[df$Gene %in% filter_genes, ]
  df
}

# Read input arguments

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {stop("Usage: Rscript imp_distribution.R /path/to/panforest/{run_id}/imp_cutoff/imp_fixed.csv /path/to/panforest/imp_cutoff/{run_id}_nodes.tsv [D-value cutoff]")} # NOT THE COINFINDER nodes_all.tsv but the output of PanForest's calculate_d.R

imp_path <- args[1]
dval_path <- args[2]
dcutoff_value <- suppressWarnings(as.numeric(args[3]))

imp_dir  <- dirname(imp_path)
out_dir  <- imp_dir
unique_id <- basename(dirname(imp_dir))

# Create symmetric importance matrix and sort scores

imp            <- read.csv(imp_path, row.names = 1, check.names = FALSE)
imp_mat_sym    <- make_symmetric_mean(as.matrix(imp))
scores         <- as.numeric(imp_mat_sym)
scores         <- scores[!is.na(scores) & scores > 0]
sorted_scores  <- sort(scores, decreasing = TRUE)

n <- length(sorted_scores)

# -------------------------
# Restricted distance-to-line elbow (ignore the pre-steep region)
# -------------------------
# Step 1: compute slopes
slopes <- diff(sorted_scores)
steep_idx <- which.min(slopes) + 1   # index of steepest drop

# Step 2: build line from steepest point to last point
all_points <- cbind(1:n, sorted_scores)
first_point <- all_points[steep_idx, , drop = FALSE]
last_point  <- all_points[n, , drop = FALSE]
line_vec <- last_point - first_point
line_vec <- line_vec / sqrt(sum(line_vec^2))

# Step 3: project all points onto that line
vec_from_first <- sweep(all_points, 2, first_point)
scalar_proj    <- vec_from_first %*% t(line_vec)
proj_point     <- scalar_proj %*% line_vec + 
                  matrix(rep(first_point, n), nrow = n, byrow = TRUE)

# Step 4: compute distances, ignoring early section
distances <- sqrt(rowSums((all_points - proj_point)^2))
distances[1:steep_idx] <- NA   # ignore shallow section

# Step 5: elbow = max distance in post-steep region
elbow_idx    <- which.max(distances)
cutoff_value <- sorted_scores[elbow_idx]

# Save cutoff_value to a text file
cutoff_file <- file.path(out_dir, paste0(unique_id, "_cutoff_value.txt"))
write(cutoff_value, file = cutoff_file)

cat("Restricted distance-based cutoff written to:", cutoff_file, "\n")

# Plot elbow value cutoff on density plot

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

stats_file <- file.path(out_dir, "hist_stats.csv")
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

# Plot elbow value on cumulative distribution plot

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

output_file_combined <- file.path(out_dir, "importance_two_panel.png")
ggsave(output_file_combined, plot = combined, width = 12, height = 6, dpi = 300)

# -------------------------------------------------
# Stop here if *_dcutoff.txt file does not exist
# -------------------------------------------------
if (is.na(dcutoff_value)) {
  cat("No D-value cutoff provided, stopping after importance score cutoff calculation.\n")
  quit(save = "no", status = 0)
}

# Combined filter for D-value cutoff and importance score cutoff
diag(imp_mat_sym) <- NA
strong_mat <- imp_mat_sym >= cutoff_value
gene_names <- rownames(strong_mat)
if (is.null(gene_names)) stop("No row names found in importance matrix - check your input CSV.")

upper_idx <- upper.tri(strong_mat, diag = FALSE)
pair_list <- which(upper_idx, arr.ind = TRUE)

# Read D-values
dvals <- read.delim(dval_path, header = TRUE, stringsAsFactors = FALSE)
if (!all(c("ID", "Result") %in% colnames(dvals))) {
  stop("D-value file must have columns: ID, Result")
}
dval_map <- setNames(dvals$Result, dvals$ID)

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

# Importance cutoff filter
genes_imp_filtered <- gene_names[rowSums(strong_mat, na.rm = TRUE) > 0]

# Both importance AND D-value cutoff filter
genes_both_filtered <- genes_imp_filtered[
  !is.na(dval_map[genes_imp_filtered]) &
  dval_map[genes_imp_filtered] >= dcutoff_value
]

# Per-gene summaries
gene_counts <- rowSums(strong_mat, na.rm = TRUE)
gene_mean_scores <- apply(imp_mat_sym, 1, function(x) {
  vals <- x[!is.na(x) & x >= cutoff_value]
  if (length(vals) == 0) return(NA)
  mean(vals)
})

# Summaries
df_counts <- make_gene_summary(gene_counts, "Partner count", gene_names, genes_both_filtered)
df_means  <- make_gene_summary(gene_mean_scores, "Mean importance score", gene_names, genes_both_filtered)

# Partner and mean importance score distribution histograms
# Distribution of partner counts (all genes)
p_counts <- ggplot(df_counts, aes(x = value)) +
  geom_histogram(binwidth = 1, fill = "#6baed6", color = "#08519c", alpha = 0.6) +
  labs(
    title = "Distribution of gene partners",
    subtitle = paste0("Cutoff = ", round(cutoff_value, 3)),
    x = "Number of partners >= cutoff",
    y = "Gene count"
  ) +
  theme_minimal()

# Distribution of mean importance scores (only genes above D-value cutoff)
p_means <- ggplot(df_means, aes(x = value)) +
  geom_histogram(fill = "#fc9272", color = "#cb181d", alpha = 0.6, bins = 30) +
  labs(
    title = "Distribution of mean importance scores (D-value filtered)",
    subtitle = paste0("D-value cutoff = ", round(dcutoff_value, 3)),
    x = "Mean importance score >= cutoff",
    y = "Gene count"
  ) +
  theme_minimal()

# Combine into two-panel figure
combined_gene_plots <- p_counts + p_means
ggsave(file.path(out_dir, "per_gene_histograms.png"),
       plot = combined_gene_plots, width = 12, height = 6, dpi = 300)

# -------------------------
# Histogram: gene-level D-values
# -------------------------

# Define breaks for both gene-level and pair-level histograms
dval_breaks <- seq(
  floor(min(df_genes$D_value, na.rm = TRUE)),
  ceiling(max(df_genes$D_value, na.rm = TRUE)),
  by = 0.25
)

# Gene-level data for importance cutoff only
df_genes_imp <- data.frame(
  Gene   = genes_imp_filtered,
  D_value = dval_map[genes_imp_filtered],
  stringsAsFactors = FALSE
)

p_hist_genes <- ggplot(df_genes_imp, aes(x = D_value)) +
  geom_histogram(
    binwidth = 0.25, boundary = dcutoff_value,
    fill = "#6baed6", color = "#08519c", alpha = 0.6
  ) +
  geom_vline(
    xintercept = dcutoff_value,
    linetype = "dashed", color = "#e31a1c", linewidth = 1
  ) +
  scale_x_continuous(
    breaks = dval_breaks,
    labels = ifelse(dval_breaks %% 0.5 == 0, dval_breaks, "")
  ) +
  labs(
    title = paste("Distribution of gene D-values -", unique_id),
    subtitle = paste0("D-value cutoff = ", round(dcutoff_value, 3)),
    x = "D-value", y = "Gene count"
  ) +
  theme_minimal()

# -------------------------
# Histogram: pair-level D-values
# -------------------------

p_hist_pairs <- ggplot(df_pairs, aes(x = d_pair)) +
  geom_histogram(
    binwidth = 0.25, boundary = dcutoff_value,
    fill = "#fdae6b", color = "#e6550d", alpha = 0.6
  ) +
  geom_vline(
    xintercept = dcutoff_value,
    linetype = "dashed", color = "#e31a1c", linewidth = 1
  ) +
  scale_x_continuous(
    breaks = dval_breaks,
    labels = ifelse(dval_breaks %% 0.5 == 0, dval_breaks, "")
  ) +
  labs(
    title = paste("Distribution of gene pairs D-values -", unique_id),
    subtitle = paste0("D-value cutoff = ", round(dcutoff_value, 3)),
    x = "D-value (min of source & target)",
    y = "Pair count"
  ) +
  theme_minimal()

# Combine into two-panel figure
combined_dval_plots <- p_hist_genes + p_hist_pairs
output_file_dvals <- file.path(out_dir, "dvalue_histograms.png")
ggsave(output_file_dvals, plot = combined_dval_plots, width = 12, height = 6, dpi = 300)

# -------------------------
# Write output data
# -------------------------

# Combine
df_out <- bind_rows(df_genes, df_pairs)

# Write to CSV
write_csv(df_out, file.path(out_dir, paste0("panforest_dvalues_", unique_id, ".csv")))

message("D-value CSV saved to: ",
        file.path(out_dir, paste0("panforest_dvalues_", unique_id, ".csv")))
