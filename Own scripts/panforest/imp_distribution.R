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

# Elbow finder accepts numeric vector or a data.frame with columns cutoff and pairs
elbow_cutoff <- function(x, na_ok = FALSE) {
  if (is.null(x) || length(x) == 0) {
    if (na_ok) return(NA_real_)
    stop("Input is empty and na_ok = FALSE")
  }

  # Normalize input
  if (is.numeric(x) && is.null(dim(x))) {
    # Case 1: numeric vector of sorted scores
    values <- x
    cutoffs <- x
  } else if (is.data.frame(x) && all(c("cutoff", "pairs") %in% names(x))) {
    # Case 2: data.frame with cutoff/pairs
    values <- x$pairs
    cutoffs <- x$cutoff
  } else {
    stop("Unsupported input type: must be numeric vector or data.frame with cutoff/pairs")
  }

  n <- length(values)

  # Step 1: compute slopes
  slopes <- diff(values)
  steep_idx <- which.min(slopes) + 1

  # Step 2: build line from steepest point to last point
  all_points <- cbind(1:n, values)
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
  distances[1:steep_idx] <- NA

  # Step 5: elbow = max distance in post-steep region
  elbow_idx    <- which.max(distances)
  imp_cutoff <- cutoffs[elbow_idx]

  return(imp_cutoff)
}

# Read input arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {stop("Usage: Rscript imp_distribution.R /path/to/panforest/{run_id}/imp_cutoff/imp_fixed.csv /path/to/panforest/{run_id}/performance.csv /path/to/panforest/imp_cutoff/{run_id}_nodes.tsv [D-value cutoff]")} # DO NOT USE THE COINFINDER nodes_all.tsv but the output of PanForest's calculate_d.R

imp_path      <- args[1]
perf_path     <- args[2]
dval_path     <- args[3]
dcutoff_value <- suppressWarnings(as.numeric(args[4]))

imp_dir   <- dirname(imp_path)
out_dir   <- imp_dir
unique_id <- basename(dirname(imp_dir))

# Parameter defaults
imp_cutoff_default <- 0.01
accuracy_cutoff    <- 0.9
f1_cutoff          <- 0.9

# Read importance matrix
imp         <- read.csv(imp_path, row.names = 1, check.names = FALSE)
imp_sym_raw <- make_symmetric_mean(as.matrix(imp))
gene_names  <- rownames(imp_sym_raw)

# Read D-values
dvals <- read.delim(dval_path, header = TRUE, stringsAsFactors = FALSE)
if (!all(c("ID", "Result") %in% colnames(dvals))) {
  stop("D-value file must have columns: ID, Result")
}
dval_map <- setNames(dvals$Result, dvals$ID)

# Apply D-value cutoff
keep_genes <- gene_names[!is.na(dval_map[gene_names]) & dval_map[gene_names] >= dcutoff_value]
lost_genes <- setdiff(gene_names, keep_genes)

imp_sym <- imp_sym_raw[keep_genes, keep_genes, drop = FALSE]

message("Initial matrix: ", length(gene_names), " genes, ",
        sum(upper.tri(imp_sym_raw)), " gene pairs")
        
message("After D-value cutoff (", dcutoff_value, "): ",
        length(keep_genes), " genes retained (",
        length(lost_genes), " lost)")

pairs_after_d <- sum(upper.tri(imp_sym))
message("Gene pairs retained after D cutoff: ", pairs_after_d)

# Importance score cutoff
imp_cutoff <- imp_cutoff_default
before_pairs <- sum(!is.na(imp_sym[upper.tri(imp_sym)]))
imp_sym[imp_sym < imp_cutoff] <- NA
after_pairs <- sum(!is.na(imp_sym[upper.tri(imp_sym)]))

message("Importance cutoff (", imp_cutoff, "): ",
        after_pairs, " gene pairs retained (",
        before_pairs - after_pairs, " lost)")

# Accuracy and F1-score cutoff
perf <- read.csv(perf_path, header = TRUE, stringsAsFactors = FALSE, na.strings = c("NA", ""))
# Expect columns: Accuracy (Ate), F1 (F1te)
perf_colnames <- colnames(perf)
if (!all(c("Ate", "F1te") %in% perf_colnames)) {
  stop("Performance file must contain columns: Ate and F1te.")
}

id_col <- perf_colnames[1]
acc_col <- "Ate"
f1_col  <- "F1te"

perf_map_acc <- setNames(as.numeric(perf[[acc_col]]), perf[[id_col]])
perf_map_f1  <- setNames(as.numeric(perf[[f1_col]]), perf[[id_col]])

gene_acc <- perf_map_acc[rownames(imp_sym)]
gene_f1  <- perf_map_f1[rownames(imp_sym)]

gene_acc <- as.numeric(gene_acc)
gene_f1  <- as.numeric(gene_f1)

pass_gene <- !is.na(gene_acc) & !is.na(gene_f1) & (gene_acc >= accuracy_cutoff) & (gene_f1 >= f1_cutoff)

n_genes_before <- nrow(imp_sym)
n_genes_pass   <- sum(pass_gene)
n_genes_fail   <- n_genes_before - n_genes_pass
message("Genes retained after per-gene accuracy/F1-score cutoffs: ", n_genes_pass, " and ", n_genes_fail, " lost)")

# For any gene that fails, set its entire row/column in imp_sym to NA
failed_genes <- rownames(imp_sym)[!pass_gene]

if (length(failed_genes) > 0) {
  imp_sym[failed_genes, ] <- NA
  imp_sym[, failed_genes] <- NA
}

pairs_remaining <- sum(!is.na(imp_sym[upper.tri(imp_sym)]))
message("Gene pairs remaining after per-gene accuracy/F1-score filtering: ", pairs_remaining)

# Prepare scores and stats
scores           <- as.numeric(imp_sym_raw)
scores           <- scores[!is.na(scores) & scores > 0]
sorted_scores    <- sort(scores, decreasing = TRUE)
n_scores <- length(sorted_scores)

# Summary stats data frame

summary_stats <- data.frame(
  Stat  = c("Min", "Q1", "Median", "Mean", "Q3", "Max", "N", "Elbow"),
  Value = c(
    min(sorted_scores),
    quantile(sorted_scores, 0.25),
    median(sorted_scores),
    mean(sorted_scores),
    quantile(sorted_scores, 0.75),
    max(sorted_scores),
    n_scores,
    imp_cutoff
  )
)

stats_vals <- setNames(summary_stats$Value, summary_stats$Stat)

stats_file <- file.path(out_dir, "hist_stats.csv")
write.csv(summary_stats, stats_file, row.names = FALSE)

# Density and CDF plots of importance score distribution
ref_lines <- data.frame(
  Stat  = c("Q1", "Median", "Q3", "Elbow"),
  Value = c(stats_vals["Q1"], stats_vals["Median"], stats_vals["Q3"], stats_vals["Elbow"])
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

# Number retained above cutoff
n_retained <- sum(sorted_scores >= imp_cutoff)

# Percentage retained
pct_above <- if (n_scores > 0) (n_retained / n_scores) * 100 else NA_real_

p_density <- ggplot(df_sorted, aes(x = score)) +
  geom_density(fill = "#6baed6", alpha = 0.4, color = "#08519c", linewidth = 1) +
  geom_vline(data = ref_lines_density,
             aes(xintercept = log_value, color = Stat, linetype = Stat),
             linewidth = 1) +
  annotate("point",
           x = log10(imp_cutoff),
           y = safe_approx(dens_obj$x, dens_obj$y, log10(imp_cutoff)),
           color = "#e31a1c", size = 3) +
  annotate("rect",
           xmin = log10(imp_cutoff), xmax = Inf, ymin = 0, ymax = Inf,
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
    title    = "Log-scaled density of symmetric gene-gene \nimportance scores",
    subtitle = paste0("Importance score cutoff at ", round(imp_cutoff, 3),
                      " (", "Gene pairs retained: ",round(pct_above, 1), "%)"),
    x        = "Importance Score (log10 scale)",
    y        = "Density"
  ) +
  theme_minimal() +
  theme(legend.position = "none")

# Plot elbow value on cumulative distribution plot

df_cdf  <- data.frame(score = sort(sorted_scores),
                      cum_prop = seq_along(sorted_scores) / length(sorted_scores))
elbow_y <- safe_approx(df_cdf$score, df_cdf$cum_prop, imp_cutoff)

p_cdf <- ggplot(df_cdf, aes(x = score, y = cum_prop)) +
  geom_line(color = "#08519c", linewidth = 1) +
  geom_vline(data = subset(ref_lines, Stat == "Elbow"),
             aes(xintercept = Value, color = Stat, linetype = Stat),
             linewidth = 1) +
  annotate("point",
           x = imp_cutoff, y = elbow_y,
           color = "#e31a1c", size = 3) +
  annotate("rect",
           xmin = imp_cutoff, xmax = Inf, ymin = 0, ymax = Inf,
           alpha = 0.1, fill = "#e31a1c") +
  geom_label_repel(
    data = data.frame(Stat = "Elbow", Value = imp_cutoff, y = elbow_y),
    aes(
      x     = Value,
      y     = y,
      label = paste0("Importance score cutoff\n", format(Value, scientific = TRUE, digits = 3)),
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
    title    = "Cumulative distribution of symmetric gene-gene \nimportance scores",
    subtitle = paste0("Cutoff score cutoff ", round(imp_cutoff, 3)),
    x        = "Importance Score",
    y        = "Cumulative Proportion"
  ) +
  theme_minimal() +
  theme(legend.position = "none")

combined <- p_density + p_cdf

output_file_combined <- file.path(out_dir, "importance_two_panel.png")
ggsave(output_file_combined, plot = combined, width = 12, height = 6, dpi = 300)

# Combined filter for D-value cutoff and importance score cutoff
diag(imp_sym) <- NA
strong_mat <- imp_sym >= imp_cutoff
gene_names <- rownames(strong_mat)
if (is.null(gene_names)) stop("No row names found in importance matrix - check: ", imp_path)

upper_idx <- upper.tri(strong_mat, diag = FALSE)
pair_list <- which(upper_idx, arr.ind = TRUE)

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

# Calculate elbow-based cutoff
if (nrow(df_pairs) > 0) {
  cutoffs <- seq(
    min(df_pairs$D_value, na.rm = TRUE),
    max(df_pairs$D_value, na.rm = TRUE),
    length.out = 200
  )

  pair_counts <- sapply(cutoffs, function(cut) {
    sum(df_pairs$D_value >= cut, na.rm = TRUE)
  })

  pair_curve <- data.frame(
    cutoff = cutoffs,
    pairs  = pair_counts
  )
} else {
  pair_curve <- data.frame(cutoff = numeric(0), pairs = numeric(0))
}

# Calculate elbow-based cutoff
elbow_dcutoff <- elbow_cutoff(pair_curve, na_ok = TRUE)

# Compute elbow pairs
pairs_elbow <- sum(df_pairs$D_value >= elbow_dcutoff, na.rm = TRUE)

# Only compute Coinfinder pairs if cutoff provided
if (!is.null(dcutoff_value) && !is.na(dcutoff_value)) {
  pairs_coinfinder <- sum(df_pairs$D_value >= dcutoff_value, na.rm = TRUE)
} else {
  pairs_coinfinder <- NA
}

# Always include the elbow cutoff
cutoff_df <- data.frame(
  Stat   = "Elbow",
  cutoff = elbow_dcutoff
)

# Add Coinfinder cutoff only if it exists
if (!is.na(dcutoff_value)) {
  cutoff_df <- rbind(
    cutoff_df,
    data.frame(
      Stat   = "Coinfinder",
      cutoff = dcutoff_value
    )
  )
}

# Build subtitle for genes
genes_imp_filtered <- gene_names[rowSums(strong_mat, na.rm = TRUE) > 0]

# Genes retained at elbow D-value cutoff
genes_elbow <- sum(
  !is.na(dval_map[genes_imp_filtered]) &
  dval_map[genes_imp_filtered] >= elbow_dcutoff
)

# Genes retained at Coinfinder cutoff (if provided)
if (!is.null(dcutoff_value) && !is.na(dcutoff_value)) {
  genes_coinfinder <- sum(
    !is.na(dval_map[genes_imp_filtered]) &
    dval_map[genes_imp_filtered] >= dcutoff_value
  )
} else {
  genes_coinfinder <- NA
}

# Build subtitle lines dynamically
subtitle_lines_genes <- c(
  paste0("Elbow method cutoff - genes retained: ", genes_elbow)
)
if (!is.na(genes_coinfinder)) {
  subtitle_lines_genes <- c(
    paste0("Q3 cutoff - genes retained: ", genes_coinfinder),
    subtitle_lines_genes
  )
}
subtitle_text_genes <- paste(subtitle_lines_genes, collapse = "\n")

# Build subtitle for pairs
subtitle_lines <- c(
  paste0("Elbow method cutoff - gene pairs retained: ", pairs_elbow)
)
if (!is.na(pairs_coinfinder)) {
  subtitle_lines <- c(
    paste0("Q3 cutoff - gene pairs retained: ", pairs_coinfinder),
    subtitle_lines
  )
}
subtitle_text_pairs <- paste(subtitle_lines, collapse = "\n")

# -------------------------
# Pairs-retained curve with elbow cutoff
# -------------------------

# Build base plot with Elbow
p_curve <- ggplot(pair_curve, aes(x = cutoff, y = pairs)) +
  geom_point(size = 1.5, color = "black") +
  geom_line(color = "#08519c", linewidth = 1) +
  geom_vline(
    data = cutoff_df,
    aes(xintercept = cutoff, color = Stat, linetype = Stat),
    linewidth = 1
  ) +
  scale_color_manual(
    values = c("Coinfinder" = "#450808", "Elbow" = "#e31a1c"),
    breaks = c("Coinfinder", "Elbow"),
    labels = c(
      paste0("Coinfinder (", signif(dcutoff_value, 3), ")"),
      paste0("Elbow method (", signif(elbow_dcutoff, 3), ")")
    ),
    drop = TRUE   # <- important: drop unused levels
  ) +
  scale_linetype_manual(
    values = c("Coinfinder" = "dotted", "Elbow" = "dashed"),
    breaks = c("Coinfinder", "Elbow"),
    labels = c(
      paste0("Coinfinder (", signif(dcutoff_value, 3), ")"),
      paste0("Elbow method (", signif(elbow_dcutoff, 3), ")")
    ),
    drop = TRUE
  ) +
  labs(
    title    = "Pairs retained vs. D-value cutoff",
    subtitle = subtitle_text_pairs,
    x        = "Lineage independence cutoff (D)",
    y        = "# of gene pairs retained",
    color    = "Cutoff type", linetype = "Cutoff type"
  ) +
  coord_cartesian(xlim = c(-5, 5)) +
  theme_minimal() +
  theme(legend.position = "right")

n_low_pairs  <- sum(df_pairs$D_value < -5, na.rm = TRUE)
n_high_pairs <- sum(df_pairs$D_value >  5, na.rm = TRUE)

if (n_low_pairs > 0) {
  p_curve <- p_curve +
    annotate("text", x = -5, y = max(pair_curve$pairs),
             label = paste0(n_low_pairs, " pairs < -5"),
             hjust = -0.1, vjust = -0.5, color = "blue")
}

if (n_high_pairs > 0) {
  p_curve <- p_curve +
    annotate("text", x = 5, y = max(pair_curve$pairs),
             label = paste0(n_high_pairs, " pairs > 5"),
             hjust = 1.1, vjust = -0.5, color = "red")
}

ggsave(file.path(out_dir, paste0("d_distribution_pairs_curve_", unique_id, ".png")),
       plot = p_curve, width = 8, height = 6, dpi = 300, bg = "white")

# Both importance AND D-value cutoff filter
genes_both_filtered <- genes_imp_filtered[
  !is.na(dval_map[genes_imp_filtered]) &
  dval_map[genes_imp_filtered] >= dcutoff_value
]

# -------------------------
# Histogram: gene-level D-values
# -------------------------

# Gene-level data for importance cutoff only
df_genes_imp <- data.frame(
  Gene   = genes_imp_filtered,
  D_value = dval_map[genes_imp_filtered],
  stringsAsFactors = FALSE
)

# Calculate outliers above and below D-value of 5 and -5
n_low  <- sum(df_genes_imp$D_value < -5, na.rm = TRUE)
n_high <- sum(df_genes_imp$D_value >  5, na.rm = TRUE)

p_hist_genes <- ggplot(df_genes_imp, aes(x = D_value)) +
  geom_histogram(
    binwidth = 0.25, boundary = dcutoff_value,
    fill = "#6baed6", color = "#08519c", alpha = 0.6
  ) +
  geom_vline(data = cutoff_df,
             aes(xintercept = cutoff, color = Stat, linetype = Stat),
             linewidth = 1) +
  scale_color_manual(values = c("Elbow" = "#e31a1c", "Coinfinder" = "#450808")
  ) +
  scale_linetype_manual(values = c("Elbow" = "dashed", "Coinfinder" = "dotted")
  ) +
  scale_x_continuous(
    breaks = seq(-5, 5, by = 1.0)   
  ) +
  labs(
    title = paste("Distribution of gene D-values -", unique_id),
    subtitle = subtitle_text_genes,
    x = "D-value", y = "Gene count",
    color    = "Cutoff type", linetype = "Cutoff type"
  ) +
  coord_cartesian(xlim = c(-5, 5)) +
  theme_minimal()

# If outliers, annotate in plot
if (n_low > 0) {
  p_hist_genes <- p_hist_genes +
    annotate(
      "text",
      x = -Inf, y = Inf,
      label = paste0(n_low, " genes < D: -5"),
      hjust = -0.1, vjust = 2,
      color = "blue"
    )
}

if (n_high > 0) {
  p_hist_genes <- p_hist_genes +
    annotate(
      "text",
      x = -Inf, y = Inf,
      label = paste0(n_high, " genes > D: 5"),
      hjust = -0.1, vjust = 3.5,
      color = "red"
    )
}

# -------------------------
# Histogram: pair-level D-values
# -------------------------

# Pair-level data for importance cutoff only
df_pairs_imp <- subset(
  df_pairs,
  !is.na(D_value) & D_value >= imp_cutoff
)

p_hist_pairs <- ggplot(df_pairs_imp, aes(x = D_value)) +
  geom_histogram(
    binwidth = 0.25, boundary = dcutoff_value,
    fill = "#fdae6b", color = "#e6550d", alpha = 0.6
  ) +
  geom_vline(data = cutoff_df,
             aes(xintercept = cutoff, color = Stat, linetype = Stat),
             linewidth = 1) +
  scale_color_manual(values = c("Elbow" = "#e31a1c", "Coinfinder" = "#450808")
  ) +
  scale_linetype_manual(values = c("Elbow" = "dashed", "Coinfinder" = "dotted")
  ) +
  scale_x_continuous(
    breaks = seq(-5, 5, by = 1.0)   
  ) +
  labs(
    title = paste("Distribution of gene pairs D-values -", unique_id),
    subtitle = subtitle_text_pairs,
    x = "D-value (min of source & target)",
    y = "Pair count",
    color    = "Cutoff type", linetype = "Cutoff type"
  ) +
  coord_cartesian(xlim = c(-5, 5)) +
  theme_minimal()
  
# If outliers, annotate in plot
if (n_low_pairs > 0) {
  p_hist_pairs <- p_hist_pairs +
    annotate(
      "text",
      x = -Inf, y = Inf,
      label = paste0(n_low_pairs, " pairs < D: -5"),
      hjust = -0.1, vjust = 2,
      color = "blue"
    )
}

if (n_high_pairs > 0) {
  p_hist_pairs <- p_hist_pairs +
    annotate(
      "text",
      x = -Inf, y = Inf,
      label = paste0(n_high_pairs, " pairs > D: 5"),
      hjust = -0.1, vjust = 3.5,
      color = "red"
    )
}

# Combine into two-panel figure
combined_dval_plots <- (p_hist_genes + p_hist_pairs) +
  plot_layout(guides = "collect") &
  theme(legend.position = "right")
output_file_dvals <- file.path(out_dir, "dvalue_histograms.png")
ggsave(output_file_dvals, plot = combined_dval_plots, width = 12, height = 6, dpi = 300)

# Per-gene summaries
gene_counts <- rowSums(strong_mat, na.rm = TRUE)
gene_mean_scores <- apply(imp_sym, 1, function(x) {
  vals <- x[!is.na(x) & x >= imp_cutoff]
  if (length(vals) == 0) return(NA)
  mean(vals)
})

genes_both_filtered <- genes_imp_filtered[!is.na(dval_map[genes_imp_filtered]) & dval_map[genes_imp_filtered] >= dcutoff_value]

# Summaries
df_counts <- data.frame(Gene = gene_names, value = gene_counts, metric = "Partner count", stringsAsFactors = FALSE) %>% filter(Gene %in% genes_both_filtered)
df_means  <- data.frame(Gene = gene_names, value = gene_mean_scores, metric = "Mean importance score", stringsAsFactors = FALSE) %>% filter(Gene %in% genes_both_filtered)

# Partner and mean importance score distribution histograms
# Distribution of partner counts (all genes)
p_counts <- ggplot(df_counts, aes(x = value)) +
  geom_histogram(binwidth = 1, fill = "#6baed6", color = "#08519c", alpha = 0.6) +
  labs(
    title = "Distribution of gene partners (D-value filtered)",
    subtitle = paste0(
      "Importance cutoff = ", round(imp_cutoff, 3),
      " | Mean gene partners = ", round(mean(df_counts$value, na.rm = TRUE), 2)
      ),
    x = "Number of partners >= cutoff",
    y = "Gene count"
  ) +
  theme_minimal()

# Distribution of mean importance scores (only genes above D-value cutoff)
p_means <- ggplot(df_means, aes(x = value)) +
  geom_histogram(fill = "#fc9272", color = "#cb181d", alpha = 0.6, bins = 30) +
  labs(
    title = "Distribution of mean importance scores (D-value filtered)",
    subtitle = paste0(
      "D-value cutoff = ", round(dcutoff_value, 3),
      " | Average mean importance score = ", round(mean(df_means$value, na.rm = TRUE), 3)
      ),
    x = "Mean importance score >= cutoff",
    y = "Gene count"
  ) +
  theme_minimal()

# Combine into two-panel figure
combined_gene_plots <- p_counts + p_means
ggsave(file.path(out_dir, "per_gene_histograms.png"),
       plot = combined_gene_plots, width = 12, height = 6, dpi = 300)

# Write output data
df_out <- bind_rows(
  df_genes %>% mutate(Gene = Gene),
  df_pairs %>% mutate(Gene = NA_character_)
)
write_csv(df_out, file.path(out_dir, paste0("panforest_dvalues_", unique_id, ".csv")))
message("D-value CSV saved to: ", file.path(out_dir, paste0("panforest_dvalues_", unique_id, ".csv")))
