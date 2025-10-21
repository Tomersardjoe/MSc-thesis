#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(ggrepel)
  library(dplyr)
  library(patchwork)
  library(readr)
})

# -------------------------
# Helpers
# -------------------------

annotate_outliers <- function(plot, n_low, n_high, label_prefix = "genes") {
  if (n_low > 0) {
    plot <- plot + annotate("text", x = -Inf, y = Inf,
                            label = paste0(n_low, " ", label_prefix, " < D: -5"),
                            hjust = -0.1, vjust = 2, color = "blue")
  }
  if (n_high > 0) {
    plot <- plot + annotate("text", x = -Inf, y = Inf,
                            label = paste0(n_high, " ", label_prefix, " > D: 5"),
                            hjust = -0.1, vjust = 3.5, color = "red")
  }
  plot
}

make_subtitle <- function(stage, entity, counts, totals, cutoffs) {
  retained <- counts[[stage]][[entity]]
  total    <- totals[[stage]][[entity]]
  pct      <- if (!is.na(total) && total > 0) round(100 * retained / total, 1) else NA
  cutoff   <- cutoffs[[stage]]

  cutoff_val <- cutoff$value
  if (is.numeric(cutoff_val) && length(cutoff_val) == 1) {
    cutoff_val <- round(cutoff_val, 2)
  } else if (is.numeric(cutoff_val) && length(cutoff_val) == 2) {
    cutoff_val <- paste(round(cutoff_val, 2), collapse = "/")
  }

  paste0(
    cutoff$label, " cutoff = ", cutoff_val,
    " | ", toupper(substr(entity,1,1)), substr(entity,2,nchar(entity)),
    " retained = ", retained, "/", total,
    if (!is.na(pct)) paste0(" (", pct, "%)") else ""
  )
}

make_hist <- function(df, subtitle, title, label_prefix, cutoff_df, dcutoff_value) {
  n_low  <- sum(df$D_value < -5, na.rm = TRUE)
  n_high <- sum(df$D_value >  5, na.rm = TRUE)

  fill_col  <- ifelse(label_prefix == "genes", "#6baed6", "#fdae6b")
  border_col <- ifelse(label_prefix == "genes", "#08519c", "#e6550d")

  p <- ggplot(df, aes(x = D_value)) +
    geom_histogram(binwidth = 0.25, boundary = dcutoff_value,
                   fill = fill_col, color = border_col, alpha = 0.6) +
    geom_vline(data = cutoff_df,
               aes(xintercept = cutoff, color = Stat, linetype = Stat),
               linewidth = 1) +
    scale_color_manual(
      values = c("Elbow"="#e31a1c","Q3"="#450808"),
      labels = c(
        paste0("Elbow (", round(elbow_dcutoff, 2), ")"),
        paste0("Q3 (", round(dcutoff_value, 2), ")")
      )
    ) +
    scale_linetype_manual(
      values = c("Elbow"="dashed","Q3"="dotted"),
      labels = c(
        paste0("Elbow (", round(elbow_dcutoff, 2), ")"),
        paste0("Q3 (", round(dcutoff_value, 2), ")")
      )
    ) +
    scale_x_continuous(breaks = seq(-5,5,1)) +
    labs(title = title, subtitle = subtitle,
         x = "D-value", y = paste0(toupper(substr(label_prefix,1,1)),
                                   substr(label_prefix,2,nchar(label_prefix)),
                                   " count"),
         color = "Cutoff type", linetype = "Cutoff type") +
    coord_cartesian(xlim = c(-5,5)) +
    theme_minimal()

  annotate_outliers(p, n_low, n_high, label_prefix)
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

# -------------------------
# Load data and compute cutoffs
# -------------------------

# Read input arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {stop("Usage: Rscript imp_distribution.R
    /path/to/panforest/{run_id}/imp_.csv 
    /path/to/panforest/{run_id}/performance.csv
    /path/to/panforest/imp_cutoff/{run_id}_nodes.tsv # DO NOT USE THE COINFINDER nodes_all.tsv but the output of PanForest's calculate_d.R
    [D-value cutoff]")}

imp_path      <- args[1]
perf_path     <- args[2]
dval_path     <- args[3]
dcutoff_value <- suppressWarnings(as.numeric(args[4]))

out_dir   <- dirname(dval_path)
unique_id <- basename(dirname(out_dir))

imp_cutoff         <- 0.01
accuracy_cutoff    <- 0.9
f1_cutoff          <- 0.9

# Read importance matrix
imp         <- read.csv(imp_path, row.names = 1, check.names = FALSE)
imp_sym_raw <- make_symmetric_mean(as.matrix(imp))
gene_names  <- rownames(imp_sym_raw)

message("Initial matrix: ", length(gene_names), " genes, ",
        sum(upper.tri(imp_sym_raw)), " gene pairs")

# Read D-values
dvals <- read.delim(dval_path, header = TRUE, stringsAsFactors = FALSE)
dval_map <- setNames(dvals$Result, dvals$ID)

# Read performance metrics
perf <- read.csv(perf_path, header = TRUE, stringsAsFactors = FALSE, na.strings = c("NA",""))
perf_map_acc <- setNames(as.numeric(perf$Ate), perf[[1]])
perf_map_f1  <- setNames(as.numeric(perf$F1te), perf[[1]])


# Raw data frames
df_genes_raw <- data.frame(Gene = names(dval_map),
                           D_value = dval_map,
                           stringsAsFactors = FALSE) %>%
  filter(!is.na(D_value))

g1 <- rownames(imp_sym_raw)[row(imp_sym_raw)[upper.tri(imp_sym_raw)]]
g2 <- rownames(imp_sym_raw)[col(imp_sym_raw)[upper.tri(imp_sym_raw)]]
d_pair_all <- pmin(dval_map[g1], dval_map[g2], na.rm = TRUE)

df_pairs_raw <- data.frame(Gene_1 = g1, Gene_2 = g2, D_value = d_pair_all,
                           stringsAsFactors = FALSE) %>%
  filter(!is.na(D_value))

# Apply D-value cutoff

keep_genes <- gene_names[!is.na(dval_map[gene_names]) & dval_map[gene_names] >= dcutoff_value]
imp_sym <- imp_sym_raw[keep_genes, keep_genes, drop = FALSE]

n_genes_after_d <- length(keep_genes)
pairs_after_d   <- sum(upper.tri(imp_sym))

message("After D-value cutoff (", dcutoff_value, "): ",
        n_genes_after_d, " genes retained and ", pairs_after_d, " gene pairs retained")

# Apply importance cutoff
before_pairs <- sum(!is.na(imp_sym[upper.tri(imp_sym)]))
imp_sym[imp_sym < imp_cutoff] <- NA
after_pairs <- sum(!is.na(imp_sym[upper.tri(imp_sym)]))

# Define cutoffs list for subtitles
cutoffs <- list(
  raw  = list(label = "D-value",     value = dcutoff_value),
  imp  = list(label = "Importance",  value = imp_cutoff),
  perf = list(label = "Accuracy/F1", value = c(accuracy_cutoff, f1_cutoff))
)

# Genes with at least one strong connection
strong_mat <- imp_sym >= imp_cutoff
genes_imp_filtered <- rownames(strong_mat)[rowSums(strong_mat, na.rm = TRUE) > 0]

# Build dataframes for genes and pairs surviving the importance score cutoff
df_genes_imp <- data.frame(Gene = rownames(strong_mat),
                           D_value = dval_map[rownames(strong_mat)],
                           stringsAsFactors = FALSE) %>%
  filter(!is.na(D_value),
         Gene %in% genes_imp_filtered)

# Pairs surviving importance cutoff
pair_list_imp <- which(upper.tri(strong_mat) & strong_mat, arr.ind = TRUE)  # TRUE where >= cutoff
df_pairs_imp <- data.frame(
  Gene_1 = rownames(strong_mat)[pair_list_imp[,1]],
  Gene_2 = rownames(strong_mat)[pair_list_imp[,2]],
  D_value = pmin(dval_map[rownames(strong_mat)[pair_list_imp[,1]]],
                 dval_map[rownames(strong_mat)[pair_list_imp[,2]]],
                 na.rm = TRUE),
  stringsAsFactors = FALSE
) %>% filter(!is.na(D_value))

message("After importance score cutoff (", imp_cutoff, "): ",
        length(genes_imp_filtered), " genes retained and ", nrow(df_pairs_imp), " pairs retained")
        
# For plotting only! Keep genes/pairs below D-value cutoff
imp_sym_all <- imp_sym_raw
imp_sym_all[imp_sym_all < imp_cutoff] <- NA
strong_mat_all <- imp_sym_all >= imp_cutoff
genes_imp_filtered_all <- rownames(strong_mat_all)[rowSums(strong_mat_all, na.rm = TRUE) > 0]

df_genes_imp_plot <- df_genes_raw %>%
  filter(Gene %in% genes_imp_filtered_all)

df_pairs_imp_plot <- df_pairs_raw %>%
  filter(Gene_1 %in% genes_imp_filtered_all & Gene_2 %in% genes_imp_filtered_all)

# Pairs-retained curve and elbow cutoff diagnostics

if (nrow(df_pairs_raw) > 0) {
  cutoff_seq <- seq(
    min(df_pairs_raw$D_value, na.rm = TRUE),
    max(df_pairs_raw$D_value, na.rm = TRUE),
    length.out = 200
  )

  pair_counts <- sapply(cutoff_seq, function(cut) {
    sum(df_pairs_raw$D_value >= cut, na.rm = TRUE)
  })

  pair_curve <- data.frame(cutoff = cutoff_seq, pairs = pair_counts)
} else {
  pair_curve <- data.frame(cutoff = numeric(0), pairs = numeric(0))
}

# Elbow cutoff on the pair curve
elbow_dcutoff <- elbow_cutoff(pair_curve, na_ok = TRUE)

# Build cutoff_df for later histograms
cutoff_df <- data.frame(Stat = "Elbow", cutoff = elbow_dcutoff)
if (!is.na(dcutoff_value)) {
  cutoff_df <- rbind(cutoff_df,
                     data.frame(Stat = "Q3", cutoff = dcutoff_value))
}

# Plot pairs-retained curve
p_curve <- ggplot(pair_curve, aes(x = cutoff, y = pairs)) +
  geom_point(size = 1.5, color = "black") +
  geom_line(color = "#08519c", linewidth = 1) +
  geom_vline(data = cutoff_df,
             aes(xintercept = cutoff, color = Stat, linetype = Stat),
             linewidth = 1) +
  scale_color_manual(
    values = c("Elbow"="#e31a1c","Q3"="#450808"),
    labels = c(
      paste0("Elbow (", round(elbow_dcutoff, 2), ")"),
      paste0("Q3 (", round(dcutoff_value, 2), ")")
    )
  ) +
  scale_linetype_manual(
    values = c("Elbow"="dashed","Q3"="dotted"),
    labels = c(
      paste0("Elbow (", round(elbow_dcutoff, 2), ")"),
      paste0("Q3 (", round(dcutoff_value, 2), ")")
    )
  ) +
  labs(title="Pairs retained vs. D-value cutoff",
       subtitle=paste0("Pairs retained after D-value cutoff: ", pairs_after_d),
       x="Lineage independence cutoff (D)", y="# of gene pairs retained",
       color="Cutoff type", linetype="Cutoff type") +
  coord_cartesian(xlim=c(-5,5)) +
  theme_minimal() +
  theme(legend.position="right")

# Annotate extreme outliers if present
n_low_pairs  <- sum(df_pairs_raw$D_value < -5, na.rm = TRUE)
n_high_pairs <- sum(df_pairs_raw$D_value >  5, na.rm = TRUE)

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

totals <- list(
  raw  = list(genes = length(gene_names),             pairs = sum(upper.tri(imp_sym_raw))),
  imp  = list(genes = length(rownames(imp_sym)),      pairs = sum(upper.tri(imp_sym))),
  perf = list(genes = length(genes_imp_filtered),     pairs = nrow(df_pairs_imp))
)


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
      label = paste0("Importance score cutoff\n", format(Value, scientific = FALSE, digits = 3)),
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

# Performance cutoff data

gene_acc <- perf_map_acc[rownames(imp_sym)]
gene_f1  <- perf_map_f1[rownames(imp_sym)]
pass_gene <- !is.na(gene_acc) & !is.na(gene_f1) &
             gene_acc >= accuracy_cutoff & gene_f1 >= f1_cutoff

n_genes_pass   <- sum(pass_gene)
pairs_remaining <- sum(!is.na(imp_sym[upper.tri(imp_sym)]))

# Final gene set after all filters
genes_both_filtered <- genes_imp_filtered[
  !is.na(dval_map[genes_imp_filtered]) &
  dval_map[genes_imp_filtered] >= dcutoff_value &
  genes_imp_filtered %in% rownames(imp_sym)[pass_gene]
]

# Build df_genes_perf / df_pairs_perf
df_genes_perf <- df_genes_imp %>% filter(Gene %in% genes_both_filtered)
df_pairs_perf <- df_pairs_imp %>% filter(Gene_1 %in% genes_both_filtered &
                                     Gene_2 %in% genes_both_filtered)
                                     
# For plotting only! Keep genes/pairs below D-value cutoff
gene_acc_all <- perf_map_acc[rownames(imp_sym_raw)]
gene_f1_all  <- perf_map_f1[rownames(imp_sym_raw)]
pass_gene_all <- !is.na(gene_acc_all) & !is.na(gene_f1_all) &
                 gene_acc_all >= accuracy_cutoff & gene_f1_all >= f1_cutoff

genes_perf_filtered_all <- rownames(imp_sym_raw)[pass_gene_all]

df_genes_perf_plot <- df_genes_raw %>%
  filter(Gene %in% genes_perf_filtered_all)

df_pairs_perf_plot <- df_pairs_raw %>%
  filter(Gene_1 %in% genes_perf_filtered_all & Gene_2 %in% genes_perf_filtered_all)
                                     
message("After accuracy/F1 cutoffs (", accuracy_cutoff, "/", f1_cutoff, "): ",
        n_genes_pass, " genes retained and ", nrow(df_pairs_perf), " pairs retained")

# -------------------------
# Organize counts and data frames
# -------------------------

counts <- list(
  raw  = list(genes = n_genes_after_d,            pairs = pairs_after_d),
  imp  = list(genes = length(genes_imp_filtered), pairs = nrow(df_pairs_imp)),
  perf = list(genes = n_genes_pass,               pairs = nrow(df_pairs_perf))
)

dfs <- list(
  raw  = list(genes = df_genes_raw,       pairs = df_pairs_raw),
  imp  = list(genes = df_genes_imp_plot,  pairs = df_pairs_imp_plot),
  perf = list(genes = df_genes_perf_plot, pairs = df_pairs_perf_plot)
)

# -------------------------
# Generate histograms for each stage
# -------------------------

for (stage in names(dfs)) {
  p_genes <- make_hist(dfs[[stage]]$genes,
                       make_subtitle(stage, "genes", counts, totals, cutoffs),
                       paste("Distribution of gene D-values -", unique_id),
                       "genes", cutoff_df, dcutoff_value)

  p_pairs <- make_hist(dfs[[stage]]$pairs,
                       make_subtitle(stage, "pairs", counts, totals, cutoffs),
                       paste("Distribution of gene pairs D-values -", unique_id),
                       "pairs", cutoff_df, dcutoff_value)

  combined <- (p_genes + p_pairs) + plot_layout(guides="collect") &
              theme(legend.position="right")

  ggsave(file.path(out_dir, paste0(stage, "_histograms.png")),
         plot=combined, width=12, height=6, dpi=300)
}

# Per-gene summaries
gene_counts <- rowSums(strong_mat, na.rm = TRUE)
gene_mean_scores <- apply(imp_sym, 1, function(x) {
  vals <- x[!is.na(x) & x >= imp_cutoff]
  if (length(vals) == 0) return(NA)
  mean(vals)
})

# Summaries
df_counts <- data.frame(
  Gene   = names(gene_counts),
  value  = gene_counts,
  metric = "Partner count",
  stringsAsFactors = FALSE
) %>% filter(Gene %in% genes_both_filtered)

df_means <- data.frame(
  Gene   = names(gene_mean_scores),
  value  = gene_mean_scores,
  metric = "Mean importance score",
  stringsAsFactors = FALSE
) %>% filter(Gene %in% genes_both_filtered)

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
  df_genes_perf %>% mutate(Gene = Gene),
  df_pairs_perf %>% mutate(Gene = NA_character_)
)
write_csv(df_out, file.path(out_dir, paste0("panforest_dvalues_", unique_id, ".csv")))
message("D-value CSV saved to: ", file.path(out_dir, paste0("panforest_dvalues_", unique_id, ".csv")))
