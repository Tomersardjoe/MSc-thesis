#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ineq)
  library(ggplot2)
  library(ggrepel)
})

select_top_generalists <- function(df, n = 5) {
  df$rank_total <- rank(-df$total_importance, ties.method = "min")   # higher is better
  df$rank_mean  <- rank(-df$mean_importance, ties.method = "min")    # higher is better
  df$rank_gini  <- rank(df$gini, ties.method = "min")                # lower is better
  df$rank_dval  <- rank(df$D_value, ties.method = "min")             # lower is better
  
  df$avg_rank <- (df$rank_total + df$rank_mean + df$rank_gini + df$rank_dval) / 4
  df[order(df$avg_rank), ][1:n, ]
}

select_top_specialists <- function(df, n = 5) {
  df$rank_total <- rank(-df$total_importance, ties.method = "min")   # higher is better
  df$rank_mean  <- rank(df$mean_importance, ties.method = "min")     # lower is better
  df$rank_gini  <- rank(-df$gini, ties.method = "min")               # higher is better
  df$rank_dval  <- rank(df$D_value, ties.method = "min")             # lower is better
  
  df$avg_rank <- (df$rank_total + df$rank_mean + df$rank_gini + df$rank_dval) / 4
  df[order(df$avg_rank), ][1:n, ]
}

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {stop("Usage: Rscript inspect_imp.R /path/to/panforest/{run_id}/imp_fixed.csv /path/to/panforest/{run_id}_nodes.tsv /path/to/panforest/{run_id}_imp_cutoff.txt /path/to/coinfinder/{run_id}_d_cutoff.txt")} # NOT THE COINFINDER nodes_all.tsv but the output of PanForest's calculate_d.R

imp_path       <- args[1]
dval_path      <- args[2]
cutoff_path    <- args[3]
d_cutoff_path  <- args[4]
imp_dir        <- dirname(imp_path)
unique_id      <- basename(dirname(imp_path))

imp_mat        <- read.csv(imp_path, row.names = 1, check.names = FALSE)
imp_mat        <- as.matrix(imp_mat)
mode(imp_mat)  <- "numeric"

row_sum  <- rowSums(imp_mat, na.rm = TRUE)   # total importance contributed by each gene
row_mean <- rowMeans(imp_mat, na.rm = TRUE)  # average importance per target

# Read in importance score cutoff value
imp_cutoff <- as.numeric(read.table(cutoff_path)[[1]])

# Create a thresholded copy: values below cutoff set to NA
imp_mat_cut <- imp_mat
imp_mat_cut[imp_mat_cut < imp_cutoff] <- NA

keep_genes <- rownames(imp_mat_cut)[rowSums(!is.na(imp_mat_cut)) > 0]

# Make sure they exist in both rows and columns
keep_genes <- intersect(keep_genes, rownames(imp_mat_cut))
keep_genes <- intersect(keep_genes, colnames(imp_mat_cut))

if (length(keep_genes) == 0) {
  warning("No genes survived cutoff")
  imp_mat_cut <- imp_mat_cut[0, 0, drop = FALSE]
} else {
  imp_mat_cut <- imp_mat_cut[keep_genes, keep_genes, drop = FALSE]
}

row_sum_cut  <- rowSums(imp_mat_cut, na.rm = TRUE)
row_mean_cut <- rowMeans(imp_mat_cut, na.rm = TRUE)

# Compute Gini for each row (~1, one or few targets, ~0, many targets)
row_gini <- apply(imp_mat, 1, function(x) {
  if (all(x == 0 | is.na(x))) return(NA_real_)
  ineq::ineq(x, type = "Gini")
})

row_gini_cut <- apply(imp_mat_cut, 1, function(x) {
  if (all(is.na(x))) return(NA_real_)
  ineq::ineq(x, type = "Gini", na.rm = TRUE)
})

# Read D-values (genes)
dvals <- read.delim(dval_path, header = TRUE, stringsAsFactors = FALSE)
if (!all(c("ID", "Result") %in% colnames(dvals))) {
  stop("D-value file must have columns: ID, Result")
}

# Adjust column headers
colnames(dvals)[colnames(dvals) == "ID"] <- "gene"
colnames(dvals)[colnames(dvals) == "Result"] <- "D_value"

# Read D-value cutoff
d_cutoff <- as.numeric(read.table(d_cutoff_path)[[1]])

row_summary <- data.frame(
  gene            = rownames(imp_mat),
  total_importance = row_sum,
  mean_importance  = row_mean,
  gini             = row_gini,
  stringsAsFactors = FALSE
)

row_summary_cut <- data.frame(
  gene             = rownames(imp_mat_cut),
  total_importance = row_sum_cut,
  mean_importance  = row_mean_cut,
  gini             = row_gini_cut,
  stringsAsFactors = FALSE
)

# Merge with D-values
row_summary     <- merge(row_summary, dvals, by = "gene", all.x = TRUE)
row_summary_cut <- merge(row_summary_cut, dvals, by = "gene", all.x = TRUE)
n_before_d <- nrow(row_summary)

# Keep only genes with D-value above cutoff
row_summary     <- subset(row_summary,     D_value > d_cutoff)
row_summary_cut <- subset(row_summary_cut, D_value > d_cutoff)
n_after_d <- nrow(row_summary) 
pct_after  <- if (n_before_d > 0) 100 * n_after_d / n_before_d else NA_real_

# Sort by total importance and then Gini
row_summary     <- row_summary[order(-row_summary$total_importance,
                                 -row_summary$gini), ]

row_summary_cut <- row_summary_cut[order(-row_summary_cut$total_importance,
                                         -row_summary_cut$gini), ]

# Save to CSV
write.csv(
  row_summary,
  file = file.path(imp_dir, paste0(unique_id, "_row_importance_gini.csv")),
  row.names = FALSE
)

# Filter for plotting
top_generalists      <- select_top_generalists(row_summary, n = 5)
top_generalists_cut  <- select_top_generalists(row_summary_cut, n = 5)
top_specialists      <- select_top_specialists(row_summary, n = 5)
top_specialists_cut  <- select_top_specialists(row_summary_cut, n = 5)

# Count how many genes in total
n_genes <- nrow(row_summary)

# Determine outline categories present in the data
outline_levels <- c()
if (nrow(top_generalists) > 0) outline_levels <- c(outline_levels, "Broad contributor")
if (nrow(top_specialists) > 0) outline_levels <- c(outline_levels, "Focused contributor")

# Identify the gene with the highest mean importance (unfiltered)
outlier_unfiltered <- row_summary[which.max(row_summary$mean_importance), , drop = FALSE]

p_unfiltered <- ggplot(row_summary, aes(x = mean_importance,
                                        y = gini,
                                        size = total_importance,
                                        color = D_value)) +
  # Background cloud
  geom_point(alpha = 0.6, shape = 16) +
  
  # Scales with explicit legend order
  scale_size_continuous(
    range = c(0.5, 4),
    guide = guide_legend(order = 1,
                         override.aes = list(shape = 16, colour = "black", fill = "black", alpha = 0.6))
  ) +
  scale_color_viridis_c(
    option = "magma", direction = -1,
    guide = guide_colorbar(order = 3)
  ) +
  scale_shape_manual(
    values = c("Broad contributor" = 21,
               "Focused contributor" = 21,
               "Highest mean importance" = 21)
  ) +
  guides(
    shape = guide_legend(
      order = 1,
      override.aes = list(
        shape  = c(21, 21, 21),
        fill   = c(NA, NA, NA),
        colour = c("black", "blue", "red"),  # black = broad, blue = focused, red = higesht mean importance
        size   = 4,
        alpha  = 1,
        stroke = 0.6
      )
    )
  ) +
  
  labs(
    title = paste("Mean importance vs Gini with D-values for", unique_id),
    subtitle = sprintf("D-cutoff = %.3f | %.1f%% survived (n = %d of %d)",
                       d_cutoff, pct_after, n_after_d, n_before_d),
    x = "Mean importance per target",
    y = "Gini coefficient",
    size = "Total importance",
    shape = "Outline"
  ) +
  theme_bw(base_size = 14) +

  # Generalists
  geom_point(data = top_generalists,
             aes(x = mean_importance, y = gini, size = total_importance, color = D_value),
             inherit.aes = FALSE, alpha = 0.9, shape = 16) +
  geom_point(data = top_generalists,
             aes(x = mean_importance, y = gini, size = total_importance),
             inherit.aes = FALSE, shape = 21, fill = NA, colour = "black", stroke = 0.6) +
  geom_text_repel(data = top_generalists,
                  aes(x = mean_importance, y = gini, label = gene),
                  inherit.aes = FALSE, size = 3, color = "black",
                  bg.color = "white", bg.r = 0.15,
                  box.padding = 0.5, point.padding = 0.3,
                  force = 1.5, segment.size = 0.2, max.overlaps = Inf) +

  # Specialists
  geom_point(data = top_specialists,
             aes(x = mean_importance, y = gini, size = total_importance, color = D_value),
             inherit.aes = FALSE, alpha = 0.9, shape = 16) +
  geom_point(data = top_specialists,
             aes(x = mean_importance, y = gini, size = total_importance),
             inherit.aes = FALSE, shape = 21, fill = NA, colour = "blue", stroke = 0.6) +
  geom_text_repel(data = top_specialists,
                  aes(x = mean_importance, y = gini, label = gene),
                  inherit.aes = FALSE, size = 3, color = "blue",
                  bg.color = "white", bg.r = 0.15,
                  box.padding = 0.5, point.padding = 0.3,
                  force = 1.5, segment.size = 0.2, max.overlaps = Inf) +

  # Outlier with max mean importance (unfiltered)
  geom_point(data = outlier_unfiltered,
             aes(x = mean_importance, y = gini, size = total_importance, color = D_value),
             inherit.aes = FALSE, shape = 16, alpha = 1) +
  geom_point(data = outlier_unfiltered,
             aes(x = mean_importance, y = gini, size = total_importance),
             inherit.aes = FALSE, shape = 21, fill = NA, colour = "red", stroke = 0.8) +
  geom_text_repel(data = outlier_unfiltered,
                  aes(x = mean_importance, y = gini, label = gene),
                  inherit.aes = FALSE, size = 3, color = "black",
                  bg.color = "white", bg.r = 0.15,
                  box.padding = 0.6, point.padding = 0.3,
                  force = 2, segment.size = 0.3, max.overlaps = Inf) +

  # Dummy legend points (invisible but trigger legend)
  geom_point(
    data = data.frame(
      highlight = c("Broad contributor", "Focused contributor", "Highest mean importance"),
      mean_importance = median(row_summary$mean_importance, na.rm = TRUE),
      gini = median(row_summary$gini, na.rm = TRUE)
    ),
    aes(x = mean_importance, y = gini, shape = highlight),
    inherit.aes = FALSE,
    alpha = 0,
    show.legend = TRUE
  )

ggsave(file.path(imp_dir, paste0(unique_id, "_mean_vs_gini.png")),
       plot = p_unfiltered, width = 8, height = 6, dpi = 300, bg = "white")

# Count how many genes above cutoff
n_survivors <- nrow(row_summary_cut)

# Determine outline categories present in the data
outline_levels <- c()
if (nrow(top_generalists_cut) > 0) outline_levels <- c(outline_levels, "Broad contributor")
if (nrow(top_specialists_cut) > 0) outline_levels <- c(outline_levels, "Focused contributor")

# Identify the gene with the highest mean importance
outlier_cut <- row_summary_cut[which.max(row_summary_cut$mean_importance), , drop = FALSE]

p_filtered <- ggplot(row_summary_cut, aes(x = mean_importance,
                                          y = gini,
                                          size = total_importance,
                                          color = D_value)) +
  # Background cloud
  geom_point(alpha = 0.6, shape = 16) +
  
  # Scales with explicit legend order
  scale_size_continuous(
    range = c(0.5, 4),
    guide = guide_legend(order = 1,
                         override.aes = list(shape = 16, colour = "black", fill = "black", alpha = 0.6))
  ) +
  scale_color_viridis_c(
    option = "magma", direction = -1,
    guide = guide_colorbar(order = 3)
  ) +
  scale_shape_manual(
    values = c("Broad contributor" = 21,
               "Focused contributor" = 21,
               "Highest mean importance" = 21)
  ) +
  guides(
    shape = guide_legend(
      order = 1,
      override.aes = list(
        shape  = c(21, 21, 21),
        fill   = c(NA, NA, NA),
        colour = c("black", "blue", "red"),  # black = broad, blue = focused, red = highest mean importance
        size   = 4,
        alpha  = 1,
        stroke = 0.6
      )
    )
  ) +
  
  labs(
    title = paste("Mean importance vs Gini with D-values for", unique_id),
    subtitle = sprintf("Importance cutoff = %.3f | %.1f%% survived (n = %d of %d)",
                       imp_cutoff, 100 * n_survivors / n_genes, n_survivors, n_genes),
    x = "Mean importance per target",
    y = "Gini coefficient",
    size = "Total importance",
    shape = "Outline" 
  ) +
  theme_bw(base_size = 14) +

  # Generalists
  geom_point(data = top_generalists_cut,
             aes(x = mean_importance, y = gini, size = total_importance, color = D_value),
             inherit.aes = FALSE, alpha = 0.9, shape = 16) +
  geom_point(data = top_generalists_cut,
             aes(x = mean_importance, y = gini, size = total_importance),
             inherit.aes = FALSE, shape = 21, fill = NA, colour = "black", stroke = 0.6) +
  geom_text_repel(data = top_generalists_cut,
                  aes(x = mean_importance, y = gini, label = gene),
                  inherit.aes = FALSE, size = 3, color = "black",
                  bg.color = "white", bg.r = 0.15,
                  box.padding = 0.5, point.padding = 0.3,
                  force = 1.5, segment.size = 0.2, max.overlaps = Inf) +

  # Specialists
  geom_point(data = top_specialists_cut,
             aes(x = mean_importance, y = gini, size = total_importance, color = D_value),
             inherit.aes = FALSE, alpha = 0.9, shape = 16) +
  geom_point(data = top_specialists_cut,
             aes(x = mean_importance, y = gini, size = total_importance),
             inherit.aes = FALSE, shape = 21, fill = NA, colour = "blue", stroke = 0.6) +
  geom_text_repel(data = top_specialists_cut,
                  aes(x = mean_importance, y = gini, label = gene),
                  inherit.aes = FALSE, size = 3, color = "blue",
                  bg.color = "white", bg.r = 0.15,
                  box.padding = 0.5, point.padding = 0.3,
                  force = 1.5, segment.size = 0.2, max.overlaps = Inf) +

  # Outlier with max mean importance
  geom_point(data = outlier_cut,
             aes(x = mean_importance, y = gini, size = total_importance, color = D_value),
             inherit.aes = FALSE, shape = 16, alpha = 1) +
  geom_point(data = outlier_cut,
             aes(x = mean_importance, y = gini, size = total_importance),
             inherit.aes = FALSE, shape = 21, fill = NA, colour = "red", stroke = 0.8) +
  geom_text_repel(data = outlier_cut,
                  aes(x = mean_importance, y = gini, label = gene),
                  inherit.aes = FALSE, size = 3, color = "black",
                  bg.color = "white", bg.r = 0.15,
                  box.padding = 0.6, point.padding = 0.3,
                  force = 2, segment.size = 0.3, max.overlaps = Inf) +

  # Dummy legend points
  geom_point(
    data = data.frame(
      highlight = c("Broad contributor", "Focused contributor", "Highest mean importance"),
      mean_importance = median(row_summary_cut$mean_importance, na.rm = TRUE),
      gini = median(row_summary_cut$gini, na.rm = TRUE)
    ),
    aes(x = mean_importance, y = gini, shape = highlight),
    inherit.aes = FALSE,
    alpha = 0,
    show.legend = TRUE
  )

ggsave(file.path(imp_dir, paste0(unique_id, "_mean_vs_gini_filtered.png")),
       plot = p_filtered, width = 8, height = 6, dpi = 300, bg = "white")

sorted_totals <- sort(row_sum, decreasing = TRUE)
top_ratio <- sorted_totals[1] / sorted_totals[2]

cat("\nTop-to-second importance ratio:", top_ratio, "\n")
