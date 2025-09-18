#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ineq)
  library(dplyr)
  library(ggplot2)
  library(ggrepel)
})

select_top_generalists <- function(df, imp_mat, n = 5) {
  # Calculate proportion of non-zero targets for each gene
  prop_nonzero <- rowMeans(imp_mat != 0, na.rm = TRUE)
  df$prop_nonzero <- prop_nonzero[match(df$gene, rownames(imp_mat))]
  
  # Rankings
  df$rank_total    <- rank(-df$total_importance, ties.method = "min")   # higher is better
  df$rank_mean     <- rank(-df$mean_importance, ties.method = "min")    # higher is better
  df$rank_gini     <- rank(df$gini, ties.method = "min")                # lower is better
  df$rank_dval     <- rank(df$D_value, ties.method = "min")             # lower is better
  df$rank_nonzero  <- rank(-df$prop_nonzero, ties.method = "min")       # higher is better
  
  # Average rank including non-zero proportion
  df$avg_rank <- (df$rank_total + df$rank_mean + df$rank_gini + df$rank_dval + df$rank_nonzero) / 5
  
  df[order(df$avg_rank), ][1:n, ]
}

select_top_specialists <- function(df, imp_mat, n = 5) {
  # Calculate proportion of non-zero targets for each gene
  prop_nonzero <- rowMeans(imp_mat != 0, na.rm = TRUE)
  df$prop_nonzero <- prop_nonzero[match(df$gene, rownames(imp_mat))]
  
  # Rankings
  df$rank_total    <- rank(-df$total_importance, ties.method = "min")   # higher is better
  df$rank_mean     <- rank(df$mean_importance, ties.method = "min")     # lower is better
  df$rank_gini     <- rank(-df$gini, ties.method = "min")               # higher is better
  df$rank_dval     <- rank(df$D_value, ties.method = "min")             # lower is better
  df$rank_nonzero  <- rank(df$prop_nonzero, ties.method = "min")        # lower is better
  
  # Average rank including non-zero proportion
  df$avg_rank <- (df$rank_total + df$rank_mean + df$rank_gini + df$rank_dval + df$rank_nonzero) / 5
  
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
top_generalists      <- select_top_generalists(row_summary,      imp_mat, n = 5)
top_generalists_cut  <- select_top_generalists(row_summary_cut,  imp_mat, n = 5)
top_specialists      <- select_top_specialists(row_summary,      imp_mat, n = 5)
top_specialists_cut  <- select_top_specialists(row_summary_cut,  imp_mat, n = 5)

# Count how many genes in total
n_genes <- nrow(row_summary)

# Determine outline categories present in the data
outline_levels <- c()
if (nrow(top_generalists) > 0) outline_levels <- c(outline_levels, "Broad contributor")
if (nrow(top_specialists) > 0) outline_levels <- c(outline_levels, "Focused contributor")

# Identify the gene with the highest mean importance (unfiltered)
outlier_unfiltered <- row_summary[which.max(row_summary$mean_importance), , drop = FALSE]

# Combine all labelled points into one data frame
label_df <- bind_rows(
  top_generalists  %>% mutate(label_color = "black", outline_color = "black"),
  top_specialists  %>% mutate(label_color = "blue",  outline_color = "blue"),
  outlier_unfiltered %>% mutate(label_color = "black", outline_color = "red")
)

p_unfiltered <- ggplot(row_summary, aes(x = mean_importance,
                                        y = gini,
                                        size = total_importance,
                                        color = D_value)) +
  geom_point(alpha = 0.6, shape = 16) +
  scale_size_continuous(range = c(0.5, 4),
                        guide = guide_legend(order = 1,
                                             override.aes = list(shape = 16, colour = "black", fill = "black", alpha = 0.6))) +
  scale_color_viridis_c(option = "magma", direction = -1,
                        guide = guide_colorbar(order = 3)) +
  scale_shape_manual(values = c("Broad contributor" = 21,
                                "Focused contributor" = 21,
                                "Highest mean importance" = 21)) +
  guides(shape = guide_legend(order = 1,
                              override.aes = list(shape  = c(21, 21, 21),
                                                  fill   = c(NA, NA, NA),
                                                  colour = c("black", "blue", "red"),
                                                  size   = 4,
                                                  alpha  = 1,
                                                  stroke = 0.6))) +
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

  # Points for each category
  geom_point(data = top_generalists,
             aes(size = total_importance, color = D_value),
             alpha = 0.9, shape = 16) +
  geom_point(data = top_generalists,
             aes(size = total_importance),
             shape = 21, fill = NA, colour = "black", stroke = 0.6) +

  geom_point(data = top_specialists,
             aes(size = total_importance, color = D_value),
             alpha = 0.9, shape = 16) +
  geom_point(data = top_specialists,
             aes(size = total_importance),
             shape = 21, fill = NA, colour = "blue", stroke = 0.6) +

  geom_point(data = outlier_unfiltered,
             aes(size = total_importance, color = D_value),
             shape = 16, alpha = 1) +
  geom_point(data = outlier_unfiltered,
             aes(size = total_importance),
             shape = 21, fill = NA, colour = "red", stroke = 0.6) +

  # Single geom_text_repel for all labels
  ggrepel::geom_text_repel(
    data = label_df,
    aes(label = gene, color = NULL),  # don't map color to D_value here
    size = 3,
    color = label_df$label_color,
    bg.color = "white", bg.r = 0.15,
    box.padding = 0.4, point.padding = 0.7,
    min.segment.length = 0.01,
    force = 1.7, segment.size = 0.2, max.overlaps = Inf
  ) +

  # Dummy legend points
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

# Save
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

# Combine all labelled points into one data frame
label_df_cut <- bind_rows(
  top_generalists_cut  %>% mutate(label_color = "black", outline_color = "black"),
  top_specialists_cut  %>% mutate(label_color = "blue",  outline_color = "blue"),
  outlier_cut %>% mutate(label_color = "black", outline_color = "red")
)

p_filtered <- ggplot(row_summary_cut, aes(x = mean_importance,
                                          y = gini,
                                          size = total_importance,
                                          color = D_value)) +
  geom_point(alpha = 0.6, shape = 16) +
  scale_size_continuous(
    range = c(0.5, 4),
    guide = guide_legend(order = 1,
                         override.aes = list(shape = 16, colour = "black", fill = "black", alpha = 0.6))
  ) +
  scale_color_viridis_c(option = "magma", direction = -1,
                        guide = guide_colorbar(order = 3)) +
  scale_shape_manual(values = c("Broad contributor" = 21,
                                "Focused contributor" = 21,
                                "Highest mean importance" = 21)) +
  guides(shape = guide_legend(order = 1,
                              override.aes = list(
                                shape  = c(21, 21, 21),
                                fill   = c(NA, NA, NA),
                                colour = c("black", "blue", "red"),
                                size   = 4,
                                alpha  = 1,
                                stroke = 0.6
                              ))) +
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

  # Points for generalists
  geom_point(data = top_generalists_cut,
             aes(size = total_importance, color = D_value),
             alpha = 0.9, shape = 16) +
  geom_point(data = top_generalists_cut,
             aes(size = total_importance),
             shape = 21, fill = NA, colour = "black", stroke = 0.6) +

  # Points for specialists
  geom_point(data = top_specialists_cut,
             aes(size = total_importance, color = D_value),
             alpha = 0.9, shape = 16) +
  geom_point(data = top_specialists_cut,
             aes(size = total_importance),
             shape = 21, fill = NA, colour = "blue", stroke = 0.6) +

  # Points for outlier
  geom_point(data = outlier_cut,
             aes(size = total_importance, color = D_value),
             shape = 16, alpha = 1) +
  geom_point(data = outlier_cut,
             aes(size = total_importance),
             shape = 21, fill = NA, colour = "red", stroke = 0.6) +

  # Single geom_text_repel for all labels
  geom_text_repel(
    data = label_df_cut,
    aes(label = gene),
    size = 3,
    color = label_df_cut$label_color,
    bg.color = "white", bg.r = 0.15,
    box.padding = 0.4, point.padding = 0.7,
    min.segment.length = 0.01,
    force = 1., segment.size = 0.2, max.overlaps = Inf
  ) +

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

# Pick one generalist and one specialist to compare
gen_gene <- top_generalists_cut$gene[1]
spec_gene <- top_specialists_cut$gene[1]

# Extract values from unthresholded matrix
gen_vals <- as.numeric(imp_mat[gen_gene, ])
spec_vals <- as.numeric(imp_mat[spec_gene, ])

# Keep zeros (important for sparsity story)
gen_vals <- gen_vals[!is.na(gen_vals)]
spec_vals <- spec_vals[!is.na(spec_vals)]

# Combine into one data frame
plot_df <- data.frame(
  Importance = c(gen_vals, spec_vals),
  Type = c(rep("Broad contributor", length(gen_vals)),
           rep("Focused contributor", length(spec_vals)))
)

zero_counts <- plot_df %>%
  group_by(Type) %>% # There is only one gene per type
  summarise(
    zeros = sum(Importance == 0),
    total = n(),
    perc_zeros = round(100 * zeros / total, 1)
  )

# Add a label column for display
zero_counts <- zero_counts %>%
  mutate(label = paste0("0-scores: ", zeros, " (", perc_zeros, "%)"))

# Log-scaled density plot
p_log <- ggplot(plot_df, aes(x = Importance, fill = Type)) +
  geom_density(alpha = 0.5, adjust = 1) +
  facet_wrap(~ Type, scales = "fixed") +
  scale_y_continuous(trans = "log1p") +  # log(1 + y)
  scale_fill_manual(values = c("Broad contributor" = "#1b9e77",
                               "Focused contributor" = "#d95f02")) +
  labs(
    title = paste("Log-transformed importance distribution for", unique_id),
    subtitle = paste(gen_gene, "(broad) vs", spec_gene, "(focused)"),
    x = "Importance score",
    y = expression(log[1](1 + density))
  ) +
  theme_minimal(base_size = 14) +
  # Add zero count labels in top-right of each facet
  geom_text(
    data = zero_counts,
    aes(
      x = max(plot_df$Importance) * 0.95,
      y = Inf,
      label = label
    ),
    inherit.aes = FALSE,
    hjust = 1, vjust = 1.5
  )

ggsave(file.path(imp_dir, paste0(unique_id, "_imp_density.png")),
       plot = p_log, width = 8, height = 6, dpi = 300, bg = "white")
