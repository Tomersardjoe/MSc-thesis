#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(ggrepel)
  library(patchwork)
  library(phangorn)
  library(dplyr)
  library(ape)
  library(ggtree)
  library(readr)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {
  stop("Usage: Rscript d_distribution_coinfinder.R /path/to/coincident_nodes_all.tsv /path/to/coincident_pairs.tsv /path/to/tree.nwk /path/to/presence_absence.csv")
}

nodes_path <- args[1]
pairs_path <- args[2]
tree_path  <- args[3]
pa_path    <- args[4]
unique_id  <- basename(dirname(nodes_path))

# Create output directory
out_dir <- file.path(dirname(nodes_path), "d_cutoff")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# -------------------------
# Read nodes.tsv (gene-level D-values)
# -------------------------
nodes <- read.table(nodes_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
nodes$Result <- as.numeric(nodes$Result)

if (all(is.na(nodes$Result))) {
    message("Exiting run ", unique_id, " - no valid D-values in", nodes_path)
    quit(save = "no", status = 1)
}

# -------------------------
# Calculate summary stats
# -------------------------
cutoff_25    <- quantile(nodes$Result, 0.25, na.rm = TRUE)
median_val   <- median(nodes$Result, na.rm = TRUE)
cutoff_value <- quantile(nodes$Result, 0.75, na.rm = TRUE)

# Save d_cutoff to a text file
cutoff_file <- file.path(out_dir, paste0(unique_id, "_d_cutoff.txt"))
write(cutoff_value, file = cutoff_file)

cat("Q3/cutoff value written to:", cutoff_file, "\n")

summary_stats <- data.frame(
  Stat  = c("Q1", "Median", "Q3/cutoff"),
  Value = c(cutoff_25, median_val, cutoff_value)
)

# calculate percentage of genes above cutoff
total_genes    <- sum(!is.na(nodes$Result))
retained_genes <- sum(nodes$Result >= cutoff_value, na.rm = TRUE)
pct_retained   <- 100 * retained_genes / total_genes

# -------------------------
# Reference lines styling
# -------------------------
stat_colors    <- c("Q1" = "#1f78b4", "Median" = "#33a02c", "Q3/cutoff" = "#e31a1c")
stat_linetypes <- c("Q1" = "dashed",   "Median" = "solid",   "Q3/cutoff" = "dotdash")

# Compute histogram bins
binwidth <- 0.25
bins <- ggplot_build(
  ggplot(nodes, aes(x = Result)) +
    geom_histogram(binwidth = binwidth)
)$data[[1]]

# Tallest bar height
y_max <- max(bins$count, na.rm = TRUE)

# For each stat, find the bin it falls into and use that bar height
ref_labs <- summary_stats %>%
  rowwise() %>%
  mutate(
    bin_idx = which(Value >= bins$xmin & Value < bins$xmax)[1],
    bar_height = bins$count[bin_idx],
    y = bar_height + 0.10 * y_max,  # 10% above that bar
    nudge_x = case_when(
      Stat == "Q1"        ~ -0.5,
      Stat == "Median"    ~  0.0,
      Stat == "Q3/cutoff" ~  0.5
    )
  ) %>%
  ungroup()

# -------------------------
# Histogram of D-value cutoff
# -------------------------
p_d_cutoff <- ggplot(nodes, aes(x = Result)) +
  geom_histogram(
    binwidth = 0.25,
    boundary = cutoff_value,
    fill = "#6baed6", color = "#08519c", alpha = 0.6
  ) +
  geom_vline(
    data = summary_stats,
    aes(xintercept = Value, color = Stat, linetype = Stat),
    linewidth = 1
  ) +
  annotate("rect",
           xmin = cutoff_value, xmax = Inf, ymin = 0, ymax = Inf,
           alpha = 0.1, fill = "#e31a1c") +
  geom_label_repel(
    data = ref_labs,
    aes(
      x     = Value,
      y     = y,
      label = paste0(Stat, "\n", round(Value, 3)),
      color = Stat
    ),
    inherit.aes   = FALSE,
    nudge_x       = ref_labs$nudge_x,
    size          = 3.5,
    max.overlaps  = Inf,
    segment.color = "grey50",
    box.padding   = 0.4,
    label.size    = 0,
    fill          = "white",
    direction     = "y",
    force         = 0,
    force_pull    = 0
  ) +
  scale_color_manual(values = stat_colors) +
  scale_linetype_manual(values = stat_linetypes) +
  labs(
    title    = paste("Distribution of gene D-values -", unique_id),
    subtitle = paste0("D-value cutoff = ", round(cutoff_value, 3),
      " | Genes retained = ", retained_genes, "/", total_genes,
      " (", round(pct_retained, 1), "%)"
    ),
    x = "D-value", y = "Gene count"
  ) +
  coord_cartesian(xlim = c(-5, 5)) +
  theme_minimal() +
  theme(legend.position = "none")

# Save plot
output_path <- file.path(out_dir, paste0("d_cutoff_", unique_id, ".png"))
ggsave(output_path, plot = p_d_cutoff, width = 12, height = 6, dpi = 300, bg = "white")

# -------------------------
# Read coincident_pairs.tsv
# -------------------------
pairs <- read.table(pairs_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)

# If no pairs, exit with status 1
if (nrow(pairs) == 0) {
  message("Exiting run ", unique_id, " - no valid pairs in", pairs_path)
  quit(save = "no", status = 1)
}

# Merge d-values for Source and Target genes
pairs <- pairs %>%
  left_join(nodes %>% select(Source = ID, d_source = Result), by = "Source") %>%
  left_join(nodes %>% select(Target = ID, d_target = Result), by = "Target")

# -------------------------
# Histogram of gene-level D-values
# -------------------------

# Calculate outliers above and below D-value of 5 and -5
n_low  <- sum(nodes$Result < -5, na.rm = TRUE)
n_high <- sum(nodes$Result >  5, na.rm = TRUE)

p_hist_genes <- ggplot(nodes, aes(x = Result)) +
  geom_histogram(
    binwidth = 0.25, boundary = cutoff_value,
    fill = "#6baed6", color = "#08519c", alpha = 0.6
  ) +
  geom_vline(
    xintercept = cutoff_value,
    linetype = "dashed", color = "#e31a1c", linewidth = 1
  ) +
    annotate("rect",
           xmin = cutoff_value, xmax = Inf, ymin = 0, ymax = Inf,
           alpha = 0.1, fill = "#e31a1c") +
  scale_x_continuous(
    breaks = seq(-5, 5, by = 1.0)
  ) +
  labs(
    title    = paste("Distribution of gene D-values -", unique_id),
    subtitle = paste0("D-value cutoff = ", round(cutoff_value, 3),
      " | Genes retained = ", retained_genes, "/", total_genes,
      " (", round(pct_retained, 1), "%)"
    ),
    x = "D-value", y = "Gene count"
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
# Histogram of pair-level D-values
# -------------------------

retained_pairs <- sum(pairs$d_source >= cutoff_value & pairs$d_target >= cutoff_value, na.rm = TRUE)
total_pairs    <- nrow(pairs)
pct_retained   <- 100 * retained_pairs / total_pairs

# Pair-level D-value = minimum of the two genes' D-values
pairs$d_pair <- pmin(pairs$d_source, pairs$d_target, na.rm = TRUE)

# Calculate outliers above and below D-value of 5 and -5 at the pair level
n_low_pairs  <- sum(pairs$d_pair < -5, na.rm = TRUE)
n_high_pairs <- sum(pairs$d_pair >  5, na.rm = TRUE)

p_hist_pairs <- ggplot(pairs, aes(x = d_pair)) +
  geom_histogram(
    binwidth = 0.25, boundary = cutoff_value,
    fill = "#fdae6b", color = "#e6550d", alpha = 0.6
  ) +
  geom_vline(
    xintercept = cutoff_value,
    linetype = "dashed", color = "#e31a1c", linewidth = 1
  ) +
  annotate("rect",
           xmin = cutoff_value, xmax = Inf, ymin = 0, ymax = Inf,
           alpha = 0.1, fill = "#e31a1c") +
  scale_x_continuous(
    breaks = seq(-5, 5, by = 1.0)
  ) +
  labs(
    title    = paste("Distribution of gene pairs D-values -", unique_id),
    subtitle = paste0("D-value cutoff = ", round(cutoff_value, 3),
      " | Gene pairs retained = ", retained_pairs, "/", total_pairs,
      " (", round(pct_retained, 1), "%)"
    ),
    x = "D-value",
    y = "Pair count"
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

# -------------------------
# Combine and save plots
# -------------------------
final_plot <- p_hist_genes + p_hist_pairs
output_path <- file.path(out_dir, paste0("d_distribution_", unique_id, ".png"))
ggsave(output_path, plot = final_plot, width = 12, height = 6, dpi = 300, bg = "white")

# -------------------------
# Histogram: number of partners per gene above D-value cutoff
# -------------------------

# Keep only pairs where both genes meet the cutoff (i.e., d_pair >= cutoff_value)
pairs_kept <- pairs %>%
  filter(!is.na(d_pair), d_pair >= cutoff_value)

# Build undirected partner list (count unique partners per gene)
partners_long <- bind_rows(
  pairs_kept %>% transmute(gene = Source, partner = Target),
  pairs_kept %>% transmute(gene = Target, partner = Source)
) %>%
  filter(!is.na(gene), !is.na(partner), gene != partner)

partners_per_gene <- partners_long %>%
  distinct(gene, partner) %>%               # unique partners only
  count(gene, name = "n_partners")          # degree per gene above cutoff

mean_partners <- mean(partners_per_gene$n_partners, na.rm = TRUE)

p_partners <- ggplot(partners_per_gene, aes(x = n_partners)) +
  geom_histogram(binwidth = 1, boundary = 0, closed = "right",
                 fill = "#6baed6", color = "#08519c", alpha = 0.6) +
  labs(
    title = paste("Distribution of gene partners >= D-value cutoff -", unique_id),
    subtitle = paste0(
      "D-value cutoff = ", round(cutoff_value, 3),
      " | Mean gene partners = ", round(mean_partners, 2)
    ),
    x = "Number of partners >= cutoff",
    y = "Gene count"
  ) +
  theme_minimal()

ggsave(file.path(out_dir, paste0("partners_hist_", unique_id, ".png")),
       plot = p_partners, width = 8, height = 6, dpi = 300, bg = "white")

# -------------------------
# Build pairs-retained curve (low → high cutoff, ≥ logic)
# -------------------------
cutoffs <- seq(min(nodes$Result, na.rm = TRUE),
               max(nodes$Result, na.rm = TRUE),
               length.out = 200)

pair_counts <- sapply(cutoffs, function(cut) {
  sum(pairs$d_source >= cut & pairs$d_target >= cut, na.rm = TRUE)
})

df_counts <- data.frame(cutoff = cutoffs, pairs = pair_counts)

# -------------------------
# Pairs-retained curve with elbow cutoff
# -------------------------
p_curve <- ggplot(df_counts, aes(x = cutoff, y = pairs)) +
  geom_point(size = 1.5, color = "black") +
  geom_line(color = "#08519c", linewidth = 1) +
    annotate("rect",
           xmin = cutoff_value, xmax = Inf, ymin = 0, ymax = Inf,
           alpha = 0.1, fill = "#e31a1c") +
  geom_vline(xintercept = cutoff_value, linetype = "dashed", color = "#e31a1c") +
  annotate("text", x = cutoff_value, y = max(df_counts$pairs), 
           label = paste0("Cut-off = ", signif(cutoff_value, 3)), 
           hjust = -0.1, vjust = 1.5, color = "#e31a1c") +
  labs(
    title = paste("Pairs retained vs. D-value cutoff -", unique_id),
    subtitle = paste0("D-value cutoff = ", round(cutoff_value, 3),
      " | Gene pairs retained = ", retained_pairs, "/", total_pairs,
      " (", round(pct_retained, 1), "%)"
    ),
    x = "Lineage independence cutoff (D)",
    y = "# of gene pairs retained"
  ) +
  theme_minimal()
  
ggsave(file.path(out_dir, paste0("d_distribution_pairs_curve_", unique_id, ".png")),
       plot = p_curve, width = 8, height = 6, dpi = 300, bg = "white")

# -------------------------
# D-value cutoff heatmap
# -------------------------

# Load tree
tree <- read.tree(tree_path)

# Load presence/absence matrix
pa_raw <- read_csv(pa_path, show_col_types = FALSE, name_repair = "minimal")

gene_ids <- pa_raw[[1]]
pa_vals  <- as.matrix(pa_raw[,-1])
mode(pa_vals) <- "numeric"

genome_ids <- colnames(pa_raw)[-1]
rownames(pa_vals) <- gene_ids
colnames(pa_vals) <- genome_ids

# Transpose so rows = genomes, cols = genes
pa_mat <- t(pa_vals)

# Align genomes to the tree tip order
keep_genomes <- intersect(rownames(pa_mat), tree$tip.label)
pa_mat <- pa_mat[keep_genomes, , drop = FALSE]
pa_mat <- pa_mat[tree$tip.label, , drop = FALSE]

# Restrict pa_mat to the genes that have D-values
ordered_genes <- intersect(nodes$ID, colnames(pa_mat))
pa_mat <- pa_mat[tree$tip.label, ordered_genes, drop = FALSE]

# Define spread score helper
get_spread_score <- function(vec) {
  sum(diff(as.numeric(vec)) != 0)
}

# Compute spread scores
spread_scores <- apply(pa_mat, 2, get_spread_score)

# Join back into nodes
nodes <- nodes %>%
  left_join(
    data.frame(ID = names(spread_scores),
               spread_score = spread_scores),
    by = "ID"
  )

# Restrict to genes near the cutoff
margin <- 0.14
near_cutoff <- nodes %>%
  filter(Result > (cutoff_value - margin), Result < (cutoff_value + margin))

# Pick one clumped gene below cutoff
gene_lowD <- near_cutoff %>%
  filter(Result < cutoff_value) %>%
  arrange(spread_score) %>%
  slice(15) %>%
  pull(ID)

# Pick one spread gene above cutoff
gene_highD <- near_cutoff %>%
  filter(Result > cutoff_value) %>%
  arrange(desc(spread_score)) %>%
  slice(1) %>%
  pull(ID)

cat("Chosen below-cutoff gene:", gene_lowD,
    "D =", nodes$Result[nodes$ID == gene_lowD], "\n")
cat("Chosen above-cutoff gene:", gene_highD,
    "D =", nodes$Result[nodes$ID == gene_highD], "\n")

# Build tidy data frames with .id column = tip labels
df_low <- data.frame(
  .id   = rownames(pa_mat),
  value = as.character(pa_mat[, gene_lowD])
)

df_high <- data.frame(
  .id   = rownames(pa_mat),
  value = as.character(pa_mat[, gene_highD])
)

# Labels with D-values
panel_low  <- paste0(gene_lowD,  " (D=", round(nodes$Result[nodes$ID == gene_lowD], 3), ")")
panel_high <- paste0(gene_highD, " (D=", round(nodes$Result[nodes$ID == gene_highD], 3), ")")

# Base tree
p_tree2 <- ggtree(tree)

# Add each gene as its own facet panel
two_gene_plot <- facet_plot(p_tree2, panel = panel_low,
                            data = df_low,
                            geom = geom_tile,
                            mapping = aes(x = 1, fill = value),
                            color = "white", size = 1.5)

two_gene_plot <- facet_plot(two_gene_plot, panel = panel_high,
                            data = df_high,
                            geom = geom_tile,
                            mapping = aes(x = 1, fill = value),
                            color = "white", size = 1.5) +
  scale_fill_manual(values = c("0" = "white", "1" = "black")) +
  labs(
    title = paste("Gene distribution patterns above and below D-value cutoff -", unique_id),
    subtitle = paste("D-value cutoff =", round(cutoff_value, 3), 
                     "| example genes below and above the D-value cutoff")
  ) +
  theme(legend.position = "none")

# Save
outfile2 <- file.path(out_dir, paste0("gene_distribution_heatmap_", unique_id, ".png"))
ggsave(outfile2, plot = two_gene_plot, width = 8, height = 6)

# Combined plot of cutoff, elbow and heatmap
combined_plot <- (p_d_cutoff | p_hist_pairs) / (p_curve | two_gene_plot)

# Save to file
outfile <- file.path(out_dir, paste0("d_combined_", unique_id, ".png"))
ggsave(outfile, combined_plot, width = 14, height = 10)

# -------------------------
# Write output data for comparative plots
# -------------------------
# Gene-level data: keep only genes above cutoff
df_genes <- nodes %>%
  filter(!is.na(Result), Result >= cutoff_value) %>%
  transmute(
    Gene   = ID,
    Gene_1 = NA_character_,
    Gene_2 = NA_character_,
    D_value = Result,
    Method  = "Coinfinder",
    Level   = "Gene"
  )

# Pair-level data: keep only pairs where both genes meet cutoff
pairs$d_pair <- pmin(pairs$d_source, pairs$d_target, na.rm = TRUE)

df_pairs <- pairs %>%
  filter(!is.na(d_pair), d_source >= cutoff_value, d_target >= cutoff_value) %>%
  transmute(
    Gene   = NA_character_,
    Gene_1 = Source,
    Gene_2 = Target,
    D_value = d_pair,
    Method  = "Coinfinder",
    Level   = "Pair"
  )

# Combine
df_out <- bind_rows(df_genes, df_pairs)

# Write to CSV
write_csv(df_out, file.path(out_dir, paste0("coinfinder_dvalues_", unique_id, ".csv")))
message("D-value CSV saved to: ", file.path(out_dir, paste0("coinfinder_dvalues_", unique_id, ".csv")))
