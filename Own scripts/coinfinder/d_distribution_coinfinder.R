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
  stop("Usage: Rscript d_distribution_coinfinder.R /path/to/coincident_nodes.tsv /path/to/coincident_pairs.tsv /path/to/tree.nwk /path/to/presence_absence.csv")
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
# Read nodes.tsv (gene-level d-values)
# -------------------------
nodes <- read.table(nodes_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
nodes$Result <- as.numeric(nodes$Result)

if (all(is.na(nodes$Result))) {
    message("Skipping run ", unique_id, " - no valid D-values in nodes file.")
    next
}

# -------------------------
# Read coincident_pairs.tsv
# -------------------------
pairs <- read.table(pairs_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)

# Merge d-values for Source and Target genes
pairs <- pairs %>%
  left_join(nodes %>% select(Source = ID, d_source = Result), by = "Source") %>%
  left_join(nodes %>% select(Target = ID, d_target = Result), by = "Target")

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
# Elbow detection on descending curve
# -------------------------
all_points <- cbind(1:nrow(df_counts), df_counts$pairs)
first_point <- all_points[1, ]
last_point  <- all_points[nrow(df_counts), ]
line_vec <- last_point - first_point
line_vec <- line_vec / sqrt(sum(line_vec^2))
vec_from_first <- sweep(all_points, 2, first_point)
scalar_proj <- vec_from_first %*% line_vec
proj_point <- scalar_proj %*% line_vec + matrix(rep(first_point, nrow(df_counts)), nrow = nrow(df_counts), byrow = TRUE)
distances <- sqrt(rowSums((all_points - proj_point)^2))
elbow_idx <- which.max(distances)
cutoff_value <- df_counts$cutoff[elbow_idx]

cutoff_file <- file.path(out_dir, paste0(unique_id, "_d_cutoff.txt"))
write.table(cutoff_value, file = cutoff_file, row.names = FALSE, col.names = FALSE, quote = FALSE)

# -------------------------
# Histogram of gene-level d-values
# -------------------------
p_hist_genes <- ggplot(nodes, aes(x = Result)) +
  geom_histogram(binwidth = 0.25, boundary = cutoff_value, fill = "#6baed6", color = "#08519c", alpha = 0.6) +
  geom_vline(xintercept = cutoff_value, linetype = "dashed", color = "#e31a1c", linewidth = 1) +
  scale_x_continuous(
    breaks = seq(
      floor(min(nodes$Result, na.rm = TRUE)),
      ceiling(max(nodes$Result, na.rm = TRUE)),
      by = 0.25
    )
  ) +
  labs(
    title = paste("Distribution of gene D-values —", unique_id),
    x = "D-value", y = "Gene count"
  ) +
  theme_minimal()

# -------------------------
# Histogram of pair-level d-values
# -------------------------

# Pair-level D-value = minimum of the two genes' D-values
pairs$d_pair <- pmin(pairs$d_source, pairs$d_target, na.rm = TRUE)

p_hist_pairs <- ggplot(pairs, aes(x = d_pair)) +
  geom_histogram(binwidth = 0.25, boundary = cutoff_value, fill = "#fdae6b", color = "#e6550d", alpha = 0.6) +
  geom_vline(xintercept = cutoff_value, linetype = "dashed", color = "#e31a1c", linewidth = 1) +
  scale_x_continuous(
    breaks = seq(
      floor(min(nodes$Result, na.rm = TRUE)),
      ceiling(max(nodes$Result, na.rm = TRUE)),
      by = 0.25
    )
  ) +
  labs(
    title = paste("Distribution of gene pairs D-values —", unique_id),
    x = "D-value (pair-level, min of source & target)",
    y = "Pair count"
  ) +
  theme_minimal()

# -------------------------
# Combine and save plots
# -------------------------
final_plot <- p_hist_genes + p_hist_pairs
output_path <- file.path(out_dir, paste0("d_distribution_", unique_id, ".png"))
ggsave(output_path, plot = final_plot, width = 12, height = 6, dpi = 300, bg = "white")

message("Plot saved as: ", output_path)
message("Recommended Coinfinder cut-off (pairs curve elbow): ", cutoff_value)

# -------------------------
# Pairs-retained curve with elbow cutoff
# -------------------------
p_curve <- ggplot(df_counts, aes(x = cutoff, y = pairs)) +
  geom_point(size = 1.5, color = "black") +
  geom_line(color = "#08519c", linewidth = 1) +
  geom_vline(xintercept = cutoff_value, linetype = "dashed", color = "#e31a1c") +
  annotate("text", x = cutoff_value, y = max(df_counts$pairs), 
           label = paste0("Cut-off = ", signif(cutoff_value, 3)), 
           hjust = -0.1, vjust = 1.5, color = "#e31a1c") +
  labs(
    title = "Pairs retained vs. D-value cutoff",
    x = "Lineage independence cutoff (D)",
    y = "# of coincident gene pairs"
  ) +
  theme_minimal()
  
ggsave(file.path(out_dir, paste0("d_distribution_pairs_curve_", unique_id, ".png")),
       plot = p_curve, width = 8, height = 6, dpi = 300, bg = "white")

# -------------------------
# Cross-check with phylogeny + clustering metric for ALL genes
# -------------------------
tree <- read.tree(tree_path)

# Read presence/absence CSV (genes as rows, taxa as columns)
pa_raw <- read.csv(pa_path, header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
gene_ids <- pa_raw[[1]]
pa_raw <- pa_raw[ , -1, drop = FALSE]

# Transpose so taxa are rows, genes are columns
pa <- as.data.frame(t(pa_raw))
colnames(pa) <- gene_ids
pa$Taxon <- rownames(pa)
rownames(pa) <- pa$Taxon
pa$Taxon <- NULL

genes_above <- nodes$ID[nodes$Result >= cutoff_value]
genes_below <- nodes$ID[nodes$Result < cutoff_value]

set.seed(123)
sample_above <- sample(genes_above, min(5, length(genes_above)))
sample_below <- sample(genes_below, min(5, length(genes_below)))

plot_gene_on_tree <- function(gene_id) {
  if (!gene_id %in% colnames(pa)) return(NULL)
  presence <- as.numeric(pa[[gene_id]])
  anno <- data.frame(label = rownames(pa), presence = presence)
  ggtree(tree) %<+% anno +
    geom_tippoint(aes(color = factor(presence)), size = 2) +
    scale_color_manual(values = c("0" = "grey70", "1" = "red")) +
    ggtitle(gene_id)
}

# Function to calculate parsimony score (lower = more clumped)
calc_parsimony <- function(gene_id) {
  if (!gene_id %in% colnames(pa)) return(NA)
  pres <- as.numeric(pa[[gene_id]])
  names(pres) <- rownames(pa)
  phyDat_obj <- phyDat(as.matrix(pres), type = "USER", levels = c(0,1))
  score <- parsimony(tree, phyDat_obj)
  return(score)
}

# Calculate metrics
metrics_above <- data.frame(
  Gene = sample_above,
  ParsimonyScore = sapply(sample_above, calc_parsimony)
)
metrics_below <- data.frame(
  Gene = sample_below,
  ParsimonyScore = sapply(sample_below, calc_parsimony)
)

# Save plots to PDF
pdf(file.path(out_dir, paste0("phylo_crosscheck_", unique_id, ".pdf")), width = 8, height = 10)
for (p in lapply(sample_above, plot_gene_on_tree)) if (!is.null(p)) print(p + ggtitle(paste("Above cutoff:", signif(cutoff_value, 3))))
for (p in lapply(sample_below, plot_gene_on_tree)) if (!is.null(p)) print(p + ggtitle(paste("Below cutoff:", signif(cutoff_value, 3))))
invisible(dev.off())

# Save metrics summary
metrics_summary <- rbind(
  data.frame(Group = "Above cutoff", metrics_above),
  data.frame(Group = "Below cutoff", metrics_below)
)
write.table(metrics_summary,
            file = file.path(out_dir, paste0("phylo_crosscheck_tree_metrics_", unique_id, ".tsv")),
            sep = "\t", row.names = FALSE, quote = FALSE)

message("Phylogeny cross-check PDF saved as: ",
        file.path(out_dir, paste0("phylo_crosscheck_", unique_id, ".pdf")))
message("Clustering metrics saved as: ",
        file.path(out_dir, paste0("phylo_crosscheck_metrics_", unique_id, ".tsv")))

# Calculate scores for ALL genes
nodes$ParsimonyScore <- sapply(nodes$ID, calc_parsimony)

# Split into groups
nodes$Group <- ifelse(nodes$Result >= cutoff_value, "Above cutoff", "Below cutoff")

# Summary stats
summary_stats <- nodes %>%
  group_by(Group) %>%
  summarise(
    Count = n(),
    Min = min(ParsimonyScore, na.rm = TRUE),
    Q1 = quantile(ParsimonyScore, 0.25, na.rm = TRUE),
    Median = median(ParsimonyScore, na.rm = TRUE),
    Mean = mean(ParsimonyScore, na.rm = TRUE),
    Q3 = quantile(ParsimonyScore, 0.75, na.rm = TRUE),
    Max = max(ParsimonyScore, na.rm = TRUE)
  )

# Save summary table
write.table(summary_stats,
            file = file.path(out_dir, paste0("phylo_crosscheck_summary_", unique_id, ".tsv")),
            sep = "\t", row.names = FALSE, quote = FALSE)

# -------------------------
# Statistical test: Wilcoxon rank-sum (Mann–Whitney U)
# -------------------------
scores_above <- nodes$ParsimonyScore[nodes$Group == "Above cutoff"]
scores_below <- nodes$ParsimonyScore[nodes$Group == "Below cutoff"]

wilcox_res <- wilcox.test(scores_above, scores_below, alternative = "two.sided", exact = FALSE)

# Save p-value to a text file
pval_path <- file.path(out_dir, paste0("phylo_crosscheck_stats_", unique_id, ".txt"))
writeLines(
  c(
    paste("Wilcoxon rank-sum test comparing parsimony scores above vs. below cutoff"),
    paste("Cutoff value:", signif(cutoff_value, 4)),
    paste("W statistic:", wilcox_res$statistic),
    paste("p-value:", signif(wilcox_res$p.value, 4))
  ),
  con = pval_path
)

message("Statistical test results saved as: ", pval_path)

# -------------------------
# Boxplot/violin plot of parsimony scores by group
# -------------------------
p_box <- ggplot(nodes, aes(x = Group, y = ParsimonyScore, fill = Group)) +
  geom_violin(trim = FALSE, alpha = 0.5) +
  geom_boxplot(width = 0.1, outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.1, size = 0.5) +
  scale_fill_manual(values = c("Above cutoff" = "#1f78b4", "Below cutoff" = "#e31a1c")) +
  labs(
    title = "Parsimony score distribution by group",
    x = "",
    y = "Parsimony score (lower = more clumped)"
  ) +
  theme_minimal()

ggsave(file.path(out_dir, paste0("phylo_crosscheck_boxplot_", unique_id, ".png")),
       plot = p_box, width = 6, height = 5, dpi = 300, bg = "white")

message("Clustering summary table saved as: ",
        file.path(out_dir, paste0("phylo_crosscheck_summary_", unique_id, ".tsv")))
message("Clustering distribution plot saved as: ",
        file.path(out_dir, paste0("phylo_crosscheck_boxplot_", unique_id, ".png")))

# -------------------------
# Write output data for comparative plots
# -------------------------

# Gene-level data
df_genes <- nodes %>%
  transmute(
    Gene = ID,
    Gene_1 = NA_character_,
    Gene_2 = NA_character_,
    D_value = Result,
    Method = "Coinfinder",
    Level = "Gene"
  ) %>%
  filter(!is.na(D_value))

# Pair-level data
pairs$d_pair <- pmin(pairs$d_source, pairs$d_target, na.rm = TRUE)

df_pairs <- pairs %>%
  transmute(
    Gene = NA_character_,
    Gene_1 = Source,
    Gene_2 = Target,
    D_value = d_pair,
    Method = "Coinfinder",
    Level = "Pair"
  ) %>%
  filter(!is.na(D_value))

# Combine
df_out <- bind_rows(df_genes, df_pairs)

# Write to CSV
write_csv(df_out, file.path(out_dir, paste0("coinfinder_dvalues_", unique_id, ".csv")))

message("D-value CSV saved to: ",
        file.path(out_dir, paste0("coinfinder_dvalues_", unique_id, ".csv")))
