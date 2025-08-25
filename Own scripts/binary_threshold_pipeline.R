#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(tidyverse)
})

get_tool_paths <- function(tool, dataset_dir, truth_file) {
  switch(
    tool,
    Coinfinder = list(
      coin = file.path(dataset_dir, "coincident_pairs.tsv"),
      truth = truth_file
    ),
    Goldfinder = list(
      assoc = file.path(dataset_dir, "simultaneous_association_significant_pairs.csv"),
      diss  = file.path(dataset_dir, "simultaneous_dissociation_significant_pairs.csv"),
      truth = truth_file
    ),
    PanForest = list(
      pan   = file.path(dataset_dir, "imp.csv"),
      truth = truth_file
    ),
    stop("Unknown tool: ", tool)
  )
}

canonical_pair <- function(a, b) {
  ifelse(a < b, paste0(a, "_", b), paste0(b, "_", a))
}

pairs_from_adj_matrix <- function(M) {
  g <- rownames(M)
  ut <- which(upper.tri(M) & M == 1, arr.ind = TRUE)
  if (nrow(ut) == 0) return(character(0))
  canonical_pair(g[ut[, 1]], g[ut[, 2]])
}

jaccard_set <- function(set_a, set_b) {
  u <- union(set_a, set_b)
  if (length(u) == 0) return(NA_real_)
  length(intersect(set_a, set_b)) / length(u)
}

compare_to_truth <- function(tool_pairs, truth_pairs) {
  tp <- length(intersect(tool_pairs, truth_pairs))
  fp <- length(setdiff(tool_pairs, truth_pairs))
  fn <- length(setdiff(truth_pairs, tool_pairs))
  precision <- if ((tp + fp) == 0) NA_real_ else tp / (tp + fp)
  recall    <- if ((tp + fn) == 0) NA_real_ else tp / (tp + fn)
  f1        <- if (is.na(precision) || is.na(recall) || (precision + recall) == 0) NA_real_ else 2 * precision * recall / (precision + recall)
  list(tp = tp, fp = fp, fn = fn, precision = precision, recall = recall, f1 = f1)
}

run_analysis_single <- function(tool, paths) {
  truth_df <- read_tsv(
    paths$truth,
    col_names = FALSE,
    col_types = cols(.default = col_character())
  )
  colnames(truth_df) <- c("Gene1", "Gene2")
  truth_pairs <- canonical_pair(truth_df$Gene1, truth_df$Gene2)

  tool_pairs <- character(0)

  if (tool == "Coinfinder") {
    coin <- read_tsv(
      paths$coin,
      col_types = cols(
        Source = col_character(),
        Target = col_character(),
        .default = col_double()
      )
    )
    coin <- select(coin, Source, Target)
    tool_pairs <- canonical_pair(coin$Source, coin$Target)

  } else if (tool == "Goldfinder") {
    gold_assoc <- if (file.exists(paths$assoc)) {
      select(read_csv(paths$assoc, col_types = cols(
        Gene_1 = col_character(),
        Gene_2 = col_character()
      )), Gene_1, Gene_2)
    } else {
      tibble(Gene_1 = character(), Gene_2 = character())
    }

    gold_diss <- if (file.exists(paths$diss)) {
      select(read_csv(paths$diss, col_types = cols(
        Gene_1 = col_character(),
        Gene_2 = col_character()
      )), Gene_1, Gene_2)
    } else {
      tibble(Gene_1 = character(), Gene_2 = character())
    }

    gold <- bind_rows(gold_assoc, gold_diss)
    tool_pairs <- canonical_pair(gold$Gene_1, gold$Gene_2)

  } else if (tool == "PanForest") {
    pan <- suppressMessages(
      read_csv(
        paths$pan,
        col_types = cols(
          .default = col_double(),
          `...1` = col_character()
        )
      )
    )
    names(pan)[1] <- "Genes"
    mat <- as.matrix(pan[,-1])
    rownames(mat) <- pan$Genes
    upper_vals <- mat[upper.tri(mat)]
    thr <- suppressWarnings(quantile(upper_vals, 0.95, na.rm = TRUE))
    M <- ifelse(mat >= thr, 1, 0)
    diag(M) <- 0
    tool_pairs <- pairs_from_adj_matrix(M)
  }

  tool_pairs <- unique(na.omit(tool_pairs))
  perf <- compare_to_truth(tool_pairs, truth_pairs)

  tibble(
    tool = tool,
    n_edges = length(tool_pairs),
    tp = perf$tp,
    fp = perf$fp,
    fn = perf$fn,
    precision = perf$precision,
    recall = perf$recall,
    f1 = perf$f1,
    edges = list(tool_pairs)
  )
}

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4) {
  stop("Usage: Rscript binary_threshold_pipeline.R <coinfinder_root> <goldfinder_root> <panforest_root> <simulation_root>")
}

tool_roots <- c(
  Coinfinder = args[1],
  Goldfinder = args[2],
  PanForest  = args[3]
)
simulation_root <- args[4]

types <- c("closed", "moderate", "open")
reps  <- 1:10

per_tool_results <- list()
pairwise_results <- list()

for (type in types) {
  for (i in reps) {
    truth_file <- file.path(simulation_root, "output", type, "grid_1", paste0("rep", i), "simulation", "pairs.txt")
    edges_by_tool <- list()
    present_tools <- character(0)

    for (tool in names(tool_roots)) {
      dataset_dir <- file.path(tool_roots[[tool]], type, paste0("rep", i))
      if (!dir.exists(dataset_dir)) next
      paths <- get_tool_paths(tool, dataset_dir, truth_file)
      if (!all(file.exists(unlist(paths)))) next

      res <- run_analysis_single(tool, paths) %>%
        mutate(type = type, replicate = i, dataset_dir = dataset_dir)

      per_tool_results[[length(per_tool_results) + 1]] <- res
      present_tools <- c(present_tools, tool)
      edges_by_tool[[tool]] <- res$edges[[1]]
    }

    if (length(present_tools) >= 2) {
      combos <- combn(present_tools, 2, simplify = FALSE)
      for (cmb in combos) {
        t1 <- cmb[1]; t2 <- cmb[2]
        e1 <- edges_by_tool[[t1]]
        e2 <- edges_by_tool[[t2]]
        j  <- jaccard_set(e1, e2)
        pairwise_results[[length(pairwise_results) + 1]] <- tibble(
          type = type,
          replicate = i,
          tool_a = t1,
          tool_b = t2,
          shared = length(intersect(e1, e2)),
          total = length(union(e1, e2)),
          jaccard = j
        )
      }
    }
  }
}

per_tool_df <- bind_rows(per_tool_results) %>%
  arrange(tool, type, replicate)

pairwise_df <- bind_rows(pairwise_results) %>%
  mutate(
    tool_a_new = pmin(tool_a, tool_b),
    tool_b_new = pmax(tool_a, tool_b)
  ) %>%
  select(-tool_a, -tool_b) %>%
  rename(tool_a = tool_a_new, tool_b = tool_b_new) %>%
  arrange(tool_a, tool_b, type, replicate)

write_csv(per_tool_df, "per_tool_metrics.csv")
write_csv(pairwise_df, "pairwise_jaccard.csv")

cat("Saved per-tool metrics to per_tool_metrics.csv\n")
cat("Saved cross-tool Jaccard to pairwise_jaccard.csv\n")
