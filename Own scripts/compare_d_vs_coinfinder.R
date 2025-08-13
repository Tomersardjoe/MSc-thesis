#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

read_coin_nodes <- function(path) {
  dt <- fread(path)
  if (nrow(dt) == 0) stop("Empty nodes.tsv: ", path)
  
  if (!all(c("ID", "Result") %in% names(dt))) {
    stop("Expected columns 'ID' and 'Result' not found in nodes.tsv. Found: ",
         paste(names(dt), collapse = ", "))
  }
  
  dt[, .(gene = as.character(ID), D_coin = as.numeric(Result))]
}

read_ours <- function(path) {
  dt <- fread(path)
  if (!("gene" %in% names(dt) && "D" %in% names(dt))) {
    stop("Our CSV must contain 'gene' and 'D' columns.")
  }
  setnames(copy(dt[, .(gene, D)]), "D", "D_ours")
}

compare_d <- function(coin_nodes, ours_csv, dataset = NA_character_, out_dir = ".", save_plot = TRUE) {
  cf <- read_coin_nodes(coin_nodes)
  us <- read_ours(ours_csv)
  m <- merge(us, cf, by = "gene", all = FALSE)
  m <- m[is.finite(D_ours) & is.finite(D_coin)]
  
  if (!nrow(m)) stop("No overlapping genes with finite D values.")
  
  pearson  <- cor(m$D_ours, m$D_coin, method = "pearson")
  spearman <- cor(m$D_ours, m$D_coin, method = "spearman")
  mae <- mean(abs(m$D_ours - m$D_coin))
  rmse <- sqrt(mean((m$D_ours - m$D_coin)^2))
  bias <- mean(m$D_ours - m$D_coin)
  
  summary <- data.table(
    Dataset = ifelse(is.na(dataset), basename(dirname(coin_nodes)), dataset),
    N_shared = nrow(m),
    Pearson_r = round(pearson, 3),
    Spearman_rho = round(spearman, 3),
    MAE = round(mae, 3),
    RMSE = round(rmse, 3),
    Mean_Bias = round(bias, 3)
  )
  
  if (save_plot) {
    if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
    plt <- ggplot(m, aes(D_coin, D_ours)) +
      geom_abline(slope = 1, intercept = 0, linetype = 2, color = "grey50") +
      geom_point(alpha = 0.6) +
      labs(x = "Coinfinder D", y = "Our D", 
           title = paste0("D-value agreement: ", summary$Dataset[1]),
           subtitle = sprintf("N=%d, r=%.3f, rho=%.3f, MAE=%.3f, RMSE=%.3f",
                              summary$N_shared, summary$Pearson_r, summary$Spearman_rho, summary$MAE, summary$RMSE)) +
      theme_minimal(base_size = 12)
    ggsave(file.path(out_dir, sprintf("dvalue_agreement_%s.png", summary$Dataset[1])),
           plt, width = 6, height = 5, dpi = 150)
  }
  
  list(summary = summary, merged = m)
}

# CLI usage: Rscript compare_d_vs_coinfinder.R <nodes.tsv> <our_dvalues.csv> [out_dir]
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  cat("\nUsage: Rscript compare_d_vs_coinfinder.R <nodes.tsv> <our_dvalues.csv> [out_dir]\n\n")
  quit(status = 1)
}
coin_nodes <- args[1]
ours_csv   <- args[2]
out_dir <- ifelse(length(args) >= 3, args[3], ".")

res <- compare_d(coin_nodes, ours_csv, dataset = basename(dirname(coin_nodes)), out_dir = out_dir, save_plot = TRUE)
fwrite(res$summary, file.path(out_dir, sprintf("dvalue_agreement_%s.csv", res$summary$Dataset[1])), quote = FALSE)
cat(paste(res$summary, collapse = ","), "\n")
