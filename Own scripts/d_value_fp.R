#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  # Core deps
  library(ape)
  library(phangorn)
  # CLI
  library(optparse)
})

# -------------------------
# Utility: tidy errors and booleans
# -------------------------

fail <- function(msg, code = 1) {
  message(sprintf("Error: %s", msg))
  quit(save = "no", status = code, runLast = FALSE)
}

parse_bool <- function(x, default = TRUE) {
  if (is.null(x) || is.na(x)) return(default)
  x <- tolower(trimws(as.character(x)))
  if (x %in% c("true", "t", "1", "yes", "y")) return(TRUE)
  if (x %in% c("false", "f", "0", "no", "n")) return(FALSE)
  default
}

# -------------------------
# Core D implementation
# -------------------------

compute_D <- function(tree, trait, n_sim = 1000, seed = 1, na_drop = TRUE) {
  set.seed(seed)
  if (!inherits(tree, "phylo")) stop("tree must be a phylo object")
  if (is.null(names(trait))) stop("trait must be a named vector with taxa names")
  # Align trait to tip order
  trait <- trait[match(tree$tip.label, names(trait))]
  names(trait) <- tree$tip.label
  
  # Optionally drop NAs and prune tree
  if (na_drop && anyNA(trait)) {
    keep <- which(!is.na(trait))
    if (length(keep) <= 1) {
      warning("Too few non-NA taxa after pruning.")
      return(list(D = NA_real_, observed = NA_real_, mean_random = NA_real_,
                  mean_brownian = NA_real_, random_vals = NA, brownian_vals = NA,
                  n_tips = length(keep)))
    }
    tree <- drop.tip(tree, setdiff(tree$tip.label, tree$tip.label[keep]))
    trait <- trait[keep]
    names(trait) <- tree$tip.label
  }
  if (anyNA(trait)) stop("Trait contains NA after alignment; set na_drop=TRUE to prune.")
  
  # Coerce to binary "0"/"1"
  trait <- as.character(trait)
  if (!all(trait %in% c("0","1"))) {
    suppressWarnings(num <- as.numeric(trait))
    if (!any(is.na(num)) && all(num %in% c(0,1))) {
      trait <- as.character(num)
    } else {
      stop("Trait must be binary 0/1 (as numeric or '0'/'1').")
    }
  }
  
  # Prevalence and invariant guard
  k1 <- sum(trait == "1")
  if (k1 == 0 || k1 == length(trait)) {
    warning("Trait is invariant (all 0 or all 1). D is undefined.")
    return(list(D = NA_real_, observed = NA_real_, mean_random = NA_real_,
                mean_brownian = NA_real_, random_vals = NA, brownian_vals = NA,
                n_tips = length(tree$tip.label)))
  }
  
  make_phyDat <- function(tr) {
    tr <- as.character(tr)
    if (!all(tr %in% c("0","1"))) stop("Internal: tr must be '0'/'1'")
    mat <- matrix(tr, ncol = 1, dimnames = list(tree$tip.label, "trait"))
    phyDat(mat, type = "USER", levels = c("0","1"))
  }
  trait_changes <- function(tr) {
    states <- make_phyDat(tr)
    as.numeric(parsimony(tree, states, method = "fitch"))
  }
  bm_threshold <- function(x, k) {
    ord <- order(x, decreasing = TRUE)
    out <- rep("0", length(x))
    if (k > 0) out[ord[seq_len(k)]] <- "1"
    names(out) <- names(x)
    out
  }
  
  # 1) Observed
  observed <- trait_changes(trait)
  # 2) Random
  random_vals <- replicate(n_sim, {
    trait_perm <- sample(trait, replace = FALSE)
    trait_changes(trait_perm)
  })
  mean_random <- mean(random_vals)
  # 3) Brownian
  brownian_vals <- replicate(n_sim, {
    cont <- rTraitCont(tree)
    cont <- cont[tree$tip.label]
    bin <- bm_threshold(cont, k = k1)
    trait_changes(bin)
  })
  mean_brownian <- mean(brownian_vals)
  # 4) Fritz & Purvis D: 0=Brownian, 1=Random
  denom <- mean_random - mean_brownian
  D <- if (abs(denom) < .Machine$double.eps) NA_real_ else (observed - mean_brownian) / denom
  
  list(
    D = D,
    observed = observed,
    mean_random = mean_random,
    mean_brownian = mean_brownian,
    random_vals = random_vals,
    brownian_vals = brownian_vals,
    n_tips = length(tree$tip.label)
  )
}

compute_D_many <- function(tree,
                           traits,
                           n_sim = 1000,
                           seed = 1,
                           na_drop = TRUE,
                           drop_monomorphic = TRUE,
                           return_distributions = FALSE) {
  stopifnot(inherits(tree, "phylo"))
  if (is.null(rownames(traits))) stop("Traits must have taxa as rownames.")
  # Ensure all tree tips are present and in order
  missing <- setdiff(tree$tip.label, rownames(traits))
  if (length(missing) > 0) {
    stop(sprintf("Traits are missing %d taxa, e.g., %s",
                 length(missing), paste(head(missing, 3), collapse = ", ")))
  }
  traits <- traits[tree$tip.label, , drop = FALSE]
  
  coerce_char01 <- function(x) {
    x_chr <- as.character(x)
    if (!all(x_chr %in% c("0","1"))) {
      suppressWarnings(num <- as.numeric(x_chr))
      if (!any(is.na(num)) && all(num %in% c(0,1))) x_chr <- as.character(num)
    }
    x_chr
  }
  is_valid_binary <- function(x) {
    x_chr <- coerce_char01(x)
    ux <- unique(stats::na.omit(x_chr))
    length(ux) == 2 && all(ux %in% c("0","1"))
  }
  is_monomorphic <- function(x) {
    x_chr <- coerce_char01(x)
    ux <- unique(stats::na.omit(x_chr))
    length(ux) <= 1
  }
  gene_seed <- function(g) {
    (seed + sum(utf8ToInt(as.character(g))) * 1009L) %% .Machine$integer.max
  }
  
  genes <- colnames(traits)
  out <- vector("list", length(genes))
  
  for (i in seq_along(genes)) {
    g <- genes[i]
    x <- traits[, g]
    note <- NA_character_
    
    if (drop_monomorphic && is_monomorphic(x)) {
      out[[i]] <- data.frame(
        gene = g, D = NA_real_, observed = NA_real_,
        mean_random = NA_real_, mean_brownian = NA_real_,
        prevalence = mean(as.numeric(coerce_char01(x)) == 1, na.rm = TRUE),
        n_tips = length(tree$tip.label),
        n_sim = n_sim,
        note = "monomorphic",
        stringsAsFactors = FALSE
      )
      next
    }
    if (!is_valid_binary(x)) {
      out[[i]] <- data.frame(
        gene = g, D = NA_real_, observed = NA_real_,
        mean_random = NA_real_, mean_brownian = NA_real_,
        prevalence = NA_real_,
        n_tips = length(tree$tip.label),
        n_sim = n_sim,
        note = "non-binary",
        stringsAsFactors = FALSE
      )
      next
    }
    
    tr_vec <- setNames(coerce_char01(x), rownames(traits))
    res <- compute_D(tree = tree,
                     trait = tr_vec,
                     n_sim = n_sim,
                     seed = gene_seed(g),
                     na_drop = na_drop)
    if (!return_distributions) {
      res$random_vals <- NULL
      res$brownian_vals <- NULL
    }
    
    out[[i]] <- data.frame(
      gene = g,
      D = res$D,
      observed = res$observed,
      mean_random = res$mean_random,
      mean_brownian = res$mean_brownian,
      prevalence = mean(as.numeric(tr_vec) == 1, na.rm = TRUE),
      n_tips = res$n_tips,
      n_sim = n_sim,
      note = ifelse(is.na(note), NA_character_, note),
      stringsAsFactors = FALSE
    )
    cat(sprintf("[%s] ✅ Finished gene: %s (%d sims)\n", Sys.time(), g, n_sim))
    flush.console()
  }
  res_df <- do.call(rbind, out)
  rownames(res_df) <- NULL
  res_df
}

# -------------------------
# I/O and harmonization
# -------------------------

read_gene_list <- function(path) {
  x <- readLines(path, warn = FALSE)
  x <- trimws(x)
  x <- x[nchar(x) > 0]
  unique(x)
}

# Prepare traits matrix with taxa as rownames and genes as columns
# Assumes presence/absence values are 0/1 (numeric or character)
prepare_traits_matrix <- function(pa_path, tree,
                                  taxa_mode = c("auto","rows","cols"),
                                  gene_col = NA,
                                  sample_cols = NULL) {
  taxa_mode <- match.arg(tolower(taxa_mode), c("auto","rows","cols"))
  raw <- utils::read.csv(pa_path, check.names = FALSE, stringsAsFactors = FALSE)
  tips <- tree$tip.label
  
  # Helper: coerce values to "0"/"1"
  coerce01_df <- function(df) {
    df[] <- lapply(df, function(x) {
      xc <- as.character(x)
      if (!all(xc %in% c("0","1"))) {
        suppressWarnings(num <- as.numeric(xc))
        if (!any(is.na(num)) && all(num %in% c(0, 1))) xc <- as.character(num)
      }
      xc
    })
    df
  }
  
  if (taxa_mode == "auto") {
    # Heuristic: more matches in columns => taxa are columns (common for PA with samples as columns)
    col_overlap <- sum(colnames(raw) %in% tips)
    row_overlap <- if (!is.null(rownames(raw))) sum(rownames(raw) %in% tips) else 0
    taxa_mode <- if (col_overlap >= row_overlap && col_overlap > 0) "cols" else "rows"
  }
  
  if (taxa_mode == "cols") {
    # Genes as rows, samples as columns
    # Identify gene ID column if present
    if (is.na(gene_col)) {
      gcand <- intersect(c("Gene", "gene", "GeneID", "locus_tag", "gene_name"), colnames(raw))
      gene_col <- if (length(gcand) >= 1) gcand[1] else NA
    }
    if (!is.na(gene_col) && gene_col %in% colnames(raw)) {
      gene_ids <- raw[[gene_col]]
      pa <- raw[, setdiff(colnames(raw), gene_col), drop = FALSE]
      rownames(pa) <- gene_ids
    } else {
      pa <- raw
      if (is.null(rownames(pa))) {
        stop("Could not identify a gene ID column; supply --gene_col when taxa_mode='cols'.")
      }
    }
    
    if (is.null(sample_cols)) {
      sample_cols <- intersect(colnames(pa), tips)
      if (length(sample_cols) == 0) {
        stop("No sample columns match tree tip labels. Provide --sample_cols or harmonize names.")
      }
    } else {
      sample_cols <- intersect(sample_cols, colnames(pa))
      if (length(sample_cols) == 0) {
        stop("Provided --sample_cols did not match any columns in presence/absence file.")
      }
    }
    
    pa <- pa[, sample_cols, drop = FALSE]
    pa <- coerce01_df(pa)
    
    # Transpose to taxa x genes
    pa_t <- t(pa)  # taxa x genes
    traits <- as.data.frame(pa_t, check.names = FALSE, stringsAsFactors = FALSE)
    rownames(traits) <- rownames(pa_t)
    colnames(traits) <- colnames(pa_t)
    
  } else {
    # taxa_mode == "rows": taxa are rows; columns are genes
    traits <- raw
    if (is.null(rownames(traits))) {
      stop("Traits must have taxa as rownames when taxa_mode='rows'.")
    }
    traits <- traits[intersect(rownames(traits), tips), , drop = FALSE]
    traits <- coerce01_df(traits)
  }
  
  # Ensure all tree tips exist in traits
  missing <- setdiff(tips, rownames(traits))
  if (length(missing) > 0) {
    stop(sprintf("Traits are missing %d taxa (e.g., %s). Harmonize names or supply --sample_cols.",
                 length(missing), paste(head(missing, 3), collapse = ", ")))
  }
  traits <- traits[tips, , drop = FALSE]
  traits
}

# -------------------------
# Pipeline runner
# -------------------------

run_d_pipeline <- function(tree_path,
                           pa_path,
                           gold_path,
                           coin_path,
                           out_prefix = "d_values",
                           out_dir = ".",
                           n_sim = 2000,
                           base_seed = 1,
                           taxa_mode = "auto",
                           gene_col = NA,
                           sample_cols = NULL,
                           na_drop = TRUE) {
  # Create output directory if it does not exist
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

  # Validate output directory and prefix
  if (is.null(out_dir) || is.na(out_dir) || !nzchar(out_dir)) {
  fail("Invalid output directory: out_dir is missing or empty.")
  }
  if (is.null(out_prefix) || is.na(out_prefix) || !nzchar(out_prefix)) {
  fail("Invalid output prefix: out_prefix is missing or empty.")
  }

  # Load inputs
  if (!file.exists(tree_path)) fail(sprintf("Tree file not found: %s", tree_path))
  if (!file.exists(pa_path)) fail(sprintf("Presence/absence file not found: %s", pa_path))
  if (!file.exists(gold_path)) fail(sprintf("Goldfinder genes file not found: %s", gold_path))
  if (!file.exists(coin_path)) fail(sprintf("Coinfinder genes file not found: %s", coin_path))
  
  tree <- read.tree(tree_path)
  if (is.null(tree$tip.label) || length(tree$tip.label) < 3) {
    fail("Tree must have at least 3 tips with labels.")
  }
  
  gf_genes <- read_gene_list(gold_path)
  cf_genes <- read_gene_list(coin_path)
  
  traits <- prepare_traits_matrix(pa_path, tree,
                                  taxa_mode = taxa_mode,
                                  gene_col = gene_col,
                                  sample_cols = sample_cols)
  
  # Filter to requested genes (keep order of provided lists)
  gf_present <- gf_genes[gf_genes %in% colnames(traits)]
  cf_present <- cf_genes[cf_genes %in% colnames(traits)]
  gf_missing <- setdiff(gf_genes, gf_present)
  cf_missing <- setdiff(cf_genes, cf_present)
  
  if (length(gf_missing)) message(sprintf("Goldfinder: %d genes not found (e.g., %s)",
                                          length(gf_missing), paste(head(gf_missing, 5), collapse = ", ")))
  if (length(cf_missing)) message(sprintf("Coinfinder: %d genes not found (e.g., %s)",
                                          length(cf_missing), paste(head(cf_missing, 5), collapse = ", ")))
  
  # Subset matrices
  traits_gf <- if (length(gf_present)) traits[, gf_present, drop = FALSE] else traits[, 0, drop = FALSE]
  traits_cf <- if (length(cf_present)) traits[, cf_present, drop = FALSE] else traits[, 0, drop = FALSE]
  
  # Compute blocks
  compute_block <- function(traits_block) {
    if (ncol(traits_block) == 0L) {
      return(data.frame(gene = character(), D = numeric(), observed = numeric(),
                        mean_random = numeric(), mean_brownian = numeric(),
                        prevalence = numeric(), n_tips = integer(), n_sim = integer(),
                        note = character(), stringsAsFactors = FALSE))
    }
    compute_D_many(tree, traits_block, n_sim = n_sim, seed = base_seed, na_drop = na_drop)
  }
  
  res_gf <- compute_block(traits_gf)
  res_cf <- compute_block(traits_cf)
  
  # Order back to input list order
  if (nrow(res_gf)) res_gf <- res_gf[match(gf_present, res_gf$gene), , drop = FALSE]
  if (nrow(res_cf)) res_cf <- res_cf[match(cf_present, res_cf$gene), , drop = FALSE]
  
  # Add source label and write outputs
  if (nrow(res_gf)) res_gf$source <- "goldfinder"
  if (nrow(res_cf)) res_cf$source <- "coinfinder"
  res_all <- rbind(res_gf, res_cf)

  utils::write.csv(res_gf, file.path(out_dir, sprintf("%s_goldfinder.csv", out_prefix)), row.names = FALSE)
  utils::write.csv(res_cf, file.path(out_dir, sprintf("%s_coinfinder.csv", out_prefix)), row.names = FALSE)
  utils::write.csv(res_all, file.path(out_dir, sprintf("%s_combined.csv", out_prefix)), row.names = FALSE)

  invisible(list(goldfinder = res_gf, coinfinder = res_cf, combined = res_all,
                 missing = list(goldfinder = gf_missing, coinfinder = cf_missing)))
}

# -------------------------
# CLI: parse args and run
# -------------------------

option_list <- list(
  make_option(c("--tree", "-t"), type = "character", help = "Path to tree file (.nwk)", metavar = "FILE"),
  make_option(c("--presence_absence", "-p"), type = "character", help = "Path to gene presence–absence file (.csv)", metavar = "FILE"),
  make_option(c("--goldfinder", "-g"), type = "character", help = "Path to Goldfinder gene list (.txt)", metavar = "FILE"),
  make_option(c("--coinfinder", "-c"), type = "character", help = "Path to Coinfinder gene list (.txt)", metavar = "FILE"),
  make_option(c("--out_prefix", "-o"), type = "character", default = "d_values", help = "Prefix for output CSVs [default %default]", metavar = "STR"),
  make_option(c("--out_dir", "-d"), type = "character", default = ".", help = "Directory to write output CSVs [default %default]", metavar = "DIR"),
  make_option(c("--n_sim", "-n"), type = "integer", default = 4000, help = "Number of simulations per gene [default %default]", metavar = "INT"),
  make_option(c("--seed", "-s"), type = "integer", default = 1, help = "Base random seed [default %default]", metavar = "INT"),
  make_option(c("--taxa_mode"), type = "character", default = "auto", help = "Taxa orientation in PA file: auto|rows|cols [default %default]", metavar = "STR"),
  make_option(c("--gene_col"), type = "character", default = NA, help = "Gene ID column in PA when taxa_mode=cols (e.g., Gene)", metavar = "STR"),
  make_option(c("--sample_cols"), type = "character", default = NA, help = "Comma-separated sample columns to use (must match tree tips)", metavar = "STR"),
  make_option(c("--na_drop"), type = "character", default = "true", help = "Prune taxa with NA per gene: true|false [default %default]", metavar = "BOOL")
)

opt_parser <- OptionParser(usage = "%prog [options]",
                           description = "Compute Fritz & Purvis D for Goldfinder and Coinfinder gene sets.",
                           option_list = option_list)
opt <- parse_args(opt_parser)

# Validate required args
required <- c("tree","presence_absence","goldfinder","coinfinder")
missing_req <- required[sapply(required, function(k) is.null(opt[[k]]))]
if (length(missing_req) > 0) {
  print_help(opt_parser)
  fail(sprintf("Missing required arguments: %s", paste(missing_req, collapse = ", ")))
}

# Parse sample_cols and na_drop
sample_cols <- NULL
if (!is.null(opt$sample_cols) && !is.na(opt$sample_cols)) {
  sample_cols <- strsplit(opt$sample_cols, ",")[[1]]
  sample_cols <- trimws(sample_cols)
  sample_cols <- sample_cols[nchar(sample_cols) > 0]
  if (length(sample_cols) == 0) sample_cols <- NULL
}

na_drop <- parse_bool(opt$na_drop, default = TRUE)

# Run
res <- tryCatch({
  run_d_pipeline(
    tree_path     = opt$tree,
    pa_path       = opt$presence_absence,
    gold_path     = opt$goldfinder,
    coin_path     = opt$coinfinder,
    out_prefix    = opt$out_prefix,
    out_dir       = opt$out_dir,
    n_sim         = as.integer(opt$n_sim),
    base_seed     = as.integer(opt$seed),
    taxa_mode     = tolower(opt$taxa_mode),
    gene_col      = if (is.na(opt$gene_col)) NA else opt$gene_col,
    sample_cols   = sample_cols,
    na_drop       = na_drop
  )
}, error = function(e) {
  fail(conditionMessage(e))
})

# Console summary
cat(sprintf("\n✅ D-values written to:\n - %s/%s_goldfinder.csv\n - %s/%s_coinfinder.csv\n - %s/%s_combined.csv\n",
            opt$out_dir, opt$out_prefix, opt$out_dir, opt$out_prefix, opt$out_dir, opt$out_prefix))

invisible(res)
