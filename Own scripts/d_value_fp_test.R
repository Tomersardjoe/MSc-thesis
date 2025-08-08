library(ape)
library(phangorn)

compute_D <- function(tree, trait, n_sim = 1000, seed = seed, na_drop = TRUE) {
  set.seed(seed)
  if (!inherits(tree, "phylo")) stop("tree must be a phylo object")
  if (is.null(names(trait))) stop("trait must be a named vector with taxa names")
  
  # Align to tip order
  trait <- trait[match(tree$tip.label, names(trait))]
  names(trait) <- tree$tip.label
  
  # Optionally drop NAs and prune tree
  if (na_drop && anyNA(trait)) {
    keep <- which(!is.na(trait))
    tree <- drop.tip(tree, setdiff(tree$tip.label, tree$tip.label[keep]))
    trait <- trait[keep]
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
                mean_brownian = NA_real_, random_vals = NA, brownian_vals = NA))
  }
  
  # Build phyDat from a matrix (most stable path)
  make_phyDat <- function(tr) {
    tr <- as.character(tr)
    if (!all(tr %in% c("0","1"))) stop("Internal: tr must be '0'/'1'")
    mat <- matrix(tr, ncol = 1,
                  dimnames = list(tree$tip.label, "trait"))
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
  
  # 2) Random permutations
  random_vals <- replicate(n_sim, {
    trait_perm <- sample(trait, replace = FALSE)
    trait_changes(trait_perm)
  })
  mean_random <- mean(random_vals)
  
  # 3) Brownian sims (preserve prevalence via thresholding to k1)
  brownian_vals <- replicate(n_sim, {
    cont <- rTraitCont(tree)
    cont <- cont[tree$tip.label]
    bin <- bm_threshold(cont, k = k1)
    trait_changes(bin)
  })
  mean_brownian <- mean(brownian_vals)
  
  # 4) D-value
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
                           seed = seed,
                           na_drop = TRUE,
                           taxa_col = NULL,           # set if species/taxon names are in a column
                           drop_monomorphic = TRUE,   # skip all-0 or all-1 traits
                           return_distributions = FALSE) {
  
  stopifnot(inherits(tree, "phylo"))
  
  # 1) Prepare trait matrix with taxa as rownames
  if (!is.null(taxa_col)) {
    stopifnot(taxa_col %in% colnames(traits))
    rn <- traits[[taxa_col]]
    traits <- traits[ , setdiff(colnames(traits), taxa_col), drop = FALSE]
    rownames(traits) <- rn
  }
  if (is.null(rownames(traits))) stop("Traits must have taxa as rownames or provide `taxa_col`.")
  
  # 2) Ensure all tree tips are present and in order
  missing <- setdiff(tree$tip.label, rownames(traits))
  if (length(missing) > 0) stop(sprintf("Traits are missing %d taxa, e.g., %s", length(missing), paste(head(missing, 3), collapse = ", ")))
  traits <- traits[tree$tip.label, , drop = FALSE]
  
  # 3) Helpers
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
    # deterministic per-gene seed from base R only
    (seed + sum(utf8ToInt(as.character(g))) * 1009L) %% .Machine$integer.max
  }
  
  # 4) Iterate per gene
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
        note = "monomorphic"
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
        note = "non-binary"
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
  }
  
  res_df <- do.call(rbind, out)
  rownames(res_df) <- NULL
  res_df
}

verify_D <- function(tree, trait_vec,
                     n_sim = 1000,
                     seed = seed,
                     tolerance = 1e-6,
                     verbose = TRUE) {
  stopifnot(inherits(tree, "phylo"))
  if (is.null(names(trait_vec))) stop("Trait vector must be named with tree$tip.label")
  
  # Align and validate
  trait_vec <- trait_vec[tree$tip.label]
  names(trait_vec) <- tree$tip.label
  if (!all(trait_vec %in% c("0", "1"))) stop("Trait must be binary '0'/'1'")
  if (length(unique(trait_vec)) <= 1) stop("Trait is monomorphic")
  
  # Helpers — mirror compute_D
  make_phyDat <- function(tr) {
    mat <- matrix(tr, ncol = 1,
                  dimnames = list(tree$tip.label, "trait"))  # force tree$tip.label
    phyDat(mat, type = "USER", levels = c("0", "1"))
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
  
  # Observed
  observed_manual <- trait_changes(trait_vec)
  
  # Use one RNG stream to mirror compute_D
  set.seed(seed)
  
  # Random permutations — identical to compute_D
  random_vals <- replicate(n_sim, {
    trait_perm <- sample(trait_vec, replace = FALSE)
    trait_changes(trait_perm)
  })
  mean_random_manual <- mean(random_vals)
  
  # Brownian sims — identical to compute_D (no reseeding)
  k1 <- sum(trait_vec == "1")
  brownian_vals <- replicate(n_sim, {
    cont <- rTraitCont(tree)
    cont <- cont[tree$tip.label]
    bin <- bm_threshold(cont, k = k1)
    trait_changes(bin)
  })
  mean_brownian_manual <- mean(brownian_vals)
  
  # D calculation — same orientation as compute_D
  denom <- mean_random_manual - mean_brownian_manual
  D_manual <- if (abs(denom) < .Machine$double.eps) NA_real_
  else (observed_manual - mean_brownian_manual) / denom
  
  # Reference from compute_D
  ref <- compute_D(tree, trait_vec, n_sim = n_sim, seed = seed)
  
  # Comparison
  delta <- function(x, y) if (any(is.na(c(x, y)))) NA else abs(x - y)
  delta_D   <- delta(ref$D, D_manual)
  delta_obs <- delta(ref$observed, observed_manual)
  delta_rand <- delta(ref$mean_random, mean_random_manual)
  delta_bm <- delta(ref$mean_brownian, mean_brownian_manual)
  
  if (verbose) {
    cat("\n--- D-Value Verification ---\n")
    cat("Observed:        ", observed_manual, "vs", ref$observed, "\n")
    cat("Mean Random:     ", round(mean_random_manual, 4), "vs", round(ref$mean_random, 4), "\n")
    cat("Mean Brownian:   ", round(mean_brownian_manual, 4), "vs", round(ref$mean_brownian, 4), "\n")
    cat("D-value (manual):", round(D_manual, 6), "vs", round(ref$D, 6), "\n")
    cat("Δ D:             ", round(delta_D, 6), "\n")
  }
  
  data.frame(
    D_ref = ref$D,
    D_manual = D_manual,
    delta_D = delta_D,
    observed_ref = ref$observed,
    observed_manual = observed_manual,
    delta_observed = delta_obs,
    mean_random_ref = ref$mean_random,
    mean_random_manual = mean_random_manual,
    delta_random = delta_rand,
    mean_brownian_ref = ref$mean_brownian,
    mean_brownian_manual = mean_brownian_manual,
    delta_brownian = delta_bm,
    consistent = ifelse(is.na(delta_D), NA, delta_D <= tolerance)
  )
}





