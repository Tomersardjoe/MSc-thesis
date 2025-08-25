#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(igraph)
})

# --- Arguments ---
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 6) {
  stop("Usage: script.R <gene_presence_absence.csv - root directory> <output directory> <coinfinder root directory> <goldfinder root directory> <panforest root directory> <dataset>")
}

gpa_root_dir        <- args[1]
output_dir          <- args[2]
coinfinder_root_dir <- args[3]
goldfinder_root_dir <- args[4]
panforest_root_dir  <- args[5]
dataset             <- args[6]

deny_dirs <- c("__pycache__")

list_nonhidden_dirs <- function(root) {
  entries <- list.files(root, full.names = TRUE, recursive = FALSE, include.dirs = TRUE, no.. = TRUE)
  if (length(entries) == 0) return(character(0))
  dirs <- entries[file.info(entries)$isdir %in% TRUE]
  keep <- !(basename(dirs) %in% deny_dirs) & !grepl("^\\.", basename(dirs))
  dirs[keep]
}

search_dirs <- function(root, dataset) {
  root <- normalizePath(root, mustWork = TRUE)
  ds_dir <- file.path(root, dataset)
  if (!dir.exists(ds_dir)) {
    stop("Dataset '", dataset, "' not found in root: ", root)
  }
  normalizePath(ds_dir, mustWork = TRUE)
}

find_exactly_one_file <- function(root, pattern, dataset) {
  dirs <- search_dirs(root, dataset)
  files <- unlist(lapply(dirs, function(d) {
    list.files(d, pattern = pattern, full.names = TRUE, recursive = TRUE)
  }), use.names = FALSE)
  files <- unique(files)
  if (length(files) != 1) {
    stop("Expected exactly one match for pattern '", pattern, "' under root/dataset: ",
         root, " / ", dataset,
         ". Found (", length(files), "): ", paste(files, collapse = " | "))
  }
  normalizePath(files, mustWork = TRUE)
}

# --- Weight mapping ---
WEIGHT_ATTR_BY_TOOL <- c(
  coinfinder = "FDR_BH",
  goldfinder = "Force",
  panforest = "Importance"
)

detect_tool <- function(file_path) {
  ext <- tolower(tools::file_ext(file_path))
  
  # Panforest
  if (ext %in% c("graphml", "xml"))
    return(list(tool = "panforest", format = "graphml"))
  
  # Goldfinder CSV (header-based)
  if (ext == "csv") {
    hdr <- names(read.csv(file_path, nrows = 1, check.names = FALSE))
    if (all(c("Node1","Node2","Force","pair_type") %in% hdr))
      return(list(tool = "goldfinder", format = "csv_goldfinder"))
    return(list(tool = "csv_unknown", format = "csv"))
  }
  
  # Coinfinder TSV (header-based)
  if (ext == "tsv") {
    hdr <- names(read.delim(file_path, nrows = 1, check.names = FALSE))
    if (all(c("Source","Target","weight","FDR_BH") %in% hdr))
      return(list(tool = "coinfinder", format = "tsv_coinfinder"))
    return(list(tool = "tsv_unknown", format = "tsv"))
  }
  
  return(list(tool = "unknown", format = ext))
}

read_native_graph <- function(file_path, info) {
  if (info$format == "graphml") {
    g <- read_graph(file_path, format = "graphml")
  } else if (info$format == "csv_goldfinder") {
    df <- read.csv(file_path, check.names = FALSE, stringsAsFactors = FALSE)
    edges <- df[, c("Node1","Node2","Force","pair_type"), drop = FALSE]
    g <- graph_from_data_frame(edges, directed = FALSE)
  } else if (info$format == "tsv_coinfinder") {
    df <- read.delim(file_path, check.names = FALSE, stringsAsFactors = FALSE)
    edges <- df[, c("Source","Target","weight",WEIGHT_ATTR_BY_TOOL["coinfinder"]), drop = FALSE]
    g <- graph_from_data_frame(edges, directed = FALSE)
  } else {
    stop("Unsupported or unknown input format for: ", file_path)
  }
  g
}

normalize_graph <- function(g, tool = NULL) {
  weight_attr <- unname(WEIGHT_ATTR_BY_TOOL[[tool]])
  if (is.null(weight_attr) || !(weight_attr %in% edge_attr_names(g))) {
    weight_attr <- NULL
  }
  attr_names <- edge_attr_names(g)
  attr_comb <- if (length(attr_names)) setNames(rep("first", length(attr_names)), attr_names) else list()
  if (!is.null(weight_attr)) attr_comb[[weight_attr]] <- "mean"
  g <- as_undirected(g, mode = "collapse", edge.attr.comb = attr_comb)
  g <- simplify(g, remove.multiple = TRUE, remove.loops = TRUE, edge.attr.comb = attr_comb)
  attr(g, "weight_attr") <- weight_attr
  g
}

harmonize_edge_weights <- function(
    g, tool, alpha = 1, concave_exp = 0.7,
    eps_p = 1e-300, winsor = c(0, 0.99),
    min_cost = 1e-9
) {
  # 1) Tool-specific strength
  s <- switch(tool,
              coinfinder = {
                p <- as.numeric(igraph::edge_attr(g, "FDR_BH"))
                -log10(pmax(p, eps_p))
              },
              goldfinder = {
                f <- as.numeric(igraph::edge_attr(g, "Force"))
                if (is.null(igraph::edge_attr(g, "pair_type"))) {
                  E(g)$association <- ifelse(is.na(f), NA_character_,
                                             ifelse(f >= 0, "co", "anti"))
                } else {
                  E(g)$association <- as.character(igraph::edge_attr(g, "pair_type"))
                }
                abs(f)
              },
              panforest = {
                imp <- as.numeric(igraph::edge_attr(g, "Importance"))
                pmax(imp, 0)
              },
              stop("Unknown tool: ", tool)
  )
  
  # 2) Handle NAs/Infs
  s[!is.finite(s)] <- NA_real_
  
  # 3) Winsorization (stabilize extremes)
  if (!is.null(winsor) && length(winsor) == 2) {
    qs <- stats::quantile(s, probs = winsor, na.rm = TRUE, names = FALSE)
    s <- pmin(pmax(s, qs[1]), qs[2])
  }
  
  # 4) ECDF percentile scaling to [0,1]
  valid <- which(!is.na(s))
  w <- rep(NA_real_, length(s))
  if (length(valid) > 1) {
    ec <- stats::ecdf(s[valid])
    w[valid] <- ec(s[valid])
  } else if (length(valid) == 1) {
    w[valid] <- 1.0
  }
  
  # 5) Tail emphasis
  if (!is.null(alpha) && alpha > 0 && alpha != 1) {
    w <- w ^ alpha
  }
  
  # 6) Attach harmonized weights and concave distance transform
  E(g)$w_harmonized <- as.numeric(w)
  E(g)$distance <- 1 - (pmin(pmax(E(g)$w_harmonized, 0), 1))^concave_exp
  
  # 7) Ensure strictly positive distances for igraph algorithms
  E(g)$distance <- pmax(E(g)$distance, min_cost)
  E(g)$weight   <- E(g)$distance
  
  g
}

restrict_to_canonical <- function(g, canonical_genes, drop_noncanonical = TRUE) {
  current <- V(g)$name
  if (is.null(current)) stop("Vertex names are missing; graphs must have gene IDs as vertex names.")
  noncanon <- setdiff(current, canonical_genes)
  if (drop_noncanonical && length(noncanon) > 0) {
    g <- induced_subgraph(g, vids = intersect(current, canonical_genes))
  }
  list(
    g = g,
    noncanonical_removed = length(noncanon)
  )
}

pad_to_canonical <- function(g, canonical_genes) {
  missing <- setdiff(canonical_genes, V(g)$name)
  if (length(missing) > 0) {
    g <- add_vertices(g, nv = length(missing), name = missing)
  }
  g
}

load_harmonize_graph <- function(file_path, canonical_genes, drop_noncanonical = TRUE,
                                 alpha = 1) {
  info <- detect_tool(file_path)
  g_raw <- read_native_graph(file_path, info)
  g_norm <- normalize_graph(g_raw, tool = info$tool)
  present_in_canonical <- sum(V(g_norm)$name %in% canonical_genes)
  coverage_ratio <- present_in_canonical / length(canonical_genes)
  r <- restrict_to_canonical(g_norm, canonical_genes, drop_noncanonical = drop_noncanonical)
  g_restricted <- r$g
  nodes_pre_padding <- vcount(g_restricted)
  g_padded <- pad_to_canonical(g_restricted, canonical_genes)
  g_h <- harmonize_edge_weights(g_padded, tool = info$tool, alpha = alpha)
  list(
    tool = info$tool,
    format = info$format,
    graph = g_h,
    coverage_ratio = coverage_ratio,
    noncanonical_removed = r$noncanonical_removed,
    nodes_pre_padding = nodes_pre_padding
  )
}

compute_metrics <- function(g, use_weights = TRUE) {
  w <- NULL
  if (use_weights) {
    # Prefer harmonized weights if available
    if ("w_harmonized" %in% edge_attr_names(g)) {
      w <- E(g)$w_harmonized
    } else {
      weight_attr <- attr(g, "weight_attr")
      if (!is.null(weight_attr) &&
          ecount(g) > 0 &&
          weight_attr %in% edge_attr_names(g)) {
        vals <- tryCatch(E(g)[[weight_attr]], error = function(e) NULL)
        if (!is.null(vals)) {
          if (!is.numeric(vals)) vals <- suppressWarnings(as.numeric(vals))
          if (length(vals) == ecount(g) && any(is.finite(vals))) {
            w <- vals
          }
        }
      }
    }
  }
  n <- vcount(g); m <- ecount(g)
  deg <- if (n > 0) degree(g, mode = "all") else numeric(0)
  
  cl <- NULL
  if (m > 0 && n > 0) {
    cl <- tryCatch(cluster_louvain(g, weights = w), error = function(e) NULL)
  }
  
  isolates <- if (n > 0) sum(deg == 0) else 0
  comp <- if (n > 0) components(g) else list(csize = integer())
  giant_frac <- if (n > 0 && length(comp$csize)) max(comp$csize) / n else NA_real_
  
  safe <- function(expr) tryCatch(expr, error = function(e) NA_real_)
  
  possible_gene_pairs <- if (n > 1) n * (n - 1) / 2 else 0
  association_rate <- if (possible_gene_pairs > 0) m / possible_gene_pairs else NA_real_
  
  module_count <- if (!is.null(cl)) length(unique(membership(cl))) else NA_real_
  comm_sizes <- sizes(cl)
  module_count <- sum(comm_sizes > 1)
  avg_genes_per_module <- if (module_count > 0) {
    sum(comm_sizes[comm_sizes > 1]) / module_count
  } else NA_real_
  avg_path_length <- if (m > 0) safe(mean_distance(g, directed = FALSE, unconnected = TRUE)) else NA_real_
  diameter_val <- if (m > 0) safe(diameter(g, directed = FALSE, unconnected = TRUE)) else NA_real_
  comp_entropy <- if (n > 0) {
    probs <- sizes(comp) / sum(sizes(comp))
    -sum(probs * log2(probs))
  } else NA_real_
  
  avg_local_clust <- if (m > 0) safe(mean(transitivity(g, type = "local", isolates = "zero"), na.rm = TRUE)) else NA_real_
  kcore_max <- if (n > 0) max(coreness(g)) else NA_real_
  
  degree_gini <- if (n > 0) ineq::Gini(deg) else NA_real_
  
  btw <- if (n > 0) betweenness(g) else numeric(0)
  btw_mean <- if (length(btw)) mean(btw) else NA_real_
  btw_max <- if (length(btw)) max(btw) else NA_real_
  eigen_var <- if (n > 0) var(eigen_centrality(g)$vector) else NA_real_
  
  comm_entropy <- if (!is.null(cl)) {
    probs <- comm_sizes / sum(comm_sizes)
    -sum(probs * log2(probs))
  } else NA_real_
  inter_edges_frac <- if (!is.null(cl)) {
    memb <- membership(cl)
    inter_edges <- E(g)[memb[ends(g, E(g))[,1]] != memb[ends(g, E(g))[,2]]]
    length(inter_edges) / gsize(g)
  } else NA_real_
  
  articulation_pts <- if (n > 0) length(articulation_points(g)) else NA_real_
  hub_tol <- if (n > 0) {
    frac_remove <- 0.05
    hubs <- order(deg, decreasing = TRUE)[1:ceiling(frac_remove * vcount(g))]
    g2 <- delete_vertices(g, hubs)
    max(sizes(components(g2))) / vcount(g)
  } else NA_real_
  
  # --- Motif enrichment: triad census ---
  if (!is_directed(g)) {
    tri <- sum(count_triangles(g)) / 3
    triples <- sum(choose(degree(g), 2))
    open <- triples - tri
    total_triples <- choose(vcount(g), 3)
    disconnected <- total_triples - open - tri
    
    triad_census_vals <- rep(NA_real_, 16)
    triad_census_vals[1]  <- disconnected # 003
    triad_census_vals[2]  <- open         # 102
    triad_census_vals[16] <- tri          # 300
  } else {
    triad_census_vals <- if (m > 0) triad_census(g) else rep(NA_real_, 16)
  }
  
  triad_names <- paste0("triad_", sprintf("%02d", seq_along(triad_census_vals)))
  names(triad_census_vals) <- triad_names
  
  cliques_size3 <- if (m > 0) length(cliques(g, min = 3, max = 3)) else NA_real_
  
  data.frame(
    nodes = n,
    edges = m,
    isolates = isolates,
    coverage_ratio = if (n > 0) 1 - isolates / n else NA_real_,
    avg_degree = if (n > 0) mean(deg) else NA_real_,
    global_clustering_coefficient = if (m > 1) safe(transitivity(g, type = "globalundirected", weights = w)) else NA_real_,
    assortativity_degree = if (m > 0) safe(assortativity_degree(g, directed = FALSE)) else NA_real_,
    modularity_louvain = if (!is.null(cl)) safe(modularity(cl)) else NA_real_,
    giant_component_frac = giant_frac,
    possible_gene_pairs = possible_gene_pairs,
    association_rate = association_rate,
    module_count = module_count,
    avg_genes_per_module = avg_genes_per_module,
    avg_path_length = avg_path_length,
    diameter = diameter_val,
    component_entropy = comp_entropy,
    avg_local_clust = avg_local_clust,
    kcore_max = kcore_max,
    degree_gini = degree_gini,
    btw_mean = btw_mean,
    btw_max = btw_max,
    eigen_var = eigen_var,
    comm_entropy = comm_entropy,
    intercomm_frac = inter_edges_frac,
    articulation_pts = articulation_pts,
    hub_tol = hub_tol,
    cliques_size3 = cliques_size3,
    as.list(triad_census_vals),
    stringsAsFactors = FALSE
  )
}

# ==== Main execution ====
gpa_path <- find_exactly_one_file(
  root    = gpa_root_dir,
  pattern = "^gene_presence_absence\\.csv$",
  dataset = dataset
)

canonical <- read.csv(gpa_path, stringsAsFactors = FALSE)
if (!("Gene" %in% names(canonical))) {
  stop("The gene_presence_absence.csv must contain a 'Gene' column.")
}
canonical_genes <- unique(canonical$Gene)

coinfinder_path <- find_exactly_one_file(
  root    = coinfinder_root_dir,
  pattern = "\\_edges.tsv$",
  dataset = dataset
)

goldfinder_path <- find_exactly_one_file(
  root    = goldfinder_root_dir,
  pattern = "cytoscape_input\\.csv$",
  dataset = dataset
)

panforest_path <- find_exactly_one_file(
  root    = panforest_root_dir,
  pattern = "\\.graphml$",
  dataset = dataset
)

paths <- c(coinfinder_path, goldfinder_path, panforest_path)

loaded <- lapply(paths, load_harmonize_graph, canonical_genes = canonical_genes, alpha = 1)
names(loaded) <- vapply(loaded, function(x) x$tool, character(1))
graphs <- lapply(loaded, function(x) x$graph)

metrics <- do.call(rbind, lapply(names(graphs), function(name) {
  m <- compute_metrics(graphs[[name]])
  m$tool <- name
  m$dataset <- dataset
  m$format <- loaded[[name]]$format
  m$nodes_pre_padding <- loaded[[name]]$nodes_pre_padding
  m$noncanonical_removed <- loaded[[name]]$noncanonical_removed
  m$coverage_ratio <- loaded[[name]]$coverage_ratio
  m
}))

metrics <- metrics[, c(
  "tool", "dataset", "format", "nodes_pre_padding",
  "nodes", "edges", "isolates", "possible_gene_pairs",
  "coverage_ratio", "association_rate", "avg_degree",
  "avg_path_length", "diameter", "component_entropy",
  "giant_component_frac",
  "global_clustering_coefficient", "avg_local_clust", "assortativity_degree", "kcore_max",
  "module_count", "avg_genes_per_module", "comm_entropy",
  "modularity_louvain", "intercomm_frac",
  "degree_gini", "btw_mean", "btw_max", "eigen_var",
  "articulation_pts", "hub_tol",
  "cliques_size3",
  "triad_01", "triad_02", "triad_16",
  "noncanonical_removed"
)]

# Mark metrics computed from harmonized weights
norm_cols <- c(
  "global_clustering_coefficient",
  "avg_path_length",
  "diameter",
  "modularity_louvain",
  "btw_mean",
  "btw_max",
  "eigen_var"
)

for (col in norm_cols) {
  if (col %in% names(metrics)) {
    names(metrics)[names(metrics) == col] <- paste0(col, "_norm")
  }
}

if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
metrics_output_path <- file.path(output_dir, paste0(dataset, "_metrics.csv"))
write.csv(metrics, metrics_output_path, row.names = FALSE)
message("Metrics written to: ", metrics_output_path)
