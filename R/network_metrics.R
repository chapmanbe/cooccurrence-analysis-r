# Network metrics and community detection.

#' Modularity Q of a partition of a weighted graph
#'
#' `Q = sum_c [L_c / m - (D_c / 2m)^2]`, where `m` is the total edge weight,
#' `L_c` the weight of edges inside community `c` (each counted once), and
#' `D_c` the total strength of its nodes. The `(D_c / 2m)^2` term carries the
#' null-model penalty for every same-community pair, edge or not, and the
#' diagonal self-terms, so a one-community partition of any graph scores 0.
#' Returns 0 for an empty or edgeless graph (igraph returns `NaN` there).
#' Computed by [igraph::modularity()].
#'
#' `adjacency` is a symmetric, non-negative weight matrix with a zero
#' diagonal; anything else raises. `assignments` is one community label per
#' row (any values).
#' @export
modularity_q <- function(adjacency, assignments) {
  n <- nrow(adjacency)
  stopifnot(length(assignments) == n)
  if (n == 0L) return(0)
  valid <- !anyNA(adjacency) && isSymmetric(unname(adjacency)) && all(adjacency >= 0) &&
    all(diag(adjacency) == 0)
  if (!valid) {
    stop("modularity_q: adjacency must be symmetric, non-negative, NA-free, with a zero diagonal",
         call. = FALSE)
  }
  if (all(adjacency == 0)) return(0)
  g <- igraph::graph_from_adjacency_matrix(adjacency, mode = "undirected", weighted = TRUE,
                                           diag = FALSE)
  igraph::modularity(g, membership = match(assignments, unique(assignments)),
                     weights = igraph::E(g)$weight)
}

#' Detect communities in a co-occurrence network
#'
#' `method = "leiden"` (default) runs [igraph::cluster_leiden()] with the
#' modularity objective, iterated to convergence; `"louvain"` runs
#' [igraph::cluster_louvain()]. Both use the edge weights and `resolution`.
#' Leiden is the default because it guarantees connected communities, which
#' Louvain can fail to give (Traag, Waltman & van Eck 2019, Sci Rep 9:5233). The Julia reference's `"label_propagation"` is not
#' ported; asking for it, or anything else, raises.
#'
#' **Divergence from the Julia reference:** Julia runs only Louvain's first
#' phase (local moving, no aggregation), which can stall below the optimum,
#' so partitions differ from Julia's by design (CLAUDE.md, "Divergences").
#'
#' Both algorithms are randomized. `seed` has no default: pass a whole number
#' for a reproducible partition (applied locally; the caller's RNG stream is
#' restored), or `NULL` to draw from the current stream. A graph with no
#' edges gets one community per vertex without calling igraph.
#'
#' Returns a list: `assignments` (integer, named by item, in `net$items`
#' order, labels `1..K` numbered by first appearance), `communities` (a list
#' named `"1"`, `"2"`, ..., each the sorted member items), `modularity`
#' ([modularity_q()] of the partition), and `n_communities`.
#' @export
detect_communities <- function(net, method = "leiden", seed, resolution = 1) {
  stopifnot(inherits(net, "cooccurrence_network"))
  if (identical(method, "label_propagation")) {
    stop(paste0("method = 'label_propagation' exists in the Julia reference but is not ",
                "ported to R: it is nondeterministic, and the analysis uses 'leiden'"),
         call. = FALSE)
  }
  .check_choice(method, c("leiden", "louvain"), "method")
  if (missing(seed)) {
    stop(paste0("detect_communities: seed is required: pass a whole number for a ",
                "reproducible partition, or NULL to draw from the current RNG stream"),
         call. = FALSE)
  }
  .check_seed(seed, "detect_communities")
  if (!(is.numeric(resolution) && length(resolution) == 1L && is.finite(resolution) &&
          resolution > 0)) {
    stop("detect_communities: resolution must be one positive number", call. = FALSE)
  }
  items <- net[["items"]]
  if (length(items) == 0L) {
    return(list(assignments = stats::setNames(integer(0), character(0)),
                communities = list(), modularity = 0, n_communities = 0L))
  }
  adj <- net[["adjacency"]]
  raw <- if (all(adj == 0)) {
    seq_along(items)
  } else {
    g <- igraph::graph_from_adjacency_matrix(adj, mode = "undirected", weighted = TRUE,
                                             diag = FALSE)
    w <- igraph::E(g)$weight
    fit <- .with_seed(seed, if (identical(method, "leiden")) {
      igraph::cluster_leiden(g, objective_function = "modularity", weights = w,
                             resolution = resolution, n_iterations = -1)
    } else {
      igraph::cluster_louvain(g, weights = w, resolution = resolution)
    })
    as.integer(igraph::membership(fit))
  }
  communities_from_assignments(net, stats::setNames(match(raw, unique(raw)), items))
}

#' A community result from a known partition
#'
#' `assignments` is an integer vector named by item, covering every vertex
#' of `net` (labels `1..K`). Returns the same structure as
#' [detect_communities()], with modularity recomputed from the graph -- for
#' a partition read back from disk rather than fitted.
#' @export
communities_from_assignments <- function(net, assignments) {
  stopifnot(inherits(net, "cooccurrence_network"))
  items <- net[["items"]]
  a <- as.integer(unname(assignments[items]))
  if (anyNA(a)) stop("assignments must name every vertex of the network", call. = FALSE)
  labels <- sort(unique(a))
  if (!identical(labels, seq_along(labels))) {
    stop("community labels must be 1..K", call. = FALSE)
  }
  communities <- lapply(labels, function(c) .sort_c(items[a == c]))
  names(communities) <- as.character(labels)
  list(assignments = stats::setNames(a, items),
       communities = communities,
       modularity = modularity_q(net[["adjacency"]], a),
       n_communities = length(labels))
}

#' Per-node network metrics
#'
#' One row per vertex, in `net$items` order: `item`, `degree`, `strength`
#' (weighted degree, 3 places), `betweenness` (normalized, on the
#' *unweighted* topology -- edge weights are association strengths, not
#' distances, so weighting would make a strong association a long path; 4
#' places), `clustering_coeff` (local, unweighted, 0 for degree < 2; 4
#' places), and `prevalence` (records holding the item).
#' @export
compute_network_metrics <- function(net) {
  stopifnot(inherits(net, "cooccurrence_network"))
  items <- net[["items"]]
  n <- length(items)
  if (n == 0L) {
    return(data.frame(item = character(0), degree = integer(0), strength = numeric(0),
                      betweenness = numeric(0), clustering_coeff = numeric(0),
                      prevalence = integer(0)))
  }
  adj <- net[["adjacency"]]
  topo <- adj != 0
  diag(topo) <- FALSE
  g <- igraph::graph_from_adjacency_matrix(topo, mode = "undirected", diag = FALSE)
  betweenness <- if (n <= 2L) rep(0, n) else
    as.numeric(igraph::betweenness(g, directed = FALSE, normalized = TRUE, weights = NULL))
  clustering <- as.numeric(igraph::transitivity(g, type = "localundirected",
                                                isolates = "zero", weights = NULL))
  data.frame(
    item = items,
    degree = as.integer(rowSums(topo)),
    strength = round_digits(rowSums(adj), 3L),
    betweenness = round_digits(betweenness, 4L),
    clustering_coeff = round_digits(clustering, 4L),
    prevalence = as.integer(unname(net[["prevalence"]][items]))
  )
}

#' Node table: metrics plus each node's community
#' @export
network_node_table <- function(net, communities) {
  metrics <- compute_network_metrics(net)
  metrics[["community"]] <- as.integer(unname(communities[["assignments"]][metrics[["item"]]]))
  if (anyNA(metrics[["community"]])) stop("a vertex has no community assignment", call. = FALSE)
  metrics
}

#' The `top_n` rows of `metrics` by a centrality column, descending
#'
#' A stable sort: tied nodes keep their `metrics` order.
#' @export
top_by_centrality <- function(metrics, centrality = "strength", top_n = 20L) {
  if (!(centrality %in% names(metrics))) {
    stop(sprintf("metrics has no '%s' column", centrality), call. = FALSE)
  }
  o <- order(metrics[[centrality]], decreasing = TRUE, method = "radix")
  out <- metrics[o[seq_len(min(top_n, nrow(metrics)))], , drop = FALSE]
  rownames(out) <- NULL
  out
}
