# Network metrics and community detection.

#' Modularity Q of a partition of a weighted graph
#'
#' `Q = sum_c [L_c / m - (D_c / 2m)^2]`, where `m` is the total edge weight,
#' `L_c` the weight of edges inside community `c` (each counted once), and
#' `D_c` the total strength of its nodes. The `(D_c / 2m)^2` term carries the
#' null-model penalty for every same-community pair, edge or not, and the
#' diagonal self-terms, so a one-community partition of any graph scores 0.
#' Returns 0 for an empty or edgeless graph.
#'
#' `adjacency` is a symmetric weight matrix; `assignments` one community
#' label per row (any integers).
#' @export
modularity_q <- function(adjacency, assignments) {
  n <- nrow(adjacency)
  stopifnot(length(assignments) == n)
  if (n == 0L) return(0)
  m2 <- sum(adjacency[upper.tri(adjacency)]) * 2
  if (m2 == 0) return(0)
  strengths <- rowSums(adjacency)
  q <- 0
  for (c in unique(assignments)) {
    members <- assignments == c
    d_c <- sum(strengths[members])
    inner <- adjacency[members, members, drop = FALSE]
    l_c <- sum(inner[upper.tri(inner)])
    q <- q + 2 * l_c / m2 - (d_c / m2)^2
  }
  q
}

#' Neighbors of vertex `v` (ascending), excluding `v` itself
#' @keywords internal
.neighbors <- function(adjacency, v) {
  nb <- which(adjacency[v, ] != 0)
  nb[nb != v]
}

#' Modularity gain from adding `node` to `community`
#'
#' `2 k_in / 2m - 2 s_tot k_i / (2m)^2`, with `k_in` the edge weight from
#' `node` into `community`, `s_tot` the community's total strength excluding
#' `node` itself, and `k_i` the node's strength; `m2` is `2m`. A move
#' `current -> target` gains `gain(target) - gain(current)`, which equals the
#' change in [modularity_q()] exactly (the `k_i^2` self-terms cancel).
#' `comm_strength[c]` is the maintained total strength of community `c`.
#' @keywords internal
.modularity_gain <- function(adjacency, node, community, assignments, strengths,
                             comm_strength, m2) {
  nb <- .neighbors(adjacency, node)
  in_comm <- nb[assignments[nb] == community]
  k_in <- sum(adjacency[node, in_comm])
  s_tot <- comm_strength[[community]]
  if (assignments[[node]] == community) s_tot <- s_tot - strengths[[node]]
  2 * k_in / m2 - 2 * (s_tot * strengths[[node]]) / (m2 * m2)
}

#' Louvain community detection, phase 1 (local moving)
#'
#' Every node starts in its own community. Nodes are visited in order
#' `1..n`; each moves to the neighboring community with the largest positive
#' modularity gain `2 k_in / 2m - 2 s_tot k_i / (2m)^2` (relative to staying,
#' with the node's own strength excluded from its current community). Sweeps
#' repeat until a full sweep moves nothing, or `max_iter` sweeps. There is no
#' aggregation phase: the graphs here are small enough that local moving is
#' the whole algorithm.
#'
#' Candidate communities are considered in ascending label order and a move
#' requires a strictly larger gain, so exact ties keep the earlier label.
#' Returns one label per node (labels are node indices, not renumbered).
#' @keywords internal
.louvain <- function(adjacency, max_iter) {
  n <- nrow(adjacency)
  assignments <- seq_len(n)
  m2 <- sum(adjacency[upper.tri(adjacency)]) * 2
  if (m2 == 0) return(assignments)
  strengths <- rowSums(adjacency)
  comm_strength <- strengths
  neighbors <- lapply(seq_len(n), function(v) .neighbors(adjacency, v))

  for (iter in seq_len(max_iter)) {
    improved <- FALSE
    for (node in seq_len(n)) {
      current <- assignments[[node]]
      candidates <- sort(unique(assignments[neighbors[[node]]]))
      loss <- .modularity_gain(adjacency, node, current, assignments, strengths,
                               comm_strength, m2)
      best <- current
      best_gain <- 0
      for (comm in candidates) {
        if (comm == current) next
        net_gain <- .modularity_gain(adjacency, node, comm, assignments, strengths,
                                     comm_strength, m2) - loss
        if (net_gain > best_gain) {
          best_gain <- net_gain
          best <- comm
        }
      }
      if (best != current) {
        comm_strength[[current]] <- comm_strength[[current]] - strengths[[node]]
        comm_strength[[best]] <- comm_strength[[best]] + strengths[[node]]
        assignments[[node]] <- best
        improved <- TRUE
      }
    }
    if (!improved) break
  }
  assignments
}

#' Detect communities in a co-occurrence network
#'
#' `method = "louvain"` is the only method; anything else raises. Community
#' labels are renumbered `1..K` in ascending order of their raw label.
#'
#' Returns a list: `assignments` (integer, named by item, in `net$items`
#' order), `communities` (a list named by community label, `"1"`, `"2"`, ...,
#' each the sorted member items), `modularity` ([modularity_q()] of the
#' partition), and `n_communities`.
#' @export
detect_communities <- function(net, method = "louvain", max_iter = 100L) {
  stopifnot(inherits(net, "cooccurrence_network"))
  .check_choice(method, "louvain", "method")
  items <- net[["items"]]
  if (length(items) == 0L) {
    return(list(assignments = stats::setNames(integer(0), character(0)),
                communities = list(), modularity = 0, n_communities = 0L))
  }
  raw <- .louvain(net[["adjacency"]], max_iter)
  labels <- sort(unique(raw))
  communities_from_assignments(net, stats::setNames(match(raw, labels), items))
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
