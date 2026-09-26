# Ported from the Julia suite (test/runtests.jl): one R case per Julia
# @testset where the behavior carries over, named after it. Cases the port
# adds beyond the Julia suite are marked "(R)".

event_df <- make_test_event_df()

#' A `cooccurrence_network` straight from a weight matrix, for graph-level
#' tests that do not need an event table.
net_from_adjacency <- function(adj, items = paste0("n", seq_len(nrow(adj)))) {
  dimnames(adj) <- list(items, items)
  structure(list(items = items, adjacency = adj, edge_data = data.frame(),
                 prevalence = stats::setNames(rep(1L, length(items)), items),
                 n_records = 1L, weight_metric = "lift", min_count = 1L, alpha = 1),
            class = "cooccurrence_network")
}

two_triangles <- function() {
  adj <- matrix(0, 6, 6)
  for (e in list(c(1, 2), c(2, 3), c(1, 3), c(4, 5), c(5, 6), c(4, 6))) {
    adj[e[[1]], e[[2]]] <- 1
    adj[e[[2]], e[[1]]] <- 1
  }
  adj[3, 4] <- 0.1
  adj[4, 3] <- 0.1
  adj
}

test_that("phi_coefficient", {
  expect_equal(phi_coefficient(matrix(c(10, 0, 0, 10), 2, byrow = TRUE)), 1)
  expect_equal(phi_coefficient(matrix(c(5, 5, 5, 5), 2, byrow = TRUE)), 0)
  expect_equal(phi_coefficient(matrix(c(0, 10, 10, 0), 2, byrow = TRUE)), -1)
  expect_identical(phi_coefficient(matrix(c(5, 0, 0, 0), 2, byrow = TRUE)), 0)
})

test_that("compute_pairwise_associations", {
  assoc <- compute_pairwise_associations(event_df)
  expect_gt(nrow(assoc), 0L)
  expect_true(all(c("item_a", "item_b", "lift", "phi", "p_adjusted") %in% names(assoc)))
  bt <- assoc[(assoc[["item_a"]] == "Item01" & assoc[["item_b"]] == "Item02") |
                (assoc[["item_a"]] == "Item02" & assoc[["item_b"]] == "Item01"), ]
  expect_identical(nrow(bt), 1L)
  expect_gt(bt[["lift"]], 1)
})

test_that("compute_pairwise_associations: columns, order, and rounding (R)", {
  assoc <- compute_pairwise_associations(event_df)
  expect_identical(names(assoc), c("item_a", "item_b", "observed", "expected", "lift", "phi",
                                   "odds_ratio", "p_value", "n_a", "n_b", "p_adjusted"))
  expect_identical(order(assoc[["item_a"]], assoc[["item_b"]], method = "radix"),
                   seq_len(nrow(assoc)))
  bt <- assoc[assoc[["item_a"]] == "Item01" & assoc[["item_b"]] == "Item02", ]
  # Item01 in 8 + 10 + 3 = 21 records, Item02 in 4 + 10 = 14, of 50.
  expect_identical(c(bt[["observed"]], bt[["n_a"]], bt[["n_b"]]), c(10L, 21L, 14L))
  expect_identical(bt[["expected"]], round_digits(21 * 14 / 50, 2L))
  expect_identical(bt[["lift"]], round_digits(10 / (21 * 14 / 50), 4L))
  expect_identical(bt[["odds_ratio"]], round_digits(10 * 25 / (11 * 4), 4L))
  expect_identical(bt[["phi"]],
                   round_digits(phi_coefficient(matrix(c(10, 11, 4, 25), 2, byrow = TRUE)), 4L))
})

test_that("build_cooccurrence_network", {
  net <- build_cooccurrence_network(event_df, min_count = 2L, alpha = 0.5)
  expect_s3_class(net, "cooccurrence_network")
  expect_identical(net[["n_records"]], 50L)
  expect_identical(nrow(net[["adjacency"]]), length(net[["items"]]))
  expect_gt(network_size(net)[["n_edges"]], 0L)
  expect_gt(length(net[["prevalence"]]), 0L)
  strict <- build_cooccurrence_network(event_df, min_count = 20L, alpha = 0.01)
  expect_lte(network_size(strict)[["n_edges"]], network_size(net)[["n_edges"]])
})

test_that("build_cooccurrence_network: vertices, weights, prevalence (R)", {
  net <- build_cooccurrence_network(event_df, weight_metric = "phi", min_count = 2L, alpha = 0.5)
  e <- network_edges(net)
  # Vertices are exactly the items in some edge, sorted.
  expect_identical(net[["items"]], sort(unique(c(e[["item_a"]], e[["item_b"]])), method = "radix"))
  ed <- net[["edge_data"]]
  m <- merge(e, ed, by = c("item_a", "item_b"))
  expect_identical(m[["weight"]], pmax(m[["phi"]], 0.001))
  expect_identical(unname(net[["prevalence"]][["Item01"]]), 21L)
  lift_net <- build_cooccurrence_network(event_df, weight_metric = "lift", min_count = 2L,
                                         alpha = 0.5)
  m <- merge(network_edges(lift_net), ed, by = c("item_a", "item_b"))
  expect_identical(m[["weight"]], m[["lift"]])
})

test_that("cooccurrence_network_from_edge_data rebuilds the same network (R)", {
  net <- build_cooccurrence_network(event_df, weight_metric = "phi", min_count = 2L, alpha = 0.5)
  again <- cooccurrence_network_from_edge_data(net[["edge_data"]], net[["prevalence"]],
                                               net[["n_records"]], "phi", 2L, 0.5)
  expect_identical(again, net)
  expect_error(cooccurrence_network_from_edge_data(net[["edge_data"]], integer(0), 50L,
                                                   "phi", 2L, 0.5), "prevalence")
})

test_that("detect_communities: louvain", {
  net <- build_cooccurrence_network(event_df, min_count = 2L, alpha = 0.5)
  comm <- detect_communities(net, method = "louvain")
  expect_identical(length(comm[["assignments"]]), length(net[["items"]]))
  expect_gte(comm[["n_communities"]], 1L)
  expect_true(comm[["modularity"]] >= 0 || comm[["n_communities"]] == 1L)
  expect_setequal(unlist(comm[["communities"]], use.names = FALSE), net[["items"]])
})

test_that("detect_communities: louvain modularity gain matches delta Q (C1, C2, T6)", {
  adj <- two_triangles()
  strengths <- rowSums(adj)
  m2 <- sum(strengths)
  comm_strength <- function(assign) vapply(1:6, function(c) sum(strengths[assign == c]), 0)
  assign <- 1:6
  cs <- comm_strength(assign)
  for (move in list(c(2L, 1L), c(5L, 4L), c(4L, 3L))) {
    node <- move[[1L]]
    target <- move[[2L]]
    gain <- CooccurrenceAnalysis:::.modularity_gain(adj, node, target, assign, strengths, cs, m2)
    loss <- CooccurrenceAnalysis:::.modularity_gain(adj, node, assign[[node]], assign,
                                                    strengths, cs, m2)
    moved <- assign
    moved[[node]] <- target
    expect_equal(gain - loss, modularity_q(adj, moved) - modularity_q(adj, assign),
                 tolerance = 1e-10)
  }
  final <- CooccurrenceAnalysis:::.louvain(adj, 100L)
  expect_true(final[[1]] == final[[2]] && final[[2]] == final[[3]])
  expect_true(final[[4]] == final[[5]] && final[[5]] == final[[6]])
  expect_true(final[[1]] != final[[4]])
  expect_gt(modularity_q(adj, final), modularity_q(adj, rep(1L, 6)))
})

test_that("detect_communities: renumbering, names, and unsupported methods (R)", {
  comm <- detect_communities(net_from_adjacency(two_triangles()))
  expect_identical(unname(comm[["assignments"]]), c(1L, 1L, 1L, 2L, 2L, 2L))
  expect_identical(names(comm[["assignments"]]), paste0("n", 1:6))
  expect_identical(comm[["communities"]],
                   list(`1` = c("n1", "n2", "n3"), `2` = c("n4", "n5", "n6")))
  # A one-community partition scores exactly zero (the diagonal/null terms).
  expect_equal(modularity_q(two_triangles(), rep(1L, 6)), 0)
  expect_error(detect_communities(net_from_adjacency(two_triangles()),
                                  method = "label_propagation"), "method")
})

test_that("communities_from_assignments reproduces a fitted partition (R)", {
  net <- net_from_adjacency(two_triangles())
  comm <- detect_communities(net)
  expect_identical(communities_from_assignments(net, comm[["assignments"]]), comm)
  expect_error(communities_from_assignments(net, comm[["assignments"]][-1]), "every vertex")
  bad <- comm[["assignments"]] + 1L
  expect_error(communities_from_assignments(net, bad), "1..K")
})

test_that("compute_network_metrics", {
  net <- build_cooccurrence_network(event_df, min_count = 2L, alpha = 0.5)
  metrics <- compute_network_metrics(net)
  expect_identical(nrow(metrics), length(net[["items"]]))
  expect_true(all(c("item", "degree", "strength", "betweenness", "clustering_coeff",
                    "prevalence") %in% names(metrics)))
  expect_true(all(metrics[["degree"]] >= 0L))
})

test_that("compute_network_metrics: known values on two bridged triangles (R)", {
  m <- compute_network_metrics(net_from_adjacency(two_triangles()))
  expect_identical(m[["degree"]], c(2L, 2L, 3L, 3L, 2L, 2L))
  expect_identical(m[["strength"]], c(2, 2, 2.1, 2.1, 2, 2))
  # Unweighted betweenness, normalized by (n-1)(n-2)/2 = 10: each bridge end
  # is interior to the 6 cross-triangle shortest paths that do not end at
  # it, so 0.6. Weighting by the 0.1 bridge must not change this.
  expect_identical(m[["betweenness"]], c(0, 0, 0.6, 0.6, 0, 0))
  expect_identical(m[["clustering_coeff"]], c(1, 1, round_digits(1 / 3, 4L),
                                              round_digits(1 / 3, 4L), 1, 1))
  path <- matrix(0, 3, 3)
  path[1, 2] <- path[2, 1] <- path[2, 3] <- path[3, 2] <- 1
  mp <- compute_network_metrics(net_from_adjacency(path))
  expect_identical(mp[["clustering_coeff"]], c(0, 0, 0))
  expect_identical(mp[["betweenness"]], c(0, 1, 0))
})

test_that("top_by_centrality is a stable descending sort (R)", {
  metrics <- data.frame(item = c("a", "b", "c", "d"), betweenness = c(0, 0.5, 0, 0.5))
  top <- top_by_centrality(metrics, "betweenness", 3L)
  expect_identical(top[["item"]], c("b", "d", "a"))
  expect_error(top_by_centrality(metrics, "strength"), "strength")
})

test_that("round_digits matches the Julia reference on ties (R)", {
  expect_identical(round_digits(0.15, 1L), 0.2)
  expect_identical(round_digits(1.015, 2L), 1.01)
  expect_identical(round_digits(c(Inf, NaN, 2.5), 0L), c(Inf, NaN, 2))
})

test_that("network_visualization: plot_cooccurrence_network", {
  net <- build_cooccurrence_network(event_df, min_count = 2L, alpha = 0.5)
  comm <- detect_communities(net)
  expect_s3_class(plot_cooccurrence_network(net, communities = comm), "ggplot")
})

test_that("network_visualization: plot_community_heatmap", {
  net <- build_cooccurrence_network(event_df, min_count = 2L, alpha = 0.5)
  comm <- detect_communities(net)
  expect_s3_class(plot_community_heatmap(net, comm), "ggplot")
})

test_that("extended_network_visualization: plot_centrality_barchart", {
  net <- build_cooccurrence_network(event_df, min_count = 1L, alpha = 1)
  metrics <- compute_network_metrics(net)
  expect_s3_class(plot_centrality_barchart(metrics), "ggplot")
})
