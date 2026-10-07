# Community detection via igraph (Leiden default, Louvain option). These
# replace the Julia-parity tests of the old phase-1-only Louvain, a
# documented divergence (CLAUDE.md, "Divergences from the Julia reference").
# Expected partitions come from fixtures that are unambiguous by
# construction, not from the Julia reference.

test_that("R-only: .with_seed seeds locally and restores the caller's stream", {
  set.seed(42)
  expected_next <- stats::runif(1)
  set.seed(42)
  a <- .with_seed(1L, stats::runif(3))
  expect_identical(stats::runif(1), expected_next)
  expect_identical(.with_seed(1L, stats::runif(3)), a)
})

test_that("R-only: .with_seed(NULL) draws from the caller's stream", {
  set.seed(3)
  a <- .with_seed(NULL, stats::runif(2))
  set.seed(3)
  expect_identical(stats::runif(2), a)
})

test_that("R-only: .with_seed leaves no .Random.seed behind in a fresh session", {
  key <- ".Random.seed"
  had <- exists(key, envir = globalenv(), inherits = FALSE)
  if (had) saved <- get(key, envir = globalenv(), inherits = FALSE)
  on.exit(if (had) assign(key, saved, envir = globalenv()))
  if (had) rm(list = key, envir = globalenv())
  .with_seed(1L, stats::runif(1))
  expect_false(exists(key, envir = globalenv(), inherits = FALSE))
})

test_that("R-only: .with_seed depends on the seed alone, not the caller's RNG kind", {
  default <- .with_seed(1L, stats::runif(3))
  old <- RNGkind("L'Ecuyer-CMRG")
  on.exit(RNGkind(old[[1L]], old[[2L]], old[[3L]]))
  expect_identical(.with_seed(1L, stats::runif(3)), default)
  # The caller's own generator is still in force afterward.
  expect_identical(RNGkind()[[1L]], "L'Ecuyer-CMRG")
})

test_that("R-only: .check_seed accepts NULL and whole numbers, raises otherwise", {
  expect_silent(.check_seed(NULL, "f"))
  expect_silent(.check_seed(7, "f"))
  expect_silent(.check_seed(-5L, "f"))
  for (bad in list(1.5, NA, NA_integer_, c(1, 2), "1", 2^31, Inf)) {
    expect_error(.check_seed(bad, "f"), "f: seed must be NULL or one whole number")
  }
})

canon <- function(m) match(m, unique(m))

test_that("R-only: detect_communities finds two bridged triangles under any seed", {
  net <- net_from_adjacency(two_triangles())
  for (method in c("leiden", "louvain")) {
    for (s in 1:5) {
      comm <- detect_communities(net, method = method, seed = s)
      expect_identical(unname(comm[["assignments"]]), c(1L, 1L, 1L, 2L, 2L, 2L))
      expect_identical(comm[["communities"]],
                       list(`1` = c("n1", "n2", "n3"), `2` = c("n4", "n5", "n6")))
      expect_identical(names(comm), c("assignments", "communities", "modularity",
                                      "n_communities"))
      expect_equal(comm[["modularity"]], modularity_q(two_triangles(), c(1, 1, 1, 2, 2, 2)))
    }
  }
})

test_that("R-only: disconnected components are separate communities", {
  adj <- two_triangles()
  adj[3, 4] <- adj[4, 3] <- 0
  for (s in 1:5) {
    expect_identical(unname(detect_communities(net_from_adjacency(adj), seed = s)[["assignments"]]),
                     c(1L, 1L, 1L, 2L, 2L, 2L))
  }
})

test_that("R-only: full Leiden/Louvain beat the reference's phase-1 local moving", {
  # 20 triangles in a ring. Phase-1-only local moving (the Julia reference,
  # and this package before the change) stops at Q = 0.70 with one community
  # per triangle; aggregation reached Q = 0.770-0.775 under every seed in 1..20. Julia's
  # `_louvain` was run on this graph: Q = 0.70, 20 communities, any max_iter.
  net <- net_from_adjacency(ring_of_cliques(20L, 3L))
  for (method in c("leiden", "louvain")) {
    for (s in 1:5) {
      comm <- detect_communities(net, method = method, seed = s)
      expect_gt(comm[["modularity"]], 0.75)
      expect_lt(comm[["n_communities"]], 20L)
    }
  }
})

test_that("R-only: labels are numbered by first appearance in item order", {
  net <- net_from_adjacency(ring_of_cliques(20L, 3L))
  for (s in 1:10) {
    a <- unname(detect_communities(net, seed = s)[["assignments"]])
    expect_identical(unique(a), seq_len(max(a)))
  }
})

test_that("R-only: same seed, same partition; seed = NULL follows the caller's stream", {
  net <- net_from_adjacency(ring_of_cliques(20L, 3L))
  expect_identical(detect_communities(net, seed = 11L), detect_communities(net, seed = 11L))
  set.seed(5)
  a <- detect_communities(net, seed = NULL)
  set.seed(5)
  expect_identical(detect_communities(net, seed = NULL), a)
})

test_that("R-only: seed is required, and bad methods and resolutions raise", {
  net <- net_from_adjacency(two_triangles())
  expect_error(detect_communities(net), "seed is required")
  expect_error(detect_communities(net, seed = 1.5), "seed must be NULL or one whole number")
  expect_error(detect_communities(net, method = "walktrap", seed = 1L), "unsupported method")
  expect_error(detect_communities(net, method = "label_propagation", seed = 1L),
               "label_propagation' exists in the Julia reference but is not ported")
  for (bad in list(0, -1, NA_real_, c(1, 2), "1")) {
    expect_error(detect_communities(net, seed = 1L, resolution = bad), "resolution")
  }
})

test_that("R-only: edgeless and empty networks need no igraph call", {
  # igraph's cluster_leiden(n_iterations = -1) errors on a graph with no edges.
  for (method in c("leiden", "louvain")) {
    comm <- detect_communities(net_from_adjacency(matrix(0, 3, 3)), method = method, seed = 1L)
    expect_identical(unname(comm[["assignments"]]), 1:3)
    expect_identical(comm[["modularity"]], 0)
  }
  empty <- net_from_adjacency(matrix(0, 0, 0), character(0))
  expect_identical(detect_communities(empty, seed = 1L)[["n_communities"]], 0L)
  expect_identical(detect_communities(empty, seed = NULL)[["communities"]], list())
})

test_that("R-only: an isolated vertex beside edges gets its own community", {
  adj <- matrix(0, 3, 3)
  adj[1, 2] <- adj[2, 1] <- 1
  comm <- detect_communities(net_from_adjacency(adj), seed = 1L)
  expect_identical(unname(comm[["assignments"]]), c(1L, 1L, 2L))
})

test_that("R-only: stratified analysis and the pipeline pass seed and method through", {
  ev <- make_test_event_df()
  res <- stratified_network_analysis(ev, exclusive_items = NULL, min_count = 1L, alpha = 1,
                                     community_method = "louvain", seed = 3L)
  for (g in names(res)) {
    expect_identical(res[[g]][["communities"]],
                     detect_communities(res[[g]][["network"]], method = "louvain", seed = 3L))
  }
  expect_error(stratified_network_analysis(ev, exclusive_items = NULL), "seed is required")
  piped <- run_network_pipeline(ev, min_count = 1L, alpha = 1, stratify_by_group = FALSE,
                                seed = 3L)
  expect_identical(piped[["communities"]], detect_communities(piped[["network"]], seed = 3L))
  expect_error(run_network_pipeline(ev, stratify_by_group = FALSE), "seed is required")
})
