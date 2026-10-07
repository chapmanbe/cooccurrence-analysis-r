# Independent cross-checks against validated R packages. These do not replace
# the Julia-derived expected values elsewhere in the suite; they add a second,
# R-side reference so a bug shared by the Julia code and its port would still
# surface. Each check skips when the other package is absent.

# Random records x items presence matrix in which every record holds at least
# two items (so no record is excluded by `min_items`) and every item appears.
crosscheck_matrix <- function(n_records = 60L, n_items = 8L, p = 0.35, seed = 20260506L) {
  set.seed(seed)
  repeat {
    m <- matrix(rbinom(n_items * n_records, 1L, p), nrow = n_items,
                dimnames = list(sprintf("it%02d", seq_len(n_items)), NULL))
    if (all(colSums(m) >= 2L) && all(rowSums(m) >= 1L)) return(m)
  }
}

events_from_matrix <- function(m) {
  hits <- which(m == 1L, arr.ind = TRUE)
  data.frame(id = hits[, 2L], item = rownames(m)[hits[, 1L]])
}

pair_key <- function(a, b) paste(pmin(a, b), pmax(a, b), sep = "|")

test_that("pairwise statistics agree with cooccur (Veech 2013) on shared quantities", {
  skip_if_not_installed("cooccur")
  m <- crosscheck_matrix()
  ours <- compute_pairwise_associations(events_from_matrix(m))
  # thresh = FALSE keeps pairs with expected co-occurrence < 1, which the
  # package under test also keeps.
  theirs <- cooccur::cooccur(m, type = "spp_site", thresh = FALSE, spp_names = TRUE)$results
  theirs$key <- pair_key(theirs$sp1_name, theirs$sp2_name)
  ours$key <- pair_key(ours$item_a, ours$item_b)

  # Pairs that never co-occur are outside this package's family by design;
  # every pair that does co-occur must appear in cooccur's table.
  expect_true(all(ours$key %in% theirs$key))
  j <- theirs[match(ours$key, theirs$key), ]

  expect_equal(ours$observed, j$obs_cooccur)
  # Per-item record counts, matched by item name rather than pair position.
  inc <- ifelse(ours$item_a == j$sp1_name, j$sp1_inc, j$sp2_inc)
  expect_equal(ours$n_a, inc)
  # `expected` is rounded to 2 places here; cooccur reports it at 1 place.
  expect_equal(ours$expected, j$exp_cooccur, tolerance = 0.06)
  # One-sided upper-tail Fisher p equals cooccur's p_gt, the probability of
  # an observed co-occurrence at least this large. cooccur rounds to 5 places,
  # so compare on the absolute scale.
  expect_lt(max(abs(ours$p_value - j$p_gt)), 6e-6)
})

test_that("pairwise p-values agree with stats::fisher.test directly", {
  m <- crosscheck_matrix(seed = 7L)
  ours <- compute_pairwise_associations(events_from_matrix(m))
  n <- ncol(m)
  ref <- vapply(seq_len(nrow(ours)), function(i) {
    a <- m[ours$item_a[[i]], ]
    b <- m[ours$item_b[[i]], ]
    stats::fisher.test(table(factor(a, 0:1), factor(b, 0:1)), alternative = "greater")$p.value
  }, numeric(1L))
  expect_equal(ours$p_value, ref, tolerance = 1e-8)
  expect_equal(ours$p_adjusted, stats::p.adjust(ref, method = "BH"), tolerance = 1e-8)
})

# A planted-partition graph: three dense blocks joined by weak noise edges.
planted_net <- function(seed = 11L, block = 8L, k = 3L) {
  set.seed(seed)
  n <- block * k
  truth <- rep(seq_len(k), each = block)
  w <- matrix(0, n, n)
  for (i in seq_len(n - 1L)) {
    for (j in (i + 1L):n) {
      p <- if (truth[[i]] == truth[[j]]) 0.8 else 0.05
      if (runif(1L) < p) w[i, j] <- w[j, i] <- runif(1L, 0.5, 2)
    }
  }
  items <- sprintf("n%02d", seq_len(n))
  dimnames(w) <- list(items, items)
  list(net = structure(list(items = items, adjacency = w), class = "cooccurrence_network"),
       truth = truth)
}

test_that("modularity_q matches igraph::modularity on weighted graphs", {
  skip_if_not_installed("igraph")
  pn <- planted_net()
  g <- igraph::graph_from_adjacency_matrix(pn$net$adjacency, mode = "undirected",
                                           weighted = TRUE, diag = FALSE)
  set.seed(3L)
  n <- length(pn$truth)
  for (labels in list(pn$truth, sample(1:4, n, replace = TRUE), rep(1L, n))) {
    expect_equal(modularity_q(pn$net$adjacency, labels),
                 igraph::modularity(g, membership = labels, weights = igraph::E(g)$weight),
                 tolerance = 1e-10)
  }
})

test_that("detect_communities recovers planted blocks and is competitive with igraph Louvain", {
  skip_if_not_installed("igraph")
  pn <- planted_net()
  ours <- detect_communities(pn$net)
  g <- igraph::graph_from_adjacency_matrix(pn$net$adjacency, mode = "undirected",
                                           weighted = TRUE, diag = FALSE)
  ref <- igraph::cluster_louvain(g)

  expect_gt(adjusted_rand(unname(ours$assignments), pn$truth), 0.9)
  expect_gt(adjusted_rand(as.integer(igraph::membership(ref)), pn$truth), 0.9)
  # Same objective, different tie-breaking: Q must be within a small margin of
  # igraph's, and the reported Q must be the true Q of the returned partition.
  expect_gte(ours$modularity, igraph::modularity(ref) - 0.02)
  expect_equal(ours$modularity,
               igraph::modularity(g, membership = unname(ours$assignments),
                                  weights = igraph::E(g)$weight),
               tolerance = 1e-10)
})
