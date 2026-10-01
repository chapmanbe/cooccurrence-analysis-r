# Ported from the Julia suite's "statistical_validation", "network
# timing_filter", "detect_communities / label_propagation", and
# "network_integration" testsets. Where the R port does not carry a feature
# (it is outside what the analysis pipeline reaches; see PORTING.md, step 3
# notes), the paired R case asserts that asking for it raises an informative
# error rather than silently falling back.

event_df <- make_test_event_df()

test_that("odds_ratio zero-cell handling (C5)", {
  # Expected values computed by running the Julia reference
  # (CooccurrenceAnalysis.odds_ratio) on the same tables.
  tbl <- function(a, b, c, d) matrix(c(a, b, c, d), 2L, byrow = TRUE)
  expect_equal(odds_ratio(tbl(10, 2, 3, 20)), 33.333333333333336, tolerance = 1e-15)
  expect_equal(odds_ratio(tbl(10, 2, 3, 20), correction = "none"), 33.333333333333336,
               tolerance = 1e-15)
  zc <- tbl(5, 0, 0, 5)
  expect_true(is.finite(odds_ratio(zc)))
  expect_identical(odds_ratio(zc), 121)
  expect_identical(odds_ratio(zc, correction = "none"), Inf)
  expect_true(is.finite(odds_ratio(tbl(0, 0, 3, 4))))
  expect_equal(odds_ratio(tbl(0, 0, 3, 4)), 1.2857142857142858, tolerance = 1e-15)
  expect_true(is.nan(odds_ratio(tbl(0, 0, 3, 4), correction = "none")))
  expect_error(odds_ratio(tbl(1, 1, 1, 1), correction = "bogus"), "correction")
  # (R) The zero-cell shapes a network table can hold (n11 >= 1 always):
  # an empty n10, n01, or n00.
  expect_equal(odds_ratio(tbl(7, 0, 4, 30)), 101.66666666666667, tolerance = 1e-15)
  expect_equal(odds_ratio(tbl(9, 3, 0, 12)), 67.85714285714286, tolerance = 1e-15)
  expect_equal(odds_ratio(tbl(4, 6, 5, 0)), 0.06293706293706294, tolerance = 1e-15)
})

test_that("pairwise odds_ratio is Haldane-corrected, never Inf (R)", {
  # Item01 + GB1: every GB1 record also holds Item01, so n01 = 0.
  assoc <- compute_pairwise_associations(event_df)
  row <- assoc[assoc[["item_a"]] == "GB1" & assoc[["item_b"]] == "Item01", ]
  expect_identical(nrow(row), 1L)
  # GB1 is a, Item01 is b: n11 = 3, n10 = 0, n01 = 18, n00 = 29.
  expect_identical(row[["odds_ratio"]], round_digits((3.5 * 29.5) / (0.5 * 18.5), 4L))
  expect_true(all(is.finite(assoc[["odds_ratio"]])))
})

test_that("test_association (Fisher: one-sided, significant for the planted pair)", {
  assoc <- compute_pairwise_associations(event_df)
  bt <- assoc[assoc[["item_a"]] == "Item01" & assoc[["item_b"]] == "Item02", ]
  expect_identical(bt[["observed"]], 10L)
  expect_gt(bt[["expected"]], 0)
  expect_gt(bt[["odds_ratio"]], 1)
  expect_lt(bt[["p_value"]], 0.05)
})

test_that("test_association with chisq: not ported, raises", {
  expect_error(compute_pairwise_associations(event_df, test = "chisq"),
               "test = 'chisq' exists in the Julia reference but is not ported")
  expect_error(build_cooccurrence_network(event_df, test = "chisq"), "not ported")
})

test_that("adjust_pvalues", {
  pvals <- c(0.001, 0.01, 0.03, 0.06, 0.5)
  adj <- bh_adjust(pvals)
  expect_identical(length(adj), 5L)
  expect_true(all(adj >= pvals))
  expect_identical(bh_adjust(numeric(0)), numeric(0))
})

test_that("adjust_pvalues on empty input: no co-occurring pairs gives an empty table (R)", {
  singles <- data.frame(id = 1:4, item = c("A", "B", "C", "A"))
  assoc <- compute_pairwise_associations(singles)
  expect_identical(nrow(assoc), 0L)
  expect_identical(names(assoc), c("item_a", "item_b", "observed", "expected", "lift", "phi",
                                   "odds_ratio", "p_value", "n_a", "n_b", "p_adjusted"))
  net <- build_cooccurrence_network(singles, min_count = 1L, alpha = 1)
  expect_identical(network_size(net)[["n_nodes"]], 0L)
})

test_that("network timing_filter: filter_event_by_timing helper -- not ported, raises", {
  for (tf in c("concurrent", "sequential")) {
    expect_error(compute_pairwise_associations(event_df, timing_filter = tf),
                 sprintf("timing_filter = '%s'.*not ported", tf))
  }
  expect_identical(compute_pairwise_associations(event_df, timing_filter = "all"),
                   compute_pairwise_associations(event_df))
})

test_that("network timing_filter: compute_pairwise_associations + timing_filter -- raises", {
  expect_error(compute_pairwise_associations(event_df, timing_filter = "concurrent"),
               "not ported")
  expect_error(compute_pairwise_associations(event_df, timing_filter = "sequential"),
               "not ported")
})

test_that("network timing_filter: build_cooccurrence_network + timing_filter -- raises", {
  expect_error(build_cooccurrence_network(event_df, min_count = 2L, alpha = 0.5,
                                          timing_filter = "concurrent"), "not ported")
  expect_error(build_cooccurrence_network(event_df, min_count = 2L, alpha = 0.5,
                                          timing_filter = "sequential"), "not ported")
})

test_that("network timing_filter: rejects bad timing_filter", {
  expect_error(compute_pairwise_associations(event_df, timing_filter = "bogus"),
               "unsupported timing_filter")
})

test_that("detect_communities: label_propagation -- not ported, raises", {
  net <- build_cooccurrence_network(event_df, min_count = 2L, alpha = 0.5)
  expect_error(detect_communities(net, method = "label_propagation"),
               "label_propagation' exists in the Julia reference but is not ported")
})

test_that("network_integration", {
  result <- run_network_pipeline(event_df, min_count = 2L, alpha = 0.5,
                                 stratify_by_group = FALSE)
  expect_s3_class(result[["network"]], "cooccurrence_network")
  expect_true(all(c("assignments", "modularity", "n_communities") %in%
                    names(result[["communities"]])))
  expect_gt(nrow(result[["metrics"]]), 0L)
  # The Julia case runs stratify_by_group = true -- also the default in both
  # languages; it needs compare_networks (not ported), so it raises either way.
  expect_error(run_network_pipeline(event_df, min_count = 2L, alpha = 0.5,
                                    stratify_by_group = TRUE), "compare_networks")
  expect_error(run_network_pipeline(event_df, min_count = 2L, alpha = 0.5), "compare_networks")
})

test_that("Louvain counts exact-gain ties (R)", {
  # A symmetric triangle: node 1's two candidate communities {2} and {3}
  # give exactly the same gain, which is where the port and the Julia
  # reference may choose differently.
  tri <- matrix(1, 3, 3) - diag(3)
  fit <- CooccurrenceAnalysis:::.louvain(tri, 100L)
  expect_gt(fit[["n_gain_ties"]], 0L)
  # A 4-cycle with equal weights: node 1's candidates {2} and {4} tie
  # exactly; R takes the smaller label.
  cyc <- matrix(0, 4, 4)
  for (e in list(c(1, 2), c(2, 3), c(3, 4), c(4, 1))) {
    cyc[e[[1]], e[[2]]] <- cyc[e[[2]], e[[1]]] <- 1
  }
  fit <- CooccurrenceAnalysis:::.louvain(cyc, 100L)
  expect_gt(fit[["n_gain_ties"]], 0L)
  expect_identical(fit[["assignments"]][[1]], fit[["assignments"]][[2]])
  # A star whose center, visited first, has three candidates with distinct
  # gains: several candidates, but no exact tie.
  star <- matrix(0, 4, 4)
  star[1, 2] <- star[2, 1] <- 1
  star[1, 3] <- star[3, 1] <- 0.5
  star[1, 4] <- star[4, 1] <- 0.25
  expect_identical(CooccurrenceAnalysis:::.louvain(star, 100L)[["n_gain_ties"]], 0L)
  net <- build_cooccurrence_network(event_df, min_count = 2L, alpha = 0.5)
  expect_true(is.integer(detect_communities(net)[["n_gain_ties"]]))
})
