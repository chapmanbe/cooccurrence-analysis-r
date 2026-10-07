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

# Expected values below come from running the Julia reference on the same
# fixture (compute_pairwise_associations / build_cooccurrence_network /
# adjust_pvalues, CooccurrenceAnalysis.jl at 2026-10-06).

test_that("test_association with chisq (two-sided Pearson, no continuity correction)", {
  assoc <- compute_pairwise_associations(event_df, test = "chisq")
  expect_identical(paste(assoc[["item_a"]], assoc[["item_b"]], sep = "|"),
                   c("GA1|Item05", "GB1|Item01", "Item01|Item02", "Item03|Item04"))
  expect_equal(assoc[["p_value"]], c(0.025595697769205516, 0.035785546721876872,
                                     0.0085578871956978522, 0.00094237240515161465),
               tolerance = 1e-12)
  expect_equal(assoc[["p_adjusted"]], c(0.034127597025607352, 0.035785546721876872,
                                        0.017115774391395704, 0.0037694896206064586),
               tolerance = 1e-12)
  net <- build_cooccurrence_network(event_df, min_count = 1L, test = "chisq")
  expect_identical(network_edges(net)[, c("item_a", "item_b")],
                   data.frame(item_a = c("GA1", "GB1", "Item01", "Item03"),
                              item_b = c("Item05", "Item01", "Item02", "Item04")))
})

test_that("R-only: chisq agrees with stats::chisq.test(correct = FALSE)", {
  # Counts large enough that chisq.test() raises no approximation warning.
  ev <- events_from_records(list(list(items = c("A", "B"), n = 40L), list(items = "A", n = 20L),
                                 list(items = "B", n = 15L), list(items = c("C", "D"), n = 30L)))
  assoc <- compute_pairwise_associations(ev, test = "chisq")
  ab <- assoc[assoc[["item_a"]] == "A" & assoc[["item_b"]] == "B", ]
  ct <- matrix(c(40, 15, 20, 30), 2)
  expect_equal(ab[["p_value"]], stats::chisq.test(ct, correct = FALSE)$p.value, tolerance = 1e-12)
})

test_that("R-only: chisq raises on a zero margin (Julia returns NaN and poisons the family)", {
  # U is in every record. Julia gives p = NaN for A-U and B-U, and then
  # q = NaN for every pair, A-B included.
  ev <- data.frame(id = c(1, 1, 2, 2, 3, 3, 3, 4, 4),
                   item = c("U", "A", "U", "B", "U", "A", "B", "U", "A"))
  expect_error(compute_pairwise_associations(ev, test = "chisq"),
               "undefined for 2 pair\\(s\\) \\(first: A-U\\).*zero margin")
  expect_identical(nrow(compute_pairwise_associations(ev)), 3L)
})

test_that("network with Bonferroni and Holm corrections", {
  bonf <- compute_pairwise_associations(event_df, correction = "bonferroni")
  expect_equal(bonf[["p_adjusted"]], c(0.15889372766080384, 0.27142857142857157,
                                       0.04176181862857771, 0.091428571428571442),
               tolerance = 1e-12)
  holm <- compute_pairwise_associations(event_df, correction = "holm")
  expect_equal(holm[["p_adjusted"]], c(0.07944686383040192, 0.07944686383040192,
                                       0.04176181862857771, 0.068571428571428589),
               tolerance = 1e-12)
  for (corr in c("bonferroni", "holm")) {
    net <- build_cooccurrence_network(event_df, min_count = 1L, correction = corr)
    expect_identical(network_edges(net)[, c("item_a", "item_b")],
                     data.frame(item_a = "Item01", item_b = "Item02"))
  }
  piped <- run_network_pipeline(event_df, min_count = 1L, stratify_by_group = FALSE,
                                correction = "holm", seed = 1L)
  expect_identical(piped[["network"]][["items"]], c("Item01", "Item02"))
})

test_that("adjust_pvalues", {
  pvals <- c(0.001, 0.01, 0.03, 0.06, 0.5)
  adj <- bh_adjust(pvals)
  expect_identical(length(adj), 5L)
  expect_true(all(adj >= pvals))
  expect_identical(bh_adjust(numeric(0)), numeric(0))
  expect_equal(.adjust_p(pvals, "bonferroni"), c(0.005, 0.05, 0.15, 0.3, 1.0))
  expect_equal(.adjust_p(pvals, "holm"), c(0.005, 0.04, 0.09, 0.12, 0.5))
  expect_error(.adjust_p(c(0.1, NaN), "holm"), "contains NA")
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
                                 stratify_by_group = FALSE, seed = 1L)
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
