# Specification-first tests for the significance family (port step 3).
#
# Written from the project's specification -- "significance testing is
# Fisher's exact with Benjamini-Hochberg correction"; an edge needs
# co-occurrence >= min_count, BH-adjusted p < alpha, and lift > 1 -- before
# any implementation existed, and without reading the Julia implementation
# of the test or the correction. Two properties the specification fixes that
# a plausible implementation could silently get wrong:
#
# 1. The BH family is the set of pairs actually tested: pairs that co-occur
#    in at least one record. Pairs that never co-occur are not tested and do
#    not enlarge the family (counting them as p = 1 inflates every q-value).
# 2. The edge filter reads the adjusted p-value, never the raw one.
#
# The test is one-sided (enrichment): the network's edges are positive
# associations, and a depleted pair must not score as significant.

pair_row <- function(ed, a, b) {
  row <- ed[ed[["item_a"]] == a & ed[["item_b"]] == b, , drop = FALSE]
  stopifnot(nrow(row) == 1L)
  row
}

test_that("only co-occurring pairs are tested, and each once with item_a < item_b", {
  ed <- compute_pairwise_associations(bh_fixture(9L, 30L))
  pairs <- paste(ed[["item_a"]], ed[["item_b"]])
  expect_setequal(pairs, c("A B", "A C", "C D", "B C", "D E"))
  expect_true(all(ed[["item_a"]] < ed[["item_b"]]))
  expect_true(all(ed[["observed"]] >= 1L))
})

test_that("p_value is the one-sided (enrichment) Fisher exact test", {
  ev <- bh_fixture(9L, 30L)
  ed <- compute_pairwise_associations(ev)
  n <- length(unique(ev[["id"]]))
  for (i in seq_len(nrow(ed))) {
    expect_equal(ed[["p_value"]][[i]],
                 fisher_greater(ed[["observed"]][[i]], ed[["n_a"]][[i]], ed[["n_b"]][[i]], n),
                 tolerance = 1e-12)
  }
  # B-C is strongly depleted: two-sided it would be highly significant.
  expect_gt(pair_row(ed, "B", "C")[["p_value"]], 0.5)
})

test_that("p_adjusted is BH over exactly the tested pairs", {
  ed <- compute_pairwise_associations(bh_fixture(9L, 30L))
  expect_equal(ed[["p_adjusted"]], stats::p.adjust(ed[["p_value"]], method = "BH"),
               tolerance = 1e-14)
  widened <- stats::p.adjust(c(ed[["p_value"]], rep(1, 5L)), method = "BH")[seq_len(nrow(ed))]
  expect_false(isTRUE(all.equal(ed[["p_adjusted"]], widened)))
})

test_that("p_adjusted is not the raw p-value", {
  ed <- compute_pairwise_associations(bh_fixture(9L, 30L))
  ab <- pair_row(ed, "A", "B")
  expect_gt(ab[["p_adjusted"]], ab[["p_value"]])
  expect_true(all(ed[["p_adjusted"]] >= ed[["p_value"]]))
  expect_true(all(ed[["p_adjusted"]] <= 1))
})

test_that("widening the family to untested pairs would lose a real edge", {
  # A-B: q ~ 0.034 over the five tested pairs, ~ 0.069 over all ten.
  net <- build_cooccurrence_network(bh_fixture(9L, 30L), weight_metric = "phi",
                                    min_count = 1L, alpha = 0.05)
  e <- network_edges(net)
  expect_identical(paste(e[["item_a"]], e[["item_b"]]), "A B")
  ab <- pair_row(net[["edge_data"]], "A", "B")
  expect_lt(ab[["p_adjusted"]], 0.05)
})

test_that("the edge filter reads q, not raw p", {
  # A-B: raw p ~ 0.028 < alpha, q ~ 0.14 > alpha -> not an edge.
  net <- build_cooccurrence_network(bh_fixture(8L, 20L), weight_metric = "phi",
                                    min_count = 1L, alpha = 0.05)
  ab <- pair_row(net[["edge_data"]], "A", "B")
  expect_lt(ab[["p_value"]], 0.05)
  expect_gte(ab[["p_adjusted"]], 0.05)
  expect_identical(nrow(network_edges(net)), 0L)
})

test_that("edges pass all three filters: count, q, and lift", {
  ev <- bh_fixture(9L, 30L)
  ed <- compute_pairwise_associations(ev)
  for (min_count in c(1L, 9L, 10L)) {
    net <- build_cooccurrence_network(ev, weight_metric = "phi",
                                      min_count = min_count, alpha = 0.05)
    keep <- ed[ed[["observed"]] >= min_count & ed[["p_adjusted"]] < 0.05 &
                 ed[["lift"]] > 1, , drop = FALSE]
    e <- network_edges(net)
    expect_setequal(paste(e[["item_a"]], e[["item_b"]]),
                    paste(keep[["item_a"]], keep[["item_b"]]))
  }
})

test_that("a depleted pair is never an edge, even at alpha = 1", {
  net <- build_cooccurrence_network(bh_fixture(9L, 30L), weight_metric = "phi",
                                    min_count = 1L, alpha = 1)
  e <- network_edges(net)
  pairs <- paste(e[["item_a"]], e[["item_b"]])
  expect_false("B C" %in% pairs)
  expect_true(all(pair_row(net[["edge_data"]], "B", "C")[["lift"]] < 1))
})

test_that("a record's duplicate events count once", {
  ev <- bh_fixture(9L, 30L)
  dup <- rbind(ev, ev[ev[["item"]] == "A", , drop = FALSE])
  expect_identical(compute_pairwise_associations(dup),
                   compute_pairwise_associations(ev))
})

test_that("unsupported tests and corrections raise rather than fall back", {
  ev <- bh_fixture(9L, 30L)
  expect_error(compute_pairwise_associations(ev, test = "chisq"), "test = 'chisq'.*not ported")
  expect_error(compute_pairwise_associations(ev, correction = "bonferroni"),
               "correction = 'bonferroni'.*not ported")
  expect_error(compute_pairwise_associations(ev, test = "bogus"), "unsupported test")
  expect_error(build_cooccurrence_network(ev, weight_metric = "jaccard"), "weight_metric")
})

test_that("missing ids or items raise rather than drop rows", {
  ev <- bh_fixture(9L, 30L)
  ev_na <- ev
  ev_na[["item"]][[1L]] <- NA
  expect_error(compute_pairwise_associations(ev_na), "NA")
  ev_na <- ev
  ev_na[["id"]][[1L]] <- NA
  expect_error(compute_pairwise_associations(ev_na), "NA")
})
