# Group-stratified networks and the caller-supplied exclusion of
# group-exclusive items. Written from the specification (the consuming
# project's AGENTS.md and port design, section 3) before the code, not
# transliterated from the Julia implementation. Cases the Julia suite also
# has are paired with it by name; cases the port adds are marked "(R)".

#' Two groups, "A" and "B". "OnlyA" is an item that can occur only in group
#' A; two records in group B carry it anyway (a recording error), one with
#' "X" and one alone. Group B is listed first in the table, so any order in
#' the result comes from sorting, not from the input.
exclusion_fixture <- function() {
  b <- events_from_records(list(
    list(items = c("X", "Y"), n = 12L),
    list(items = "X", n = 6L),
    list(items = "Y", n = 6L),
    list(items = "Z", n = 10L),
    list(items = c("OnlyA", "X"), n = 1L),
    list(items = "OnlyA", n = 1L)
  ))
  b[["Group"]] <- "B"
  a <- events_from_records(list(
    list(items = c("OnlyA", "X"), n = 10L),
    list(items = "OnlyA", n = 5L),
    list(items = "X", n = 5L),
    list(items = "Y", n = 8L)
  ))
  a[["id"]] <- a[["id"]] + 1000L
  a[["Group"]] <- "A"
  rbind(b, a)
}

excl <- list(A = "OnlyA")

test_that("stratified_network_analysis: named by group, in sorted order (R)", {
  ev <- exclusion_fixture()
  res <- stratified_network_analysis(seed = 1L, ev, exclusive_items = excl, min_count = 1L)
  expect_identical(names(res), c("A", "B"))
  for (g in names(res)) {
    expect_named(res[[g]], c("network", "communities", "metrics", "excluded_items",
                             "n_excluded_events"))
    expect_s3_class(res[[g]][["network"]], "cooccurrence_network")
  }
})

test_that("each stratum records what the exclusion removed from it (R)", {
  ev <- exclusion_fixture()
  res <- stratified_network_analysis(seed = 1L, ev, exclusive_items = excl, min_count = 1L)
  expect_identical(res[["B"]][["excluded_items"]], "OnlyA")
  expect_identical(res[["B"]][["n_excluded_events"]], 2L)
  expect_identical(res[["A"]][["excluded_items"]], character(0))
  expect_identical(res[["A"]][["n_excluded_events"]], 0L)
  none <- stratified_network_analysis(seed = 1L, ev, exclusive_items = NULL, min_count = 1L)
  expect_identical(none[["B"]][["excluded_items"]], character(0))
  expect_identical(none[["B"]][["n_excluded_events"]], 0L)
})

test_that("stratified_network_analysis: order is C-locale byte order (R)", {
  ev <- exclusion_fixture()
  ev[["Group"]] <- ifelse(ev[["Group"]] == "A", "b", "B")
  res <- stratified_network_analysis(ev, exclusive_items = list(b = "OnlyA"), min_count = 1L,
                                     seed = 1L)
  expect_identical(names(res), c("B", "b"))
})

test_that("exclusive_items has no default: forgetting it errors (R)", {
  ev <- exclusion_fixture()
  expect_error(stratified_network_analysis(seed = 1L, ev, min_count = 1L), "exclusive_items")
})

test_that("exclusive_items = NULL is the explicit 'no exclusion' (R)", {
  ev <- exclusion_fixture()
  res <- stratified_network_analysis(seed = 1L, ev, exclusive_items = NULL, min_count = 1L)
  expect_true("OnlyA" %in% res[["B"]][["network"]][["edge_data"]][["item_a"]])
  expect_identical(res[["B"]][["network"]][["n_records"]], 36L)
})

test_that("each stratum is the network of that group's events, exclusions removed", {
  ev <- exclusion_fixture()
  res <- stratified_network_analysis(seed = 1L, ev, exclusive_items = excl, min_count = 1L,
                                     weight_metric = "phi")
  b_events <- ev[ev[["Group"]] == "B" & ev[["item"]] != "OnlyA", c("id", "item")]
  expect_equal(res[["B"]][["network"]],
               build_cooccurrence_network(b_events, weight_metric = "phi", min_count = 1L))
  a_events <- ev[ev[["Group"]] == "A", c("id", "item")]
  expect_equal(res[["A"]][["network"]],
               build_cooccurrence_network(a_events, weight_metric = "phi", min_count = 1L))
})

test_that("an excluded item is dropped from the other group, not scored there (R)", {
  ev <- exclusion_fixture()
  res <- stratified_network_analysis(seed = 1L, ev, exclusive_items = excl, min_count = 1L)
  b <- res[["B"]][["network"]]
  expect_false("OnlyA" %in% c(b[["edge_data"]][["item_a"]], b[["edge_data"]][["item_b"]]))
  expect_false("OnlyA" %in% names(b[["prevalence"]]))
  # The record holding OnlyA alone has no event left in group B, so it is not
  # a record of that stratum; the OnlyA + X record keeps its X.
  expect_identical(b[["n_records"]], 35L)
  # The group the item belongs to keeps it.
  a <- res[["A"]][["network"]]
  expect_true("OnlyA" %in% a[["edge_data"]][["item_a"]])
  expect_identical(a[["n_records"]], 28L)
})

test_that("known answer: X-Y in group B with and without the exclusion (R)", {
  ev <- exclusion_fixture()
  row <- function(res) {
    ed <- res[["B"]][["network"]][["edge_data"]]
    ed[ed[["item_a"]] == "X" & ed[["item_b"]] == "Y", , drop = FALSE]
  }
  with_excl <- row(stratified_network_analysis(ev, exclusive_items = excl, min_count = 1L,
                                               seed = 1L))
  without <- row(stratified_network_analysis(seed = 1L, ev, exclusive_items = NULL, min_count = 1L))
  # Group B: 12 X+Y, 6 X, 6 Y, 10 Z, 1 OnlyA+X, 1 OnlyA. Excluded: N = 35
  # (the OnlyA-only record is gone), n_X = 19, n_Y = 18.
  expect_equal(with_excl[["p_value"]], fisher_greater(12, 19, 18, 35))
  expect_equal(without[["p_value"]], fisher_greater(12, 19, 18, 36))
  expect_false(isTRUE(all.equal(with_excl[["p_value"]], without[["p_value"]])))
})

test_that("vocabulary guard: no exclusive item in the data errors (R)", {
  ev <- exclusion_fixture()
  wrong <- list(A = c("Only-A", "onlya"))
  expect_error(stratified_network_analysis(seed = 1L, ev, exclusive_items = wrong, min_count = 1L),
               "none of the exclusive_items")
})

test_that("vocabulary guard: some exclusive items not in the data warns, naming them (R)", {
  ev <- exclusion_fixture()
  partial <- list(A = c("OnlyA", "Ghost"))
  expect_warning(stratified_network_analysis(ev, exclusive_items = partial, min_count = 1L,
                                             seed = 1L),
                 "Ghost")
})

test_that("exclusive_items keyed by a group that is not in the data errors (R)", {
  ev <- exclusion_fixture()
  expect_error(stratified_network_analysis(seed = 1L, ev, exclusive_items = list(a = "OnlyA"),
                                           min_count = 1L),
               "not a value of")
  expect_error(stratified_network_analysis(seed = 1L, ev, exclusive_items = list("OnlyA"),
                                           min_count = 1L),
               "named")
})

test_that("an item listed as exclusive to two groups errors (R)", {
  ev <- exclusion_fixture()
  expect_error(stratified_network_analysis(seed = 1L, ev, exclusive_items = list(A = "X", B = "X"),
                                           min_count = 1L),
               "more than one group")
})

test_that("a missing or NA group column errors rather than dropping rows (R)", {
  ev <- exclusion_fixture()
  expect_error(stratified_network_analysis(ev, exclusive_items = excl, group_col = "Sex",
                                           seed = 1L),
               "no 'Sex' column")
  ev[["Group"]][[3L]] <- NA
  expect_error(stratified_network_analysis(seed = 1L, ev, exclusive_items = excl), "NA")
})

test_that("stratified_network_analysis on the shared fixture", {
  # Julia: FIXTURE_EXCLUSIVE; GA1 occurs only in group A, GB1 only in B.
  res <- stratified_network_analysis(seed = 1L, make_test_event_df(),
                                     exclusive_items = list(A = "GA1", B = "GB1"),
                                     min_count = 1L, alpha = 1)
  expect_identical(names(res), c("A", "B"))
  expect_s3_class(res[["A"]][["network"]], "cooccurrence_network")
  expect_s3_class(res[["B"]][["network"]], "cooccurrence_network")
  tested <- function(g) {
    ed <- res[[g]][["network"]][["edge_data"]]
    c(ed[["item_a"]], ed[["item_b"]])
  }
  expect_false("GB1" %in% tested("A"))
  expect_false("GA1" %in% tested("B"))
  expect_true("GB1" %in% tested("B"))
  # Item01 + Item02 is planted in group B only; GA1 + Item05 in group A only.
  # Tested pairs, not edges: within its own group neither pair beats its
  # expectation (lift < 1), because the group-only background is small.
  pairs <- function(g) {
    e <- res[[g]][["network"]][["edge_data"]]
    paste(e[["item_a"]], e[["item_b"]])
  }
  expect_true("Item01 Item02" %in% pairs("B"))
  expect_false("Item01 Item02" %in% pairs("A"))
  expect_true("GA1 Item05" %in% pairs("A"))
  for (g in names(res)) {
    m <- res[[g]][["metrics"]]
    expect_identical(m[["item"]], res[[g]][["network"]][["items"]])
  }
})

test_that("network stratification over 3 groups", {
  # Julia's three-group fixture: per group, 8 S1+S2 records plus its own
  # items; one group-A record carries B's exclusive item B1.
  rows <- list()
  add <- function(g, items, n) {
    for (k in seq_len(n)) rows[[length(rows) + 1L]] <<- list(g = g, items = items)
  }
  for (g in c("A", "B", "C")) {
    add(g, c("S1", "S2"), 8L)
    add(g, paste0(g, c("1", "2")), 4L)
  }
  add("A", c("A1", "B1"), 1L)
  ev <- do.call(rbind, lapply(seq_along(rows), function(i) {
    data.frame(id = i, item = rows[[i]][["items"]], Group = rows[[i]][["g"]])
  }))
  excl3 <- list(A = c("A1", "A2"), B = c("B1", "B2"), C = c("C1", "C2"))
  res <- stratified_network_analysis(ev, exclusive_items = excl3, min_count = 1L, alpha = 1,
                                     seed = 1L)
  expect_identical(names(res), c("A", "B", "C"))
  for (g in names(res)) expect_s3_class(res[[g]][["network"]], "cooccurrence_network")
  expect_setequal(res[["A"]][["excluded_items"]], c("B1", "B2", "C1", "C2"))
  expect_identical(res[["A"]][["n_excluded_events"]], 1L)
  expect_false("B1" %in% res[["A"]][["network"]][["edge_data"]][["item_b"]])
})

test_that("run_network_pipeline: stratify_by_group = TRUE still raises (R)", {
  expect_error(run_network_pipeline(make_test_event_df()), "compare_networks")
})

test_that("plot_network_comparison takes two named networks (R)", {
  ev <- exclusion_fixture()
  res <- stratified_network_analysis(ev, exclusive_items = excl, min_count = 1L, alpha = 1,
                                     seed = 1L)
  p <- plot_network_comparison(res)
  expect_s3_class(p[["plot"]], "ggplot")
  # Read back from the built plot: panel order and edge colour as drawn.
  expect_identical(p[["panels"]][["group"]], c("A", "B"))
  expect_identical(p[["panels"]][["panel"]], 1:2)
  expect_identical(p[["panels"]][["edge_colour"]], c("#B3B3B3", "#B3B3B3"))
  rev_p <- plot_network_comparison(res[c("B", "A")])
  expect_identical(rev_p[["panels"]][["group"]], c("B", "A"))
  expect_identical(rev_p[["panels"]][["panel"]], 1:2)
  expect_error(plot_network_comparison(unname(res)), "named")
  expect_error(plot_network_comparison(res["A"]), "two")
})

#' Two triangles joined by one weak edge: Louvain finds two communities.
two_triangle_stratum <- function() {
  items <- paste0("n", 1:6)
  adj <- matrix(0, 6, 6, dimnames = list(items, items))
  for (e in list(c(1, 2), c(2, 3), c(1, 3), c(4, 5), c(5, 6), c(4, 6))) {
    adj[e[[1]], e[[2]]] <- adj[e[[2]], e[[1]]] <- 1
  }
  adj[3, 4] <- adj[4, 3] <- 0.1
  net <- structure(list(items = items, adjacency = adj, edge_data = data.frame(),
                        prevalence = stats::setNames(rep(1L, 6), items), n_records = 1L,
                        weight_metric = "lift", min_count = 1L, alpha = 1),
                   class = "cooccurrence_network")
  list(network = net, communities = detect_communities(net, seed = 1L))
}

test_that("plot_network_comparison draws as Julia: edge alpha, node colour by community (R)", {
  s <- two_triangle_stratum()
  expect_identical(s[["communities"]][["n_communities"]], 2L)
  p <- plot_network_comparison(list(A = s, B = s))[["plot"]]
  seg <- ggplot2::layer_data(p, 1L)
  expect_true(all(seg[["alpha"]] == 0.6))
  pts <- ggplot2::layer_data(p, 2L)
  # Julia: palette[((c - 1) % 7) + 1] of Makie's wong_colors().
  wong <- c("#0072B2", "#E69F00", "#009E73", "#CC79A7", "#56B4E9", "#D55E00", "#F0E442")
  comm <- s[["communities"]][["assignments"]][s[["network"]][["items"]]]
  expected <- rep(wong[(comm - 1L) %% 7L + 1L], 2L)
  expect_identical(toupper(pts[["colour"]]), expected)
  expect_length(unique(pts[["colour"]]), 2L)
})

test_that("community colours recycle Makie's seven wong colours by label (R)", {
  expect_identical(CooccurrenceAnalysis:::.community_palette,
                   c("#0072B2", "#E69F00", "#009E73", "#CC79A7", "#56B4E9", "#D55E00",
                     "#F0E442"))
})
