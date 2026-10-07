# Synthetic fixtures only. Every fixture is hand-built so its answer can be
# worked out on paper (or with an independent call such as stats::fisher.test)
# rather than read back from the code under test.

#' Event table (`id`, `item`) from a list of record templates.
#'
#' `records` is a list of `list(items = <character>, n = <count>)`; each
#' template becomes `n` records holding exactly those items. Ids are assigned
#' in order, so the result is deterministic.
events_from_records <- function(records) {
  rows <- list()
  next_id <- 1L
  for (r in records) {
    for (k in seq_len(r[["n"]])) {
      rows[[length(rows) + 1L]] <- data.frame(id = next_id, item = r[["items"]])
      next_id <- next_id + 1L
    }
  }
  do.call(rbind, rows)
}

#' The BH-family fixture. Items A..E; five of the ten possible pairs ever
#' co-occur, so five pairs are untested. `ab` sets how often A and B
#' co-occur and `solo_e` how many records hold E alone, which together set
#' how significant A-B is:
#'
#' - `bh_fixture(9, 30)`: A-B has one-sided p ~ 0.0069. BH over the five
#'   tested pairs gives q ~ 0.034 (an edge at alpha = 0.05); BH over all ten
#'   pairs, the five untested scored p = 1, gives q ~ 0.069 (no edge).
#' - `bh_fixture(8, 20)`: A-B has p ~ 0.028 but q ~ 0.14 -- significant
#'   only if raw p-values are (wrongly) used in place of q-values.
#'
#' B-C co-occurs once against an expectation of ~9: strongly depleted, with
#' a tiny two-sided p but a one-sided (enrichment) p near 1.
bh_fixture <- function(ab, solo_e) {
  events_from_records(list(
    list(items = c("A", "B"), n = ab),
    list(items = c("A", "C"), n = 4L),
    list(items = c("C", "D"), n = 3L),
    list(items = c("B", "C"), n = 1L),
    list(items = c("D", "E"), n = 2L),
    list(items = "A", n = 10L),
    list(items = "B", n = 10L),
    list(items = "C", n = 30L),
    list(items = "D", n = 10L),
    list(items = "E", n = solo_e)
  ))
}

#' One-sided (enrichment) Fisher p-value for a pair, computed independently
#' of the package from the 2x2 table.
fisher_greater <- function(n11, n_a, n_b, n) {
  ct <- matrix(c(n11, n_a - n11, n_b - n11, n - n_a - n_b + n11), 2L, byrow = TRUE)
  stats::fisher.test(ct, alternative = "greater")[["p.value"]]
}

#' Adjusted Rand index of two labelings, written out from the definition
#' (as the Julia suite's `_test_ari`), independent of the package.
adjusted_rand <- function(a, b) {
  tab <- table(a, b)
  c2 <- function(x) x * (x - 1) / 2
  sn <- sum(c2(tab))
  sa <- sum(c2(rowSums(tab)))
  sb <- sum(c2(colSums(tab)))
  nc2 <- c2(length(a))
  if (nc2 == 0) return(1)
  expected <- sa * sb / nc2
  denom <- 0.5 * (sa + sb) - expected
  if (denom == 0) return(1)
  (sn - expected) / denom
}

#' The Julia suite's shared fixture (`make_test_event_df()` in
#' test/runtests.jl), reduced to the columns these tests read. 50 records:
#' 30 single-item background records (group B: 8 Item01, 4 Item02, 3 Item03;
#' group A: 8 GA1, 4 Item05, 3 Item03); 10 Item01 + Item02 (B, planted strong
#' association); 5 GA1 + Item05 (A); 3 Item01 + GB1 (B); 2 Item03 + Item04
#' (one in each group). Records are numbered in that order, one id each.
make_test_event_df <- function() {
  rec <- function(items, n, group) list(items = items, n = n, group = group)
  templates <- list(
    rec("Item01", 8L, "B"), rec("Item02", 4L, "B"), rec("Item03", 3L, "B"),
    rec("GA1", 8L, "A"), rec("Item05", 4L, "A"), rec("Item03", 3L, "A"),
    rec(c("Item01", "Item02"), 10L, "B"),
    rec(c("GA1", "Item05"), 5L, "A"),
    rec(c("Item01", "GB1"), 3L, "B"),
    rec(c("Item03", "Item04"), 1L, "A"),
    rec(c("Item03", "Item04"), 1L, "B")
  )
  ev <- events_from_records(templates)
  ev[["Group"]] <- unlist(lapply(templates, function(t) {
    rep(t[["group"]], t[["n"]] * length(t[["items"]]))
  }))
  ev
}
