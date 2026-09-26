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
