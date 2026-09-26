# Network construction from item co-occurrence data.
#
# An event table has one row per (record, item) event: columns `id` (the
# record) and `item`. A record may hold an item more than once; it counts
# once. Every count here is a count of records.

#' Phi coefficient of a 2x2 contingency table
#'
#' Layout, as a matrix with rows "has A / no A" and columns "has B / no B":
#' `ct[1, 1] = n11`, `ct[1, 2] = n10`, `ct[2, 1] = n01`, `ct[2, 2] = n00`.
#' Ranges from -1 to 1; returns 0 when a margin is empty (zero denominator).
#' @export
phi_coefficient <- function(ct) {
  stopifnot(is.matrix(ct), identical(dim(ct), c(2L, 2L)))
  .phi(ct[1L, 1L], ct[1L, 2L], ct[2L, 1L], ct[2L, 2L])
}

#' Vectorized phi from the four cells, in double precision throughout
#' @keywords internal
.phi <- function(n11, n10, n01, n00) {
  n11 <- as.numeric(n11)
  n10 <- as.numeric(n10)
  n01 <- as.numeric(n01)
  n00 <- as.numeric(n00)
  denom <- sqrt((n11 + n10) * (n01 + n00) * (n11 + n01) * (n10 + n00))
  out <- (n11 * n00 - n10 * n01) / denom
  out[denom == 0] <- 0
  out
}

#' Odds ratio `(a d) / (b c)` of a 2x2 table `[a b; c d]`
#'
#' With `correction = "haldane"` (the default, as in the Julia reference),
#' the Haldane-Anscombe correction adds 0.5 to all four cells whenever any
#' cell is zero, so an empty cell gives a large but finite ratio instead of
#' `Inf` or `NaN`. `correction = "none"` returns the raw ratio, which may be
#' `Inf` (empty `b c`) or `NaN` (0/0). Any other correction raises.
#' @export
odds_ratio <- function(ct, correction = "haldane") {
  stopifnot(is.matrix(ct), identical(dim(ct), c(2L, 2L)))
  .odds_ratio(ct[1L, 1L], ct[1L, 2L], ct[2L, 1L], ct[2L, 2L], correction)
}

#' Vectorized odds ratio from the four cells, in double precision
#' @keywords internal
.odds_ratio <- function(a, b, c, d, correction) {
  .check_choice(correction, c("haldane", "none"), "correction")
  a <- as.numeric(a)
  b <- as.numeric(b)
  c <- as.numeric(c)
  d <- as.numeric(d)
  raw <- (a * d) / (b * c)
  if (identical(correction, "none")) return(raw)
  zero <- a == 0 | b == 0 | c == 0 | d == 0
  raw[zero] <- ((a[zero] + 0.5) * (d[zero] + 0.5)) / ((b[zero] + 0.5) * (c[zero] + 0.5))
  raw
}

#' Benjamini-Hochberg adjusted p-values over one family of tests
#'
#' The family is exactly the vector passed in: every element is one test
#' that was performed. Do not pad it with p = 1 for tests that were never
#' run -- that enlarges the family and inflates every q-value. Raises on
#' missing or out-of-range p-values rather than letting `p.adjust()` shrink
#' the family around an `NA`. An empty family gives an empty result.
#' @export
bh_adjust <- function(p) {
  if (!is.numeric(p)) stop("bh_adjust: p must be numeric", call. = FALSE)
  if (anyNA(p)) stop("bh_adjust: p contains NA; every tested pair needs a p-value", call. = FALSE)
  if (any(p < 0 | p > 1)) stop("bh_adjust: p-values must lie in [0, 1]", call. = FALSE)
  stats::p.adjust(p, method = "BH")
}

#' Validate test / correction / timing options, naming what is not ported
#'
#' The Julia reference also offers a chi-square test, Bonferroni and Holm
#' corrections, and `timing_filter = "concurrent"` / `"sequential"`. The
#' analysis never uses them (the v2 network fit passes Fisher + BH and the
#' default timing), so they are not ported; asking for one raises with a
#' message that says so rather than falling back to a supported option.
#' @keywords internal
.check_options <- function(test, correction, timing_filter) {
  not_ported <- function(arg, value) {
    stop(sprintf(paste0("%s = '%s' exists in the Julia reference but is not ported to R: ",
                        "the analysis pipeline never uses it"), arg, value), call. = FALSE)
  }
  if (identical(test, "chisq")) not_ported("test", test)
  if (isTRUE(correction %in% c("bonferroni", "holm"))) not_ported("correction", correction)
  if (isTRUE(timing_filter %in% c("concurrent", "sequential"))) {
    not_ported("timing_filter", timing_filter)
  }
  .check_choice(test, "fisher", "test")
  .check_choice(correction, "bh", "correction")
  .check_choice(timing_filter, "all", "timing_filter")
}

#' Record-level counts behind every pairwise statistic
#'
#' Returns `items` (sorted, C-locale), `n_item` (records holding each item,
#' in `items` order), `n_records`, and `pair_counts`, a `K x K` integer
#' matrix whose upper triangle holds the number of records holding both
#' items. Raises on `NA` ids or items: dropping them would silently shrink
#' every count.
#' @keywords internal
.record_counts <- function(event_df) {
  for (col in c("id", "item")) {
    if (!(col %in% names(event_df))) {
      stop(sprintf("event table has no '%s' column (columns: %s)", col,
                   paste(names(event_df), collapse = ", ")), call. = FALSE)
    }
    if (anyNA(event_df[[col]])) {
      stop(sprintf("event table column '%s' contains NA; every event needs one", col),
           call. = FALSE)
    }
  }
  item_chr <- as.character(event_df[["item"]])
  items <- .sort_c(unique(item_chr))
  k <- length(items)
  code <- match(item_chr, items)
  ids <- event_df[["id"]]
  unique_ids <- unique(ids)
  rec <- match(ids, unique_ids)
  n_records <- length(unique_ids)

  # One row per (record, item): a repeated event counts once.
  keep <- !duplicated(as.numeric(rec) * (k + 1) + code)
  rec <- rec[keep]
  code <- code[keep]

  n_item <- tabulate(code, nbins = k)
  size <- tabulate(rec, nbins = n_records)

  pair_counts <- matrix(0L, k, k, dimnames = list(items, items))
  multi <- size[rec] >= 2L
  if (any(multi)) {
    rec_m <- rec[multi]
    code_m <- code[multi]
    o <- order(rec_m, code_m, method = "radix")
    rec_m <- rec_m[o]
    code_m <- code_m[o]
    sizes <- size[rec_m]
    # Records with the same number of items form a matrix (one row each,
    # items ascending); every column pair i < j is one within-record pair.
    for (s in sort(unique(sizes))) {
      m <- matrix(code_m[sizes == s], ncol = s, byrow = TRUE)
      for (i in seq_len(s - 1L)) {
        for (j in (i + 1L):s) {
          idx <- (m[, j] - 1L) * k + m[, i]
          pair_counts <- pair_counts + matrix(tabulate(idx, nbins = k * k), k, k)
        }
      }
    }
  }
  list(items = items, n_item = n_item, n_records = n_records, pair_counts = pair_counts)
}

#' Association statistics for every pair of items that co-occur
#'
#' One row per unordered pair `item_a < item_b` (C-locale order) that
#' co-occurs in at least one record, ordered by `item_a` then `item_b`.
#' Pairs that never co-occur are not tested and have no row; they are not
#' part of the multiple-testing family.
#'
#' Columns: `item_a`, `item_b`, `observed` (records holding both),
#' `expected` (`n_a * n_b / N`, rounded to 2 places), `lift` (`observed /
#' expected`, 4 places), `phi` (4 places), `odds_ratio` ([odds_ratio()] with its
#' default Haldane-Anscombe correction, so always finite; 4 places),
#' `p_value` (one-sided Fisher exact test for enrichment, unrounded), `n_a`,
#' `n_b` (records holding each item), and `p_adjusted` (Benjamini-Hochberg
#' over the rows of this table). `N` is the number of distinct records.
#'
#' `test` must be `"fisher"`, `correction` `"bh"`, and `timing_filter`
#' `"all"`; the Julia reference's other options are not ported and raise.
#' @export
compute_pairwise_associations <- function(event_df, test = "fisher", correction = "bh",
                                          timing_filter = "all") {
  .check_options(test, correction, timing_filter)
  counts <- .record_counts(event_df)
  .associations_from_counts(counts)
}

#' @keywords internal
.associations_from_counts <- function(counts) {
  items <- counts[["items"]]
  pc <- counts[["pair_counts"]]
  hit <- which(upper.tri(pc) & pc > 0L, arr.ind = TRUE)
  hit <- hit[order(hit[, 1L], hit[, 2L], method = "radix"), , drop = FALSE]
  a <- hit[, 1L]
  b <- hit[, 2L]

  n <- as.numeric(counts[["n_records"]])
  n11 <- pc[hit]
  n_a <- counts[["n_item"]][a]
  n_b <- counts[["n_item"]][b]
  n10 <- as.numeric(n_a) - n11
  n01 <- as.numeric(n_b) - n11
  n00 <- n - n_a - n_b + n11

  expected <- as.numeric(n_a) * as.numeric(n_b) / n
  p_value <- stats::phyper(n11 - 1, n_a, n - n_a, n_b, lower.tail = FALSE)

  data.frame(
    item_a = items[a],
    item_b = items[b],
    observed = as.integer(n11),
    expected = round_digits(expected, 2L),
    lift = round_digits(n11 / expected, 4L),
    phi = round_digits(.phi(n11, n10, n01, n00), 4L),
    odds_ratio = round_digits(.odds_ratio(n11, n10, n01, n00, "haldane"), 4L),
    p_value = p_value,
    n_a = as.integer(n_a),
    n_b = as.integer(n_b),
    p_adjusted = bh_adjust(p_value)
  )
}

#' Build a weighted co-occurrence network from an event table
#'
#' Edges are the pairs that pass all three filters: `observed >= min_count`,
#' `p_adjusted < alpha`, and `lift > 1`. Vertices are the items in at least
#' one edge, sorted (isolated items carry no network information). Edge
#' weights are `lift`, or for `weight_metric = "phi"` the phi coefficient
#' floored at 0.001.
#'
#' Returns a `cooccurrence_network`: `items`, `adjacency` (symmetric weight
#' matrix, zero where there is no edge), `edge_data` (every tested pair, see
#' [compute_pairwise_associations()]), `prevalence` (records holding each
#' vertex item, named), `n_records`, and the parameters that produced it.
#' @export
build_cooccurrence_network <- function(event_df, weight_metric = "lift", min_count = 30L,
                                       alpha = 0.05, test = "fisher", correction = "bh",
                                       timing_filter = "all") {
  .check_choice(weight_metric, c("lift", "phi"), "weight_metric")
  .check_options(test, correction, timing_filter)
  counts <- .record_counts(event_df)
  edge_data <- .associations_from_counts(counts)
  prevalence <- stats::setNames(counts[["n_item"]], counts[["items"]])
  cooccurrence_network_from_edge_data(edge_data, prevalence = prevalence,
                                      n_records = counts[["n_records"]],
                                      weight_metric = weight_metric,
                                      min_count = min_count, alpha = alpha)
}

#' Rebuild a co-occurrence network from its pairwise table
#'
#' The graph is a pure function of `edge_data` and the filter parameters, so
#' a network written to disk as tables can be reconstructed exactly.
#' `prevalence` is a named integer vector covering at least every vertex
#' item; a vertex missing from it raises.
#' @export
cooccurrence_network_from_edge_data <- function(edge_data, prevalence, n_records,
                                                weight_metric, min_count, alpha) {
  .check_choice(weight_metric, c("lift", "phi"), "weight_metric")
  stopifnot(length(min_count) == 1L, length(alpha) == 1L, length(n_records) == 1L)
  sig <- edge_data[["observed"]] >= min_count &
    edge_data[["p_adjusted"]] < alpha &
    edge_data[["lift"]] > 1
  if (anyNA(sig)) stop("edge_data has NA in observed, p_adjusted, or lift", call. = FALSE)
  significant <- edge_data[sig, , drop = FALSE]

  items <- .sort_c(unique(c(significant[["item_a"]], significant[["item_b"]])))
  n <- length(items)
  adjacency <- matrix(0, n, n, dimnames = list(items, items))
  w <- if (identical(weight_metric, "phi")) {
    pmax(significant[["phi"]], 0.001)
  } else {
    significant[["lift"]]
  }
  i <- match(significant[["item_a"]], items)
  j <- match(significant[["item_b"]], items)
  adjacency[cbind(i, j)] <- w
  adjacency[cbind(j, i)] <- w

  missing_prev <- setdiff(items, names(prevalence))
  if (length(missing_prev) > 0L) {
    stop(sprintf("prevalence has no entry for vertex item(s): %s",
                 paste(missing_prev, collapse = ", ")), call. = FALSE)
  }
  structure(list(
    items = items,
    adjacency = adjacency,
    edge_data = edge_data,
    prevalence = prevalence[items],
    n_records = as.integer(n_records),
    weight_metric = weight_metric,
    min_count = min_count,
    alpha = alpha
  ), class = "cooccurrence_network")
}

#' The network's edges as a table
#'
#' One row per edge, `item_a < item_b`, ordered by `item_a` then `item_b`,
#' with its `weight`.
#' @export
network_edges <- function(net) {
  stopifnot(inherits(net, "cooccurrence_network"))
  adj <- net[["adjacency"]]
  hit <- which(upper.tri(adj) & adj != 0, arr.ind = TRUE)
  hit <- hit[order(hit[, 1L], hit[, 2L], method = "radix"), , drop = FALSE]
  data.frame(item_a = net[["items"]][hit[, 1L]],
             item_b = net[["items"]][hit[, 2L]],
             weight = adj[hit])
}

#' Number of vertices and edges
#' @export
network_size <- function(net) {
  stopifnot(inherits(net, "cooccurrence_network"))
  adj <- net[["adjacency"]]
  c(n_nodes = length(net[["items"]]), n_edges = sum(upper.tri(adj) & adj != 0))
}
