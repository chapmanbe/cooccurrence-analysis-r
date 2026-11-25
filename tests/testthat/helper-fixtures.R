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

#' Each record's distinct items, read straight from the event table
#' (independent of the package). Records with fewer than `min_items` distinct
#' items are dropped, as the mining functions drop them.
record_sets <- function(ev, min_items = 1L) {
  sets <- lapply(split(as.character(ev[["item"]]), ev[["id"]]), unique)
  unname(sets[lengths(sets) >= min_items])
}

#' Number of records in `sets` holding every one of `items`.
n_holding <- function(sets, items) {
  sum(vapply(sets, function(s) all(items %in% s), logical(1)))
}

#' Every pair of items that co-occurs in at least one record, as sorted
#' "a|b" keys in C-locale order.
cooccurring_pairs <- function(sets) {
  keys <- unlist(lapply(sets, function(s) {
    if (length(s) < 2L) return(character(0))
    pr <- utils::combn(sort(s, method = "radix"), 2L)
    paste(pr[1L, ], pr[2L, ], sep = "|")
  }))
  sort(unique(keys), method = "radix")
}

#' The 2x2 table for items `a` and `b` over `sets`:
#' rows a present/absent, columns b present/absent.
pair_table <- function(sets, a, b) {
  n <- length(sets)
  n11 <- n_holding(sets, c(a, b))
  n_a <- n_holding(sets, a)
  n_b <- n_holding(sets, b)
  matrix(c(n11, n_a - n11, n_b - n11, n - n_a - n_b + n11), 2L, byrow = TRUE)
}

#' Pearson chi-square p-value (1 df, no continuity correction) from the
#' definition, sum((O - E)^2 / E), without stats::chisq.test().
pearson_p <- function(tab) {
  e <- outer(rowSums(tab), colSums(tab)) / sum(tab)
  stats::pchisq(sum((tab - e)^2 / e), df = 1, lower.tail = FALSE)
}

#' Every itemset of 2 to `max_len` items held by at least `min_n` of `sets`,
#' enumerated from the records themselves, as sorted "a+b" keys with counts.
frequent_sets <- function(sets, min_n, max_len = 4L) {
  keys <- unlist(lapply(sets, function(s) {
    s <- sort(s, method = "radix")
    sizes <- seq_len(min(length(s), max_len))
    unlist(lapply(sizes[sizes >= 2L], function(k) {
      apply(utils::combn(s, k), 2L, paste, collapse = "+")
    }))
  }))
  counts <- table(keys)
  counts <- counts[counts >= min_n]
  stats::setNames(as.integer(counts), names(counts))
}

#' Every rule LHS => RHS (one RHS item) from the itemsets of `frequent_sets()`,
#' with N, Support, Confidence, Coverage and Lift computed from record counts.
rules_from_records <- function(sets, min_n, min_conf, max_len = 4L) {
  fs <- frequent_sets(sets, min_n, max_len)
  n <- length(sets)
  rows <- lapply(names(fs), function(key) {
    items <- strsplit(key, "+", fixed = TRUE)[[1L]]
    lapply(items, function(rhs) {
      lhs <- setdiff(items, rhs)
      n_lhs <- n_holding(sets, lhs)
      conf <- fs[[key]] / n_lhs
      data.frame(rule = paste(paste(lhs, collapse = "+"), rhs, sep = ">"), N = fs[[key]],
                 Support = fs[[key]] / n, Confidence = conf, Coverage = n_lhs / n,
                 Lift = conf / (n_holding(sets, rhs) / n))
    })
  })
  out <- do.call(rbind, unlist(rows, recursive = FALSE))
  out <- out[out[["Confidence"]] >= min_conf, ]
  out <- out[order(out[["rule"]], method = "radix"), ]
  rownames(out) <- NULL
  out
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

#' Ring of `n_cliques` complete graphs of `k` nodes, each joined to the next
#' by one unit edge. Phase-1-only local moving stops at one community per
#' clique; full Louvain/Leiden aggregates neighbors and scores higher.
ring_of_cliques <- function(n_cliques, k) {
  n <- n_cliques * k
  adj <- matrix(0, n, n)
  for (c in seq_len(n_cliques) - 1L) {
    idx <- c * k + seq_len(k)
    adj[idx, idx] <- 1
    nxt <- (c + 1L) %% n_cliques * k + 1L
    adj[c * k + k, nxt] <- adj[nxt, c * k + k] <- 1
  }
  diag(adj) <- 0
  adj
}
