# Frequent itemsets and association rules. R port of the Julia reference's
# mining.jl. Enumeration is delegated to arules (Borgelt's Apriori); counts,
# thresholds and metrics are computed here from the incidence matrix, so they
# follow the reference's definitions exactly rather than arules' floating-point
# thresholds.

#' Records x items logical matrix for mining: one row per record (`id`)
#' holding at least `min_items` distinct items, columns the items those
#' records hold (C-locale order). The reference's `build_transactions`.
#' @keywords internal
.transactions <- function(event_df, min_items) {
  .check_events(event_df)
  item_chr <- as.character(event_df[["item"]])
  items <- .sort_c(unique(item_chr))
  ids <- event_df[["id"]]
  uid <- unique(ids)
  inc <- .incidence(match(ids, uid), match(item_chr, items), length(uid), length(items))
  inc <- inc[Matrix::rowSums(inc) >= min_items, , drop = FALSE]
  held <- which(Matrix::colSums(inc) > 0)
  x <- as.matrix(inc[, held, drop = FALSE]) > 0
  dimnames(x) <- list(NULL, items[held])
  x
}

#' Minimum record count for a support threshold. `min_count`, when given, is
#' used as the count itself. The Julia reference converts it to a fraction
#' and RuleMiner back with `ceil(min_support * n)`, which can land one above
#' it (min_count = 7 of n = 100 gives 8); this port does not round-trip.
#' Without it, `ceil(min_support * n)` as RuleMiner computes it.
#' @keywords internal
.min_n <- function(min_support, min_count, n) {
  if (!is.null(min_count)) {
    if (!(is.numeric(min_count) && length(min_count) == 1L && !is.na(min_count) &&
            min_count >= 1 && min_count == round(min_count))) {
      stop("min_count must be NULL or one whole number >= 1", call. = FALSE)
    }
    return(as.integer(min_count))
  }
  if (!(is.numeric(min_support) && length(min_support) == 1L && !is.na(min_support) &&
          min_support > 0 && min_support <= 1)) {
    stop("min_support must be one number in (0, 1]", call. = FALSE)
  }
  max(1L, as.integer(ceiling(min_support * n)))
}

.check_mining_args <- function(max_length, min_items) {
  for (arg in c("max_length", "min_items")) {
    v <- get(arg)
    if (!(is.numeric(v) && length(v) == 1L && !is.na(v) && v >= 2 && v == round(v))) {
      stop(sprintf("%s must be one whole number >= 2", arg), call. = FALSE)
    }
  }
}

.require_arules <- function(caller) {
  if (!requireNamespace("arules", quietly = TRUE)) {
    stop(sprintf("%s needs the 'arules' package: install.packages(\"arules\")", caller),
         call. = FALSE)
  }
}

# Records holding every item of each itemset (a list of character vectors).
.itemset_counts <- function(x, itemsets) {
  vapply(itemsets, function(s) sum(rowSums(x[, s, drop = FALSE]) == length(s)), integer(1L))
}

# Run arules::apriori with the support threshold half a record below `min_n`,
# so every itemset with at least `min_n` records is enumerated whatever
# rounding arules applies; callers then filter on exact counts.
.apriori <- function(x, min_n, parameter) {
  trans <- methods::as(x, "transactions")
  parameter[["support"]] <- (min_n - 0.5) / nrow(x)
  arules::apriori(trans, parameter = parameter, control = list(verbose = FALSE))
}

.itemset_key <- function(sets) vapply(sets, paste, character(1L), collapse = "\r")

#' Frequent itemsets of two or more items
#'
#' Records are the ids holding at least `min_items` distinct items (the
#' reference's `build_transactions`). An itemset is frequent when at least
#' `min_n` of them hold it, where `min_n` is `min_count` if given and
#' otherwise `ceiling(min_support * N)`; both thresholds are inclusive.
#' Itemsets run from 2 to `max_length` items.
#'
#' Returns a data.frame with columns `Itemset` (a list of character vectors,
#' items in C-locale order), `Support` (`N` over the number of records), `N`,
#' and `Length`, ordered by `N` descending, then `Length`, then itemset.
#' Needs the arules package.
#' @export
mine_frequent_itemsets <- function(event_df, min_support = 0.005, min_count = NULL,
                                   max_length = 4L, min_items = 2L) {
  .check_mining_args(max_length, min_items)
  .require_arules("mine_frequent_itemsets")
  x <- .transactions(event_df, min_items)
  empty <- data.frame(Support = numeric(0), N = integer(0), Length = integer(0))
  empty[["Itemset"]] <- list()
  empty <- empty[, c("Itemset", "Support", "N", "Length")]
  if (nrow(x) == 0L) return(empty)
  min_n <- .min_n(min_support, min_count, nrow(x))
  found <- .apriori(x, min_n, list(minlen = 2L, maxlen = as.integer(max_length),
                                   target = "frequent itemsets"))
  sets <- lapply(arules::LIST(arules::items(found), decode = TRUE), .sort_c)
  if (length(sets) == 0L) return(empty)
  n_hold <- .itemset_counts(x, sets)
  keep <- n_hold >= min_n
  sets <- sets[keep]
  n_hold <- n_hold[keep]
  len <- lengths(sets)
  o <- order(-n_hold, len, .itemset_key(sets), method = "radix")
  out <- data.frame(Support = n_hold[o] / nrow(x), N = n_hold[o], Length = as.integer(len[o]))
  out[["Itemset"]] <- unname(sets[o])
  rownames(out) <- NULL
  out[, c("Itemset", "Support", "N", "Length")]
}

#' Association rules with a non-empty left-hand side
#'
#' Rules `LHS => RHS` with one RHS item, built from itemsets of 2 to
#' `max_length` items held by at least `min_n` records (as in
#' [mine_frequent_itemsets()]), keeping those with `Confidence >=
#' min_confidence`. Metrics follow the Julia reference (RuleMiner):
#' `Support = N / n`, `Confidence = N / cov`, `Coverage = cov / n` (the
#' share of records holding the LHS), and `Lift = (N / n) / ((cov / n) *
#' (n_rhs / n))`, where `N` holds LHS and RHS together, `cov` the LHS, and
#' `n_rhs` the RHS.
#'
#' Returns a data.frame with columns `LHS` (a list of character vectors, C-locale
#' order), `RHS`, `Support`, `Confidence`, `Coverage`, `Lift`, `N`, `Length`
#' (LHS size + 1), ordered by `N` descending, then `LHS`, then `RHS`. Needs
#' the arules package.
#' @export
mine_association_rules <- function(event_df, min_support = 0.005, min_confidence = 0.1,
                                   min_count = NULL, max_length = 4L, min_items = 2L) {
  .check_mining_args(max_length, min_items)
  if (!(is.numeric(min_confidence) && length(min_confidence) == 1L && !is.na(min_confidence) &&
          min_confidence >= 0 && min_confidence <= 1)) {
    stop("min_confidence must be one number in [0, 1]", call. = FALSE)
  }
  .require_arules("mine_association_rules")
  x <- .transactions(event_df, min_items)
  empty <- data.frame(RHS = character(0), Support = numeric(0), Confidence = numeric(0),
                      Coverage = numeric(0), Lift = numeric(0), N = integer(0),
                      Length = integer(0))
  empty[["LHS"]] <- list()
  cols <- c("LHS", "RHS", "Support", "Confidence", "Coverage", "Lift", "N", "Length")
  empty <- empty[, cols]
  if (nrow(x) == 0L) return(empty)
  n <- nrow(x)
  min_n <- .min_n(min_support, min_count, n)
  # Confidence is filtered exactly below; arules gets a threshold just under it.
  found <- .apriori(x, min_n, list(confidence = max(0, min_confidence - 1e-9), minlen = 2L,
                                   maxlen = as.integer(max_length), target = "rules"))
  if (length(found) == 0L) return(empty)
  lhs <- lapply(arules::LIST(arules::lhs(found), decode = TRUE), .sort_c)
  rhs <- vapply(arules::LIST(arules::rhs(found), decode = TRUE), `[[`, character(1L), 1L)
  n_both <- .itemset_counts(x, Map(c, lhs, rhs))
  cov <- .itemset_counts(x, lhs)
  n_rhs <- as.integer(colSums(x)[rhs])
  conf <- n_both / cov
  keep <- n_both >= min_n & conf >= min_confidence
  if (!any(keep)) return(empty)
  lhs <- lhs[keep]
  rhs <- rhs[keep]
  n_both <- n_both[keep]
  cov <- cov[keep]
  n_rhs <- n_rhs[keep]
  conf <- conf[keep]
  o <- order(-n_both, .itemset_key(lhs), rhs, method = "radix")
  out <- data.frame(
    RHS = unname(rhs[o]),
    Support = n_both[o] / n,
    Confidence = conf[o],
    Coverage = cov[o] / n,
    Lift = (n_both[o] / n) / ((cov[o] / n) * (n_rhs[o] / n)),
    N = n_both[o],
    Length = lengths(lhs)[o] + 1L
  )
  out[["LHS"]] <- unname(lhs[o])
  rownames(out) <- NULL
  out[, cols]
}
