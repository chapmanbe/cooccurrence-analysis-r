# Multi-seed HDP fits: one fit per seed, the best one by ELBO, and how much
# the fits agree -- pairwise ARI between record assignments, and whether
# each of the best fit's clusters recurs in the others. R-only: the Julia
# reference reports a single seed's fit. `adjusted_rand_index()` is the one
# exception; it ports the reference's `_adjusted_rand_index`
# (network_group_stratification.jl).

#' Adjusted Rand index between two partitions of the same elements
#'
#' Port of the reference's `_adjusted_rand_index`, edge cases included:
#' 0 when there are fewer than two elements, 1 when both partitions are
#' trivial (the denominator vanishes). Labels may be any atomic values.
#' Pair counts are doubles, not integers, so the products stay finite at a
#' million elements (the reference's `Int` product overflows near there).
#' Raises on unequal lengths or `NA` labels.
#' @export
adjusted_rand_index <- function(a, b) {
  n <- length(a)
  if (n != length(b)) stop("Label vectors must have same length", call. = FALSE)
  if (anyNA(a) || anyNA(b)) stop("labels must not be NA", call. = FALSE)
  if (n < 2L) return(0)
  ia <- match(a, unique(a))
  ib <- match(b, unique(b))
  na <- max(ia)
  nij <- tabulate(ia + na * (ib - 1L), nbins = na * max(ib))
  c2 <- function(x) {
    x <- as.numeric(x)
    sum(x * (x - 1) / 2)
  }
  sum_nij <- c2(nij)
  sum_a <- c2(tabulate(ia))
  sum_b <- c2(tabulate(ib))
  expected <- sum_a * sum_b / c2(n)
  denom <- 0.5 * (sum_a + sum_b) - expected
  if (denom == 0) return(1)
  (sum_nij - expected) / denom
}

#' Solve a linear sum assignment problem
#'
#' The minimum-total-cost matching of rows to columns of `cost` (a finite
#' numeric matrix), by the Hungarian method with potentials, O(n^2 m) for
#' n rows and m >= n columns. Each row gets a distinct column; when there
#' are more rows than columns, the rows left over get `NA`. Among equal-cost
#' optima the one returned is deterministic but otherwise unspecified.
#'
#' Returns an integer vector, one entry per row: its column.
#' @export
solve_lsap <- function(cost) {
  if (!is.matrix(cost) || !is.numeric(cost) || any(!is.finite(cost))) {
    stop("cost must be a finite numeric matrix", call. = FALSE)
  }
  n <- nrow(cost)
  m <- ncol(cost)
  if (n == 0L || m == 0L) return(rep(NA_integer_, n))
  if (n > m) {
    col_of_row <- rep(NA_integer_, n)
    col_of_row[.hungarian(t(cost))] <- seq_len(m)
    return(col_of_row)
  }
  .hungarian(cost)
}

# The Hungarian method for n <= m (rows <= columns): the column of each row.
# Every array carries a sentinel in position 1 (index 0 in the textbook
# statement), so row i is u[i + 1], and column j is v[j + 1] and p[j + 1].
.hungarian <- function(a) {
  n <- nrow(a)
  m <- ncol(a)
  u <- numeric(n + 1L)
  v <- numeric(m + 1L)
  p <- integer(m + 1L) # p[j + 1]: the row matched to column j (0: none yet)
  way <- integer(m + 1L)
  for (i in seq_len(n)) {
    p[[1L]] <- i
    j0 <- 0L
    minv <- rep(Inf, m + 1L)
    used <- rep(FALSE, m + 1L)
    repeat {
      used[[j0 + 1L]] <- TRUE
      i0 <- p[[j0 + 1L]]
      delta <- Inf
      j1 <- 0L
      for (j in seq_len(m)) {
        if (!used[[j + 1L]]) {
          cur <- a[i0, j] - u[[i0 + 1L]] - v[[j + 1L]]
          if (cur < minv[[j + 1L]]) {
            minv[[j + 1L]] <- cur
            way[[j + 1L]] <- j0
          }
          if (minv[[j + 1L]] < delta) {
            delta <- minv[[j + 1L]]
            j1 <- j
          }
        }
      }
      for (j in 0:m) {
        if (used[[j + 1L]]) {
          r <- p[[j + 1L]] + 1L
          u[[r]] <- u[[r]] + delta
          v[[j + 1L]] <- v[[j + 1L]] - delta
        } else {
          minv[[j + 1L]] <- minv[[j + 1L]] - delta
        }
      }
      j0 <- j1
      if (p[[j0 + 1L]] == 0L) break
    }
    repeat {
      j1 <- way[[j0 + 1L]]
      p[[j0 + 1L]] <- p[[j1 + 1L]]
      j0 <- j1
      if (j0 == 0L) break
    }
  }
  out <- integer(n)
  for (j in seq_len(m)) if (p[[j + 1L]] != 0L) out[[p[[j + 1L]]]] <- j
  out
}

#' Fit once per seed
#'
#' For each seed in `seeds`, seeds R's stream with `set.seed(seed)` under
#' pinned kinds (Mersenne-Twister, Inversion, Rejection; recorded as the
#' result's `rng_kind` attribute) and calls `fit_fn()`, which must
#' return an `hdp_clustering_result` (e.g. a closure over
#' [hdp_clustering()]). A fit's `responsibilities` (records x `k_max`) are
#' dropped to bound memory; its `assignments` and `record_assignments` stay.
#' Each fit depends on its seed alone, so `cores > 1` (forked workers via
#' `parallel::mclapply`; not on Windows) returns the same fits as `cores = 1`.
#'
#' Returns a list of fits named by seed, in the order of `seeds`, with
#' attribute `rng_kind`. Raises on
#' a seed that is `NA`, not a whole number, or repeated; on a fit of the
#' wrong class; and on any worker's error, naming its seed.
#' @export
hdp_fit_seeds <- function(seeds, fit_fn, cores = 1L) {
  if (!is.numeric(seeds) || length(seeds) == 0L || anyNA(seeds) || any(seeds != round(seeds))) {
    stop("seeds must be a non-empty vector of whole numbers", call. = FALSE)
  }
  if (anyDuplicated(seeds)) stop("seeds must not repeat", call. = FALSE)
  if (!is.function(fit_fn)) stop("fit_fn must be a function", call. = FALSE)
  seeds <- as.integer(seeds)
  # Seeding under pinned kinds changes the session's kinds; give the
  # caller's back (and its stream position, if it had one).
  old_kind <- RNGkind()
  had_seed <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
  if (had_seed) old_seed <- get(".Random.seed", envir = globalenv(), inherits = FALSE)
  on.exit({
    do.call(RNGkind, as.list(old_kind))
    # nolint next: object_name_linter. .Random.seed is R's own name.
    if (had_seed) assign(".Random.seed", old_seed, envir = globalenv())
  }, add = TRUE)
  one <- function(seed) {
    set.seed(seed, kind = .HDP_RNG_KIND[["kind"]], normal.kind = .HDP_RNG_KIND[["normal.kind"]],
             sample.kind = .HDP_RNG_KIND[["sample.kind"]])
    fit <- fit_fn()
    if (!inherits(fit, "hdp_clustering_result")) {
      stop(sprintf("seed %d: fit_fn() did not return an hdp_clustering_result", seed),
           call. = FALSE)
    }
    fit[["hdp"]][["responsibilities"]] <- NULL
    fit
  }
  # An error becomes a value, so every seed is tried and every failure named
  # (and mclapply has no worker error to warn about); it is raised below.
  guarded <- function(seed) {
    tryCatch(one(seed), error = function(e) {
      structure(list(conditionMessage(e)), class = "hdp_seed_error")
    })
  }
  fits <- if (cores > 1L) {
    parallel::mclapply(seeds, guarded, mc.cores = cores, mc.preschedule = FALSE)
  } else {
    lapply(seeds, guarded)
  }
  # A worker killed outright comes back NULL or as a try-error.
  failed <- vapply(fits, function(f) {
    is.null(f) || inherits(f, "try-error") || inherits(f, "hdp_seed_error")
  }, logical(1))
  if (any(failed)) {
    msgs <- vapply(fits[failed], function(f) {
      if (is.null(f)) "no result" else as.character(f[[1L]])
    }, character(1))
    stop(sprintf("fit failed for seed(s) %s: %s", paste(seeds[failed], collapse = ", "),
                 paste(unique(trimws(msgs)), collapse = "; ")), call. = FALSE)
  }
  fits <- stats::setNames(fits, as.character(seeds))
  attr(fits, "rng_kind") <- .HDP_RNG_KIND
  fits
}

# The RNG kinds every per-seed fit is seeded under, named so a change of R's
# defaults cannot move a seed's stream.
# nolint next: object_name_linter. A constant.
.HDP_RNG_KIND <- c(kind = "Mersenne-Twister", normal.kind = "Inversion",
                   sample.kind = "Rejection")

#' The best fit by final ELBO
#'
#' The index of the lowest seed among the fits whose `elbo` is within
#' `rel_tol * |max(elbo)|` of the largest. The tolerance treats ELBOs that
#' differ only by floating-point summation order as tied, so the choice
#' cannot hinge on the last bits. Raises on a non-finite ELBO, on `NA` or
#' repeated seeds, and on lengths that disagree.
#' @export
hdp_best_elbo <- function(elbo, seeds, rel_tol = 1e-9) {
  if (length(elbo) != length(seeds) || length(elbo) == 0L) {
    stop("elbo and seeds must be non-empty and of equal length", call. = FALSE)
  }
  if (!is.numeric(elbo) || any(!is.finite(elbo))) stop("every ELBO must be finite", call. = FALSE)
  if (anyNA(seeds) || anyDuplicated(seeds)) stop("seeds must be distinct, not NA", call. = FALSE)
  if (!is.numeric(rel_tol) || length(rel_tol) != 1L || is.na(rel_tol) || rel_tol < 0) {
    stop("rel_tol must be one non-negative number", call. = FALSE)
  }
  best <- max(elbo)
  top <- which(elbo >= best - rel_tol * abs(best))
  top[[which.min(seeds[top])]]
}

#' Pairwise adjusted Rand index between fits' record assignments
#'
#' `fits` is a named list of `hdp_clustering_result`s on the same records
#' (as [hdp_fit_seeds()] returns). Raises unless every fit's
#' `record_assignments` holds the same `id` and `group` in the same order,
#' since assignments are compared position by position.
#'
#' Returns a data frame with one row per unordered pair: `seed_a`, `seed_b`
#' (the list's names), `ari`.
#' @export
hdp_pairwise_ari <- function(fits) {
  if (length(fits) < 2L || is.null(names(fits))) {
    stop("fits must be a named list of at least two fits", call. = FALSE)
  }
  ref <- fits[[1L]][["record_assignments"]]
  for (nm in names(fits)) {
    ra <- fits[[nm]][["record_assignments"]]
    if (!identical(ra[["id"]], ref[["id"]]) || !identical(ra[["group"]], ref[["group"]])) {
      stop(sprintf("fit %s is not on the same records, in the same order, as fit %s",
                   nm, names(fits)[[1L]]), call. = FALSE)
    }
  }
  pairs <- utils::combn(length(fits), 2L)
  ari <- vapply(seq_len(ncol(pairs)), function(p) {
    adjusted_rand_index(fits[[pairs[1L, p]]][["record_assignments"]][["assignment"]],
                        fits[[pairs[2L, p]]][["record_assignments"]][["assignment"]])
  }, numeric(1))
  data.frame(seed_a = names(fits)[pairs[1L, ]], seed_b = names(fits)[pairs[2L, ]], ari = ari)
}

# Cosine similarity between every row of `a` and every row of `b`.
.cosine_matrix <- function(a, b) {
  na <- sqrt(rowSums(a * a))
  nb <- sqrt(rowSums(b * b))
  if (any(na == 0) || any(nb == 0)) stop("a profile is all zeros", call. = FALSE)
  (a %*% t(b)) / outer(na, nb)
}

#' Match one fit's clusters to a reference fit's
#'
#' The optimal one-to-one matching of clusters `ref_clusters` of `ref` to
#' clusters `clusters` of `fit` (both `hdp_clustering_result`s over the same
#' items in the same order), minimizing the total of 1 - cosine similarity
#' between their item-probability profiles (`theta` rows), by
#' [solve_lsap()]. Cosine rather than an L1 distance, because profiles are
#' sparse: two clusters with different dominant items differ in a handful
#' of entries, which an L1 distance averaged over items dilutes toward zero.
#'
#' Returns a data frame, one row per reference cluster: `ref_cluster`,
#' `cluster` and `cosine` (both `NA` when `fit` has fewer clusters to offer
#' than the reference).
#' @export
hdp_match_clusters <- function(ref, fit, ref_clusters, clusters) {
  if (!identical(ref[["item_names"]], fit[["item_names"]])) {
    stop("the two fits are not over the same items in the same order", call. = FALSE)
  }
  out <- data.frame(ref_cluster = as.integer(ref_clusters),
                    cluster = rep(NA_integer_, length(ref_clusters)),
                    cosine = rep(NA_real_, length(ref_clusters)))
  if (length(ref_clusters) == 0L || length(clusters) == 0L) return(out)
  sim <- .cosine_matrix(ref[["hdp"]][["theta"]][ref_clusters, , drop = FALSE],
                        fit[["hdp"]][["theta"]][clusters, , drop = FALSE])
  col <- solve_lsap(1 - sim)
  hit <- !is.na(col)
  out[["cluster"]][hit] <- as.integer(clusters[col[hit]])
  out[["cosine"]][hit] <- sim[cbind(which(hit), col[hit])]
  out
}

#' How often each of a reference fit's clusters recurs in the other fits
#'
#' For every other fit in `fits`, matches its active clusters to the
#' reference fit's active clusters ([hdp_match_clusters()]). A reference
#' cluster *recurs* in that fit when its match has cosine similarity of at
#' least `min_cosine`; its category *agrees* when it recurs and the matched
#' cluster has the same category.
#'
#' - `fits`: a named list of `hdp_clustering_result`s; `reference`: the name
#'   (or index) of the reference fit among them.
#' - `active`: a list parallel to `fits`, each fit's active cluster indices;
#'   `NULL` uses [hdp_active_clusters()] at `beta_threshold`.
#' - `categories`: a list parallel to `fits`, each a character vector with
#'   one category per cluster (length `k_max`); `NULL` uses
#'   [hdp_cluster_categorization()]'s `category` at `beta_threshold`.
#'
#' Returns a data frame, one row per active reference cluster: `cluster`,
#' `category`, `n_other` (fits compared), `n_recur`, `recurrence_share`,
#' `n_category_agree`, `category_agreement_share` (both shares out of all
#' `n_other` fits compared, not out of the fits in which it recurs), and
#' `median_cosine` of its match over the fits in which it was matched at
#' all (`NA` if none).
#' @export
hdp_cluster_recurrence <- function(fits, reference, active = NULL, categories = NULL,
                                   beta_threshold = 0.01, min_cosine = 0.9) {
  if (is.null(names(fits)) || length(fits) < 2L) {
    stop("fits must be a named list of at least two fits", call. = FALSE)
  }
  ref_i <- if (is.character(reference)) match(reference, names(fits)) else as.integer(reference)
  if (length(ref_i) != 1L || is.na(ref_i) || ref_i < 1L || ref_i > length(fits)) {
    stop("reference must name or index one fit in fits", call. = FALSE)
  }
  if (!is.numeric(min_cosine) || length(min_cosine) != 1L || is.na(min_cosine) ||
        min_cosine < -1 || min_cosine > 1) {
    stop("min_cosine must be one number in [-1, 1]", call. = FALSE)
  }
  if (is.null(active)) {
    active <- lapply(fits, hdp_active_clusters, beta_threshold = beta_threshold)
  }
  if (is.null(categories)) {
    categories <- lapply(fits, function(f) {
      hdp_cluster_categorization(f, beta_threshold = beta_threshold)[["category"]]
    })
  }
  if (length(active) != length(fits) || length(categories) != length(fits)) {
    stop("active and categories must be parallel to fits", call. = FALSE)
  }
  for (i in seq_along(fits)) {
    if (length(categories[[i]]) != fits[[i]][["hdp"]][["k_max"]] || anyNA(categories[[i]])) {
      stop(sprintf("categories for fit %s must hold one category per cluster",
                   names(fits)[[i]]), call. = FALSE)
    }
  }
  ref <- fits[[ref_i]]
  ref_clusters <- as.integer(active[[ref_i]])
  ref_cats <- categories[[ref_i]][ref_clusters]
  others <- setdiff(seq_along(fits), ref_i)
  k <- length(ref_clusters)
  recur <- matrix(FALSE, k, length(others))
  agree <- matrix(FALSE, k, length(others))
  cosine <- matrix(NA_real_, k, length(others))
  for (o in seq_along(others)) {
    i <- others[[o]]
    m <- hdp_match_clusters(ref, fits[[i]], ref_clusters, as.integer(active[[i]]))
    cosine[, o] <- m[["cosine"]]
    recur[, o] <- !is.na(m[["cosine"]]) & m[["cosine"]] >= min_cosine
    same <- categories[[i]][m[["cluster"]]] == ref_cats
    agree[, o] <- recur[, o] & !is.na(same) & same
  }
  n_other <- length(others)
  med <- vapply(seq_len(k), function(r) {
    x <- cosine[r, ]
    if (all(is.na(x))) NA_real_ else stats::median(x, na.rm = TRUE)
  }, numeric(1))
  data.frame(cluster = ref_clusters, category = ref_cats,
             n_other = rep(n_other, k), n_recur = as.integer(rowSums(recur)),
             recurrence_share = rowSums(recur) / n_other,
             n_category_agree = as.integer(rowSums(agree)),
             category_agreement_share = rowSums(agree) / n_other,
             median_cosine = med)
}
