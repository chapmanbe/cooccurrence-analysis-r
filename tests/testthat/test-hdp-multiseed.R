# Multi-seed HDP: the seed runner, best-ELBO selection, pairwise ARI, the
# LSAP solver and cluster recurrence. R-only except adjusted_rand_index(),
# whose expected values come from running the reference's
# `_adjusted_rand_index` (network_group_stratification.jl).

# Every permutation of 1..n, one per row (brute force for small n).
all_perms <- function(n) {
  if (n == 1L) return(matrix(1L, 1L, 1L))
  rest <- all_perms(n - 1L)
  do.call(rbind, lapply(seq_len(n), function(first) {
    cbind(first, matrix(setdiff(seq_len(n), first)[rest], nrow(rest)))
  }))
}

# The minimum total cost over every injective row -> column map (rows <= cols).
brute_force_lsap <- function(cost) {
  n <- nrow(cost)
  m <- ncol(cost)
  best <- Inf
  for (cols in utils::combn(m, n, simplify = FALSE)) {
    perms <- all_perms(n)
    for (r in seq_len(nrow(perms))) {
      best <- min(best, sum(cost[cbind(seq_len(n), cols[perms[r, ]])]))
    }
  }
  best
}

# 300 records in three item blocks, groups A (first 150) and B, each block
# item present with probability 0.9: a smaller copy of test-hdp.R's
# blocky_events(), which other test files cannot see.
seed_events <- function() {
  blocks <- list(c("a1", "a2", "a3"), c("b1", "b2", "b3"), c("c1", "c2", "c3"))
  parts <- lapply(1:300, function(pid) {
    blk <- blocks[[(pid %% 3L) + 1L]]
    keep <- stats::runif(3L) < 0.9
    data.frame(id = rep(pid, sum(keep)), item = blk[keep],
               Group = rep(if (pid <= 150L) "A" else "B", sum(keep)))
  })
  do.call(rbind, parts)
}

# A stand-in fit carrying only what the matching reads: items and theta.
fake_fit <- function(theta, items = paste0("i", seq_len(ncol(theta)))) {
  structure(list(hdp = list(theta = theta, k_max = nrow(theta)), item_names = items),
            class = "hdp_clustering_result")
}

test_that("adjusted_rand_index matches the Julia reference", {
  # Values printed by CooccurrenceAnalysis._adjusted_rand_index (Julia).
  expect_equal(adjusted_rand_index(c(1, 1, 2, 2, 3, 3), c(1, 1, 2, 3, 3, 3)),
               0.4444444444444444, tolerance = 1e-15)
  expect_equal(adjusted_rand_index(c(1, 1, 1, 2, 2, 2, 3, 3), c(2, 2, 1, 1, 1, 3, 3, 3)),
               0.23809523809523808, tolerance = 1e-15)
  expect_identical(adjusted_rand_index(1:4, c(1, 1, 1, 1)), 0)
  expect_identical(adjusted_rand_index(c(1, 1, 1), c(2, 2, 2)), 1)
  expect_equal(adjusted_rand_index(c(1, 1, 2, 2, 1, 2, 3, 1, 2, 3),
                                   c(3, 3, 1, 1, 3, 2, 2, 3, 1, 2)),
               0.7232472324723247, tolerance = 1e-15)
})

test_that("adjusted_rand_index: R-only edge cases", {
  a <- c(1, 1, 2, 2, 3, 3, 3)
  expect_identical(adjusted_rand_index(a, c(9, 9, 4, 4, 7, 7, 7)), 1)   # relabeling
  expect_identical(adjusted_rand_index(c("x", "x", "y"), c(2, 2, 5)), 1)
  expect_identical(adjusted_rand_index(1, 1), 0)
  expect_error(adjusted_rand_index(1:3, 1:2), "same length")
  expect_error(adjusted_rand_index(c(1, NA), c(1, 2)), "NA")
})

test_that("solve_lsap: known answers by brute force over every permutation", {
  set.seed(11)
  for (n in 1:6) {
    for (rep in 1:5) {
      # Small integer costs make ties common, which is where solvers go wrong.
      cost <- matrix(if (rep %% 2L == 0L) sample(0:3, n * n, TRUE) else stats::runif(n * n),
                     n, n)
      col <- solve_lsap(cost)
      expect_identical(sort(col), seq_len(n))
      expect_equal(sum(cost[cbind(seq_len(n), col)]), brute_force_lsap(cost), tolerance = 1e-12)
    }
  }
})

test_that("solve_lsap: known answers where greedy matching is not optimal", {
  # Greedy takes the cheapest cell (1,1) = 0.1 and is left with 0.9: total
  # 1.0. The optimum crosses over: 0.2 + 0.15 = 0.35.
  square <- matrix(c(0.1, 0.15, 0.2, 0.9), 2, 2)
  expect_identical(solve_lsap(square), c(2L, 1L))
  # The same trap with a third, worse column on offer.
  wide <- matrix(c(0.1, 0.15, 0.2, 0.9, 0.95, 0.99), 2, 3)
  expect_identical(solve_lsap(wide), c(2L, 1L))
  # And transposed: three rows, two columns; the third row goes unmatched.
  expect_identical(solve_lsap(t(wide)), c(2L, 1L, NA))
  # Three by three: greedy takes (1,1) = 0, then (2,2) = 5, then (3,3) = 9,
  # total 14; the optimum is the anti-diagonal, 1 + 5 + 1 = 7.
  three <- matrix(c(0, 1, 1,
                    1, 5, 9,
                    1, 9, 9), 3, 3, byrow = TRUE)
  expect_identical(solve_lsap(three), c(3L, 2L, 1L))
})

test_that("solve_lsap: rectangular costs, either way round", {
  set.seed(12)
  for (rep in 1:5) {
    wide <- matrix(stats::runif(3 * 5), 3, 5)
    col <- solve_lsap(wide)
    expect_false(anyDuplicated(col) > 0L)
    expect_equal(sum(wide[cbind(1:3, col)]), brute_force_lsap(wide), tolerance = 1e-12)
    tall <- t(wide)
    col_t <- solve_lsap(tall)
    expect_identical(sum(is.na(col_t)), 2L)
    hit <- which(!is.na(col_t))
    expect_identical(sort(col_t[hit]), 1:3)
    expect_equal(sum(tall[cbind(hit, col_t[hit])]), brute_force_lsap(wide), tolerance = 1e-12)
  }
  expect_identical(solve_lsap(matrix(numeric(0), 0, 3)), integer(0))
  expect_error(solve_lsap(matrix(c(1, Inf, 2, 3), 2)), "finite")
})

test_that("hdp_best_elbo picks the maximum and breaks ties toward the lowest seed", {
  expect_identical(hdp_best_elbo(c(-10, -3, -7), c(2026, 2027, 2028)), 2L)
  # A tie listed out of seed order: the lower seed wins, not the first listed.
  expect_identical(hdp_best_elbo(c(-5, -5, -7), c(2030, 2026, 2028)), 2L)
  expect_identical(hdp_best_elbo(c(-1, -2, -1), c(2026, 2027, 2028)), 1L)
  expect_error(hdp_best_elbo(c(-1, NaN), c(1, 2)), "finite")
  expect_error(hdp_best_elbo(c(-1, -2), c(1, 1)), "distinct")
  expect_error(hdp_best_elbo(c(-1, -2), 1), "equal length")
})

test_that("hdp_best_elbo treats ELBOs within the relative tolerance as tied", {
  # 5e-10 relative below the maximum: inside the default 1e-9, so a tie,
  # and the lower seed (2026, listed second) wins over the larger ELBO.
  near <- c(-1e6, -1e6 * (1 + 5e-10), -1.1e6)
  expect_identical(hdp_best_elbo(near, c(2030, 2026, 2028)), 2L)
  # 1e-8 relative below: outside the tolerance, so the maximum wins.
  far <- c(-1e6, -1e6 * (1 + 1e-8), -1.1e6)
  expect_identical(hdp_best_elbo(far, c(2030, 2026, 2028)), 1L)
  # rel_tol = 0 is the exact rule.
  expect_identical(hdp_best_elbo(near, c(2030, 2026, 2028), rel_tol = 0), 1L)
  expect_error(hdp_best_elbo(near, c(1, 2, 3), rel_tol = -1), "rel_tol")
})

test_that("hdp_fit_seeds: one fit per seed, reproducible, the same in parallel", {
  set.seed(5)
  ev <- seed_events()
  fit_fn <- function() {
    hdp_clustering(ev, group_col = "Group", k_max = 5L, n_init = 1L, max_iter = 40L)
  }
  fits <- hdp_fit_seeds(c(7, 3, 9), fit_fn)
  expect_identical(names(fits), c("7", "3", "9"))
  expect_true(all(vapply(fits, inherits, logical(1), "hdp_clustering_result")))
  expect_false("responsibilities" %in% names(fits[[1L]][["hdp"]]))
  expect_identical(length(fits[[1L]][["hdp"]][["assignments"]]),
                   nrow(fits[[1L]][["record_assignments"]]))
  # A fit depends on its seed alone, not on what ran before it.
  again <- hdp_fit_seeds(c(9, 7), fit_fn)
  expect_identical(again[["7"]], fits[["7"]])
  expect_identical(again[["9"]], fits[["9"]])
  if (.Platform[["OS.type"]] == "unix") {
    expect_identical(hdp_fit_seeds(c(7, 3, 9), fit_fn, cores = 2L), fits)
  }
})

test_that("hdp_fit_seeds pins the RNG kinds, whatever the caller's", {
  draw <- function() {
    structure(list(hdp = list(u = stats::runif(1), n = stats::rnorm(1))),
              class = "hdp_clustering_result")
  }
  old <- RNGkind()
  on.exit(do.call(RNGkind, as.list(old)), add = TRUE)
  RNGkind("Mersenne-Twister", "Inversion", "Rejection")
  pinned <- hdp_fit_seeds(c(1, 2), draw)
  RNGkind("L'Ecuyer-CMRG", "Box-Muller")
  expect_identical(hdp_fit_seeds(c(1, 2), draw), pinned)
  # The caller's kinds survive the call.
  expect_identical(RNGkind()[1:2], c("L'Ecuyer-CMRG", "Box-Muller"))
  expect_identical(attr(pinned, "rng_kind"),
                   c(kind = "Mersenne-Twister", normal.kind = "Inversion",
                     sample.kind = "Rejection"))
})

test_that("hdp_fit_seeds raises rather than dropping a seed", {
  expect_error(hdp_fit_seeds(c(1, 1), function() NULL), "repeat")
  expect_error(hdp_fit_seeds(c(1, NA), function() NULL), "whole numbers")
  expect_error(hdp_fit_seeds(1.5, function() NULL), "whole numbers")
  expect_error(hdp_fit_seeds(4, function() list()), "seed 4")
  bad <- function() if (stats::runif(1) >= 0) stop("boom")
  expect_error(hdp_fit_seeds(c(1, 2), bad), "seed\\(s\\) 1, 2: boom")
  if (.Platform[["OS.type"]] == "unix") {
    expect_error(hdp_fit_seeds(c(1, 2), bad, cores = 2L), "seed\\(s\\) 1, 2: boom")
  }
})

test_that("hdp_pairwise_ari: every unordered pair, records aligned", {
  set.seed(5)
  ev <- seed_events()
  fits <- hdp_fit_seeds(c(1, 2, 3), function() {
    hdp_clustering(ev, group_col = "Group", k_max = 5L, n_init = 1L, max_iter = 40L)
  })
  ari <- hdp_pairwise_ari(fits)
  expect_identical(ari[["seed_a"]], c("1", "1", "2"))
  expect_identical(ari[["seed_b"]], c("2", "3", "3"))
  expect_true(all(ari[["ari"]] <= 1 + 1e-12))
  a1 <- fits[[1L]][["record_assignments"]][["assignment"]]
  a2 <- fits[[2L]][["record_assignments"]][["assignment"]]
  expect_identical(ari[["ari"]][[1L]], adjusted_rand_index(a1, a2))
  same <- list(a = fits[[1L]], b = fits[[1L]])
  expect_identical(hdp_pairwise_ari(same)[["ari"]], 1)
  shuffled <- fits[[2L]]
  ra <- shuffled[["record_assignments"]]
  shuffled[["record_assignments"]] <- ra[rev(seq_len(nrow(ra))), ]
  expect_error(hdp_pairwise_ari(list(a = fits[[1L]], b = shuffled)), "same records")
  expect_error(hdp_pairwise_ari(fits[1L]), "at least two")
})

test_that("hdp_match_clusters recovers a permutation of the same profiles", {
  theta <- rbind(c(0.9, 0.8, 0.0, 0.0, 0.0),
                 c(0.0, 0.0, 0.9, 0.7, 0.0),
                 c(0.1, 0.0, 0.0, 0.0, 0.9))
  ref <- fake_fit(theta)
  perm <- c(3L, 1L, 2L)
  other <- fake_fit(theta[perm, ])
  m <- hdp_match_clusters(ref, other, 1:3, 1:3)
  expect_identical(m[["cluster"]], match(1:3, perm))
  expect_equal(m[["cosine"]], c(1, 1, 1))
  # Fewer clusters on offer: the worst-fitting reference cluster goes unmatched.
  m2 <- hdp_match_clusters(ref, other, 1:3, c(2L, 3L))
  expect_identical(m2[["cluster"]], c(2L, 3L, NA))
  expect_error(hdp_match_clusters(ref, fake_fit(theta, letters[1:5]), 1:3, 1:3), "same items")
})

test_that("hdp_match_clusters is optimal where greedy matching is not", {
  # Cosines: ref 1 ~ other 1 is 0.981, the best single pair, but taking it
  # leaves ref 2 ~ other 2 at 0. The optimum crosses: 0.703 + 0.832.
  ref <- fake_fit(rbind(c(0.9, 0.9, 0.0), c(0.9, 0.0, 0.0)))
  other <- fake_fit(rbind(c(0.9, 0.6, 0.0), c(0.0, 0.9, 0.1)))
  m <- hdp_match_clusters(ref, other, 1:2, 1:2)
  expect_identical(m[["cluster"]], c(2L, 1L))
  expect_equal(m[["cosine"]], c(0.81 / sqrt(1.62 * 0.82), 0.81 / sqrt(0.81 * 1.17)),
               tolerance = 1e-12)
})

test_that("hdp_cluster_recurrence counts recurrence and category agreement", {
  a <- c(0.9, 0.8, 0.0, 0.0, 0.0)
  b <- c(0.0, 0.0, 0.9, 0.7, 0.0)
  cc <- c(0.0, 0.0, 0.0, 0.1, 0.9)
  far <- c(0.5, 0.0, 0.5, 0.0, 0.5)
  fits <- list(
    ref = fake_fit(rbind(a, b, cc, far)),       # cluster 4 inactive
    s2 = fake_fit(rbind(b, a, cc)),             # all three recur; cc's category differs
    s3 = fake_fit(rbind(a, far, c(0.3, 0.3, 0.3, 0.3, 0.3)))  # only a recurs
  )
  active <- list(1:3, 1:3, 1:3)
  categories <- list(c("u", "x", "y", "negligible"), c("x", "u", "u"), c("x", "y", "y"))
  rec <- hdp_cluster_recurrence(fits, "ref", active = active, categories = categories,
                                min_cosine = 0.9)
  expect_identical(rec[["cluster"]], 1:3)
  expect_identical(rec[["category"]], c("u", "x", "y"))
  expect_identical(rec[["n_other"]], c(2L, 2L, 2L))
  expect_identical(rec[["n_recur"]], c(2L, 1L, 1L))
  expect_identical(rec[["recurrence_share"]], c(1, 0.5, 0.5))
  # a: u in s2 (agrees), x in s3 (disagrees); b: x in s2 (agrees); cc: u in s2.
  expect_identical(rec[["n_category_agree"]], c(1L, 1L, 0L))
  expect_identical(rec[["category_agreement_share"]], c(0.5, 0.5, 0))
  expect_equal(rec[["median_cosine"]][[1L]], 1)
  # The threshold is what decides recurrence.
  lax <- hdp_cluster_recurrence(fits, 1L, active = active, categories = categories,
                                min_cosine = -1)
  expect_identical(lax[["n_recur"]], c(2L, 2L, 2L))
  expect_error(hdp_cluster_recurrence(fits, "nope", active = active, categories = categories),
               "reference")
  expect_error(hdp_cluster_recurrence(fits, 1L, active = active,
                                      categories = list("u", "x", "y")), "one category")
})
