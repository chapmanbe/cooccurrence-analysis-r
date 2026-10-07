# HDP Bernoulli mixture. One case per Julia @testset in the reference's
# test/runtests.jl ("empirical_bayes_priors", "HDP Bernoulli mixture",
# "HDP item-domain wrappers", "HDP native 3-group categorization"), named
# after it, plus R-only cases marked as such. Julia's StableRNG seeds do not
# carry over; each fixture here is re-drawn from R's stream.

# Planted 3-cluster data, as the reference's "3-cluster planted recovery":
# A (items 1-2) in both groups, B (items 3-4) group 1 only, C (items 5-6)
# group 2 only.
planted_matrix <- function(n_per = 130L, d = 8L) {
  x <- matrix(FALSE, 3L * n_per, d)
  truth <- rep(1:3, each = n_per)
  grp <- c(ifelse(seq_len(n_per) <= n_per %/% 2L, 1L, 2L), rep(1L, n_per), rep(2L, n_per))
  a <- seq_len(n_per)
  b <- n_per + a
  cc <- 2L * n_per + a
  x[a, 1] <- stats::runif(n_per) < 0.85
  x[a, 2] <- stats::runif(n_per) < 0.80
  x[b, 3] <- stats::runif(n_per) < 0.85
  x[b, 4] <- stats::runif(n_per) < 0.80
  x[cc, 5] <- stats::runif(n_per) < 0.85
  x[cc, 6] <- stats::runif(n_per) < 0.80
  list(x = x, truth = truth, group = grp)
}

# The reference's C9-regression fixture: 1200 records in three item blocks,
# groups A (first 600) and B, each block item present with probability 0.9.
blocky_events <- function() {
  blocks <- list(c("a1", "a2", "a3"), c("b1", "b2", "b3"), c("c1", "c2", "c3"))
  parts <- lapply(1:1200, function(pid) {
    blk <- blocks[[(pid %% 3L) + 1L]]
    keep <- stats::runif(3L) < 0.9
    data.frame(id = rep(pid, sum(keep)), item = blk[keep],
               Group = rep(if (pid <= 600L) "A" else "B", sum(keep)))
  })
  do.call(rbind, parts)
}

# The reference's "planted labels" fixture: S1/S2 in 40 records of each
# group (ids 1-80), A1/A2 in 60 group-A records, B1/B2 in 60 group-B ones.
shared_and_private_events <- function() {
  rec <- function(items, n) list(items = items, n = n)
  ev <- events_from_records(list(rec(c("S1", "S2"), 40L), rec(c("S1", "S2"), 40L),
                                 rec(c("A1", "A2"), 60L), rec(c("B1", "B2"), 60L)))
  ev[["Group"]] <- rep(c("A", "B", "A", "B"), times = c(80L, 80L, 120L, 120L))
  ev
}

test_that("empirical_bayes_priors", {
  x <- matrix(FALSE, 100L, 2L)
  x[1:50, 1] <- TRUE
  x[1:10, 2] <- TRUE
  pr <- empirical_bayes_priors(x, concentration = 10, floor = 1)
  expect_length(pr[["alpha"]], 2L)
  expect_length(pr[["beta"]], 2L)
  expect_equal(unname(pr[["alpha"]]), c(6, 2))
  expect_equal(unname(pr[["beta"]]), c(6, 10))
  zero <- empirical_bayes_priors(matrix(FALSE, 50L, 1L))
  expect_equal(unname(zero[["alpha"]]), 1)
  expect_equal(unname(zero[["beta"]]), 11)
  # R-only: no rows raises rather than dividing by zero.
  expect_error(empirical_bayes_priors(matrix(FALSE, 0L, 2L)), "empty data")
})

test_that("3-cluster planted recovery", {
  set.seed(2026)
  pm <- planted_matrix()
  set.seed(99)
  res <- fit_hdp_bernoulli_mixture(pm[["x"]], pm[["group"]], 2L, k_max = 8L, prior = "flat",
                                   n_init = 6L)
  expect_s3_class(res, "hdp_bernoulli_result")
  expect_identical(res[["k_max"]], 8L)
  expect_identical(res[["n_groups"]], 2L)
  expect_identical(dim(res[["theta"]]), c(8L, 8L))
  expect_identical(dim(res[["responsibilities"]]), c(390L, 8L))
  expect_length(res[["assignments"]], 390L)
  expect_length(res[["group"]], 390L)
  expect_identical(dim(res[["pi_mean"]]), c(2L, 8L))
  expect_length(res[["beta_mean"]], 8L)

  expect_true(all(res[["theta"]] >= 0 & res[["theta"]] <= 1))
  expect_equal(sum(res[["beta_mean"]]), 1, tolerance = 0.05)
  expect_true(all(abs(rowSums(res[["pi_mean"]]) - 1) <= 0.1))
  expect_true(all(abs(rowSums(res[["responsibilities"]]) - 1) <= 1e-6))
  expect_gt(adjusted_rand(res[["assignments"]], pm[["truth"]]), 0.8)
  expect_true(res[["effective_k"]] %in% 2:5)
  expect_true(is.finite(res[["elbo"]]))
})

test_that("K_max truncation insensitivity", {
  set.seed(42)
  x <- matrix(FALSE, 180L, 6L)
  grp <- c(rep(1L, 90L), rep(2L, 90L))
  x[1:60, 1] <- stats::runif(60) < 0.8
  x[1:60, 2] <- stats::runif(60) < 0.7
  x[61:120, 3] <- stats::runif(60) < 0.8
  x[61:120, 4] <- stats::runif(60) < 0.7
  x[121:180, 5] <- stats::runif(60) < 0.8
  x[121:180, 6] <- stats::runif(60) < 0.7
  set.seed(1)
  r6 <- fit_hdp_bernoulli_mixture(x, grp, 2L, k_max = 6L, n_init = 2L)
  set.seed(1)
  r12 <- fit_hdp_bernoulli_mixture(x, grp, 2L, k_max = 12L, n_init = 2L)
  expect_lte(abs(r6[["effective_k"]] - r12[["effective_k"]]), 2L)
  expect_gte(r6[["effective_k"]], 2L)
  expect_gte(r12[["effective_k"]], 2L)
})

test_that("empirical_bayes prior", {
  set.seed(7)
  x <- matrix(FALSE, 120L, 8L)
  x[1:60, 1] <- stats::runif(60) < 0.8
  x[61:120, 2] <- stats::runif(60) < 0.8
  res <- fit_hdp_bernoulli_mixture(x, rep(1:2, each = 60L), 2L, k_max = 6L,
                                   prior = "empirical_bayes", concentration = 5, n_init = 2L)
  expect_s3_class(res, "hdp_bernoulli_result")
  expect_true(all(res[["theta"]] >= 0 & res[["theta"]] <= 1))
  # R-only: an unknown prior raises instead of falling back to flat.
  expect_error(fit_hdp_bernoulli_mixture(x, rep(1:2, each = 60L), 2L, prior = "jeffreys"),
               "unsupported prior")
})

test_that("error on empty group", {
  x <- matrix(FALSE, 20L, 4L)
  x[1:10, 1] <- TRUE
  expect_error(fit_hdp_bernoulli_mixture(x, rep(1L, 20L), 2L, k_max = 4L), "empty group")
})

test_that("single group reduces gracefully", {
  set.seed(3)
  x <- matrix(FALSE, 80L, 4L)
  x[1:40, 1] <- stats::runif(40) < 0.8
  x[41:80, 2] <- stats::runif(40) < 0.8
  res <- fit_hdp_bernoulli_mixture(x, rep(1L, 80L), 1L, k_max = 6L, n_init = 2L)
  expect_s3_class(res, "hdp_bernoulli_result")
  expect_identical(res[["n_groups"]], 1L)
  expect_identical(dim(res[["pi_mean"]]), c(1L, 6L))
})

test_that("N < K_max warns and runs", {
  x <- matrix(FALSE, 5L, 3L)
  x[1:2, 1] <- TRUE
  set.seed(1)
  expect_warning(fit_hdp_bernoulli_mixture(x, c(1L, 1L, 2L, 2L, 2L), 2L, k_max = 8L,
                                           n_init = 1L),
                 "N=5 < k_max=8")
})

test_that("hdp_clustering on event_df", {
  set.seed(42)
  res <- hdp_clustering(make_test_event_df(), group_col = "Group", k_max = 4L, n_init = 2L)
  expect_s3_class(res, "hdp_clustering_result")
  expect_s3_class(res[["hdp"]], "hdp_bernoulli_result")
  expect_identical(res[["hdp"]][["n_groups"]], 2L)
  expect_gt(length(res[["item_names"]]), 0L)
  expect_length(res[["n_records_per_group"]], 2L)
  expect_true(all(res[["n_records_per_group"]] > 0L))
  expect_identical(sum(res[["n_records_per_group"]]), nrow(res[["record_assignments"]]))
  cp <- res[["class_profiles"]]
  expect_identical(nrow(cp), 4L)
  expect_true(all(c("cluster", "beta_weight") %in% names(cp)))
  expect_equal(sum(cp[["beta_weight"]]), 1, tolerance = 0.05)
  gp <- res[["group_profiles"]]
  expect_identical(nrow(gp), 4L)
  expect_true("cluster" %in% names(gp))
  expect_identical(names(gp)[-1L], c("pi_A", "pi_B"))
  ra <- res[["record_assignments"]]
  expect_true(all(c("id", "group", "assignment") %in% names(ra)))
  expect_true(all(ra[["assignment"]] >= 1L & ra[["assignment"]] <= 4L))
  # R-only, from the fixture's arithmetic: multi-item records are 10 + 3 + 1
  # in B and 5 + 1 in A; named by group, not position.
  expect_identical(res[["n_records_per_group"]], c(A = 6L, B = 14L))
})

test_that("hdp_cluster_categorization structure", {
  set.seed(77)
  res <- hdp_clustering(make_test_event_df(), k_max = 4L, n_init = 2L)
  cat_df <- hdp_cluster_categorization(res)
  expect_s3_class(cat_df, "data.frame")
  expect_identical(names(cat_df), c("cluster", "beta_weight", "pi_A", "pi_B", "category"))
  expect_identical(nrow(cat_df), 4L)
  valid <- c("universal", "negligible", "both_present_but_unequal", "group_a_only",
             "group_b_only")
  expect_true(all(cat_df[["category"]] %in% valid))
})

test_that("hdp_cluster_categorization planted labels (C3, T2)", {
  ev <- shared_and_private_events()
  set.seed(2026)
  res <- hdp_clustering(ev, k_max = 6L, n_init = 4L)
  cats <- hdp_cluster_categorization(res)[["category"]]
  expect_true("universal" %in% cats)
  expect_true("group_a_only" %in% cats)
  expect_true("group_b_only" %in% cats)
  expect_false(any(c("a_only", "b_only") %in% cats))
})

test_that("timing_filter propagates", {
  set.seed(5)
  res <- hdp_clustering(make_test_event_df(), k_max = 4L, timing_filter = "all", n_init = 1L)
  expect_s3_class(res, "hdp_clustering_result")
  # R-only: the reference's other timing filters are not ported, and say so.
  expect_error(hdp_clustering(make_test_event_df(), timing_filter = "concurrent"),
               "not ported")
})

test_that("error on unknown group_by column", {
  expect_error(hdp_clustering(make_test_event_df(), group_col = "NoSuchColumn"),
               "NoSuchColumn not found")
})

test_that("CAVI runs past the plateau window and separates blocks (C9 regression)", {
  set.seed(20260729)
  blocky <- blocky_events()
  set.seed(7)
  res <- hdp_clustering(blocky, k_max = 8L, n_init = 1L, max_iter = 300L, tol = 1e-5)
  h <- res[["hdp"]]
  expect_gt(h[["n_iter"]], 10L)
  active <- which(h[["beta_mean"]] > 0.01)
  expect_gte(length(active), 3L)
  dominant <- unique(vapply(active, function(k) res[["item_names"]][[which.max(h[["theta"]][k, ])]],
                            character(1)))
  expect_gte(length(dominant), 3L)
  expect_gt(max(h[["theta"]][active, ]), 0.6)
})

test_that("HDP native 3-group categorization", {
  private <- lapply(c("A", "B", "C"), function(g) list(items = paste0(g, 1:2), n = 14L))
  shared <- rep(list(list(items = c("S1", "S2"), n = 8L)), 3L)
  leak <- list(list(items = c("A1", "B1"), n = 1L))
  ev <- events_from_records(c(private, shared, leak))
  ev[["Group"]] <- rep(c("A", "B", "C", "A", "B", "C", "A"),
                       times = c(28L, 28L, 28L, 16L, 16L, 16L, 2L))
  set.seed(2026)
  res <- hdp_clustering(ev, k_max = 6L, n_init = 3L)
  expect_identical(res[["hdp"]][["n_groups"]], 3L)
  cats <- hdp_cluster_categorization(res)[["category"]]
  expect_true(any(c("universal", "both_present_but_unequal") %in% cats))
  expect_true(any(c("group_a_only", "group_b_only", "group_c_only") %in% cats))
  expect_false(any(c("a_only", "b_only", "c_only") %in% cats))
})

# ── R-only cases ─────────────────────────────────────────────────────────────

test_that("R-only: the stopping rule needs five consecutive small ELBO steps", {
  set.seed(20260729)
  blocky <- blocky_events()
  set.seed(7)
  h <- hdp_clustering(blocky, k_max = 8L, n_init = 1L)[["hdp"]]
  expect_true(h[["converged"]])
  tr <- h[["elbo_trace"]]
  expect_length(tr, h[["n_iter"]])
  small <- abs(diff(tr)) < 1e-5 * (1 + abs(utils::head(tr, -1L)))
  # The last five steps are small, and no earlier run of five is.
  expect_true(all(utils::tail(small, 5L)))
  runs <- stats::filter(as.numeric(small), rep(1, 5L), sides = 1L)
  expect_identical(which(runs == 5)[[1L]], length(small))
  # Stopping on "no improvement over the best ELBO" would have ended where
  # the ELBO first dipped below an earlier peak; that point lies before the
  # plateau on this fixture, as the regression above requires.
  expect_identical(h[["elbo"]], tr[[length(tr)]])
})

test_that("R-only: the stopping rule is on per-iteration change, not on the best ELBO", {
  # A trace shaped like real data's: an early peak (-800), a dip, then a
  # steady climb in steps of 0.1 -- well above tol * (1 + |ELBO|) = 0.0085 --
  # that stays below the early peak for many iterations.
  climb <- -850 + 0.1 * (1:40)
  trace <- c(-1000, -900, -800, -850, climb)
  for (n in 5:length(trace)) {
    expect_false(.hdp_converged(trace[seq_len(n)], tol = 1e-5), label = sprintf("n = %d", n))
  }
  # "No improvement over the best seen" would have stopped five steps into
  # the climb, with the atoms still moving. A plateau stops it: five steps
  # each below the tolerance.
  plateau <- utils::tail(climb, 1L) + 1e-4 * (1:5)
  expect_false(.hdp_converged(c(trace, plateau[1:4]), tol = 1e-5))
  expect_true(.hdp_converged(c(trace, plateau), tol = 1e-5))
  # One large step resets the count.
  expect_false(.hdp_converged(c(trace, plateau[1:3], -700, -700 + 1e-4), tol = 1e-5))
  # Too short a trace, or a non-finite value in the window, never converges.
  expect_false(.hdp_converged(c(-1, -1, -1, -1, -1), tol = 1e-5))
  expect_false(.hdp_converged(c(-1, NaN, -1, -1, -1, -1), tol = 1e-5))
  expect_true(.hdp_converged(c(-1, -1, -1, -1, -1, -1), tol = 1e-5))
})

test_that("R-only: max_iter caps an unconverged fit", {
  set.seed(20260729)
  blocky <- blocky_events()
  set.seed(7)
  h <- hdp_clustering(blocky, k_max = 8L, n_init = 1L, max_iter = 3L)[["hdp"]]
  expect_false(h[["converged"]])
  expect_identical(h[["n_iter"]], 3L)
})

test_that("R-only: a joint fit shares atoms across groups", {
  ev <- shared_and_private_events()
  set.seed(2026)
  res <- hdp_clustering(ev, k_max = 6L, n_init = 4L)
  ra <- res[["record_assignments"]]
  shared <- ra[["id"]] <= 80L
  # Both groups' S1/S2 records land in one cluster: one atom, not one per group.
  expect_length(unique(ra[["assignment"]][shared]), 1L)
  expect_identical(sort(unique(ra[["group"]][shared])), c("A", "B"))
})

test_that("R-only: record_item_matrix keeps multi-item records per group", {
  ev <- data.frame(id = c(1, 1, 1, 2, 3, 3, 4, 4),
                   item = c("x", "y", "x", "x", "y", "z", "x", "y"),
                   Group = c("B", "B", "B", "B", "A", "A", "A", "A"))
  m <- record_item_matrix(ev)
  expect_identical(m[["group_labels"]], c("A", "B"))
  expect_identical(m[["ids"]], c(3, 4, 1))
  expect_identical(m[["group"]], c(1L, 1L, 2L))
  expect_identical(colnames(m[["x"]]), c("x", "y", "z"))
  expect_identical(unname(m[["x"]][3L, ]), c(TRUE, TRUE, FALSE))  # a repeated x counts once
  expect_error(record_item_matrix(transform(ev, item = replace(item, 2L, NA))), "NA")
})

test_that("R-only: an item named like a profile column raises", {
  ev <- events_from_records(list(list(items = c("cluster", "y"), n = 5L)))
  ev[["Group"]] <- "A"
  expect_error(hdp_clustering(ev, k_max = 2L, n_init = 1L), "collide")
})

test_that("R-only: Beta KL and stick means match independent formulas", {
  # -KL(Beta(a, b) || Beta(a0, b0)) by numerical integration.
  kl_num <- function(a, b, a0, b0) {
    f <- function(t) {
      log_ratio <- stats::dbeta(t, a, b, log = TRUE) - stats::dbeta(t, a0, b0, log = TRUE)
      stats::dbeta(t, a, b) * log_ratio
    }
    -stats::integrate(f, 0, 1)[["value"]]
  }
  expect_equal(.neg_kl_beta(2, 5, 1, 1), kl_num(2, 5, 1, 1), tolerance = 1e-6)
  expect_equal(.neg_kl_beta(3.5, 1.2, 2, 4), kl_num(3.5, 1.2, 2, 4), tolerance = 1e-6)
  # Stick means: v = (0.5, 0.5, last) -> (0.5, 0.25, 0.25).
  expect_equal(.stick_means(c(1, 1, 1), c(1, 1, 1e-6)), c(0.5, 0.25, 0.25))
  expect_identical(hdp_effective_k(c(0.5, 0.3, 0.15, 0.05)), 3L)
  expect_identical(hdp_effective_k(c(0.4, 0.4), threshold = 0.95), 2L)
})
