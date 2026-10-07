# R12: the HDP is hand-written (no R package fits it), so cross-check it on the
# one case an R package does fit: one group, three well-separated Bernoulli
# components. A small-gamma HDP should recover the truth as well as a flat
# K = 3 maximum-likelihood fit from flexmix. Tolerances come from five seeds
# (HDP ARI 0.975-0.995, equal to flexmix to the digit; worst matched item
# probability error 0.064), with margin.

adjusted_rand <- function(a, b) {
  pairs <- function(n) n * (n - 1) / 2
  tab <- table(a, b)
  sum_cells <- sum(pairs(tab))
  sum_rows <- sum(pairs(rowSums(tab)))
  sum_cols <- sum(pairs(colSums(tab)))
  expected <- sum_rows * sum_cols / pairs(sum(tab))
  (sum_cells - expected) / ((sum_rows + sum_cols) / 2 - expected)
}

# For each true component, the largest item-probability error of the closest
# fitted component.
worst_matched_error <- function(fitted, truth) {
  vapply(seq_len(nrow(truth)), function(i) {
    min(apply(fitted, 1L, function(r) max(abs(r - truth[i, ]))))
  }, numeric(1))
}

test_that("HDP recovers three Bernoulli components as well as flexmix (R12)", {
  skip_if_not_installed("flexmix")
  set.seed(20261007)
  n <- 600L
  d <- 12L
  truth <- matrix(0.05, 3L, d)
  truth[1L, 1:4] <- truth[2L, 5:8] <- truth[3L, 9:12] <- 0.85
  z <- sample(3L, n, replace = TRUE, prob = c(0.4, 0.35, 0.25))
  x <- matrix(runif(n * d) < truth[z, ], n, d)

  hdp <- fit_hdp_bernoulli_mixture(x, rep(1L, n), 1L, k_max = 10L, gamma = 0.5, n_init = 5L)
  flex <- flexmix::flexmix(x + 0 ~ 1, k = 3L, model = flexmix::FLXMCmvbinary(),
                           control = list(minprior = 0))

  expect_equal(hdp$effective_k, 3L)
  ari_hdp <- adjusted_rand(z, hdp$assignments)
  ari_flex <- adjusted_rand(z, flexmix::clusters(flex))
  expect_gt(ari_hdp, 0.9)
  expect_gte(ari_hdp, ari_flex - 0.02)

  active <- sort(unique(hdp$assignments))
  expect_length(active, 3L)
  expect_lt(max(worst_matched_error(hdp$theta[active, , drop = FALSE], truth)), 0.10)
  expect_lt(max(worst_matched_error(t(flexmix::parameters(flex)), truth)), 0.10)
})
