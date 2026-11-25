# The known-answer fixture, shared by test-hdp-known-answers.R and
# test-hdp-elbo.R.

# Restated from tools/hdp_known_answers.jl: 60 records x 9 items in three
# blocks; block items present unless (i * j) %% q == 0, others present when
# (i + 2j) %% 11 == 0. Groups: records 1-30, 31-60.
ka_x <- function(n = 60L, q = 5L) {
  x <- matrix(FALSE, n, 9L)
  for (i in seq_len(n)) {
    for (j in 1:9) {
      in_block <- (j - 1L) %/% 3L == i %% 3L
      x[i, j] <- if (in_block) (i * j) %% q != 0L else (i + 2L * j) %% 11L == 0L
    }
  }
  x
}
ka_init <- function(n = 60L, k = 6L, p = 17L) {
  outer(seq_len(n), seq_len(k), function(i, kk) ((i * 7 + kk * 13) %% p + 1) / (p + 1))
}
ka_group <- rep(1:2, each = 30L)

ka_fit <- function(case, max_iter = 300L) {
  fit_hdp_bernoulli_mixture(ka_x(q = case[["q"]]), ka_group, 2L, k_max = 6L,
                            prior = case[["prior"]], alpha = case[["alpha"]],
                            gamma = case[["gamma"]], n_init = 1L, max_iter = max_iter,
                            tol = 1e-5, init = ka_init(p = case[["p"]]))
}
