# R-only: the HDP ELBO against elbo_from_definition() (helper-hdp-elbo.R),
# an evaluator written from the model's definition rather than from
# R/hdp_clustering.R. The known-answer traces (test-hdp-known-answers.R) come
# from the Julia package alone; here each of their iterations is recomputed
# independently, so an ELBO bug shared by Julia and R cannot pass.

test_that("R-only: the textbook Beta KL agrees with numerical integration", {
  cases <- rbind(c(2, 3, 1, 1), c(5.5, 1.2, 1.3, 0.7), c(12, 40, 3, 9), c(1, 1, 2, 2),
                 c(30.2, 2.4, 1, 0.7))
  for (i in seq_len(nrow(cases))) {
    p <- cases[i, ]
    expect_equal(kl_beta_textbook(p[[1L]], p[[2L]], p[[3L]], p[[4L]]),
                 kl_beta_numeric(p[[1L]], p[[2L]], p[[3L]], p[[4L]]), tolerance = 1e-8)
  }
  expect_identical(kl_beta_textbook(2, 3, 2, 3), 0)
})

test_that("R-only: mean_weights gives stick-breaking weights that sum to 1", {
  w <- mean_weights(c(1, 1, 1), c(1, 1, 1e-6))
  expect_equal(w, c(0.5, 0.25, 0.25))
  expect_equal(sum(mean_weights(c(3, 1, 2, 1), c(1, 4, 2, 1e-6))), 1)
})

# Every .hdp_elbo() call of a fit: its arguments (the variational state) and
# its value. The call order is the iteration order.
capture_elbo_states <- function(fit_call) {
  states <- list()
  original <- .hdp_elbo
  testthat::local_mocked_bindings(.hdp_elbo = function(...) {
    value <- original(...)
    states[[length(states) + 1L]] <<- list(args = list(...), value = value)
    value
  })
  fit <- fit_call()
  list(fit = fit, states = states)
}

# .hdp_elbo(llp, r, theta_alpha, theta_beta, beta_alpha, beta_beta, beta_mean,
#           pi_alpha, pi_beta, alpha_prior, beta_prior, alpha, gamma) is called
# positionally. llp and beta_mean are the package's own derived values, so the
# evaluator takes neither: it recomputes both from x and the stick shapes.
independent_elbo <- function(state, x, group, alpha, gamma) {
  a <- state[["args"]]
  elbo_from_definition(x, group, r = a[[2L]], theta_a = a[[3L]], theta_b = a[[4L]],
                       v_a = a[[5L]], v_b = a[[6L]], w_a = a[[8L]], w_b = a[[9L]],
                       a0 = a[[10L]], b0 = a[[11L]], alpha = alpha, gamma = gamma)
}

for (name in names(hdp_known_answers)) {
  label <- sprintf("R-only: every ELBO of known answer %s, recomputed from the definition", name)
  test_that(label, {
    case <- hdp_known_answers[[name]]
    x <- ka_x(q = case[["q"]])
    run <- capture_elbo_states(function() ka_fit(case))
    expect_identical(length(run[["states"]]), run[["fit"]][["n_iter"]])
    mine <- vapply(run[["states"]], independent_elbo, numeric(1), x = x, group = ka_group,
                   alpha = case[["alpha"]], gamma = case[["gamma"]])
    pkg <- vapply(run[["states"]], function(s) s[["value"]], numeric(1))
    expect_equal(mine, pkg, tolerance = 1e-10)
    # The Julia trace, now with a second source.
    expect_equal(mine, case[["elbo_trace"]], tolerance = 1e-9)
  })
}
