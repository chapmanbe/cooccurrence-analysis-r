# An evidence lower bound for the HDP Bernoulli mixture, written out from the
# model's definition, independent of R/hdp_clustering.R, so the package's
# ELBO (and the Julia traces it reproduces) can be checked against it.
#
# Model, truncated at K components:
#   global sticks   v_k ~ Beta(1, gamma), k < K; v_K = 1
#   global weights  b_k = v_k prod_{l<k} (1 - v_l)
#   group sticks    w_jk ~ Beta(alpha b_k, alpha sum_{l>k} b_l), k < K
#   group weights   pi_jk = w_jk prod_{l<k} (1 - w_jl)
#   atoms           theta_kd ~ Beta(a0_d, b0_d)
#   assignment      z_i ~ Categorical(pi_{g(i)})
#   data            x_id ~ Bernoulli(theta_{z_i, d})
# Variational family: q(z_i) = Categorical(r_i), q(theta_kd) = Beta, q(v_k) =
# Beta, q(w_jk) = Beta. Two conventions are the package's, stated here so the
# check is of the arithmetic, not of the choice:
#   - the group-stick prior uses the global weights at their variational mean
#     (a plug-in for E_q over b, which has no closed form), each shape floored
#     at 1e-6;
#   - the last group stick keeps its q = Beta(a, 1e-6) in E_q[log pi], rather
#     than being fixed at 1 (the difference is of order 1e-6 per record).
#
# ELBO = sum_i E_q[log p(x_i | z_i, theta)] + E_q[log p(z_i | pi)] + H[q(z_i)]
#        - sum KL(q(theta) || p(theta)) - sum KL(q(w) || p(w)) - sum KL(q(v) || p(v)).

#' KL(Beta(a, b) || Beta(a0, b0)) in the textbook form
#' log B(a0, b0) - log B(a, b) + (a - a0) psi(a) + (b - b0) psi(b)
#'   + (a0 - a + b0 - b) psi(a + b).
kl_beta_textbook <- function(a, b, a0, b0) {
  lbeta(a0, b0) - lbeta(a, b) + (a - a0) * digamma(a) + (b - b0) * digamma(b) +
    (a0 - a + b0 - b) * digamma(a + b)
}

#' The same KL by numerical integration of q log(q / p), to check the form above.
kl_beta_numeric <- function(a, b, a0, b0) {
  f <- function(t) {
    q <- stats::dbeta(t, a, b)
    ifelse(q > 0, q * (stats::dbeta(t, a, b, log = TRUE) - stats::dbeta(t, a0, b0, log = TRUE)), 0)
  }
  stats::integrate(f, 0, 1, rel.tol = 1e-10)[["value"]]
}

#' E_q[log w] and E_q[log (1 - w)] for w ~ Beta(a, b).
e_log_beta <- function(a, b) {
  list(log = digamma(a) - digamma(a + b), log1m = digamma(b) - digamma(a + b))
}

#' Mean stick-breaking weights from stick shapes (independence of the sticks
#' makes the mean of the product the product of the means); the last weight
#' is what the first K - 1 sticks leave.
mean_weights <- function(a, b) {
  k <- length(a)
  v <- a[-k] / (a[-k] + b[-k])
  c(v * cumprod(c(1, 1 - v))[-k], prod(1 - v))
}

#' The ELBO at one variational state. `x` is the logical record-by-item
#' matrix, `group` each record's group index, `r` the N x K responsibilities;
#' theta_*, v_*, w_* are the q shapes (K x D, length K, G x K); a0, b0 the
#' per-item prior shapes.
elbo_from_definition <- function(x, group, r, theta_a, theta_b, v_a, v_b, w_a, w_b,
                                 a0, b0, alpha, gamma) {
  n <- nrow(x)
  k <- ncol(r)
  n_groups <- nrow(w_a)
  th <- e_log_beta(theta_a, theta_b)

  # E_q[log pi_jk] = E log w_jk + sum_{l<k} E log (1 - w_jl).
  e_log_pi <- matrix(0, n_groups, k)
  for (j in seq_len(n_groups)) {
    w <- e_log_beta(w_a[j, ], w_b[j, ])
    for (kk in seq_len(k)) {
      e_log_pi[j, kk] <- w[["log"]][[kk]] + sum(w[["log1m"]][seq_len(kk - 1L)])
    }
  }

  # Data, assignment and entropy terms, one record and component at a time.
  per_record <- numeric(n)
  for (i in seq_len(n)) {
    xi <- x[i, ]
    total <- 0
    for (kk in seq_len(k)) {
      rik <- r[i, kk]
      if (rik == 0) next
      log_lik <- sum(ifelse(xi, th[["log"]][kk, ], th[["log1m"]][kk, ]))
      total <- total + rik * (log_lik + e_log_pi[group[[i]], kk] - log(rik))
    }
    per_record[[i]] <- total
  }

  kl_theta <- 0
  for (kk in seq_len(k)) kl_theta <- kl_theta + sum(kl_beta_textbook(theta_a[kk, ], theta_b[kk, ],
                                                                     a0, b0))

  head <- seq_len(k - 1L)
  b_mean <- mean_weights(v_a, v_b)
  b_tail <- vapply(head, function(kk) sum(b_mean[(kk + 1L):k]), numeric(1))
  prior_a <- pmax(alpha * b_mean[head], 1e-6)
  prior_b <- pmax(alpha * b_tail, 1e-6)
  kl_w <- 0
  for (j in seq_len(n_groups)) {
    kl_w <- kl_w + sum(kl_beta_textbook(w_a[j, head], w_b[j, head], prior_a, prior_b))
  }

  kl_v <- sum(kl_beta_textbook(v_a[head], v_b[head], 1, gamma))

  sum(per_record) - kl_theta - kl_w - kl_v
}
