# HDP Bernoulli mixture: batch CAVI with truncated stick-breaking.
#
# R port of the Julia reference's hdp_clustering.jl. Model: Teh, Jordan, Beal
# & Blei 2006 (hierarchical Dirichlet process). Algorithm: in the spirit of
# Wang, Paisley & Blei 2011 section 3.1 batch CAVI, adapted to Bernoulli
# observations, using the single-level shared-atom stick-breaking
# representation: record assignment z_i selects directly from corpus-level
# atom k. Bernoulli atoms with empirical-Bayes Beta priors follow Ye, Zhang &
# Nie 2018.
#
# Every update below is the reference's, line for line; only the random
# stream differs (R's Mersenne-Twister, not Julia's), so a seed does not
# carry across. The global-stick update is the reference's documented
# approximation (see .hdp_update_global_sticks) and is deliberately kept.

#' Empirical-Bayes Beta priors per item (Ye et al. 2018)
#'
#' For a logical record-by-item matrix `x` with item prevalence `p_d`,
#' returns `list(alpha = concentration * p + floor, beta = concentration *
#' (1 - p) + floor)`. `floor` keeps the prior proper for an item that never
#' occurs. Raises on a matrix with no rows.
#' @export
empirical_bayes_priors <- function(x, concentration = 10, floor = 1) {
  if (!is.matrix(x) || !is.logical(x)) stop("x must be a logical matrix", call. = FALSE)
  if (anyNA(x)) stop("x contains NA", call. = FALSE)
  n <- nrow(x)
  if (n == 0L) stop("cannot compute empirical Bayes priors on empty data", call. = FALSE)
  p <- colSums(x) / n
  list(alpha = concentration * p + floor, beta = concentration * (1 - p) + floor)
}

# -KL(Beta(a, b) || Beta(a0, b0)), elementwise.
.neg_kl_beta <- function(a, b, a0, b0) {
  e_log_x <- digamma(a) - digamma(a + b)
  e_log_1mx <- digamma(b) - digamma(a + b)
  lgamma(a) + lgamma(b) - lgamma(a + b) - lgamma(a0) - lgamma(b0) + lgamma(a0 + b0) -
    (a - a0) * e_log_x - (b - b0) * e_log_1mx
}

# tail[k] = sum(x[(k + 1):K]), tail[K] = 0, accumulated from the end as the
# reference does.
.tail_sums <- function(x) c(rev(cumsum(rev(x)))[-1L], 0)

# E_q[log pi_{j,k}] for each row of stick parameters (groups x K).
.e_log_pi <- function(pi_alpha, pi_beta) {
  t1 <- digamma(pi_alpha) - digamma(pi_alpha + pi_beta)
  t2 <- digamma(pi_beta) - digamma(pi_alpha + pi_beta)
  out <- t1
  cum <- numeric(nrow(t1))
  for (kk in seq_len(ncol(t1))) {
    out[, kk] <- t1[, kk] + cum
    cum <- cum + t2[, kk]
  }
  out
}

# Stick-breaking means: v_k = a_k / (a_k + b_k); the last stick takes the rest.
.stick_means <- function(a, b) {
  k <- length(a)
  out <- numeric(k)
  log_remain <- 0
  for (kk in seq_len(k - 1L)) {
    vk <- a[[kk]] / (a[[kk]] + b[[kk]])
    out[[kk]] <- exp(log(max(vk, 1e-300)) + log_remain)
    log_remain <- log_remain + log(max(1 - vk, 1e-300))
  }
  out[[k]] <- exp(log_remain)
  out
}

.stick_means_rows <- function(a, b) {
  out <- a
  for (j in seq_len(nrow(a))) out[j, ] <- .stick_means(a[j, ], b[j, ])
  out
}

#' Effective number of components
#'
#' The fewest components whose largest weights sum to at least `threshold`;
#' all of them if none do.
#' @export
hdp_effective_k <- function(beta_mean, threshold = 0.95) {
  hit <- which(cumsum(sort(beta_mean, decreasing = TRUE)) >= threshold)
  if (length(hit) == 0L) length(beta_mean) else hit[[1L]]
}

# The record-by-item matrix as the positions of its TRUE cells. Records
# hold a handful of items out of dozens, so the two products CAVI needs
# every iteration are sums over those cells rather than dense products.
.sparse_binary <- function(x) {
  list(rows = lapply(seq_len(ncol(x)), function(d) which(x[, d])), n = nrow(x), d = ncol(x))
}

# X %*% t(w) for a K x D matrix w: N x K, row i the sum of w's columns at
# record i's items (zero for a record with none).
.sx_times <- function(sx, w) {
  out <- matrix(0, sx[["n"]], nrow(w))
  for (k in seq_len(nrow(w))) {
    v <- numeric(sx[["n"]])
    for (d in seq_len(sx[["d"]])) {
      rows <- sx[["rows"]][[d]]
      v[rows] <- v[rows] + w[k, d]
    }
    out[, k] <- v
  }
  out
}

# crossprod(r, X) for an N x K matrix r: K x D, column d the sum of r's rows
# at the records holding item d.
.sx_crossprod <- function(r, sx) {
  out <- matrix(0, ncol(r), sx[["d"]])
  for (d in seq_len(sx[["d"]])) {
    rows <- sx[["rows"]][[d]]
    if (length(rows) > 0L) out[, d] <- colSums(r[rows, , drop = FALSE])
  }
  out
}

# Row-wise maximum of a matrix, without apply().
.row_max <- function(m) {
  out <- m[, 1L]
  for (kk in seq_len(ncol(m))[-1L]) out <- pmax(out, m[, kk])
  out
}

# N x K expected log-likelihood E_q[log p(x_i | theta_k)], plus the
# per-group prior term E_q[log pi_{g(i),k}].
.hdp_log_lik <- function(xf, group, theta_alpha, theta_beta, e_log_pi) {
  e_log_th <- digamma(theta_alpha) - digamma(theta_alpha + theta_beta)
  e_log_1m <- digamma(theta_beta) - digamma(theta_alpha + theta_beta)
  ll <- .sx_times(xf, e_log_th - e_log_1m)
  ll <- ll + rep(rowSums(e_log_1m), each = xf[["n"]])
  ll + e_log_pi[group, , drop = FALSE]
}

# E-step: responsibilities from the .hdp_log_lik() matrix, normalized by
# row-wise log-sum-exp.
.hdp_e_step <- function(log_r) {
  r <- exp(log_r - .row_max(log_r))
  r / rowSums(r)
}

# Atoms: q(theta_{k,d}) = Beta(alpha_d + sum_i r_ik x_id, beta_d + sum_i r_ik (1 - x_id)).
.hdp_update_atoms <- function(xf, r, alpha_prior, beta_prior) {
  k <- ncol(r)
  wx <- .sx_crossprod(r, xf)
  w1mx <- colSums(r) - wx
  list(alpha = wx + rep(alpha_prior, each = k), beta = w1mx + rep(beta_prior, each = k))
}

# Global sticks. APPROXIMATION (the reference's, documented and kept): this
# pools the raw soft counts across groups (N_k = sum_i r_ik) rather than the
# group-level table counts a full Chinese-restaurant-franchise update would
# use, so the ELBO is a valid lower bound but not guaranteed monotone.
.hdp_update_global_sticks <- function(r, gamma) {
  k <- ncol(r)
  n_k <- unname(colSums(r))
  tail <- .tail_sums(n_k)
  list(alpha = c(1 + n_k[-k], 1), beta = c(pmax(gamma + tail[-k], 1e-6), 1e-6))
}

# Per-group sticks, coupled to the global weights.
.hdp_update_group_sticks <- function(r, group, n_groups, beta_mean, alpha) {
  n_jk <- rowsum(r, group, reorder = TRUE)
  if (nrow(n_jk) != n_groups) stop("a group has no records", call. = FALSE)
  dimnames(n_jk) <- NULL
  beta_tail <- .tail_sums(beta_mean)
  n_tail <- n_jk
  for (j in seq_len(n_groups)) n_tail[j, ] <- .tail_sums(n_jk[j, ])
  list(alpha = pmax(sweep(n_jk, 2L, alpha * beta_mean, "+"), 1e-6),
       beta = pmax(sweep(n_tail, 2L, alpha * beta_tail, "+"), 1e-6))
}

# `llp` is .hdp_log_lik() at the current theta and per-group sticks -- the
# same matrix the next E-step starts from, so the caller computes it once.
.hdp_elbo <- function(llp, r, theta_alpha, theta_beta, beta_alpha, beta_beta,
                      beta_mean, pi_alpha, pi_beta, alpha_prior, beta_prior,
                      alpha, gamma) {
  k <- ncol(r)
  keep <- r > 1e-300
  l_data <- sum(r[keep] * (llp[keep] - log(r[keep])))

  l_theta <- sum(.neg_kl_beta(theta_alpha, theta_beta, rep(alpha_prior, each = k),
                              rep(beta_prior, each = k)))

  beta_tail <- .tail_sums(beta_mean)
  head <- seq_len(k - 1L)
  a0 <- pmax(alpha * beta_mean[head], 1e-6)
  b0 <- pmax(alpha * beta_tail[head], 1e-6)
  l_pi <- 0
  for (j in seq_len(nrow(pi_alpha))) {
    l_pi <- l_pi + sum(.neg_kl_beta(pi_alpha[j, head], pi_beta[j, head], a0, b0))
  }
  l_beta <- sum(.neg_kl_beta(beta_alpha[head], beta_beta[head], 1, gamma))
  l_data + l_theta + l_pi + l_beta
}

# Consecutive small ELBO steps required to stop. The test is on the
# per-iteration delta, NOT on improvement over the best ELBO seen: the
# objective is non-monotone (see the global-stick note), so "no improvement
# over the best" fills the window within a few iterations while the atoms
# still sit on the marginal prevalences. That was a real bug in the
# reference (a 153,700-record fit stopped at iteration 5 with 3
# undifferentiated clusters instead of 6), invisible on small fixtures.
.PLATEAU_WINDOW <- 5L # nolint: object_name_linter. A constant, as in the reference.

# Has CAVI converged, given the ELBO after each iteration so far? TRUE when
# each of the last `window` iterations changed the ELBO by less than
# `tol * (1 + |previous ELBO|)` -- the reference's counter of consecutive
# small steps, which a step of any size resets. The first iteration has no
# previous value and never counts. A pure function so the rule itself is
# tested on hand-made traces, including the non-monotone ones real data
# produce and small fixtures do not.
.hdp_converged <- function(trace, tol, window = .PLATEAU_WINDOW) {
  n <- length(trace)
  if (n < window + 1L) return(FALSE)
  prev <- trace[(n - window):(n - 1L)]
  step <- abs(trace[(n - window + 1L):n] - prev)
  isTRUE(all(is.finite(prev)) && all(step < tol * (1 + abs(prev))))
}

# One CAVI run from a random start: `init` (N x K, positive, rows normalized
# here as the reference normalizes its uniform draws), or, when NULL, a
# uniform draw from R's global stream.
.run_hdp_cavi <- function(xf, group, n_groups, k_max, alpha, gamma, alpha_prior, beta_prior,
                          max_iter, tol, effective_k_threshold, init = NULL) {
  n <- xf[["n"]]
  r <- if (is.null(init)) matrix(stats::runif(n * k_max), n, k_max) else init
  r <- r / rowSums(r)
  atoms <- .hdp_update_atoms(xf, r, alpha_prior, beta_prior)
  theta_alpha <- atoms[["alpha"]]
  theta_beta <- atoms[["beta"]]
  beta_a <- rep(1, k_max)
  beta_b <- c(rep(1, k_max - 1L), 1e-6)   # the last stick absorbs the rest
  pi_a <- matrix(1, n_groups, k_max)
  pi_b <- matrix(1, n_groups, k_max)
  pi_b[, k_max] <- 1e-6
  e_log_pi <- .e_log_pi(pi_a, pi_b)

  converged <- FALSE
  n_iter <- max_iter
  trace <- numeric(0)
  llp <- .hdp_log_lik(xf, group, theta_alpha, theta_beta, e_log_pi)
  for (iter in seq_len(max_iter)) {
    r <- .hdp_e_step(llp)
    atoms <- .hdp_update_atoms(xf, r, alpha_prior, beta_prior)
    theta_alpha <- atoms[["alpha"]]
    theta_beta <- atoms[["beta"]]
    sticks <- .hdp_update_global_sticks(r, gamma)
    beta_a <- sticks[["alpha"]]
    beta_b <- sticks[["beta"]]
    beta_mean <- .stick_means(beta_a, beta_b)
    gs <- .hdp_update_group_sticks(r, group, n_groups, beta_mean, alpha)
    pi_a <- gs[["alpha"]]
    pi_b <- gs[["beta"]]
    e_log_pi <- .e_log_pi(pi_a, pi_b)

    llp <- .hdp_log_lik(xf, group, theta_alpha, theta_beta, e_log_pi)
    elbo <- .hdp_elbo(llp, r, theta_alpha, theta_beta, beta_a, beta_b, beta_mean,
                      pi_a, pi_b, alpha_prior, beta_prior, alpha, gamma)
    trace <- c(trace, elbo)
    if (.hdp_converged(trace, tol)) {
      converged <- TRUE
      n_iter <- iter
      break
    }
  }

  beta_mean <- .stick_means(beta_a, beta_b)
  pi_mean <- .stick_means_rows(pi_a, pi_b)
  # Components by global weight, descending; ties keep their order.
  ord <- order(-beta_mean, method = "radix")
  theta_alpha <- theta_alpha[ord, , drop = FALSE]
  theta_beta <- theta_beta[ord, , drop = FALSE]
  r <- r[, ord, drop = FALSE]
  beta_mean <- beta_mean[ord]
  structure(list(
    k_max = k_max, n_groups = n_groups,
    group_labels = paste0("Group_", seq_len(n_groups)),
    theta = theta_alpha / (theta_alpha + theta_beta),
    theta_alpha = theta_alpha, theta_beta = theta_beta,
    beta_alpha = beta_a[ord], beta_beta = beta_b[ord], beta_mean = beta_mean,
    pi_mean = pi_mean[, ord, drop = FALSE], pi_alpha = pi_a[, ord, drop = FALSE],
    pi_beta = pi_b[, ord, drop = FALSE],
    responsibilities = r,
    assignments = max.col(r, ties.method = "first"),
    group = group,
    elbo = trace[[length(trace)]], elbo_trace = trace, n_iter = n_iter, converged = converged,
    effective_k = hdp_effective_k(beta_mean, threshold = effective_k_threshold)
  ), class = "hdp_bernoulli_result")
}

#' Fit an HDP Bernoulli mixture by CAVI
#'
#' `x` is an N x D logical record-by-item matrix and `group` an integer
#' vector of length N with values in `1:n_groups`. Truncated stick-breaking
#' at `k_max` components; closed-form Beta-Bernoulli updates.
#'
#' - `alpha`: per-group DP concentration (larger: groups closer to global).
#' - `gamma`: global DP concentration (smaller: fewer effective clusters).
#' - `prior`: `"flat"` (Beta(`alpha_prior`, `beta_prior`) on every atom) or
#'   `"empirical_bayes"` ([empirical_bayes_priors()] with `concentration`,
#'   `floor`).
#' - `max_iter`, `tol`: CAVI stops after 5 consecutive iterations whose ELBO
#'   change is below `tol * (1 + |previous ELBO|)`, or at `max_iter`.
#' - `n_init`: random restarts from R's global random stream (seed it first);
#'   the restart with the highest ELBO is returned.
#' - `init`: an explicit starting point instead of a random one -- an N x
#'   `k_max` matrix of positive values, normalized by row into the initial
#'   responsibilities (the reference's start is a uniform draw, normalized
#'   the same way). Requires `n_init = 1`. Given the same matrix, the fit is
#'   deterministic and comparable value for value with the reference.
#'
#' Returns an `hdp_bernoulli_result` list: `k_max`, `n_groups`,
#' `group_labels`, `theta` (K x D posterior means), `theta_alpha`,
#' `theta_beta`, `beta_alpha`, `beta_beta`, `beta_mean` (global weights),
#' `pi_mean`, `pi_alpha`, `pi_beta` (groups x K), `responsibilities` (N x
#' K), `assignments` (argmax, first on ties), `group`, `elbo`, `elbo_trace`,
#' `n_iter`, `converged`, `effective_k`. Components are sorted by
#' `beta_mean`, descending.
#'
#' Raises on an empty group, on bad shapes, and if every restart ends with a
#' non-finite ELBO; warns when N < `k_max`.
#' @export
fit_hdp_bernoulli_mixture <- function(x, group, n_groups, k_max = 20L, alpha = 1, gamma = 1,
                                      prior = "flat", alpha_prior = 1, beta_prior = 1,
                                      concentration = 10, floor = 1, max_iter = 300L,
                                      tol = 1e-5, n_init = 3L, effective_k_threshold = 0.95,
                                      init = NULL) {
  if (!is.matrix(x) || !is.logical(x) || anyNA(x)) {
    stop("x must be a logical matrix without NA", call. = FALSE)
  }
  .check_choice(prior, c("flat", "empirical_bayes"), "prior")
  n <- nrow(x)
  d <- ncol(x)
  if (length(group) != n) stop("group length must equal the number of rows in x", call. = FALSE)
  if (n_groups < 1L) stop("n_groups must be at least 1", call. = FALSE)
  if (anyNA(group) || any(group != round(group)) || any(group < 1L | group > n_groups)) {
    stop("group values must be whole numbers in 1:n_groups", call. = FALSE)
  }
  if (k_max < 2L) stop("k_max must be at least 2", call. = FALSE)
  if (n_init < 1L) stop("n_init must be at least 1", call. = FALSE)
  if (max_iter < 1L) stop("max_iter must be at least 1", call. = FALSE)
  if (!is.null(init)) {
    if (n_init != 1L) stop("an explicit init needs n_init = 1", call. = FALSE)
    if (!is.matrix(init) || !is.numeric(init) || !identical(dim(init), c(n, as.integer(k_max))) ||
          anyNA(init) || any(init <= 0)) {
      stop(sprintf("init must be a positive numeric %d x %d matrix", n, as.integer(k_max)),
           call. = FALSE)
    }
    dimnames(init) <- NULL
  }
  group <- as.integer(group)
  counts <- tabulate(group, nbins = n_groups)
  if (any(counts == 0L)) {
    stop(sprintf("empty group detected: group(s) %s have zero records",
                 paste(which(counts == 0L), collapse = ", ")), call. = FALSE)
  }
  if (n < k_max) {
    warning(sprintf("N=%d < k_max=%d; some components will be empty", n, k_max), call. = FALSE)
  }

  if (prior == "empirical_bayes") {
    pr <- empirical_bayes_priors(x, concentration = concentration, floor = floor)
    a_prior <- unname(pr[["alpha"]])
    b_prior <- unname(pr[["beta"]])
  } else {
    a_prior <- rep(alpha_prior, d)
    b_prior <- rep(beta_prior, d)
  }

  xf <- .sparse_binary(x)
  best <- NULL
  best_elbo <- -Inf
  for (i in seq_len(n_init)) {
    res <- .run_hdp_cavi(xf, group, n_groups, as.integer(k_max), alpha, gamma, a_prior, b_prior,
                         as.integer(max_iter), tol, effective_k_threshold, init = init)
    # A non-finite ELBO never wins: NaN > best is NA, not TRUE.
    if (isTRUE(res[["elbo"]] > best_elbo)) {
      best_elbo <- res[["elbo"]]
      best <- res
    }
  }
  if (is.null(best)) {
    stop(sprintf(paste0("HDP CAVI failed: all %d restart(s) produced a non-finite ELBO. ",
                        "Try more records, a smaller k_max, or another seed."), n_init),
         call. = FALSE)
  }
  best
}
