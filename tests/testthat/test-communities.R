# Community detection via igraph (Leiden default, Louvain option). These
# replace the Julia-parity tests of the old phase-1-only Louvain, a
# documented divergence (CLAUDE.md, "Divergences from the Julia reference").
# Expected partitions come from fixtures that are unambiguous by
# construction, not from the Julia reference.

test_that("R-only: .with_seed seeds locally and restores the caller's stream", {
  set.seed(42)
  expected_next <- stats::runif(1)
  set.seed(42)
  a <- .with_seed(1L, stats::runif(3))
  expect_identical(stats::runif(1), expected_next)
  expect_identical(.with_seed(1L, stats::runif(3)), a)
})

test_that("R-only: .with_seed(NULL) draws from the caller's stream", {
  set.seed(3)
  a <- .with_seed(NULL, stats::runif(2))
  set.seed(3)
  expect_identical(stats::runif(2), a)
})

test_that("R-only: .with_seed leaves no .Random.seed behind in a fresh session", {
  key <- ".Random.seed"
  had <- exists(key, envir = globalenv(), inherits = FALSE)
  if (had) saved <- get(key, envir = globalenv(), inherits = FALSE)
  on.exit(if (had) assign(key, saved, envir = globalenv()))
  if (had) rm(list = key, envir = globalenv())
  .with_seed(1L, stats::runif(1))
  expect_false(exists(key, envir = globalenv(), inherits = FALSE))
})

test_that("R-only: .check_seed accepts NULL and whole numbers, raises otherwise", {
  expect_silent(.check_seed(NULL, "f"))
  expect_silent(.check_seed(7, "f"))
  expect_silent(.check_seed(-5L, "f"))
  for (bad in list(1.5, NA, NA_integer_, c(1, 2), "1", 2^31, Inf)) {
    expect_error(.check_seed(bad, "f"), "f: seed must be NULL or one whole number")
  }
})
