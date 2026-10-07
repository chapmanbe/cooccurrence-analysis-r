#' Round to `digits` decimal places the way the Julia reference does
#'
#' `round(x * 10^digits) / 10^digits`: binary floating-point scaling, then
#' R's `round()` at zero digits, which breaks ties to even exactly as Julia's
#' `round(x; digits = d)` does. R's own `round(x, digits = d)` for `d > 0`
#' uses a different, "decimally correct" algorithm and disagrees with the
#' reference on ties (e.g. `round(0.15, 1)`), so the published tables could
#' not be reproduced with it. Non-finite values pass through unchanged.
#'
#' Vectorized over `x`; `digits` must be a single non-negative whole number.
#' @export
round_digits <- function(x, digits) {
  stopifnot(is.numeric(x), length(digits) == 1L, digits >= 0, digits == round(digits))
  scale <- 10^digits
  out <- round(x * scale) / scale
  nonfinite <- !is.finite(x)
  out[nonfinite] <- x[nonfinite]
  out
}

#' Stop unless `value` is one of `allowed`, naming the argument
#' @keywords internal
.check_choice <- function(value, allowed, arg) {
  if (!(is.character(value) && length(value) == 1L && value %in% allowed)) {
    stop(sprintf("unsupported %s '%s'; expected one of: %s", arg,
                 paste(format(value), collapse = ","), paste(allowed, collapse = ", ")),
         call. = FALSE)
  }
  invisible(value)
}

#' Sparse record-by-item incidence matrix
#'
#' `rec` and `code` give each event's record (`1..n_rec`) and item
#' (`1..n_item`). Returns an `n_rec x n_item` `dgCMatrix` holding 1 where the
#' record holds the item; an event repeated within a record counts once.
#' Numeric rather than logical so `crossprod()` returns counts, not a pattern.
#' @keywords internal
.incidence <- function(rec, code, n_rec, n_item) {
  keep <- !duplicated(as.numeric(rec) * (n_item + 1) + code)
  Matrix::sparseMatrix(i = rec[keep], j = code[keep], x = 1, dims = c(n_rec, n_item))
}

#' Sort character values by byte (C-locale) order, whatever the session locale
#' @keywords internal
.sort_c <- function(x) sort(x, method = "radix")

#' Raise unless `seed` is NULL or one whole number `set.seed()` accepts
#' @keywords internal
.check_seed <- function(seed, caller) {
  ok <- is.null(seed) ||
    (is.numeric(seed) && length(seed) == 1L && is.finite(seed) && seed == round(seed) &&
       abs(seed) <= .Machine$integer.max)
  if (!ok) {
    stop(sprintf(paste0("%s: seed must be NULL or one whole number within the integer range ",
                        "(NULL draws from the current RNG stream)"), caller), call. = FALSE)
  }
  invisible(seed)
}

#' Evaluate `expr` with the RNG seeded to `seed`, then restore the caller's
#' stream (or remove `.Random.seed` if the caller had none), so a seeded call
#' never shifts the caller's later random draws. The generator is pinned to
#' R's defaults, so the same seed gives the same draws whatever `RNGkind()` the
#' caller has set. `seed = NULL` evaluates
#' `expr` on the caller's stream, which it advances.
#' @keywords internal
.with_seed <- function(seed, expr) {
  if (is.null(seed)) return(expr)
  env <- globalenv()
  key <- ".Random.seed"
  had <- exists(key, envir = env, inherits = FALSE)
  if (had) saved <- get(key, envir = env, inherits = FALSE)
  on.exit({
    if (had) {
      assign(key, saved, envir = env)
    } else if (exists(key, envir = env, inherits = FALSE)) {
      rm(list = key, envir = env)
    }
  })
  # Pin the generator too: a seed alone gives different draws under another
  # RNGkind (L'Ecuyer under parallel, say). Restoring .Random.seed below also
  # restores the caller's kind, which that vector encodes.
  set.seed(seed, kind = "Mersenne-Twister", normal.kind = "Inversion",
           sample.kind = "Rejection")
  expr
}
