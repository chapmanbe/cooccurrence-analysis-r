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
