# Sourced by testthat before any test file. The same silent-failure guards
# the consuming analysis project applies to every script and test: partial
# `$` / argument / attribute matching warns, and collation is pinned to "C"
# so sort order and factor levels do not depend on the machine. This package
# is standalone, so it carries its own copy rather than sourcing the
# consumer's.
options(
  warnPartialMatchDollar = TRUE,
  warnPartialMatchArgs = TRUE,
  warnPartialMatchAttr = TRUE,
  stringsAsFactors = FALSE
)
invisible(Sys.setlocale("LC_COLLATE", "C"))
