# R-only: package metadata the code depends on.

test_that("R-only: every pkg:: namespace used in R/ is declared in DESCRIPTION", {
  root <- testthat::test_path("..", "..")
  skip_if_not(file.exists(file.path(root, "DESCRIPTION")), "package source not available")
  desc <- read.dcf(file.path(root, "DESCRIPTION"), fields = c("Imports", "Suggests"))
  declared <- trimws(sub("\\(.*", "", unlist(strsplit(desc[!is.na(desc)], ","))))
  code <- unlist(lapply(list.files(file.path(root, "R"), "\\.R$", full.names = TRUE), readLines))
  code <- sub("#.*", "", code)
  used <- unique(unlist(regmatches(code, gregexpr("[A-Za-z][A-Za-z0-9.]*(?=:::?)", code,
                                                  perl = TRUE))))
  expect_identical(setdiff(used, c(declared, "base")), character(0))
})
