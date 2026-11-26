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

test_that("R-only: code and tests use American spelling", {
  root <- testthat::test_path("..", "..")
  skip_if_not(file.exists(file.path(root, "DESCRIPTION")), "package source not available")
  paths <- c(list.files(file.path(root, "R"), "\\.R$", full.names = TRUE),
             list.files(file.path(root, "tests", "testthat"), "\\.R$", full.names = TRUE))
  # ggplot2 names its built-data column "colour"; that spelling is its, not ours.
  british <- paste0("colour|labell(ed|ing)|modell|behaviour|favour|grey|centre|programme|",
                    "isation|towards|amongst|whilst")
  this_test <- "test-package-metadata.R"
  hits <- unlist(lapply(paths[basename(paths) != this_test], function(f) {
    lines <- gsub('[[]{2}"colour"[]]{2}', "", readLines(f))
    i <- grep(british, lines, ignore.case = TRUE)
    if (length(i) > 0L) paste0(basename(f), ":", i) else character(0)
  }))
  expect_identical(hits, character(0))
})
