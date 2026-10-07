# Ported from the Julia suite's "mining" testsets. Expected values come from
# running the Julia reference (build_transactions(event_df; min_items = 2) then
# mine_frequent_itemsets / mine_association_rules, CooccurrenceAnalysis.jl at
# 2026-10-06) on the shared fixture: 20 multi-item records.

skip_if_not_installed("arules")
event_df <- make_test_event_df()

set_keys <- function(sets) vapply(sets, paste, character(1L), collapse = "+")

test_that("mine_frequent_itemsets", {
  fi <- mine_frequent_itemsets(event_df, min_support = 0.1)
  expect_identical(names(fi), c("Itemset", "Support", "N", "Length"))
  expect_true(all(fi[["Length"]] >= 2L))
  expect_identical(set_keys(fi[["Itemset"]]),
                   c("Item01+Item02", "GA1+Item05", "GB1+Item01", "Item03+Item04"))
  expect_identical(fi[["N"]], c(10L, 5L, 3L, 2L))
  expect_equal(fi[["Support"]], c(0.5, 0.25, 0.15, 0.1))
})

test_that("mine_frequent_itemsets with min_count: thresholds are inclusive", {
  # Julia keeps GA1+Item05 (N = 5) at min_count = 5 and at min_support = 0.25.
  for (fi in list(mine_frequent_itemsets(event_df, min_count = 5L),
                  mine_frequent_itemsets(event_df, min_support = 0.25))) {
    expect_identical(set_keys(fi[["Itemset"]]), c("Item01+Item02", "GA1+Item05"))
  }
  expect_identical(nrow(mine_frequent_itemsets(event_df, min_count = 3L)), 3L)
  expect_identical(nrow(mine_frequent_itemsets(event_df, min_count = 2L)), 4L)
})

test_that("R-only: min_count is used as a count, without Julia's fraction round trip", {
  # 100 records, 7 holding X+Y. Julia passes min_count = 7 as 7/100 and
  # RuleMiner takes ceil(7/100 * 100) = 8, dropping X+Y; here it is kept.
  ev <- events_from_records(list(list(items = c("X", "Y"), n = 7L),
                                 list(items = c("P", "Q"), n = 93L)))
  fi <- mine_frequent_itemsets(ev, min_count = 7L)
  expect_true("X+Y" %in% set_keys(fi[["Itemset"]]))
  expect_identical(nrow(mine_frequent_itemsets(ev, min_count = 8L)), 1L)
})

test_that("mine_association_rules", {
  rules <- mine_association_rules(event_df, min_support = 0.1, min_confidence = 0.3)
  expect_identical(names(rules), c("LHS", "RHS", "Support", "Confidence", "Coverage", "Lift",
                                   "N", "Length"))
  expect_true(all(lengths(rules[["LHS"]]) >= 1L))
  # Julia: 7 rules; Item01 => GB1 (confidence 3/13) falls below 0.3.
  expect_identical(nrow(rules), 7L)
  expect_false(any(set_keys(rules[["LHS"]]) == "Item01" & rules[["RHS"]] == "GB1"))
})

test_that("mine_association_rules metrics match the reference exactly", {
  rules <- mine_association_rules(event_df, min_count = 2L, min_confidence = 0, max_length = 2L)
  got <- rules[order(set_keys(rules[["LHS"]]), rules[["RHS"]]), ]
  expect_identical(paste(set_keys(got[["LHS"]]), got[["RHS"]], sep = ">"),
                   c("GA1>Item05", "GB1>Item01", "Item01>GB1", "Item01>Item02",
                     "Item02>Item01", "Item03>Item04", "Item04>Item03", "Item05>GA1"))
  expect_identical(got[["N"]], c(5L, 3L, 3L, 10L, 10L, 2L, 2L, 5L))
  expect_identical(got[["Length"]], rep(2L, 8L))
  expect_equal(got[["Confidence"]], c(1, 1, 0.23076923076923078, 0.76923076923076927, 1, 1, 1, 1),
               tolerance = 1e-15)
  expect_equal(got[["Coverage"]], c(0.25, 0.15, 0.65, 0.65, 0.5, 0.1, 0.1, 0.25), tolerance = 1e-15)
  expect_equal(got[["Lift"]], c(4, rep(1.5384615384615383, 4), 9.9999999999999982,
                                9.9999999999999982, 4), tolerance = 1e-15)
  expect_equal(got[["Support"]], c(0.25, 0.15, 0.15, 0.5, 0.5, 0.1, 0.1, 0.25), tolerance = 1e-15)

  # Julia at min_count = 5, min_confidence = 0.5: the four rules on the N >= 5 pairs.
  strong <- mine_association_rules(event_df, min_count = 5L, min_confidence = 0.5)
  expect_identical(sort(paste(set_keys(strong[["LHS"]]), strong[["RHS"]], sep = ">"),
                        method = "radix"),
                   c("GA1>Item05", "Item01>Item02", "Item02>Item01", "Item05>GA1"))
})

test_that("R-only: rules of three items are counted from the records, not from arules", {
  ev <- events_from_records(list(list(items = c("A", "B", "C"), n = 6L),
                                 list(items = c("A", "B"), n = 4L),
                                 list(items = c("C", "D"), n = 10L)))
  rules <- mine_association_rules(ev, min_count = 6L, min_confidence = 0)
  ab_c <- rules[set_keys(rules[["LHS"]]) == "A+B" & rules[["RHS"]] == "C", ]
  expect_identical(ab_c[["N"]], 6L)
  expect_identical(ab_c[["Length"]], 3L)
  expect_equal(ab_c[["Confidence"]], 6 / 10)
  expect_equal(ab_c[["Lift"]], (6 / 20) / ((10 / 20) * (16 / 20)))
  expect_false(any(lengths(rules[["LHS"]]) == 0L))
})

test_that("empty transactions", {
  singles <- data.frame(id = 1:3, item = c("A", "B", "C"))
  fi <- mine_frequent_itemsets(singles)
  expect_identical(nrow(fi), 0L)
  expect_identical(names(fi), c("Itemset", "Support", "N", "Length"))
  ar <- mine_association_rules(singles)
  expect_identical(nrow(ar), 0L)
  expect_identical(names(ar), c("LHS", "RHS", "Support", "Confidence", "Coverage", "Lift",
                                "N", "Length"))
})

test_that("R-only: invalid mining arguments raise", {
  expect_error(mine_frequent_itemsets(event_df, min_count = 2.5), "min_count")
  expect_error(mine_frequent_itemsets(event_df, min_support = 0), "min_support")
  expect_error(mine_frequent_itemsets(event_df, max_length = 1L), "max_length")
  expect_error(mine_association_rules(event_df, min_confidence = 1.5), "min_confidence")
  expect_error(mine_frequent_itemsets(transform(event_df, id = replace(id, 1L, NA))),
               "contains NA")
})
