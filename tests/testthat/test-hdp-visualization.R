# HDP plots. One case per Julia @testset under "HDP visualization", plus
# R-only cases on what the butterfly draws.

hdp_fixture <- function() {
  set.seed(42)
  hdp_clustering(make_test_event_df(), k_max = 4L, n_init = 2L)
}

test_that("plot_hdp_stick_weights", {
  expect_s3_class(plot_hdp_stick_weights(hdp_fixture()), "ggplot")
})

test_that("plot_hdp_class_profiles", {
  expect_s3_class(plot_hdp_class_profiles(hdp_fixture()), "ggplot")
})

test_that("plot_hdp_sharing_heatmap", {
  expect_s3_class(plot_hdp_sharing_heatmap(hdp_fixture()), "ggplot")
})

test_that("plot_hdp_cluster_butterfly", {
  out <- plot_hdp_cluster_butterfly(hdp_fixture(), left = "A", right = "B")
  expect_s3_class(out[["plot"]], "ggplot")
})

test_that("plot_hdp_cluster_butterfly requires 2 groups", {
  ev <- make_test_event_df()
  set.seed(7)
  single <- hdp_clustering(ev[ev[["Group"]] == "B", ], k_max = 3L, n_init = 1L)
  expect_error(plot_hdp_cluster_butterfly(single, left = "B", right = "B"), "exactly 2 groups")
})

test_that("all three respect beta_threshold", {
  res <- hdp_fixture()
  expect_error(plot_hdp_stick_weights(res, beta_threshold = 2), "No active clusters")
  expect_error(plot_hdp_class_profiles(res, beta_threshold = 2), "No active clusters")
  expect_error(plot_hdp_sharing_heatmap(res, beta_threshold = 2), "No active clusters")
})

test_that("R-only: the butterfly puts each group on the side named for it", {
  res <- hdp_fixture()
  ab <- plot_hdp_cluster_butterfly(res, left = "A", right = "B")[["bars"]]
  ba <- plot_hdp_cluster_butterfly(res, left = "B", right = "A")[["bars"]]
  expect_identical(unique(ab[["side"]][ab[["group"]] == "A"]), "left")
  expect_identical(unique(ab[["side"]][ab[["group"]] == "B"]), "right")
  expect_identical(unique(ba[["side"]][ba[["group"]] == "A"]), "right")
  # The bar lengths are the groups' weights, whichever side they are drawn on.
  h <- res[["hdp"]]
  a_rows <- ab[ab[["group"]] == "A", ]
  expect_equal(a_rows[["pi"]], h[["pi_mean"]][1L, a_rows[["cluster"]]])
  # Colour follows the category, not the side.
  expect_identical(ab[["fill"]][order(ab[["group"]], ab[["cluster"]])],
                   ba[["fill"]][order(ba[["group"]], ba[["cluster"]])])
})

test_that("R-only: the butterfly refuses positions or unknown names", {
  res <- hdp_fixture()
  expect_error(plot_hdp_cluster_butterfly(res), "left")
  expect_error(plot_hdp_cluster_butterfly(res, left = "A", right = "A"), "must name")
  expect_error(plot_hdp_cluster_butterfly(res, left = "a", right = "B"), "must name")
})

test_that("R-only: profile item order groups items by their peak cluster", {
  theta <- rbind(c(0.9, 0.01, 0.3, 0.7), c(0.1, 0.02, 0.8, 0.2))
  # Items 1 and 4 peak in cluster 1 (0.9 then 0.7); item 3 in cluster 2;
  # item 2 never reaches the threshold.
  expect_identical(.filter_and_order_items(theta, 0.05), c(1L, 4L, 3L))
})
