#' Whole-population network pipeline: build, detect communities, measure
#'
#' Returns a list: `network`, `communities`, `metrics`. The Julia reference's
#' `run_network_pipeline` also builds per-group networks by default
#' (`stratify_by_group = true`) and compares every pair of them with
#' `compare_networks`. The per-group networks are ported
#' ([stratified_network_analysis()]); `compare_networks` is not, since no
#' pipeline rule calls it. So the default here matches Julia's (`TRUE`) and
#' raises, rather than silently returning a result without the comparisons;
#' pass `stratify_by_group = FALSE` for the pooled network, and call
#' [stratified_network_analysis()] for the per-group ones. `test` and
#' `correction` go to [build_cooccurrence_network()]; `community_method` and `seed` to
#' [detect_communities()].
#' @export
run_network_pipeline <- function(event_df, weight_metric = "lift", min_count = 30L,
                                 alpha = 0.05, community_method = "leiden",
                                 stratify_by_group = TRUE, test = "fisher", correction = "bh",
                                 seed) {
  if (!identical(stratify_by_group, FALSE)) {
    stop(paste0("stratify_by_group = TRUE needs compare_networks, which exists in the ",
                "Julia reference but is not ported to R; pass stratify_by_group = FALSE ",
                "and call stratified_network_analysis() for the per-group networks"),
         call. = FALSE)
  }
  if (missing(seed)) {
    stop(paste0("run_network_pipeline: seed is required: pass a whole number for a ",
                "reproducible partition, or NULL to draw from the current RNG stream"),
         call. = FALSE)
  }
  net <- build_cooccurrence_network(event_df, weight_metric = weight_metric,
                                    min_count = min_count, alpha = alpha, test = test,
                                    correction = correction)
  list(network = net,
       communities = detect_communities(net, method = community_method, seed = seed),
       metrics = compute_network_metrics(net))
}
