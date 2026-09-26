#' Whole-population network pipeline: build, detect communities, measure
#'
#' Returns a list: `network`, `communities`, `metrics`. The Julia reference's
#' `run_network_pipeline` also builds per-group networks by default
#' (`stratify_by_group = true`); that half arrives with the stratified
#' network (port step 4), so here `stratify_by_group` must be `FALSE` and
#' `TRUE` raises rather than silently returning only the pooled network.
#' @export
run_network_pipeline <- function(event_df, weight_metric = "lift", min_count = 30L,
                                 alpha = 0.05, community_method = "louvain",
                                 stratify_by_group = FALSE) {
  if (!identical(stratify_by_group, FALSE)) {
    stop(paste0("stratify_by_group = TRUE is not ported yet: group-stratified networks ",
                "arrive at port step 4; pass stratify_by_group = FALSE"), call. = FALSE)
  }
  net <- build_cooccurrence_network(event_df, weight_metric = weight_metric,
                                    min_count = min_count, alpha = alpha)
  list(network = net,
       communities = detect_communities(net, method = community_method),
       metrics = compute_network_metrics(net))
}
