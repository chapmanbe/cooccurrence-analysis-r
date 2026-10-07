# Group-stratified co-occurrence networks.
#
# One network per value of a group column. Some items can occur in only one
# group (the caller knows which; this package does not). Within every other
# group such an item is impossible, and leaving its stray events in would
# make it a vocabulary item there -- scored as absent in every record that
# lacks it. The caller therefore supplies `exclusive_items`, and this file
# makes each way of getting that argument wrong raise instead of passing
# silently: it has no default, its keys must be group values, and its items
# must be items of the data.

#' Validate `exclusive_items` against the data's groups and items
#'
#' Returns the list unchanged, or raises / warns:
#' - not a list with every element named by a group value: error;
#' - an item listed under more than one group: error;
#' - no listed item occurs anywhere in the data: error (a list built for a
#'   different item vocabulary matches nothing, which would switch the
#'   exclusion off without a sound);
#' - some listed items do not occur: a warning naming them.
#' @keywords internal
.check_exclusive_items <- function(exclusive_items, groups, items, group_col) {
  if (is.null(exclusive_items)) return(NULL)
  if (!is.list(exclusive_items) || length(exclusive_items) == 0L) {
    stop("exclusive_items must be NULL (no exclusion) or a non-empty named list ",
         "of group -> items", call. = FALSE)
  }
  keys <- names(exclusive_items)
  if (is.null(keys) || anyNA(keys) || any(!nzchar(keys)) || anyDuplicated(keys)) {
    stop("exclusive_items must be named by group value, each name once", call. = FALSE)
  }
  unknown <- setdiff(keys, groups)
  if (length(unknown) > 0L) {
    stop(sprintf("exclusive_items key(s) %s not a value of '%s' (values: %s)",
                 paste0("'", unknown, "'", collapse = ", "), group_col,
                 paste(groups, collapse = ", ")), call. = FALSE)
  }
  for (k in keys) {
    v <- exclusive_items[[k]]
    if (!is.character(v) || anyNA(v)) {
      stop(sprintf("exclusive_items[['%s']] must be a character vector without NA", k),
           call. = FALSE)
    }
  }
  listed <- unlist(lapply(exclusive_items, unique), use.names = FALSE)
  shared <- unique(listed[duplicated(listed)])
  if (length(shared) > 0L) {
    stop(sprintf("item(s) listed as exclusive to more than one group: %s",
                 paste(shared, collapse = ", ")), call. = FALSE)
  }
  absent <- setdiff(listed, items)
  if (length(absent) == length(listed)) {
    stop(sprintf(paste0("none of the exclusive_items occurs in the data's item vocabulary ",
                        "(listed: %s); the list was probably built for another vocabulary"),
                 paste(listed, collapse = ", ")), call. = FALSE)
  }
  if (length(absent) > 0L) {
    warning(sprintf(paste0("exclusive_items not in the data's item vocabulary, so they ",
                           "exclude nothing: %s"), paste(absent, collapse = ", ")),
            call. = FALSE)
  }
  exclusive_items
}

#' Co-occurrence networks per group, excluding other groups' exclusive items
#'
#' Splits `event_df` (columns `id`, `item`, and `group_col`) by group value
#' and fits, per group, [build_cooccurrence_network()], [detect_communities()]
#' and [compute_network_metrics()].
#'
#' `exclusive_items` has **no default**: pass `NULL` to say there is nothing
#' to exclude, or a list named by group value whose elements are the items
#' that can occur only in that group. In the stratum of every *other* group,
#' events of those items are dropped before counting: the item is not part of
#' that stratum's vocabulary, and a record left with no events is not one of
#' its records. See `.check_exclusive_items()` for what raises.
#'
#' Returns a list named by group value, in C-locale sorted order, each
#' element `list(network, communities, metrics, excluded_items,
#' n_excluded_events)`: `excluded_items` is the sorted list of items this
#' stratum was set to drop (the other groups' exclusive items), and
#' `n_excluded_events` how many of its events they matched and removed. The
#' list lets a caller check that the right exclusion was requested even on
#' data where it removes nothing; only the count reflects the data. Read the
#' result by name
#' (`result[["A"]]`), never by position. A missing group column, or `NA` in
#' it, raises: dropping those rows would shrink every count silently.
#'
#' `seed` is passed to [detect_communities()] unchanged for every group, so a
#' group's partition does not depend on which groups precede it.
#' @export
stratified_network_analysis <- function(event_df, exclusive_items, group_col = "Group",
                                        weight_metric = "lift", min_count = 30L,
                                        alpha = 0.05, test = "fisher", correction = "bh",
                                        community_method = "leiden", timing_filter = "all", seed) {
  if (missing(exclusive_items)) {
    stop(paste0("exclusive_items is required: pass NULL for no exclusion, or a list ",
                "named by group of the items only that group can hold"), call. = FALSE)
  }
  if (missing(seed)) {
    stop(paste0("stratified_network_analysis: seed is required: pass a whole number for ",
                "reproducible partitions (every group uses it), or NULL to draw from the ",
                "current RNG stream"), call. = FALSE)
  }
  .check_seed(seed, "stratified_network_analysis")
  .check_options(test, correction, timing_filter)
  if (!(group_col %in% names(event_df))) {
    stop(sprintf("event table has no '%s' column (columns: %s)", group_col,
                 paste(names(event_df), collapse = ", ")), call. = FALSE)
  }
  g <- event_df[[group_col]]
  if (anyNA(g)) {
    stop(sprintf("group column '%s' contains NA; every event needs a group", group_col),
         call. = FALSE)
  }
  g <- as.character(g)
  groups <- .sort_c(unique(g))
  items_chr <- as.character(event_df[["item"]])
  exclusive_items <- .check_exclusive_items(exclusive_items, groups,
                                            unique(items_chr), group_col)

  out <- lapply(groups, function(grp) {
    others <- setdiff(names(exclusive_items), grp)
    dropped <- .sort_c(as.character(unique(unlist(exclusive_items[others],
                                                  use.names = FALSE))))
    in_group <- g == grp
    excluded <- in_group & items_chr %in% dropped
    stratum <- event_df[in_group & !excluded, c("id", "item"), drop = FALSE]
    rownames(stratum) <- NULL
    net <- build_cooccurrence_network(stratum, weight_metric = weight_metric,
                                      min_count = min_count, alpha = alpha,
                                      test = test, correction = correction,
                                      timing_filter = timing_filter)
    list(network = net,
         communities = detect_communities(net, method = community_method, seed = seed),
         metrics = compute_network_metrics(net),
         excluded_items = dropped,
         n_excluded_events = sum(excluded))
  })
  names(out) <- groups
  out
}
