# Item-domain wrappers for the HDP Bernoulli mixture. R port of the Julia
# reference's hdp_wrappers.jl.

#' Per-group record-by-item matrix
#'
#' Within each group (sorted, C-locale order), one row per record (`id`)
#' holding at least `min_items` distinct items; columns are the sorted union
#' of items over the kept records. A record is identified by its id within
#' its group, as the reference builds each group's transactions separately.
#' Rows come group by group, records in order of first appearance.
#'
#' Returns a list: `x` (logical matrix, items as column names), `ids`,
#' `group` (integer index into `group_labels`), `group_labels`. Raises on
#' `NA` in `id`, `item` or the group column.
#' @export
record_item_matrix <- function(event_df, group_col = "Group", min_items = 2L) {
  for (col in c("id", "item", group_col)) {
    if (!(col %in% names(event_df))) {
      stop(sprintf("Column %s not found in event_df", col), call. = FALSE)
    }
    if (anyNA(event_df[[col]])) {
      stop(sprintf("event_df column '%s' contains NA; every event needs one", col),
           call. = FALSE)
    }
  }
  grp_chr <- as.character(event_df[[group_col]])
  group_labels <- .sort_c(unique(grp_chr))
  if (length(group_labels) == 0L) stop(sprintf("No groups found in column %s", group_col),
                                       call. = FALSE)
  g <- match(grp_chr, group_labels)
  ids <- event_df[["id"]]
  uid <- unique(ids)
  id_code <- match(ids, uid)
  item_chr <- as.character(event_df[["item"]])

  # A record is (group, id). Rows: group-major, then first appearance.
  rec_key <- (g - 1) * length(uid) + id_code
  rec_levels <- unique(rec_key)
  rec_levels <- rec_levels[order(ceiling(rec_levels / length(uid)), method = "radix")]
  rec <- match(rec_key, rec_levels)

  items_all <- .sort_c(unique(item_chr))
  code <- match(item_chr, items_all)
  inc <- .incidence(rec, code, length(rec_levels), length(items_all))
  kept_rec <- which(Matrix::rowSums(inc) >= min_items)
  if (length(kept_rec) == 0L) {
    return(list(x = matrix(FALSE, 0L, 0L), ids = ids[0L], group = integer(0),
                group_labels = group_labels))
  }
  inc <- inc[kept_rec, , drop = FALSE]
  held <- which(Matrix::colSums(inc) > 0)
  x <- as.matrix(inc[, held, drop = FALSE]) > 0
  dimnames(x) <- list(NULL, items_all[held])

  first_event <- match(kept_rec, rec)
  list(x = x, ids = ids[first_event], group = g[first_event], group_labels = group_labels)
}

#' Fit an HDP Bernoulli mixture jointly across groups
#'
#' One joint fit over a shared atom pool with per-group mixing weights, so it
#' needs no exclusive-item filtering: a group-exclusive item simply has low
#' theta in clusters dominated by the other groups. Groups are the sorted
#' values of `group_col`; records are those with at least `min_items`
#' distinct items (see [record_item_matrix()]). Fit arguments pass to
#' [fit_hdp_bernoulli_mixture()]; restarts draw from R's global random
#' stream, so seed it first -- or pass `init`, an explicit start keyed by
#' record: a data frame with `group`, `id`, and `k_max` numeric columns, one
#' row per fitted record (raises otherwise; needs `n_init = 1`).
#' `timing_filter` other than `"all"` is not
#' ported and raises.
#'
#' Returns an `hdp_clustering_result` list:
#' - `hdp`: the [fit_hdp_bernoulli_mixture()] result, with `group_labels`
#'   set to the group values.
#' - `item_names`: the items, in column order of `theta`.
#' - `class_profiles`: `cluster`, `beta_weight`, then one column per item
#'   (theta), all rounded to 4 places.
#' - `group_profiles`: `cluster`, then `pi_<group>` per group, rounded to 4.
#' - `record_assignments`: `id`, `group`, `assignment`, `max_responsibility`.
#' - `n_records_per_group`: named by group.
#'
#' Raises if an item name collides with `cluster`, `beta_weight` or a
#' `pi_<group>` column, or if a group has no multi-item records.
#' @export
hdp_clustering <- function(event_df, group_col = "Group", min_items = 2L, k_max = 20L,
                           alpha = 1, gamma = 1, prior = "flat", alpha_prior = 1,
                           beta_prior = 1, concentration = 10, floor = 1, max_iter = 300L,
                           tol = 1e-5, n_init = 3L, timing_filter = "all", init = NULL) {
  if (!(group_col %in% names(event_df))) {
    stop(sprintf("Column %s not found in event_df", group_col), call. = FALSE)
  }
  if (!identical(timing_filter, "all")) {
    stop(paste0("timing_filter = '", format(timing_filter), "' exists in the Julia reference ",
                "but is not ported to R; only 'all' is"), call. = FALSE)
  }
  m <- record_item_matrix(event_df, group_col = group_col, min_items = min_items)
  labels <- m[["group_labels"]]
  if (nrow(m[["x"]]) == 0L) stop("No multi-item records found after filtering", call. = FALSE)
  items <- colnames(m[["x"]])
  reserved <- c("cluster", "beta_weight", paste0("pi_", labels))
  clash <- intersect(items, reserved)
  if (length(clash) > 0L) {
    stop(sprintf(paste0("Item name(s) %s collide with reserved profile metadata columns ",
                        "(cluster, beta_weight, pi_<group>); rename the item(s)."),
                 paste(clash, collapse = ", ")), call. = FALSE)
  }
  per_group <- stats::setNames(tabulate(m[["group"]], nbins = length(labels)), labels)
  if (any(per_group == 0L)) {
    stop(sprintf("Group %s has no multi-item records",
                 paste(labels[per_group == 0L], collapse = ", ")), call. = FALSE)
  }

  x <- m[["x"]]
  dimnames(x) <- NULL
  init_m <- if (is.null(init)) NULL else .init_for_records(init, m, k_max)
  hdp <- fit_hdp_bernoulli_mixture(x, m[["group"]], length(labels), k_max = k_max,
                                   alpha = alpha, gamma = gamma, prior = prior,
                                   alpha_prior = alpha_prior, beta_prior = beta_prior,
                                   concentration = concentration, floor = floor,
                                   max_iter = max_iter, tol = tol, n_init = n_init,
                                   init = init_m)
  hdp[["group_labels"]] <- labels

  k <- hdp[["k_max"]]
  profiles <- data.frame(cluster = seq_len(k), beta_weight = round_digits(hdp[["beta_mean"]], 4))
  theta <- round_digits(hdp[["theta"]], 4)
  for (d in seq_along(items)) profiles[[items[[d]]]] <- theta[, d]

  gprofiles <- data.frame(cluster = seq_len(k))
  for (j in seq_along(labels)) {
    gprofiles[[paste0("pi_", labels[[j]])]] <- round_digits(hdp[["pi_mean"]][j, ], 4)
  }

  r <- hdp[["responsibilities"]]
  assignments <- data.frame(id = m[["ids"]], group = labels[m[["group"]]],
                            assignment = hdp[["assignments"]],
                            max_responsibility = r[cbind(seq_len(nrow(r)), hdp[["assignments"]])])

  structure(list(hdp = hdp, item_names = items, class_profiles = profiles,
                 group_profiles = gprofiles, record_assignments = assignments,
                 n_records_per_group = per_group),
            class = "hdp_clustering_result")
}

# An explicit starting point keyed by record: `init` has `group`, `id`, then
# `k_max` numeric columns (one per component, in order). Returns it as a
# matrix in the row order of `m` (a record_item_matrix() result). Raises
# unless every record has exactly one row, so a start can never be applied
# to the wrong record or silently padded.
.init_for_records <- function(init, m, k_max) {
  if (!is.data.frame(init) || !all(c("group", "id") %in% names(init))) {
    stop("init must be a data frame with 'group' and 'id' columns", call. = FALSE)
  }
  comp <- setdiff(names(init), c("group", "id"))
  if (length(comp) != k_max) {
    stop(sprintf("init has %d component columns; k_max is %d", length(comp), k_max),
         call. = FALSE)
  }
  key <- function(g, id) paste(g, sprintf("%.17g", as.numeric(id)), sep = "\r")
  init_key <- key(as.character(init[["group"]]), init[["id"]])
  if (anyDuplicated(init_key)) stop("init has a record more than once", call. = FALSE)
  rows <- match(key(m[["group_labels"]][m[["group"]]], m[["ids"]]), init_key)
  if (anyNA(rows) || nrow(init) != length(rows)) {
    stop(sprintf("init does not cover exactly the %d fitted records (%d missing, %d rows given)",
                 length(rows), sum(is.na(rows)), nrow(init)), call. = FALSE)
  }
  out <- as.matrix(init[rows, comp, drop = FALSE])
  dimnames(out) <- NULL
  out
}

#' Categorize each HDP cluster by its presence across groups
#'
#' A cluster is present in group j when `pi_mean[j, k] > presence_threshold`.
#' Categories: `"negligible"` (`beta_mean <= beta_threshold`, or present in
#' no group); `"universal"` (present in every group, max/min weight ratio
#' below `universal_ratio`); `"both_present_but_unequal"` (present in every
#' group, ratio at or above it); otherwise named for the present group(s),
#' lowercased: `"group_<a>_only"`, `"group_<a>_and_<b>_only"`.
#'
#' Returns a data frame: `cluster`, `beta_weight`, `pi_<group>` per group
#' (rounded to 4 places; the tests use unrounded values), `category`.
#' @export
hdp_cluster_categorization <- function(result, presence_threshold = 0.05,
                                       beta_threshold = 0.01, universal_ratio = 2,
                                       ratio_floor = 1e-8) {
  if (!inherits(result, "hdp_clustering_result")) {
    stop("result must be an hdp_clustering_result", call. = FALSE)
  }
  hdp <- result[["hdp"]]
  labels <- hdp[["group_labels"]]
  pi_mean <- hdp[["pi_mean"]]
  beta_mean <- hdp[["beta_mean"]]
  category <- vapply(seq_len(hdp[["k_max"]]), function(k) {
    pis <- pi_mean[, k]
    present <- pis > presence_threshold
    if (beta_mean[[k]] <= beta_threshold || !any(present)) return("negligible")
    if (all(present)) {
      ratio <- max(pis) / max(min(pis), ratio_floor)
      return(if (ratio < universal_ratio) "universal" else "both_present_but_unequal")
    }
    paste0("group_", paste(tolower(labels[present]), collapse = "_and_"), "_only")
  }, character(1))
  out <- data.frame(cluster = seq_len(hdp[["k_max"]]), beta_weight = round_digits(beta_mean, 4))
  for (j in seq_along(labels)) out[[paste0("pi_", labels[[j]])]] <- round_digits(pi_mean[j, ], 4)
  out[["category"]] <- category
  out
}

#' Active clusters: those whose global weight exceeds `beta_threshold`
#' @export
hdp_active_clusters <- function(result, beta_threshold = 0.01) {
  which(result[["hdp"]][["beta_mean"]] > beta_threshold)
}
