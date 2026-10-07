# HDP plots (ggplot2). R port of the HDP half of the Julia reference's
# bayesian_visualization.jl. Each returns a ggplot object; saving is the
# caller's job. Figures are judged by eye -- the data behind each one is
# what gets compared -- so layout follows the reference only loosely (e.g.
# the stick-weight figure carries each cluster's category and top items in
# its axis labels rather than in a second panel).

.check_hdp_result <- function(result) {
  if (!inherits(result, "hdp_clustering_result")) {
    stop("result must be an hdp_clustering_result", call. = FALSE)
  }
  invisible(result)
}

.active_or_stop <- function(result, beta_threshold, inclusive = FALSE) {
  b <- result[["hdp"]][["beta_mean"]]
  active <- if (inclusive) which(b >= beta_threshold) else which(b > beta_threshold)
  if (length(active) == 0L) {
    stop(sprintf("No active clusters above beta_threshold=%s", format(beta_threshold)),
         call. = FALSE)
  }
  active
}

# The n items with the highest theta in cluster k, highest first (ties by
# item order).
.top_items <- function(result, k, n) {
  th <- result[["hdp"]][["theta"]][k, ]
  idx <- order(-th, method = "radix")[seq_len(min(n, length(th)))]
  result[["item_names"]][idx]
}

#' Items to show in a cluster-profile heatmap, and their order
#'
#' Keeps items whose largest theta over the given clusters (rows of
#' `theta`) is at least `threshold`, ordered by the cluster where each peaks,
#' then by that peak, descending. Returns the kept column indices.
#' @keywords internal
.filter_and_order_items <- function(theta, threshold) {
  peak <- vapply(seq_len(ncol(theta)), function(d) max(theta[, d]), numeric(1))
  keep <- which(peak >= threshold)
  if (length(keep) == 0L) return(integer(0))
  peak_cluster <- vapply(keep, function(d) which.max(theta[, d]), integer(1))
  keep[order(peak_cluster, -peak[keep], method = "radix")]
}

#' Global and per-group stick weights of the active clusters
#'
#' One bar per group (`pi[j, k]`, in `group_labels` order) and one for the
#' global weight (`beta[k]`, gray) at each active cluster (`beta > beta_threshold`).
#' Axis labels carry each cluster's category and two highest-theta items.
#' Raises if no cluster is active.
#' @export
plot_hdp_stick_weights <- function(result, beta_threshold = 0.01) {
  .check_hdp_result(result)
  hdp <- result[["hdp"]]
  active <- .active_or_stop(result, beta_threshold)
  labels <- hdp[["group_labels"]]
  cats <- hdp_cluster_categorization(result, beta_threshold = beta_threshold)[["category"]]
  axis_labels <- vapply(active, function(k) {
    sprintf("C%d\n%s\n%s", k, gsub("_", " ", cats[[k]]),
            paste(.top_items(result, k, 2L), collapse = "\n"))
  }, character(1))
  series <- c(labels, "Global β")
  per_group <- lapply(seq_along(labels), function(j) {
    data.frame(cluster = active, series = labels[[j]], weight = hdp[["pi_mean"]][j, active])
  })
  global <- data.frame(cluster = active, series = "Global β",
                       weight = hdp[["beta_mean"]][active])
  df <- do.call(rbind, c(per_group, list(global)))
  df[["series"]] <- factor(df[["series"]], levels = series)
  df[["cluster"]] <- factor(df[["cluster"]], levels = active)
  hue <- (seq_along(labels) - 1L) %% length(.community_palette) + 1L
  fill <- stats::setNames(c(.community_palette[hue], "gray40"), series)
  ggplot2::ggplot(df, ggplot2::aes(x = .data$cluster, y = .data$weight, fill = .data$series)) +
    ggplot2::geom_col(position = ggplot2::position_dodge(width = 0.8), width = 0.8) +
    ggplot2::scale_fill_manual(values = fill, name = NULL) +
    ggplot2::scale_x_discrete(labels = stats::setNames(axis_labels, active)) +
    ggplot2::labs(title = "HDP Stick Weights: Global β and Per-Group π",
                  x = "Cluster", y = "Weight") +
    ggplot2::theme_minimal(base_size = 10) +
    ggplot2::theme(legend.position = "top")
}

#' Heatmap of posterior mean theta for the active clusters
#'
#' Rows: active clusters (`beta > beta_threshold`), labelled with category
#' and global weight, the heaviest at the bottom as in the reference.
#' Columns: items whose theta reaches `prob_threshold` in some active
#' cluster, grouped by the cluster where they peak. Cells with theta >= 0.15
#' carry the value. Raises if no cluster is active or no item qualifies.
#' @export
plot_hdp_class_profiles <- function(result, beta_threshold = 0.01, prob_threshold = 0.05,
                                    item_label = "Item") {
  .check_hdp_result(result)
  hdp <- result[["hdp"]]
  active <- .active_or_stop(result, beta_threshold)
  theta <- hdp[["theta"]][active, , drop = FALSE]
  cols <- .filter_and_order_items(theta, prob_threshold)
  if (length(cols) == 0L) {
    stop(sprintf("No items above prob_threshold=%s", format(prob_threshold)), call. = FALSE)
  }
  cats <- hdp_cluster_categorization(result, beta_threshold = beta_threshold)[["category"]]
  row_labels <- sprintf("C%d [%s]  β=%s", active, gsub("_", " ", cats[active]),
                        format(round_digits(hdp[["beta_mean"]][active], 3)))
  items <- result[["item_names"]][cols]
  df <- data.frame(row = rep(seq_along(active), times = length(cols)),
                   item = factor(rep(items, each = length(active)), levels = items),
                   theta = as.vector(theta[, cols, drop = FALSE]))
  df[["label"]] <- ifelse(df[["theta"]] >= 0.15, format(round_digits(df[["theta"]], 2)), "")
  df[["text_colour"]] <- ifelse(df[["theta"]] > 0.6, "white", "black")
  ggplot2::ggplot(df, ggplot2::aes(x = .data$item, y = .data$row)) +
    ggplot2::geom_tile(ggplot2::aes(fill = .data$theta)) +
    ggplot2::geom_text(ggplot2::aes(label = .data$label, colour = .data$text_colour),
                       size = 2.2) +
    ggplot2::scale_colour_identity() +
    ggplot2::scale_fill_distiller(palette = "YlOrRd", direction = 1, limits = c(0, 1),
                                  name = "E_q[θ[k,d]]") +
    ggplot2::scale_y_continuous(breaks = seq_along(active), labels = row_labels,
                                expand = c(0, 0)) +
    ggplot2::labs(title = "HDP Cluster Profiles — Posterior Mean θ[k,d]",
                  x = item_label, y = "Cluster") +
    ggplot2::theme_minimal(base_size = 9) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 60, hjust = 1))
}

#' Cluster-by-group heatmap of the per-group mixing weights
#'
#' Rows: active clusters (`beta > beta_threshold`) labelled with their three
#' highest-theta items and their category; columns: groups in
#' `group_labels` order. Cell text is `pi[j, k]` to 3 places. Raises if no
#' cluster is active.
#' @export
plot_hdp_sharing_heatmap <- function(result, beta_threshold = 0.01, group_label = "Group") {
  .check_hdp_result(result)
  hdp <- result[["hdp"]]
  active <- .active_or_stop(result, beta_threshold)
  labels <- hdp[["group_labels"]]
  cats <- hdp_cluster_categorization(result, beta_threshold = beta_threshold)[["category"]]
  short <- gsub("_", " ", sub("both_present_but_unequal", "unequal", cats[active], fixed = TRUE))
  row_labels <- vapply(seq_along(active), function(i) {
    top <- paste(.top_items(result, active[[i]], 3L), collapse = ", ")
    sprintf("C%d: %s  [%s]", active[[i]], top, short[[i]])
  }, character(1))
  pi_act <- hdp[["pi_mean"]][, active, drop = FALSE]
  top <- max(pi_act) + 1e-6
  df <- data.frame(row = rep(seq_along(active), each = length(labels)),
                   group = factor(rep(labels, times = length(active)), levels = labels),
                   pi = as.vector(pi_act))
  df[["label"]] <- format(round_digits(df[["pi"]], 3))
  df[["text_colour"]] <- ifelse(df[["pi"]] > 0.3 * max(pi_act), "white", "black")
  ggplot2::ggplot(df, ggplot2::aes(x = .data$group, y = .data$row)) +
    ggplot2::geom_tile(ggplot2::aes(fill = .data$pi)) +
    ggplot2::geom_text(ggplot2::aes(label = .data$label, colour = .data$text_colour), size = 3) +
    ggplot2::scale_colour_identity() +
    ggplot2::scale_fill_distiller(palette = "YlOrRd", direction = 1, limits = c(0, top),
                                  name = "π[j,k]") +
    ggplot2::scale_y_continuous(breaks = seq_along(active), labels = row_labels,
                                expand = c(0, 0)) +
    ggplot2::labs(title = "HDP Cross-Strata Cluster Sharing (π[j,k])",
                  x = group_label, y = "Cluster") +
    ggplot2::theme_minimal(base_size = 9)
}

#' Butterfly chart of two groups' mixing weights
#'
#' For each active cluster (`beta >= beta_threshold`, as in the reference),
#' `left`'s weight `pi` extends left and `right`'s extends right; clusters
#' run top to bottom by global weight. `left` and `right` are group *names*,
#' required, and must be the result's two groups: nothing here falls back to
#' a position. Bars are coloured by category: universal blue, the first
#' group's (in `group_labels` order) only-cluster orange, the second's
#' green, both-present-but-unequal reddish purple, negligible gray.
#'
#' Returns a list: `plot`, and `bars`, one row per drawn bar as read back
#' from `ggplot2::ggplot_build()` -- `cluster`, `group`, `side` (`"left"` or
#' `"right"`, from the bar's drawn extent), `pi`, `fill` -- so a caller
#' records what the figure shows. Raises unless the result has exactly two
#' groups.
#' @export
plot_hdp_cluster_butterfly <- function(result, left, right, beta_threshold = 0.01) {
  .check_hdp_result(result)
  hdp <- result[["hdp"]]
  labels <- hdp[["group_labels"]]
  if (length(labels) != 2L) {
    stop(sprintf("Butterfly chart requires exactly 2 groups; got %d", length(labels)),
         call. = FALSE)
  }
  if (!(is.character(left) && is.character(right) && length(left) == 1L &&
          length(right) == 1L && setequal(c(left, right), labels) && left != right)) {
    stop(sprintf("left and right must name the two groups (%s); got %s, %s",
                 paste(labels, collapse = ", "), format(left), format(right)), call. = FALSE)
  }
  active <- .active_or_stop(result, beta_threshold, inclusive = TRUE)
  active <- active[order(-hdp[["beta_mean"]][active], method = "radix")]
  cats <- hdp_cluster_categorization(result, beta_threshold = beta_threshold)[["category"]]
  pal <- .community_palette
  cat_fill <- c(universal = pal[[1L]], negligible = "#D3D3D3",
                both_present_but_unequal = pal[[4L]])
  cat_fill[[paste0("group_", tolower(labels[[1L]]), "_only")]] <- pal[[2L]]
  cat_fill[[paste0("group_", tolower(labels[[2L]]), "_only")]] <- pal[[3L]]
  fill <- unname(cat_fill[cats[active]])
  if (anyNA(fill)) stop("a cluster category has no colour", call. = FALSE)

  ys <- rev(seq_along(active))   # the heaviest cluster at the top
  pi_of <- function(g) hdp[["pi_mean"]][match(g, labels), active]
  left_df <- data.frame(cluster = active, group = left, y = ys, x = -pi_of(left), fill = fill)
  right_df <- data.frame(cluster = active, group = right, y = ys, x = pi_of(right), fill = fill)
  df <- rbind(left_df, right_df)
  xmax <- max(abs(df[["x"]])) * 1.15
  ticks <- seq(0, xmax * 0.9, length.out = 5L)
  breaks <- c(-rev(ticks[-1L]), 0, ticks[-1L])
  y_labels <- sprintf("C%d  β=%s", active, format(round_digits(hdp[["beta_mean"]][active], 3)))
  p <- ggplot2::ggplot(df) +
    ggplot2::geom_rect(ggplot2::aes(xmin = pmin(0, .data$x), xmax = pmax(0, .data$x),
                                    ymin = .data$y - 0.4, ymax = .data$y + 0.4,
                                    fill = .data$fill),
                       alpha = 0.85, colour = "black", linewidth = 0.2) +
    ggplot2::geom_vline(xintercept = 0, colour = "black") +
    ggplot2::annotate("text", x = -xmax * 0.95, y = length(active) + 0.6, label = left,
                      hjust = 0, fontface = "bold") +
    ggplot2::annotate("text", x = xmax * 0.95, y = length(active) + 0.6, label = right,
                      hjust = 1, fontface = "bold") +
    ggplot2::scale_fill_identity() +
    ggplot2::scale_x_continuous(limits = c(-xmax, xmax), breaks = breaks,
                                labels = format(round_digits(abs(breaks), 3))) +
    ggplot2::scale_y_continuous(breaks = ys, labels = y_labels) +
    ggplot2::labs(title = sprintf("HDP Cluster Mixing Weights: %s ← | → %s", left, right),
                  x = "Mixing weight π", y = NULL) +
    ggplot2::theme_minimal(base_size = 9)
  list(plot = p, bars = .drawn_bars(p, df))
}

# What the butterfly draws, read back from the built plot: each bar's side
# from its drawn extent, its fill as drawn. `df` supplies cluster and group
# for the rows, in the order the rect layer received them.
.drawn_bars <- function(p, df) {
  built <- ggplot2::ggplot_build(p)
  is_rect <- vapply(p[["layers"]], function(l) inherits(l[["geom"]], "GeomRect"), logical(1))
  if (sum(is_rect) != 1L) stop("expected exactly one bar layer", call. = FALSE)
  bars <- built[["data"]][[which(is_rect)]]
  if (nrow(bars) != nrow(df)) stop("the bar layer did not draw one bar per row", call. = FALSE)
  side <- ifelse(bars[["xmax"]] <= 0 & bars[["xmin"]] < 0, "left",
                 ifelse(bars[["xmin"]] >= 0 & bars[["xmax"]] > 0, "right", "none"))
  data.frame(cluster = df[["cluster"]], group = df[["group"]], side = side,
             pi = bars[["xmax"]] - bars[["xmin"]], fill = bars[["fill"]])
}
