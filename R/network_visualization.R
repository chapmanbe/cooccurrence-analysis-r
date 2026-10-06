# Network plots (ggplot2). Each returns a ggplot object; saving is the
# caller's job. Figures are judged by eye -- the data behind each one is
# what gets compared -- so these aim for legibility, not pixel parity with
# the Julia (Makie) versions.

# Makie's `wong_colors()`, exactly as the Julia reference uses it: seven
# color-blind-safe hues (Okabe-Ito without black), recycled by community
# label -- community c gets hue ((c - 1) mod 7) + 1.
.community_palette <- c("#0072B2", "#E69F00", "#009E73", "#CC79A7",
                        "#56B4E9", "#D55E00", "#F0E442")

.minmax_scale <- function(x, lo, hi) {
  if (length(x) == 0L) return(numeric(0))
  r <- range(x)
  if (r[[2L]] == r[[1L]]) return(rep((lo + hi) / 2, length(x)))
  lo + (x - r[[1L]]) / (r[[2L]] - r[[1L]]) * (hi - lo)
}

#' Plot a co-occurrence network
#'
#' Node area grows with prevalence, node color is community membership (when
#' `communities` is given), edge width grows with weight. The layout is
#' Kamada-Kawai on the unweighted topology, which igraph computes from a
#' deterministic starting layout, so the figure is reproducible.
#' @export
plot_cooccurrence_network <- function(net, communities = NULL,
                                      title = "Item Co-occurrence Network",
                                      node_size_range = c(3, 12),
                                      edge_width_range = c(0.4, 3)) {
  stopifnot(inherits(net, "cooccurrence_network"))
  items <- net[["items"]]
  if (length(items) == 0L) stop("cannot plot an empty network", call. = FALSE)
  topo <- net[["adjacency"]] != 0
  g <- igraph::graph_from_adjacency_matrix(topo, mode = "undirected", diag = FALSE)
  xy <- igraph::layout_with_kk(g)

  comm <- if (is.null(communities)) rep(1L, length(items)) else
    as.integer(unname(communities[["assignments"]][items]))
  nodes <- data.frame(
    item = items, x = xy[, 1L], y = xy[, 2L],
    size = .minmax_scale(as.numeric(net[["prevalence"]][items]),
                         node_size_range[[1L]], node_size_range[[2L]]),
    community = factor(comm, levels = sort(unique(comm)))
  )
  e <- network_edges(net)
  ia <- match(e[["item_a"]], items)
  ib <- match(e[["item_b"]], items)
  edges <- data.frame(x = xy[ia, 1L], y = xy[ia, 2L], xend = xy[ib, 1L], yend = xy[ib, 2L],
                      width = .minmax_scale(e[["weight"]], edge_width_range[[1L]],
                                            edge_width_range[[2L]]))
  palette <- rep_len(.community_palette, nlevels(nodes[["community"]]))
  names(palette) <- levels(nodes[["community"]])

  ggplot2::ggplot() +
    ggplot2::geom_segment(data = edges,
                          ggplot2::aes(x = .data$x, y = .data$y,
                                       xend = .data$xend, yend = .data$yend,
                                       linewidth = .data$width),
                          color = "gray55", alpha = 0.7) +
    ggplot2::geom_point(data = nodes,
                        ggplot2::aes(x = .data$x, y = .data$y, size = .data$size,
                                     color = .data$community)) +
    ggplot2::geom_text(data = nodes, ggplot2::aes(x = .data$x, y = .data$y, label = .data$item),
                       size = 3, vjust = -1.4) +
    ggplot2::scale_linewidth_identity() +
    ggplot2::scale_size_identity() +
    ggplot2::scale_color_manual(values = palette, name = "Community") +
    ggplot2::guides(color = ggplot2::guide_legend(override.aes = list(size = 4))) +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = 0.15)) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = 0.12)) +
    ggplot2::labs(title = title) +
    ggplot2::theme_void() +
    ggplot2::theme(plot.title = ggplot2::element_text(size = 16, hjust = 0.5))
}

#' Heatmap of edge weights, rows and columns ordered by community
#'
#' Items are ordered by community, then by name; dashed lines mark the
#' community boundaries.
#' @export
plot_community_heatmap <- function(net, communities, weight_label = net[["weight_metric"]],
                                   title = "Co-occurrence Strength by Community",
                                   axis_label = "Item") {
  stopifnot(inherits(net, "cooccurrence_network"))
  items <- net[["items"]]
  n <- length(items)
  if (n == 0L) stop("cannot plot a heatmap for an empty network", call. = FALSE)
  comm <- as.integer(unname(communities[["assignments"]][items]))
  o <- order(comm, items, method = "radix")
  ordered <- items[o]
  mat <- net[["adjacency"]][ordered, ordered, drop = FALSE]
  cells <- data.frame(
    x = factor(rep(ordered, each = n), levels = ordered),
    y = factor(rep(ordered, times = n), levels = ordered),
    weight = as.vector(t(mat))
  )
  bounds <- which(diff(comm[o]) != 0) + 0.5
  ggplot2::ggplot(cells, ggplot2::aes(x = .data$x, y = .data$y, fill = .data$weight)) +
    ggplot2::geom_tile() +
    ggplot2::geom_hline(yintercept = bounds, linetype = "dashed", linewidth = 0.5) +
    ggplot2::geom_vline(xintercept = bounds, linetype = "dashed", linewidth = 0.5) +
    ggplot2::scale_fill_gradient(low = "#FFFFCC", high = "#BD0026", name = weight_label) +
    ggplot2::labs(title = title, x = axis_label, y = axis_label) +
    ggplot2::coord_equal() +
    ggplot2::theme_minimal() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 60, hjust = 1, size = 7),
                   axis.text.y = ggplot2::element_text(size = 7),
                   panel.grid = ggplot2::element_blank())
}

#' Horizontal bar chart of the top nodes by a centrality metric
#'
#' Input is [compute_network_metrics()] output; the bars are
#' [top_by_centrality()], largest at the top.
#' @export
plot_centrality_barchart <- function(metrics, centrality = "strength", top_n = 20L) {
  if (nrow(metrics) == 0L) stop("cannot plot an empty metrics table", call. = FALSE)
  top <- top_by_centrality(metrics, centrality, top_n)
  top[["item"]] <- factor(top[["item"]], levels = rev(top[["item"]]))
  ggplot2::ggplot(top, ggplot2::aes(x = .data[[centrality]], y = .data$item)) +
    ggplot2::geom_col(fill = "steelblue") +
    ggplot2::labs(title = paste("Network Centrality:", centrality), x = centrality, y = NULL) +
    ggplot2::theme_minimal()
}

#' Two group networks side by side
#'
#' `strata` is a list of exactly two elements named by group, each holding
#' `network` and `communities` -- the shape [stratified_network_analysis()]
#' returns. Panels follow the list's name order, so a caller passing that
#' result gets its sorted order, and one wanting another order reorders by
#' name (`strata[c("B", "A")]`), never by position. Drawn as the Julia
#' reference draws it: each panel laid out on its own (here Kamada-Kawai on
#' the unweighted topology), node area by prevalence, node colour by the
#' group's own communities, edges gray (`gray70`, alpha 0.6) with width by
#' weight rescaled within the panel. No element carries a per-group colour.
#' An empty network gets an empty panel.
#'
#' Returns a list: `plot` (the ggplot) and `panels`, one row per panel read
#' back from `ggplot2::ggplot_build()` of that plot -- `group`, `panel`
#' (1 = left) and `edge_colour` (hex, `NA` for a panel without edges) -- so
#' a caller records what the figure draws, not what was asked for.
#' @export
plot_network_comparison <- function(strata, item_label = "items") {
  if (!is.list(strata) || is.null(names(strata)) || any(!nzchar(names(strata)))) {
    stop("strata must be a list named by group", call. = FALSE)
  }
  if (length(strata) != 2L) {
    stop(sprintf("plot_network_comparison compares two groups; got %d", length(strata)),
         call. = FALSE)
  }
  groups <- names(strata)
  node_parts <- list()
  edge_parts <- list()
  titles <- character(0)
  for (grp in groups) {
    net <- strata[[grp]][["network"]]
    stopifnot(inherits(net, "cooccurrence_network"))
    items <- net[["items"]]
    size <- network_size(net)
    titles[[grp]] <- sprintf("%s Network (%d %s, %d edges)", grp, size[["n_nodes"]],
                             item_label, size[["n_edges"]])
    if (length(items) == 0L) next
    g <- igraph::graph_from_adjacency_matrix(net[["adjacency"]] != 0, mode = "undirected",
                                             diag = FALSE)
    xy <- igraph::layout_with_kk(g)
    comm <- as.integer(unname(strata[[grp]][["communities"]][["assignments"]][items]))
    if (anyNA(comm)) stop(sprintf("a vertex of '%s' has no community", grp), call. = FALSE)
    size_scaled <- .minmax_scale(as.numeric(net[["prevalence"]][items]), 3, 10)
    node_parts[[grp]] <- data.frame(group = grp, item = items, x = xy[, 1L], y = xy[, 2L],
                                    size = size_scaled, community = comm)
    e <- network_edges(net)
    ia <- match(e[["item_a"]], items)
    ib <- match(e[["item_b"]], items)
    edge_parts[[grp]] <- data.frame(group = rep(grp, nrow(e)), x = xy[ia, 1L], y = xy[ia, 2L],
                                    xend = xy[ib, 1L], yend = xy[ib, 2L],
                                    width = .minmax_scale(e[["weight"]], 0.4, 3))
  }
  empty_nodes <- data.frame(group = character(0), item = character(0), x = numeric(0),
                            y = numeric(0), size = numeric(0), community = integer(0))
  empty_edges <- data.frame(group = character(0), x = numeric(0), y = numeric(0),
                            xend = numeric(0), yend = numeric(0), width = numeric(0))
  nodes <- do.call(rbind, c(list(empty_nodes), unname(node_parts)))
  edges <- do.call(rbind, c(list(empty_edges), unname(edge_parts)))
  nodes[["group"]] <- factor(nodes[["group"]], levels = groups)
  edges[["group"]] <- factor(edges[["group"]], levels = groups)
  levels_c <- sort(unique(nodes[["community"]]))
  nodes[["community"]] <- factor(nodes[["community"]], levels = levels_c)
  hue <- (levels_c - 1L) %% length(.community_palette) + 1L
  palette <- stats::setNames(.community_palette[hue], levels_c)

  p <- ggplot2::ggplot() +
    ggplot2::geom_segment(data = edges,
                          ggplot2::aes(x = .data$x, y = .data$y,
                                       xend = .data$xend, yend = .data$yend,
                                       linewidth = .data$width),
                          color = "gray70", alpha = 0.6) +
    ggplot2::geom_point(data = nodes,
                        ggplot2::aes(x = .data$x, y = .data$y, size = .data$size,
                                     color = .data$community)) +
    ggplot2::geom_text(data = nodes, ggplot2::aes(x = .data$x, y = .data$y, label = .data$item),
                       size = 2.4, vjust = -1.3) +
    ggplot2::facet_wrap(ggplot2::vars(.data$group), nrow = 1L, scales = "free", drop = FALSE,
                        labeller = ggplot2::as_labeller(titles)) +
    ggplot2::scale_linewidth_identity() +
    ggplot2::scale_size_identity() +
    ggplot2::scale_color_manual(values = palette, name = "Community", drop = FALSE) +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = 0.15)) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = 0.12)) +
    ggplot2::theme_void() +
    ggplot2::theme(strip.text = ggplot2::element_text(size = 13,
                                                      margin = ggplot2::margin(b = 6)))
  list(plot = p, panels = .drawn_panels(p))
}

#' Panel order and edge colour as the built plot draws them
#'
#' Reads `ggplot2::ggplot_build(p)`: the facet layout gives each group's
#' column (the plot is one row), and the segment layer's built data give the
#' colour each panel's edges are drawn in, as hex. Raises if a panel's edges
#' are drawn in more than one colour, or if the plot has no single facet
#' column named `group` or no single segment layer.
#' @keywords internal
.drawn_panels <- function(p) {
  built <- ggplot2::ggplot_build(p)
  layout <- built[["layout"]][["layout"]]
  if (!("group" %in% names(layout)) || any(layout[["ROW"]] != 1L)) {
    stop("expected a one-row plot faceted by 'group'", call. = FALSE)
  }
  is_segment <- vapply(p[["layers"]], function(l) inherits(l[["geom"]], "GeomSegment"),
                       logical(1))
  if (sum(is_segment) != 1L) stop("expected exactly one segment layer", call. = FALSE)
  seg <- built[["data"]][[which(is_segment)]]
  edge_colour <- vapply(seq_len(nrow(layout)), function(i) {
    in_panel <- as.integer(seg[["PANEL"]]) == as.integer(layout[["PANEL"]][[i]])
    cols <- unique(seg[["colour"]][in_panel])
    if (length(cols) == 0L) return(NA_character_)
    if (length(cols) > 1L) {
      stop(sprintf("panel '%s' draws edges in more than one colour", layout[["group"]][[i]]),
           call. = FALSE)
    }
    rgb <- grDevices::col2rgb(cols)
    grDevices::rgb(rgb[1L, 1L], rgb[2L, 1L], rgb[3L, 1L], maxColorValue = 255)
  }, character(1))
  out <- data.frame(group = as.character(layout[["group"]]),
                    panel = as.integer(layout[["COL"]]),
                    edge_colour = edge_colour)
  out <- out[order(out[["panel"]]), , drop = FALSE]
  rownames(out) <- NULL
  out
}
