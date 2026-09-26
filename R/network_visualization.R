# Network plots (ggplot2). Each returns a ggplot object; saving is the
# caller's job. Figures are judged by eye -- the data behind each one is
# what gets compared -- so these aim for legibility, not pixel parity with
# the Julia (Makie) versions.

# Okabe-Ito / Makie "wong" palette: color-blind safe, eight hues, recycled.
.community_palette <- c("#0072B2", "#E69F00", "#009E73", "#CC79A7",
                        "#56B4E9", "#D55E00", "#F0E442", "#000000")

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
