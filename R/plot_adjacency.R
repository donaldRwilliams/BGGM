#' Plot Adjacency Matrix
#'
#' @param x adjacency matrix
#' @param layout network layout (\link[sna]{gplot.layout})
#' @param node_labels node labels
#' @param node_labels_color node labels color
#' @param node_groups node group indicator
#' @param node_outer_size node border size
#' @param node_inner_size node size
#' @param alpha edge transparency
#' @param txt_size node text size
#' @param ... additional arguments (\link[GGally]{ggnet2})
#' @importFrom GGally ggnet2
#' @importFrom ggplot2 ggtitle
#' @importFrom network network.vertex.names<- set.edge.value set.edge.attribute %e% %v%<-
#' @importFrom sna gplot.layout.circle
#' @return object of class \code{ggplot}
#' @export
#'
#' @examples
#'
#' # data
#'Y <- rsa[,-ncol(rsa)]
#'
#'# fit model
#'fit_explore <- explore(Y, iter = 1000)
#'
#'# select the graph (edge set E)
#'E <- select(fit_explore)
#'
#'plot_adjacency(E$Adj_10, node_groups = BGGM:::rsa_labels)
plot_adjacency <-  function(x,  layout = "circle",
                            node_labels = NULL,
                            node_labels_color = "black",
                            node_groups = NULL,
                            node_outer_size = 12,
                            node_inner_size = 11,
                            alpha = 0.50, txt_size = 8,
                            ... ){

p <- ncol(x)
net <- network::network(x, directed = FALSE)

# default labels
if (is.null(node_labels)) {
  network::network.vertex.names(net) <- 1:p
  # custom labels
} else {
  # check label length
  if (isFALSE(length(node_labels) == p)) {
    stop("labels must be of length p (number of nodes)")
  }

  network::network.vertex.names(net) <- node_labels
}


if (is.null(node_groups)) {
  plt_null <- ggnet2(
    net = net, edge.alpha = alpha,
    mode = layout, node.size = node_outer_size,
    node.color = "black",
    label = TRUE
  ) +
    geom_point(color = "white",
               size = node_inner_size,
               alpha = 1) +
    geom_text(aes(label = label),
              color = node_labels_color,
              size = txt_size)

  plt$plt_null <- plt_null

} else{
  if (isFALSE(length(node_groups) == p)) {
    stop("labels must be of length p (number of nodes)")
  }

  net %v% "group" <- node_groups

  plt_null <- ggnet2(
    net = net,edge.alpha = alpha,
    mode = layout, node.alpha = 0.5,
    node.size = node_outer_size,
    node.color = "group",
    label = TRUE,
    ...
  ) +
    geom_point(aes(color = color),
               size = node_inner_size,
               alpha = 1) +
    geom_text(aes(label = label),
              color = node_labels_color,
              size = txt_size)


}
plt_null
}

