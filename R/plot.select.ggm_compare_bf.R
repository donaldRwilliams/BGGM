#' Plot Edge Differences and Equivalence Between GGMs (\code{ggm_compare_bf })
#'
#' @param x object of class \code{select.ggm_compare_bf}
#' @param mat_type "adj" or "BF"
#' @param type \code{heatmap} or \code{network} plot (network is currently the only option)
#' @param lower_tri heatmap lower triangular elements
#' @param layout layout of network plot. see igraph for options
#' @param node_outer border for node
#' @param node_inner node size (see notes)
#' @param node_text_size  node text size
#' @param labels node labels. Default is 1:p, but can change to colnames(Y)
#'
#' @return
#' @export
#'
#' @examples
plot.select.ggm_compare_bf <- function(x, mat_type = "adj", type = "network", lower_tri = TRUE,
                                       layout = "circle", node_outer = 10,
                                       node_inner = 9,
                                       node_text_size = 8, labels = NULL){

  if(mat_type == "BF"){


    plt1 <- BGGM:::net_plot(x$BF_10_adj,
                            layout = layout,
                            mat_type = "partials",
                            node_outer = node_outer,
                            node_inner = node_inner,
                            node_text_size = node_text_size,
                            labels = labels)

    plt2 <- BGGM:::net_plot((x$BF_01_adj),
                            layout = layout,
                            mat_type = "partials",
                            node_outer = node_outer,
                            node_inner = node_inner,
                            node_text_size = node_text_size,
                            labels = labels)


    temp_undecided <- ifelse( x$adj_10 + x$adj_01 == 0, 1, 0)

    diag(temp_undecided) <- 0

    plt3 <- BGGM:::net_plot(temp_undecided,
                            layout = layout,
                            mat_type = "adj",
                            node_outer = node_outer,
                            node_inner = node_inner,
                            node_text_size = node_text_size,
                            labels = labels)








    plt <- list(plot_nonzero = plt1, plot_zero = plt2, plt_undecided = plt3)
  }


  if(mat_type == "adj"){


    plt1 <- BGGM:::net_plot(x$adj_10,
                            layout = layout,
                            mat_type = "adj",
                            node_outer = node_outer,
                            node_inner = node_inner,
                            node_text_size = node_text_size,
                            labels = labels)

    plt2 <- BGGM:::net_plot(x$adj_01,
                            layout = layout,
                            mat_type = "adj",
                            node_outer = node_outer,
                            node_inner = node_inner,
                            node_text_size = node_text_size,
                            labels = labels)


    temp_undecided <- ifelse( x$adj_10 + x$adj_01 == 0, 1, 0)

    diag(temp_undecided) <- 0

    plt3 <- BGGM:::net_plot(temp_undecided,
                            layout = layout,
                            mat_type = "adj",
                            node_outer = node_outer,
                            node_inner = node_inner,
                            node_text_size = node_text_size,
                            labels = labels)








    plt <- list(plot_nonzero = plt1, plot_zero = plt2, plt_undecided = plt3)
  }
  plt


}
