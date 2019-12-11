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
#' @return objects of class \code{ggplot}
#' @export
#'
#' @examples
#' # data (null is true)
#' Y1 <- MASS::mvrnorm(500, mu = rep(0,16),
#'                     Sigma = Sigma = BGGM::ptsd_cor2,
#'                     empirical = T)
#' Y2 <- MASS::mvrnorm(500, mu = rep(0,16),
#'                     Sigma = BGGM::ptsd_cor3,
#'                     empirical = T)
#' # fit model
#' fit <- ggm_compare_bf(Y1, Y2,
#'                      prior_sd = 0.22,
#'                      iter = 5000,
#'                      cores = 4)
#'# select E
#'sel <- select(fit, BF_cut = 3)
#'# plot
#'plot(sel, mat_type = "BF", type = "heatmap")

plot.select.ggm_compare_bf <- function(x, mat_type = "adj",
                                       type = "heatmap",
                                       lower_tri = TRUE,
                                       layout = "circle",
                                       node_outer = 10,
                                       node_inner = 9,
                                       node_text_size = 8, labels = NULL){

  if(type == "network"){

    if(mat_type == "BF"){
      # log BF
      mat_temp1 <- log(x$BF_10_adj)
      mat_temp1 <- ifelse(is.infinite(mat_temp1), 0, mat_temp1)

      plt1 <- BGGM:::net_plot(mat_temp1,
                            layout = layout,
                            mat_type = "partials",
                            node_outer = node_outer,
                            node_inner = node_inner,
                            node_text_size = node_text_size,
                            labels = labels)

      # log BF
      mat_temp2 <- log(x$BF_01_adj)
      mat_temp2 <- ifelse(is.infinite(mat_temp2), 0, mat_temp2)

    plt2 <- BGGM:::net_plot(mat_temp2,
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

    plt <- list(plot_nonzero = plt1,
                plot_zero = plt2,
                plt_undecided = plt3)
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

    plt <- list(plot_nonzero = plt1,
                plot_zero = plt2,
                plt_undecided = plt3)
  }
}
  if(type == "heatmap"){

    cutoff <- x$BF

    melt_01 <- reshape::melt(x$BF_01)
    melt_10 <- reshape::melt(x$BF_10)

    melt_01 <- subset(melt_01, value > cutoff)
    melt_01$X1 <- as.factor(melt_01$X1)

    max_01 <- log(max(melt_01$value))
    melt_10 <- subset(melt_10, value > cutoff)

    melt_10$X1 <- as.factor(melt_10$X1)
    max_10 <- log(max(melt_10$value))

      plt_null <- ggplot(melt_01, aes(x = as.factor(X2),
                                      y = as.factor(X1),
                                      fill = log(value))) +
        geom_tile() +
        scale_fill_gradientn(colours = c("white", "yellow", "red"),
                             values = c(0, .1, 1),
                             limits = c(log(cutoff), max_01),
                             labels = round(c(log(cutoff),    max_01)),
                             breaks =  c(log(cutoff),  max_01),
                             name = "BF 01") +

        xlab("Null Hypothesis Matrix") + ylab(" ") +
        ylim(rev(levels(melt_01$X1)))+
        theme_bw(base_size = 12) +
        theme(panel.grid = element_blank())  +
        theme(panel.grid = element_blank(),
              legend.title=element_text(size=9),
              legend.text=element_text(size=9))

      plt_alt <- ggplot(melt_10, aes(x = as.factor(X2),
                                     y = as.factor(X1),
                                     fill = round(log(value)))) +
        geom_tile() +
        scale_fill_gradientn(colours = c("white", "lightblue", "purple"),
                             values = c(0, .1, 1), limits = round(c(log(cutoff), max_10)),
                             labels = round(c(log(cutoff),   max_10)),
                             breaks =  round(c(log(cutoff),  max_10)),
                             name = "BF 10") +
        xlab("Alternative Hypothesis Matrix") + ylab(" ") +
        ylim(rev(levels(melt_10$X1)))+
        theme_bw(base_size = 12) +
        theme(panel.grid = element_blank(),
              legend.title=element_text(size=9),
              legend.text=element_text(size=9))

      plt <- list(plt_null = plt_null, plt_alt = plt_alt)
      }
  plt
}
