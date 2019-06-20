#' Plot Edge Differences and Equivalence Between GGMs (\code{ggm_compare_estimate})
#'
#' @param x object of class \code{select.ggm_compare_estimate}
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
#' Y1 <- MASS::mvrnorm(5000, rep(0, 16), Sigma = diag(16))
#' Y2 <- MASS::mvrnorm(5000, rep(0, 16), Sigma = diag(16))
#' Y3 <- MASS::mvrnorm(5000, rep(0, 16), Sigma = diag(16))
#"
#' ggm_est <- ggm_compare_estimate(Y1, Y2, Y3,  iter = 5000, ci_width = .95)
#'
#'sel <- select(ggm_est, type = "ci")
#'
#'plot(x, labels = letters[1:16])


plot.select.ggm_compare_estimate <- function(x, type = "network", lower_tri = TRUE,
                                             layout = "circle", node_outer = 10,
                                             node_inner = 9,
                                             node_text_size = 8, labels = NULL){





  if(!is.null(x$rope)){
    plt1 <- list()
    plt2 <- list()
    plt3 <- list()


    for(i in 1:length(x$mats_null)){




      plt1[[i]] <- BGGM:::net_plot(x$mats_diff[[i]][[1]],
                                   layout = layout,
                                   mat_type = "adj",
                                   node_outer = node_outer,
                                   node_inner = node_inner,
                                   node_text_size = node_text_size,
                                   labels = labels)

      plt2[[i]] <- BGGM:::net_plot(x$mats_null[[i]][[1]],
                                   layout = layout,
                                   mat_type = "adj",
                                   node_outer = node_outer,
                                   node_inner = node_inner,
                                   node_text_size = node_text_size,
                                   labels = labels)


      temp_undecided <- ifelse( x$mats_diff[[i]][[1]] + x$mats_null[[i]][[1]] == 0, 1, 0)
      diag(temp_undecided) <- 0

      plt3[[i]] <- BGGM:::net_plot(temp_undecided,
                                   layout = layout,
                                   mat_type = "adj",
                                   node_outer = node_outer,
                                   node_inner = node_inner,
                                   node_text_size = node_text_size,
                                   labels = labels)




    }

    temp_name <- unlist(lapply(1:length(x$mats_diff), function(z) names( x$mats_diff[[z]])))
    names(plt1) <- temp_name
    names(plt2) <- temp_name
    names(plt3) <- temp_name

    plt <- list(plot_nonzero = plt1, plot_zero = plt2, plt_undecided = plt3)

  }

  if(x$type == "ci"){
    plt1 <- list()

    for(i in 1:length(x$mats_diff)){




      plt1[[i]] <- BGGM:::net_plot(x$mats_diff[[i]][[1]],
                                   layout = layout,
                                   mat_type = "adj",
                                   node_outer = node_outer,
                                   node_inner = node_inner,
                                   node_text_size = node_text_size,
                                   labels = labels)
    }
    temp_name <- unlist(lapply(1:length(x$mats_diff), function(z) names( x$mats_diff[[z]])))
    names(plt1) <- temp_name
    plt <- list(plot_nonzero = plt1)
  }

  plt

}
