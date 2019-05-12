#' Plot Estimation Based Graphs
#'
#' @param x object of class \code{select.estimate}
#' @param type \code{heatmap} or \code{network} plot
#' @param lower_tri heatmap lower triangular elements
#' @param layout layout of network plot. see igraph for options
#' @param node_outer border for node
#' @param node_inner node size (see notes)
#' @param node_text_size  node text size
#'
#' @return
#' @export
#'
#' @note
#' @examples
plot.select.estimate <- function(x, type,
                                 lower_tri = TRUE,
                                 layout = "circle",
                                 node_outer = 10,
                                 node_inner = 9,
                                 node_text_size = 8){

  if(type == "heatmap"){

    message("type = network arguments are ignored")

    if(is.null(x$rope)){

      if(isTRUE(lower_tri)){

        diag(x$partials) <- 1

        partials <- na.omit(reshape::melt(BGGM:::get_lower_tri(x$partials)))

        partials$value <- ifelse(partials$X1 == partials$X2, 0, partials$value)

        } else{

        partials <- reshape::melt(x$partials)

      }

    partials$X1 <- as.factor(partials$X1)

    limits <- max(abs(partials$value))

    plt <- ggplot(partials, aes(x = as.factor(X2),
                                y = as.factor(X1),
                                fill = value)) +
      geom_tile() +
      scale_fill_gradientn(colours = c("brown3", "white", "palegreen3"),
                           name ="", limits = c(-limits, limits)) +
      theme_bw(base_size = 12) +
      theme(panel.grid = element_blank()) +
      scale_y_discrete(limits= rev(levels(partials$X1)), expand = c(0, 0)) +
      scale_x_discrete(expand = c(0,0))     +
      ylab("") +
      xlab("")
  }


  if(!is.null(x$rope)){

    if(isTRUE(lower_tri)){
      diag(x$partials_non_zero) <- 1
      partials <- na.omit(reshape::melt(BGGM:::get_lower_tri(x$partials_non_zero)))
      partials$value <- ifelse(partials$X1 == partials$X2, 0, partials$value)

    } else{

      partials <- reshape::melt(x$partials_non_zero)

      }
    partials$X1 <- as.factor(partials$X1)

    limits <- max(abs(partials$value))

    plt1 <- ggplot(partials, aes(x = as.factor(X2),
                                 y = as.factor(X1),
                                 fill = value)) +
      geom_tile(color = "white") +
       scale_fill_gradientn(colours = c("brown3", "white", "palegreen3"), name =
       "", limits = c(-limits, limits)) +
      theme_bw(base_size = 12) +
      theme(panel.grid = element_blank(), plot.background = element_rect(fill  = "white")) +
      scale_y_discrete(limits= rev(levels(partials$X1)), expand = c(0, 0)) +
      scale_x_discrete(expand = c(0,0))     +
      ylab("") +
      xlab("")


    if(isTRUE(lower_tri)){

      zeros  <- na.omit(reshape::melt(BGGM:::get_lower_tri(x$adjacency_zero)))
    } else{

      zeros <- reshape::melt(x$adjacency_zero)

    }

    # zeros$X1 <- as.factor(zeros$X1)
    # zeros$value <- ifelse(x$in_rope > x$prob, 1, 0)

   plt2 <- ggplot(zeros, aes(x = as.factor(X2),
                              y = as.factor(X1),
                              fill = as.factor(value))) +
      geom_tile(show.legend = F) +
      scale_fill_manual(values = c("white", "black"), name = "") +
      theme_bw(base_size = 12) +
      theme(panel.grid = element_blank()) +
      scale_y_discrete(limits= rev(levels(partials$X1)), expand = c(0, 0)) +
      scale_x_discrete(expand = c(0,0))     +
      ylab("") +
      xlab("")


    plt <- list(plot_nonzero = plt1, plot_zero = plt2)

    }
}

  if(type == "network"){
    if(is.null(x$rope)){
   plt <- BGGM:::net_plot(x$partials,
                    layout = layout,
                    mat_type = "partials",
                    node_outer = node_outer,
                    node_inner = node_inner,
                    node_text_size = node_text_size)


  }

  if(!is.null(x$rope)){

    plt1 <- BGGM:::net_plot(x$partials_non_zero,
                            layout = layout,
                            mat_type = "partials",
                            node_outer = node_outer,
                            node_inner = node_inner,
                            node_text_size = node_text_size)
    plt2 <- BGGM:::net_plot(x,
                            layout = layout,
                            mat_type = "adj",
                            node_outer = node_outer,
                            node_inner = node_inner,
                            node_text_size = node_text_size)



  plt <- list(plot_nonzero = plt1, plot_zero = plt2)

  }

}

  plt





}
