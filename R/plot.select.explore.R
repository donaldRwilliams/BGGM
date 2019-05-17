#' Plot Exploratory (hypothesis testing) Based Graphs
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
plot.select.explore <- function(x, type,
                                lower_tri = TRUE,
                                layout = "circle",
                                node_outer = 10,
                                node_inner = 9,
                                node_text_size = 8){

  # heat map
  if(type == "heatmap"){

    # ignore network arguments
    message("type = network arguments are ignored")

    # two sided hypothesis testing
    if(x$alternative == "two.sided"){

      # plot lower triangle
      if(isTRUE(lower_tri)){

        diag(x$partials_non_zero) <- 1

        partials <- na.omit(reshape::melt(BGGM:::get_lower_tri(x$partials_non_zero)))

        partials$value <- ifelse(partials$X1 == partials$X2, 0, partials$value)

        # plot lower and upper
      } else{

        partials <- reshape::melt(x$partials_non_zero)

      }

      partials$X1 <- as.factor(partials$X1)

      # max partials for legend
      limits <- max(abs(partials$value))

      # plot 1 (non-zero)
      plt1 <- ggplot(partials, aes(x = as.factor(X2),
                                   y = as.factor(X1),
                                   fill = value)) +
        geom_tile() +
        scale_fill_gradientn(colours = c("brown3", "white", "palegreen3"),
                             name ="",
                             limits = c(-limits, limits)) +
        theme_bw(base_size = 12) +
        theme(panel.grid = element_blank()) +
        scale_y_discrete(limits= rev(levels(partials$X1)),
                         expand = c(0, 0)) +
        scale_x_discrete(expand = c(0,0))  +
        ylab("") +
        xlab("")

      # plot lower triangle
      if(isTRUE(lower_tri)){

        zeros  <- na.omit(reshape::melt(BGGM:::get_lower_tri(x$Adj_01)))

        # plot lower and upper
      } else{

        zeros <- reshape::melt(x$Adj_01)

      }

      # plot 2 (zero)
      plt2 <- ggplot(zeros, aes(x = as.factor(X2),
                                y = as.factor(X1),
                                fill = as.factor(value))) +
        geom_tile(show.legend = F) +
        scale_fill_manual(values = c("white", "black"),
                          name = "") +
        theme_bw(base_size = 12) +
        theme(panel.grid = element_blank()) +
        scale_y_discrete(limits= rev(levels(partials$X1)),
                         expand = c(0, 0)) +
        scale_x_discrete(expand = c(0,0))     +
        ylab("") +
        xlab("")

      # returned plots
      plt <- list(plot_nonzero = plt1,
                  plot_zero = plt2)


    }

    # greater than 0
    if(x$alternative == "greater"){

      # lower triangle
      if(isTRUE(lower_tri)){

        diag(x$partials_positive) <- 1

        partials <- na.omit(reshape::melt(BGGM:::get_lower_tri(x$partials_positive)))

        partials$value <- ifelse(partials$X1 == partials$X2, 0, partials$value)

        # lower and upper
      } else{

        partials <- reshape::melt(x$partials_positive)

      }

      partials$X1 <- as.factor(partials$X1)

      # max for legend
      limits <- max(abs(partials$value))

      # plot 1 (rho > 0)
      plt1 <- ggplot(partials, aes(x = as.factor(X2),
                                   y = as.factor(X1),
                                   fill = value)) +
        geom_tile() +
        scale_fill_gradientn(colours = c("brown3", "white", "palegreen3"),
                             name ="",
                             limits = c(-limits, limits)) +
        theme_bw(base_size = 12) +
        theme(panel.grid = element_blank()) +
        scale_y_discrete(limits= rev(levels(partials$X1)),
                         expand = c(0, 0)) +
        scale_x_discrete(expand = c(0,0))     +
        ylab("") +
        xlab("")

      # lower triangle
      if(isTRUE(lower_tri)){

        zeros  <- na.omit(reshape::melt(BGGM:::get_lower_tri(x$Adj_01)))
        # lower and upper
      } else{

        zeros <- reshape::melt(x$Adj_01)

      }
      # plot 2 (rho = 0)
      plt2 <- ggplot(zeros, aes(x = as.factor(X2),
                                y = as.factor(X1),
                                fill = as.factor(value))) +
        geom_tile(show.legend = F) +
        scale_fill_manual(values = c("white", "black"),
                          name = "") +
        theme_bw(base_size = 12) +
        theme(panel.grid = element_blank()) +
        scale_y_discrete(limits= rev(levels(partials$X1)),
                         expand = c(0, 0)) +
        scale_x_discrete(expand = c(0,0))     +
        ylab("") +
        xlab("")

      # returned plots
      plt <- list(plot_nonzero = plt1, plot_zero = plt2)

    }

    # less than zero
    if(x$alternative == "less"){

      # lower tri
      if(isTRUE(lower_tri)){

        diag(x$partials_negative) <- 1

        partials <- na.omit(reshape::melt(BGGM:::get_lower_tri(x$partials_negative)))

        partials$value <- ifelse(partials$X1 == partials$X2, 0, partials$value)

        # lower and upper tri
      } else{

        partials <- reshape::melt(x$partials_negative)

      }

      partials$X1 <- as.factor(partials$X1)

      # max for legend
      limits <- max(abs(partials$value))

      # plot 1 (rho < 0)
      plt1 <- ggplot(partials, aes(x = as.factor(X2),
                                   y = as.factor(X1),
                                   fill = value)) +
        geom_tile() +
        scale_fill_gradientn(colours = c("brown3", "white", "palegreen3"),
                             name ="",
                             limits = c(-limits, limits)) +
        theme_bw(base_size = 12) +
        theme(panel.grid = element_blank()) +
        scale_y_discrete(limits= rev(levels(partials$X1)),
                         expand = c(0, 0)) +
        scale_x_discrete(expand = c(0,0))     +
        ylab("") +
        xlab("")

      # lower tri
      if(isTRUE(lower_tri)){

        zeros  <- na.omit(reshape::melt(BGGM:::get_lower_tri(x$Adj_01)))

        # lower and upper tri
      } else{

        zeros <- reshape::melt(x$Adj_20)

      }

      # plot 2 (rho = 0)
      plt2 <- ggplot(zeros, aes(x = as.factor(X2),
                                y = as.factor(X1),
                                fill = as.factor(value))) +
        geom_tile(show.legend = F) +
        scale_fill_manual(values = c("white", "black"),
                          name = "") +
        theme_bw(base_size = 12) +
        theme(panel.grid = element_blank()) +
        scale_y_discrete(limits= rev(levels(partials$X1)),
                         expand = c(0, 0)) +
        scale_x_discrete(expand = c(0,0))     +
        ylab("") +
        xlab("")

      # returned plots
      plt <- list(plot_nonzero = plt1, plot_zero = plt2)

    }

    # exhaustive hyp testing
    if(x$alternative == "exhaustive"){

      # lower tri
      if(isTRUE(lower_tri)){

        pos_partial <- x$pos_mat * x$pcor_mat

        diag(pos_partial) <- 1

        partials <- na.omit(reshape::melt(BGGM:::get_lower_tri(pos_partial)))

        partials$value <- ifelse(partials$X1 == partials$X2, 0, partials$value)

        # lower and upper tri
      } else{

        pos_partial <- x$pos_mat * x$pcor_mat

        partials <- reshape::melt(pos_partial)

      }


      partials$X1 <- as.factor(partials$X1)

      # max for legend
      limits <- max(abs(partials$value))

      # plot 1 (rho > 0)
      plt1 <- ggplot(partials, aes(x = as.factor(X2),
                                   y = as.factor(X1),
                                   fill = value)) +
        geom_tile() +
        scale_fill_gradientn(colours = c("brown3", "white", "palegreen3"),
                             name ="",
                             limits = c(-limits, limits)) +
        theme_bw(base_size = 12) +
        theme(panel.grid = element_blank()) +
        scale_y_discrete(limits= rev(levels(partials$X1)),
                         expand = c(0, 0)) +
        scale_x_discrete(expand = c(0,0))     +
        ylab("") +
        xlab("")

      # lower tri
      if(isTRUE(lower_tri)){

        neg_partial <- x$neg_mat * x$pcor_mat

        diag(neg_partial) <- 1

        partials <- na.omit(reshape::melt(BGGM:::get_lower_tri(neg_partial)))

        partials$value <- ifelse(partials$X1 == partials$X2, 0, partials$value)

      } else{

        neg_partial <- x$neg_mat * x$pcor_mat

        partials <- reshape::melt(neg_partial)

      }


      partials$X1 <- as.factor(partials$X1)

      # max for legend
      limits <- max(abs(partials$value))

      # plot 2 (rho < 0)
      plt2 <- ggplot(partials, aes(x = as.factor(X2),
                                   y = as.factor(X1),
                                   fill = value)) +
        geom_tile() +
        scale_fill_gradientn(colours = c("brown3", "white", "palegreen3"),
                             name ="",
                             limits = c(-limits, limits)) +
        theme_bw(base_size = 12) +
        theme(panel.grid = element_blank()) +
        scale_y_discrete(limits= rev(levels(partials$X1)),
                         expand = c(0, 0)) +
        scale_x_discrete(expand = c(0,0))     +
        ylab("") +
        xlab("")

      # lower tri
      if(isTRUE(lower_tri)){

        adj_zero <- ifelse(x$null_mat > x$prob, 1, 0)

        zeros  <- na.omit(reshape::melt(BGGM:::get_lower_tri(adj_zero)))

        # lower and upper tri
      } else{

        adj_zero <- ifelse(x$null_mat > x$prob, 1, 0)

        zeros <- reshape::melt(adj_zero)

      }

      # plot 3 (rho = 0)
      plt3 <- ggplot(zeros, aes(x = as.factor(X2),
                                y = as.factor(X1),
                                fill = as.factor(value))) +
        geom_tile(show.legend = F) +
        scale_fill_manual(values = c("white", "black"),
                          name = "") +
        theme_bw(base_size = 12) +
        theme(panel.grid = element_blank()) +
        scale_y_discrete(limits= rev(levels(partials$X1)),
                         expand = c(0, 0)) +
        scale_x_discrete(expand = c(0,0))     +
        ylab("") +
        xlab("")

      # plots
      plt <- list(plot_H1 = plt1, plot_H2 = plt2, plot_H0 = plt3)

    }
  }

  if(type == "network")  {

    if(x$alternative == "two.sided"){

      if(sum(x$partials_positive[upper.tri(x$partials_positive)]) == 0 ){

        plt1 <- NA

      } else {


        plt1 <- BGGM:::net_plot(x$partials_non_zero,
                                layout = layout,
                                mat_type = "partials",
                                node_outer = node_outer,
                                node_inner = node_inner,
                                node_text_size = node_text_size)
      }

      if(sum(x$Adj_01[upper.tri(x$Adj_01)]) == 0){

        plt2 <- NA

      }
      plt2 <- BGGM:::net_plot(x$Adj_01,
                              layout = layout,
                              mat_type = "adj",
                              node_outer = node_outer,
                              node_inner = node_inner,
                              node_text_size = node_text_size)



      plt <- list(plot_nonzero = plt1, plot_zero = plt2)
    }


    if(x$alternative == "greater"){

      if(sum(x$partials_positive[upper.tri(x$partials_positive)]) == 0 ){
        plt1 <- NA

      } else {
        plt1 <- BGGM:::net_plot(x$partials_positive,
                                layout = layout,
                                mat_type = "partials",
                                node_outer = node_outer,
                                node_inner = node_inner,
                                node_text_size = node_text_size)
      }

      if(sum(x$Adj_01[upper.tri(x$Adj_01)] ) == 0 ) {

        plt2 <- NA

      } else {
        plt2 <- BGGM:::net_plot(x$Adj_01,
                                layout = layout,
                                mat_type = "adj",
                                node_outer = node_outer,
                                node_inner = node_inner,
                                node_text_size = node_text_size)
      }


      plt <- list(plot_nonzero = plt1, plot_zero = plt2)




    }
    if(x$alternative == "less"){

      if( sum(x$partials_negative[upper.tri(x$partials_negative)] ) == 0 ) {
        plt1 <- NA

      } else {
        plt1 <- BGGM:::net_plot(x$partials_negative,
                                layout = layout,
                                mat_type = "partials",
                                node_outer = node_outer,
                                node_inner = node_inner,
                                node_text_size = node_text_size)
      }

      if( sum(x$Adj_01[upper.tri(x$Adj_01)]) == 0 ){
        plt2 <- NA

      } else {
        plt2 <- BGGM:::net_plot(x$Adj_01,
                                layout = layout,
                                mat_type = "adj",
                                node_outer = node_outer,
                                node_inner = node_inner,
                                node_text_size = node_text_size)

      }

      plt <- list(plot_nonzero = plt1, plot_zero = plt2)

    }


    if(x$alternative == "exhaustive"){

      if(sum(x$pos_mat[upper.tri(x$pos_mat)]) == 0  ){
        plt1 <- NA

      } else{

        plt1 <- BGGM:::net_plot(x$pos_mat * x$pcor_mat,
                                layout = layout,
                                mat_type = "partials",
                                node_outer = node_outer,
                                node_inner = node_inner,
                                node_text_size = node_text_size) +
          ggtitle( bquote("   p("~H[1]~"|"~bold(Y)~")" >= .(x$prob)))
      }

      if(sum(x$neg_mat[upper.tri(x$neg_mat)]) == 0){
        plt2 <- NA

      } else{

        plt2 <- BGGM:::net_plot(x$neg_mat * x$pcor_mat,
                                layout = layout,
                                mat_type = "partials",
                                node_outer = node_outer,
                                node_inner = node_inner,
                                node_text_size = node_text_size) +
          ggtitle(bquote("   p("~H[2]~"|"~bold(Y)~")" >= .(x$prob)))
      }

      adj_zero <- ifelse(x$null_mat > x$prob, 1, 0)

      if(sum(adj_zero[upper.tri(adj_zero)])   == 0) {
        plt3 <- NA

      } else{
        plt3 <-  BGGM:::net_plot(adj_zero,
                                 layout = layout,
                                 mat_type = "adj",
                                 node_outer = node_outer,
                                 node_inner = node_inner,
                                 node_text_size = node_text_size) +
          ggtitle( bquote("   p("~H[0]~"|"~bold(Y)~")" >= .(x$prob)))
      }
      plt <- list(plot_H1 = plt1, plot_H2 = plt2, plot_H0 = plt3)
    }
  }
  plt
}

