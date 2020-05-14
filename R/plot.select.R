#' Network Plot for \code{select} Objects
#'
#' Visualize the conditional (in)dependence structure.
#'
#' @param x An object of class \code{select}.
#'
#' @param layout Character string. Which graph layout (defaults is \code{circle}) ?.
#'                See \link[sna]{gplot.layout}.
#'
#' @param pos_col Character string. Color for the positive edges (defaults to \code{green}).
#'
#' @param neg_col Character string.  Color for the negative edges (defaults to \code{green}).
#'
#' @param node_size Numeric. The size of the nodes (defaults to \code{10}).
#'
#' @param edge_magnify Numeric. A value that is multiplied by the edge weights. This can increase (> 1) or
#'                     derease (< 1) the line widths.
#'
#' @param groups A character string of length \emph{p} (the number of nodes in the model).
#'               This indicates groups of nodes that should be the same color
#'               (e.g., "clusters" or "communities").
#'
#' @param palette A character string sepcifying the palette for the \code{groups}.
#'                (default is \code{Set3}). See \href{http://www.cookbook-r.com/Graphs/Colors_(ggplot2)/}{palette options here}.
#'
#' @param ... Additional options passed to \link[GGally]{ggnet2}
#'
#' @importFrom GGally ggnet2
#'
#' @importFrom network network.vertex.names<- set.edge.value set.edge.attribute %e% %v%<- network
#'
#' @importFrom sna gplot.layout.circle
#'
#' @return An object (or list of objects) of class \code{ggplot}
#' that can then be further customized.
#'
#'
#' @examples
#' \donttest{
#' #########################
#' ### example 1: one ggm ##
#' #########################
#'
#' # data
#' Y <- bfi[,1:25]
#'
#' # estimate
#' fit <- estimate(Y, iter = 250)
#'
#' # "communities"
#' comm <- substring(colnames(Y), 1, 1)
#'
#' # edge set
#' E <- select(fit)
#'
#' # plot edge set
#' plt_E <- plot(E, edge_magnify = 5,
#'               palette = "Set1",
#'               groups = comm)
#'
#'
#' #############################
#' ### example 2: ggm compare ##
#' #############################
#' # compare males vs. females
#'
#' # data
#' Y <- bfi[,1:26]
#'
#' Ym <- subset(Y, gender == 1,
#'              select = -gender)
#'
#' Yf <- subset(Y, gender == 2,
#'               select = -gender)
#'
#' # estimate
#' fit <- ggm_compare_estimate(Ym, Yf, iter = 250)
#'
#' # "communities"
#' comm <- substring(colnames(Ym), 1, 1)
#'
#' # edge set
#' E <- select(fit)
#'
#' # plot edge set
#' plt_E <- plot(E, edge_magnify = 5,
#'               palette = "Set1",
#'               groups = comm)
#'
#'
#'}
#'
#' @export

plot.select <- function(x,
                        layout = "circle",
                        pos_col = "#009E73",
                        neg_col = "#D55E00",
                        node_size = 10,
                        edge_magnify = 1,
                        groups = NULL,
                        palette = "Set3",
                        ...){

  # select estimate
  if(is(x, "select.estimate")){

    cn <- colnames(x$object$Y)

    p <- ncol(x$pcor_adj)

    diag(x$pcor_adj) <- 0

    net <- network::network(x$pcor_adj)

    if(is.null(cn) ) {
      cn <- 1:p
    }

    # edge weights
    network::set.edge.value(x = net, attrname = "weights",
                            value = x$pcor_adj)

    # edge weights absolute
    network::set.edge.value(x = net, attrname = "abs_weights",
                            value = abs(x$pcor_adj) * edge_magnify)


    # edge colors
    network::set.edge.attribute(x = net, attrname = "edge_color",
                                value = ifelse(net %e% "weights" < 0,
                                               neg_col,
                                               pos_col))

    e <- abs(as.numeric(x$pcor_adj))

    if(is.null(groups)){

      plt <- ggnet2(net, edge.alpha = e[e != 0] / max(e),
                    edge.size = "abs_weights",
                    edge.color = "edge_color",
                    node.size = 1,
                    mode = layout) +
        geom_point(color = "black",
                   size = node_size+1) +
        geom_point(size = node_size, color = "white")  +
        guides(color = FALSE) +
        geom_text(label = cn)


    } else {

      net %v% "group" <- groups

      suppressMessages(
      plt <-  ggnet2(net, edge.alpha = e[e != 0] / max(e),
                     edge.size = "abs_weights",
                     edge.color = "edge_color",
                     node.color = "group",
                     node.size = 1,
                     mode = layout) +
        geom_point(aes(color = groups),
                   size = node_size + 2,
                   alpha = 0.2) +
        geom_point(aes(color = groups),
                   size = node_size,
                   alpha = 1) +
        guides(colour = guide_legend(override.aes = list(size=node_size))) +
        theme(legend.title = element_blank()) +
        scale_color_brewer(palette = palette) +
        geom_text(label = cn)
      )
    }

   list(plt = plt)

    # end: select.estimate
  } else if (is(x, "select.explore")){

    if(x$alternative == "two.sided" |

       x$alternative == "greater" |

       x$alternative == "less"){


      cn <- colnames(x$object$Y)

      p <- ncol(x$pcor_mat_zero)

      diag(x$pcor_mat_zero) <- 0

      if(x$alternative == "two.sided"){
        Adj_alt <- x$Adj_10
        Adj_null <- x$Adj_01
      }

      if(x$alternative == "greater" |
         x$alternative == "less"){

        warning(paste0("interpret the conditional indepedence structure cautiously, as the Bayes factor\n",
                       "is a measure of 'relative' evidence. In this case, ",
                       x$alternative, " than zero was compared\n",
                       "to a null model. This does not consider the opposite direction."
        ))

        Adj_alt <- x$Adj_20
        Adj_null <- x$Adj_02

      }

      ambiguous <- matrix(1, p, p) - diag(p) -  Adj_alt - Adj_null

      net_alt <- network::network(x$pcor_mat_zero)

      net_null <- network::network(Adj_null)

      net_ambigous <- network::network(ambiguous)

      if(is.null(cn) ) {
        cn <- 1:p
      }



      # edge weights
      network::set.edge.value(x = net_alt, attrname = "weights",
                              value = x$pcor_mat_zero)

      # edge weights absolute
      network::set.edge.value(x = net_alt, attrname = "abs_weights",
                              value = abs(x$pcor_mat_zero) * edge_magnify)


      # edge colors
      network::set.edge.attribute(x = net_alt, attrname = "edge_color",
                                  value = ifelse(net_alt %e% "weights" < 0,
                                                 neg_col,
                                                 pos_col))


      e <- abs(as.numeric( x$pcor_mat_zero))

      if(is.null(groups)){

        plt_alt <- ggnet2(
          net_alt,
          edge.alpha = e[e != 0] / max(e),
          edge.size = "abs_weights",
          edge.color = "edge_color",
          node.size = 1,
          mode = layout
        ) +
          geom_point(color = "black",
                     size = node_size + 1) +
          geom_point(size = node_size, color = "white")  +
          guides(color = FALSE) +
          geom_text(label = cn)


        plt_null <- ggnet2(net_null,
                           node.size = 1,
                           mode = layout) +
          geom_point(color = "black",
                     size = node_size + 1) +
          geom_point(size = node_size, color = "white")  +
          guides(color = FALSE) +
          geom_text(label = cn)



        plt_ambiguous <- ggnet2(net_ambigous,
                                node.size = 1,
                                mode = layout) +
          geom_point(color = "black",
                     size = node_size + 1) +
          geom_point(size = node_size, color = "white")  +
          guides(color = FALSE) +
          geom_text(label = cn)

        list(plt_alt = plt_alt,
             plt_null = plt_null,
             plt_ambiguous = plt_ambiguous)

      } else {

        net_alt %v% "group" <- groups
        net_null %v% "group" <- groups
        net_ambigous %v% "group" <- groups

       suppressMessages(
        plt_alt <- ggnet2(net_alt, edge.alpha = e[e != 0] / max(e),
               edge.size = "abs_weights",
               edge.color = "edge_color",
               node.color = "group",
               node.size = 1,
               mode = layout) +
          geom_point(aes(color = groups),
                     size = node_size + 2,
                     alpha = 0.2) +
          geom_point(aes(color = groups),
                     size = node_size,
                     alpha = 1) +
          guides(colour = guide_legend(override.aes = list(size=node_size))) +
          theme(legend.title = element_blank()) +
          scale_color_brewer(palette = palette) +
          geom_text(label = cn)
       )

       suppressMessages(

        plt_null <- ggnet2(net_null,
                           node.size = 1,
                           node.color = "group",
                           mode = layout) +
          geom_point(aes(color = groups),
                     size = node_size + 2,
                     alpha = 0.2) +
          geom_point(aes(color = groups),
                     size = node_size,
                     alpha = 1) +
          guides(colour = guide_legend(override.aes = list(size=node_size))) +
          theme(legend.title = element_blank()) +
          scale_color_brewer(palette = palette) +
          geom_text(label = cn)
        )


       suppressMessages(
       plt_ambiguous <-  ggnet2(net_ambigous,
                                 node.size = 1,
                                 node.color = "group",
                                 mode = layout) +
          geom_point(aes(color = groups),
                     size = node_size + 2,
                     alpha = 0.2) +
          geom_point(aes(color = groups),
                     size = node_size,
                     alpha = 1) +
          guides(colour = guide_legend(override.aes = list(size=node_size))) +
          theme(legend.title = element_blank()) +
          scale_color_brewer(palette = palette) +
          geom_text(label = cn)
       )


       list(H1_plt = plt_alt,
            H0_plt = plt_null,
            ambiguous_plt = plt_ambiguous)

      } # end groups


    } else if(x$alternative == "exhaustive"){

      cn <- colnames(x$object$Y)

      p <- ncol(x$pcor_mat)

      diag(x$pcor_mat) <- 0

      ambiguous <- ifelse((x$neg_mat + x$pos_mat + x$null_mat) == 0, 1, 0)

      net_pos <- network::network(x$pos_mat * x$pcor_mat)
      net_neg <- network::network(x$neg_mat * x$pcor_mat)
      net_null <- network::network(x$null_mat)
      net_ambigous <- network::network(ambiguous)

      if(is.null(cn) ) {
        cn <- 1:p
      }

      # positive
      # edge weights
      network::set.edge.value(x = net_pos, attrname = "weights",
                              value = x$pos_mat * x$pcor_mat)

      # edge weights absolute
      network::set.edge.value(x = net_pos, attrname = "abs_weights",
                              value = abs(x$pos_mat * x$pcor_mat) * edge_magnify)

      # edge colors
      network::set.edge.attribute(x = net_pos, attrname = "edge_color",
                                  value = ifelse(net_pos %e% "weights" < 0,
                                                 neg_col,
                                                 pos_col))

      # negative
      # edge weights
      network::set.edge.value(x = net_neg, attrname = "weights",
                              value = x$neg_mat * x$pcor_mat)

      # edge weights absolute
      network::set.edge.value(x = net_neg, attrname = "abs_weights",
                              value = abs(x$neg_mat * x$pcor_mat ) * edge_magnify)


      # edge colors
      network::set.edge.attribute(x = net_neg, attrname = "edge_color",
                                  value = ifelse(net_neg %e% "weights" < 0,
                                                 neg_col,
                                                 pos_col))

      if(is.null(groups)){

        e <- abs(as.numeric( x$pcor_mat * x$pos_mat))

        plt_pos <- ggnet2(
          net_pos,
          edge.alpha = e[e != 0] / max(e),
          edge.size = "abs_weights",
          edge.color = "edge_color",
          node.size = 1,
          mode = layout
        ) +
          geom_point(color = "black",
                     size = node_size + 1) +
          geom_point(size = node_size, color = "white")  +
          guides(color = FALSE) +
          geom_text(label = cn)


        e <- abs(as.numeric( x$pcor_mat * x$neg_mat))

        plt_neg <- ggnet2(net_neg,
                          node.size = 1,
                          edge.alpha = e[e != 0] / max(e),
                          edge.size = "abs_weights",
                          edge.color = "edge_color",
                          mode = layout) +
          geom_point(color = "black",
                     size = node_size + 1) +
          geom_point(size = node_size, color = "white")  +
          guides(color = FALSE) +
          geom_text(label = cn)



        plt_null <- ggnet2(net_null,
                           node.size = 1,
                           mode = layout) +
          geom_point(color = "black",
                     size = node_size + 1) +
          geom_point(size = node_size, color = "white")  +
          guides(color = FALSE) +
          geom_text(label = cn)

        plt_ambiguous <- ggnet2(net_ambigous,
                                node.size = 1,
                                mode = layout) +
          geom_point(color = "black",
                     size = node_size + 1) +
          geom_point(size = node_size, color = "white")  +
          guides(color = FALSE) +
          geom_text(label = cn)

        list(plt_pos = plt_pos,
             plt_neg = plt_neg,
             plt_null = plt_null,
             plt_ambiguous = plt_ambiguous)

      } else {

        net_pos %v% "group" <- groups
        net_neg %v% "group" <- groups
        net_null %v% "group" <- groups
        net_ambigous %v% "group" <- groups


        e <- abs(as.numeric( x$pcor_mat * x$pos_mat))

        suppressMessages(
        plt_pos <- ggnet2(
          net_pos,
          edge.alpha = e[e != 0] / max(e),
          edge.size = "abs_weights",
          edge.color = "edge_color",
          node.size = 1,
          node.color = "group",
          mode = layout
        ) +
          geom_point(aes(color = groups),
                     size = node_size + 2,
                     alpha = 0.2) +
          geom_point(aes(color = groups),
                     size = node_size,
                     alpha = 1) +
          guides(colour = guide_legend(override.aes = list(size=node_size))) +
          theme(legend.title = element_blank()) +
          scale_color_brewer(palette = palette) +
          geom_text(label = cn)
        )


        e <- abs(as.numeric( x$pcor_mat * x$neg_mat))


        suppressMessages(
        plt_neg <- ggnet2(net_neg,
                          node.size = 1,
                          edge.alpha = e[e != 0] / max(e),
                          edge.size = "abs_weights",
                          edge.color = "edge_color",
                          node.color = "group",
                          mode = layout) +
          geom_point(aes(color = groups),
                     size = node_size + 2,
                     alpha = 0.2) +
          geom_point(aes(color = groups),
                     size = node_size,
                     alpha = 1) +
          guides(colour = guide_legend(override.aes = list(size=node_size))) +
          theme(legend.title = element_blank()) +
          scale_color_brewer(palette = palette) +
          geom_text(label = cn)
        )


       suppressMessages(
        plt_null <- ggnet2(net_null,
                           node.size = 1,
                           mode = layout,
                           node.color = "group") +
          geom_point(aes(color = groups),
                     size = node_size + 2,
                     alpha = 0.2) +
          geom_point(aes(color = groups),
                     size = node_size,
                     alpha = 1) +
          guides(colour = guide_legend(override.aes = list(size=node_size))) +
          theme(legend.title = element_blank()) +
          scale_color_brewer(palette = palette) +
          geom_text(label = cn)
)

       suppressMessages(
       plt_ambiguous <- ggnet2(net_ambigous,
                                node.size = 1,
                                mode = layout,
                                node.color =  "group") +
          geom_point(aes(color = groups),
                     size = node_size,
                     alpha = 1) +
          guides(colour = guide_legend(override.aes = list(size=node_size))) +
          theme(legend.title = element_blank()) +
          scale_color_brewer(palette = palette) +
          geom_text(label = cn)
       )

        list(plt_pos = plt_pos,
             plt_neg = plt_neg,
             plt_null = plt_null,
             plt_ambiguous = plt_ambiguous)

      } # end groups
    }
  } else if(is(x, "select.ggm_compare_estimate")){


    cn <- colnames(x$object$info$dat[[1]])

    p <- ncol(x$pcor_adj[[1]])

    comparisons <- length(x$pcor_adj)

    if(is.null(cn) ) {
      cn <- 1:p
    }



  lapply(1:comparisons, function(z){

            net <- network::network(x$pcor_adj[[z]])

            # edge weights
            network::set.edge.value(x = net, attrname = "weights",
                                    value = x$pcor_adj[[z]])

            # edge weights absolute
            network::set.edge.value(x = net, attrname = "abs_weights",
                                    value = abs(x$pcor_adj[[z]]) * edge_magnify)


            # edge colors
            network::set.edge.attribute(x = net, attrname = "edge_color",
                                        value = ifelse(net %e% "weights" < 0,
                                                       neg_col,
                                                       pos_col))

            diag(x$pcor_adj[[z]]) <- 0
            e <- abs(as.numeric(x$pcor_adj[[z]]))

            if(is.null(groups)){

              ggnet2(net, edge.alpha = e[e != 0] / max(e),
                            edge.size = "abs_weights",
                            edge.color = "edge_color",
                            node.size = 1,
                            mode = layout) +
                geom_point(color = "black",
                           size = node_size+1) +
                geom_point(size = node_size, color = "white")  +
                guides(color = FALSE) +
                geom_text(label = cn) +
                ggtitle(names(x$object$diff)[z])


            } else {

              net %v% "group" <- groups

              suppressMessages(
                 ggnet2(net, edge.alpha = e[e != 0] / max(e),
                               edge.size = "abs_weights",
                               edge.color = "edge_color",
                               node.color = "group",
                               node.size = 1,
                               mode = layout) +
                  geom_point(aes(color = groups),
                             size = node_size + 2,
                             alpha = 0.2) +
                  geom_point(aes(color = groups),
                             size = node_size,
                             alpha = 1) +
                  guides(colour = guide_legend(override.aes = list(size=node_size))) +
                  theme(legend.title = element_blank()) +
                  scale_color_brewer(palette = palette) +
                  geom_text(label = cn) +
                  ggtitle(names(x$object$diff[z]))
              )
            }
        })
  } else if(is(x, "select.ggm_compare_bf")){


    if(x$post_prob == 0.50){


      cn <- colnames(x$object$info$dat[[1]])



      p <- ncol(x$pcor_adj[[1]])

     if(is.null(cn) ){
        cn <- 1:p
      }



      if(length(x$info$dat) == 2){


        net_alt <- network::network(x$adj_10 * x$pcor_mat)
        net_null <- network::network(x$adj_01)


        # edge weights
        network::set.edge.value(x = net_alt, attrname = "weights",
                                value = x$adj_10 * x$pcor_mat)

        # edge weights absolute
        network::set.edge.value(x = net_alt, attrname = "abs_weights",
                                value = abs(x$adj_10 * x$pcor_mat) * edge_magnify)

        # edge colors
        network::set.edge.attribute(x = net_alt, attrname = "edge_color",
                                    value = ifelse(net_alt %e% "weights" < 0,
                                                   neg_col,
                                                   pos_col))


        if(is.null(groups)){

          e <- abs(as.numeric( x$pcor_mat * x$adj_10))

          plt_alt <- ggnet2(
            net_alt,
            edge.alpha = e[e != 0] / max(e),
            edge.size = "abs_weights",
            edge.color = "edge_color",
            node.size = 1,
            mode = layout,...
          ) +
            geom_point(color = "black",
                       size = node_size + 1) +
            geom_point(size = node_size, color = "white")  +
            guides(color = FALSE) +
            geom_text(label = cn)

          plt_null <- ggnet2(net_null,
                             node.size = 1,
                             mode = layout) +
            geom_point(color = "black",
                       size = node_size + 1) +
            geom_point(size = node_size, color = "white")  +
            guides(color = FALSE) +
            geom_text(label = cn)



          list(plt_alt = plt_alt,
               plt_null = plt_null)

        } else {

          net_alt %v% "group" <- groups
          net_null %v% "group" <- groups


          e <- abs(as.numeric( x$pcor_mat * x$adj_10))

          suppressMessages(
          plt_alt <- ggnet2(
            net_alt,
            edge.alpha = e[e != 0] / max(e),
            edge.size = "abs_weights",
            edge.color = "edge_color",
            node.color = "group",
            node.size = 1,
            mode = layout
          ) +
            geom_point(aes(color = groups),
                     size = node_size + 2,
                     alpha = 0.2) +
            geom_point(aes(color = groups),
                       size = node_size,
                       alpha = 1) +
            guides(colour = guide_legend(override.aes = list(size=node_size))) +
            theme(legend.title = element_blank()) +
            scale_color_brewer(palette = palette) +
            geom_text(label = cn)
          )

        suppressMessages(
          plt_null <- ggnet2(net_null,
                             node.size = 1,
                             mode = layout,
                             node.color = "group") +
            geom_point(color = "black",
                       size = node_size + 1) +
            geom_point(size = node_size, color = "white")  +
            guides(color = FALSE) +
            geom_text(label = cn) +
            geom_point(aes(color = groups),
                       size = node_size + 2,
                       alpha = 0.2) +
            geom_point(aes(color = groups),
                       size = node_size,
                       alpha = 1) +
            guides(colour = guide_legend(override.aes = list(size=node_size))) +
            theme(legend.title = element_blank()) +
            scale_color_brewer(palette = palette) +
            geom_text(label = cn)
          )


        list(plt_alt = plt_alt,
             plt_null = plt_null)

        } # end clusters

        }  else {

          net_alt <- network::network(x$adj_10)
          net_null <- network::network(x$adj_01)

          if(is.null(groups)){

            plt_alt <- ggnet2(net_alt,
                             node.size = 1,
                             mode = layout) +
            geom_point(color = "black",
                       size = node_size + 1) +
            geom_point(size = node_size, color = "white")  +
            guides(color = FALSE) +
            geom_text(label = cn)

            plt_null <- ggnet2(net_null,
                             node.size = 1,
                             mode = layout) +
            geom_point(color = "black",
                       size = node_size + 1) +
            geom_point(size = node_size, color = "white")  +
            guides(color = FALSE) +
            geom_text(label = cn)



          list(plt_alt = plt_alt,
               plt_null = plt_null)

        } else {

          net_alt %v% "group" <- groups
          net_null %v% "group" <- groups



          suppressMessages(
            plt_alt <- ggnet2(net_alt,
                               node.size = 1,
                               mode = layout,
                               node.color = "group") +
              geom_point(color = "black",
                         size = node_size + 1) +
              geom_point(size = node_size, color = "white")  +
              guides(color = FALSE) +
              geom_text(label = cn) +
              geom_point(aes(color = groups),
                         size = node_size + 2,
                         alpha = 0.2) +
              geom_point(aes(color = groups),
                         size = node_size,
                         alpha = 1) +
              guides(colour = guide_legend(override.aes = list(size=node_size))) +
              theme(legend.title = element_blank()) +
              scale_color_brewer(palette = palette) +
              geom_text(label = cn)
          )

          suppressMessages(
            plt_null <- ggnet2(net_null,
                               node.size = 1,
                               mode = layout,
                               node.color = "group") +
              geom_point(color = "black",
                         size = node_size + 1) +
              geom_point(size = node_size, color = "white")  +
              guides(color = FALSE) +
              geom_text(label = cn) +
              geom_point(aes(color = groups),
                         size = node_size + 2,
                         alpha = 0.2) +
              geom_point(aes(color = groups),
                         size = node_size,
                         alpha = 1) +
              guides(colour = guide_legend(override.aes = list(size=node_size))) +
              theme(legend.title = element_blank()) +
              scale_color_brewer(palette = palette) +
              geom_text(label = cn)
          )


          list(plt_alt = plt_alt,
               plt_null = plt_null)

      } # end of clusters

    } # end two groups

      # more than 0.50
  } else {


    cn <- colnames(x$object$info$dat[[1]])



    p <- ncol(x$BF_10)

    if(is.null(cn) ){
      cn <- 1:p
    }



    if(length(x$info$dat) == 2){

      Adj_alt <- x$adj_10
      Adj_null <- x$adj_01

      ambiguous <- matrix(1, p, p) - diag(p) -  Adj_alt - Adj_null

      net_alt <- network::network(x$pcor_mat_10 * Adj_alt)
      net_null <- network::network(Adj_null)
      net_ambigous <- network(ambiguous)



      # edge weights
      network::set.edge.value(x = net_alt, attrname = "weights",
                              value = x$pcor_mat_10 * Adj_alt)

      # edge weights absolute
      network::set.edge.value(x = net_alt, attrname = "abs_weights",
                              value = abs(x$pcor_mat_10 * Adj_alt) * edge_magnify)

      # edge colors
      network::set.edge.attribute(x = net_alt, attrname = "edge_color",
                                  value = ifelse(net_alt %e% "weights" < 0,
                                                 neg_col,
                                                 pos_col))

      e <- abs(as.numeric(  x$pcor_mat_10 * x$adj_10))

      if(is.null(groups)){

        plt_alt <- ggnet2(
          net_alt,
          edge.alpha = e[e != 0] / max(e),
          edge.size = "abs_weights",
          edge.color = "edge_color",
          node.size = 1,
          mode = layout
        ) +
          geom_point(color = "black",
                     size = node_size + 1) +
          geom_point(size = node_size, color = "white")  +
          guides(color = FALSE) +
          geom_text(label = cn)


        plt_null <- ggnet2(net_null,
                           node.size = 1,
                           mode = layout) +
          geom_point(color = "black",
                     size = node_size + 1) +
          geom_point(size = node_size, color = "white")  +
          guides(color = FALSE) +
          geom_text(label = cn)



        plt_ambiguous <- ggnet2(net_ambigous,
                                node.size = 1,
                                mode = layout) +
          geom_point(color = "black",
                     size = node_size + 1) +
          geom_point(size = node_size, color = "white")  +
          guides(color = FALSE) +
          geom_text(label = cn)

        list(plt_alt = plt_alt,
             plt_null = plt_null,
             plt_ambiguous = plt_ambiguous)

      } else {

        net_alt %v% "group" <- groups
        net_null %v% "group" <- groups
        net_ambigous %v% "group" <- groups

        suppressMessages(
          plt_alt <- ggnet2(net_alt, edge.alpha = e[e != 0] / max(e),
                            edge.size = "abs_weights",
                            edge.color = "edge_color",
                            node.color = "group",
                            node.size = 1,
                            mode = layout) +
            geom_point(aes(color = groups),
                       size = node_size + 2,
                       alpha = 0.2) +
            geom_point(aes(color = groups),
                       size = node_size,
                       alpha = 1) +
            guides(colour = guide_legend(override.aes = list(size=node_size))) +
            theme(legend.title = element_blank()) +
            scale_color_brewer(palette = palette) +
            geom_text(label = cn)
        )

        suppressMessages(

          plt_null <- ggnet2(net_null,
                             node.size = 1,
                             node.color = "group",
                             mode = layout) +
            geom_point(aes(color = groups),
                       size = node_size + 2,
                       alpha = 0.2) +
            geom_point(aes(color = groups),
                       size = node_size,
                       alpha = 1) +
            guides(colour = guide_legend(override.aes = list(size=node_size))) +
            theme(legend.title = element_blank()) +
            scale_color_brewer(palette = palette) +
            geom_text(label = cn)
        )


        suppressMessages(
          plt_ambiguous <-  ggnet2(net_ambigous,
                                   node.size = 1,
                                   node.color = "group",
                                   mode = layout) +
            geom_point(aes(color = groups),
                       size = node_size + 2,
                       alpha = 0.2) +
            geom_point(aes(color = groups),
                       size = node_size,
                       alpha = 1) +
            guides(colour = guide_legend(override.aes = list(size=node_size))) +
            theme(legend.title = element_blank()) +
            scale_color_brewer(palette = palette) +
            geom_text(label = cn)
        )


        list(H1_plt = plt_alt,
             H0_plt = plt_null,
             ambiguous_plt = plt_ambiguous)

        } # end cluster

      # end: two groups
  } else {

      Adj_alt <- x$adj_10
      Adj_null <- x$adj_01

      ambiguous <- matrix(1, p, p) - diag(p) -  Adj_alt - Adj_null

      net_alt <- network::network(Adj_alt)
      net_null <- network::network(Adj_null)
      net_ambigous <- network::network(ambiguous)

      if(is.null(groups)){

        plt_alt <- ggnet2(net_alt,
                           node.size = 1,
                           mode = layout) +
          geom_point(color = "black",
                     size = node_size + 1) +
          geom_point(size = node_size, color = "white")  +
          guides(color = FALSE) +
          geom_text(label = cn)


        plt_null <- ggnet2(net_null,
                           node.size = 1,
                           mode = layout) +
          geom_point(color = "black",
                     size = node_size + 1) +
          geom_point(size = node_size, color = "white")  +
          guides(color = FALSE) +
          geom_text(label = cn)

        plt_ambiguous <-  ggnet2(net_ambigous,
                                 node.size = 1,
                                 mode = layout) +
          geom_point(color = "black",
                     size = node_size + 1) +
          geom_point(size = node_size, color = "white")  +
          guides(color = FALSE) +
          geom_text(label = cn)


        list(plt_alt = plt_alt,
             plt_null = plt_null,
             plt_ambiguous = plt_ambiguous)

      } else {

        net_alt %v% "group" <- groups
        net_null %v% "group" <- groups
        net_ambigous %v% "group" <- groups

        suppressMessages(
          plt_alt <- ggnet2(net_alt,
                            node.size = 1,
                            node.color = "group",
                            mode = layout) +
            geom_point(aes(color = groups),
                       size = node_size + 2,
                       alpha = 0.2) +
            geom_point(aes(color = groups),
                       size = node_size,
                       alpha = 1) +
            guides(colour = guide_legend(override.aes = list(size=node_size))) +
            theme(legend.title = element_blank()) +
            scale_color_brewer(palette = palette) +
            geom_text(label = cn)
        )

        suppressMessages(

          plt_null <- ggnet2(net_null,
                             node.size = 1,
                             node.color = "group",
                             mode = layout) +
            geom_point(aes(color = groups),
                       size = node_size + 2,
                       alpha = 0.2) +
            geom_point(aes(color = groups),
                       size = node_size,
                       alpha = 1) +
            guides(colour = guide_legend(override.aes = list(size=node_size))) +
            theme(legend.title = element_blank()) +
            scale_color_brewer(palette = palette) +
            geom_text(label = cn)
        )

        suppressMessages(

          plt_ambiguous <-  ggnet2(net_ambigous,
                                   node.size = 1,
                                   node.color = "group",
                                   mode = layout) +
            geom_point(aes(color = groups),
                       size = node_size + 2,
                       alpha = 0.2) +
            geom_point(aes(color = groups),
                       size = node_size,
                       alpha = 1) +
            guides(colour = guide_legend(override.aes = list(size=node_size))) +
            theme(legend.title = element_blank()) +
            scale_color_brewer(palette = palette) +
            geom_text(label = cn)

          )

        list(H1_plt = plt_alt,
             H0_plt = plt_null,
             ambiguous_plt = plt_ambiguous)

        } # end cluster

      } # end: more than 2 groups

    } # end not 0.50.

  } # end select.ggm_compare.explore

}

