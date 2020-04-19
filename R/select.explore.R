#' Select Graphical Structure with the Bayes Factor
#' @name select.explore
#' @description This allows for not only estimating the conditional dependence structure, that is non-zero edges, but also the conditional \strong{in}dependence
#' structure (evidence for no relation).
#'
#' @param object object of class \code{explore.default}
#' @param BF_cut evidentiary threshold
#' @param alternative type of hypothesis (see notes)
#' @param ... currently not used
#'
#' @note The \code{alternative} can be either \code{greater}, \code{less}, \code{two.sided}, or \code{exhaustive}. The argument \code{hyp_prob} is used only
#' when \code{alternative = hypothesis}. \code{greater} and \code{less} test directional hypotheses, and thus, the graphical structure will
#' only included edges in that direction (i.e., positive or negative). \code{two.sided} is the customary approach, and test for the presence or
#' absence of an edge. \code{exhaustive} tests negative vs. positive  vs. zero. Here \code{hyp_prob} is the posterior probability threshold for the respective
#' hypotheses.
#'
#'
#' @return list of class \code{select.explore}:
#'
#' \code{alternative = "two.sided"}:
#' \itemize{
#'  \item \code{partials_non_zero} selected partial correlation matrix
#'  \item \code{pcor_mat} partial correlation matrix (non set to zero)
#'  \item \code{pcor_sd} partial correlation standard deviations
#'  \item \code{Adj_10} adjacency matrix for the selected edges (in favor of the \code{alternative})
#'  \item \code{Adj_01} adjacency matrix for the null hypothesis  (conditional independence)
#'  \item \code{BF_10} Bayes factors for \code{alternative}
#'  \item \code{BF_01} Bayes factors for the null hypothesis
#'  \item \code{BF_cut} evidentiary threshold
#'  \item \code{alternative} \code{"two.sided"}
#'  \item \code{code} \code{match.call()}
#' }
#'
#' \code{alternative = "greater"}:
#' \itemize{
#'  \item \code{partials_positive} selected partial correlation matrix
#'  \item \code{pcor_mat} partial correlation matrix (none set to zero)
#'  \item \code{pcor_sd} partial correlation standard deviations
#'  \item \code{Adj_20} adjacency matrix for the selected edges (in favor of the \code{alternative})
#'  \item \code{Adj_01} adjacency matrix for the null hypothesis  (conditional independence)
#'  \item \code{BF_20} Bayes factors for \code{alternative}
#'  \item \code{BF_01} Bayes factors for the null hypothesis
#'  \item \code{BF_cut} evidentiary threshold
#'  \item \code{alternative} \code{"greater"}
#'  \item \code{code} \code{match.call()}
#' }
#'
#' \code{alternative = "less"}:
#' \itemize{
#'  \item \code{partials_negative} selected partial correlation matrix
#'  \item \code{pcor_mat} partial correlation matrix (none set to zero)
#'  \item \code{pcor_sd} partial correlation standard deviations
#'  \item \code{Adj_20} adjacency matrix for the selected edges (in favor of the \code{alternative})
#'  \item \code{Adj_01} adjacency matrix for the null hypothesis  (conditional independence)
#'  \item \code{BF_20} Bayes factors for \code{alternative}
#'  \item \code{BF_01} Bayes factors for the null hypothesis
#'  \item \code{BF_cut} evidentiary threshold
#'  \item \code{alternative} \code{"less"}
#'  \item \code{code} \code{match.call()}
#' }
#'
#' \code{alternative = "exhaustive"}
#' \itemize{
#' \item \code{post_prob} data.frame with posterior probabilities for each edge
#' \item \code{neg_mat} adjacency matrix for negative edges
#' \item \code{post_mat} adjacency matrix for positive edges
#' \item \code{null_mat} adjacency matrix for zero  (conditional independence)
#' \item \code{"alternative"} "exhaustive"
#' \item \code{pcor_mat} partial correlation matrix (non set to zero)
#' \item \code{pcor_sd} partial correlation standard deviations
#' \item \code{code} \code{match.call()}
#' \item \code{prob} \code{hyp_prob}
#' }
#'
#' @examples
#' \donttest{
#' # p = 10
#' Y <- BGGM::bfi[,1:10]
#'
#' # sample posterior
#' fit <- explore(Y, iter = 5000)
#'
#' # select E
#' E <- select(fit, BF_cut = 3)
#'
#' # summarize
#' summary(E)
#'
#' # non-zero edges
#' E$partials_non_zero
#'
#' # adjacency matrix
#' E$Adj_10
#'
#' # null adjacency matrix
#' E$Adj_01
#' }
#' @export
select.explore <- function(object,
                           BF_cut = 3,
                           alternative = "two.sided", ...){
  # rename
  x <- object

  # hyp probability
  hyp_prob <- BF_cut / (BF_cut + 1)

  # posterior samples
  posterior_samples <- do.call(rbind.data.frame,
                               lapply(1:x$cores, function(z)
                                 x$samples[[z]]$fisher_z_post))

  # prior samples
  prior_samples <- unlist(do.call(rbind.data.frame,
                                  lapply(1:x$cores, function(z)
                                  x$samples[[z]]$fisher_z_prior)))

  # two sided testing
  if(alternative == "two.sided"){

    # matrices for storage
    BF_10_mat <- BF_01_mat <- matrix(0, x$p, x$p)

    # prior density
    prior_dens <- dnorm(0, mean(prior_samples), sd(prior_samples))

    # BF for alternative
    BF_10 <- apply(posterior_samples, MARGIN = 2,
                   FUN = function(z){ prior_dens / dnorm(0, mean(z), sd(z)) })

    # BF for null
    BF_01 <- 1 / BF_10

    # alternative BF mat
    BF_10_mat[upper.tri(BF_10_mat)] <- BF_10[1:x$edge]
    BF_10_mat <- symmteric_mat(BF_10_mat)

    # null BF mat
    BF_01_mat <- 1 / BF_10_mat
    diag(BF_01_mat) <- 0

    # selected edges (alternative)
    Adj_10 <- ifelse(BF_10_mat > BF_cut, 1, 0)

    # selected edges (null)
    Adj_01 <- ifelse(BF_10_mat < 1 / BF_cut, 1, 0)
    diag(Adj_01) <- 0

    # returned object
    returned_object = list(partials_non_zero = x$parcors_mat * Adj_10,
                           pcor_mat = round(x$parcors_mat,3),
                           pcor_sd = round(x$parcors_sd,3),
                           Adj_10 = Adj_10,
                           Adj_01 = Adj_01,
                           BF_10 = BF_10_mat,
                           BF_01 = BF_01_mat,
                           BF_cut = BF_cut,
                           alternative = alternative,
                           call = match.call())
  }

  # one sided (greater)
  if(alternative == "greater"){

    # prior density
    prior_dens <- dnorm(0, mean(prior_samples), sd(prior_samples))

    # matrices for storage
    BF_20_mat <- BF_01_mat <- matrix(0, x$p, x$p)

    # BF for alternative
    BF_10 <- apply(posterior_samples, MARGIN = 2,
                   FUN = function(z){ prior_dens / dnorm(0, mean(z), sd(z)) })

    # density greater than zero
    dens_greater <-  apply(posterior_samples, MARGIN = 2,
                           FUN = function(z){ (1  - pnorm(0, mean(z), sd(z))) * 2 })

    # one sided BF
    BF_20 <- dens_greater * (BF_10)

    # BF null
    BF_01 <- 1 / BF_10

    # alternative BF mat
    BF_20_mat[upper.tri(BF_20_mat)] <- BF_20[1:x$edge]
    BF_20_mat <- symmteric_mat(BF_20_mat)

    # null BF mat
    BF_01_mat[upper.tri(BF_01_mat)] <- BF_01[1:x$edge]
    BF_01_mat <- symmteric_mat(BF_01_mat)
    diag(BF_01_mat) <- 0

    # selected edges (alternative)
    Adj_20 <- ifelse(BF_20_mat > BF_cut, 1, 0)

    # selected edges (null)
    Adj_01 <- ifelse(BF_01_mat > BF_cut, 1, 0)
    diag(Adj_01) <- 0

    # returned object
    returned_object = list(partials_positive = x$parcors_mat * Adj_20,
                           pcor_mat = round(x$parcors_mat,3),
                           pcor_sd = round(x$parcors_sd,3),
                           Adj_01 = Adj_01,
                           Adj_20 =  Adj_20,
                           BF_20 = BF_20_mat,
                           BF_01 = BF_01_mat,
                           BF_cut = BF_cut,
                           alternative = alternative,
                           call = match.call())


  }

  # one sided (less)
  if(alternative == "less")  {

    # prior density
    prior_dens <- dnorm(0, mean(prior_samples), sd(prior_samples))

    # matrices for storage
    BF_20_mat <- BF_01_mat <- matrix(0, x$p, x$p)

    # BF for alternative
    BF_10 <- apply(posterior_samples, MARGIN = 2,
                   FUN = function(z){ prior_dens / dnorm(0, mean(z), sd(z)) })

    # density less than zero
    dens_less <-  apply(posterior_samples, MARGIN = 2,
                        FUN = function(z){ (pnorm(0, mean(z), sd(z))) * 2 })

    # one sided BF
    BF_20 <- dens_less * (BF_10)

    # BF null
    BF_01 <- 1 / BF_10

    # alternative BF mat
    BF_20_mat[upper.tri(BF_20_mat)] <- BF_20[1:x$edge]
    BF_20_mat <- symmteric_mat(BF_20_mat)

    # null BF mat
    BF_01_mat[upper.tri(BF_01_mat)] <- BF_01[1:x$edge]
    BF_01_mat <- symmteric_mat(BF_01_mat)
    diag(BF_01_mat) <- 0

    # selected edges (alternative)
    Adj_20 <- ifelse(BF_20_mat > BF_cut, 1, 0)

    # selected edges (null)
    Adj_01 <- ifelse(BF_01_mat > BF_cut, 1, 0)
    diag(Adj_01) <- 0

    # returned object
    returned_object = list(partials_negative = x$parcors_mat * Adj_20,
                           pcor_mat = round(x$parcors_mat,3),
                           pcor_sd = round(x$parcors_sd, 3),
                           Adj_01 = Adj_01,
                           Adj_20 = Adj_20,
                           BF_20 = BF_20_mat,
                           BF_01 = BF_01_mat,
                           BF_cut = BF_cut,
                           alternative = alternative,
                           call = match.call())

  }

  # exhaustive testing
  if(alternative == "exhaustive"){

    if(is.null(hyp_prob)){
      stop("posterior probability must be specificed \n for exhaustive hypothesis testing")

    }

    # matrices for storage
    mat_names <- null_mat <- pos_mat <- neg_mat <- matrix(0, x$p, x$p)

    # matrix names
    mat_names[] <-  unlist(lapply(1:x$p,
                                  FUN = function(z) paste(1:x$p, z, sep = "--")))

    # prior density
    prior_dens <- dnorm(0, mean(prior_samples), sd(prior_samples))

    # two sided BF
    BF_10 <- apply(posterior_samples, MARGIN = 2,
                   FUN = function(z){ prior_dens / dnorm(0, mean(z), sd(z)) })

    # density less than zero
    dens_less <-  apply(posterior_samples, MARGIN = 2,
                        FUN = function(z){ (pnorm(0, mean(z), sd(z))) * 2 })

    # BF less than zero
    BF_less <- dens_less * (BF_10)

    # BF null
    BF_01 <- 1 / BF_10

    # density greater than zero
    dens_greater <-  apply(posterior_samples, MARGIN = 2,
                           FUN = function(z){ (1  - pnorm(0, mean(z), sd(z))) * 2 })

    # BF greater than zero
    BF_greater <- dens_greater * (BF_10)

    # combine BFs
    BF_mat <- data.frame(BF_01 = BF_01,
                         BF_greater = BF_greater,
                         BF_less = BF_less)

    # posterior hypothesis probabilities
    post_prob <-  data.frame(t(apply(BF_mat, 1, FUN = function(x) { round(x / sum(x),3) })))

    # column names
    colnames(post_prob) <- c("prob_zero", "prob_greater", "prob_less")

    # no rownames
    row.names(post_prob) <- c()

    # round
    post_prob <-  round(post_prob, 3)

    # combine probabilities and edge names
    post_prob <- cbind.data.frame(edge = mat_names[upper.tri(mat_names)], post_prob)[1:x$edge,]

    # selected edges (null)
    null_mat[upper.tri(null_mat)] <- ifelse(post_prob$prob_zero > hyp_prob, 1, 0)
    null_mat <- symmteric_mat(null_mat)

    # selected edges (positive)
    pos_mat[upper.tri(pos_mat)] <-  ifelse(post_prob$prob_greater >  hyp_prob, 1, 0)
    pos_mat <- symmteric_mat(pos_mat)

    # selected edges (negative)
    neg_mat[upper.tri(neg_mat)]  <-  ifelse(post_prob$prob_less > hyp_prob, 1 , 0)
    neg_mat <- symmteric_mat(neg_mat)

    # negative edges
    returned_object <- list(post_prob = post_prob,
                            neg_mat = neg_mat,
                            pos_mat = pos_mat,
                            null_mat = null_mat,
                            alternative = alternative,
                            pcor_mat = round(x$parcors_mat,3),
                            pcor_sd = round(x$parcors_sd, 3),
                            call = match.call(),
                            prob = hyp_prob)
    }

  class(returned_object) <- c("BGGM",
                              "explore",
                              "select.explore")
  returned_object
}


#' Plot \code{select.explore} Network
#'
#' @param x object of class \code{select.explore}
#' @param layout network layout (\link[sna]{gplot.layout})
#' @param edge_colors color theme for positive and negative edges
#' @param node_labels node labels
#' @param node_labels_color node labels color
#' @param node_groups node group indicator
#' @param node_outer_size node border size
#' @param node_inner_size node size
#' @param alpha edge transparency
#' @param txt_size node text size
#' @param edge_multiplier constant to change edge width (egde * edge_multiplier)
#' @param ... additional arguments (\link[GGally]{ggnet2})
#'
#' @importFrom GGally ggnet2
#' @importFrom ggplot2 ggtitle
#' @importFrom network network.vertex.names<- set.edge.value set.edge.attribute %e% %v%<-
#' @importFrom sna gplot.layout.circle
#' @return object of class \code{ggplot}
#' @export
#' @examples
#' \donttest{
#' Y <- BGGM::bfi[1:500, 1:20]
#'
#' # fit model
#' fit_explore <- explore(Y)
#'
#' # select the graph (edge set E)
#'
#' E <- select(fit_explore)
#'
#'
#' # plot
#' plt <- plot(E,
#'             node_labels = letters[1:20],
#'             node_labels_color = "white",
#'             node_groups = rep(c("1", "2", "3", "4"), each = 5),
#'             edge_colors = "classic", txt_size = 8,
#'             alpha = 0.5, palette = "Set2")
#'}
plot.select.explore <-
  function(x,
           layout = "circle",
           edge_colors = "classic",
           node_labels = NULL,
           node_labels_color = "black",
           node_groups = NULL,
           node_outer_size = 12,
           node_inner_size = 11,
           alpha = 0.50, txt_size = 8,
           edge_multiplier = 1,
           ...) {
    label <- NULL
    color <- NULL
    # number of nodes

    if(x$alternative != "exhaustive"){
    names(x)[1] <- "partials_non_zero"
    p <- ncol(x$partials_non_zero)
    # network
    if(x$alternative == "greater" | x$alternative == "less" ){
      net <- network::network(x$Adj_20)
    }else{
    net <- network::network(x$Adj_10)
}

    # default labels
    if (is.null(node_labels)) {
      network::network.vertex.names(net) <- 1:p
      # custom labels
    } else {
      # check label length
      if (length(node_labels) != p) {
        stop("labels must be of length p (number of nodes)")
      }

      network::network.vertex.names(net) <- node_labels
    }
    # set edge weights
    network::set.edge.value(
      x = net,
      attrname =  "weights",
      value = x$partials_non_zero
    )
    #
    network::set.edge.value(
      x = net,
      attrname =  "abs_weights",
      value = abs(x$partials_non_zero) * edge_multiplier
    )
    if (edge_colors == "classic") {
      network::set.edge.attribute(
        x = net,
        attrname = "edge_color",
        value = ifelse(net %e% "weights" < 0,
                       "brown3", "palegreen3")
      )
    } else if (edge_colors == "color_blind") {
      network::set.edge.attribute(
        x = net,
        attrname = "edge_color",
        value = ifelse(net %e% "weights" < 0,
                       "#009E73", "#D55E00")
      )
    } else if (edge_colors == "vivid") {
      network::set.edge.attribute(
        x = net,
        attrname = "edge_color",
        value = ifelse(net %e% "weights" < 0,
                       "darkorange1", "darkorchid4")
      )

    }

    if (is.null(node_groups)) {
      plt <- ggnet2(
        net = net, edge.alpha = alpha,
        mode = layout,
        node.size = node_outer_size,
        node.color = "black",
        edge.color = "edge_color",
        edge.size = "abs_weights",
        label = TRUE) +
        geom_point(color = "white",
                   size = node_inner_size,
                   alpha = 1) +
        geom_text(aes(label = label),
                  color = node_labels_color,
                  size = txt_size)

      plt <- list(plt = plt)

    } else {
      if (length(node_groups) != p) {
        stop("labels must be of length p (number of nodes)")
      }

      net %v% "group" <- node_groups

      plt <- ggnet2(
        net = net,edge.alpha = alpha,
        mode = layout, node.size = node_outer_size,
        node.color = "group", node.alpha = 0.5,
        edge.color = "edge_color",
        edge.size = "abs_weights",
        label = TRUE,
        ...
      ) +
        # geom_point(aes(color = color),
        #            size = node_outer_size,
        #            alpha = 0.5) +
        geom_point(aes(color = color),
                   size = node_inner_size,
                   alpha = 1) +
        geom_text(aes(label = label),
                  color = node_labels_color,
                  size = txt_size)

      plt <- list(plt = plt)
    }

    if (is.null(x$rope)) {
      # network
      net <- network::network(x$Adj_01, directed = FALSE)

      # default labels
      if (is.null(node_labels)) {
        network::network.vertex.names(net) <- 1:p
        # custom labels
      } else {
        # check label length
        if (length(node_labels) != p) {
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
          # geom_point(color = "black",
          #            size = node_outer_size,
          #            alpha = 1) +
          geom_point(color = "white",
                     size = node_inner_size,
                     alpha = 1) +
          geom_text(aes(label = label),
                    color = node_labels_color,
                    size = txt_size)

        plt$plt_null <- plt_null

      } else{
        if (length(node_groups) != p) {
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

        plt$plt_null <- plt_null
      }
    }
  } else{


      p <-  ncol(x$neg_mat)

      if (is.null(node_groups)) {

        # positive
        w_pos <- x$pos_mat * x$pcor_mat
        net <- network::network(x$pos_mat, directed = FALSE)

      # default labels
      if (is.null(node_labels)) {
        network::network.vertex.names(net) <- 1:p
        # custom labels
      } else {
        # check label length
        if (length(node_labels) != p) {
          stop("labels must be of length p (number of nodes)")
        }

        network::network.vertex.names(net) <- node_labels
      }
        network::set.edge.value(x = net,
                                attrname =  "weights",
                                value = w_pos)

        network::set.edge.value(x = net,
                                attrname =  "abs_weights",
                                value = abs(w_pos) * edge_multiplier)

      if (edge_colors == "classic") {
        network::set.edge.attribute(
          x = net,
          attrname = "edge_color",
          value = ifelse(net %e% "weights" < 0,
                         "brown3", "palegreen3")
        )
      } else if (edge_colors == "color_blind") {
        network::set.edge.attribute(
          x = net,
          attrname = "edge_color",
          value = ifelse(net %e% "weights" < 0,
                         "#009E73", "#D55E00")
        )
      } else if (edge_colors == "vivid") {
        network::set.edge.attribute(
          x = net,
          attrname = "edge_color",
          value = ifelse(net %e% "weights" < 0,
                         "darkorange1", "darkorchid4")
        )

      }

        plt_pos <- ggnet2(
          net = net,edge.alpha = alpha,
          mode = layout,
          node.size = node_outer_size,
          node.color = "black",
          edge.color = "edge_color",
          edge.size = "abs_weights",
          label = TRUE
        ) +
          # geom_point(color = "black",
          #            size = node_outer_size,
          #            alpha = 1) +
          geom_point(color = "white",
                     size = node_inner_size,
                     alpha = 1) +
          geom_text(aes(label = label),
                    color = node_labels_color,
                    size = txt_size)

        # negative
        w_neg <- x$neg_mat * x$pcor_mat
        net <- network::network(x$neg_mat, directed = FALSE)

        # default labels
        if (is.null(node_labels)) {
          network::network.vertex.names(net) <- 1:p
          # custom labels
        } else {
          # check label length
          if (length(node_labels) != p) {
            stop("labels must be of length p (number of nodes)")
          }

          network::network.vertex.names(net) <- node_labels
        }


        network::set.edge.value(
          x = net,
          attrname =  "weights",
          value = w_neg
        )
        network::set.edge.value(
          x = net,
          attrname =  "abs_weights",
          value = abs(w_neg) * edge_multiplier
        )

        if (edge_colors == "classic") {
          network::set.edge.attribute(
            x = net,
            attrname = "edge_color",
            value = ifelse(net %e% "weights" < 0,
                           "brown3", "palegreen3")
          )
        } else if (edge_colors == "color_blind") {
          network::set.edge.attribute(
            x = net,
            attrname = "edge_color",
            value = ifelse(net %e% "weights" < 0,
                           "#009E73", "#D55E00")
          )
        } else if (edge_colors == "vivid") {
          network::set.edge.attribute(
            x = net,
            attrname = "edge_color",
            value = ifelse(net %e% "weights" < 0,
                           "darkorange1", "darkorchid4")
          )

        }



        plt_neg <- ggnet2(
          net = net, edge.alpha = alpha,
          mode = layout, node.size = node_outer_size,
          node.color = "black",
          edge.color = "edge_color",
          edge.size = "abs_weights",
          label = TRUE
        ) +
          # geom_point(color = "black",
          #            size = node_outer_size,
          #            alpha = 1) +
          geom_point(color = "white",
                     size = node_inner_size,
                     alpha = 1) +
          geom_text(aes(label = label),
                    color = node_labels_color,
                    size = txt_size)



        # null
        net <- network::network(x$null_mat, directed = FALSE)

        # default labels
        if (is.null(node_labels)) {
          network::network.vertex.names(net) <- 1:p
          # custom labels
        } else {
          # check label length
          if (length(node_labels) != p) {
            stop("labels must be of length p (number of nodes)")
          }

          network::network.vertex.names(net) <- node_labels
        }





        plt_null <- ggnet2(
          net = net, edge.alpha = alpha,
          mode = layout, node.size = node_outer_size,
          node.color = "white",
          label = TRUE
        ) +
          # geom_point(color = "black",
          #            size = node_outer_size,
          #            alpha = 1) +
          geom_point(color = "white",
                     size = node_inner_size,
                     alpha = 1) +
          geom_text(aes(label = label),
                    color = node_labels_color,
                    size = txt_size)

        plt <- list(plt_pos = plt_pos,
                    plt_neg = plt_neg,
                    plt_null = plt_null)

        } else if(!is.null(node_groups))  {



          if (length(node_groups) != p) {
            stop("labels must be of length p (number of nodes)")
          }



          # positive
          w_pos <- x$pos_mat * x$pcor_mat
          net <- network::network(x$pos_mat, directed = FALSE)

          net %v% "group" <- node_groups


          # add positive with node groups
          network::set.edge.value(
            x = net,
            attrname =  "weights",
            value = w_pos
          )
          network::set.edge.value(
            x = net,
            attrname =  "abs_weights",
            value = abs(w_pos) * edge_multiplier
          )

          if (edge_colors == "classic") {
            network::set.edge.attribute(
              x = net,
              attrname = "edge_color",
              value = ifelse(net %e% "weights" < 0,
                             "brown3", "palegreen3")
            )
          } else if (edge_colors == "color_blind") {
            network::set.edge.attribute(
              x = net,
              attrname = "edge_color",
              value = ifelse(net %e% "weights" < 0,
                             "#009E73", "#D55E00")
            )
          } else if (edge_colors == "vivid") {
            network::set.edge.attribute(
              x = net,
              attrname = "edge_color",
              value = ifelse(net %e% "weights" < 0,
                             "darkorange1", "darkorchid4")
            )

          }



          plt_pos <- ggnet2(
            net = net, edge.alpha = alpha,
            mode = layout, node.size = node_outer_size,
            node.color = "group",
            edge.color = "edge_color",
            edge.size = "abs_weights",
            label = TRUE,...
          ) +
            geom_point(aes(color = color),
                       size = node_outer_size,
                       alpha = 0.5) +
            geom_point(aes(color = color),
                       size = node_inner_size,
                       alpha = 1) +
            geom_text(aes(label = label),
                      color = node_labels_color,
                      size = txt_size)


          # negative
          w_neg <- x$neg_mat * x$pcor_mat
          net <- network::network(x$neg_mat, directed = FALSE)

          net %v% "group" <- node_groups


          # add positive with node groups






          network::set.edge.value(
            x = net,
            attrname =  "weights",
            value = w_neg
          )
          network::set.edge.value(
            x = net,
            attrname =  "abs_weights",
            value = abs(w_neg) * edge_multiplier
          )

          if (edge_colors == "classic") {
            network::set.edge.attribute(
              x = net,
              attrname = "edge_color",
              value = ifelse(net %e% "weights" < 0,
                             "brown3", "palegreen3")
            )
          } else if (edge_colors == "color_blind") {
            network::set.edge.attribute(
              x = net,
              attrname = "edge_color",
              value = ifelse(net %e% "weights" < 0,
                             "#009E73", "#D55E00")
            )
          } else if (edge_colors == "vivid") {
            network::set.edge.attribute(
              x = net,
              attrname = "edge_color",
              value = ifelse(net %e% "weights" < 0,
                             "darkorange1", "darkorchid4")
            )

          }

          plt_neg <- ggnet2(
            net = net, edge.alpha = alpha,
            mode = layout,
            node.size = node_outer_size,
            node.color = "group",
            edge.color = "edge_color",
            edge.size = "abs_weights",
            label = TRUE,...
          ) +
            geom_point(aes(color = color),
                       size = node_outer_size,
                       alpha = 0.5) +
            geom_point(aes(color = color),
                       size = node_inner_size,
                       alpha = 1) +
            geom_text(aes(label = label),
                      color = node_labels_color,
                      size = txt_size)

          # null
          net <- network::network(x$null_mat, directed = FALSE)

          net %v% "group" <- node_groups

          if (edge_colors == "classic") {
            network::set.edge.attribute(
              x = net,
              attrname = "edge_color",
              value = ifelse(net %e% "weights" < 0,
                             "brown3", "palegreen3")
            )
          } else if (edge_colors == "color_blind") {
            network::set.edge.attribute(
              x = net,
              attrname = "edge_color",
              value = ifelse(net %e% "weights" < 0,
                             "#009E73", "#D55E00")
            )
          } else if (edge_colors == "vivid") {
            network::set.edge.attribute(
              x = net,
              attrname = "edge_color",
              value = ifelse(net %e% "weights" < 0,
                             "darkorange1", "darkorchid4")
            )

          }

          plt_null <- ggnet2(
            net = net, edge.alpha = alpha,
            mode = layout, node.size = node_outer_size,
            node.color = "group",
            label = TRUE,...
          ) +
            geom_point(aes(color = color),
                       size = node_outer_size,
                       alpha = 0.5) +
            geom_point(aes(color = color),
                       size = node_inner_size,
                       alpha = 1) +
            geom_text(aes(label = label),
                      color = node_labels_color,
                      size = txt_size)

          plt <- list(plt_pos = plt_pos,
                      plt_neg = plt_neg,
                      plt_null = plt_null)
        }

}
    return(plt)
}
