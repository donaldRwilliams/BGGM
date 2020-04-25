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
                           alternative = "two.sided",
                           ...){
  # rename
  x <- object

  # hyp probability
  hyp_prob <- BF_cut / (BF_cut + 1)


  if(x$type == "continuous"){

    # posterior samples
    post_samp <- x$post_samp

    # prior samples
    prior_samp <- x$prior_samp

    } else if (x$type == "binary" |
               x$type == "ordinal" |
               x$type == "mixed"){


    posterior_samples <- x$samples$fisher_z_post

    prior_samples <- x$samples$fisher_z_prior

    }

  # two sided testing
  if(alternative == "two.sided"){

    # posterior
    post_sd <- apply(post_samp$fisher_z[,,(51:x$iter)], 1:2, sd)
    post_mean  <- x$pcor_mat
    post_dens <- dnorm(0, post_mean, post_sd )

    # prior
    prior_sd <- apply(prior_samp$fisher_z[,,(51:x$iter)], 1:2, sd)
    prior_dens <- dnorm(0, 0, mean(prior_sd))

    # BF
    BF_10_mat <- prior_dens / post_dens
    BF_01_mat <- 1 / BF_10_mat
    diag(BF_01_mat) <- 0
    diag(BF_10_mat) <- 0

    # selected: alternative
    Adj_10 <- ifelse(BF_10_mat > BF_cut, 1, 0)

    # selected: null
    Adj_01 <- ifelse(BF_10_mat < 1 / BF_cut, 1, 0)
    diag(Adj_01) <- 0
    diag(Adj_10) <- 0

    # returned object
    returned_object = list(pcor_mat_zero = post_mean * Adj_10,
                           pcor_mat = round(post_mean, 3),
                           pcor_sd = round(post_sd, 3),
                           Adj_10 = Adj_10,
                           Adj_01 = Adj_01,
                           BF_10 = BF_10_mat,
                           BF_01 = BF_01_mat,
                           BF_cut = BF_cut,
                           alternative = alternative,
                           call = match.call(),
                           type = x$type)

    # one sided greater
    } else if(alternative == "greater"){

      # posterior
      post_sd <- apply(post_samp$fisher_z[,,(51:x$iter)], 1:2, sd)
      post_mean  <- x$pcor_mat
      post_dens <- dnorm(0, post_mean, post_sd )

      # prior
      prior_sd <- apply(prior_samp$fisher_z[,,(51:x$iter)], 1:2, sd)
      prior_dens <- dnorm(0, 0, mean(prior_sd))

      # BF (two sided)
      BF_10_mat <- prior_dens / post_dens

      # BF one sided
      BF_20_mat <-  BF_10_mat * ((1 - pnorm(0, post_mean, post_sd)) * 2)

      # BF null
      BF_02_mat <- 1 / BF_20_mat
      diag(BF_02_mat) <- 0
      diag(BF_20_mat) <- 0

      # selected edges (alternative)
      Adj_20 <- ifelse(BF_20_mat > BF_cut, 1, 0)

      # selected edges (null)
      Adj_02 <- ifelse(BF_02_mat > BF_cut, 1, 0)
      diag(Adj_02) <- 0
      diag(Adj_20) <- 0

      # returned object
      returned_object = list(
        pcor_mat_zero = post_mean * Adj_20,
        pcor_mat = round(post_mean, 3),
        pcor_sd = round(post_sd, 3),
        Adj_20 = Adj_20,
        Adj_02 = Adj_02,
        BF_20 = BF_20_mat,
        BF_02 = BF_02_mat,
        BF_cut = BF_cut,
        alternative = alternative,
        call = match.call(),
        type = x$type
      )

    # one side less
    } else if(alternative == "less")  {

      # posterior
      post_sd <- apply(post_samp$fisher_z[,,(51:x$iter)], 1:2, sd)
      post_mean  <- x$pcor_mat
      post_dens <- dnorm(0, post_mean, post_sd )

      # prior
      prior_sd <- apply(prior_samp$fisher_z[,,(51:x$iter)], 1:2, sd)
      prior_dens <- dnorm(0, 0, mean(prior_sd))

      # BF (two sided)
      BF_10_mat <- prior_dens / post_dens

      # BF one sided
      BF_20_mat <- BF_10_mat * (pnorm(0, post_mean, post_sd) * 2)

      # BF null
      BF_02_mat <- 1 / BF_20_mat
      diag(BF_02_mat) <- 0
      diag(BF_20_mat) <- 0

      # selected edges (alternative)
      Adj_20 <- ifelse(BF_20_mat > BF_cut, 1, 0)

      # selected edges (null)
      Adj_02 <- ifelse(BF_02_mat > BF_cut, 1, 0)
      diag(Adj_02) <- 0
      diag(Adj_20) <- 0

      # returned object
      returned_object = list(
        pcor_mat_zero = post_mean * Adj_20,
        pcor_mat = round(post_mean, 3),
        pcor_sd = round(post_sd, 3),
        Adj_20 = Adj_20,
        Adj_02 = Adj_02,
        BF_20 = BF_20_mat,
        BF_02 = BF_02_mat,
        BF_cut = BF_cut,
        alternative = alternative,
        call = match.call(),
        type = x$type
      )

      # exhaustive testing
      } else if (alternative == "exhaustive")

        if(alternative == "exhaustive"){

          if(is.null(hyp_prob)){

            stop("posterior probability must be specificed \n for exhaustive hypothesis testing")

            }


          mat_names <-  matrix(0, x$p, x$p)

          # matrix names
          mat_names[] <-  unlist(lapply(1:x$p,
                                        FUN = function(z) paste(1:x$p, z, sep = "--")))

          # posterior
          post_sd <- apply(post_samp$fisher_z[,,(51:x$iter)], 1:2, sd)
          post_mean  <- x$pcor_mat
          post_dens <- dnorm(0, post_mean, post_sd)

          # prior
          prior_sd <- apply(prior_samp$fisher_z[,,(51:x$iter)], 1:2, sd)
          prior_dens <- dnorm(0, 0, mean(prior_sd))

          # BF (two sided)
          BF_10_mat <- prior_dens / post_dens

          # BF less
          BF_less <- BF_10_mat  * (pnorm(0, post_mean, post_sd) * 2)
          BF_greater <-  BF_10_mat * ((1 - pnorm(0, post_mean, post_sd)) * 2)

          # BF null
          BF_null <- 1 / BF_10_mat

          prob_null <-  BF_null / (BF_null + BF_greater + BF_less)
          prob_greater <-  BF_greater / (BF_null + BF_greater + BF_less)
          prob_less <-  BF_less / (BF_null + BF_greater + BF_less)

          prob_mat <-  prob_null + prob_greater + prob_less

          prob_dat = data.frame(edge = mat_names[upper.tri(mat_names)],
                                prob_zero = prob_null[upper.tri(prob_null)],
                                prob_greater = prob_greater[upper.tri(prob_greater)],
                                prob_less = prob_less[upper.tri(prob_less)])

          # no rownames
          row.names(prob_dat) <- c()

          # selected  (null)
          null_mat <- ifelse(prob_null > hyp_prob, 1, 0)

          # selected (positive)
          pos_mat <-  ifelse(prob_greater > hyp_prob, 1, 0)


          # selected  (negative)
          neg_mat <-  ifelse(prob_less > hyp_prob, 1, 0)

          # negative edges
          returned_object <- list(
            post_prob = prob_dat,
            neg_mat = neg_mat,
            pos_mat = pos_mat,
            null_mat = null_mat,
            alternative = alternative,
            pcor_mat = round(post_mean, 3),
            pcor_sd = round(post_sd, 3),
            call = match.call(),
            prob = hyp_prob,
            type = x$type
          )
        } else {

          stop("alternative not supported. see documentation")
    }

  class(returned_object) <- c("BGGM",
                              "select.explore",
                              "explore")
  returned_object
}




#' Summary Method for \code{select.explore} Objects
#'
#' @param object object of class \code{select.explore}.
#'
#' @return a data frame including the posterior mean, standard deviation,
#' and posterior hypothesis probabilities for each relation.
#' @export
summary.select.explore <- function(object){

  x <- object
  p <- ncol(x$pcor_mat)
  mat_names <- matrix(0, p, p)

  mat_names[] <-  unlist(lapply(1:p,
                                FUN = function(z) paste(1:p, z, sep = "--")))

  if(x$alternative == "two.sided"){

    post_mean <- x$pcor_mat[upper.tri(x$pcor_mat)]
    post_sd <-  x$pcor_sd[upper.tri(x$pcor_sd)]
    prob_H1 <- x$BF_10[upper.tri(x$BF_10)] / (x$BF_10[upper.tri(x$BF_10)] + 1)
    prob_H0 <- 1 - prob_H1
    summ <-  data.frame(
      Relation = mat_names[upper.tri(mat_names)],
      Post.mean = post_mean,
      Post.sd = post_sd,
      Pr.H0 = round(prob_H0, 3),
      Pr.H1 = round(prob_H1, 3)
    )

  } else if (x$alternative == "greater"){

    post_mean <- x$pcor_mat[upper.tri(x$pcor_mat)]
    post_sd <-  x$pcor_sd[upper.tri(x$pcor_sd)]
    prob_H1 <- x$BF_20[upper.tri(x$BF_20)] / (x$BF_20[upper.tri(x$BF_20)] + 1)
    prob_H0 <- 1 - prob_H1
    summ <-  data.frame(
      Relation = mat_names[upper.tri(mat_names)],
      Post.mean = post_mean,
      Post.sd = post_sd,
      Pr.H0 = round(prob_H0, 3),
      Pr.H1 = round(prob_H1, 3)
    )



  } else if (x$alternative == "less" | x$alternative == "greater"){

    post_mean <- x$pcor_mat[upper.tri(x$pcor_mat)]
    post_sd <-  x$pcor_sd[upper.tri(x$pcor_sd)]
    prob_H1 <- x$BF_20[upper.tri(x$BF_20)] / (x$BF_20[upper.tri(x$BF_20)] + 1)
    prob_H0 <- 1 - prob_H1
    summ <-  data.frame(
      Relation = mat_names[upper.tri(mat_names)],
      Post.mean = post_mean,
      Post.sd = post_sd,
      Pr.H0 = round(prob_H0, 3),
      Pr.H1 = round(prob_H1, 3)
    )



  } else {

    summ <- cbind.data.frame( x$post_prob[,1],
                              x$pcor_mat[upper.tri(x$pcor_mat)],
                              x$pcor_sd[upper.tri(x$pcor_sd)],
                              round(x$post_prob[,2:4], 3))

    colnames(summ) <- c("Relation",
                        "Post.mean",
                        "Post.sd",
                        "Pr.H0",
                        "Pr.H1",
                        "Pr.H2")


  }

  returned_object <- list(summary = summ, object = object)

  class(returned_object) <- c("BGGM",
                              "explore", "select.explore",
                              "summary")
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
