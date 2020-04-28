#' GGM Selection: Estimation
#' @name select.estimate
#'
#' @description Select the edge set for \code{estimate} objects. The graph is determined with
#' the posterior distribution, in particular credible interval exclusion of zero. This
#' corresponds to a directional posterior probability \insertCite{Marsman2017a}{BGGM}. For example,
#' the probability (conditional on the model) of a positive edge is at least 97.5% when the 95%
#' credible interval excludes zero.
#'
#' @param object object of class \code{estimate.default}
#' @param cred credible interval width used for the decision rule
#'
#' @param alternative a character string specifying the alternative hypothesis,
#'                    must be one of "two.sided", "greater" (default) or "less".
#'                    See note for futher details.

#' @param ... not currently used
#'
#' @references
#' \insertAllCited{}
#'
#' @return An object of class \code{select.estimate}:
#'
#' \itemize{
#' \item \code{pcor_adj} selected partial correlation matrix
#' \item \code{adj} adjacency matrix for the selected edges
#' \item \code{ci} credible interval width
#' \item \code{object} object of class \code{estimate}
#' }
#'

#'
#'
#' @note
#'
#' \strong{Alternative Default:}
#'
#' This package was built for the social-behavioral sciences in particular. In these applications, there is
#' strong theory that expects \emp{all} effects to be positive. This is known as a "positive manifold" and
#' this notion has a rich tradition in psychometrics. Hence, by default this expecation is included into graph selection
#' (\code{alternative = "greater"}). This results in the estimted structure including only positive edges. Further
#' details can be found at the blog "Dealing with Negative (Red) Edges in Psychological Networks: Frequentist Edition"
#' (\href{https://donaldrwilliams.github.io/2020/03/29/dealing-with-negative-red-edges-in-psychological-networks-frequentist-edition/}{link})
#'
#' @examples
#'
#' # Analytic = TRUE
#'# p = 5
#' Y <- BGGM::bfi[,1:5]
#'
#' # analytic solution
#' fit_analytic <- estimate(Y, analytic = TRUE)
#'
#' # select E
#' E <- select(fit_analytic, ci_width = 0.95)
#'
#' # non-zero partial correlations
#' E$partials_non_zero
#'
#' # adjacency matrix
#' E$adjacency_non_zero
#' @export
select.estimate <- function(object,
                            cred = 0.95,
                            alternative = "greater"){

  if(isFALSE(object$analytic)){

    pcors <- object$post_samp$pcors[,,51:(fit$iter +50)]

    if(alternative == "two.sided"){

      lb <- (1 - cred) / 2
      ub <- 1 - lb

      adj <- ifelse(apply(pcors, 1:2, quantile, lb) < 0 &
                  apply(pcors, 1:2, quantile, ub) > 0, 0, 1)

      } else if(alternative == "greater") {

        lb <- (1 - cred)
        adj <- ifelse(apply(pcors, 1:2, quantile, lb) > 0, 1, 0)

        } else {

          ub <- cred
          adj <- ifelse(apply(pcors, 1:2, quantile, ub) < 0, 1, 0)

          }

    # analytic
    } else {

      if(alternative == "two.sided"){

        lb <- (1 - cred) / 2
        ub <- 1 - lb

        z_stat <- abs(object$analytic_fit$inv_map /  sqrt(object$analytic_fit$inv_var))

        adj <- ifelse(z_stat >  qnorm(ub), 1, 0)

      } else if (alternative == "greater") {

        ub <- 1 - cred

        z_stat <- (-object$analytic_fit$inv_map) /  sqrt(object$analytic_fit$inv_var)

        adj <- ifelse( z_stat > qnorm(ub, lower.tail = FALSE), 1, 0)



      } else if(alternative == "less"){


        ub <- 1 - cred

        z_stat <- (object$analytic_fit$inv_map) /  sqrt(object$analytic_fit$inv_var)

        adj <- ifelse(z_stat > qnorm(ub, lower.tail = FALSE), 1, 0)

      }





}

  pcor_adj <- adj * object$pcor_mat

  returned_object <- list(pcor_adj = pcor_adj,
                           adj = adj,
                           alternative = alternative,
                           cred = cred,
                           object = object)

class(returned_object) <- c("BGGM",
                             "select.estimate",
                             "estimate")
  returned_object
}

#' @title S3 select method
#' @name select
#' @description S3 select method
#' @param object object of class \code{estimate}, \code{explore}, or ..
#' @param ... not currently used
#' @return \code{select} works with the following methods:
#' \itemize{
#' \item \code{\link{select.estimate}}
#' \item \code{\link{select.explore}}
#' \item \code{\link{select.ggm_compare_estimate}}
#' }
#' @export
select <- function(object,...){
  UseMethod("select", object)
}


print_select_estimate <- function(x, ...){
  object <- x
  p <- ncol(object$pcor_adj)
  cat("BGGM: Bayesian Gaussian Graphical Models \n")
  cat("--- \n")
  cat("Type:", object$object$type, "\n")
  cat("Analytic:", object$object$analytic, "\n")
  cat("Formula:", paste(as.character(fit$formula), collapse = " "), "\n")
  cat("Posterior Samples:", object$object$iter, "\n")
  cat("Credible Interval:",  gsub("*0.","", formatC( round(object$cred, 4), format='f', digits=2)), "% \n")
  cat("--- \n")
  cat("Call: \n")
  print(object$object$call)
  cat("--- \n")
  cat("Selected:\n\n")
  mat <- object$pcor_adj
  colnames(mat) <- 1:p
  row.names(mat) <- 1:p
  print(round(mat, 3))
  cat("--- \n")
}






















#' Plot \code{select.estimate} Network
#'
#' @param x object of class \code{select.estimate}
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
#'
#' @examples
#' Y <- BGGM::bfi[, 1:20]
#'
#' # analytic approach (sample by setting analytic = FALSE)
#' fit_analytic <- estimate(Y, analytic = TRUE)
#'
#' # select the graph (edge set E)
#' E <- select(fit_analytic, ci_width = 0.95)
#'
#' # plot
#' plt <- plot(E,
#'             node_labels = letters[1:20],
#'             node_labels_color = "white",
#'             node_groups = rep(c("1", "2", "3", "4"), each = 5),
#'             edge_colors = "classic",
#'             alpha = 0.5, palette = "Set2")


plot.select.estimate <-
  function(x,
           layout = "circle",
           edge_colors = "classic",
           node_labels = NULL,
           node_labels_color = "black",
           node_groups = NULL,
           node_outer_size = 12,
           node_inner_size = 11,
           alpha = 0.50,
           txt_size = 8,
           edge_multiplier = 1,
           ...) {
    label <- NULL
    color <- NULL
    # number of nodes
    p <- ncol(x$adjacency_non_zero)

    # network
    net <- network::network(x$adjacency_non_zero)


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
        geom_point(aes(color = color),
                   size = node_inner_size,
                   alpha = 1) +
        geom_text(aes(label = label),
                  color = node_labels_color,
                  size = txt_size)

      plt <- list(plt = plt)
    }

    if (!is.null(x$rope)) {
      # network
      net <- network::network(x$adjacency_zero, directed = FALSE)

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
    return(plt)
  }


