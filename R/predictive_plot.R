#' predictive plot
#'
#' @description plots measures of predictive accuracy (in and out of sample)
#'
#' @param x object of class \code{bayes_R2} or \code{predictive error}
#' @param size size of the points
#' @param color color of points
#'
#'
#' @return
#' @export
#'
#'@note The error bars correspond to what was specificed in the function \code{bayes_R2}
#'
#' @examples
#'
#' # X is the data matrix of dimensions n by p
#' fit <- bayes_estimate(X)
#'
#' # select the graphical structure
#' select <- estimate_select(fit, ci_width = .95)
#'
#' # compute Bayesian R2
#' r2 <-  bayes_R2(fit, selected = select$adjacency_mat, ci_width = .95, samples = 1000)
#'
#' plt <- predictive_plot(x = r2, size = 3, color = color)
#'
#' # plot can be changed further with ggplot--e.g., changing theme
#'
#' plt + theme_classic()
#'
#'
#'
predictive_plot <- function(x, size, color){


    # column for nodes
    x$summary_error$node <- 1:x$p

    # add ordered levels
    x$summary_error$node <- factor(x$summary_error$node,
                                levels = order(x$summary_error$post_mean),
                                labels = order(x$summary_error$post_mean))





    if(x$measure == "bayes_R2"){
      # column for nodes
      x$summary_error$node <- 1:x$p

      # add ordered levels
      x$summary_error$node <- factor(x$summary_error$node,
                                     levels = order(x$summary_error$post_mean),
                                     labels = order(x$summary_error$post_mean))

      plt <- ggplot(x$summary_error, aes(x = node,
                                           y = post_mean)) +
        geom_errorbar(aes(ymin = `2.5%`,
                          ymax =`97.5%`), width = .1) +
        geom_point(size = size,
                   color = color) +
        coord_flip() +
        xlab("Node") +
        ylab("Bayesian R2")
      }
    if(x$measure == "kl"){
      # column for nodes
      x$summary_error$node <- 1:x$p

      # add ordered levels
      x$summary_error$node <- factor(x$summary_error$node,
                                     levels = order(x$summary_error$post_mean),
                                     labels = order(x$summary_error$post_mean))
      plt <- ggplot(x$summary_error, aes(x = reorder(node, desc(node)),
                                         y = post_mean)) +
        geom_errorbar(aes(ymin = `2.5%`,
                          ymax =`97.5%`), width = .1) +
        geom_point(size = size,
                   color = color) +
        coord_flip() +
        xlab("Node") +
        ylab("KL-Divergence")

    }

    if(x$measure == "mse"){
      # column for nodes
      x$summary_error$node <- 1:x$p

      # add ordered levels
      x$summary_error$node <- factor(x$summary_error$node,
                                     levels = rev(order(x$summary_error$post_mean)),
                                     labels = rev(order(x$summary_error$post_mean)))


      plt <- ggplot(x$summary_error, aes(x = reorder(node, desc(node)),
                                         y = post_mean)) +
        geom_errorbar(aes(ymin = `2.5%`,
                          ymax =`97.5%`), width = .1) +
        geom_point(size = size,
                   color = color) +
        coord_flip() +
        xlab("Node") +
        ylab("Mean Squared Error")

    }

    if(x$measure == "mae"){

      # column for nodes
      x$summary_error$node <- 1:x$p

      # add ordered levels
      x$summary_error$node <- factor(x$summary_error$node,
                                     levels = rev(order(x$summary_error$post_mean)),
                                     labels = rev(order(x$summary_error$post_mean)))


      plt <- ggplot(x$summary_error, aes(x = node,
                                         y = post_mean)) +
        geom_errorbar(aes(ymin = `2.5%`,
                          ymax =`97.5%`), width = .1) +
        geom_point(size = size,
                   color = color) +
        coord_flip() +
        xlab("Node") +
        ylab("Mean Absolute Error")
    }
   plt
}

