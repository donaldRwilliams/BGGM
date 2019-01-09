#' Title
#'
#' @param x object of class \code{bayes_R2} or \code{predictive error}
#' @param size size of the points
#' @param color color of points
#'
#' @return
#' @export
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
predictive_plot <- function(x, size, color){

  if(class(x) == "bayes_R2"){
    # column for nodes
    x$summary_r2$node <- 1:x$p

    # add ordered levels
    x$summary_r2$node <- factor(x$summary_r2$node,
                                levels = order(x$summary_r2$post_mean),
                                labels = order(x$summary_r2$post_mean))

    plt <- ggplot(x$summary_r2, aes(x = node,
                                    y = post_mean)) +
      geom_point(size = size,
                 color = color) +
      geom_errorbar(aes(ymin = `2.5%`,
                        ymax =`97.5%`), width = .1) +
      coord_flip() +
      ylab("Bayesian R2") +
      xlab("Node")
  }
  plt
}
