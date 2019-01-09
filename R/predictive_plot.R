#' Title
#'
#' @param x
#' @param size
#' @param color
#'
#' @return
#' @export
#'
#' @examples
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
