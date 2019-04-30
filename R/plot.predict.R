#' Predict plot
#'
#' @param x object of class "predict"
#' @param size geom_point size
#' @param color geom_point color
#' @param width geom_errorbar width
#'
#' @return ggplot object
#' @export
#'
#' @examples
#'
#' dat <- BGGM::ptsd
#' fit <- estimate(dat, samples)
#'
#' select graph
#' selected <- select(fit, ci_width = 0.95)$adjacency_mat
#'
#' error <- predict(fit,
#'               selected = selected,
#'               measure = "MSE",
#'               test_data  = NULL,
#'               ci_width = 0.95,
#'               samples = 1000)
#'
#' plt <- plot(error, size = 5, color = "red", width = 1)
#'
#' # plot can be changed further with ggplot--e.g., changing theme
#'
#' plt + theme_classic()
plot.predict <- function(x, size = 2, color = "red", width = .1){

# column for nodes
x$summary_error$node <- 1:x$p

# add ordered levels
x$summary_error$node <- factor(x$summary_error$node,
                               levels = order(x$summary_error$post_mean),
                               labels = order(x$summary_error$post_mean))
if(x$measure == "R2"){
  # column for nodes
  x$summary_error$node <- 1:x$p

  # add ordered levels
  x$summary_error$node <- factor(x$summary_error$node,
                                 levels = order(x$summary_error$post_mean),
                                 labels = order(x$summary_error$post_mean))

  plt <- ggplot(x$summary_error, aes(x = node,
                                     y = post_mean)) +
    geom_errorbar(aes(ymin =  x$summary_error[,3],
                      ymax = x$summary_error[,4]), width = width) +
    geom_point(size = size,
               color = color) +
    coord_flip() +
    xlab("Node") +
    ylab("Bayesian R2") +
    theme_classic()
}
if(x$measure == "MSE"){
  # column for nodes
  x$summary_error$node <- 1:x$p

  # add ordered levels
  x$summary_error$node <- factor(x$summary_error$node,
                                 levels = rev(order(x$summary_error$post_mean)),
                                 labels = rev(order(x$summary_error$post_mean)))


  plt <- ggplot(x$summary_error, aes(x = reorder(node, desc(node)),
                                     y = post_mean)) +
    geom_errorbar(aes(ymin = `2.5%`,
                      ymax =`97.5%`), width = width) +
    geom_point(size = size,
               color = color) +
    coord_flip() +
    xlab("Node") +
    ylab("Mean Squared Error")

}
plt
}



