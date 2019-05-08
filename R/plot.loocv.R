#' Title
#'
#' @param x
#' @param size
#' @param color
#' @param width
#'
#' @return
#' @export
#'
#' @examples
plot.loocv <- function(x, size = 2, color = "blue", width = 0.01){

  dat <- x$returned_object

  dat$node <- factor(dat$node,

                     levels = rev(dat$node[order(dat$loo)]),

                     labels = rev(dat$node[order(dat$loo)]))
  plt <- ggplot(dat) +
         geom_linerange(aes(x = node,
                            ymin = 0,
                            ymax = loo), color = "grey75") +
         geom_point(aes(x = dat$node,
                        y = loo),
                    size = size,
                    color = color) +
    xlab("Node")

  if(is.numeric(x$samples)){
    plt <- ggplot(dat, aes(x = node,
                           y = loo)) +
      geom_errorbar(aes(ymax = loo + loo_se,
                        ymin = loo - loo_se),
                    width = width) +
      geom_point(aes(x = dat$node,
                     y = loo),
                 size = size,
                 color = color) +
      xlab("Node")
    }

  plt
}
