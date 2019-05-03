#' Title
#'
#' @param x
#' @param ci_width
#' @param size
#' @param width
#'
#' @return
#' @export
#'
#' @examples
plot.estimate <- function(x, ci_width, size, width){
  low <- ((1 - ci_width) / 2)
  up <- 1 - low
  if(!isFALSE(x$analytic)){
    stop("posterior samples required--i.e., analytic = TRUE")
  }
  mat_mu <- mat_low <- mat_up <- mat_names <- matrix(0, nrow = x$p, ncol = x$p)
  parcors <- x$posterior_samples[,  grep("pcors", colnames(x$posterior_samples))]
  mat_mu[] <-   colMeans(parcors)
  mat_low[] <- apply(parcors, 2, quantile, low)
  mat_up[]  <-  apply(parcors, 2, quantile, up)

  mat_names[] <-  unlist(lapply(1:x$p, function(y) paste(y, 1:x$p, sep = "--")))

  res  <- data.frame(edge = mat_names[upper.tri(mat_names)],
                     mu = mat_mu[upper.tri(mat_mu)],
                     low = mat_low[upper.tri(mat_low)],
                     up = mat_up[upper.tri(mat_mu)])


  res$sig <- ifelse(res$low< 0 & res$up > 0, 0, 1)
  res$edge2 <- factor(res$edge, levels = res$edge[order(res$mu)], labels = res$edge[order(res$mu)])

  plt <- ggplot(res, aes(y = mu, x = edge2))  +

    geom_errorbar(aes(ymax = up, ymin = low), color = "grey75", width = width) +
    geom_point(aes(color = as.factor(sig)),show.legend = F, size = size) +
    theme_bw() +
    geom_hline( yintercept = 0, linetype = "dotted") +
    theme(panel.grid = element_blank()) +
    coord_flip() +
    ylab("Edge (Partial Correlation)") +
    xlab("")

  return(plt)



}
