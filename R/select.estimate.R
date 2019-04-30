#' Title
#'
#' @param x object of class \code{bayes_estimate}
#' @param ci_width credible interval width used for the decision rule
#'
#' @return the selected partial correlation matrix (the conditional dependence structure) and
#'  the adjacency matrix is returned
#' @export
#'
#' @examples
select.estimate  <- function(x, ci_width){

  # check object class
  if(class(x) !=  "estimate"){
    stop("Must be an object class bayes_estimate")
  }

  # ensure ci_width is allowed
  if(ci_width >= 1 | ci_width <= 0){
    stop("ci_width must be between 0 and 1")
  }

  # matrices for selected edges and the adjacency matrix
  adjacency_mat <- matrix(0, x$p, x$p)

  pcor_samples <- x$posterior_samples[,  grep("pcors", colnames(x$posterior_samples))]

  adjacency_mat[] <- apply(pcor_samples, 2, ci_helper, ci_width)

  returned_object <- list(pcor_selected = x$parcors_mat * adjacency_mat,
                          adjacency_mat = adjacency_mat)

  class(returned_object) <- "estimate"
  returned_object
}
