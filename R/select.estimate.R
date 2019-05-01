#' Select graphical structure
#'
#' @param x object of class \code{bayes_estimate}
#' @param ci_width credible interval width used for the decision rule
#' @param rope
#' @param prob
#'
#' @return partial  selected partial correlations
#' @return adjacency selected adjacency matrix
#' @export
#'
#' @examples
select.estimate <- function(x, ci_width = 0.95, rope = NULL, prob = 0.975){

  # check object class
  if(class(x) !=  "estimate"){
    stop("Must be an object class bayes_estimate")
  }

  # ensure ci_width is allowed
  if(ci_width >= 1 | ci_width <= 0){
    stop("ci_width must be between 0 and 1")
  }
  pcor_samples <- x$posterior_samples[,  grep("pcors", colnames(x$posterior_samples))]
  if(is.null(rope)){
    # matrices for selected edges and the adjacency matrix
    adjacency_mat <- matrix(0, x$p, x$p)



    adjacency_mat[] <- apply(pcor_samples, 2, BGGM:::ci_helper, ci_width)

    returned_object <- list(partials = x$parcors_mat * adjacency_mat,
                            adjacency = adjacency_mat, call = match.call(), ci = ci_width, rope = rope)
  }

  if(is.numeric(rope)){
    if(!is.numeric(prob)){
      stop("prob must be specificed (0 - 1) when rope = TRUE")
    }
    in_rope <- apply(pcor_samples, 2, BGGM:::rope_helper,  .1)
    out_rope <- 1 - in_rope
    nonzero_mat <- zero_mat <- matrix(0, x$p, x$p)
    nonzero_mat[] <- ifelse(out_rope >= prob, 1, 0)
    zero_mat[] <- ifelse(in_rope >= prob, 1, 0)

    returned_object <- list(partials_non_zero = x$parcors_mat * nonzero_mat,
                            adjacency_non_zero = nonzero_mat,
                            partial_zero = x$parcors_mat * zero_mat,
                            adjacency_zero = zero_mat,
                            call = match.call(), rope = rope, prob = prob)
  }


  class(returned_object) <- "select.estimate"
  returned_object$call <- match.call()
  returned_object
}