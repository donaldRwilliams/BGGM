#' Partial Correlations to Correlations
#'
#' @name pcor_to_cor
#'
#' @description Convert the partial correlation matrices into correlation matrices.
#'
#' @param object an object of class \code{estimate} or \code{explore}
#'
#' @param iter numeric. How many iterations (i.e., posterior samples) should be used?
#'             The default uses all of the samples, but note that this can take a long
#'             time with large matrices.
#'
#' @return
#'
#' \itemize{
#'
#' \item \code{R} An array including the correlation matrices
#'               (of dimensions \emph{p} by \emph{p} by \emph{iter})
#'
#' \item \code{R_mean} Posterior mean of the correlations (of dimensions \emph{p} by \emph{p})
#' }
#'
#' @note
#' The 'default' prior distributions are specified for partial correlations in particular. This
#' means that the implied prior distribution will not be the same for the correlations.
#'

#' @export
pcor_to_cor <- function(object, iter = NULL){

  if(!is(object, "default")){

    stop("Class not supported. Must but an 'estimate' or 'explore' object.")

  }

  post_samps <- -object$post_samp$pcors

  dims <- dim(post_samps)

  if(!is.null(iter)){

    if((dims[3] - 50) < iter){

    warning("Iterations do not exist (too large). Using all iterations in the object.")

    iter <- dims[3] - 50
  }

  } else {

    iter <- dims[3] - 50

}

  p <- dims[1]

  object <- post_samps[,,-c(1:50)]

  object <- post_samps[,,1:iter]

  # call c ++ for speed
  returned_object <- .Call("_BGGM_pcor_to_cor_internal",
                           PACKAGE = "BGGM",
                           x = object, p = p)

  returned_object

}
