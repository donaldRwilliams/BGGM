#' @title Constrained Posterior Distribution
#'
#' @description Compute the posterior distribution
#'              with off-diagonal elements of the precision matrix constrained
#'              to zero.
#'
#' @param object An object of class \code{estimate} or \code{explore}
#'
#' @param adj A \code{p} by \code{p} adjacency matrix. The zero entries denote the
#'            elements that should be contrained to zero.
#'
#' @param iter Number of iterations (posterior samples; defaults to 5000).
#'
#'
#' @param progress Logical. Should a progress bar be included (defaults to \code{TRUE}) ?
#'
#' @return An object of class \code{contrained}, including
#'
#' \itemize{
#'
#' \item \code{precision_mean} The posterior mean for the precision matrix.
#'
#' \item \code{pcor_mean} The posterior mean for the precision matrix.
#'
#' \item \code{precision_samps} A 3d array of dimension \code{p} by \code{p} by \code{iter}
#'                              including the sampled precision matrices.
#'
#'  \item \code{pcor_samps} A 3d array of dimension \code{p} by \code{p} by \code{iter}
#'                              including sampled partial correlations matrices.
#' }
#'
#' @note
#' The \code{object} includes posterior samples from the full graph. This function
#' then corrects those samples to reflect the estimated structure of zeros. Note
#' that this should be used with caution, as this is essentially treating the graph
#' as known (but it is not).
#'
#' @examples
#' \donttest{
#'
#' # data
#' Y <- bfi[,1:10]
#'
#' # sample posterior
#' fit <- estimate(Y, iter = 100)
#'
#' # select graph
#' sel <- select(fit)
#'
#' # constrained posterior
#' post <- constrained_posterior(object = fit,
#'                               adj = sel$adj,
#'                               iter = 100,
#'                               progress = FALSE)
#'
#' }
#' @export
constrained_posterior <- function(object,
                                  adj,
                                  iter = 1000,
                                  progress = TRUE){

  if (!any(class(object) %in% c("estimate", "explore"))) {
    stop("object must be of class 'estimate' or 'explore'")
  }

  if (object$iter < iter) {
    stop("iter exceeds iter in the object")
  }

  # ensure diagonal is 1
  diag(adj) <- 1

  # number of nodes
  p <- object$p

  # covert to correlations
  cors <- pcor_to_cor(object, iter = iter)

  if(progress){
    pb <- txtProgressBar(min = 0, max = iter, style = 3)

  }

  # precision matrix
  Theta_samps <- vapply(1:iter, function(s) {

    Theta <- BGGM:::hft_algorithm(Sigma = cors$R[,,s],
                                  adj = adj, tol = 0.00001, max_iter = 10)$Theta

    if(progress){
      utils::setTxtProgressBar(pb, s)
    }

    Theta

  }, FUN.VALUE =  matrix(0,p,p)

  )

  # partial correlations
  pcor_samps <- vapply(1:iter, function(x)
    -cov2cor(Theta_samps[,,x]) + diag(2, p),
    FUN.VALUE =  matrix(0, p, p)
  )

  Theta_mu <- apply(Theta_samps, 1:2, mean)
  pcor_mu <- apply(pcor_samps, 1:2, mean)

  colnames(pcor_mu) <- colnames( object$Y )
  row.names(pcor_mu) <- colnames( object$Y )
  returned_object <- list(
    precision_mean = Theta_mu,
    pcor_mean = pcor_mu,
    precision_samps = Theta_samps,
    pcor_samps = pcor_samps
  )

  class(returned_object) <- c("BGGM",
                              "constrained")

  return(returned_object)

}


print_constrained <- function(x, ...){
  cat("BGGM: Bayesian Gaussian Graphical Models\n")
  cat("Constrained posterior\n")
  cat("---\n")
  cat("Estimates: \n\n")
  print(round(x$pcor_mean, 3))
}
