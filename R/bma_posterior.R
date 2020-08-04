#' @title  Bayesian Model Averaged Posterior Distribution
#'
#' @description Draw samples from the posterior distribution according to the
#'              most probable graphs (weighted by their respective posterior model probabilities)
#'
#' @param object An object of class \code{ggm_search}
#'
#' @param param Character string. What parameter should be computed ? The default is
#'              \code{param = "pcor"} which provides the partial correlations. The other
#'              option is the precision matrix (i.e., \code{precision}).
#'
#' @param iter Number of iterations (posterior samples; defaults to 5000).
#'
#' @param progress Logical. Should a progress bar be included (defaults to \code{TRUE}) ?
#'
#' @return
#' \itemize{
#'
#' \item \code{bma_mean} The mean of the partial correlation or precision matrix
#'                       (\code{p} by \code{p} matrix).
#'
#' \item \code{samples} 3d array of dimensions \code{p} by \code{p} by \code{iter}
#' including the samples partial correlation or precision matrices.
#'
#' }
#'
#'
#' @examples
#' \donttest{
#' # data
#' Y <- ptsd
#'
#' # fit model
#' fit <- ggm_search(Y)
#'
#' # bma
#' bma <- bma_posterior(fit, iter = 100)
#'
#' }
#'
#' @export
bma_posterior <- function(object,
                          param = "pcor",
                          iter = 5000,
                          progress = TRUE){

  if(!is(object, "ggm_search")){
    stop("class not supported. Must be 'ggm_search'")
  }

  n <- object$n

  p <- ncol(object$adj)

  I_p <- diag(p)

  scatter <- object$S * (n - 1)

  approx_marg_ll <- object$approx_marg_ll

  graphs <- object$adj_path[,,which(duplicated(approx_marg_ll) == 0)]

  graphs <- graphs[,,-1]

  graphs_n <- dim(graphs)[3]

  probs <- object$probs

  if(isTRUE(progress)){
    pb <- utils::txtProgressBar(min = 0, max = iter, style = 3)
  }


  samples <- vapply(1:iter, function(s){

    graph_s <- sample(1:graphs_n, 1,replace = FALSE, probs)

    Sigma <- solve(rWishart(1, n + p - 1, solve(scatter + I_p * p))[,,1])

    Theta <- hft_algorithm(Sigma =  Sigma,
                           graphs[, , graph_s],
                           tol = 0.0001,
                           max_iter = 10)$Theta

    if(isTRUE(progress)){
      utils::setTxtProgressBar(pb, s)
    }

    if(param == "pcor"){

      -cov2cor( Theta ) + diag(2, p)

    } else {


      Theta
    }

  }, FUN.VALUE = matrix(0, p, p))

  bma_mean <- apply(samples, 1:2, mean)

  if( is.null( colnames(scatter) ) ){
    colnames(bma_mean) <- 1:p
    row.names(bma_mean ) <- 1:p
  } else {
    colnames(bma_mean) <- colnames(scatter)
    row.names(bma_mean) <- colnames(scatter)
  }

  returned_object <- list(bma_mean = bma_mean,
                          samples = samples)

  class(returned_object) <- c("BGGM",
                              "bma_posterior")

  return(returned_object)

}

print_bma <- function(x,...){
  cat("BGGM: Bayesian Gaussian Graphical Models \n")
  cat("--- \n")
  cat("Bayesian Model Averaged Graph:\n\n")
  print(round(x$bma_mean, 3))
}
