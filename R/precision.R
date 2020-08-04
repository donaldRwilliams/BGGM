#' @title Precision Matrix Posterior Distribution
#'
#' @description Transform the sampled correlation matrices to
#' precision matrices (i.e., inverse covariance matrices).
#'
#' @param object An object of class \code{estimate}.
#'
#' @param progress Logical. Should a progress bar be included (defaults to \code{TRUE}) ?
#'
#' @note The estimated precision matrix is the inverse of the \strong{correlation} matrix.
#'
#' @return
#'
#' \itemize{
#'
#' \item \code{precision_mean} The mean of the precision matrix (\code{p} by \code{p} matrix).
#'
#' \item \code{precision} 3d array of dimensions \code{p} by \code{p} by \code{iter}
#' including \strong{unconstrained} (i.e., from th full graph)
#' precision matrices.
#'
#'
#' }
#' @examples
#' \donttest{
#' # data
#' Y <- ptsd
#'
#' # fit model
#' fit <- estimate(Y)
#'
#' # precision matrix
#' Theta <- precision(fit)
#'
#' }
#'
#' @export
precision <- function(object,
                      progress = TRUE){

  if(is(object,"estimate") & is(object,"default")){

    iter <- object$iter

    p <- object$p

    cors <- pcor_to_cor(object)$R

    if(isTRUE(progress)){
      pb <- utils::txtProgressBar(min = 0, max = iter, style = 3)
    }

  precision <- vapply(1:iter, function(s){

    Theta <- solve(cors[,,s])

    if(isTRUE(progress)){
      utils::setTxtProgressBar(pb, s)
    }
    Theta

    }, FUN.VALUE = matrix(0, 20, 20))


  } else {

    stop("class not currently supported")

    }

  precision_mean = apply(precision, 1:2, mean)

  returned_object <- list(precision_mean = precision_mean,
                          precision = precision)

  class(returned_object) <- c("BGGM",
                              "precision")

  return(returned_object)

}


print_precision <- function(x,...) {
  mat <- x$precision_mean
  p <- ncol(mat)
  colnames(mat) <- 1:p
  row.names(mat) <- 1:p
  cat("BGGM: Bayesian Gaussian Graphical Models \n")
  cat("--- \n")
  cat("Estimate:\n\n")
  print(round(mat, 3))
  }
