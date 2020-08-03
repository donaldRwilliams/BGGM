#' @title Precision Matrix Posterior Distribution
#'
#' @description Transform the sampled correlation matrices to
#' precision matrices (i.e., inverse covariance matrices).
#'
#' @param object An object of class \code{estimate}.
#'
#' @param progress Logical. Should a progress bar be included (defaults to \code{TRUE}) ?
#'
#' @return A 3d array of dimensions \code{p} by \code{p} by \code{iter} including
#' **unconstrained** (i.e., from th full graph) precision matrices.
#'
#'
#' @export
#'
#' @examples
#' \donttest{
#' # data
#' Y <- ptsd
#'
#' # fit model
#' fit <- estimate(Y)
#'
#' # precisio matrix
#' Theta <- precision(fit)
#'
#' }
#'
precision <- function(object, progress = TRUE){

  if(is(x,"estimate") & is(x,"default")){

    iter <- x$iter

    p <- x$p

    cors <- pcor_to_cor(x)$R

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


  return(precision)

}
