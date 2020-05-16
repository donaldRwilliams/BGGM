#' @title Extract the Partial Correlation Matrix
#'
#' @param object A model estimated with \strong{BGGM}. All classes are supported, assuming
#' there is matrix to extract.
#'
#' @param ... Currently ignored.
#'
#' @return The estimated partial correlation matrix.
#'
#' @examples
#' \donttest{
#' # note: iter = 250 for demonstrative purposes
#'
#' # data
#' Y <- ptsd[,1:5] + 1
#'
#' # ordinal
#' fit <- estimate(Y, type = "ordinal",
#'                 iter = 250)
#'
#' pcor_mat(fit)
#' }
#' @export
pcor_mat <- function(object, ...){

  cn <- colnames(object$Y)

  colnames(object$pcor_mat) <- cn

  rownames(object$pcor_mat) <- cn

  if(is.null(cn)) {
    cn <- 1:object$p
  }

  round(object$pcor_mat, 3)

}
