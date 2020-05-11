#' Compute Custom Network Statistics
#'
#' This function allows for computing custom network statistics for
#' weighted adjacency matrices (partial correlations). The statistics are computed for
#' each of the sampled matrices, resulting in a distribution.
#'
#' @param object An object of class \code{estimate}.
#'
#' @param FUN A custom function for computing the statistic.
#'
#' @param select Logical. Should the graph be selected ? The default is currently \code{FALSE}.
#'
#' @param cred Numeric. Credible interval between 0 and 1  (default is 0.95) that is used for selecting the graph.
#'
#' @param extract A character specifying what, if anything, should be extracted from the function. For example,
#'                to extract the object "r", set \code{extract = "r"} and internally \code{$r} will be added
#'                to the end of the function.
#'
#' @param ... Currently ignored
#'
#' @return An object defined by \code{FUN}.
#'
#' @details
#'
#' The user has complete control of this function. Hence, care must be taken as to what \code{FUN} returns
#' and in what format.
#'
#' @export
#'
#' @examples
#' \donttest{
#'
#' ####################################
#' ###### example 1: assortment #######
#' ####################################
#' library(assortnet)
#'
#' Y <- BGGM::bfi[,1:10]
#' membership <- c(rep("a", 5), rep("c", 5))
#'
#' # fit model
#' fit <- estimate(Y = Y,
#'                 analytic = FALSE,
#'                 iter = 1000)
#'
#' # list of columns belowinging in each group
#' e.g., first 5 are "a", last 5 are "c"
#'
#' membership <- c(rep("a", 5), rep("c", 5))
#'
#'
#' net_stat <- roll_your_own(object = fit,
#'                           FUN = assortment.discrete,
#'                           types = membership, weighted = TRUE,
#'                           SE = FALSE, M = 1,
#'                           extract = "r")
#'
#' hist(net_stat)
#' }
#'
#'
roll_your_own <- function(object,
                          FUN,
                          select = FALSE,
                          cred = 0.95,
                          extract = NULL, ...){

  if(! all( c("estimate", "default") %in% class(fit)) ){
    stop("class must be 'estimate'")
  }

  if(!isFALSE(select)){

    sel <- select(fit, cred = cred)
    adj <- sel$adj

  } else {

    p <- ncol(object$pcor_mat)
    adj <- matrix(1, p, p)

  }

  iter <- object$iter

  pcors <- object$post_samp$pcors[, , 51:(iter + 50)]

  results <- sapply(1:iter, function(x) {

    pcors_s <- pcors[, , x] * adj

    est <- FUN(pcors_s, ...)

    if (!is.null(extract)) {
      est <- eval(parse(text = paste0("est$", extract)))
    }

    est
  })

  results
}

