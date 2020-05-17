#' Extract the Weighted Adjacency Matrix
#'
#' @description Extract the weighted adjacency matrix (posterior mean) from
#' \code{\link{estimate}}, \code{\link{explore}}, \code{\link{ggm_compare_estimate}},
#' and \code{\link{ggm_compare_explore}} objects.
#'
#' @param object A model estimated with \strong{BGGM}. All classes are supported, assuming
#' there is matrix to be extracted.
#'
#' @param ... Currently ignored.
#'
#' @return The weighted adjacency matrix (partial correlation matrix with zeros).
#'
#' @examples
#' \donttest{
#' # note: iter = 250 for demonstrative purposes
#' Y <- bfi
#'
#' # estimate
#' fit <- estimate(Y, iter = 250)
#'
#' # select graph
#' E <- select(fit)
#'
#' # extract weighted adj matrix
#' weighted_adj_mat(E)
#'
#' }
#' @export
weighted_adj_mat <- function(object, ...){


  if(all(c("select.estimate", "estimate") %in% class(object))){

    weighted_adj_mat <- round(object$pcor_adj, 3)
    weighted_adj_mat

  } else if(all(c("select.estimate", "estimate") %in% class(object))){

    weighted_adj_mat <- round(object$pcor_mat_zero, 3)
    weighted_adj_mat

  } else if(all(c("select.ggm_compare_estimate", "estimate") %in% class(object))){

    contrasts <- length(object$pcor_adj)

    weighted_adj_mat <- lapply(1:length(contrasts), function(x) round(object$pcor_adj[[x]], 3))
    names(weighted_adj_mat) <- names(object$object$pcor_mats)
    weighted_adj_mat

  } else if(c("select.ggm_compare_bf") %in% class(object)){

    if(object$object$groups > 2){

      stop("weigthed adjacency only available for two groups")
    }

    weighted_adj_mat <- round(object$pcor_mat_10, 3)
    weighted_adj_mat

  } else {

    stop("weighted adjacency matrix not found.")

  }

}



