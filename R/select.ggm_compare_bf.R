#' @title Graph selection for \code{explore} Objects
#'
#' @description Provides the selected graph (of differences) based on the Bayes factor
#' \inserCite{williams2020comparing}{BGGM}.
#'
#' @param post_prob Numeric. Posterior `inclusion` probability (defaults to 0.50)
#'                  for including an edge.
#'
#' @param BF_cut Numeric. Threshold for including an edge (defaults to 3).
#'
#' @param ... Currently ignored.
#'
#' @return The returned object of class \code{select.ggm_compare_explore} contains
#' a lot of information that is used for printing and plotting the results.
#' For users of \strong{BGGM}, the following are the useful objects:
#'
#'
#' \itemize{
#'
#' \item \code{adj_10} Adjacency matrix for which there was evidence for a difference.
#'
#' \item \code{adj_10} Adjacency matrix for which there was evidence for a null relation
#'
#' \item \code{pcor_mat_10} Selected partial correlation matrix (weighted adjacency; only for two groups).
#'
#' }
#'
#' @seealso \code{\link{explore}} and \code{\link{ggm_compare_explore}} for several examples.
#'
#' @examples
#' \donttest{
#'
#' ##################
#' ### example 1: ###
#' ##################
#'
#' # fit model
#' fit <- ggm_compare_explore(Ymale, Yfemale,
#'                            iter = 250,
#'                            type = "continuous")
#'
#'
#' E <- select(fit, post_prob = 0.50)
#'
#' }
#'
#' @export
select.ggm_compare_explore <- function(object,
                                       post_prob = 0.50,
                                       BF_cut = NULL,...){

  x <- object

  if(is.null(BF_cut)){

    BF_cut = (post_prob) / (1 - post_prob)

    } else {

      BF_cut <- post_prob

  }

  BF_10 <- 1/ x$BF_01

  diag(BF_10) <- 0

  adj_10 <- ifelse(BF_10 > BF_cut, 1, 0)

  adj_01 <- ifelse(x$BF_01 > BF_cut, 1, 0)

  BF_01_adj <- adj_01 * x$BF_01

  BF_10_adj <- adj_10 * BF_10

  pcor_mat <- matrix(0, x$p, x$p)

  if(x$groups == 2){

    pcor_mat <- adj_10 * x$pcor_diff

    }

  returned_object <- list(BF_10 = BF_10,
                          BF_01 = x$BF_01,
                          BF_01_adj = BF_01_adj,
                          BF_10_adj = BF_10_adj,
                          adj_10 = adj_10,
                          adj_01 = adj_01,
                          call = match.call(),
                          p = ncol(BF_10),
                          iter = x$iter,
                          info = x$info,
                          post_prob = post_prob,
                          BF = BF,
                          pcor_mat_10 = pcor_mat,
                          object = object)

  class(returned_object) <- c("BGGM",
                              "explore",
                              "select",
                              "select.ggm_compare_bf")
  returned_object

}
