#' Select Graphical Structure with the Bayes Factor
#'
#' @description Compare GGMs with the Bayes factor. This method allows for
#' assessing (relative) evidence for edge equality or edges differences across any number of groups. Further, confirmatory hypothesis testing
#' can be used to test predictions or expectations regarding difference or similiarities in different groups (e.g., male vs. female).
#'
#'
#' @param x object of class \code{select.ggm_compare_bf}
#' @param BF_cut evidentiary threshold
#'
#' @note The test provides relative evidence for whether all groups have the same edge strength (for each edge in the model). It is possibel to test whether
#' many groups are all the same. In this case, if there is evidence for the alternative (not equal), information is not immediately available for which group differs
#' from the others. Thus pairwsie group contrasts could be used after testing the equality of more than two groups.
#'
#' @return list of class \code{select.ggm_compare_bf}
#'
#' \itemize{
#' \item \code{BF_10} Bayes factors for the alternative ("not equal")
#' \item \code{BF_01} Bayes factors for the null hypothesis
#' \item \code{BF_10_adj} Bayes factor adjacency matrix for the alternative ("not equal")
#' \item \code{BF_01_adj} Bayes factor adjacency matrix for the null hypothesis
#' \item \code{adj_10} adjacency matrix for the selected edges (in favor of the "not equal")
#' \item \code{adj_01} adjacency matrix for the selected edges (in favor of the null hypothesis)
#' }
#'
#' @export
#'
#' @examples
#'
#' Y1 <- MASS::mvrnorm(500, rep(0, 16), Sigma = diag(16))
#' Y2 <- MASS::mvrnorm(500, rep(0, 16), Sigma = diag(16))
#' Y3 <- MASS::mvrnorm(500, rep(0, 16), Sigma = diag(16))
#'
#' ggm_bf <- ggm_compare_bf(Y1, Y2, Y3,
#'                         prior_sd = .5,
#'                         iter = 5000,
#'                         cores = 2)
#'# select with BF_cut = 3
#' ggm_bf_sel <- select(ggm_bf, BF_cut = 3)
#'
#' # summary
#' summary(ggm_bf_sel, type = "adj")

select.ggm_compare_bf <- function(x, BF_cut = 3){

  BF_10 <- 1/ x$BF_01

  diag(BF_10) <- 0

  adj_10 <- ifelse(BF_10 > BF_cut, 1, 0)

  adj_01 <- ifelse(x$BF_01 > BF_cut, 1, 0)

  BF_01_adj <- adj_01 * x$BF_01

  BF_10_adj <- adj_10 * BF_10

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
                          BF = BF_cut)

  class(returned_object) <- "select.ggm_compare_bf"
  returned_object

}
