#' Select Partial Correlattion Differences for \code{ggm_compare_bf} Objects
#'
#' @param object object of class \code{ggm_compare_bf}.
#'
#' @param post_prob numeric. Posterior `inclusion` probability.  The default is set to 0.50.
#'
#' @param ... not currently used.
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
#' \item \code{pcor_mat_10} partial correlation matrix for the alternative ("not equal")
#' }
#'
#' @export
#'
#' @examples
#' \donttest{
#' Y1 <- MASS::mvrnorm(500, rep(0, 16), Sigma = diag(16))
#' Y2 <- MASS::mvrnorm(500, rep(0, 16), Sigma = diag(16))
#' Y3 <- MASS::mvrnorm(500, rep(0, 16), Sigma = diag(16))
#'
#' ggm_bf <- ggm_compare_bf(Y1, Y2, Y3,
#'                         prior_sd = .20,
#'                         iter = 500,
#'                         cores = 2)
#'# select with BF_cut = 3
#' ggm_bf_sel <- select(ggm_bf, BF_cut = 3)
#'
#' # summary
#' summary(ggm_bf_sel)
#' }

select.ggm_compare_explore <- function(object,
                                       post_prob = 0.50,
                                       BF = NULL,...){

  x <- object

  if(is.null(BF)){

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
