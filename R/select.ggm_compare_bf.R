#' @title Graph selection for \code{ggm_compare_explore} Objects
#'
#' @description Provides the selected graph (of differences) based on the Bayes factor
#' \insertCite{williams2020comparing}{BGGM}.
#'
#' @param object An object of class \code{ggm_compare_explore}.
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
#' # data
#' Y <- bfi
#'
#' # males and females
#' Ymale <- subset(Y, gender == 1,
#'                    select = -c(gender,
#'                                education))[,1:10]
#'
#' Yfemale <- subset(Y, gender == 2,
#'                      select = -c(gender,
#'                                  education))[,1:10]
#'
#' # fit model
#' fit <- ggm_compare_explore(Ymale, Yfemale,
#'                            iter = 250,
#'                            type = "continuous",
#'                            progress = FALSE)
#'
#'
#' E <- select(fit, post_prob = 0.50)
#'
#' }
#'
#' @export
select.ggm_compare_explore <- function(object,
                                       BF_cut = 3,
                                       ...){

  # change to x
  x <- object

  if(is.null(BF_cut)){

    BF_cut <- (post_prob) / (1 - post_prob)

    } else {

      BF_cut <- BF_cut

  }

  # post
  post_prob <- BF_cut / (BF_cut + 1)

  # BF
  BF_10 <- 1 / x$BF_01

  # BF mat diagonal
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


print_select_ggm_compare_bf <- function(x,...){

  groups <- x$object$groups
  p <- x$p
  cat("BGGM: Bayesian Gaussian Graphical Models \n")
  cat("--- \n")
  cat("Type:",  x$object$type, "\n")
  # number of iterations
  cat("Posterior Samples:", x$object$iter, "\n")
  # number of observations
  cat("Observations (n):\n")
  groups <- length(x$object$info$dat)
  for(i in 1:groups){
    cat("  Group", paste( i, ":", sep = "") , x$object$info$dat_info$n[[i]], "\n")
  }
  # number of variables
  cat("Variables (p):", x$object$p, "\n")
  # number of edges
  cat("Relations:", .5 * (x$object$p * (x$object$p-1)), "\n")
  cat("Delta:", x$object$delta, "\n")
  cat("--- \n")
  cat("Call: \n")
  print(x$object$call)
  cat("--- \n")
  cat("Hypotheses:\n")
  cat("H0:", paste0("rho_g", 1:groups, collapse = " = "), "\n")
  cat("H1:", paste0("rho_g", 1:groups, collapse = " - "), " = 0\n")
  cat("--- \n\n")
  if(groups ==2){
    cat("Partial Correlations:\n\n")
    colnames(x$pcor_mat_10) <- 1:p
    row.names(x$pcor_mat_10) <- 1:p
    print(round(x$pcor_mat_10, 2))
    cat("--- \n")
  }
  cat("Adjacency:\n\n")
  colnames(x$adj_10) <- 1:p
  row.names(x$adj_10) <- 1:p
  print(round(x$adj_10, 2))

}

