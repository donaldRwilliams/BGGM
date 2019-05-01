#' estimate
#'
#' @description Samples from the posterior with the Wishart prior distribution.
#' @param x data matrix
#' @param samples number of posterior samples
#'
#' @return  An object of class \code{estimate} that includes posterior samples, partial correlation matrix, etc. Note
#' this object is then used with other functions for model selection, prediction, regression coverion, etc.
#' For Bayesian hypothesis testing see \code{explore} and \code{confirm}.
#'
#' @export
#'
#' @examples
#' fit <- estimate(X)
#' fit

estimate.default  <- function(x, samples = 5000, ...){

  # remove the NAs
  X <- na.omit(as.matrix(x))

  # mean center the data
  X <- scale(X, scale = F)


  # number of observations
  n <- nrow(X)

  # number of columns
  p <- ncol(X)

  # number of columns inv + pcor
  cols_samps <- p^2 + p^2

  S <- t(X) %*% X


  # store posterior samples
  df_samps <- matrix(0, nrow = samples, ncol = cols_samps)

  for(i in 1:samples){
    # draw directly from Wishart
    inv_mat <- rWishart(1, n-1, solve(S))[,,1]

    # compute partial correlations
    pcor_mat <-   -1 * cov2cor(inv_mat)

    # into the $i$th row
    df_samps[i,1:cols_samps] <- c(as.numeric(inv_mat), as.numeric(pcor_mat))

  }

  # name the columns
  inv_names <- unlist(lapply(1:p, function(x)  BGGM:::samps_inv_helper(x, p)))
  pcor_names <-unlist(lapply(1:p, function(x)  BGGM:::samps_pcor_helper(x, p)))

  colnames(df_samps) <- c(inv_names, pcor_names)

  parcor_mat <- matrix(0, ncol = ncol(x), ncol(x))
  inv_mat <-  parcor_mat

  parcor_mat[] <- colMeans(df_samps[,  grep("pcors", colnames(df_samps))])
  diag(parcor_mat) <- 0

  inv_mat[]   <- colMeans(df_samps[,  grep("cov_inv", colnames(df_samps))])

  returned_object  <- list(parcors_mat = parcor_mat,
                         inv_mat = inv_mat,
                         posterior_samples = as.data.frame(df_samps),
                         p = ncol(x),
                         dat = X,
                         iter = samples, call = match.call())

class(returned_object) <- "estimate"

return(returned_object)
}
