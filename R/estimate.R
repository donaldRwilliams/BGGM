#' Gaussian Graphical Models with Credible Intervals or the Region of Practical Equivalence
#'
#' @description Learn the conditional (in)dependence structure with credible intervals or the region of practical equivalence.
#' For the former, there is an analytic solution avaiable, whereas for the latter, samples are efficiently drawn from the posterior
#' distribution.
#'
#' @param x data matrix (\emph{n} $\times$  \emph{p}).
#' @param iter number of posterior samples
#' @param analytic analytic solution. see notes for futher details.
#'
#' @return  list of class \code{estimate}:
#'
#' \code{analytic = TRUE}:
#' \itemize{
#' \item \code{fit} list of analytic solution estimates
#' \itemize{
#' \item \code{inv_mu} inverse covariance matrix (mean)
#' \item \code{inv_var} inverse covariance matrix (variance)
#' \item \code{partial} partial correlation matrix
#' }
#' \item \code{analytic} TRUE
#' \item \code{call} match.call()
#' \item \code{dat} data matrix
#' \item \code{p} number of variables
#' }
#'
#' \code{analytic = FALSE}:
#' \itemize{
#' \item \code{parcors_mat} partial correlation matrix
#' \item \code{inv_mat} inverse covariance matrix
#' \item \code{posterior samples} posterior samples for partial correlations and inverse covariance matrix
#' \item \code{p} number of variables
#' \item \code{dat} data matrix
#' \item \code{iter} number of posterior samples
#' \item \code{call} match.call()
#' \item \code{analytic} FALSE
#' }
#'
#' @export
#'
#'
#' @note The default is to draws samples from the posterior distribution (\code{analytic = FALSE}). The samples are required for computing edge differences,
#' Bayesian R2, etc. If the goal is to *only* determined the non-zero effects, this can be accomplished by setting \code{analytic = TRUE}. This is accomplished
#' by estimating the posterior mean and variance, from which the credible intervals can computed. Note also sampling is very fast--i.e., less than 1 second
#' with p = 25, n = 2500 and 5,000 samples. There is one function that makes use of the analytic solution. Namely, \code{loocv} computes node-wise leave-one-out
#' error (also analytically).
#'
#'see \code{methods("estimate")}

#' @examples
#'
#' # p = 5
#' Y <- BGGM::bfi[,1:20]
#'
#' # analytic approach (sample by setting analytic = F)
#' fit_analytic <- estimate(Y, analytic = T)
#'
#' # select the graph (edge set E)
#' E <- select(fit_analytic, ci_width = 0.95)
#'
#' # plot
#' plot(E, type = "network")

estimate.default  <- function(x, iter = 5000, analytic = FALSE){

  # remove the NAs
  X <- na.omit(as.matrix(x))

  # mean center the data
  X <- scale(X, scale = T)


  # number of observations
  n <- nrow(X)

  # number of columns
  p <- ncol(X)

  # number of columns inv + pcor
  cols_samps <- p^2 + p^2

  S <- t(X) %*% X


  if(isFALSE(analytic)){
  # store posterior samples
  df_samps <- matrix(0, nrow = iter, ncol = cols_samps)

  for(i in 1:iter){
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
                         iter = iter, call = match.call(),
                         analytic = analytic)
  } else{


    fit <-  BGGM:::analytic_solve(X)

    returned_object <- list(fit = fit, analytic = analytic, call = match.call(), data = X, p = ncol(X))

    }

class(returned_object) <- "estimate"

return(returned_object)
}


