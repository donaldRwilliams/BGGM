#' Gaussian Graphical Models with Exploratory Bayesian Hypothesis Testing
#'
#' @description Learn the conditional (in)dependence structure with the Bayes factor computed from the matrix-F prior distribution. It is
#' possible to test for only positive or negative edges, as well as two sided hypothesis testing (which is the customary approach). Further
#' there is also an exhaustive option that provides the posterior probability of the null, greater than zero, and less than zero.
#'
#' @param X data matrix  (\emph{n} $\times$  \emph{p}).
#' @param iter number of posterior samples.
#' @param cores number of cores for parallel computing. The default is 2, but this can be adjusted.
#' @param prior_sd hypothesized standard deviation of the prior distribution.
#'
#' @return list of class \code{explore}:
#' \itemize{
#' \item \code{parcors_mat} partial correlation matrix
#' \item \code{parcors_sd} partial correlation standard deviations
#' \item \code{samples} list of prior and posterior samples
#' \itemize{
#'  \item \code{fisher_z_post} Fisher z transformed posterior distributions (partial correlations)
#'  \item  \code{pcor_post} partial correlation posterior distribtions (not tranformed)
#'  \item \code{inv_cov_post} inverse covariance matrix posterior distribtions
#'  \item \code{pcor_prior} partial correlation prior distribution
#'  \item \code{fisher_z_prior}  Fisher z transformed prior distributions (partial correlations)
#'  }
#' \item \code{delta} hyperparameter of matrix-F prior distribution (corresponds to \code{prior_sd})
#' \item \code{iter} number of posterior and prior samples
#' \item \code{dat} data matrix
#' \item \code{call} match.call()
#' \item \code{p} number of variables
#' \item \code{cores} number of cores
#' \item \code{edge} number of estimated edges
#' }
#' @note After sampling from the posterior distribution, use \code{select} to determine the edge set and \code{plot} for visualizing the
#' edge set. see \code{methods(class = "explore")}
#'
#' @export
#'
#' @examples
#'
#' # p = 5
#' Y <- BGGM::bfi[,1:5]
#'
#' # fit model
#' fit_bf <- explore(Y, prior_sd = 0.5,
#'                  iter = 5000,
#'                  cores = 2)
#'
#'
#'# two.sided hypothesis testing
#'E <- select(fit_bf, BF_cut = 3,
#'             alternative = "two.sided")
#'
#'# 'network' plot
#'plot(E, type = "network")

explore.default <- function(X, prior_sd, iter = 5000,  cores = 2){

  delta <- BGGM:::delta_solve(prior_sd)

  X <- na.omit(X)

  p <- ncol(X)

  parcors_mat <- parcors_sd <- matrix(0, p, p)

  edges <- 0.5 * (p * (p -1))

  samples <- BGGM:::sampling(X, nu = 10000,
                             delta = delta,
                             n_samples = iter,
                             cores = cores)




  posterior_samples <- do.call(rbind.data.frame,
                               lapply(1:cores, function(x)  samples[[x]]$pcor_post))


  parcors_mat[upper.tri(parcors_mat)] <- colMeans(posterior_samples)[1:edges]

  pacors_mat <- BGGM:::symmteric_mat(parcors_mat)

  parcors_sd[upper.tri(parcors_sd)] <- apply(posterior_samples, 2,sd)[1:edges]

  pacors_sd <- BGGM:::symmteric_mat(parcors_sd)

  returned_object <- list(parcors_mat = pacors_mat,
                          parcors_sd = parcors_sd,
                          samples = samples,
                          delta = delta,
                          iter = iter,
                          dat = X,
                          call = match.call(),
                          p = p,
                          cores = cores,
                          edge = edges)

  class(returned_object) <- "explore"

  return(returned_object)

}



