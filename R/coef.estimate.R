#' Precision Matrix to Multiple Regression
#' @name coef.estimate
#' @description There is a direct correspondence between the covariance matrix and multiple regression. In the case of GGMs, it is possible
#' to estimate the edge set with multiple regression (i.e., neighborhood selection). In *BGGM*, the precision matrix is first sampled from, and then
#' each draws is converted to the corresponding coefficients and error variances. This results in a posterior distribution. This function can be used
#' to perform Bayesian multiple regression.
#'
#' @param object object of class \code{estimate} (analytic = F)
#' @param node which node to summarize (i.e., the outcome)
#' @param cred credible interval used in the summary output
#' @param iter number of samples used in the conversion.
#' @param ... e.g., \code{digits}
#'
#' @return list of class \code{coef.estimate}:
#'
#' list \code{inv_2_beta}:
#' \itemize{
#'  \item \code{betas} posterior samples for the regression coefficients
#'  \item \code{sigma} posterior samples for sigma (residual sd)
#'  \item \code{call} \code{match.call()}
#' }
#'
#' data frame \code{summary_inv_2_beta}:
#' \itemize{
#' \item summary of regression coefficients
#' }
#'
#'
#' \code{call} \code{match.call()}
#'
#' @examples
#' # p = 10
#' Y <- BGGM::bfi[,1:10]
#'
#' # sample posterior
#' fit <- estimate(Y, iter = 5000)
#'
#' # precision to regression
#' coefficients(fit, node = 1, cred = 0.95)
#' @export
coef.estimate <- function(object, node = 1, cred = 0.95, iter = 500, ...){

  # check for samples
  if(isTRUE(object$analytic)) stop("posterior samples are required (analytic = F)")

  # inverse to beta
  inv_2_beta <- beta_summary(object, node = node,
                             ci_width = cred, samples = iter)

  # summary regression coefficients
  summary_inv_2_beta <- inv_2_beta[[1]][[1]][,1:5]

  returned_object <- list(inv_2_beta = inv_2_beta,
                          summary_inv_2_beta = summary_inv_2_beta,
                          call = match.call(), data = object$dat,
                          iter = iter, node = node, cred = cred)
  class(returned_object) <- c("BGGM", "estimate",  "coef")
  returned_object
}


