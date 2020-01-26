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
  class(returned_object) <- "coef.estimate"
  returned_object
}

#' @name print.coef.estimate
#' @title  Print method for \code{coef.estimate} objects
#' @param x An object of class \code{coef.estimate}
#' @param ... currently ignored
#'
#' @seealso  \code{\link{coef.estimate}}
#'
#' @export
print.coef.estimate <- function(x,...){

  res_sigma <- x$inv_2_beta$sigma

  lb <- (1 - x$cred) / 2

  ub <- 1 - lb

  cred_int <- stats::quantile(res_sigma, prob = c(lb, ub))

  res_sigma_summ <- data.frame(Estimate = mean(res_sigma),
                               Est.Error = sd(res_sigma),
                               t(cred_int))

  # R2
  ypred <- t(apply(as.matrix(x$inv_2_beta$betas)[1:x$iter,], 1,
                   function(z)  z %*% t(as.matrix(x$data[,- x$node]))))

  r2 <- R2_helper(ypred, x$data[,x$node], ci_width = x$cred)
  cred_in <- stats::quantile(r2$R2, prob = c(lb, ub))

  res_r2_summ <- data.frame(Estimate = mean(r2$R2), Est.Error = sd(r2$R2), t(cred_in))

  colnames(res_sigma_summ) <- c("Estimate", "Est.Error", "CrI.lb", "CrI.ub")

  colnames(res_r2_summ) <- c("Estimate", "Est.Error", "CrI.lb", "CrI.ub")

  colnames(x$summary_inv_2_beta) <- c("Node", "Estimate",
                                      "Est.Error", "Cred.lb",
                                      "Cred.ub")
  cat("BGGM: Bayesian Gaussian Graphical Models \n")
  cat("--- \n")
  cat("Type: Inverse to Regression \n")
  cat("Credible Interval:",  gsub("*0.","", formatC( round(x$cred, 4), format='f', digits=2)), "% \n")
  cat("Node:", x$node, "\n")
  cat("--- \n")
  cat("Call: \n")
  print(x$call)
  cat("--- \n")
  cat("Coefficients: \n \n")
  summary_inv_2_beta <- data.frame(x$summary_inv_2_beta,
                                   check.names = F)
  print(summary_inv_2_beta, row.names = FALSE, ...)
  cat("--- \n")
  cat("Sigma:\n\n")
  print(round(res_sigma_summ, 3), row.names = FALSE, ...)
  cat("--- \n")
  cat("Bayesian R2:\n\n")
  print(round(res_r2_summ, 3), row.names = FALSE, ...)
  cat("--- \n")
}
