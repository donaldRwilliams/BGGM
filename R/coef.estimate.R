#' Precision Matrix to Multiple Regression
#'
#' @description There is a direct correspondence between the covariance matrix and multiple regression. In the case of GGMs, it is possible
#' to estimate the edge set with multiple regression (i.e., neighborhood selection). In *BGGM*, the precision matrix is first sampled from, and then
#' each draws is converted to the corresponding coefficients and error variances. This results in a posterior distribution. This function can be used
#' to perform Bayesian multiple regression.
#'
#' @param object object of class \code{estimate} (analytic = F)
#' @param node which node to summarize (i.e., the outcome)
#' @param CrI credible interval used in the summary output
#' @param iter number of samples used in the conversion.
#' @param ... e.g., \code{digits}
#'
#' @return summary of the multiple regression coefficients for \code{node}
#' @export
#'
#' @examples
#'
#' # p = 10
#' Y <- BGGM::bfi[,1:10]
#'
#' # sample posterior
#' fit <- estimate(Y, iter = 5000)
#'
#' # precision to regression
#' coefficients(fit, node = 1, ci_width = 0.95)
coef.estimate <- function(object, node = 1, CrI = 0.95, iter = 1000, ...){

  # check for samples
  if(isTRUE(object$analytic)) stop("posterior samples are required (analytic = F)")

  # inverse to beta
  inv_2_beta <- beta_summary(object, node = node,
                             ci_width = CrI, samples = iter)

  # summary regression coefficients
  summary_inv_2_beta <- inv_2_beta[[1]][[1]][,1:5]

  # rename columns
  colnames(summary_inv_2_beta) <- c("Node", "Estimate", "Est.Error",
                                    paste(c("lb.", "ub."),  gsub("0.", "", CrI), "%", sep = ""))
  # call beta summary
  call_inv_2_beta <- inv_2_beta[[1]][[2]]

  cat("BGGM: Bayesian Gaussian Graphical Models \n")
  cat("--- \n")
  cat("Type: Inverse to Regression \n")
  cat("--- \n")
  cat("Call: \n")
  print(call_inv_2_beta)
  cat("--- \n")
  cat("Estimates: \n \n")
  suummary_inv_2_beta <- data.frame(summary_inv_2_beta, check.names = F)
  print(summary_inv_2_beta, row.names = FALSE, ...)
  cat("--- \n")
}


