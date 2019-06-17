#' Precision Matrix to Multiple Regression
#'
#' @description There is a direct correspondence between the covariance matrix and multiple regression. In the case of GGMs, it is possible
#' to estimate the edge set with multiple regression (i.e., neighborhood selection). In *BGGM*, the precision matrix is first sampled from, and then
#' each draws is converted to the corresponding coefficients and error variances. This results in a posterior distribution. This function can be used
#' to perform Bayesian multiple regression.
#'
#' @param fit object of class \code{estimate} (analytic = F)
#' @param node which node to summarize (i.e., the outcome)
#' @param ci_width credible interval width used in the summary output
#' @param samples number of samples used in the conversion.
#' @param ... e.g., \code{digits}
#'
#' @return summary of the multiple regresion coefficients for \code{node}
#' @export
#'
#' @examples
#'
#' # p = 10
#' Y <- BGGM::bfi[,1:10]
#'
#' # sample posterior
#' fit <- estimate(Y, samples = 5000)
#'
#' # precision to regression
#' coefficients(fit, node = 1, ci_width = 0.95)

coef.estimate <- function(fit, node, ci_width,  samples = 1000, ...){

  if(isTRUE(fit$analytic)) stop("posterior samples are required (analytic = F)")

  test <- BGGM:::beta_summary(fit, node = node, ci_width = ci_width, samples = samples)

  sums <- list()

  for(i in 1:max(node)){

    sums[[i]] <- test[[i]][[1]]
    names(sums)[[i]] <-  paste("Predicting node", node[i])
  }
  cat("BGGM: Bayesian Gaussian Graphical Models \n")
  cat("--- \n")
  cat("Type: Inverse to Regression \n")
  cat("--- \n")
  cat("Call: \n")
  print(test$call)
  cat("--- \n")
  cat("Estimates: \n \n")
  names(sums) <- ""
  sums <- data.frame( sums, check.names = F )
  colnames(sums)[2] <- "post_mean"
  print(sums[,1:5], row.names = FALSE, ...)
  cat("--- \n")
}
