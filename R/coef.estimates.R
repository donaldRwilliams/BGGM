#' Title
#'
#' @param fit
#' @param node
#' @param ci_width
#' @param samples
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
coef.estimate <- function(fit, node, ci_width,  samples = 1000, ...){
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
