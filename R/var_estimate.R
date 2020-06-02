#' VAR: Estimation
#'
#' @param Y Matrix (or data frame) of dimensions \emph{n} (observations) by  \emph{p} (variables).
#' @param rho_sd Numeric. Scale of the prior distribution for the partial correlations,
#' approximately the standard deviation of a beta distribution
#' (defaults to 0.50).
#' @param beta_sd Numeric. Standard deviation of the prior distribution for the regression coefficients
#'        (defaults to 1).
#'
#' @param iter Number of iterations (posterior samples; defaults to 5000).
#'
#' @return An object of class \code{VAR_estimate}
#' @export
VAR_estimate <- function(Y, rho_sd = 0.50,
                         beta_sd = 1,
                         iter = 5000) {

  # number of nodes
  p <- ncol(Y)

  # Y lagged: add NA
  Y_lag <-   rbind(NA, Y)

  # Y lagged: column names
  colnames(Y_lag) <- paste0(colnames(Y), ".l1")

  # combine all data
  Y_all <- na.omit(cbind.data.frame(rbind(Y, NA), Y_lag))

  # nodes in GGM
  Y <- scale(as.matrix(Y_all[,1:p]), scale = TRUE)

  # predictors (lagged effects)
  X <- scale(as.matrix(Y_all[,(p+1):(p*2)]), scale = TRUE)

  # delta: rho ~ beta(delta/2, delta/2)
  delta <- BGGM:::delta_solve(rho_sd)

  # prior variance
  beta_var <- beta_sd^2

  message(paste0("BGGM: Posterior Sampling "))

  fit <-.Call(
      "_BGGM_var",
      Y =  as.matrix(Y),
      X = as.matrix(X),
      delta = delta,
      epsilon = 0.001,
      beta_prior = diag(p) * (1 / beta_var),
      iter = iter + 50,
      start = solve(cor(Y)),
      progress = TRUE
    )
  message("BGGM: Finished")

  returned_object <- list(fit = fit,
                          iter = iter,
                          p = p,
                          Y = Y,
                          X = X)

  class(returned_object) <- c("BGGM",
                              "VAR",
                              "VAR_estimate")
  return(returned_object)
}
