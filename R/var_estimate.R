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

  p <- ncol(Y)

  Y_lag <-   rbind(NA, Y)

  colnames(Y_lag) <- paste0(colnames(Y), ".l1")

  Y_all <- na.omit(cbind.data.frame(rbind(Y, NA), Y_lag))

  Y <- scale(as.matrix(Y_all[,1:p]), scale = TRUE)

  X <- scale(as.matrix(Y_all[,(p+1):(p*2)]), scale = TRUE)

  delta <- BGGM:::delta_solve(rho_sd)

  beta_var <- beta_sd^2

  fit <- BGGM:::var(Y =  as.matrix( Y ),
                    X = as.matrix( X ),
                    delta = delta,
                    epsilon = 0.001,
                    beta_prior = diag(p) * (1 / beta_var),
                    iter = iter + 50,
                    start = diag(p),
                    progress = TRUE)

  returned_object <- fit
  class(returned_object) <- c("BGGM", "VAR", "VAR_estimate")
  returned_object

}
