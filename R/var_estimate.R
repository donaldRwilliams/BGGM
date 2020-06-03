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
var_estimate <- function(Y, rho_sd = 0.50,
                         beta_sd = 1,
                         iter = 5000) {

  # number of nodes
  p <- ncol(Y)

  # number of obs
  n <- nrow(Y)

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
                          n = n,
                          Y = Y,
                          X = X,
                          call = match.call())

  class(returned_object) <- c("BGGM",
                              "var_estimate",
                              "default")
  return(returned_object)
}


print_var_estimate <- function(x, ...){
  cat("BGGM: Bayesian Gaussian Graphical Models \n")
  cat("--- \n")
  cat("Vector Autoregressive Model (VAR) \n")
  cat("--- \n")
  # number of iterations
  cat("Posterior Samples:", x$iter, "\n")
  # number of observations
  cat("Observations (n):", x$n,"\n")
  # number of variables
  cat("Nodes (p):", x$p, "\n")
  # number of edges
  cat("--- \n")
  cat("Call: \n")
  print(x$call)
  cat("--- \n")
  cat("Date:", date(), "\n")
}


#' @title Summary Method for \code{var_estimate} Objects
#'
#' @name summary.var_estimate
#'
#' @description Summarize the posterior distribution of each partial correlation
#' and regression coefficient with the posterior mean, standard deviation, and
#' credible intervals.
#'
#' @param object An object of class \code{var_estimate}
#'
#' @param cred Numeric. The credible interval width for summarizing the posterior
#' distributions (defaults to 0.95; must be between 0 and 1).
#'
#' @param ... Currently ignored.
#'
#' @seealso \code{\link{var_estimate}}
#'
#' @return A dataframe containing the summarized posterior distributions,
#' including both the partial correlations and the regression coefficients.
#' @export
summary.var_estimate <- function(object,
                                 cred = 0.95,
                                 ...){

  # nodes
  p <- object$p

  # identity matrix
  I_p <- diag(p)

  # lower bound
  lb <- (1 - cred) / 2

  # upper bound
  ub <- 1 - lb

  # column names
  cn <-  colnames(object$Y)

  if(is.null(cn)){

    mat_names <- sapply(1:p , function(x) paste(1:p, x, sep = "--"))[upper.tri(I_p)]

  } else {

    mat_names <-  sapply(cn , function(x) paste(cn, x, sep = "--"))[upper.tri(I_p)]

  }

  pcor_mean <- round(
    apply(object$fit$pcors[, , 51:(object$iter + 50)], 1:2, mean),
    3)[upper.tri(I_p)]

  pcor_sd  <- round(
    apply(object$fit$pcors[,, 51:(object$iter + 50) ], 1:2, sd),
    digits = 3)[upper.tri(I_p)]

  pcor_lb <- round(
    apply( object$fit$pcors[,, 51:(object$iter + 50) ], 1:2, quantile, lb),
    digits = 3)[upper.tri(I_p)]

  pcor_ub <- round(
    apply(object$fit$pcors[,, 51:(object$iter + 50) ], 1:2, quantile, ub),
    digits =  3)[upper.tri(I_p)]

  beta_mean <- round(
    apply(object$fit$beta[,, 51:(object$iter + 50) ], 1:2, mean),
    digits = 3)

  beta_sd  <- round(
    apply(object$fit$beta[,, 51:(object$iter + 50) ], 1:2, sd),
    digits = 3)

  beta_lb <- round(
    apply( object$fit$beta[,, 51:(object$iter + 50) ], 1:2, quantile, lb),
    digits =  3)

  beta_ub <- round(
    apply(object$fit$beta[,, 51:(object$iter + 50) ], 1:2, quantile, ub),
    digits = 3)

  pcor_results <-
    data.frame(
      relation = mat_names,
      post_mean =  pcor_mean,
      post_sd = pcor_sd,
      post_lb = pcor_lb,
      post_ub = pcor_ub
    )

  colnames(pcor_results) <- c(
    "Relation",
    "Post.mean",
    "Post.sd",
    "Cred.lb",
    "Cred.ub")

  beta_results <-
    lapply(1:p, function (x) {
      res_p <- data.frame(
        relation = colnames(object$X),
        post_mean =  beta_mean[, x],
        post_sd = beta_sd[, x],
        post_lb = beta_lb[, x],
        post_ub = beta_ub[, x]
      )

      colnames(res_p) <- c("Relation",
                           "Post.mean",
                           "Post.sd",
                           "Cred.lb",
                           "Cred.ub")
      res_p
    })

  names(beta_results) <- colnames(object$Y)

  returned_object <- list(pcor_results = pcor_results,
                          beta_results = beta_results)

  class(returned_object) <- c("BGGM",
                              "var_estimate",
                              "summary.var_estimate")
  return(returned_object)


}

