#' @title GGMs with Credible Intervals or the Region of Practical Equivalence
#'
#' @description Estimate the conditional (in)dependence structure with credible intervals or the region of practical equivalence.
#' For the former, there is an analytic solution available, whereas for the latter, samples are efficiently drawn from the posterior
#' distribution.
#'
#' @param Y data matrix (\emph{n} by  \emph{p}).
#' @param iter number of posterior samples
#' @param analytic analytic solution. see notes for further details.
#' @param ... not used
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
#' # p = 20
#' Y <- BGGM::bfi[, 1:5]
#'
#' # analytic approach (sample by setting analytic = FALSE)
#' fit_analytic <- estimate(Y, analytic = TRUE)
#'
#' # select the graph (edge set E)
#' E <- select(fit_analytic, ci_width = 0.95)
#'
#' @export
estimate.default  <- function(Y, iter = 5000,
                              analytic = FALSE, ...){

  x <- Y
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

  # scatter matrix
  S <- t(X) %*% X

  # sample from Wishart
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
  inv_names <- unlist(lapply(1:p, function(x)  samps_inv_helper(x, p)))
  pcor_names <-unlist(lapply(1:p, function(x)  samps_pcor_helper(x, p)))

  # name columns
  colnames(df_samps) <- c(inv_names, pcor_names)

  # matrix for storage
  parcor_mat <-  inv_mat <- matrix(0, ncol = p, p)

  # posterior means (partials)
  parcor_mat[] <- colMeans(df_samps[,  grep("pcors", colnames(df_samps))])
  diag(parcor_mat) <- 0

  # posterior means (inverse)
  inv_mat[]   <- colMeans(df_samps[,  grep("cov_inv", colnames(df_samps))])

  returned_object  <- list(parcors_mat = parcor_mat, inv_mat = inv_mat,
                           posterior_samples = as.data.frame(df_samps),
                           p = ncol(x), dat = X, iter = iter,
                           call = match.call(), analytic = analytic)
  # analytic solution
  } else {

    # analytic
    fit <-  analytic_solve(X)

    returned_object <- list(fit = fit,
                            analytic = analytic,
                            call = match.call(),
                            data = X,
                            p = ncol(X))
    }

  class(returned_object) <- "estimate"
  return(returned_object)

}

#' @title S3 estimate method
#' @name estimate
#' @param ... currently not used
#'
#' @description S3 estimate method
#' @seealso \code{\link{estimate.default}}
#' @export
estimate <- function(...) {
  UseMethod("estimate")
}

#' @name print.estimate
#' @title  Print method for \code{estimate.default} objects
#'
#' @param x An object of class \code{estimate}
#' @param ... currently ignored
#' @seealso \code{\link{select.estimate}}
#' @examples
#' # data
#' Y <- BGGM::bfi[, 1:5]
#' # analytic approach (sample by setting analytic = FALSE)
#' fit <- estimate(Y, analytic = TRUE)
#' fit
#' @export
print.estimate <- function(x, ...){
  cat("BGGM: Bayesian Gaussian Graphical Models \n")
  cat("--- \n")
  # analytic == TRUE
  if(!isFALSE( x$analytic)){
    cat("Type: Estimation (Analytic Solution) \n")
  }
  # analytic  == FALSE
  if(isFALSE( x$analytic)){
    cat("Type: Estimation (Sampling) \n")
  }
  # number of iterations
  cat("Posterior Samples:", x$iter, "\n")
  # number of observations
  cat("Observations (n):", nrow(x$dat), "\n")
  # number of variables
  cat("Variables (p):", x$p, "\n")
  # number of edges
  cat("Edges:", .5 * (x$p * (x$p-1)), "\n")
  cat("--- \n")
  cat("Call: \n")
  print(x$call)
  cat("--- \n")
  cat("Date:", date(), "\n")
}


#' @name summary.estimate
#' @title Summary method for \code{estimate.default} objects
#'
#' @param object An object of class \code{estimate}
#' @seealso \code{\link{select.estimate}}
#' @param ... currently ignored
#' @param cred credible interval width
#' @return A list containing the summarized posterior distributions
#' # data
#' Y <- BGGM::bfi[, 1:5]
#' # analytic approach (sample by setting analytic = FALSE)
#' fit <- estimate(Y, analytic = TRUE)
#' summary(fit)
#' @export
summary.estimate <- function(object, cred = 0.95, ...) {

  if (isTRUE(object$analytic)) {
    returned_object <- list(object = object)
    } else{

    lb <- (1 - cred) / 2
    ub <- 1 - lb

    name_temp <- matrix(0, object$p, object$p)

    edge_names <- unlist(lapply(1:object$p , function(x)
      paste(1:object$p, x, sep = "--")))

    name_temp[] <- edge_names
    up_tri <- name_temp[upper.tri(name_temp)]

    pcor_samples <-
      object$posterior_samples[,  grep("pcors", colnames(object$posterior_samples))]

    colnames(pcor_samples) <- edge_names
    pcor_upper <- pcor_samples[, up_tri]

    ci <- apply(
      pcor_upper,
      MARGIN = 2,
      FUN = function(x) {
        quantile(x, probs = c(lb, ub))
      }
    )
    diff_mu <-
      apply(pcor_upper, MARGIN = 2, mean)

    diff_sd <-
      apply(pcor_upper, MARGIN = 2, sd)

    dat_results <-
      data.frame(
        edge = name_temp[upper.tri(name_temp)],
        post_mean =  round(diff_mu, 3),
        post_sd = round(diff_sd, 3),
        ci = round(t(ci), 3)
      )

    colnames(dat_results) <- c(
      "Edge",
      "Estimate",
      "Est.Error",
      paste("lb.", gsub(
        "*0.", "", formatC(round(cred, 4), format = 'f', digits = 2)
      ), "%", sep = ""),
      paste("ub.", gsub(
        "*0.", "", formatC(round(cred, 4), format = 'f', digits = 2)
      ), "%", sep = "")
    )

    returned_object <- list(dat_results = dat_results,
                            object = object,
                            pcor_samples = pcor_samples)
  }

  class(returned_object) <- "summary.estimate"
  returned_object
}

#' @title Summary method for \code{summary.estimate} objects
#' @name print.summary.estimate
#'
#' @param x An object of class \code{summary.estimate}
#' @param ... currently ignored
#' @seealso \code{\link{summary.estimate}}
#' @export
print.summary.estimate <- function(x, ...) {
  # analytic == TRUE
  if (isTRUE(x$object$analytic)) {
    cat(print(x$object), "\n")
    cat("note: posterior summary not available for analytic solution")
  } else {
    cat("BGGM: Bayesian Gaussian Graphical Models \n")
    cat("--- \n")
    # number of iterations
    cat("Posterior Samples:", x$object$iter, "\n")
    # number of observations
    cat("Observations (n):", nrow(x$object$dat), "\n")
    # number of variables
    cat("Variables (p):", x$p, "\n")
    # number of edges
    cat("Edges:", .5 * (x$object$p * (x$object$p - 1)), "\n")
    cat("--- \n")
    cat("Call: \n")
    print(x$object$call)
    cat("--- \n")
    cat("Estimates:\n\n")
    print(x$dat_results, row.names = F)
    cat("--- \n")

  }
}
