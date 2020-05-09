#' @title  Plot: Prior Distribution
#'
#' @description Visualize the implied prior distribution for the partial correlations. This is
#'              particularly useful for the Bayesian hypothesis testing methods.
#'
#' @name plot_prior
#'
#' @description Visualize the prior distribution used for the
#' Bayesian hypothesis testing methods. This is governed by \code{prior_sd}.
#'
#' @param delta numeric. Hyperparameter that governs the scale of
#'              the prior distribution. Larger values are more informative
#'              (i.e., narrower).
#'
#'  @param iter number of iterations (posterior samples; defaults to 5000).
#'
#' @return a \code{ggplot} object.
#' @export
#'
#' @examples

plot_prior <- function(prior_sd = 0.2,
                       iter = 5000){


  # matrix dimensions for prior
  Y_dummy <- matrix(rnorm(10 * 3),
                    nrow = 10, ncol = 3)

  delta <- BGGM:::delta_solve(prior_sd)

  # sample prior
  prior_samp <- .Call('_BGGM_sample_prior',
                      PACKAGE = 'BGGM',
                      Y = Y_dummy,
                      iter = iter,
                      delta = delta,
                      epsilon = 0.01,
                      prior_only = 1,
                      explore = 1)

  qplot(prior_samp$pcors[1,2,], geom = "density") +
    xlab("Implied Prior Distribution")

  }
