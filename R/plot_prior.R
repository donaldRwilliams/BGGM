#' @title  Plot: Prior Distribution
#'
#' @description Visualize the implied prior distribution for the partial correlations. This is
#'              particularly useful for the Bayesian hypothesis testing methods.
#'
#' @name plot_prior
#'
#'
#' @param prior_sd Scale of the prior distribution, approximately the standard deviation
#'                 of a beta distribution (defaults to 0.25).
#'
#' @param iter Number of iterations (prior samples; defaults to 5000).
#'
#' @return A \code{ggplot} object.
#'
#' @examples
#' plot_prior(prior_sd = 0.25, iter = 5000)
#'
#' @export
plot_prior <- function(prior_sd = 0.2, iter = 5000){


  # matrix dimensions for prior
  Y_dummy <- matrix(rnorm(10 * 3),
                    nrow = 10, ncol = 3)

  delta <- delta_solve(prior_sd)

  # sample prior
  prior_samp <- .Call('_BGGM_sample_prior',
                      PACKAGE = 'BGGM',
                      Y = Y_dummy,
                      iter = iter,
                      delta = delta,
                      epsilon = 0.01,
                      prior_only = 1,
                      explore = 1,
                      progress = FALSE
                      )

  qplot(prior_samp$pcors[1,2,], geom = "density") +
    xlab("Implied Prior Distribution")

  }
