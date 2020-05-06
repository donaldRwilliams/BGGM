#' Prior Distribution Plot
#'
#' @plot_prior
#'
#' @description Visualize the prior distribution used for the
#' Bayesian hypothesis testing methods. This is governed by \code{prior_sd}.
#' @param prior_sd hypothesized standard deviation for the edges or partial correlations.
#' @param samples number of samples drawn from the prior distribution
#'
#' @return \code{ggplot object}
#' @export
#'
#' @examples
#' plot_prior(prior_sd = c(0.1, 0.2, 0.3), samples = 1000)
plot_prior <- function(prior_sd = 0.2, samples = 10000){

  prior_samples <- list()

  # loop through rho_sd
  for(i in 1:length(prior_sd)){

    # sample from the prior distribution
    prior_samples[[i]] <- (rbeta(samples, BGGM:::delta_solve(prior_sd[i]) / 2,
                                 BGGM:::delta_solve(prior_sd[i]) / 2) - 0.5) * 2

  }

  # name prior samples for plotting
  names(prior_samples) <- prior_sd

  # unlist prior samples
  samps <-  unlist(prior_samples)

  # data frame for plotting
  dat <- data.frame(SD =  as.factor(rep(prior_sd, each = samples)),
                    samples = samps)

  # plot
  plt <- ggplot(dat) +
    stat_density(aes(x= samples,
                     group = SD,
                     linetype = SD),
                 geom = "line",
                 position = "identity",
                 adjust = 1.5) +
    theme_bw() +
    xlab(expression("Prior Distribution   "  *rho[i][j]))

  plt

}
