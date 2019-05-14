#' Hypothesis Plot
#'
#' @description This function allows for plotting prior distributions and prior distribution differences.
#'  This can assist in selecting a hypothesized effect size, which is used for computing the the Bayes factor. Further,
#'  if a fitted object is supplied (\code{fit}), this allows for plotting the chosen prior and a selected posterior distribution (\code{edge}).
#'  This shows how the Bayes factor is computed (ratio of denisty evaluated at zero).
#'
#' @param rho_sd standard deviation of the prior distribution (see notes)
#' @param difference \code{TRUE} plots the difference between prior distributions (see notes)
#' @param fit option object of class \code{explore}
#' @param edge option for plotting selected posterior distribtion (only if a fitted model is provided)
#' @param size point size for prior and posteriro plot (only needed if a fitted model is provided)
#' @param samples number of samples drawn from the prior distirbition
#'
#' @return a \code{ggplot2} object
#' @export
#'
#' @examples
#'
#'# plot prior distributions
#'plot_priors <- hypothesis_plot(c(0.25, 0.5))
#'plot_prior
#'
#'#plot prior distribution differences
#' plot_priors <- hypothesis_plot(c(0.25, 0.5))
#' plot_prior
#'
#' # plot posterior and prior
#' Y <- BGGM::ptsd[,1:10]
#' fit <- explore(Y, prior_sd = .5, iter = 5000)
#'
#' plot_post_prior <- hypothesis_plot(fit = fit, edge = "2--3", size = 4, difference = FALSE)
#' plot_post_prior
#'
hypothesis_plot <- function(rho_sd = NULL,
                            difference = FALSE,
                            fit = NULL,
                            edge = NULL,
                            samples = 10000,
                            size = NULL){


  if(isFALSE(is.element(c("ggplot2"), installed.packages()))) stop("please install ggplot2")

  if(isFALSE(difference) && is.null(fit)){
    delta <- list()
    prior_samples <- list()

    for(i in 1:length(rho_sd)){

      prior_samples[[i]] <- (rbeta(samples, BGGM:::delta_solve(rho_sd[i]) / 2, BGGM:::delta_solve(rho_sd[i]) / 2) - 0.5) * 2

    }

    names(prior_samples) <- rho_sd
    samps <-  t(do.call(rbind, prior_samples))
    samps <-  unlist(prior_samples)

    dat <- data.frame(SD =  as.factor(rep(rho_sd, each = samples)), samples = samps)
    plt <- ggplot(dat) + stat_density(aes(x= samples,
                                          group = SD,
                                          linetype = SD),
                                      geom = "line",
                                      position = "identity",
                                      adjust = 1.5) +
      theme_bw() +
      xlab(expression("Prior Distribution   "  *rho[i][j]))

  }
  if(isTRUE(difference)){

    delta <- list()
    prior_samples <- list()

    for(i in 1:length(rho_sd)){
      # delta[[i]] <- BGGM:::delta_solve(rho_sd)
      one <- (rbeta(samples, BGGM:::delta_solve(rho_sd[i]) / 2, BGGM:::delta_solve(rho_sd[i]) / 2) - 0.5) * 2
      two <- (rbeta(samples, BGGM:::delta_solve(rho_sd[i]) / 2, BGGM:::delta_solve(rho_sd[i]) / 2) - 0.5) * 2
      prior_samples[[i]] <-  one - two

    }

    names(prior_samples) <- rho_sd
    samps <-  t(do.call(rbind, prior_samples))
    samps <-  unlist(prior_samples)

    dat <- data.frame(SD =  as.factor(rep(rho_sd, each = samples)), samples = samps)
    plt <- ggplot(dat) + stat_density(aes(x= samples,
                                          group = SD,
                                          linetype = SD),
                                      geom = "line",
                                      position = "identity",
                                      adjust = 1.5) +
      theme_bw() +
      xlab(expression("Prior Distribution (Difference) "  *rho[1][2]* "-"*rho[1][3]))
  }


  if(!is.null(fit)){
    if(class(fit) == "explore"){
      if(is.null(edge))stop("please specifcy one edge--e.g., 1--2")
      if(is.null(size))stop("size must be specified--e.g., 3")

      mat_names <- matrix(0, fit$p, fit$p)

      mat_names[] <-  unlist(lapply(1:fit$p, function(z) paste(1:fit$p, z, sep = "--")))


      post_samples <- do.call(rbind, lapply(1:fit$cores, function(x)  fit$samples[[x]]$pcor_post))

      colnames(post_samples) <-   c(mat_names[upper.tri(mat_names)], mat_names[lower.tri(mat_names)])

      edge_posterior <- post_samples[,edge]

      samples <- length(edge_posterior)
      prior_samps <-  (rbeta(samples, fit$delta / 2, fit$delta / 2) - 0.5) * 2

      dat <- data.frame(type = rep(c("posterior", "prior"), each = samples ),
                        samples = c(edge_posterior, prior_samps))


      temp_dens1 <- density(edge_posterior, n = 10000)



      dens_post <- approx(temp_dens1$x, temp_dens1$y, xout = c(0))$y
      if(is.na(dens_post)) dens_post <- 0
      dens_prior <-  dbeta(.5, fit$delta / 2, fit$delta / 2) * .5

      plt <-  ggplot(dat, aes(x = samples,
                              group = type,
                              linetype = type)) +
        stat_density(geom = "line",
                     position = "identity",
                     adjust = 1.5) +
        annotate(geom = "point", x = 0, y = dens_post, color = "red", size = size) +
        annotate(geom = "point", x = 0, y = dens_prior, color = "blue", size = size) +
        scale_linetype(name = "") +
        theme_bw()
    }
  }
  plt
}
