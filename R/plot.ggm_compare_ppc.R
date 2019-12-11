#' Plot Posterior Predictive GGM Comparison
#'
#' @description Plot \href{https://cran.r-project.org/web/packages/ggridges/vignettes/introduction.html}{ggridges} for
#' the GGM comparison with posterior predictive KL-divergence. The plots contain the predictive distribution, assuming group equality, as well as
#' the observed KL-divergence. Further, the predictive distributions are conveniently colored to infer whether the null of group equality
#' should be rejected. This is accomplished by having the critical region, corresponding to a desired 'significance' level, shaded in red.
#' Thus, if the observed value is in the red region, this suggests the null hypothesis of group equality should be rejected.
#'
#' @param x object of class \code{ggm_compare_ppc}
#' @param critical 'significance' level
#' @param col_noncritical fill color of non critical region
#' @param col_critical  fill color of critical region (e.g., \code{critical = 0.05})
#' @param point_size point size for the observed KL-divergence
#' @param log log transformation. useful for small values and skewed predictive distributions
#'
#' @references
#' Williams, D. R., Rast, P., Pericchi, L. R., & Mulder, J. (2019). Comparing Gaussian Graphical
#' Models with the Posterior Predictive Distribution and Bayesian Model Selection. \href{https://psyarxiv.com/yt386}{pre print}
#'
#' @return one object of class \code{ggplot} when \code{type = "global"}. One object for each pairwise contrast when \code{type = "nodewise"}
#'
#' @note This method is Bayesian, as it relies on the posterior predictive distribution. That said, there are clear parallels to frequentist testing-e.g.,
#' assuming group equality and critical regions. Most importantly, this method CANNOT provide evidence for the null hypothesis. Thus it can only reject
#' the underlying assumption of group equality. For gaining (relative) evidence for the null hypothesis see X.
#'
#' @export
#'
#' @importFrom ggridges stat_density_ridges
#'
#' @examples
#' # assume group equality
#' Y1 <- MASS::mvrnorm(500, rep(0, 16), Sigma = diag(16))
#' Y2 <- MASS::mvrnorm(500, rep(0, 16), Sigma = diag(16))
#' Y3 <- MASS::mvrnorm(500, rep(0, 16), Sigma = diag(16))
#'
#'# global
#'ggm_ppc  <- ggm_compare_ppc(Y1, Y2, Y3, type = "global", iter = 50)
#'
#'# plot
#' plot(ggm_ppc, ) +
#' theme_minimal(base_size = 14)  +
#' theme(legend.position = "none",
#' panel.grid.minor.x = element_blank())
#'
#'# nodewise
#'ggm_ppc  <- ggm_compare_ppc(Y1, Y2, Y3, type = "global", iter = 5000)
#'
#' plot(ggm_ppc, log = T)[[1]] +
#' theme_minimal(base_size = 14)  +
#' theme(legend.position = "none",  panel.grid.minor.x = element_blank())
plot.ggm_compare_ppc <- function(x,
                                 critical = 0.05,
                                 col_noncritical = "#84e184A0",
                                 col_critical = "red",
                                 point_size = 2,
                                 log = FALSE){

  # change object name
  fit <- x



  if(fit$type == "global"){

    # number of contrasts
    k <- length(fit$pvalue)

    ppc <- unlist(fit$predictive_risk)

    dat_ppc <- data.frame(ppc = ppc,
                          contrast = rep( gsub(x = fit$names,pattern =  "_", replacement = ""),
                                          each = fit$iter) )

    dat_obs <- data.frame(contrast =  gsub(x = fit$names,pattern =  "_",
                                           replacement = ""),
                          ppc = unlist(fit$obs_jsd))

    if(isFALSE(log)){

      plt <- ggplot(dat_ppc, aes(x = ppc,
                                 y = contrast,
                                 fill = factor(..quantile..))) +

        stat_density_ridges(geom = "density_ridges_gradient",
                            calc_ecdf = TRUE, alpha = 0.5,
                            quantiles = c(0.025, 1 - (critical))) +
        scale_fill_manual( values = c(col_noncritical, col_noncritical, col_critical)) +
        theme(legend.position = "none") +
        xlab("Predictive Check") +
        ylab("Contrast") +
        geom_point(inherit.aes = F,
                   data = dat_obs,
                   aes(x = ppc,
                       y = contrast),
                   size = point_size) +
        scale_y_discrete(limits = rev(levels(dat_obs$contrast)))

    }

    if(isTRUE(log)){

      plt <- ggplot(dat_ppc, aes(x = log(ppc),
                                 y = contrast,
                                 fill = factor(..quantile..))) +
        stat_density_ridges(geom = "density_ridges_gradient",
                            calc_ecdf = TRUE,
                            quantiles = c(0.025, 1 - (critical ))) +
        scale_fill_manual( values = c(col_noncritical, col_noncritical, col_critical)) +
        theme(legend.position = "none") +
        xlab("Predictive Check") +
        ylab("Contrast") +
        geom_point(inherit.aes = F,
                   data = dat_obs,
                   aes(x = log(ppc),
                       y = contrast),
                   size = point_size) +
        scale_y_discrete(limits = rev(levels(dat_obs$contrast)))

    }

  }
  if(fit$type == "nodewise" ){


    if(isTRUE(log)){
      plt <- list()
      k <- length(fit$names)

      for(i in 1:k){
        dat_obs <- data.frame(ppc = unlist(fit$obs_jsd[[i]]),
                              node = 1:fit$info$dat_info$p[1])


        test <- reshape::melt( do.call(rbind, fit$predictive_risk[[i]]))

        test$node <- factor(test$X2, levels = rev(1:fit$info$dat_info$p[1]), labels = rev(1:fit$info$dat_info$p[1]))
        dat_obs$node <- factor(dat_obs$node, levels = rev(1:fit$info$dat_info$p[1]), labels = rev(1:fit$info$dat_info$p[1]))

        plt[[i]] <- ggplot(test, aes(x = log(value),
                                     y = node,
                                     fill = factor(..quantile..))) +
          stat_density_ridges(geom = "density_ridges_gradient",rel_min_height = 0.01,
                              calc_ecdf = TRUE,
                              quantiles = c(0.025, 1 - (critical))) +
          scale_fill_manual( values = c(col_noncritical, col_noncritical, col_critical)) +

          geom_point(data = dat_obs, inherit.aes = F,  aes(x = log(ppc), y = node), size = point_size) +
          theme(legend.position = "none") +
          xlab("Predictive Check") +
          ylab("Node") +
          ggtitle( gsub(x =fit$names[[i]],pattern =  "_",
                        replacement = ""))
      }
    }
    if(isFALSE(log)){

      plt <- list()
      k <- length(fit$names)

      for(i in 1:k){
        dat_obs <- data.frame(ppc = unlist(fit$obs_jsd[[i]]),
                              node = 1:fit$info$dat_info$p[1])


        test <- reshape::melt( do.call(rbind, fit$predictive_risk[[i]]))

        test$node <- factor(test$X2, levels = rev(1:fit$info$dat_info$p[1]), labels = rev(1:fit$info$dat_info$p[1]))
        dat_obs$node <- factor(dat_obs$node, levels = rev(1:fit$info$dat_info$p[1]), labels = rev(1:fit$info$dat_info$p[1]))

        plt[[i]] <- ggplot(test, aes(x = value,
                                     y = node,
                                     fill = factor(..quantile..))) +
          stat_density_ridges(geom = "density_ridges_gradient",
                              calc_ecdf = TRUE,
                              quantiles = c(0.025, 1 - (critical))) +
          scale_fill_manual( values = c(col_noncritical, col_noncritical, col_critical)) +
          # scale_y_discrete(limits = rev(levels(as.factor(test$X2)))) +
          geom_point(data = dat_obs, inherit.aes = F,  aes(x = ppc, y = node)) +
          # scale_y_discrete(limits = rev(levels(as.factor(dat_obs$node)))) +
          theme(legend.position = "none") +
          xlab("Predictive Check") +
          ylab("Node") +
          ggtitle( gsub(x =fit$names[[i]],pattern =  "_",
                        replacement = ""))
      }


    }
  }
  plt
}

