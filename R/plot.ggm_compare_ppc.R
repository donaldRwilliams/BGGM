#' Plot GGM Comparison with Predictive DIstribution
#'
#' @param x object of class \code{ggm_compare_ppc}
#' @param critical significance level
#' @param col_noncritical fill of non critical region
#' @param col_critical  fill of critical region (e.g., \code{critical = 0.05})
#' @param point_size point size for the observed KL-divergence
#' @param log log transformation. Useful for small values and skewed predictive distributions
#'
#' @return object of class \code{ggplot}
#' @export
#'
#' @examples
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
                              node = 1:fit$info$dat_info$p[i])


        test <- reshape::melt( do.call(rbind, fit$predictive_risk[[i]]))

        test$node <- factor(test$X2, levels = rev(1:fit$info$dat_info$p[1]), labels = rev(1:fit$info$dat_info$p[1]))
        dat_obs$node <- factor(dat_obs$node, levels = rev(1:fit$info$dat_info$p[1]), labels = rev(1:fit$info$dat_info$p[1]))

        plt[[i]] <- ggplot(test, aes(x = log(value),
                                     y = node,
                                     fill = factor(..quantile..))) +
          stat_density_ridges(geom = "density_ridges_gradient",
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

