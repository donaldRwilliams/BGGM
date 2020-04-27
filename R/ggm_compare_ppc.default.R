#' GGM Compare: Posterior Predictive Distribution
#'
#' @name ggm_compare_ppc
#'
#' @description Compare GGMs with the posterior predictive distribution. The method assume group equality, and the predictive check
#' allows for testing whether that assumption should be modified--i.e., the GGMs are actually different. The current test statistic available is
#' Kullback-Leibler divergence, which in this case, can be understood as a likelihood ratio for multivariate normal distributions. There are
#' two options: (1) 'global' and (2) 'nodewise.' The former tests the entire GGM, whereas the latter allows for testing specific nodes (variables)
#' in the model.
#'
#' @param ... data matrices (\emph{n} by  \emph{p}). Requires at least two.
#' @param type \code{type = "global"} for testing the entire precision matrix. \code{type = "nodewise"} for testing each node (i.e., variable)
#' @param iter number of replicated data sets. default is \code{iter = 1000}.
#' @param cores number of cores for parallel computing. The default is 2, but this can be adjusted.
#'
#'
#' @return list of class \code{ggm_compare_ppc}:
#'
#' \itemize{
#' \item \code{pvalue} posterior predictive p-values
#' \item \code{obs_jsd} observed symmetric KL divergence (Jensen-Shannon divergence)
#' \item \code{predictive_risk} list of predictive distributions
#' \item  \code{info} list of information about the data matrices
#'
#' \itemize{
#' \item \code{dat} list containing the data matrices
#' \item \code{dat_info} sample size for each data matrix
#' \item \code{pairwise} matrix of pairwise combinations
#' }
#'
#' \item \code{names} contrast names (e.g.., Yg1  vs Yg2)
#' \item \code{iter}  number of posterior samples
#' \item \code{type} "global"
#' \item \code{call} match.call()
#' }
#'
#'

#'
#' @note This method is Bayesian, as it relies on the posterior predictive distribution.
#' That said, there are clear parallels to frequentist testing-e.g., assuming group
#' equality and critical regions. Most importantly, this method CANNOT provide evidence
#' for the null hypothesis. Thus it can only reject the underlying assumption of group equality.
#' For gaining (relative) evidence for the null hypothesis see..
#'
#' see methods(class = "ggm_compare_ppc")
#'
#' @examples
#' \donttest{
#' # Assume null is true
#' Y1 <- MASS::mvrnorm(500, rep(0, 16), Sigma = diag(16))
#' Y2 <- MASS::mvrnorm(500, rep(0, 16), Sigma = diag(16))
#' Y3 <- MASS::mvrnorm(500, rep(0, 16), Sigma = diag(16))
#'
#' # global predictive check (iter = 50 for demonstrative purposes)
#' ggm_ppc <- ggm_compare_ppc(Y1, Y2, Y3,
#'                            type = "global", iter = 50)
#' summary(ggm_ppc)
#'
#' plot(ggm_ppc)
#'
#' # nodewise
#' ggm_ppc  <- ggm_compare_ppc(Y1, Y2, Y3, type = "nodewise", iter = 50)
#'
#' plot(ggm_ppc, log = TRUE)
#' }
#' @export
ggm_compare_ppc <- function(..., type = "global", iter = 5000, cores = 1){

  if(type == "global"){

    info <- Y_combine(...)

    groups <- length(info$dat)

    if(groups < 2){
      stop("must have (at least) two groups")
    }

    Y_G <- rbind(info$dat)

    n_total <- sum(info$dat_info$n)

    Y_G <- scale(do.call(rbind, info$dat), scale = T)

    # inverse scatter matrix
    S_G <- solve(t(Y_G) %*% Y_G)

    # M_0 posterior (group equality)
    post <- rWishart(iter, n_total - 1, S_G)

    p <- info$dat_info$p[1]

    predictive_risk <- list()

    obs_jsd <- list()

    nms <- list()

    for(i in 1:nrow(info$pairwise)){

      n1 <- info$dat_info$n[info$pairwise[i,1]]

      n2 <- info$dat_info$n[info$pairwise[i,2]]

      ppc  <-  list(unlist(lapply(1:iter, function(x) Mo_risk_help(x, post, n1, n2, p))))

      predictive_risk[[i]] <- ppc

      y1 <- info$dat[info$pairwise[i,1]][[1]]

      y2 <- info$dat[info$pairwise[i,2]][[1]]

      obs <- 0.5 * KL(unbiased_cov(y1), unbiased_cov(y2)) +
             0.5 * KL(unbiased_cov(y2), unbiased_cov(y1))

      obs_jsd[[i]] <- obs

      nms[[i]] <- paste("Y_g", info$pairwise[i,], sep = "", collapse = " vs ")

      }

    results_ppc <- do.call(cbind.data.frame, predictive_risk)

    pvalues <- unlist(lapply(1:nrow(info$pairwise), function(x)  mean(na.omit(results_ppc[,x])  > obs_jsd[[x]])))

    returned_object <- list(pvalue = pvalues,
                            obs_jsd = obs_jsd,
                            predictive_risk = predictive_risk,
                            info = info,
                            names = nms,
                            iter = iter,
                            type = type,
                            call = match.call())
    }


 if(type == "nodewise"){

   info <- Y_combine(...)

   groups <- length(info$dat)

   if(groups < 2){
     stop("must have two groups")
   }


   Y_G <- rbind(info$dat)

   n_total <- sum(info$dat_info$n)

   Y_G <- scale(do.call(rbind, info$dat), scale = T)

   # inverse scatter matrix
   S_G <- solve(t(Y_G) %*% Y_G)

   # M_0 posterior (group equality)
   post <- rWishart(iter, n_total - 1, S_G)

   p <- info$dat_info$p[1]

   predictive_risk <- list()

   obs_jsd <- list()

   cl <- parallel::makeCluster(cores)

   doParallel::registerDoParallel(cl)

   nms <- list()

    for(i in 1:nrow(info$pairwise)){

     n1 <- info$dat_info$n[info$pairwise[i,1]]

     n2 <- info$dat_info$n[info$pairwise[i,2]]

     ppc <- parallel::parLapply(cl=cl, X = 1:iter,  function(x)
                                Mo_risk_help_node(x, post, n1, n2, p))

     predictive_risk[[i]] <- ppc
   }


   for(i in 1:nrow(info$pairwise)){

     y1 <- info$dat[info$pairwise[i,1]][[1]]

     y2 <- info$dat[info$pairwise[i,2]][[1]]

     obs <- lapply(1:info$dat_info$p[1], function(x) node_jsd_help(x, y1, y2))

     nms[[i]] <- paste("Y_g", info$pairwise[i,], sep = "", collapse = " vs ")

     obs_jsd[[i]] <- obs

   }
   parallel::stopCluster(cl)
   pvalue <- list()

   for(i in 1:nrow(info$pairwise)){

     ppc_i <- do.call(rbind, predictive_risk[[i]])

     obs_i <- obs_jsd[[i]]

     pvalues <- lapply(1:info$dat_info$p[1], function(x)  mean(ppc_i[,x]  > obs_i[x]))

     pvalue[[i]] <-  pvalues

   }
 # parallel::stopCluster(cl)
 returned_object <-  list(pvalue = pvalue,
                          obs_jsd = obs_jsd,
                          predictive_risk = predictive_risk,
                          info = info,
                          names = nms,
                          iter = iter,
                          type = type,
                          call = match.call())
}

class(returned_object) <- c("BGGM", "estimate", "ggm_compare_ppc")
returned_object

}




#' Plot \code{ggm_compare_ppc} Objects
#'
#' @description Plot \href{https://CRAN.R-project.org/package=ggridges/vignettes/introduction.html}{ggridges} for
#' the GGM comparison with posterior predictive KL-divergence. The plots contain the predictive distribution, assuming group equality, as well as
#' the observed KL-divergence. Further, the predictive distributions are conveniently colored to infer whether the null of group equality
#' should be rejected. This is accomplished by having the critical region, corresponding to a desired 'significance' level, shaded in red.
#' Thus, if the observed value is in the red region, this suggests the null hypothesis of group equality should be rejected.
#'
#' @param x object of class \code{ggm_compare_ppc}
#' @param critical 'significance' level
#' @param col_noncritical fill color of the non critical region
#' @param col_critical  fill color of the critical region (e.g., \code{critical = 0.05})
#' @param point_size point size for the observed KL-divergence
#' @param log log transformation. useful for small values and skewed predictive distributions
#' @param ... currently ignored
#' @references
#' Williams, D. R., Rast, P., Pericchi, L. R., & Mulder, J. (2019). Comparing Gaussian Graphical
#' Models with the Posterior Predictive Distribution and Bayesian Model Selection. \href{https://psyarxiv.com/yt386}{pre print}
#'
#' @return one object of class \code{ggplot} when \code{type = "global"}. One object for each pairwise contrast when \code{type = "nodewise"}
#'
#' @note This method is Bayesian, as it relies on the posterior predictive distribution. That said, there are clear parallels to frequentist testing-e.g.,
#' assuming group equality and critical regions. Most importantly, this method CANNOT provide evidence for the null hypothesis. Thus it can only reject
#' the underlying assumption of group equality..
#'

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
#' ggm_ppc  <- ggm_compare_ppc(Y1, Y2, Y3, type = "global", iter = 50)
#'
#' # plot
#' plot(ggm_ppc)
#'
#'\donttest{
#'# nodewise
#'ggm_ppc  <- ggm_compare_ppc(Y1, Y2, Y3, type = "nodewise", iter = 50)
#'
#' plot(ggm_ppc, log = TRUE)
#' }
#' @export
plot.ggm_compare_ppc <- function(x,
                                 critical = 0.05,
                                 col_noncritical = "#84e184A0",
                                 col_critical = "red",
                                 point_size = 2,
                                 log = FALSE,...){

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

        test$node <- factor(test$X2,
                            levels = rev(1:fit$info$dat_info$p[1]),
                            labels = rev(1:fit$info$dat_info$p[1]))

        dat_obs$node <- factor(dat_obs$node, levels = rev(1:fit$info$dat_info$p[1]),
                               labels = rev(1:fit$info$dat_info$p[1]))

        plt[[i]] <- ggplot(test, aes(x = value,
                                     y = node,
                                     fill = factor(..quantile..))) +
          stat_density_ridges(geom = "density_ridges_gradient",
                              calc_ecdf = TRUE,
                              quantiles = c(0.025, 1 - (critical))) +
          scale_fill_manual( values = c(col_noncritical, col_noncritical, col_critical)) +
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

