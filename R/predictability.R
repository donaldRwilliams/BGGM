#' Nodewise Predicatability
#'
#' @param object object of class \code{estimate}
#' @param cred credible interval width used for the decision rule
#' @param iter iterations used for computing R2
#' @param cores number of cores for parallel computing
#' @param ... currently ignored
#'
#' @return
#' @export
#'
#' @examples
predictability <- function(object,
                           cred = 0.95,
                           iter = 1000,
                           cores = 2,...){

  p <- object$p

  cl <- parallel::makeCluster(cores)

  doSNOW::registerDoSNOW(cl)

  adj <- BGGM::select(object, cred = cred)$adjacency_non_zero

  betas <- BGGM:::inverse_2_beta(object, samples = iter)

  # parallel::clusterExport(cl, varlist = c("R2_ppc"))

  pred <- parallel::parLapply(cl = cl,
                              X = 1:p, function(x) R2_ppc(fit = object,
                                                          betas = betas,
                                                          adj = adj,
                                                          which_one = x, sims = iter))

  class(pred) <- "predictability"
  pred

}




#' Assess Predictability of \code{Predictability} Objects
#'
#' @param ... object(s) of class \code{predictability}
#' @param contrast
#' @param prob
#'
#' @return
#' @export
#'
#' @examples
assess_predictability <-  function(..., contrast = "exhaustive",
                                                  prob = 0.95){


  temp <- list(...)
  cred <- 0.95
  lb <- (1 - cred) / 2
  ub <- 1 - lb
  groups <- length(temp)

  if(groups != 2){
    stop("currently only comparing groups nodewise is possible")
  }

  g1 <- temp[[1]]
  g2 <- temp[[2]]

  p <- length(g1)

  diff <- lapply(1:p, function(x) {temp <- g1[[x]] - g2[[x]];
  pp <- post_prob(temp);
  ci <- quantile(pp, probs = c(lb, ub));
  list(diff = temp, pp = pp, ci = ci)
  })
  class(diff) <- "assess_predictability"
  diff
}

#' Plot \code{assess_predictability} Objects
#'
#' @param x object of class \code{assess_predictability}
#' @param ...
#' @importFrom ggridges stat_density_ridges
#' @return
#' @export
#'
#' @examples
plot.assess_predictability <- function(x,...){

  p <- length(x)

  diffs <- lapply(1:p, function(z)  x[[z]]$diff)

  dat <- reshape2::melt(diffs)
  dat$L1 <- as.factor(dat$L1)

  dat$L1 <- factor(dat$L1, labels = order(tapply(dat$value, dat$L1, mean)),
                   levels = order(tapply(dat$value, dat$L1, mean)))




  plt <- ggplot(dat, aes(x=value, y=as.factor(L1), fill=factor(..quantile..))) +
    geom_vline(xintercept = 0, linetype = 'dotted') +
    stat_density_ridges( rel_min_height = 0.01, geom = "density_ridges_gradient", calc_ecdf = TRUE, quantiles = c(0.025, 0.975)) +
    scale_fill_manual(
      name = "Probability", values = c("#FF0000A0", "#A0A0A0A0", "#0000FFA0"),
      labels = c("(0, 0.025]", "(0.025, 0.975]", "(0.975, 1]")) +
    theme_classic() +
    theme(legend.position = "none",
          panel.grid.major = element_line(color = "grey95"))
  plt
}


