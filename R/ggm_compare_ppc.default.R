#' Compare GGMs with the Posterior Predictive Distribution
#'
#' @description Compare GGMs with the posterior predictive distribution. The method assume group equality, and the predictive check
#' allows for testing whether that assumption should be modified--i.e., the GGMs are actually different. The current test statistic avaiable is
#' Kullback–Leibler divergence, which in this case, can be understood as a likelihood ratio for multiariate normal distributions. There are
#' two options: (1) 'global' and (2) 'nodewise.' The former tests the entire GGM, whereas the latter allows for testing specific nodes (variables)
#' in the model.
#'
#' @param ... data matrices (\emph{n} \times  \emph{p}). Requires at least two.
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
#' @export
#'
#' @note This method is Bayesian, as it relies on the posterior predictive distribution.
#' That said, there are clear parallels to frequentist testing-e.g., assuming group
#' equality and critical regions. Most importanly, this method CANNOT provide evidence
#' for the null hypothesis. Thus it can only reject the unerlying assumption of group equality.
#' For gaining (relative) evidence for the null hypothesis see \link[BGGM:ggm_compare_bf.default]{ggm_compare_bf}.
#'
#' see methods(class = "ggm_compare_ppc")
#'
#' @references
#' Williams, D. R., Rast, P., Pericchi, L. R., & Mulder, J. (2019). Comparing Gaussian Graphical
#' Models with the Posterior Predictive Distribution and Bayesian Model Selection. \href{https://psyarxiv.com/yt386}{pre print}
#'
#' @examples
#'
#' Assume null is true
#' Y1 <- MASS::mvrnorm(500, rep(0, 16), Sigma = diag(16))
#' Y2 <- MASS::mvrnorm(500, rep(0, 16), Sigma = diag(16))
#' Y3 <- MASS::mvrnorm(500, rep(0, 16), Sigma = diag(16))
#'
#' # global predictive check
#' ggm_ppc <- ggm_compare_ppc(Y1, Y2, Y3,
#'                            type = "global", iter = 5000)
#' summary(ggm_ppc)
#'
#' # nodewise predictive check
#' ggm_ppc <- ggm_compare_ppc(Y1, Y2,
#'                            type = "nodewise", iter = 5000)
#' summary(ggm_ppc)
#'
ggm_compare_ppc.default <- function(..., type = "global", iter = 5000, cores = 2){

  if(type == "global"){

    info <- BGGM:::Y_combine(...)

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

    nms <- list()

    for(i in 1:nrow(info$pairwise)){

      n1 <- info$dat_info$n[info$pairwise[i,1]]

      n2 <- info$dat_info$n[info$pairwise[i,2]]

      ppc  <-  list(unlist(lapply(1:iter, function(x) BGGM:::Mo_risk_help(x, post, n1, n2, p))))



      predictive_risk[[i]] <- ppc

      y1 <- info$dat[info$pairwise[i,1]][[1]]

      y2 <- info$dat[info$pairwise[i,2]][[1]]

      # y1 <- cov2cor(y1)
      # y2 <- cov2cor(y2)
      obs <- 0.5 * BGGM:::KL(BGGM:::unbiased_cov(y1), BGGM:::unbiased_cov(y2)) +
             0.5 * BGGM:::KL(BGGM:::unbiased_cov(y2), BGGM:::unbiased_cov(y1))

      obs_jsd[[i]] <- obs

      nms[[i]] <- paste("Y_g", info$pairwise[i,], sep = "", collapse = " vs ")

      }

    results_ppc <- do.call(cbind.data.frame, predictive_risk)

    pvalues <- unlist(lapply(1:nrow(info$pairwise), function(x)  mean(results_ppc[,x]  > obs_jsd[[x]])))

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

   info <- BGGM:::Y_combine(...)

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

   doSNOW::registerDoSNOW(cl)

   nms <- list()

    for(i in 1:nrow(info$pairwise)){

     n1 <- info$dat_info$n[info$pairwise[i,1]]

     n2 <- info$dat_info$n[info$pairwise[i,2]]

     ppc <- parallel::parLapply(cl=cl, X = 1:iter,  function(x)
                            BGGM:::Mo_risk_help_node(x, post, n1, n2, p))

     predictive_risk[[i]] <- ppc
   }


   for(i in 1:nrow(info$pairwise)){

     y1 <- info$dat[info$pairwise[i,1]][[1]]

     y2 <- info$dat[info$pairwise[i,2]][[1]]

     obs <- lapply(1:info$dat_info$p[1], function(x) node_jsd_help(x, y1, y2))

     nms[[i]] <- paste("Y_g", info$pairwise[i,], sep = "", collapse = " vs ")

     obs_jsd[[i]] <- obs

   }

   pvalue <- list()

   for(i in 1:nrow(info$pairwise)){

     ppc_i <- do.call(rbind, predictive_risk[[i]])

     obs_i <- obs_jsd[[i]]

     pvalues <- lapply(1:info$dat_info$p[1], function(x)  mean(ppc_i[,x]  > obs_i[x]))

     pvalue[[i]] <-  pvalues

   }

 returned_object <-  list(pvalue = pvalue,
                          obs_jsd = obs_jsd,
                          predictive_risk = predictive_risk,
                          info = info,
                          names = nms,
                          iter = iter,
                          type = type,
                          call = match.call())
}

class(returned_object) <- "ggm_compare_ppc"
returned_object

}
