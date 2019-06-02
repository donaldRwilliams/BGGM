#' Compare GGMs with Posterior Predictive KL-Divergence
#'
#' @param ... data matrices (\emph{n} \times  \emph{p}). Requires at least two.
#' @param type \code{type = "global"} for testing the entire precision matrix. \code{type = "nodewise"} for testing each node (i.e., variable)
#' @param iter number of replicated data sets. default is \code{iter = 1000}.
#' @param cores Number of cores for parallel computing. The defaul is 2, but this can be adjusted.
#'
#' @return object of class \code{ggm_compare_ppc}
#' @export
#'
#' @examples
#'
#' # Assume null is true
#' Y1 <- MASS::mvrnorm(500, rep(0, 16), Sigma = BGGM::ptsd_cor1)
#' Y2 <- MASS::mvrnorm(500, rep(0, 16), Sigma = BGGM::ptsd_cor1)
#'
#'
#' ggm_ppc <- ggm_compare_ppc(Y1, Y2, iter = 1000)
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
