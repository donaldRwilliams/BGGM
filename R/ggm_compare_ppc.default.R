#' Compare GGMs with the Posterior Predictive Distribution
#' @name ggm_compare_ppc.default
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
#'
#' @examples
#'
#' # Assume null is true
#' Y1 <- MASS::mvrnorm(500, rep(0, 16), Sigma = diag(16))
#' Y2 <- MASS::mvrnorm(500, rep(0, 16), Sigma = diag(16))
#' Y3 <- MASS::mvrnorm(500, rep(0, 16), Sigma = diag(16))
#'
#' # global predictive check (iter = 50 for demonstrative purposes)
#' ggm_ppc <- ggm_compare_ppc(Y1, Y2, Y3,
#'                            type = "global", iter = 50)
#' summary(ggm_ppc)
#' @export
ggm_compare_ppc.default <- function(..., type = "global", iter = 5000, cores = 2){

  if(type == "global"){

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
   parallel::clusterExport(cl, varlist = c("Mo_risk_help_node",
                                           "node_jsd_help",
                                           "beta_helper",
                                           "kl_func"))

   doSNOW::registerDoSNOW(cl)

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

class(returned_object) <- "ggm_compare_ppc"
returned_object

}


#' @title S3 estimate method
#' @name ggm_compare_ppc
#' @param ... currently not used
#'
#' @description S3 estimate method
#' @seealso \code{\link{ggm_compare_ppc.default}}
#' @export
ggm_compare_ppc <- function(...) {
  UseMethod("ggm_compare_ppc")
}


#' @name print.ggm_compare_ppc
#' @title  Print method for \code{ggm_compare_ppc.default} objects
#'
#' @param x An object of class \code{ggm_compare_ppc}
#' @param ... currently ignored
#' @export
print.ggm_compare_ppc <- function(x,...){
  cat("BGGM: Bayesian Gaussian Graphical Models \n")
  cat("--- \n")
  if(x$type == "nodewise"){
    cat("Type: GGM Compare (Nodewise Predictive Check) \n")
  } else{
    cat("Type: GGM Compare (Global Predictive Check) \n")
  }
  p <- x$info$dat_info$p[1]
  cat("Posterior Samples:", x$iter, "\n")

  groups <- length(x$info$dat)
  for (i in 1:groups) {
    cat("  Group",
        paste(i, ":", sep = "") ,
        x$info$dat_info$n[[i]],
        "\n")
  }
  # cat("Observations (total):", sum(x$info$dat_info$n), "\n")
  cat("Variables (p):", p, "\n")
  cat("Edges:", .5 * (p * (p-1)), "\n")
  cat("--- \n")
  cat("Call: \n")
  print(x$call)
  cat("--- \n")
  cat("Date:", date(), "\n")
}


#' @name summary.ggm_compare_ppc
#' @title Summary method for \code{ggm_compare_ppc.default} objects
#'
#' @param object An object of class \code{ggm_compare_ppc}
#' @param ... currently ignored
#' @seealso \code{\link{ggm_compare_ppc.default}}
#' @return A list containing the summarized results
#' @export
summary.ggm_compare_ppc <- function(object, ...){

  p <- object$info$dat_info$p[1]

  if (object$type == "global") {
    results <- data.frame(
      contrast = do.call(rbind, object$names),
      KLD =  do.call(rbind, object$obs_jsd),
      p_value = object$pvalue
    )
  }
  if (object$type == "nodewise") {
    results <- list()
    for (i in 1:length(object$obs_jsd)) {
      results[[i]] <-
        data.frame(
          node = 1:p ,
          KLD =  round(do.call(rbind, object$obs_jsd[[i]]), 3),
          p_value = unlist(object$pvalue[[i]])
        )
      names(results)[[i]] <- object$names[[i]]
    }
  }

  returned_object <- list(results = results,
                          object = object)
  class(returned_object) <- "summary.ggm_compare_ppc"
  return(returned_object)
}


#' @name print.summary.ggm_compare_ppc
#' @title Summary method for \code{summary.ggm_compare_ppc} objects
#' @param x An object of class \code{summary.ggm_compare_ppc}
#' @param ... currently ignored
#' @seealso \code{\link{summary.ggm_compare_ppc}}
#' @export
print.summary.ggm_compare_ppc <- function(x, ...){
  cat("BGGM: Bayesian Gaussian Graphical Models \n")
  cat("--- \n")
  if(x$object$type == "nodewise"){
    cat("Type: GGM Comparison (Nodewise Predictive Check) \n")
  } else{
    cat("Type: GGM Comparison (Global Predictive Check) \n")
  }
  p <- x$object$info$dat_info$p[1]
  cat("Posterior Samples:", x$object$iter, "\n")

  groups <- length(x$object$info$dat)
  for (i in 1:groups) {
    cat("  Group",
        paste(i, ":", sep = "") ,
        x$object$info$dat_info$n[[i]],
        "\n")
  }
  cat("Variables (p):", p, "\n")
  cat("Edges:", .5 * (p * (p-1)), "\n")
  cat("--- \n")
  cat("Call: \n")
  print(x$object$call)
  cat("--- \n")
  cat("Estimates: \n \n")
  if(x$object$type == "global"){
    print(x$results, right = T, row.names = F,...)
    cat("--- \n")
    cat("note: \np_value = p(T(Y_rep) > T(y)|Y)\nKLD = (symmetric) Kullback-Leibler divergence")
  }
  if(x$object$type == "nodewise"){
    for(i in 1:length(x$object$obs_jsd)){
      cat(do.call(rbind, x$object$names)[[i]], "\n")
      print(   x$results[[i]],  row.names = F, ...)
      cat("\n")
    }

    cat("--- \n")
    cat("note: \np_value = p(T(Y_rep) > T(y)|Y)\nKLD = (symmetric) Kullback-Leibler divergence")
  }
}
