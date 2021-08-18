ggm_search <- function(x, n = NULL,
                       method = "mc3",
                       prior_prob = 0.1,
                       iter = 5000,
                       stop_early = 1000,
                       bma_mean = TRUE,
                       seed = 1,
                       progress = TRUE, ...){

  set.seed(seed)

  x <- Y

  if (base::isSymmetric(as.matrix(x))) {
    S <- x
  } else {
    S <- cor(x)
    n <- nrow(x)
  }

  p <- ncol(S)

  if(method == "mc3"){

    if(is.null(stop_early)){

      stop_early <- iter

    }

    pcors <- -cov2cor( solve(S) ) + diag(2, p)

    # test full vs missing one edge
    BF_01 <- exp(-0.5 *  (tstat(r = pcors, n = n, p - 2)^2 - log(n)))
    BF_10 <- 1/BF_01

    adj_start <- ifelse(BF_10 > 1, 1, 0)

    old <- hft_algorithm(Sigma = S,
                         adj = adj_start,
                         tol = 1e-10,
                         max_iter = 100)

    bic_old <- bic_fast(Theta = old$Theta,
                        S = S,
                        n = n,
                        prior_prob = prior_prob)

    if(isTRUE(progress)){
      message(paste0("BGGM: Sampling Graphs"))
    }

    fit <- .Call('_BGGM_search',
                 PACKAGE = 'BGGM',
                  S = S,
                  iter = iter,
                  old_bic = bic_old,
                  start_adj = adj_start,
                  n = n,
                  gamma = prior_prob,
                  stop_early = stop_early,
                  progress = progress)

    if(isTRUE(progress)){
      message("BGGM: Finished")
    }

    # accepted
    acc <- fit$acc

    # first matrix (starting values)
    fit$adj[,,1] <- adj_start

    # approximate marginal likelihood
    approx_marg_ll <- fit$bics

    # starting bic
    approx_marg_ll[1] <- bic_old

    if(!is.null(stop_early)){
      approx_marg_ll <- approx_marg_ll[which(approx_marg_ll != 0)]
      fit$adj <- fit$adj[,, which(approx_marg_ll != 0)]
    }

    adj_path <- fit$adj

    selected <- which.min(approx_marg_ll)

    if(acc == 0){
      adj <- fit$adj[,,1]
    } else {
      adj <-  fit$adj[,,selected]
    }

    # BFs vs mpm
    delta <- approx_marg_ll - min(approx_marg_ll)

    probs <- exp(-0.5 * delta) / sum( exp(-0.5 * delta) )

    Theta_map <- hft_algorithm(
      Sigma = S,
      adj = adj,
      tol = 1e-10,
      max_iter = 100
    )

      pcor_adj <- -cov2cor(Theta_map$Theta) + diag(2, p)

  }

  if(bma_mean & acc > 0){

    graph_ids <- which(duplicated(approx_marg_ll) == 0)[-1]

    delta <- (approx_marg_ll[graph_ids] - min(approx_marg_ll[graph_ids])) * (6 / (2*sqrt(2*n)))

    probs <- exp(- 0.5 * delta) / sum(exp(- 0.5 * delta))

    graphs <- adj_path[,,graph_ids]

    Theta_bma <- lapply(1:length(probs), function(x){

      hft_algorithm(Sigma = S,
                    adj = graphs[,,x],
                    tol =1e-10,
                    max_iter = 1000)$Theta * probs[x]
      })

    Theta_bma <- Reduce("+", Theta_bma)
    pcor_bma <- -cov2cor(Theta_bma) + diag(2, p)

  } else {

    Theta_bma <- NULL
    pcor_bma <- NULL

  }

  returned_object <- list(pcor_adj = pcor_adj,
                          Theta_map = Theta_map,
                          Theta_bma = Theta_bma,
                          pcor_bma = pcor_bma,
                          adj = adj,
                          adj_start = adj_start,
                          probs = probs,
                          approx_marg_ll = approx_marg_ll,
                          selected = selected,
                          BF_start = BF_10,
                          adj_path = adj_path,
                          acc = acc,
                          S = S,
                          n = n)

  #rm(.Random.seed, envir=.GlobalEnv)

  class(returned_object) <- c("BGGM",
                              "ggm_search")

  return( returned_object )
}


print_ggm_search <- function(x, ...){

  cat("BGGM: Bayesian Gaussian Graphical Models \n")
  cat("--- \n")

  if(x$acc == 0){
    mat <- x$pcor_adj
    p <- ncol(mat)

    if(is.null( colnames(x$S))){
      colnames(mat) <- 1:p
      row.names(mat) <- 1:p

    } else {
      colnames(mat) <- colnames(x$S)
      row.names(mat) <- colnames(x$S)
    }

    cat("Most Probable Graph:\n\n")
    print(round(mat, 3))

  } else {
    mat <- x$pcor_bma
    p <- ncol(mat)

    if(is.null( colnames(x$S))){
      colnames(mat) <- 1:p
      row.names(mat) <- 1:p

    } else {
      colnames(mat) <- colnames(x$S)
      row.names(mat) <- colnames(x$S)
    }
    cat("Bayesian Model Averaged Graph:\n\n")
    print(round(mat, 3))

  }
}
