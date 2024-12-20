##' Perform Bayesian Graph Search and Optional Model Averaging
##'
##' The `ggm_search` function performs a Bayesian graph search to identify the 
##' most probable graph structure (MAP solution) using the Metropolis-Hastings 
##' algorithm. It also computes an optional Bayesian Model Averaged (BMA) solution 
##' across the graph structures sampled during the search. 
##'
##' This function is ideal for exploring the graph space and obtaining an initial 
##' estimate of the graph structure or adjacency matrix. 
##'
##' To refine the results or compute posterior distributions of graph parameters 
##' (e.g., partial correlations), use the \code{\link{bma_posterior}} function, 
##' which builds on the output of `ggm_search` to account for parameter uncertainty.
##'
##' @return A list containing the MAP graph structure, BMA solution (if specified),
##'         and posterior probabilities of the sampled graphs.
##'
##' @seealso \code{\link{bma_posterior}}
##'  
##' @param x Data, either raw data or covariance matrix 
##' @param n For x = covariance matrix, provide number of observations
##' @param method mc3 defaults to MH sampling
##' @param prior_prob Prior prbability of sparseness. 
##' @param iter Number of iterations
##' #@param burn_in Burn in. Defaults to iter/2
##' @param stop_early Default to 1000. Stop MH algorithm if proposals keep being rejected (stopping by default after 1000 rejections).
##' @param bma_mean Compute Bayesian Model Averaged solution
##' @param seed Set seed. Current default is to set R's random seed.
##' @param progress Show progress bar, defaults to TRUE
##' @param ... Not currently in use
##' @author Donny Williams and Philippe Rast
ggm_search <- function(x, n = NULL,
                       method = "mc3",
                       prior_prob = 0.3,
                       iter = 5000,
                       #burn_in = NULL, 
                       stop_early = 1000,
                       bma_mean = TRUE,
                       seed = NULL,
                       progress = TRUE, ...){

  set.seed(seed)
  ## Random seed unless user provided
  if(!is.null(seed) ) {
    set.seed(seed)
  }

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
 
    ## Create starting adjacency matrix, based on the BF_10
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

    ## Add a burnin unless defined by user
    burn_in = 0 # Drop this line once we found solution to creeping BIC in ggm_search MH algo
    if(is.null(burn_in)) {
      burn_in <- round(iter/2) 
    }
   
    # approximate marginal likelihood
    approx_marg_ll <- fit$bics

    # starting bic
    approx_marg_ll[1] <- bic_old

    if(!is.null(stop_early)){
      approx_marg_ll <- approx_marg_ll[which(approx_marg_ll != 0)]
      fit$adj <- fit$adj[,, which(approx_marg_ll != 0)]
    }

    adj_path <- fit$adj

    # Find position of smallest bic
    selected <-  which.min(approx_marg_ll)

    ## acc are accepted proposals
    if(acc == 0){
      adj <- fit$adj[,,1]
    } else {
      adj <-  fit$adj[,,selected]
    }

    ## BFs vs Most Probably Model (mpm)
    ## Comute delta of all models compared to best model
    delta <- approx_marg_ll - min(approx_marg_ll)

    ## Convert differences in marginal log-likelihoods into posterior probabilities for each model
    probs <- exp(-0.5 * delta) / sum( exp(-0.5 * delta) )

    ## MPM:
    Theta_map <- hft_algorithm(
      Sigma = S,
      adj = adj,
      tol = 1e-10,
      max_iter = 1000
    )

    ## Partial Correlation for the MPM model
    pcor_adj <- -cov2cor(Theta_map$Theta) + diag(2, p)

  }

    ## ## Keep samples after burn-in
    ## valid_indices <- (burn_in + 1):length(approx_marg_ll)  
    ## ## Filter BICs and adjacency matrices
    ## approx_marg_ll <- approx_marg_ll[valid_indices]
    ## adj_path <- adj_path[,, valid_indices]

  
  if(bma_mean & acc > 0){

    graph_ids <-  which(duplicated(approx_marg_ll) == 0)[-1]

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



