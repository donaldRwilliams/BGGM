#' @title Compute Posterior Distributions from Graph Search Results
#'
#' The `bma_posterior` function samples posterior distributions of graph 
#' parameters (e.g., partial correlations or precision matrices) based on the 
#' graph structures sampled during a Bayesian graph search performed by 
#' \code{\link{ggm_search}}.
#'
#' This function incorporates uncertainty in both graph structure and parameter 
#' estimation, providing Bayesian Model Averaged (BMA) parameter estimates. 
#'
#' Use `bma_posterior` when detailed posterior inference on graph parameters 
#' is needed, or to refine results obtained from `ggm_search`. 
#'
#' @return A list containing posterior samples and the Bayesian Model Averaged 
#'         parameter estimates.
#'
#' @seealso \code{\link{ggm_search}}
#' @param object 
#' @param param 
#' @param iter 
#' @param progress 
bma_posterior <- function(object,
                          param = "pcor",
                          iter = 5000,
                          progress = TRUE){

  if(!is(object, "ggm_search")){
    stop("class not supported. Must be 'ggm_search'")
  }

  n <- object$n

  p <- ncol(object$adj)

  I_p <- diag(p)

  scatter <- object$S * (n - 1)

  approx_marg_ll <- object$approx_marg_ll

  graphs <- object$adj_path[,,which(duplicated(approx_marg_ll) == 0)]

  ## check if there are more than one graph
  if( length(which(duplicated(approx_marg_ll) == 0)) > 1 ) {
    graphs <- graphs[,,-1]
    graphs_n <- dim(graphs)[3]
    probs <- object$probs
  } else {
    graphs_n <- 1
    probs <- object$probs[1]
  }

  

  if(isTRUE(progress)){
    pb <- utils::txtProgressBar(min = 0, max = iter, style = 3)
  }


  samples <- vapply(1:iter, function(s){

    graph_s <- sample(1:graphs_n, 1,replace = FALSE, probs)

    Sigma <- solve(rWishart(1, n + p - 1, solve(scatter + I_p * p))[,,1])

    ## Check dimension of graphs
    if (length(dim(graphs)) == 3) {
      selected_graph <- graphs[, , graph_s]
    } else if (length(dim(graphs)) == 2) {
      selected_graph <- graphs
    } 
    
    Theta <- hft_algorithm(Sigma =  Sigma,
                          selected_graph,
                           tol = 0.0001,
                           max_iter = 10)$Theta

    if(isTRUE(progress)){
      utils::setTxtProgressBar(pb, s)
    }

    if(param == "pcor"){

      -cov2cor( Theta ) + diag(2, p)

    } else {


      Theta
    }

  }, FUN.VALUE = matrix(0, p, p))

  bma_mean <- apply(samples, 1:2, mean)

  if( is.null( colnames(scatter) ) ){
    colnames(bma_mean) <- 1:p
    row.names(bma_mean ) <- 1:p
  } else {
    colnames(bma_mean) <- colnames(scatter)
    row.names(bma_mean) <- colnames(scatter)
  }

  returned_object <- list(bma_mean = bma_mean,
                          samples = samples)

  class(returned_object) <- c("BGGM",
                              "bma_posterior")

  return(returned_object)

}

print_bma <- function(x,...){
  cat("BGGM: Bayesian Gaussian Graphical Models \n")
  cat("--- \n")
  cat("Bayesian Model Averaged Graph:\n\n")
  print(round(x$bma_mean, 3))
}
