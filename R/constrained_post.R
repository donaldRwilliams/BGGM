constrained_posterior <- function(object,
                                  adj,
                                  iter = 5000,
                                  progress = TRUE){

  if (!any(class(object) %in% c("estimate", "explore"))) {
    stop("object must be of class 'estimate' or 'explore'")
  }

  if (object$iter < iter) {
    stop("iter exceeds iter in the object")
  }

  # ensure diagonal is 1
  diag(adj) <- 1

  # number of nodes
  p <- object$p

  # covert to correlations
  cors <- pcor_to_cor(object, iter = iter)

  if(progress){
    pb <- txtProgressBar(min = 0, max = iter, style = 3)

  }

  # precision matrix
  Theta_samps <- vapply(1:iter, function(s) {

    Theta <- BGGM:::hft_algorithm(Sigma = cors$R[,,s],
                  adj = adj, tol = 0.00001, max_iter = 10)$Theta

    if(progress){
      utils::setTxtProgressBar(pb, s)
    }

    Theta

    }, FUN.VALUE =  matrix(0,p,p)

    )

  # partial correlations
  pcor_samps <- vapply(1:iter, function(x)
     -cov2cor(samps[,,x]) + diag(2, p),
     FUN.VALUE =  matrix(0, p, p)
     )

   Theta_mu <- apply(Theta_samps, 1:2, mean)
   pcor_mu <- apply(pcor_samps, 1:2, mean)

   colnames(pcor_mu) <- colnames( object$Y )
   row.names(pcor_mu) <- colnames( object$Y )
   returned_object <- list(
     Theta_mu = Theta_mu,
     pcor_mu = pcor_mu,
     Theta_samps = Theta_samps,
     pcor_samps = pcor_samps
   )

   class(returned_object) <- c("BGGM",
                               "constrained")

   returned_object

  }


print_constrained <- function(x, ...){
  cat("BGGM: Bayesian Gaussian Graphical Models\n")
  cat("Constrained posterior\n")
  cat("---\n")
  cat("Estimates: \n\n")
  print(round(x$pcor_mu, 3))
}

