#' @title Constrained Posterior Distribution
#'
#' @description Compute the posterior distribution
#'              with off-diagonal elements of the precision matrix constrained
#'              to zero.
#'
#' @param object An object of class \code{estimate} or \code{explore}
#'
#' @param adj A \code{p} by \code{p} adjacency matrix. The zero entries denote the
#'            elements that should be constrained to zero.
#'
#' @param method Character string. Which method should be used ? Defaults to
#' the "direct sampler" (i.e., \code{method = "direct"}) described in
#' \insertCite{@page 122, section 2.4,  @lenkoski2013direct;textual}{BGGM}. The other
#' option is a Metropolis-Hastings algorithm (\code{MH}).
#' See details.
#'
#' @param iter Number of iterations (posterior samples; defaults to 5000).
#'
#'
#' @param progress Logical. Should a progress bar be included (defaults to \code{TRUE}) ?
#'
#' @param ... Currently ignored.
#'
#' @references
#' \insertAllCited{}
#'
#'
#' @return An object of class \code{contrained}, including
#'
#' \itemize{
#'
#' \item \code{precision_mean} The posterior mean for the precision matrix.
#'
#' \item \code{pcor_mean} The posterior mean for the precision matrix.
#'
#' \item \code{precision_samps} A 3d array of dimension \code{p} by \code{p} by \code{iter}
#'                              including the sampled precision matrices.
#'
#'  \item \code{pcor_samps} A 3d array of dimension \code{p} by \code{p} by \code{iter}
#'                              including sampled partial correlations matrices.
#' }
#'
#'
#' @examples
#' \donttest{
#'
#' # data
#' Y <- bfi[,1:10]
#'
#' # sample posterior
#' fit <- estimate(Y, iter = 100)
#'
#' # select graph
#' sel <- select(fit)
#'
#' # constrained posterior
#' post <- constrained_posterior(object = fit,
#'                               adj = sel$adj,
#'                               iter = 100,
#'                               progress = FALSE)
#'
#' }
#' @export
constrained_posterior <- function(object,
                                  adj,
                                  method = "direct",
                                  iter = 5000,
                                  progress = TRUE,
                                  ...){

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

  if(progress){

    message(paste0("BGGM: Constrained Posterior"))

    }

  if(method == "direct"){

    # covert to correlations
    cors <- pcor_to_cor(object, iter = iter)

    diag(adj) <- 1

   Theta_samps <- .Call("_BGGM_contrained_helper",
                        cors = cors$R,
                        adj = adj,
                        iter = iter,
                        progress = progress)

   mh_object <- NULL

  }

  if(method == "MH"){

    prior_sd <- object$prior_sd

    Y <- object$Y

    p <- object$p

    n <- object$n

    prec <- precision(object, progress = FALSE)

    iter <- object$iter

    samples <- matrix(prec$precision[ ,, 1:iter][lower.tri(diag(p), diag = TRUE)],
                      nrow = iter,
                      ncol = p*(p-1)*0.5 + p,
                      byrow = TRUE)

    samps1 <- rw_helper(adj = adj,
                        means_post = colMeans(samples),
                        cov_post = cov(samples))


    # burning in: 500; thin: 20
    iter_comb <- iter * 20 + 500

    mh_object <- .Call("_BGGM_fast_g_matrix_F",
                         PACKAGE = "BGGM",
                         Y = Y,
                         adj = adj,
                         mu_samples = samps1$mean_post,
                         cov_samples = samps1$cov_post,
                         iter = iter * 20 + 500,
                         p = p,
                         N = n,
                         prior_sd = prior_sd,
                         kappa1 = 0.05,
                         progrss = TRUE)

    Theta_samps <- mh_object$Theta_G[,,seq(501, iter_comb, 20)]
    mh_object$Theta_G <- NULL
    }

  if(isTRUE(progress)){

    message("BGGM: Finished")

  }


  # partial correlations
  pcor_samps <- vapply(1:iter, function(x)
    -cov2cor(Theta_samps[,,x]) + diag(2, p),
    FUN.VALUE =  matrix(0, p, p)
  )

  Theta_mu <- mean_array(Theta_samps)
  pcor_mu <-  mean_array(pcor_samps)

  colnames(pcor_mu) <- colnames(object$Y)
  row.names(pcor_mu) <- colnames( object$Y )

  returned_object <- list(
    precision_mean = Theta_mu,
    pcor_mean = pcor_mu,
    precision_samps = Theta_samps,
    pcor_samps = pcor_samps,
    mh_object = mh_object
  )

  class(returned_object) <- c("BGGM",
                              "constrained")

  return(returned_object)

}

print_constrained <- function(x, ...){
  cat("BGGM: Bayesian Gaussian Graphical Models\n")
  cat("Constrained posterior\n")
  cat("---\n")
  cat("Estimates: \n\n")
  print(round(x$pcor_mean, 3))
}



rw_helper <- function(adj, means_post, cov_post){
  diag(adj) <- 1
  p <- ncol(adj)
  selectlower <- lower.tri(diag(p), diag = TRUE)
  which0 <- which(adj[selectlower]==0)
  which1 <- which(adj[selectlower]==1)
  means1_post <- c(means_post[which1] +
                     cov_post[which1,which0] %*%
                     solve(cov_post[which0,which0]) %*%
                     (-means_post[which0]))

  cov1_post <- cov_post[which1,which1] -
    cov_post[which1,which0] %*%
    solve(cov_post[which0,which0]) %*%
    cov_post[which0,which1]

  list(mean_post = means1_post, cov_post = cov1_post)
}
