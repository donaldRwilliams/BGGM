#' @title GGMs with Exploratory Bayesian Hypothesis Testing
#' @name explore
#' @description Learn the conditional (in)dependence structure with the Bayes factor computed from the matrix-F prior distribution. It is
#' possible to test for only positive or negative edges, as well as two sided hypothesis testing (which is the customary approach). Further
#' there is also an exhaustive option that provides the posterior probability of the null, greater than zero, and less than zero.
#'
#' @param Y data matrix  (\emph{n} by  \emph{p}).
#' @param iter number of posterior samples.
#' @param cores number of cores for parallel computing. The default is 2, but this can be changed.
#' @param prior_sd hypothesized standard deviation of the prior distribution.
#' @param ... currently not used
#'
#' @return list of class \code{explore}:
#' \itemize{
#' \item \code{parcors_mat} partial correlation matrix
#' \item \code{parcors_sd} partial correlation standard deviations
#' \item \code{samples} list of prior and posterior samples
#' \itemize{
#'  \item \code{fisher_z_post} Fisher z transformed posterior distributions (partial correlations)
#'  \item  \code{pcor_post} partial correlation posterior distributions (not transformed)
#'  \item \code{inv_cov_post} inverse covariance matrix posterior distribution
#'  \item \code{pcor_prior} partial correlation prior distribution
#'  \item \code{fisher_z_prior}  Fisher z transformed prior distributions (partial correlations)
#'  }
#' \item \code{delta} hyperparameter of matrix-F prior distribution (corresponds to \code{prior_sd})
#' \item \code{iter} number of posterior and prior samples
#' \item \code{dat} data matrix
#' \item \code{call} match.call()
#' \item \code{p} number of variables
#' \item \code{cores} number of cores
#' \item \code{edge} number of estimated edges
#' }
#' @note After sampling from the posterior distribution, use \code{select} to determine the edge set and \code{plot} for visualizing the
#' edge set. see \code{methods(class = "explore")}
#' @examples
#' # p = 10
#' Y <- BGGM::bfi[,1:10]
#'
#' # sample posterior
#' fit <- explore(Y, iter = 5000)
#'
#' # select E
#' E <- select(fit, BF_cut = 3)
#'
#' # summarize
#' summary(E)
#'
#' # non-zero edges
#' E$partials_non_zero
#'
#' # adjacency matrix
#' E$Adj_10
#'
#' # null adjacency matrix
#' E$Adj_01
#' @export
explore <- function(Y, prior_sd = 0.25,
                            iter = 5000, cores = 2,...){
  # rename
  X <- Y

  # delta parameter
  delta <- delta_solve(prior_sd)

  # remove NAs
  X <- na.omit(X)

  # number of columns
  p <- ncol(X)

  # matrices for storage
  parcors_mat <- parcors_sd <- matrix(0, p, p)

  # number of edges
  edges <- 0.5 * (p * (p -1))

  # sample from posterior
  samples <- sampling(X, nu = 10000,
                      delta = delta, n_samples = iter,
                      cores = cores)

  # extract posterior samples (each "chain")
  posterior_samples <- do.call(rbind.data.frame,
                               lapply(1:cores, function(x)
                                 samples[[x]]$pcor_post))

  # posterior mean
  parcors_mat[upper.tri(parcors_mat)] <- colMeans(posterior_samples)[1:edges]
  pacors_mat <- symmteric_mat(parcors_mat)

  # posterior sd
  parcors_sd[upper.tri(parcors_sd)] <- apply(posterior_samples, 2, sd)[1:edges]
  pacors_sd <- symmteric_mat(parcors_sd)

  # returned object
  returned_object <- list(parcors_mat = pacors_mat,
                          parcors_sd = parcors_sd,
                          samples = samples,
                          delta = delta,
                          iter = iter,
                          dat = X,
                          call = match.call(),
                          p = p,
                          cores = cores,
                          edge = edges)

  class(returned_object) <- "explore"
  return(returned_object)

}

#' @name print.explore
#' @title  Print method for \code{explore} objects
#' @param x An object of class \code{explore}
#' @param ... currently ignored
#' @seealso \code{\link{explore}}
#' @export
print.explore <- function(x,...){
  cat("BGGM: Bayesian Gaussian Graphical Models \n")
  cat("--- \n")
  cat("Type: Hypothesis Testing (Exploratory) \n")
  cat("Posterior Samples:", x$iter, "\n")
  cat("Observations (n):", nrow(x$dat), "\n")
  cat("Variables (p):", x$p, "\n")
  cat("Edges:", .5 * (x$p * (x$p-1)), "\n")
  cat("Delta:", x$delta, "\n")
  cat("--- \n")
  cat("Call: \n")
  print(x$call)
  cat("--- \n")
  cat("Date:", date(), "\n")
}


#' @name summary.explore
#' @title  Summary method for \code{explore} objects
#' @param object An object of class \code{explore}
#' @seealso \code{\link{explore}}
#' @param ... currently ignored
#' @export
summary.explore <- function(object,...){
dat_results <- select(object,
                BF_cut = 0,
                alternative = "exhaustive")

returned_object <- list(dat_results = dat_results,
                        object = object)
class(returned_object) <- "summary.explore"
return(returned_object)
}

#' @title Summary method for \code{summary.explore} objects
#' @name print.summary.explore
#'
#' @param x An object of class \code{summary.explore}
#' @param ... currently ignored
#' @seealso \code{\link{summary.explore}}
#' @export
print.summary.explore <- function(x,...){
   summary(x$dat_results, summarize = TRUE)
}

#' Plot \code{summary.explore}
#'
#' @param x an object of class \code{summary.explore}
#' @param ... currently ignored
#'
#' @return an object of class \code{ggplot}
#' @export

plot.summary.explore <- function(x,...){

  not_h0 <-  1 - (x$dat_results$post_prob$prob_zero)
  h0 <- (x$dat_results$post_prob$prob_zero)

  test_dat <- data.frame(Edge = x$dat_results$post_prob$edge, probability = c(not_h0 - 0.5))

  dat_temp <- test_dat[order(test_dat$probability, decreasing = T),]
  dat_temp$Edge <-
    factor(dat_temp$Edge,
           levels = rev(dat_temp$Edge),
           labels = rev(dat_temp$Edge))


  ggplot(dat_temp, aes(x = Edge, y = probability)) +
    geom_linerange(aes(x = Edge, ymin = 0, ymax = probability),
                   position = position_dodge(width = 1)) +

    scale_y_continuous(limits = c(-0.5, 0.5),
                       labels = c(1, 0.5, 0, 0.5, 1)) +
    geom_point(aes(x = Edge, y = probability),
               position = position_dodge(width = 1)) +

    ylab(
      expression(
        italic(H)[0] * symbol(' \254 ') * "Posterior Probability " * symbol('\256 ') *
          "'not " * italic(H)[0] * "'"
      )
    ) +
    geom_point(color = "blue", size = 1) +
    geom_hline(yintercept = 0, linetype = "dotted") +
    coord_flip()
}
