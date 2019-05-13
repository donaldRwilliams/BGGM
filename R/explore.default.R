#' Title
#'
#' @param X
#' @param delta
#' @param iter
#' @param cores
#'
#' @return
#' @export
#'
#' @examples
explore.default <- function(X, prior_sd, iter = 5000,  cores = 2){

  delta <- BGGM:::delta_solve(prior_sd)
  X <- na.omit(X)
  p <- ncol(X)
  parcors_mat <- parcors_sd <- matrix(0, p, p)
  edges <- 0.5 * (p * (p -1))

  samples <- BGGM:::sampling(X, nu = 10000,
                             delta = delta,
                             n_samples = iter,
                             cores = cores)




  posterior_samples <- do.call(rbind.data.frame,
                               lapply(1:cores, function(x)  samples[[x]]$pcor_post))


  parcors_mat[upper.tri(parcors_mat)] <- colMeans(posterior_samples)[1:edges]
  pacors_mat <- BGGM:::symmteric_mat(parcors_mat)

  parcors_sd[upper.tri(parcors_sd)] <- apply(posterior_samples, 2,sd)[1:edges]

  pacors_sd <- BGGM:::symmteric_mat(parcors_sd)


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



