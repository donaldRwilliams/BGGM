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
explore.default <- function(X, delta, iter = 5000,  cores = 2){

  X <- na.omit(X)
  p <- ncol(X)
  parcors_mat <- matrix(0, p, p)
  edges <- 0.5 * (p * (p -1))

  samples <- BGGM:::sampling(X, nu = 10000,
                             delta = delta,
                             n_samples = iter,
                             cores = cores)




  posterior_samples <- do.call(rbind.data.frame,
                               lapply(1:cores, function(x)  samples[[x]]$pcor_post))


  parcors_mat[upper.tri(parcors_mat)] <- colMeans(posterior_samples)[1:edges]
  pacors_mat <- BGGM:::symmteric_mat(parcors_mat)



  returned_object <- list(parcors_mat = pacors_mat,
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



