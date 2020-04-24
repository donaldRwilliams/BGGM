#' Extract posterior samples
#'
#' @param x A R object of class \code{estimate} or \code{explore}
#' @param ... currently ignored.
#'
#' @return A matrix of posterior samples (corresponding to the upper-triangular elements)
#'
#' @examples
#' Y <- BGGM::ptsd
#' fit <- estimate(Y)
#' posterior_samples(fit)
#' @export
posterior_samples <- function(x, ...){

  if(is(x, class2 = "explore")){

    if(!is(x, class2 = "default")){
      stop("must be of class default (the object from explore(Y))")
      }

    p <- x$p
    posterior_samples <- do.call(rbind.data.frame,
                                 lapply(1:fit$cores, function(z)
                                  fit$samples[[z]]$fisher_z_post))[1:x$edge]

    posterior_samples <- apply(posterior_samples, 2, BGGM:::fisher_r2z)
    mat_name_num <- sapply(1:p, function(x) paste0(1:p, "--", x, sep = ""))
    colnames(posterior_samples) <-  mat_name_num[upper.tri(mat_name_num)]
    # end explore
  } else if(is(x, "estimate")){

    if(!is(x, class2 = "default")){
      stop("must be of class default (the object from explore(Y))")
    }


    p <- x$p
    posterior_samples <- x$posterior_samples[,  (p^2+1):(p^2+p^2) ]

    mat_name_num <- sapply(1:p, function(x) paste0(1:p, "--", x, sep = ""))
    colnames(posterior_samples) <-  matrix(mat_name_num)
    posterior_samples <- posterior_samples[,mat_name_num[upper.tri(mat_name_num)]]
    } else {
      stop("class not supported")
    }
   posterior_samples

}
