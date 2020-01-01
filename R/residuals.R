#' Residuals for \code{estimate} Objects
#' @name residuals.estimate
#' @param object object of class \code{estimate}
#' @param iter iters used to compute the residuals
#' @param cred credible interval used for summarizing
#' @param summary summarize the posterior samples (Default is \code{TRUE}).
#'
#' @param ... currently ignored
#'
#' @return 3D array of dimensions n (observations),
#'         4 (posterior summary),
#'         p (number of nodes)
#' @export
#'
#' @examples
#' # data
#' Y <- subset(tas, gender == "M")[,-ncol(tas)]
#'
#' # fit model
#' fit <- estimate(Y)
#'
#' # diagnostic plot
#' residuals(fit, iter = 25)
residuals.estimate <- function(object, iter = 500,
                               cred = 0.95,
                               summary = TRUE,
                               ...){


  if(object$iter < iter){
    stop("iter cannot be larger than the number of samples used for fitting")
  }

  dat <- object$dat
  p <- object$p
  n <- nrow(dat)
  lb <- (1 - cred) / 2
  ub <- 1 - lb

  resid_array <- array(0, dim = c(n, 4, p ),
                       dimnames = list(1:n,
                                       c("Post.mean", "Post.sd", "Cred.lb", "Cred.ub"),
                                       paste0("node_", 1:p)))

  betas <- inverse_2_beta(object, samples = iter)

  for(i in 1:p){
    beta <- betas$betas[[i]]
    pred <- t(sapply(1:iter, function(s) {dat[,-i] %*% t(beta[s,]) - dat[,i]}))
    resid_array[,1,i] <- apply(pred, 2, mean)
    resid_array[,2,i] <- apply(pred, 2, sd)
    resid_array[,3:4,i] <- t(apply(pred, 2, quantile, prob = c(lb, ub)))
  }

  resid_array

}
