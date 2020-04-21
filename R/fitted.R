#' Fitted Values for \code{estimate} Objects
#'
#' @inheritParams residuals.estimate
#' @return 3D array of dimensions n (observations),
#'         4 (posterior summary),
#'         p (number of nodes)
#' @export
#'
#' @examples
#' \donttest{
#' # data
#' Y <- subset(tas, gender == "M")[,-ncol(tas)]
#'
#' # fit model
#' fit <- estimate(Y)
#'
#' # fitted values
#' fitted(fit, iter = 25)
#'}
fitted.estimate <- function(object, iter = 500,
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


  betas <- inverse_2_beta(object, samples = iter)

  if(isTRUE(summary)){

  fitted_array <- array(0, dim = c(n, 4, p ),
                        dimnames = list(1:n,
                                        c("Post.mean",
                                          "Post.sd",
                                          "Cred.lb",
                                          "Cred.ub"),
                                        paste0("node_", 1:p)))



  for(i in 1:p){
    beta <- betas$betas[[i]]
    pred <- t(sapply(1:iter, function(s) {dat[,-i] %*% t(beta[s,])}))
    fitted_array[,1,i] <- apply(pred, 2, mean)
    fitted_array[,2,i] <- apply(pred, 2, sd)
    fitted_array[,3:4,i] <- t(apply(pred, 2, quantile, prob = c(lb, ub)))
  }

  returned_object <-  fitted_array
  class(returned_object) <- c("BGGM", "estimate", "fitted")

  } else {

    pred <- list()

    for(i in 1:p){

      beta <- betas$betas[[i]]
      pred[[i]] <- t(sapply(1:iter, function(s) {dat[,-i] %*% t(beta[s,])}))

      }

    returned_object <- list(pred = pred, dat = dat)
    class(returned_object) <- c("estimate", "fitted")

    }

return(returned_object)

}

