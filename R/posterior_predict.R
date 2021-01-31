#' Posterior Predictive Distribution
#'
#' @description Draw samples from the posterior predictive distribution.
#'
#' @param object An object of class \code{estimate} or \code{explore}
#'
#' @param iter Numeric. Number of samples from the predictive distribution
#'
#' @param progress Logical. Should a progress bar be included (defaults to \code{TRUE})
#'
#' @return A 3D array containing the predicted datasets
#'
#' @note Currently only implemented for \code{type = "mixed"} and \code{type = "binary"}.
#'       Note the term mixed is confusing, in that it can be used with only, say, ordinal data. In this case,
#'       reestimate the model with \code{type = "mixed"} until all data types are supported.
#'
#' @export
#'
#' @examples
#' \donttest{
#' Y <- gss
#'
#' fit <- estimate(as.matrix(Y),
#'                 impute = TRUE,
#'                iter = 150, type = "mixed")
#'
#' yrep <- posterior_predict(fit, iter = 100)
#' }
posterior_predict <- function(object,
                              iter = 1000,
                              progress = TRUE){


  if(!any(class(object) %in% c("estimate", "explore"))) {
    stop("object must be of class 'estimate' or 'explore'.")
  }
  if(!object$type %in% c("binary", "mixed", "ordinal")){
    stop("type must be 'mixed' or 'binary'")
  }

  Y <- object$Y
  cors <- pcor_to_cor(object, iter = iter)
  cors <- cors$R

  n <- object$n
  p <- object$p

  predicted <- array(0, c(n, p, iter))


  if(isTRUE(progress)){
    pb <- utils::txtProgressBar(min = 0, max = iter, style = 3)
  }


  if(object$type == "mixed"){

    for(s in 1:iter){
    cors_s <- cors[,,s]
    ypred_s <- mvnrnd(n, rep(0, p), cors_s)

    predicted[,,s] <- sapply(1:p, function(j) {
      quantile(na.omit(Y[, j]), pnorm(ypred_s[, j], 0,
                                      sqrt(cors_s[j, j])),
               na.rm = TRUE, type = 1)
    }

    )

    if(isTRUE(progress)){
      utils::setTxtProgressBar(pb, s)
    }
  }

  } else if(object$type == "binary"){

    betas <- t(object$post_samp$beta[,,-c(1:50)])

    for(s in 1:iter){
      cors_s <- cors[,,s]
      yrep_s <- mvnrnd(n = n, betas[s,], cors_s)
      predicted[,,s] <- apply(yrep_s, 2, function(x) {ifelse(x > 0, 1, 0) })

      if(isTRUE(progress)){
        utils::setTxtProgressBar(pb, s)
      }
    }
  } else if(object$type == "ordinal"){
    betas <- t(object$post_samp$beta[,,-c(1:50)])
    thresh <- object$post_samp$thresh[-c(1:50),,]
    K <- ncol(thresh) - 1
    temp <- matrix(0, n, p)

    # nasty loop
    # todo: write in c++
    for(s in 1:iter){
      cors_s <- cors[,,s]
      yrep_s <- mvnrnd(n = n, betas[s,], cors_s)
      for(j in 1:p){
        for(n in 1:n){
          for(i in 1:K){
            if(yrep_s[n,j] > thresh[s,i,j] & yrep_s[n,j] < thresh[s,i+1, j]){
              temp[n, j] <- i
            }
          }
        }
      }
      predicted[,,s] <- temp

      if(isTRUE(progress)){
        utils::setTxtProgressBar(pb, s)
      }
    }

}

  dimnames(predicted)[[2]]  <- colnames(Y)
  class(predicted) <- c("array", "posterior_predict")
  return(predicted)
}
