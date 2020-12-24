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
#' @note Currently only implemented for \code{type = "mixed"}. Note the term mixed is confusing, in that
#' it can be used with only, say, ordinal data. In this case, re-estimate the model with
#' \code{type = "mixed"}
#' until all data types are supported.
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
  if(object$type != "mixed"){
    stop("type must be 'mixed'.")
  }
  Y <- object$Y
  cors <- pcor_to_cor(object, iter = iter + 1)
  cors <- cors$R[,,-1]

  n <- object$n
  p <- object$p

  predicted <- array(0, c(n, p, iter))


  if(isTRUE(progress)){
    pb <- utils::txtProgressBar(min = 0, max = iter, style = 3)
  }

  for(s in 1:iter){
    cors_s <- cors[,,s]
    ypred_s <- mvnrnd(object$n, rep(0, p), cors_s)

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
  dimnames(predicted)[[2]]  <- colnames(Y)
  class(predicted) <- c("array", "posterior_predict")
  return(predicted)
}
