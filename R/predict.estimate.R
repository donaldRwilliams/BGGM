#' Prediction from an \code{estimate} Object
#' @name predict.estimate
#' @param object object of class \code{estimate}
#' @param iter iters used to compute the residuals
#' @param cred credible interval used for summarizing
#' @param newdata an optional data frame for obtaining predictions (e.g., on test data)
#' @param summary summarize the posterior samples (Default is \code{TRUE}).
#' Setting it to \code{FALSE} can be used to then compute performance metrics.
#' @param ... currently ignored
#'
#' @return \code{summary = TRUE}: 3D array of dimensions n (observations),
#'         4 (posterior summary),
#'         p (number of nodes). \code{summary = FALSE}:
#'         list containing predictions for each variable
#' @examples
#' \donttest{
#' # data
#' Y <- subset(tas, gender == "M")[,-ncol(tas)]
#'
#' # fit model
#' fit <- estimate(Y)
#'
#' # predict
#' predict(fit, iter = 25)
#' }
#' @export
predict.estimate <- function(object, iter = 500, cred = 0.95,
                             newdata = NULL, summary = TRUE,...){

  if(object$iter < iter){
    stop("iter cannot be larger than the number of samples used for fitting")
  }
  p <- object$p
  if(is.null(newdata)){
  dat <- object$dat
  n <- nrow(dat)
  } else {

    dat <- newdata
    n <- nrow(dat)
    if(ncol(dat) != p){
      stop("new data not allowed (number of variables differs)")
    }
  }



  lb <- (1 - cred) / 2
  ub <- 1 - lb

  predict_array <- array(0, dim = c(n, 4, p ),
                        dimnames = list(1:n,
                                        c("Post.mean", "Post.sd",
                                          "Cred.lb", "Cred.ub"),
                                        paste0("node_", 1:p)))

  betas <- inverse_2_beta(object, samples = iter)

  if(isTRUE(summary)){
  for(i in 1:p){
    beta <- betas$betas[[i]]
    pred <- t(sapply(1:iter, function(s) {dat[,-i] %*% t(beta[s,])}))
    predict_array[,1,i] <- apply(pred, 2, mean)
    predict_array[,2,i] <- apply(pred, 2, sd)
    predict_array[,3:4,i] <- t(apply(pred, 2, quantile, prob = c(lb, ub)))
  }

  returned_object <-  predict_array

  class(returned_object) <- "predict.estimate"

  } else {

    pred <- list()

    for(i in 1:p){

      beta <- betas$betas[[i]]
      pred[[i]] <- t(sapply(1:iter, function(s) {dat[,-i] %*% t(beta[s,])}))

      }

  returned_object <- list(pred = pred, dat = dat)
  class(returned_object) <- "predict.estimate"

  }
returned_object
}


#' Print Method for \code{predict.estimate} Objects
#'
#' @param x object of class \code{predict.estimate}
#' @param ... currently ignored
#' @export
print.predict.estimate <- function(x,...){
  if(length(x) != 2){
    class(x) <- ""
    x <- round(x, 3)
    print(x)
  } else {
    cat("'summary = FALSE' not printed. See object contents")
  }
}

