#' Model Predictions for \code{estimate} Objects
#'
#' @name predict.estimate
#'
#' @param object object of class \code{estimate}
#'
#' @param iter number of posterior samples (defaults to all in the object).
#'
#' @param cred credible interval used for summarizing
#'
#' @param newdata an optional data frame for obtaining predictions (e.g., on test data)
#'
#' @param summary summarize the posterior samples (defaults to \code{TRUE}).
#'
#' @param progress Logical. Should a progress bar be included (defaults to \code{TRUE}) ?
#'
#' @param ... currently ignored
#'
#' @return \code{summary = TRUE}: 3D array of dimensions n (observations),
#'         4 (posterior summary),
#'         p (number of nodes). \code{summary = FALSE}:
#'         list containing predictions for each variable
#'
#'
#' @examples
#' \donttest{
#' # # data
#' Y <- ptsd
#'
#' fit <- estimate(Y, iter = 250,
#'                 progress = FALSE)
#'
#' pred <- predict(fit,
#'                 progress = FALSE)
#' }
#'
#' @export
predict.estimate <- function(object,
                           newdata = NULL,
                           summary = TRUE,
                           cred = 0.95,
                           iter = NULL,
                           progress =  TRUE,
                           ...){

  # lower bound
  lb <- (1 - cred) / 2

  # uppder bound
  ub <- 1 - lb


  if(is.null(iter)){

    iter <- object$iter

  }

  # correlations
  cors <- pcor_to_cor(object, iter = iter)$R

  if(object$type == "continuous") {

    if(is.null(newdata)){

      # data matrix
      Y <- object$Y

      # nodes
      p <- ncol(Y)

      # observations
      n <- nrow(Y)

    } else {

      # scale
      Y <- scale(newdata, scale = F)

      # nodes
      p <- ncol(Y)

      # observations
      n <- nrow(Y)

    }

  } else{

    stop("type not currently supported. must be continuous")
}
  if(object$p != p){

    stop(paste0("the number of nodes in the newdata does",
                "not match the number of nodes in the object"))

  }

  if(isTRUE(progress)){
    pb <- utils::txtProgressBar(min = 0, max = p, style = 3)
  }

  # yhats
  yhats <- lapply(1:p, function(x) {

    yhat_p <- .Call("_BGGM_pred_helper_latent",
                    Y = Y[,-x],
                    XX = cors[-x, -x,],
                    Xy = cors[x, -x,],
                    quantiles = c(lb, ub),
                    n = n,
                    iter = iter
    )

    if(isTRUE(progress)){
      utils::setTxtProgressBar(pb, x)
    }

    yhat_p

  })


  # node names
  cn <- colnames(object$Y)

  # check for column names
  if(is.null(cn)) {
    cn <- 1:p
  }

  fitted_array <- array(0, dim = c(n, 4, p))

  dimnames(fitted_array)[[2]] <- c("Post.mean",
                                   "Post.sd",
                                   "Cred.lb",
                                   "Cred.ub")


  dimnames(fitted_array)[[3]] <- cn

  if(isTRUE(summary)){

    for(i in 1:p){

      fitted_array[,,i] <- cbind(t(as.matrix(yhats[[i]]$yhat_mean)),
                                 t(as.matrix(yhats[[i]]$yhat_sd)),
                                 t(yhats[[i]]$yhat_quantiles))
    }

  } else {

    fitted_array <- array(0, dim = c(iter, n, p))

    dimnames(fitted_array)[[3]] <- cn

    for(i in 1:p){

      fitted_array[,,i] <- as.matrix(yhats[[i]]$yhat)

    }
  }

  return(fitted_array)

}




#' Model Predictions for \code{explore} Objects
#'
#' @name predict.explore
#'
#' @param object object of class \code{explore}
#'
#' @param iter number of posterior samples (defaults to all in the object).
#'
#' @param cred credible interval used for summarizing
#'
#' @param newdata an optional data frame for obtaining predictions (e.g., on test data)
#'
#' @param summary summarize the posterior samples (defaults to \code{TRUE}).
#'
#' @param progress Logical. Should a progress bar be included (defaults to \code{TRUE}) ?
#'
#' @param ... currently ignored
#'
#' @return \code{summary = TRUE}: 3D array of dimensions n (observations),
#'         4 (posterior summary),
#'         p (number of nodes). \code{summary = FALSE}:
#'         list containing predictions for each variable
#'
#' @examples
#' \donttest{
#' # data
#' Y <- ptsd
#'
#' # fit model
#' fit <- explore(Y, iter = 250,
#'                progress = FALSE)
#'
#' # predict
#' pred <- predict(fit,
#'                 progress = FALSE)
#'
#' }
#' @export
predict.explore <- function(object,
                             newdata = NULL,
                             summary = TRUE,
                             cred = 0.95,
                             iter = NULL,
                             progress = TRUE,
                             ...){

  # lower bound
  lb <- (1 - cred) / 2

  # uppder bound
  ub <- 1 - lb


  if(is.null(iter)){

    iter <- object$iter

  }

  # correlations
  cors <- pcor_to_cor(object, iter = iter)$R

  if(object$type == "continuous") {

    if(is.null(newdata)){

      # data matrix
      Y <- object$Y

      # nodes
      p <- ncol(Y)

      # observations
      n <- nrow(Y)

    } else {

      # scale
      Y <- scale(newdata, scale = F)

      # nodes
      p <- ncol(Y)

      # observations
      n <- nrow(Y)

    }

  } else{

    stop("type not currently supported. must be continuous")
  }



  if(object$p != p){

    stop(paste0("the number of nodes in the newdata does",
                "not match the number of nodes in the object"))

  }

  if(isTRUE(progress)){
    pb <- utils::txtProgressBar(min = 0, max = p, style = 3)
  }

  # yhats
  yhats <- lapply(1:p, function(x) {

    yhat_p <- .Call("_BGGM_pred_helper_latent",
                    Y = Y[,-x],
                    XX = cors[-x, -x,],
                    Xy = cors[x, -x,],
                    quantiles = c(lb, ub),
                    n = n,
                    iter = iter
    )

    if(isTRUE(progress)){
      utils::setTxtProgressBar(pb, x)
    }

    yhat_p

  })


  # node names
  cn <- colnames(object$Y)

  # check for column names
  if(is.null(cn)) {
    cn <- 1:p
  }

  fitted_array <- array(0, dim = c(n, 4, p))

  dimnames(fitted_array)[[2]] <- c("Post.mean",
                                   "Post.sd",
                                   "Cred.lb",
                                   "Cred.ub")


  dimnames(fitted_array)[[3]] <- cn

  if(isTRUE(summary)){

    for(i in 1:p){

      fitted_array[,,i] <- cbind(t(as.matrix(yhats[[i]]$yhat_mean)),
                                 t(as.matrix(yhats[[i]]$yhat_sd)),
                                 t(yhats[[i]]$yhat_quantiles))
    }

  } else {

    fitted_array <- array(0, dim = c(iter, n, p))

    dimnames(fitted_array)[[3]] <- cn

    for(i in 1:p){

      fitted_array[,,i] <- as.matrix(yhats[[i]]$yhat)

    }
  }

  return(fitted_array)

}


#' Model Predictions for \code{var_estimate} Objects
#'
#' @name predict.var_estimate
#'
#' @param object object of class \code{var_estimate}
#'
#' @param summary summarize the posterior samples (defaults to \code{TRUE}).
#'
#' @param cred credible interval used for summarizing
#'
#' @param iter number of posterior samples (defaults to all in the object).
#'
#' @param progress Logical. Should a progress bar be included (defaults to \code{TRUE}) ?
#'
#' @param ... Currently ignored
#'
#' @return The predicted values for each regression model.
#'
#' @examples
#' \donttest{
#' # data
#' Y <- subset(ifit, id == 1)[,-1]
#'
#' # fit model with alias (var_estimate also works)
#' fit <- var_estimate(Y, progress = FALSE)
#'
#' # fitted values
#' pred <- predict(fit, progress = FALSE)
#'
#' # predicted values (1st outcome)
#' pred[,,1]
#'
#' }
#' @export
predict.var_estimate <- function(object,
                                 summary = TRUE,
                                 cred = 0.95,
                                 iter = NULL,
                                 progress = TRUE,
                                 ...){


  # lower bound
  lb <- (1 - cred) / 2

  # uppder bound
  ub <- 1 - lb


    # data matrix
    X <- object$X
    n <- nrow(X)


  if(is.null(iter)){

    iter <- object$iter

  }

  p <- object$p

  post_names <- sapply(1:p, function(x) paste0(
    colnames(object$Y)[x], "_",  colnames(object$X))
  )

  post_samps <- posterior_samples(object)

  if(isTRUE(progress)){
    pb <- utils::txtProgressBar(min = 0, max = p, style = 3)
  }

  yhats <- lapply(1:p, function(x){

    yhat_p <- post_samps[, post_names[,x]] %*% t(X)


    if(isTRUE(progress)){
      utils::setTxtProgressBar(pb, x)
    }

    yhat_p

  })

  if(isTRUE(summary)){



    fitted_array <- array(0, dim = c(n, 4, p))

    dimnames(fitted_array)[[2]] <- c("Post.mean",
                                     "Post.sd",
                                     "Cred.lb",
                                     "Cred.ub")

    dimnames(fitted_array)[[3]] <- colnames(object$Y)


    for(i in 1:p){

      fitted_array[,,i] <- cbind(colMeans(yhats[[i]]),
                                 apply(yhats[[i]], 2, sd),
                                 apply(yhats[[i]], 2, quantile, lb),
                                 apply(yhats[[i]], 2, quantile, ub)
      )
    }


  } else {

    fitted_array <- array(0, dim = c(iter, n, p))

    dimnames(fitted_array)[[3]] <- colnames(object$Y)

    for(i in 1:p){

      fitted_array[,,i] <- t(as.matrix(yhats[[i]]))

    }

  }

  return(fitted_array)

}

