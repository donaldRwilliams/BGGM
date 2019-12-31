#' Mean Squared Error
#' @description Compute mean squared error for either the observed data or future data. The former
#' is computed by plugging in the observed y (the predicted node). The latter is computed from
#' replicated data sets (posterior predictive y), which results in the posterior predictive
#' mean squared error.  Both provide a measure of uncertainty,  as the error is computed from
#' the posterior samples. However, the posterior predictive approach fully captures uncertainty.
#'
#' @name mse
#' @param object object of class \code{post.pred} or  \code{predict.estimate}
#' @param ... currently ignored
#'
#' @return object of class \code{metric}
#' @export
#'
#' @examples
#' # data
#' Y <- subset(tas, gender == "M")[,-ncol(tas)]
#'
#' # fit model
#' fit <- estimate(Y)
#'
#' # predict (note summary = FALSE)
#' pred <- predict(fit, summary = FALSE)
#'
#' mse(pred)
mse <- function(object, ...){

  # data
  dat <- object$dat

  # number of variables
  p <- ncol(object$dat)

  # check summary is false
  if(length(object) != 2){
    stop("summary must be set to false")
  }

  # predictions
  pred <- object$pred

  if(class(object) == "predict.estimate"){

    type <- class(object)

    # scores
    scores <- lapply(1:p, function(x) colMeans((t(pred[[x]]) - dat[,x])^2))

    } else if (class(object) == "post.pred"){

      type <- class(object)

      scores <- lapply(1:p, function(x) colMeans((t(pred[[x]]) - dat[,x])^2))

    } else {

      stop("object class not supported (must be predict.estimate or post.pred)")

      }

  # returned object
  returned_object <- list(scores = scores,
                          type = type,
                          metric = "mse")

  class(returned_object) <- "metric"

  return(returned_object)

  }

#' Mean Absolute Error for \code{predict.estimate} Objects
#' @name mae
#' @inheritParams mse
#' @return
#' @export
#'
#' @examples
#' # data
#' Y <- subset(tas, gender == "M")[,-ncol(tas)]
#'
#' # fit model
#' fit <- estimate(Y)
#'
#' # predict (note summary = FALSE)
#' pred <- predict(fit, summary = FALSE)
#'
#' mae(pred)
mae <- function(object, ...){

  # data
  dat <- object$dat

  # number of variables
  p <- ncol(object$dat)

  # predictions
  pred <- object$pred


  if(class(object) == "predict.estimate"){

    type <- class(object)

    # scores
    scores <- lapply(1:p, function(x) colMeans(abs(t(pred[[x]]) - dat[,x])))

  } else if (class(object) == "post.pred"){

    type <- class(object)

    # scores
    scores <- lapply(1:p, function(x) colMeans(abs(t(pred[[x]]) - dat[,x])))

  } else {

    stop("object class not supported (must be predict.estimate or post.pred)")

  }

  returned_object <- list(scores = scores,
                          metric = "mae",
                          type = type)

  class(returned_object) <- "metric"
  return(returned_object)
}

#' Root Mean Square Error for \code{predict.estimate} Objects
#' @name rmse
#' @inheritParams mse
#' @return
#' @export
#'
#' @examples
#' # data
#' Y <- subset(tas, gender == "M")[,-ncol(tas)]
#'
#' # fit model
#' fit <- estimate(Y)
#'
#' # predict (note summary = FALSE)
#' pred <- predict(fit, summary = FALSE)
#'
#' rmse(pred)
rmse <- function(object, ...){

  # data
  dat <- object$dat

  # number of variables
  p <- ncol(object$dat)

  # predictions
  pred <- object$pred

  if(class(object) == "predict.estimate"){

    type <- class(object)

    # scores
    scores <- lapply(1:p, function(x) sqrt(colMeans((t(pred[[x]]) - dat[,x])^2)))

  } else if (class(object) == "post.pred"){

    type <- class(object)

    # scores
    scores <- lapply(1:p, function(x) sqrt(colMeans((t(pred[[x]]) - dat[,x])^2)))

  } else {

    stop("object class not supported (must be predict.estimate or post.pred)")

  }

  # returned object
  returned_object <- list(scores = scores,
                          metric = "rmse",
                          type = type)
  class(returned_object) <- "metric"
  return(returned_object)

}

#' Mean Absolute Percentage Error for \code{predict.estimate} Objects
#' @name mape
#' @inheritParams mse
#'
#' @return
#' @export
#'
#' @examples
#' # data
#' Y <- subset(tas, gender == "M")[,-ncol(tas)]
#'
#' # fit model
#' fit <- estimate(Y)
#'
#' # predict (note summary = FALSE)
#' pred <- predict(fit, summary = FALSE)
#'
#' mape(pred)
mape <- function(object, ...){

  # data
  dat <- object$dat

  # number of variables
  p <- ncol(object$dat)

  # predictions
  pred <- object$pred


  if(class(object) == "predict.estimate"){

    type <- class(object)

    # scores
    scores <- lapply(1:p, function(x) colMeans((abs((dat[,x] - t(pred[[x]]))/dat[,x]))))

  } else if (class(object) == "post.pred"){

    type <- class(object)

    # scores
    scores <- lapply(1:p, function(x) colMeans((abs((dat[,x] - t(pred[[x]]))/dat[,x]))))

  } else {

    stop("object class not supported (must be predict.estimate or post.pred)")

  }

  # returned object
  returned_object <- list(scores = scores,
                          metric = "mape",
                          type = type)

  class(returned_object) <- "metric"
  return(returned_object)

}


#' Summary Method for \code{metric} Objects
#'
#' @param object object of class \code{metric}
#' @param cred credible interval
#' @param ... currently ignored
#'
#' @return
#' @export
summary.metric <- function(object, cred = 0.95, ...){

  cred <- 0.95

  lb <- (1 - cred) / 2

  ub <- 1 - lb

  p <- length(object$scores)

  iter <- length(object$scores[[1]])
  dat_summ <- data.frame(Node = 1:p,
                         Post.mean  = sapply(object$scores, mean),
                         Post.sd = sapply(object$scores, sd),
                         Cred = t(sapply(object$scores,
                                         quantile,
                                         c(lb, ub))))
  returned_object <- list(summary = dat_summ,
                          metric = object$metric,
                          type = object$type,
                          iter = iter,
                          cred = cred)
  class(returned_object) <- c("summary.metric", "data.frame")
  returned_object
}

#' Print Method for \code{summary.metric} Object
#'
#' @param x object of class \code{summary.metric}
#' @param digits digits used to round the values
#' @param ... currently ignored
#'
#' @return
#' @export
print.summary.metric <- function(x, digits = 2,...){
  cat("BGGM: Bayesian Gaussian Graphical Models \n")
  cat("--- \n")
  cat("Metric:", x$metric, "\n")
  cat("Type:", x$type, "\n")
  cat("Credible Interval:", x$cred, "\n")
  cat("--- \n")
  cat("Estimates:\n\n")
  dat <- x$summary
  colnames(dat) <- c(colnames(dat)[1:3], "Cred.lb", "Cred.ub")
  print(as.data.frame( sapply(dat , round, digits)),
        row.names = FALSE)
}

#' Print Method for \code{metric} Objects
#'
#' @param x object of class \code{metric}
#' @param ... currently ignored
#'
#' @return
#' @export
print.metric <- function(x,...){
  print(summary(x))
}




