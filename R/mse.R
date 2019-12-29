#' Mean Square Error for \code{predict.estimate} Objects
#' @name mse
#' @param object object of class \code{predict.estimate}
#' @param ... currently ignored
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
#' mse(pred)
mse <- function(object, ...){

  # data
  dat <- object$dat

  # number of variables
  p <- ncol(object$dat)

  # predictions
  pred <- object$pred

  # scores
  scores <- lapply(1:p, function(x) colMeans((t(pred[[x]]) - dat[,x])^2))

  # returned object
  returned_object <- list(scores = scores, metric = "mse")
  class(returned_object) <- "metric"
  return(returned_object)

  }
