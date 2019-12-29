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

  # scores
  scores <- lapply(1:p, function(x) sqrt(colMeans((t(pred[[x]]) - dat[,x])^2)))

  # returned object
  returned_object <- list(scores = scores, metric = "rmse")
  class(returned_object) <- "metric"
  return(returned_object)

}
