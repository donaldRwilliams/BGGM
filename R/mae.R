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

  # scores
  scores <- lapply(1:p, function(x) colMeans(abs(t(pred[[x]]) - dat[,x])))


  returned_object <- list(scores = scores, metric = "mae")
  class(returned_object) <- "metric"

  return(returned_object)

}
