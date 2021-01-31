#' Predicted Probabilities
#'
#' Compute the predicted probabilities for discrete data, with the possibility
#' of conditional predictive probabilities (i.e., at fixed values of other nodes)
#'
#' @param object An object of class \code{posterior_predict}
#'
#' @param outcome Character string. Node for which the probabilities are computed.
#'
#' @param Y Matrix (or data frame) of dimensions \emph{n} (observations) by  \emph{p} (variables).
#'        This must include the column names.
#'
#' @param ... Compute conditional probabilities by specifying a column name in \code{Y}
#'            (besides the \code{outcome}) and a fixed value. This can include
#'            any number of nodes. See example below. Leave this blank to compute
#'            unconditional probabilities for \code{outcome}.
#'
#' @return A list containing a matrix with the computed probabilities
#'        (a row for each predictive sample and a column for each category).
#'
#'
#' @note There are no checks that the conditional probability exists, i.e., suppose
#'       you wish to condition on, say, B3 = 2 and B4 = 1, yet there is no instance in
#'       which B3 is 2 AND B4 is 1. This will result in an uninformative error.
#'
#' @export
#'
#' @examples
#' \donttest{
#' Y <- ptsd
#' fit <- estimate(as.matrix(Y), iter = 150, type = "mixed")
#'
#' pred <- posterior_predict(fit, iter = 100)
#'
#' prob <- predicted_probability(pred,
#'                               Y = Y,
#'                               outcome = "B3",
#'                               B4 = 0,
#'                               B5 = 0)
#'
#' }
predicted_probability <- function(object,  outcome, Y, ...){

  # note: test is for checking the subsetting of the 3d arrays
  unique_values <- sort(unique(Y[,outcome]))
  K <- length(unique_values)
  iter <- dim(object)[3]
  collect <- matrix(0, iter, K)

  if(!is(object, "posterior_predict")){
    stop("must be of class 'posterior_predict'")
  }

  dots <- list(...)
  test <- list()

  if(length(dots) == 0){

    for(i in 1:iter){
      collect[i,] <- sapply(1:K, function(x) sum(object[,outcome,i] == unique_values[x])  )
    }



  }  else {
    text_eval <- sapply(1:length(dots), function(x) {
      paste("sub_set[,", x, "] == ", dots[[x]])

    })

    for(i in 1:iter){

      sub_set <- as.matrix(object[,names(dots),i])

      conditional <- eval(parse(text =
                                  paste("object[which(", paste(unlist(text_eval),
                                                               collapse = "  & "), "),,", i,"]")
      ))

      test[[i]] <- conditional

      collect[i, ] <- sapply(1:K, function(x) sum(conditional[,outcome]
                                                  == unique_values[x]))
    }
  }
  collect <- t(apply(collect, 1,function(x){x / sum(x)}))
  colnames(collect) <-  unique_values
  returned_object <- list(collect = collect, sub_sets = test)
  return(returned_object)

}
