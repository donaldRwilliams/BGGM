#' @title Multivariate Normal Imputation
#'
#' @description Impute values, assuming a multivariate normal distribution, with the posterior
#' predictive distribution.
#'
#' @param Y Matrix (or data frame) of dimensions \emph{n} (observations) by  \emph{p} (variables).
#'
#' @param type Character string. Which type of data for \code{Y} ? The options include \code{continuous},
#' \code{binary}, \code{ordinal}, or \code{mixed}. Note that mixed can be used for data with only
#' ordinal variables. See the note for further details.
#'
#' @param iter Number of iterations (posterior samples; defaults to 5000).
#'
#' @param progress Logical. Should a progress bar be included (defaults to \code{TRUE}) ?
#'
#' @return An object of class \code{mvn_imputation}
#'
#' @examples
#' \donttest{
#'
#' d <- scale(na.omit(subset(ifit, id ==1)[,-1]))[,1:7]
#'
#' Y <- MASS::mvrnorm(5000, rep(0, 7), cor(d))
#'
#' # for checking
#' Ymain <- Y
#'
#' indices <- which( matrix(0, 5000, 7) == 0, arr.ind = TRUE)
#' na_indices <- indices[sample(5:nrow(indices), size = 1000, replace = F),]
#' na_indices <- na_indices[na_indices[,1] > 1,]
#'
#' Y[na_indices] <- NA
#' Y_miss <- ifelse(is.na(Y), 1, 0)
#'
#' # true values
#' true <- unlist(sapply(1:7, function(x)  Ymain[ which(Y_miss[,x] == 1),x] ))
#'
#' # impute
#' fit_missing <- mvn_imputation(Y, progress = FALSE)
#'
#' # correlation
#' cor( true, colMeans(fit_missing$ppc_missing))
#'
#' }
#' @export
mvn_imputation <- function(Y, type = "continuous",
                           iter = 1000,
                           progress = TRUE){

  p <- ncol(Y)

  Y_miss <- ifelse(is.na(Y), 1, 0)

  if(isTRUE(progress)){
    message(paste0("BGGM: Imputing"))
  }

  if(type == "continuous"){

    # impute means
    for(i in 1:p){
      Y[which(is.na( Y[,i])) ,i] <- mean(na.omit(Y[,i]))

    }

    fit <-.Call(
      "_BGGM_missing_gaussian",
      Y =  as.matrix(Y),
      Y_miss = as.matrix(Y_miss),
      Sigma = cov(Y),
      iter__missing = iter,
      progress = progress
      )

  }

  if(isTRUE(progress)){
    message("BGGM: Finished")
  }

  if(any( rowSums(Y_miss) == p)){
    warning("entire row missing. imputed values based on posterior predicted means.")
  }

  returned_object <- fit
  class(returned_object) <- "mvn_imputation"
  returned_object
}
