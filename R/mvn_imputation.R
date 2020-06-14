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
#' @param iter Number of iterations (posterior samples; defaults to 1000).
#'
#' @param progress Logical. Should a progress bar be included (defaults to \code{TRUE}) ?
#'
#' @return An object of class \code{mvn_imputation}
#'
#' @examples
#' \donttest{
#' # obs
#' n <- 5000
#'
#' # n missing
#' n_missing <- 1000
#'
#' # variables
#' p <- 16
#'
#' # data
#' Y <- MASS::mvrnorm(n, rep(0, p), ptsd_cor1)
#'
#' # for checking
#' Ymain <- Y
#'
#' # all possible indices
#' indices <- which(matrix(0, n, p) == 0,
#'                  arr.ind = TRUE)
#'
#' # random sample of 1000 missing values
#' na_indices <- indices[sample(5:nrow(indices),
#'                              size = n_missing,
#'                              replace = FALSE),]
#'
#' # fill with NA
#' Y[na_indices] <- NA
#'
#' # missing = 1
#' Y_miss <- ifelse(is.na(Y), 1, 0)
#'
#' # true values (to check)
#' true <- unlist(sapply(1:p, function(x)
#'   Ymain[which(Y_miss[,x] == 1),x] ))
#'
#' # impute
#' fit_missing <- mvn_imputation(Y, progress = FALSE, iter = 250)
#'
#' # plot
#' plot(x =  true,
#'      y = colMeans(fit_missing$ppc_missing),
#'      main = "BGGM: Imputation",
#'      xlab = "Actual",
#'      ylab = "Posterior Mean")
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
