#' @title Multivariate Normal Imputation
#'
#' @description Impute values, assuming a  multivariate normal distribution, with the posterior
#' predictive distribution. For binary, ordinal, and mixed (a combination of discrete and continuous)
#' data, the values are first imputed for the latent data and then converted to the original scale.
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
#' @param save_all Logical. Should each imputed dataset be stored
#' (defaults to \code{FALSE} which saves the imputed missing values) ?
#'
#' @references
#' \insertAllCited{}
#'
#' @return An object of class \code{mvn_imputation}:
#'
#'\itemize{
#'
#'\item \code{Y} The last imputed dataset.
#'
#'\item \code{ppd_missing} A matrix of dimensions \code{iter} by the number of missing values.
#'
#'\item \code{ppd_mean} A vector including the means of the posterior predictive distribution for
#' the missing values.
#'
#'\item \code{Y_all} An 3D array with \code{iter} matrices of dimensions \emph{n} by \emph{p}
#'(\code{NULL} when \code{save_all = FALSE}).
#'
#'}
#'
#' @details
#' Missing values are imputed with the approach described in \insertCite{hoff2009first;textual}{BGGM}.
#' The basic idea is to impute the missing values with the respective posterior pedictive distribution,
#' given the observed data, as the model is being estimated. Note that the default is \code{TRUE},
#' but this ignored when there are no missing values. If set to \code{FALSE}, and there are missing
#' values, list-wise deletion is performed with \code{na.omit}.
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
#'         Ymain[which(Y_miss[,x] == 1),x] ))
#'
#' # impute
#' fit_missing <- mvn_imputation(Y, progress = FALSE, iter = 250)
#'
#' print(fit_missing, n_rows = 20)
#'
#'
#' # plot
#' plot(x =  true,
#'      y = fit_missing$ppd_mean,
#'      main = "BGGM: Imputation",
#'      xlab = "Actual",
#'      ylab = "Posterior Mean")
#' }
#' @export
mvn_imputation <- function(Y, type = "continuous",
                           iter = 1000,
                           progress = TRUE,
                           save_all = FALSE){
  p <- ncol(Y)

  Y_miss <- ifelse(is.na(Y), 1, 0)

  if(sum(Y_miss) == 0){
    stop("no missing values detected")
    }

  missing_location <- unlist(sapply(1:p, function(x) paste0(which(Y_miss[,x] ==1), "--", x)))

  if(isTRUE(progress)){
    message(paste0("BGGM: Imputing"))
  }

  if(type == "continuous"){

    # impute means
    for(i in 1:p){
      Y[which(is.na(Y[,i])) ,i] <- mean(na.omit(Y[,i]))
      }

    # scale data
    Y <- scale(Y)

    fit <-.Call(
      "_BGGM_missing_gaussian",
      Y =  as.matrix(Y),
      Y_miss = as.matrix(Y_miss),
      Sigma = cov(Y),
      iter_missing = iter,
      progress = progress,
      store_all = save_all
      )
    }

  if(isTRUE(progress)){
    message("BGGM: Finished")
  }

  if(!save_all){
    fit$Y_all <- NULL
  }

  if(any(rowSums(Y_miss) == p)){
    warning("entire row missing. imputed values based on posterior predicted means")
  }

  fit$ppd_summary <- cbind(missing_location,
                           as.data.frame(round( fit$ppd_summary, 3)))

  colnames(fit$ppd_summary) <- c("Value",
                                 "Post.mean",
                                 "Post.sd",
                                 "Cred.lb",
                                 "Cred.ub")

  returned_object <- fit
  class(returned_object) <- c("BGGM", "mvn_imputation")
  returned_object
}


print_mvn_impute <- function(x, n_rows = NULL, ...) {
  cat("BGGM: Bayesian Gaussian Graphical Models \n")
  cat("--- \n")
  cat("Multivariate Normal Imputation\n")
  cat("--- \n")
  cat("Estimates:\n\n")
  dat <- x$ppd_summary
  if(is.null(n_rows)) {
    print(dat, row.names = FALSE)
  } else{
    print(dat[1:n_rows, ], row.names = FALSE)
  }
}
