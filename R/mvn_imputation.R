#' @title Obtain Imputed Datasets
#'
#' @description Impute missing values, assuming a  multivariate normal distribution, with the posterior
#' predictive distribution. For binary, ordinal, and mixed (a combination of discrete and continuous)
#' data, the values are first imputed for the latent data and then converted to the original scale.
#'
#' @param Y Matrix (or data frame) of dimensions \emph{n} (observations) by  \emph{p} (variables).
#'
#' @param type Character string. Which type of data for \code{Y} ? The options include \code{continuous},
#' \code{binary}, \code{ordinal}, or \code{mixed}. Note that mixed can be used for data with only
#' ordinal variables. See the note for further details.
#'
#' @param lambda Numeric. A regularization parameter, which defaults to p + 2. A larger value results
#'               in more shrinkage.
#'
#' @param mixed_type Numeric vector. An indicator of length \emph{p} for which variables should be treated as ranks.
#' (1 for rank and 0 to assume the observed marginal distribution).
#' The default is currently to treat all integer variables as ranks when
#' \code{type = "mixed"} and \code{NULL} otherwise. See note for further details.
#'
#' @param iter Number of iterations (posterior samples; defaults to 1000).
#'
#' @param progress Logical. Should a progress bar be included (defaults to \code{TRUE}) ?
#'
#'
#' @references
#' \insertAllCited{}
#'
#' @return An object of class \code{mvn_imputation}:
#'
#'\itemize{
#'
#'\item \code{imputed_datasets} An array including the imputed datasets.
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
#' fit_missing <- impute_data(Y, progress = FALSE, iter = 250)
#'
#' # impute
#' fit_missing <- impute_data(Y,
#'                            progress = TRUE,
#'                            iter = 250)
#'
#' # plot
#' plot(BGGM:::mean_array( fit_missing$imputed_datasets)[na_indices], Ymain[na_indices])
#'
#' }
#' @export
impute_data <- function(Y,
                        type = "continuous",
                        lambda = NULL,
                        mixed_type = NULL,
                        iter = 1000,
                        progress = TRUE){

  if(!type %in% c("continuous", "mixed")){

    stop(paste0("currently only 'continuous' and 'mixed' data are supported."))
  }

  p <- ncol(Y)

  if(is.null(lambda)){
    lambda <- p + 2
  }

  if(is.null(mixed_type)){
    idx <- rep(1, p)
  } else {
    idx <- mixed_type
  }
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

    fit <-.Call(
      "_BGGM_missing_gaussian",
      Y =  as.matrix(Y),
      Y_miss = as.matrix(Y_miss),
      Sigma = cov(Y),
      iter_missing = iter,
      progress = progress,
      store_all = TRUE,
      lambda = lambda
      )

    names(fit) <- "imputed_datasets"

    } else if(type == "mixed"){

      rank_help <- rank_helper(Y)
      rank_help$levels[na_indices] <- NA
      rank_help$z0_start[is.na(rank_help$z0_start)] <- rnorm(sum(Y_missing))

     fit <- .Call(
        "_BGGM_missing_copula_data",
        Y = Y,
        Y_missing = Y_miss,
        z0_start = rank_help$z0_start,
        Sigma_start = cov(rank_help$z0_start),
        levels = rank_help$levels,
        iter_missing = iter,
        progress_impute = TRUE,
        K = rank_help$K,
        idx = idx,
        lambda = lambda
      )

     names(fit) <- "imputed_datasets"

  }

  if(isTRUE(progress)){
    message("BGGM: Finished")
  }

  returned_object <- fit
  class(returned_object) <- c("BGGM", "mvn_imputation")
  returned_object
}


print_mvn_impute <- function(x, ...) {
  cat("BGGM: Bayesian Gaussian Graphical Models \n")
  cat("--- \n")
  cat("Multivariate Normal Imputation\n")
  cat("--- \n")
  cat(date())
}
