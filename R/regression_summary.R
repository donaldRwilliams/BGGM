#' Summarary Method for Multivariate or Univarate Regression
#'
#' @param object An object of class \code{estimate}
#'
#' @param cred Numeric. The credible interval width for summarizing the posterior
#' distributions (defaults to 0.95; must be between 0 and 1).
#'
#' @param ... Currently ignored
#'
#' @return A list of length \emph{p} including the
#'         summaries for each regression.
#'
#' @export
#'
#' @examples
#' \donttest{
#' # note: iter = 250 for demonstrative purposes
#'
#' # data
#' # Y <- bfi
#'
#' Y <- subset(Y, select = c("E5", "N5",
#'                           "gender", "education"))
#'
#'
#' fit_mv_ordinal <- estimate(Y, formula = ~ gender + as.factor(education),
#'                            type = "ordinal",
#'                            iter = 250)
#'
#' fit_mv_ordinal
#'}
regression_summary <- function(object, cred = 0.95, ...){


  if(!all(c("estimate", "default") %in% class(object))){

    stop("class not supported. must be an estimate object")

    }

  lb <- (1-cred)/2
  ub <- 1 - lb

  iter <- object$iter

  beta <- object$post_samp$beta[,,51:(iter + 50)]

  dims <- dim(beta)[1:2]

  post_mean <- apply(beta, 1:2, mean)
  post_sd <- apply(beta, 1:2, sd)
  post_lb <- apply(beta, 1:2, quantile, lb)
  post_ub <- apply(beta, 1:2, quantile, ub)


  outcomes <- dims[2]

  summ <- list()

  for(i in 1:outcomes){


    summ[[i]] <- round(data.frame(Post.mean = post_mean[,i],
                                  Post.sd = post_mean[,i],
                                  Cred.lb = post_lb[,i],
                                  Cred.ub = post_ub[,i] ), 3)

    rownames(  summ[[i]]) <- colnames(object$X)

  }

  # check colnames
  cn <- colnames(object$Y)

  if(is.null(cn)){
    cn <- 1:outcomes
  }

  # colnames
  names(summ) <- cn

  # correlation
  cors <- pcor_to_cor(object)$R

  # residual correlation mean
  cor_mean <- apply(cors, 1:2, mean)

  colnames(cor_mean) <- cn
  rownames(cor_mean) <- cn

  object$post_samp <- NULL

  returned_object <- list(reg_summary = summ,
                          resid_cor = cor_mean,
                          object = object)

  class(returned_object) <- c("BGGM",
                              "regression_summary")
  returned_object
}



print_regression_summary <- function(x, ...){
  cat("BGGM: Bayesian Gaussian Graphical Models \n")
  cat("--- \n")
  cat("Type:",  x$object$type, "\n")
  cat("Formula:", paste(as.character(x$object$formula), collapse = " "), "\n")
  cat("--- \n")
  outcomes <- length(x$reg_summary)
  cat("Coefficients: \n \n")

  for(i in 1:outcomes){
    cat(names(x$reg_summary)[i], "\n")
    print(x$reg_summary$E5)
    cat("--- \n")
  }

cat("Residual Correlation Matrix: \n")
  print(round(x$resid_cor, 3))
  cat("--- \n")
}
