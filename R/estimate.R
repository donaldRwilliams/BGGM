#' @title GGM: Estimation
#'
#' @description Estimate the conditional (in)dependence with either an analytic solution or efficiently
#' sampling from the posterior distribution. These methods were introduced in \insertCite{Williams2019;textual}{BGGM}.
#' The graph is then selected with \code{\link{select.estimate}}, with either directional posterior probabilities
#' \insertCite{Marsman2017a}{BGGM}, credible intervals, or a region of practical equivalence \insertCite{Kruschke2017}{BGGM}.
#' Bayesian hypothesis testing is implemented in \code{\link{explore}} and \code{\link{confirm}} \insertCite{Williams2019_bf}{BGGM}.
#'
#' @name estimate
#' @param Y  matrix (or data frame) of dimensions \emph{n} (observations) by  \emph{p} (variables).
#'
#' @param formula an object of class \code{\link[stats]{formula}}. This allows for including
#' control variables in the model (i.e., \code{~ gender}). See the note for further details.
#'
#' @param type character string. Which type of data for \strong{Y} ? The options include \code{continuous},
#' \code{binary}, \code{ordinal}, or \code{mixed}. See the note for further details.
#'
#' @param mixed_type numeric vector. An indicator of length p for which varibles should be treated as ranks.
#' (1 for rank and 0 to assume normality). The default is currently to treat all integer variables as ranks
#' when \code{type = "mixed"} and \code{NULL} otherwise. See note for further details.
#'
#' Treating continous data as ranks is computationally expensive and can be avoided by assuming normality.
#'
#' @param iter number of iterations (posterior samples; defaults to 5000).
#'
#' @param analytic logical. Should the analytic solution be computed (default is \code{FALSE})?
#'
#' @param ... currently ignored.
#'
#' @references
#' \insertAllCited{}
#'
#'
#' @return  list of class \code{estimate}:
#'
#' \code{analytic = TRUE}:
#' \itemize{
#' \item \code{fit} list of analytic solution estimates
#' \itemize{
#' \item \code{inv_mu} inverse covariance matrix (mean)
#' \item \code{inv_var} inverse covariance matrix (variance)
#' \item \code{partial} partial correlation matrix
#' }
#' \item \code{analytic} TRUE
#' \item \code{call} match.call()
#' \item \code{dat} data matrix
#' \item \code{p} number of variables
#' }
#'
#' \code{analytic = FALSE}:
#' \itemize{
#' \item \code{parcors_mat} partial correlation matrix
#' \item \code{inv_mat} inverse covariance matrix
#' \item \code{posterior samples} posterior samples for partial correlations and inverse covariance matrix
#' \item \code{p} number of variables
#' \item \code{dat} data matrix
#' \item \code{iter} number of posterior samples
#' \item \code{call} match.call()
#' \item \code{analytic} FALSE
#' }
#'
#'
#' @note The default is to draws samples from the posterior distribution (\code{analytic = FALSE}). The samples are
#' required for computing edge differences (see \code{\link{ggm_compare_estimate}}), Bayesian R2 introduced in
#' \insertCite{gelman_r2_2019;textual}{BGGM} (see \code{\link{bayes_R2}}), etc. If the goal is to *only* determine
#' the non-zero effects, this can be accomplished by setting \code{analytic = TRUE}. Note also sampling is very fast--i.e.,
#' less than 1 second with p = 25, n = 2500 and 5,000 samples.
#'
#' These methods are inherently Bayesian. This also means there is a close correspondence to "frequentist" methods. The prior
#' distribution for \code{analytic = TRUE} is set to mimic the unbiased estimate of the sample based
#' covariance matrix. When \code{analytic = FALSE}, samples are efficiently drawn from the posterior distribution.
#' The advantage compared to frequentist methods is that a measure of uncertainty is readily available. This allows for
#' seamlessly comparing partial correlation (see \code{\link{ggm_compare_estimate}}) and computing Bayesian R2
#' (see \code{\link{test.R2}}). Further, the posterior probability of a null region (see \code{\link{select.estimate}})
#' provides an estimate of the conditional independence structure (null effects).
#'
#'
#'\strong{Controlling for Variables}:
#'
#' When controlling for variables, it is assumed that \code{Y} includes \emph{only}
#' the nodes in the GGM and the control variables. Internally, \code{only} the predictors
#' that are included in \code{formula} are removed from \code{Y}. An example is provided below.
#'
#' \strong{Dealing with Errors}:
#'
#' An error is most likely to arise when \code{type = "ordinal"}. The are two common errors (although still rare).
#' The first is due to ampling the thresholds, especially when the data is
#' heavily skewed. This can result in an ill-defined matrix. If this occurs, we recommend to first try
#' decreasing \code{prior_sd} (i.e., a more informative prior). If that does not work, then  change the
#' data type to \code{type = mixed}. This should work without a problem.
#'
#' The second error is from the categories. For example, if the error staes that the index is out of bounds, this
#' means there are zeros in the data. The first category must be \code{1}. This can easily be addressed by adding 1 to the
#' data matrix.
#'
#'
#' \strong{Interpretation of conditional (in)dependence models for latent data}:
#'
#' See \code{\link{BGGM}} for details about the interpretation of GGMs based on latent data
#' (i.e, all data types besides \code{"continuous"})
#'
#'
#'
#' @examples
#' # p = 5
#' Y <- BGGM::bfi[, 1:5]
#'
#' # analytic approach (sample by setting analytic = FALSE)
#' fit_analytic <- estimate(Y, analytic = TRUE)
#'
#' # select the graph (edge set E)
#' E <- select(fit_analytic, ci_width = 0.95)
#'
#' @export
estimate  <- function(Y,
                      formula = NULL,
                      type = "continuous",
                      mixed_type = NULL,
                      analytic = FALSE,
                      prior_sd = 0.25,
                      iter = 5000,
                      ...){

  # delta rho ~ beta(delta/2, delta/2)
  delta <- delta_solve(prior_sd)

  # sample posterior
  if(!analytic){

    message(paste0("BGGM: Posterior Sampling ", ...))

    # continuous
    if(type == "continuous"){

      # no control
      if(is.null(formula)){

        # na omit
        Y <- as.matrix(na.omit(Y))

        # scale Y
        Y <- scale(Y, scale = F)

        # nodes
        p <- ncol(Y)

        n <- nrow(Y)

        # posterior sample
        post_samp <- .Call(
          '_BGGM_Theta_continuous',
          PACKAGE = 'BGGM',
          Y = Y,
          iter = iter + 50,
          delta = delta,
          epsilon = 0.1,
          prior_only = 0,
          explore = 1
          )

        # control for variables
        } else {

          control_info <- remove_predictors_helper(list(as.data.frame(Y)),
                                                   formula = formula)

          # data
          Y <- as.matrix(scale(control_info$Y_groups[[1]], scale = F))

          # nodes
          p <- ncol(Y)

          # observations
          n <- nrow(Y)

          # model matrix
          X <- as.matrix(control_info$model_matrices[[1]])

        # posterior sample
        post_samp <- .Call(
          "_BGGM_mv_continuous",
          Y = Y,
          X = X,
          delta = delta,
          epsilon = 0.1,
          iter = iter + 50
        )
      # end control
      }

      # binary
    } else if (type == "binary") {

      # intercept only
      if (is.null(formula)) {

          # data
          Y <- as.matrix(na.omit(Y))

          # obervations
          n <- nrow(Y)

          # nodes
          p <- ncol(Y)

          X <- matrix(1, n, 1)

        } else {

          control_info <- remove_predictors_helper(list(as.data.frame(Y)),
                                                   formula = formula)

          # data
          Y <-  as.matrix(control_info$Y_groups[[1]])

          # observations
          n <- nrow(Y)

          # nodes
          p <- ncol(Y)

          # model matrix
          X <- as.matrix(control_info$model_matrices[[1]])

          }

      # posterior sample
      post_samp <-  .Call(
        "_BGGM_mv_binary",
        Y = Y,
        X = X,
        delta = delta,
        epsilon = 0.1,
        iter = iter + 50,
        beta_prior = 0.0001,
        cutpoints = c(-Inf, 0, Inf)
      )

      # ordinal
      } else if(type == "ordinal"){

      # intercept only
      if(is.null(formula)){

        # data
        Y <- as.matrix(na.omit(Y))

        # obervations
        n <- nrow(Y)

        # nodes
        p <- ncol(Y)

        # intercept only
        X <- matrix(1, n, 1)

        } else {

          control_info <- remove_predictors_helper(list(as.data.frame(Y)),
                                                   formula = formula)

          # data
          Y <-  as.matrix(control_info$Y_groups[[1]])

          # observations
          n <- nrow(Y)

          # nodes
          p <- ncol(Y)

          # model matrix
          X <- as.matrix(control_info$model_matrices[[1]])

          }

        # categories
        K <- max(apply(Y, 2, function(x) { length(unique(x))   } ))

        # call c ++
        post_samp <- .Call(
          "_BGGM_mv_ordinal_albert",
          Y = Y,
          X = X,
          iter = iter + 50,
          delta = delta,
          epsilon = 0.1,
          K = K
        )

        } else if(type == "mixed"){

          # no control variables allowed
          if(!is.null(formula)){

            warning("formula ignored for mixed data at this time")

            control_info <- remove_predictors_helper(list(as.data.frame(Y)),
                                                     formula = formula)

            # data
            Y <-  as.matrix(control_info$Y_groups[[1]])

            formula <- NULL

          } else {

             Y <- na.omit(Y)

            }

          # default for ranks
          if(is.null(mixed_type)) {

            idx = colMeans(round(Y) == Y)
            idx = ifelse(idx == 1, 1, 0)

            # user defined
            } else {

              idx = mixed_type

              }

          # observations
          n <- nrow(Y)

          # nodes
          p <- ncol(Y)

          # rank following hoff (2008)
          rank_vars <- rank_helper(Y)

      post_samp <- .Call(
        "_BGGM_copula",
        z0_start = rank_vars$z0_start,
        levels = rank_vars$levels,
        K = rank_vars$K,
        Sigma_start = rank_vars$Sigma_start,
        iter = iter + 50,
        delta = delta,
        epsilon = 0.1,
        idx = idx
      )
      } else {

        stop("'type' not supported: must be continuous, binary, ordinal, or mixed.")
      }

    message("BGGM: Finished")

    pcor_mat <- round(apply(post_samp$pcors[,,51:(iter + 50)], 1:2, mean), 3)
    pcor_sd <- round(apply(post_samp$pcors[,,51:(iter + 50)], 1:2, sd), 3)

    results <- list(pcor_mat = pcor_mat,
                            pcor_sd = pcor_sd,
                            analytic = analytic,
                            formula = formula,
                            post_samp = post_samp,
                            type = type,
                            iter = iter,
                            Y = Y,
                            call = match.call(),
                            p = p,
                            n = n)

    #  analytic
    } else {

      if(type != "continuous"){

        warning("analytic solution only available for 'type = continuous'")
        type <- "continuous"

      }
      if(!is.null(formula)){

        stop("formula note permitted with the analytic solution")

      }

    Y <- na.omit(Y)

    # observations
    n <- nrow(Y)

    formula <- NULL

    analytic_fit <- analytic_solve(Y)

    results <- list(pcor_mat = analytic_fit$pcor_mat,
                    analytic_fit = analytic_fit,
                    analytic = analytic,
                    formula = formula,
                    type = type,
                    iter = iter,
                    Y = Y,
                    call = match.call(),
                    p = p,
                    n = n)

    } # end analytic

  returned_object <- results
  class(returned_object) <- c("BGGM", "estimate", "default")
  return(returned_object)

  }

#' @name summary.estimate
#' @title Summary method for \code{estimate.default} objects
#'
#' @param object an object of class \code{estimate}
#'
#' @param col_names logical. Should the summary include the column names (default is \code{TRUE})?
#'                  Setting to \code{FALSE} includes the column numbers (e.g., \code{1--2}).
#'
#' @param cred credible interval width
#' @param ... currently ignored


#' @seealso \code{\link{select.estimate}}
#' @return a list containing the summarized posterior distributions
#' # data
#' Y <- BGGM::bfi[, 1:5]
#' # analytic approach (sample by setting analytic = FALSE)
#' fit <- estimate(Y, analytic = TRUE)
#' summary(fit)
#' @export
summary.estimate <- function(object,
                             col_names = TRUE,
                             cred = 0.95, ...) {

  # nodes
  p <- object$p

  # identity matrix
  I_p <- diag(p)

  # lower bound
  lb <- (1 - cred) / 2

  # upper bound
  ub <- 1 - lb

  # column names
  cn <-  colnames(object$Y)


  if(is.null(cn) | isFALSE(col_names)){

    mat_names <- sapply(1:p , function(x) paste(1:p, x, sep = "--"))[upper.tri(I_p)]

  } else {

    mat_names <-  sapply(cn , function(x) paste(cn, x, sep = "--"))[upper.tri(I_p)]

  }


  if(isFALSE(object$analytic)){

    post_mean <- round(apply( object$post_samp$pcors[,, 51:(object$iter + 50) ], 1:2, mean), 3)[upper.tri(I_p)]
    post_sd  <- round(apply( object$post_samp$pcors[,, 51:(object$iter + 50) ], 1:2, sd), 3)[upper.tri(I_p)]
    post_lb <- round(apply( object$post_samp$pcors[,, 51:(object$iter + 50) ], 1:2, quantile, lb), 3)[upper.tri(I_p)]
    post_ub <- round(apply( object$post_samp$pcors[,, 51:(object$iter + 50) ], 1:2, quantile, ub), 3)[upper.tri(I_p)]



  dat_results <-
    data.frame(
      relation = mat_names,
      post_mean =  post_mean,
      post_sd = post_sd,
      post_lb = post_lb,
      post_ub = post_ub
    )

  colnames(dat_results) <- c(
    "Relation",
    "Post.mean",
    "Post.sd",
    "Cred.lb",
    "Cred.ub")

  } else {

    dat_results <-
      data.frame(
        relation = mat_names,
        post_mean =  object$pcor_mat[upper.tri(I_p)]
      )

    colnames(dat_results) <- c(
      "Relation",
      "Post.mean")

  }


  returned_object <- list(dat_results = dat_results,
                          object = object)


class(returned_object) <- c("BGGM", "estimate",
                            "summary_estimate",
                            "summary.estimate")
returned_object
}


print_summary_estimate <- function(x, ...) {
    cat("BGGM: Bayesian Gaussian Graphical Models \n")
    cat("--- \n")
    cat("Type:",  x$object$type, "\n")
    cat("Analytic:", x$object$analytic, "\n")
    cat("Formula:", paste(as.character(fit$formula), collapse = " "), "\n")
    # number of iterations
    cat("Posterior Samples:", x$object$iter, "\n")
    # number of observations
    cat("Observations (n):\n")
    # number of variables
    cat("Nodes (p):", x$object$p, "\n")
    # number of edges
    cat("Relations:", .5 * (x$object$p * (x$object$p - 1)), "\n")
    cat("--- \n")
    cat("Call: \n")
    print(x$object$call)
    cat("--- \n")
    cat("Estimates:\n")
    print(x$dat_results, row.names = F)
    cat("--- \n")
}


print_estimate <- function(x, ...){
  cat("BGGM: Bayesian Gaussian Graphical Models \n")
  cat("--- \n")
  cat("Type:",  x$type, "\n")
  cat("Analytic:", x$analytic, "\n")
  cat("Formula:", paste(as.character(fit$formula), collapse = " "), "\n")
  # number of iterations
  cat("Posterior Samples:", x$iter, "\n")
  # number of observations
  cat("Observations (n):\n")
  # number of variables
  cat("Nodes (p):", x$p, "\n")
  # number of edges
  cat("Relations:", .5 * (x$p * (x$p-1)), "\n")
  cat("--- \n")
  cat("Call: \n")
  print(x$call)
  cat("--- \n")
  cat("Date:", date(), "\n")
}



#' Plot \code{summary.estimate} Objects
#'
#' @param x an object of class \code{summary.estimate}
#' @param color color of error bar
#' @param width width of error bar cap
#' @param ... currently ignored
#'
#' @return an object of class \code{ggplot}
#' @export
plot.summary.estimate <- function(x, color = "black",
                                  size = 2,
                                  width = 0, ...){

    dat_temp <- x$dat_results[order(x$dat_results$Post.mean,
                                         decreasing = F), ]

    dat_temp$Relation <-
      factor(dat_temp$Relation,
             levels = dat_temp$Relation,
             labels = dat_temp$Relation)


   if(isFALSE(x$object$analytic)){
    ggplot(dat_temp,
           aes(x = Relation,
               y = Post.mean)) +

      geom_errorbar(aes(ymax = dat_temp[, 4],
                        ymin = dat_temp[, 5]),
                    width = width,
                    color = color) +
      geom_point(size = size) +
      xlab("Index") +
      theme(axis.text.x = element_text(
        angle = 90,
        vjust = 0.5,
        hjust = 1
      ))
   } else {

     ggplot(dat_temp,
            aes(x = Relation,
                y = Post.mean)) +
       geom_point(size = size) +
       xlab("Index") +
       theme(axis.text.x = element_text(
         angle = 90,
         vjust = 0.5,
         hjust = 1
       ))
     }
}
