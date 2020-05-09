#' @title GGM: Exploratory Bayesian Hypothesis Testing
#'
#' @name explore
#'
#' @description Learn the conditional (in)dependence structure with the Bayes factor using the matrix-F prior distribution
#' \insertCite{Mulder2018}{BGGM}. It is
#' possible to test for only positive or negative edges, as well as two sided hypothesis testing (which is the customary approach). Further
#' there is also an exhaustive option that provides the posterior probability of the null, greater than zero, and less than zero.
#' These methods were introduced in \insertCite{Williams2019_bf;textual}{BGGM}.
#'
#' @param Y  matrix (or data frame) of dimensions \emph{n} (observations) by  \emph{p} (variables).
#'
#' @param formula an object of class \code{\link[stats]{formula}}. This allows for including
#' control variables in the model (i.e., \code{~ gender}).
#'
#' @param type character string. Which type of data for \strong{Y} ? The options include \code{continuous},
#' \code{binary}, \code{ordinal}, or \code{mixed} (semi-parametric copula). See the note for further details.
#'
#' @param mixed_type numeric vector. An indicator of length p for which varibles should be treated as ranks.
#' (1 for rank and 0 to assume normality). The default is currently (dev version) to treat all integer variables
#' as ranks when \code{type = "mixed"} and \code{NULL} otherwise. See note for further details.
#'
#' @param iter number of iterations (posterior samples; defaults to 5000).
#'
#' @param prior_sd hypothesized standard deviation for the edges or partial correlations
#'
#' @param ... currently ignored.
#'
#' @references
#' \insertAllCited{}
#'
#' @return list of class \code{explore}:
#' \itemize{
#' \item \code{parcors_mat} partial correlation matrix
#' \item \code{parcors_sd} partial correlation standard deviations
#' \item \code{samples} list of prior and posterior samples
#' \itemize{
#'  \item \code{fisher_z_post} Fisher z transformed posterior distributions (partial correlations)
#'  \item  \code{pcor_post} partial correlation posterior distributions (not transformed)
#'  \item \code{inv_cov_post} inverse covariance matrix posterior distribution
#'  \item \code{pcor_prior} partial correlation prior distribution
#'  \item \code{fisher_z_prior}  Fisher z transformed prior distributions (partial correlations)
#'  }
#' \item \code{delta} hyperparameter of matrix-F prior distribution (corresponds to \code{prior_sd})
#' \item \code{iter} number of posterior and prior samples
#' \item \code{dat} data matrix
#' \item \code{call} match.call()
#' \item \code{p} number of variables
#' \item \code{cores} number of cores
#' \item \code{edge} number of estimated edges
#' }
#'
#'
#'
#'
#' @note After sampling from the posterior distribution, use \code{select} to determine the edge set and \code{plot} for visualizing the
#' edge set. see \code{methods(class = "explore")}
#'
#' \strong{Interpretation of conditional (in)dependence models for latent data:}
#'
#' A  tetrachoric correlation (binary data) is a special case of a polychoric correlation (ordinal data). Both relations are
#' between "theorized normally distributed continuous latent variables"
#' (\href{https://en.wikipedia.org/wiki/Polychoric_correlation}{Wikipedia})
#' In both instances, the correpsonding partial correlation between observed variables is conditioned
#' on the remaining variables in the \emph{latent} space. This implies that interpration is much the same as
#' for continuous data, but with respect to latent variables. We refer interested reader to
#' \insertCite{@page 2364, section 2.2, in  @webb2008bayesian;textual}{BGGM}.
#'
#' @examples
#' # p = 10
#' Y <- BGGM::bfi[,1:10]
#'
#' # sample posterior
#' fit <- explore(Y, iter = 500)
#'
#' # select E
#' E <- select(fit, BF_cut = 3)
#'
#' # summarize
#' summary(E)
#'
#' # non-zero edges
#' E$partials_non_zero
#'
#' # adjacency matrix
#' E$Adj_10
#'
#' # null adjacency matrix
#' E$Adj_01
#' @export
explore <- function(Y,
                    formula = NULL,
                    type = "continuous",
                    mixed_type = NULL,
                    analytic = FALSE,
                    prior_sd = 0.25,
                    iter = 5000,...){

  # delta parameter
  delta <- delta_solve(prior_sd)

  if(!analytic){

    message(paste0("BGGM: Posterior Sampling ", ...))

    # continuous
    if (type == "continuous") {

      # no control
      if (is.null(formula)) {

        # na omit
        Y <- as.matrix(na.omit(Y))

        # scale Y
        Y <- scale(Y, scale = F)

        X <- NULL

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

        X <- NULL

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

      } # end control

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

        formula <- ~ 1

      } else {

        control_info <- remove_predictors_helper(list(as.data.frame(Y)),
                                                 formula = formula)

        # data
        Y <-  as.matrix(control_info$Y_groups[[1]])

        X <- NULL

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
        beta_prior = 0.1,
        cutpoints = c(-Inf, 0, Inf)
      )

      # ordinal
    } else if (type == "ordinal") {

      # intercept only
      if(is.null(formula)){

        # data
        Y <- as.matrix(na.omit(Y))

        X <- NULL

        # obervations
        n <- nrow(Y)

        # nodes
        p <- ncol(Y)

        # intercept only
        X <- matrix(1, n, 1)

        formula <- ~ 1

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
    K <- max(apply(Y, 2, function(x) { length(unique(x)) } ))

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

      X = NULL
      }

    # observations
    n <- nrow(Y)

    # nodes
    p <- ncol(Y)

    # default for ranks
    if(is.null(mixed_type)) {

      idx = colMeans(round(Y) == Y)
      idx = ifelse(idx == 1, 1, 0)

      # user defined
    } else {

      idx = mixed_type
    }

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
    # finished sampling
    message("BGGM: Finished")

    # matrix dimensions for prior
    Y_dummy <- matrix(rnorm( 10 * 3 ),
                    nrow = 10, ncol = 3)

    # sample prior
    prior_samp <- .Call('_BGGM_sample_prior',
                      PACKAGE = 'BGGM',
                      Y = Y_dummy,
                      iter = iter,
                      delta = delta,
                      epsilon = 0.1,
                      prior_only = 1,
                      explore = 1)

    # compute post.mean
    pcor_mat <- round(apply(post_samp$pcors[,,51:(iter + 50)], 1:2, mean), 3)

    # compute post.sd
    pcor_sd <- round(apply(post_samp$pcors[,,51:(iter + 50)], 1:2, sd), 3)


  returned_object <- list(
    pcor_mat = pcor_mat,
    pcor_sd = pcor_sd,
    analytic = analytic,
    formula = formula,
    post_samp = post_samp,
    prior_samp = prior_samp,
    type = type,
    iter = iter,
    Y = Y,
    call = match.call(),
    p = p,
    n = n,
    X = X
  )

  }else {

    stop("analytic solution not currently available")

    }

  returned_object

  class(returned_object) <- c("BGGM",
                              "explore",
                              "default")
  return(returned_object)
}




#' @name summary.explore
#'
#' @title Summary method for \code{explore.default} objects
#'
#' @param object an object of class \code{estimate}
#'
#' @param col_names logical. Should the summary include the column names (default is \code{TRUE})?
#'                  Setting to \code{FALSE} includes the column numbers (e.g., \code{1--2}).
#'
#' @param ... currently ignored
#'
#' @seealso \code{\link{select.explore}}
#'
#' @return a list containing the summarized posterior distributions
#' # data
#' Y <- BGGM::bfi[, 1:5]
#' # analytic approach (sample by setting analytic = FALSE)
#' fit <- estimate(Y, analytic = TRUE)
#' summary(fit)
#' @export
summary.explore <- function(object,
                             col_names = TRUE) {

  # nodes
  p <- object$p

  # identity matrix
  I_p <- diag(p)

  # column names
  cn <-  colnames(object$Y)


  if(is.null(cn) | isFALSE(col_names)){

    mat_names <- sapply(1:p , function(x) paste(1:p, x, sep = "--"))[upper.tri(I_p)]

  } else {

    mat_names <-  sapply(cn , function(x) paste(cn, x, sep = "--"))[upper.tri(I_p)]

  }


  if(isFALSE(object$analytic)){

    post_mean <- round(object$pcor_mat, 3)[upper.tri(I_p)]
    post_sd  <- round(object$pcor_sd, 3)[upper.tri(I_p)]

    dat_results <-
      data.frame(
        relation = mat_names,
        post_mean =  post_mean,
        post_sd = post_sd
      )

    colnames(dat_results) <- c(
      "Relation",
      "Post.mean",
      "Post.sd")

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


  class(returned_object) <- c("BGGM", "explore",
                              "summary_explore",
                              "summary.explore")
  returned_object
}



print_explore <- function(x,...){
  cat("BGGM: Bayesian Gaussian Graphical Models \n")
  cat("--- \n")
  cat("Type:",  x$type, "\n")
  cat("Analytic:", x$analytic, "\n")
  cat("Formula:", paste(as.character(fit$formula), collapse = " "), "\n")
  # number of iterations
  cat("Posterior Samples:", x$iter, "\n")
  # number of observations
  cat("Observations (n):", x$n,  "\n")
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


#' Plot \code{summary.explore}
#'
#' @param x an object of class \code{summary.explore}
#' @param ... currently ignored
#'
#' @return an object of class \code{ggplot}
#' @export
plot.summary.explore <- function(x, color = "black",
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

      geom_errorbar(aes(ymax = Post.mean + dat_temp[, 3],
                        ymin = Post.mean -  dat_temp[, 3]),
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
