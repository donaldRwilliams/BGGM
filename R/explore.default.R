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
#' @param data an optional data frame, list or environment (or an object coercible by \code{\link[base]{as.data.frame}})
#' to a data frame containing the variables in \code{formula}. This is required when controlling for variables.
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
                    data = NULL,
                    type = "continuous",
                    mixed_type = NULL,
                    analytic = FALSE,
                    prior_sd = 0.25,
                    iter = 5000,...){

  # delta parameter
  delta <- delta_solve(prior_sd)

  # remove NAs
  Y <- na.omit(Y)

  # number of columns
  p <- ncol(Y)

  # number of rows
  n <- nrow(Y)

  if(!analytic){

  message("BGGM: Posterior Sampling")
  # continuous
  if(type == "continuous"){

    # scale Y
    Y <- scale(Y, scale = T)

    # no control
    if(is.null(formula)){

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

      # model matrix
      X <- model.matrix(formula, data)

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
  } else if (type == "binary"){


    # intercept only
    if(is.null(formula)){
      X <- matrix(1, n, 1)
      # model matrix
    } else {

      X <- model.matrix(formula, data)

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
  } else if(type == "ordinal"){

    # intercept only
    if(is.null(formula)){
      X <- matrix(1, n, 1)

    } else {

      # model matrix
      X <- model.matrix(formula, data)

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
      formula <- NULL
    }

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
  message("BGGM: Finished")

  # matrix dimensions for prior
  Y_dummy <- matrix(rnorm( 10 * 3 ),
                    nrow = 10, ncol = 3)

  prior_samp <- .Call('_BGGM_Theta_continuous',
                      PACKAGE = 'BGGM',
                      Y = Y_dummy,
                      iter = iter,
                      delta = delta,
                      epsilon = 0.1,
                      prior_only = 1,
                      explore = 1)


  pcor_mat <- round(apply(post_samp$pcors[,,51:(iter + 50)], 1:2, mean), 3)
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
    n = n
  )
  }else {
    stop("analytic solution not currently available")
    }

  returned_object
  class(returned_object) <- c("BGGM", "default",
  "explore")
  return(returned_object)
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



#' @name summary.explore
#' @title  Summary method for \code{explore} objects
#' @param object An object of class \code{explore}
#' @seealso \code{\link{explore}}
#' @param ... currently ignored
#' @export
summary.explore <- function(object,...){

  dat_results <- select(object,
                      BF_cut = 0,
                      alternative =
                      "exhaustive")

returned_object <- list(dat_results = dat_results,
                        object = object)

class(returned_object) <- c("BGGM",
                           "explore",
                           "summary.explore")
return(returned_object)
}


#' Plot \code{summary.explore}
#'
#' @param x an object of class \code{summary.explore}
#' @param ... currently ignored
#'
#' @return an object of class \code{ggplot}
#' @export

plot.summary.explore <- function(x,...){

  not_h0 <-  1 - (x$dat_results$post_prob$prob_zero)
  h0 <- (x$dat_results$post_prob$prob_zero)

  test_dat <- data.frame(Edge = x$dat_results$post_prob$edge, probability = c(not_h0 - 0.5))

  dat_temp <- test_dat[order(test_dat$probability, decreasing = T),]
  dat_temp$Edge <-
    factor(dat_temp$Edge,
           levels = rev(dat_temp$Edge),
           labels = rev(dat_temp$Edge))


  ggplot(dat_temp, aes(x = Edge, y = probability)) +
    geom_linerange(aes(x = Edge, ymin = 0, ymax = probability),
                   position = position_dodge(width = 1)) +

    scale_y_continuous(limits = c(-0.5, 0.5),
                       labels = c(1, 0.5, 0, 0.5, 1)) +
    geom_point(aes(x = Edge, y = probability),
               position = position_dodge(width = 1)) +

    ylab(
      expression(
        italic(H)[0] * symbol(' \254 ') * "Posterior Probability " * symbol('\256 ') *
          "'not " * italic(H)[0] * "'"
      )
    ) +
    geom_point(color = "blue", size = 1) +
    geom_hline(yintercept = 0, linetype = "dotted") +
    coord_flip()
}
