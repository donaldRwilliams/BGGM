#' @title GGM: Exploratory Hypothesis Testing
#'
#' @name explore
#'
#' @description Learn the conditional (in)dependence structure with the Bayes factor using the matrix-F
#' prior distribution \insertCite{Mulder2018}{BGGM}. These methods were introduced in
#' \insertCite{Williams2019_bf;textual}{BGGM}. The graph is selected with \code{\link{select.explore}} and
#' then plotted with \code{\link{plot.select}}.
#'
#' @param Y  Matrix (or data frame) of dimensions \emph{n} (observations) by  \emph{p} (variables).
#'
#' @param formula An object of class \code{\link[stats]{formula}}. This allows for including
#' control variables in the model (i.e., \code{~ gender}).
#'
#' @param type Character string. Which type of data for \code{Y} ? The options include \code{continuous},
#' \code{binary}, \code{ordinal}, or \code{mixed} (semi-parametric copula). See the note for further details.
#'
#' @param mixed_type Numeric vector. An indicator of length p for which varibles should be treated as ranks.
#' (1 for rank and 0 to assume normality). The default is to treat all integer variables as ranks
#' when \code{type = "mixed"} and \code{NULL} otherwise. See note for further details.
#'
#' @param analytic Logical. Should the analytic solution be computed (default is \code{FALSE})?
#'                 (currently not implemented)
#'
#' @param prior_sd Scale of the prior distribution, approximately the standard deviation
#'                 of a beta distribution (defaults to sqrt(1/12)).
#'
#' @param iter Number of iterations (posterior samples; defaults to 5000).
#'
#' @param impute Logicial. Should the missing values (\code{NA})
#'               be imputed during model fitting (defaults to \code{TRUE}) ?
#'
#' @param progress Logical. Should a progress bar be included (defaults to \code{TRUE}) ?
#'
#' @param seed An integer for the random seed.
#'
#' @param ... Currently ignored (leave empty).
#'
#' @references
#' \insertAllCited{}
#'
#' @return The returned object of class \code{explore} contains a lot of information that
#'         is used for printing and plotting the results. For users of \strong{BGGM}, the following
#'         are the useful objects:
#'
#' \itemize{
#'
#' \item \code{pcor_mat} partial correltion matrix (posterior mean).
#'
#' \item \code{post_samp} an object containing the posterior samples.
#'
#' }
#'
#' @details
#'
#' \strong{Controlling for Variables}:
#'
#' When controlling for variables, it is assumed that \code{Y} includes \emph{only}
#' the nodes in the GGM and the control variables. Internally, \code{only} the predictors
#' that are included in \code{formula} are removed from \code{Y}. This is not behavior of, say,
#' \code{\link{lm}}, but was adopted to ensure  users do not have to write out each variable that
#' should be included in the GGM. An example is provided below.
#'
#' \strong{Mixed Type}:
#'
#'  The term "mixed" is somewhat of a misnomer, because the method can be used for data including \emph{only}
#'  continuous or \emph{only} discrete variables. This is based on the ranked likelihood which requires sampling
#'  the ranks for each variable (i.e., the data is not merely transformed to ranks). This is computationally
#'  expensive when there are many levels. For example, with continuous data, there are as many ranks
#'  as data points!
#'
#'  The option \code{mixed_type} allows the user to determine  which variable should be treated as ranks
#'  and the "emprical" distribution is used otherwise. This is accomplished by specifying an indicator
#'  vector of length \emph{p}. A one indicates to use the ranks, whereas a zero indicates to "ignore"
#'  that variable. By default all integer variables are handled as ranks.
#'
#' \strong{Dealing with Errors}:
#'
#' An error is most likely to arise when \code{type = "ordinal"}. The are two common errors (although still rare):
#'
#' \itemize{
#'
#' \item The first is due to sampling the thresholds, especially when the data is heavily skewed.
#'       This can result in an ill-defined matrix. If this occurs, we recommend to first try
#'       decreasing \code{prior_sd} (i.e., a more informative prior). If that does not work, then
#'       change the data type to \code{type = mixed} which then estimates a copula GGM
#'       (this method can be used for data containing \strong{only} ordinal variable). This should
#'       work without a problem.
#'
#' \item  The second is due to how the ordinal data are categorized. For example, if the error states
#'        that the index is out of bounds, this indicates that the first category is a zero. This is not allowed, as
#'        the first category must be one. This is addressed by adding one (e.g., \code{Y + 1}) to the data matrix.
#'
#' }
#'
#'
#' \strong{Imputing Missing Values}:
#'
#' Missing values are imputed with the approach described in \insertCite{hoff2009first;textual}{BGGM}.
#' The basic idea is to impute the missing values with the respective posterior pedictive distribution,
#' given the observed data, as the model is being estimated. Note that the default is \code{TRUE},
#' but this ignored when there are no missing values. If set to \code{FALSE}, and there are missing
#' values, list-wise deletion is performed with \code{na.omit}.
#'
#' @note
#'
#' \strong{Posterior Uncertainty}:
#'
#' A key feature of \bold{BGGM} is that there is a posterior distribution for each partial correlation.
#' This readily allows for visiualizing uncertainty in the estimates. This feature works
#' with all data types and is accomplished by plotting the summary of the \code{explore} object
#' (i.e., \code{plot(summary(fit))}). Note that in contrast to \code{estimate} (credible intervals),
#' the posterior standard deviation is plotted for \code{explore} objects.
#'
#'
#' \strong{"Default" Prior}:
#'
#'  In Bayesian statistics, a default Bayes factor needs to have several properties. I refer
#'  interested users to \insertCite{@section 2.2 in @dablander2020default;textual}{BGGM}. In
#'  \insertCite{Williams2019_bf;textual}{BGGM}, some of these propteries were investigated including
#'  model selection consistency. That said, we would not consider this a "default" (or "automatic")
#'  Bayes factor and thus we encourage users to perform sensitivity analyses by varying
#'  the scale of the prior distribution.
#'
#'  Furthermore, it is important to note there is no "correct" prior and, also, there is no need
#'  to entertain the possibility of a "true" model. Rather, the Bayes factor can be interpreted as
#'  which hypothesis best (\strong{relative} to each other) predicts the observed data
#'  \insertCite{@Section 3.2 in @Kass1995}{BGGM}.
#'
#' \strong{Interpretation of Conditional (In)dependence Models for Latent Data}:
#'
#' See \code{\link{BGGM-package}} for details about interpreting GGMs based on latent data
#' (i.e, all data types besides \code{"continuous"})
#'
#'
#' @examples
#' \donttest{
#' # note: iter = 250 for demonstrative purposes
#'
#' ###########################
#' ### example 1:  binary ####
#' ###########################
#' Y <- women_math[1:500,]
#'
#' # fit model
#' fit <- explore(Y, type = "binary",
#'                 iter = 250,
#'                 progress = FALSE)
#'
#' # summarize the partial correlations
#' summ <- summary(fit)
#'
#' # plot the summary
#' plt_summ <- plot(summary(fit))
#'
#' # select the graph
#' E <- select(fit)
#'
#' # plot the selected graph
#' plt_E <- plot(E)
#'
#' plt_E$plt_alt
#'}
#' @export
explore <- function(Y,
                    formula = NULL,
                    type = "continuous",
                    mixed_type = NULL,
                    analytic = FALSE,
                    prior_sd = sqrt(1/12),
                    iter = 5000,
                    progress = TRUE,
                    impute = FALSE,
                    seed = NULL, ...){

  # temporary warning until missing data is fully implemented
  if(!type %in% c("continuous", "mixed")){

    if(impute){

      warning(paste0("imputation during model fitting is\n",
                   "currently only implemented for 'continuous'
                   and 'mixed' data."))
    }
  }

  set.seed(seed)
  ## Random seed unless user provided
  if(!is.null(seed) ) {
    set.seed(seed)
  }

  
  dot_dot_dot <- list(...)

  eps <- 0.01

  # delta parameter
  delta <- delta_solve(prior_sd)

  if(!analytic){

    if(isTRUE(progress)){

      message(paste0("BGGM: Posterior Sampling ", ...))

    }

    # continuous
    if (type == "continuous") {

      # no control
      if (is.null(formula)) {

        # nodes
        p <- ncol(Y)

        if(!impute){

          # na omit
          Y <- as.matrix(na.omit(Y))

          Y_miss <- Y

        } else {

          Y_miss <- ifelse(is.na(Y), 1, 0)

          if(sum(Y_miss) == 0){
            impute <- FALSE
          }

          # impute means
          for(i in 1:p){
            Y[which(is.na(Y[,i])),i] <- mean(na.omit(Y[,i]))
          }

        }

        # scale Y
        Y <- scale(Y, scale = F)

        # NULL design matrix
        X <- NULL

        # observations
        n <- nrow(Y)

        # starting values
        start <- solve(cov(Y))

        # posterior sample
        post_samp <- .Call(
          '_BGGM_Theta_continuous',
          PACKAGE = 'BGGM',
          Y = Y,
          iter = iter + 50,
          delta = delta,
          epsilon = eps,
          prior_only = 0,
          explore = 1,
          start = start,
          progress = progress,
          impute = impute,
          Y_miss = Y_miss
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

        start <- solve(cov(Y))

        # posterior sample
        post_samp <- .Call(
          "_BGGM_mv_continuous",
          Y = Y,
          X = X,
          delta = delta,
          epsilon = eps,
          iter = iter + 50,
          start = start,
          progress = progress
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

        start <- solve(cov(Y))

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

        start <- solve(cov(Y))

      }

      # posterior sample
      post_samp <-  .Call(
        "_BGGM_mv_binary",
        Y = Y,
        X = X,
        delta = delta,
        epsilon = 0.01,
        iter = iter + 50,
        beta_prior = 0.1,
        cutpoints = c(-Inf, 0, Inf),
        start = start,
        progress = progress
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

        # start
        start <- solve(cov(Y))

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

        # start
        start <- solve(cov(Y))

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
      epsilon = 0.01,
      K = K,
      start = start,
      progress = progress
    )

  } else if(type == "mixed"){
    X <- NULL
    # no control variables allowed
    if(!is.null(formula)){

      warning("formula ignored for mixed data at this time")

      control_info <- remove_predictors_helper(list(as.data.frame(Y)),
                                               formula = formula)
      # data
      Y <-  as.matrix(control_info$Y_groups[[1]])

      formula <- NULL

    }

    # observations
    n <- nrow(Y)

    # nodes
    p <- ncol(Y)

    # default for ranks
    if(is.null(mixed_type)) {

      # idx = rep(1, ncol(Y))
      # idx = ifelse(idx == 1, 1, 0)

      idx = rep(1, ncol(Y))

      # user defined
    } else {

      idx = mixed_type
    }

    # rank following hoff (2008)
    rank_vars <- rank_helper(Y)

    if(impute){

      Y_missing <- ifelse(is.na(Y), 1, 0)

      rank_vars$z0_start[is.na(rank_vars$z0_start)] <- rnorm(sum(Y_missing))

      post_samp <- .Call(
        "_BGGM_missing_copula",
        Y = Y,
        Y_missing = Y_missing,
        z0_start = rank_vars$z0_start,
        Sigma_start = cov(rank_vars$z0_start),
        levels = rank_vars$levels,
        iter_missing = iter + 50,
        progress_impute = TRUE,
        K = rank_vars$K,
        idx = idx,
        epsilon = 0.01,
        delta = delta
      )

    } else {

    post_samp <- .Call(
      "_BGGM_copula",
      z0_start = rank_vars$z0_start,
      levels = rank_vars$levels,
      K = rank_vars$K,
      Sigma_start = rank_vars$Sigma_start,
      iter = iter + 50,
      delta = delta,
      epsilon = 0.01,
      idx = idx,
      progress = progress
    )

    }
  } else {

    stop("'type' not supported: must be continuous, binary, ordinal, or mixed.")

    }

    ## matrix dimensions for prior
    ## Old:
    ## Y_dummy <- matrix(rnorm( 10 * 3 ),
    ##                  nrow = 10, ncol = 3)

    ## Replaced with:
    ## 10 times as many rows as columns
    n_row = ncol(Y) * 10
    Y_dummy <- matrix(rnorm( n_row * ncol(Y) ),
                      nrow = n_row, ncol = ncol(Y))
    ## Probably not necessary to scale up Y dim as k=3 was good enough approx.

    if(isTRUE(progress)){

      message(paste0("BGGM: Prior Sampling ", ...))

    }

    # sample prior
    prior_samp <- .Call('_BGGM_sample_prior',
                      PACKAGE = 'BGGM',
                      Y = Y_dummy,
                      iter = iter + 50,
                      delta = delta,
                      epsilon = eps,
                      prior_only = 1,
                      explore = 0,  ## with explore = 0,  k takes number of Y columns instead of k = 3 
                      progress = progress)

    if(isTRUE(progress)){

      message("BGGM: Finished")

    }

    # # compute post.mean
    pcor_mat <- post_samp$pcor_mat

    returned_object <- list(
      pcor_mat = pcor_mat,
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
      X = X,
      eps = eps,
      ppd_mean = post_samp$ppd_mean
    )

  } else {

    stop("analytic solution not currently available")

    }

  # removed per CRAN (8/12/21)
  #.Random.seed <<- old

  returned_object

  class(returned_object) <- c("BGGM",
                              "explore",
                              "default")
  return(returned_object)
}



#' @title Summary Method for \code{explore.default} Objects
#'
#' @description Summarize the posterior distribution for each partial correlation
#' with the posterior mean and standard deviation.
#'
#' @name summary.explore
#'
#' @param object An object of class \code{estimate}
#'
#' @param col_names Logical. Should the summary include the column names (default is \code{TRUE})?
#'                  Setting to \code{FALSE} includes the column numbers (e.g., \code{1--2}).
#'
#' @param ... Currently ignored
#'
#' @seealso \code{\link{select.explore}}
#'
#' @return A dataframe containing the summarized posterior distributions.
#'
#' @examples
#' \donttest{
#' # note: iter = 250 for demonstrative purposes
#'
#' Y <- ptsd[,1:5]
#'
#' fit <- explore(Y, iter = 250,
#'                progress = FALSE)
#'
#' summ <- summary(fit)
#'
#' summ
#' }
#' @export
summary.explore <- function(object,
                             col_names = TRUE, ...) {

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
    post_sd  <- round(apply(object$post_samp$pcors[,, 51:(object$iter + 50) ], 1:2, sd), 3)[upper.tri(I_p)]

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
  cat("Formula:", paste(as.character(x$formula), collapse = " "), "\n")
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


#' @title Plot \code{summary.explore} Objects
#'
#' @description Visualize the posterior distributions for each partial correlation.
#'
#' @name plot.summary.explore
#'
#'
#' @param x An object of class \code{summary.explore}
#'
#' @param size Numeric. The size for the points (defaults to \code{2}).
#'
#' @param color Character string. The color for the error bars.
#' (defaults to \code{"black"}).
#'
#' @param width Numeric. The width of error bar ends (defaults to \code{0} ).
#'
#' @param ... Currently ignored
#'
#' @return A \code{ggplot} object
#'
#' @seealso \code{\link{explore}}
#'
#' @examples
#' \donttest{
#' # note: iter = 250 for demonstrative purposes
#'
#' Y <- ptsd[,1:5]
#'
#' fit <- explore(Y, iter = 250,
#'                progress = FALSE)
#'
#' plt <- plot(summary(fit))
#'
#' plt
#' }
#' @export
plot.summary.explore <- function(x,
                                 color = "black",
                                  size = 2,
                                  width = 0,
                                 ...){

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
