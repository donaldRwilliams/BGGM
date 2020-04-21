#' @title Bayesian Estimation of Gaussian Graphical Models
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
#' control variables in the model (i.e., \code{~ gender}).
#'
#' @param data an optional data frame, list or environment (or an object coercible by \code{\link[base]{as.data.frame}})
#' to a data frame containing the variables in \code{formula}. This is required when controlling for variables.
#'
#' @param type character string. Which type of data for \strong{Y} ? The options include \code{continuous},
#' \code{binary}, \code{ordinal}, or \code{mixed}. See the note for further details.
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
#'
#' @note The default is to draws samples from the posterior distribution (\code{analytic = FALSE}). The samples are
#' required for computing edge differences (see \code{\link{ggm_compare_estimate}}), Bayesian R2 introduced in
#'  \insertCite{gelman_r2_2019;textual}{BGGM} (see \code{\link{bayes_R2}}), etc. If the goal is to *only* determine the non-zero effects, this can be accomplished by setting \code{analytic = TRUE}. Note also sampling is
#' very fast--i.e., less than 1 second with p = 25, n = 2500 and 5,000 samples.
#'
#' These methods are inherently Bayesian. This also means there is a close correspondence to "frequentist" methods or, say,
#' maximum likelihood. In fact, the prior distribution is set to mimic the unbiased estimate of the sample based
#' covariance matrix. This allows for efficiently drawing samples from the posterior distribution. The advantage compared
#' to frequentist methods is that a measure of uncertainty is readily available. This allows for
#' seamlessly computing partial correlation (see \code{\link{ggm_compare_estimate}}) and Bayesian R2  (see \code{\link{test.R2}}) differences.
#' Further, the posterior probability of a null region (see \code{\link{select.estimate}}) can be computed and this provides
#' an estimate of the conditional independence structure (null effects).
#'
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
#'
#' see \code{methods("estimate")}
#' @examples
#' # p = 20
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
                      data = NULL,
                      type = "continuous",
                      iter = 5000,
                      analytic = FALSE, ep = 0.001,...){



  if(type == "continuous"){
  # change to x
  x <- Y

  # remove the NAs
  X <- na.omit(as.matrix(x))

  # mean center the data
  X <- scale(X, scale = T)


  # number of observations
  n <- nrow(X)

  # number of columns
  p <- ncol(X)


  if(is.null(formula)){

  control <- "no_control"

  # scatter matrix
  S <- t(X) %*% X

  # sample from Wishart
  if(isFALSE(analytic)){

    # number of columns inv + pcor
    cols_samps <- p^2 + p^2

    # name the columns
    inv_names <- unlist(lapply(1:p, function(x)  samps_inv_helper(x, p)))
    pcor_names <-unlist(lapply(1:p, function(x)  samps_pcor_helper(x, p)))

    # store posterior samples
    df_samps <- matrix(0, nrow = iter, ncol = cols_samps)

  for(i in 1:iter){

    # draw directly from Wishart
    inv_mat <- rWishart(1, n-1, solve(S))[,,1]

    # compute partial correlations
    pcor_mat <-   -1 * cov2cor(inv_mat)

    # into the $i$th row
    df_samps[i,1:cols_samps] <- c(as.numeric(inv_mat), as.numeric(pcor_mat))

  }

  # name columns
  colnames(df_samps) <- c(inv_names, pcor_names)

  # matrix for storage
  pcor_mat <-  inv_mat <- matrix(0, ncol = p, p)

  # posterior means (partials)
  pcor_mat[] <- colMeans(df_samps[,  grep("pcors", colnames(df_samps))])
  diag(pcor_mat) <- 0

  # posterior means (inverse)
  inv_mat[]   <- colMeans(df_samps[,  grep("cov_inv", colnames(df_samps))])

  returned_object  <- list(pcor_mat = pcor_mat,
                           inv_mat = inv_mat,
                           posterior_samples = as.data.frame(df_samps),
                           p = p, dat = X,
                           iter = iter,
                           call = match.call(),
                           analytic = analytic,
                           control = control,
                           type = type)


  # analytic solution
  } else {

    # analytic
    fit <-  analytic_solve(X)

    returned_object <- list(fit = fit,
                            analytic = analytic,
                            call = match.call(),
                            dat = X,
                            p = p,
                            control = control,
                            type = type)


 }

  } else {

    control <- "control"

    # matrix for storage
    pcor_mat <-  inv_mat <- matrix(0, ncol = p, p)

    # name the columns
    inv_names <- unlist(lapply(1:p, function(x)  samps_inv_helper(x, p)))
    pcor_names <-unlist(lapply(1:p, function(x)  samps_pcor_helper(x, p)))

    X_pred <- model.matrix(formula, as.data.frame(data))

    fit_mvn <- mvn_continuous(X, X_pred,
                              delta = 20,
                                      epsilon = 0.0001,
                                      iter = iter + 50)



    inv_cov <- matrix(as.numeric( fit_mvn$Theta[,,51:(iter+50)]),
                   nrow = iter, ncol = p^2, byrow = TRUE )

    pcors <- matrix(as.numeric( fit_mvn$pcors[,,51:(iter+50)]),
                      nrow = iter, ncol = p^2, byrow = TRUE)
    #
    df_samps <- cbind(inv_cov, pcors)
    #
    colnames(df_samps) <- c(inv_names, pcor_names)

    # posterior means (partials)
    pcor_mat[] <- colMeans(df_samps[,  grep("pcors", colnames(df_samps))])
    diag(pcor_mat) <- 0

    # posterior means (inverse)
    inv_mat[]   <- colMeans(df_samps[,  grep("cov_inv", colnames(df_samps))])

    returned_object  <- list(pcor_mat = pcor_mat,
                             inv_mat = inv_mat,
                             posterior_samples = as.data.frame(df_samps),
                             p = p, dat = X,
                             iter = iter,
                             call = match.call(),
                             analytic = analytic,
                             betas = fit_mvn$beta,
                             coef_names = colnames(X_pred),
                             control = control,
                             type = type)
    }

  } else if(type == "binary"){


    if(is.null(formula)){

      control <- "no_control"

      p <- ncol(Y)

      # matrix for storage
      pcor_mat <-  inv_mat <- matrix(0, ncol = p, p)

      # name the columns
      inv_names <- unlist(lapply(1:p, function(x)  samps_inv_helper(x, p)))
      pcor_names <-unlist(lapply(1:p, function(x)  samps_pcor_helper(x, p)))

      X_pred <- model.matrix(~1, data = as.data.frame( Y))

      fit_mvn <- mvn_binary(Y, X_pred,
                              delta = 20,
                              epsilon = ep,
                              iter = iter + 50,
                          beta_prior = 0.0001,
                          cutpoints = c(-Inf, 0, Inf))



    inv_cov <- matrix(as.numeric( fit_mvn$Theta[,,51:(iter+50)]),
                      nrow = iter, ncol = p^2, byrow = TRUE )

    pcors <- matrix(as.numeric( fit_mvn$pcors[,,51:(iter+50)]),
                    nrow = iter, ncol = p^2, byrow = TRUE)
    #
    df_samps <- cbind(inv_cov, pcors)
    #
    colnames(df_samps) <- c(inv_names, pcor_names)

    # posterior means (partials)
    pcor_mat[] <- colMeans(df_samps[,  grep("pcors", colnames(df_samps))])
    diag(pcor_mat) <- 0

    # posterior means (inverse)
    inv_mat[]   <- colMeans(df_samps[,  grep("cov_inv", colnames(df_samps))])

    returned_object  <- list(pcor_mat = pcor_mat,
                             inv_mat = inv_mat,
                             posterior_samples = as.data.frame(df_samps),
                             p = p, dat = Y,
                             iter = iter,
                             call = match.call(),
                             analytic = analytic,
                             betas = fit_mvn$beta,
                             coef_names = colnames(X_pred),
                             control = control,
                             type = type)

    } # end no control

# end of binary
  } else if(type == "ordinal") {



    if(is.null(formula)){

      control <- "no_control"

      p <- ncol(Y)

      # matrix for storage
      pcor_mat <-  inv_mat <- matrix(0, ncol = p, p)

      # name the columns
      inv_names <- unlist(lapply(1:p, function(x)  samps_inv_helper(x, p)))
      pcor_names <-unlist(lapply(1:p, function(x)  samps_pcor_helper(x, p)))

      X_pred <- model.matrix(~1, data = as.data.frame( Y))

      fit_mvn <- mvn_ordinal(Y, X_pred,
                            delta = 20,
                            epsilon = ep,
                            iter = iter + 50,
                            MH = 0.001)



      inv_cov <- matrix(as.numeric( fit_mvn$Theta[,,51:(iter+50)]),
                        nrow = iter, ncol = p^2, byrow = TRUE )

      pcors <- matrix(as.numeric( fit_mvn$pcors[,,51:(iter+50)]),
                      nrow = iter, ncol = p^2, byrow = TRUE)
      #
      df_samps <- cbind(inv_cov, pcors)
      #
      colnames(df_samps) <- c(inv_names, pcor_names)

      # posterior means (partials)
      pcor_mat[] <- colMeans(df_samps[,  grep("pcors", colnames(df_samps))])
      diag(pcor_mat) <- 0

      # posterior means (inverse)
      inv_mat[]   <- colMeans(df_samps[,  grep("cov_inv", colnames(df_samps))])

      returned_object  <- list(pcor_mat = pcor_mat,
                               inv_mat = inv_mat,
                               posterior_samples = as.data.frame(df_samps),
                               p = p, dat = Y,
                               iter = iter,
                               call = match.call(),
                               analytic = analytic,
                               betas = fit_mvn$beta,
                               coef_names = colnames(X_pred),
                               control = control,
                               type = type)

    } # end no control


}



  class(returned_object) <- c("BGGM", "estimate", "default")
  return(returned_object)


}

#' @name summary.estimate
#' @title Summary method for \code{estimate.default} objects
#'
#' @param object An object of class \code{estimate}
#' @seealso \code{\link{select.estimate}}
#' @param ... currently ignored
#' @param cred credible interval width
#' @return A list containing the summarized posterior distributions
#' # data
#' Y <- BGGM::bfi[, 1:5]
#' # analytic approach (sample by setting analytic = FALSE)
#' fit <- estimate(Y, analytic = TRUE)
#' summary(fit)
#' @export
summary.estimate <- function(object, cred = 0.95, ...) {

  if (isTRUE(object$analytic)) {
    returned_object <- list(object = object)

    } else {

    lb <- (1 - cred) / 2
    ub <- 1 - lb

    name_temp <- matrix(0, object$p, object$p)

    edge_names <- unlist(lapply(1:object$p , function(x)
      paste(1:object$p, x, sep = "--")))

    name_temp[] <- edge_names
    up_tri <- name_temp[upper.tri(name_temp)]

    pcor_samples <-
      object$posterior_samples[,  grep("pcors", colnames(object$posterior_samples))]

    colnames(pcor_samples) <- edge_names
    pcor_upper <- pcor_samples[, up_tri]

    ci <- apply(
      pcor_upper,
      MARGIN = 2,
      FUN = function(x) {
        quantile(x, probs = c(lb, ub))
      }
    )
    diff_mu <-
      apply(pcor_upper, MARGIN = 2, mean)

    diff_sd <-
      apply(pcor_upper, MARGIN = 2, sd)

    dat_results <-
      data.frame(
        edge = name_temp[upper.tri(name_temp)],
        post_mean =  round(diff_mu, 3),
        post_sd = round(diff_sd, 3),
        ci = round(t(ci), 3)
      )

    colnames(dat_results) <- c(
      "Edge",
      "Estimate",
      "Est.Error",
      "Cred.lb", "Cred.ub")

    returned_object <- list(dat_results = dat_results,
                            object = object,
                            pcor_samples = pcor_samples)
  }

  class(returned_object) <- c("BGGM", "estimate", "summary.estimate")
  returned_object
}



#' Plot \code{summary.estimate}
#'
#' @param x an object of class \code{summary.estimate}
#' @param color color of error bar
#' @param width width of error bar cap
#' @param ... currently ignored
#'
#' @return an object of class \code{ggplot}
#' @export
plot.summary.estimate <- function(x, color = "black", width = 0,...){

  dat_temp <- x$dat_results[order(x$dat_results$Estimate,
                                  decreasing = F), ]

  dat_temp$Edge <-
    factor(dat_temp$Edge,
           levels = dat_temp$Edge,
           labels = dat_temp$Edge)

  dat_temp$selected <-
    as.factor(ifelse(dat_temp[, 4] < 0 & dat_temp[, 5] > 0, 0, 1))

  plt <- ggplot(dat_temp,
                aes(x = Edge,
                    y = Estimate,
                    color = selected)) +

    geom_errorbar(aes(ymax = dat_temp[, 4],
                      ymin = dat_temp[, 5]),
                  width = width,
                  color = color) +
    geom_point() +
    xlab("Index")

  return(plt)
 }

