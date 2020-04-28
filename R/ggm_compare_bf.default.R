#' GGM Compare: Exploratory Hypothesis Testing
#'
#' @description Compare Gaussian graphical models with (exploratory) Bayesian hypothesis testing using the matrix-F prior
#' distribution \insertCite{Mulder2018}{BGGM}. A test for each partial correlation in the model for any number of groups.
#' This provides evidence for the null hypothesis of no difference and the alternative hypothesis
#' of difference. With more than two groups, the test is for \emph{all} groups simultaneously (i.e., the relation is the same or different in all groups).
#' This method was introduced in \insertCite{williams2020comparing;textual}{BGGM}. For confirmatory hypothesis testing see
#' \code{confirm_groups}.
#'
#' @param ... at least two matrices (or data frame) of dimensions \emph{n} (observations) by  \emph{p} (variables).
#'
#' @param formula an object of class \code{\link[stats]{formula}}. This allows for including
#' control variables in the model (i.e., \code{~ gender}).
#'
#' @param data an optional data frame, list or environment (or an object coercible by \code{\link[base]{as.data.frame}})
#' to a data frame containing the variables in \code{formula}. This is required when controlling for variables.
#'
#' @param prior_sd The scale of the prior distribution (centered at zero), in reference to a beta distribtuion.
#' The `default` is 0.25. See note for further details.
#'
#' @param type character string. Which type of data for \strong{Y} ? The options include \code{continuous},
#' \code{binary}, or \code{ordinal}. See the note for further details.
#'
#' @param mixed_type numeric vector. An indicator of length p for which varibles should be treated as ranks.
#' (1 for rank and 0 to assume normality). The default is currently (dev version) to treat all integer variables
#' as ranks when \code{type = "mixed"} and \code{NULL} otherwise. See note for further details.
#'
#' @param analytic logical. Should the analytic solution be computed (default is \code{FALSE}) ? See note for details.
#'
#' @param iter number of iterations (posterior samples; defaults to 5000).
#'
#' @param seed The seed for random number generation (default set to \code{1}).
#'
#' @references
#' \insertAllCited{}
#'
#'
#' @return list of class \code{ggm_compare_bf}
#' \itemize{
#' \item \code{BF_01}  Bayes factors in favor of the null hypothesis.
#' \item \code{pcor_diff} partial correlation difference (for two groups only)
#' \item \code{p} number of variables
#' \item  \code{info} list of information about the data matrices
#'
#' \itemize{
#' \item \code{dat} list containing the data matrices
#' \item \code{dat_info} sample size for each data matrix
#' \item \code{pairwise} matrix of pairwise combinations
#' }
#' \item \code{iter} number of posterior and prior samples
#' \item \code{call} match.call()
#' \item \code{delta} hyperparameter of matrix-F prior distribution (corresponds to \code{prior_sd})
#'
#' \item \code{groups} number of groups
#' \item \code{post_samp} matrix of posterior samples
#' }
#'
#' @note
#'
#' \strong{More than Two Groups:}
#'
#'  It is possible to compare any number of group. This does \emph{not} provide pairwise differences. Rather, the test
#'  is for all group, say, that each partial correlation is the same in each. This is analagous to an ANOVA--an
#'  `omnibus` test that does not indicate which groups are different. A follow up test would be required for this purpose.
#'
#' \strong{A Note on Defaults:}
#'
#' The `default` for \code{prior_sd} is 0.25. This can and should be changed to reflect the hypothesized difference.
#' If differences are expected to be small \code{prior_sd} should be set to a smaller value (e.g., 0.15).
#' This might raise concerns of the `correct` prior.` Note, however, that the interpretation remains unchanged: which hypothesis better
#' predicts the observed data. \insertCite{@p. 777, in  @Kass1995}{BGGM}
#'
#' If the desired inference is \emph{not} (relative) evidence between models, the method in \code{\link{ggm_compare_estimate}}
#' (for each partial correlation ) or  \code{\link{ggm_compare_ppc}} (a global test) can be used.
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
#' \strong{Mixed Data:}
#'
#'
#' @export
#'
#'
#' @examples
#' \donttest{
#'# assume null is true
#' Y1 <- MASS::mvrnorm(500, rep(0, 16), Sigma = diag(16))
#' Y2 <- MASS::mvrnorm(500, rep(0, 16), Sigma = diag(16))
#' Y3 <- MASS::mvrnorm(500, rep(0, 16), Sigma = diag(16))
#'
#' # fit model
#' bf_ggm <- ggm_compare_bf(Y1, Y2, Y3,
#'                          prior_sd = .2,
#'                          iter = 500, cores = 2)
#'
#' # select graph
#' sel <- select(bf_ggm, BF_cut = 3)
#'
#' # plot
#' plot(sel)
#'
#' }
ggm_compare_explore <- function(...,
                           formula = NULL,
                           data = NULL,
                           type = "continuous",
                           mixed_type = NULL,
                           analytic = FALSE,
                           prior_sd = 0.20,
                           iter = 5000,
                           seed = 1){


  # set seed
  set.seed(seed)

  # group info (e.g., n, p, etc)
  info <- BGGM:::Y_combine(...)

  # number of variables
  p <- info$dat_info$p[1]

  # number of observation
  n = info$dat_info$n[1]

  # groups
  groups <- length(info$dat)

  # check at least two groups
  if(groups < 2){
    stop("must have (at least) two groups")
    }

  # matrix-F hyperparameter
  delta <- BGGM:::delta_solve(prior_sd)

  if(!analytic){

  if(type == "continuous"){

      if(is.null(formula)){

      # sample posterior
  post_samp <- lapply(1:groups, function(x) {

    Y <- as.matrix(scale(info$dat[[x]], scale = FALSE))

    .Call(
      '_BGGM_Theta_continuous',
      PACKAGE = 'BGGM',
      Y = Y,
      iter = iter + 50,
      delta = delta,
      epsilon = 0.1,
      prior_only = 0,
      explore = 1
    )
  })

  # control for variables
  } else {

    # model matrix
    X <- model.matrix(formula, data)

    # posterior sample
    post_samp <- lapply(1:groups, function(x) {

      # data Y_gj
      Y <- as.matrix(scale(info$dat[[x]], scale = F))

      .Call(
        "_BGGM_mv_continuous",
        Y = Y,
        X = X,
        delta = delta,
        epsilon = 0.1,
        iter = iter + 50
      )
    })
  }

  } else if (type == "binary") {


    # intercept only
    if(is.null(formula)){

      X <- matrix(1, n, 1)

      # model matrix
    } else {

      X <- model.matrix(formula, data)

    }


    # posterior sample
    post_samp <- lapply(1:groups, function(x) {

      # Y_gj
      Y <- as.matrix(info$dat[[x]])

      .Call(
        "_BGGM_mv_binary",
        Y = Y,
        X = X,
        delta = delta,
        epsilon = 0.1,
        iter = iter + 50,
        beta_prior = 0.0001,
        cutpoints = c(-Inf, 0, Inf)
      )
    })

  } else if(type == "ordinal"){

    # intercept only
    if(is.null(formula)){

      X <- matrix(1, n, 1)

    } else {

      # model matrix
      X <- model.matrix(formula, data)

    }

    # posterior sample
    post_samp <- lapply(1:groups, function(x) {

      # Y_gj
      Y <- as.matrix(info$dat[[x]])

      # categories
      K <- max(apply(Y, 2, function(x) { length(unique(x))   } ))

      # call c ++
      .Call("_BGGM_mv_ordinal_albert",
            Y = Y,
            X = X,
            iter = iter + 50,
            delta = delta,
            epsilon = 0.1,
            K = K)
    })

    } else if (type == "mixed") {

    # no control variables allowed
    if(!is.null(formula)){

      warning("formula ignored for mixed data at this time")
      formula <- NULL
    }

    # posterior samples
    post_samp <- lapply(1:groups, function(x) {

      # Y_gj
      Y <- as.matrix(info$dat[[x]])

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

      .Call("_BGGM_copula",
            z0_start = rank_vars$z0_start,
            levels = rank_vars$levels,
            K = rank_vars$K,
            Sigma_start = rank_vars$Sigma_start,
            iter = iter + 50,
            delta = delta,
            epsilon = 0.1,
            idx = idx)
      })

    } else {

      stop("'type' not supported: must be continuous, binary, ordinal, or mixed.")

      }




  # matrix dimensions for prior
  Y_dummy <- matrix(rnorm( groups * (groups+1) ),
                    nrow = groups+1, ncol = groups)

  # sample prior
  prior_samp <- lapply(1:groups, function(x) {

    # Y <- info$dat[[x]]

    .Call(
      '_BGGM_Theta_continuous',
      PACKAGE = 'BGGM',
      Y = Y_dummy,
      iter = 10000,
      delta = delta,
      epsilon = 0.1,
      prior_only = 1,
      explore = 0
    )$fisher_z
  })





  # store pcor diff
  pcor_diff <- BF_01_mat <- matrix(0, p, p)

  # upper triangular elements
  indices <- which(upper.tri(diag(p)), arr.ind = TRUE )

  # make contrast matrices
  ## words for compatability
  groups_as_words <- numbers2words(1:groups)

  ## hypotheses
  hyp <- paste(groups_as_words, sep = " ", collapse = "=")

  ## `framed` hypotheses
  framed <- framer(hyp)

  ## contrast matrices
  mats <- create_matrices(framed = framed,
                          varnames = groups_as_words)


  # loop through upper triangular
  for(i in seq_len(nrow(indices))){

    rho_ij <- indices[i,]

    # start
    post_group <-  post_samp[[1]]$fisher_z[ rho_ij[1], rho_ij[2], (51:(iter + 50))]
    prior_group <-  prior_samp[[1]][ 1, 2,]

    # combined groups
    for(j in 2:(groups)){
      post_group <-  cbind(post_group,  post_samp[[j]]$fisher_z[ rho_ij[1], rho_ij[2], (51:(iter + 50))])
      prior_group <-  cbind(prior_group,  prior_samp[[j]][1, 2,])
    }

    # posterior covariance
    cov_post <- cov(post_group)

    # prior covariance
    cov_prior <- cov(prior_group)

    # posterior mean
    post_mean <- colMeans(post_group)

    # tranformed posterior
    mu_post <- mats$R_e %*% post_mean
    s_post <- mats$R_e %*% cov_post %*% t(mats$R_e)

    # transformed prior
    mu_prior <- mats$R_e %*% rep(0, groups)
    s_prior <- mats$R_e %*% cov_prior %*% t(mats$R_e)

    # bayes factor
    log_BF <- mvtnorm::dmvnorm(x = t(mats$r_e),
                                 mean = mu_post,
                                 sigma = s_post,
                                 log = TRUE) -
                mvtnorm::dmvnorm(x = t(mats$r_e),
                                 mean = mu_prior,
                                 sigma = s_prior,
                                 log = TRUE)

    BF_01_mat[ rho_ij[1], rho_ij[2] ] <- exp(log_BF)

    if(groups == 2){
      pcor_diff[ rho_ij[1], rho_ij[2] ] <-  (z2r(post_mean)[1] - z2r(post_mean)[2])
        }
  }

  BF_01 <- BGGM:::symmteric_mat(BF_01_mat)
  pcor_diff <- BGGM:::symmteric_mat(pcor_diff)

  returned_object <- list(BF_01 = BF_01,
                          info = info,
                          iter = iter,
                          prior_sd = prior_sd,
                          call = match.call(),
                          delta = delta,
                          groups = groups,
                          pcor_diff = pcor_diff,
                          post_samp = post_samp,
                          type = type,
                          p = p)

  # analytic solution
  } else {

    stop("analytic not currently implemented")

  }



    class(returned_object) <- c("BGGM",
                                "ggm_compare_explore",
                                "explore")
    returned_object
}



print_summary_ggm_compare_bf <- function(x, ...){
  groups <- x$object$groups
  cat("BGGM: Bayesian Gaussian Graphical Models \n")
  cat("--- \n")
  cat("Type:",  x$object$type, "\n")
  # number of iterations
  cat("Posterior Samples:", x$object$iter, "\n")
  # number of observations
  cat("Observations (n):\n")
  groups <- length(x$object$info$dat)
  for(i in 1:groups){
    cat("  Group", paste( i, ":", sep = "") , x$object$info$dat_info$n[[i]], "\n")
  }
  # number of variables
  cat("Variables (p):", x$object$p, "\n")
  # number of edges
  cat("Relations:", .5 * (x$object$p * (x$object$p-1)), "\n")
  cat("Delta:", x$object$delta, "\n")
  cat("--- \n")
  cat("Call: \n")
  print(x$object$call)
  cat("--- \n")
  cat("Hypotheses:\n")
  cat("H0:", paste0("rho_g", 1:groups, collapse = " = "), "\n")
  cat("H1:", paste0("rho_g", 1:groups, collapse = " - "), " = 0\n")
  cat("--- \n\n")

  print(x$results, right = FALSE, row.names = FALSE)
  cat("--- \n")
}



print_ggm_compare_bf <- function(x, ...){
  cat("BGGM: Bayesian Gaussian Graphical Models \n")
  cat("--- \n")
  cat("Type:",  x$type, "\n")
  # number of iterations
  cat("Posterior Samples:", x$iter, "\n")
  # number of observations
  cat("Observations (n):\n")
  groups <- length(x$info$dat)
  for(i in 1:groups){
    cat("  Group", paste( i, ":", sep = "") , x$info$dat_info$n[[i]], "\n")
  }
  # number of variables
  cat("Variables (p):", x$p, "\n")
  # number of edges
  cat("Relations:", .5 * (x$p * (x$p-1)), "\n")
  cat("Delta:", x$delta, "\n")
  cat("--- \n")
  cat("Call: \n")
  print(x$call)
  cat("--- \n")
  cat("Hypotheses:\n")
  cat("H0:", paste0("rho_g", 1:groups, collapse = " = "), "\n")
  cat("H1:", paste0("rho_g", 1:groups, collapse = " - "), " = 0\n")
  cat("--- \n")
  cat("Date:", date(), "\n")
}




#' Summary Method for \code{ggm_compare_explore} Objects
#'
#' @description Summarize the posterior hypothesis probabilities
#'
#' @param object object of class \code{ggm_compare_explore}
#' @param col_names logical. Should the summary include the column names (default is \code{TRUE})?
#'                  Setting to \code{FALSE} includes the column numbers (e.g., \code{1--2}).
#'
#' @param ... currently ignored
#'
#' @return
#' @export
#'
#' @examples
summary.ggm_compare_explore <- function(object,
                                        col_names = TRUE,
                                        ...){

  # nodes
  p <- object$p

  # identity matrix
  I_p <- diag(p)

  # prob null
  prob_H0 <-   round(object$BF_01 / (object$BF_01 + 1), 3)

  # prob h1
  prob_H1 <-   round(1 - prob_H0, 3)

  # column names
  cn <-  colnames(object$info$dat[[1]])


  if(!isTRUE(col_names) | is.null(cn)){

    mat_names <- sapply(1:p , function(x) paste(1:p, x, sep = "--"))[upper.tri(I_p)]

  } else {


    mat_names <-  sapply(cn , function(x) paste(cn, x, sep = "--"))[upper.tri(I_p)]

  }



  if(object$groups == 2){
  post_mean <- round(object$pcor_diff[upper.tri(I_p)], 3)

  post_sd <- round(apply(object$post_samp[[1]]$pcors -
                          object$post_samp[[2]]$pcors, 1:2, sd)[upper.tri(I_p)], 3)


  results <- data.frame(Relation = mat_names,
                                Post.mean = post_mean,
                                Post.sd = post_sd,
                                Pr.H0 = prob_H0[upper.tri(I_p)],
                                Pr.H1 = prob_H1[upper.tri(I_p)])


  } else {

    results <- data.frame(Relation = mat_names,
                          Pr.H0 = prob_H0[upper.tri(I_p)],
                          Pr.H1 = prob_H1[upper.tri(I_p)])


  }
  returned_object <- list(results = results,
                          object = object)

  class(returned_object) <- c("BGGM",
                              "ggm_compare_explore",
                              "summary.ggm_compare_explore",
                              "explore")
  returned_object
}



#' GGM Compare: Plot \code{ggm_compare_explore] Object}
#'
#' @description Visualize posterior hypothesis probabilities.
#'
#' @param x object of class \code{ggm_compare_explore}
#' @param size numeric. Size of the points.
#' @param color string. Color of the points.
#'
#' @return
#' @export
#'
plot.summary.ggm_compare_explore <- function(x, size = 2, color = "black"){


  dat_temp <- x$results[order(x$results$Pr.H1,
                              decreasing = F), ]

  dat_temp$Relation <-
    factor(dat_temp$Relation,
           levels = dat_temp$Relation,
           labels = dat_temp$Relation)


  ggplot(dat_temp,
         aes(x = Relation,
             y = Pr.H1)) +
    geom_point(size = size, color = color) +

    theme(axis.text.x = element_text(
      angle = 90,
      vjust = 0.5,
      hjust = 1
    )) +
    coord_flip()

}


