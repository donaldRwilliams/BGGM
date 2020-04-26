#' Compare GGMs: Exploratory Bayesian Hypothesis Testing
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
#' @param prior_sd The scale of the prior distribution (centered at zero), in reference to a beta distribtuion.
#' The `default` is 0.25. See note for further details.
#'
#' @param iter number of iterations (posterior samples; defaults to 5000).
#'
#' @param seed The seed for random number generation (default set to \code{1}).
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
#'
#'
#' \strong{A Note on Defaults:}
#'
#' The `default` for \code{prior_sd} is 0.25. This can and should be changed to reflect the hypothesized difference.
#' If differences are expected to be small the \code{prior_sd} should be set to a smaller value (e.g., 0.15).
#' This might raise concerns of the `correct` prior.` Note, however, that the interpretation remains unchanged: which hypothesis better
#' predicts the observed data. \insertCite{@p. 777, in  @Kass1995}{BGGM}
#'
#' If the desired inference is \emph{not} (relative) evidence between models, the method in \code{\link{ggm_compare_estimate}}
#' (for each partial correlation ) or  \code{\link{ggm_compare_ppc}} (a global test) can be used.
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
ggm_compare_bf <- function(...,
                           prior_sd = 0.20,
                           iter = 5000,
                           seed = 1){

  # group info (e.g., n, p, etc)
  info <- BGGM:::Y_combine(...)

  # number of variables
  p <- info$dat_info$p[1]

  # groups
  groups <- length(info$dat)

  # check at least two groups
  if(groups < 2){
    stop("must have (at least) two groups")
    }

  # matrix-F hyperparameter
  delta <- BGGM:::delta_solve(prior_sd)

  # sample posterior
  post_samp <- lapply(1:groups, function(x) {

    Y <- info$dat[[x]]

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
                          p = p)

    class(returned_object) <- c("BGGM",
                                "ggm_compare_bf",
                                "explore")
    returned_object
}

