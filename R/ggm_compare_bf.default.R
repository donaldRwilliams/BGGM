#' Compare GGMs with Bayesian Hypothesis Testing
#'
#' @description Compare GGMs with the Bayes factor. This method allows for
#' assessing (relative) evidence for edge equality or edges differences across any number of groups. Further, confirmatory hypothesis testing
#' can be used to test predictions or expectations regarding difference or similarities in different groups (e.g., male vs. female).
#'
#' @param ... data matrices (\emph{n} by  \emph{p}). Requires at least two.
#' @param prior_sd hypothesized standard deviation for the edges or partial correlations.
#' @param iter number of posterior samples.
#' @param cores number of cores for parallel computing. The default is 2, but this can be adjusted.
#'
#' @return list of class \code{ggm_compare_bf}
#' \itemize{
#' \item \code{BF_01} Bayes factors in favor of the null hypothesis
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
#' \item \code{post_samps} matrix of posterior samples
#' \item \code{prior_samps} matrix of prior samples
#' }
#'
#' @note After fitting, use \code{select} to determine which partial correlations were different or the same (i.e., evidence for the null hypothesis of
#' equality). This assumes \code{hypothesis = NULL}. If a hypothesis is tested, then use \code{summary} which provides
#' information including the Bayes factors and posterior probabilities for each hypothesis.
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

