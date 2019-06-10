#' Compare GGMs with Bayesian Hypothesis Testing
#'
#' @description Comparing any number of GGMs with Bayesian hypothesis testing (i.e., with the Bayes factor). This method allows for
#' assessing (relative) evidence for edge equality or edges differences across any number of groups.
#'
#' @param ... data matrices (\emph{n} \times  \emph{p}). Requires at least two.
#' @param prior_sd hypothesized standard deviation for the edges or partial correlations.
#' @param iter number of posterior samples.
#' @param cores number of cores for parallel computing. The default is 2, but this can be adjusted.
#'
#' @return an object of class \code{ggm_compare_bf}
#' @export
#'
#'
#' @examples
#'
#'# assume null is true
#' Y1 <- MASS::mvrnorm(500, rep(0, 16), Sigma = diag(16))
#' Y2 <- MASS::mvrnorm(500, rep(0, 16), Sigma = diag(16))
#' Y3 <- MASS::mvrnorm(500, rep(0, 16), Sigma = diag(16))
#'
#' bf_ggm <- ggm_compare_bf(Y1, Y2, Y3, prior_sd = .5, iter = 2500, cores = 2)
#' select( bf_ggm, BF_cut = 3)

ggm_compare_bf.default <- function(..., prior_sd, iter = 25000, cores = 2){

  info <- BGGM:::Y_combine(...)

  groups <- length(info$dat)

  delta <- BGGM:::delta_solve(prior_sd)

  fits <- lapply(1:groups, function(x)   BGGM:::sampling(info$dat[[x]], nu = 10000,
                                                         delta = delta,
                                                         n_samples = iter,
                                                         cores = cores))

  # p
  p <- info$dat_info$p[1]

  edges <- 0.5 * (p * (p -1))

  prior_samps <- list()
  post_samps <- list()

  for(i in 1:groups){

    prior_samps[[i]] <- do.call(rbind.data.frame, lapply(1:cores, function(x)   fits[[i]][[x]]$fisher_z_prior[,1:edges]))

    post_samps[[i]] <- do.call(rbind.data.frame, lapply(1:cores, function(x)   fits[[i]][[x]]$fisher_z_post[,1:edges]))

  }

  post_string <- list()
  prior_string <- list()

  for(i in 1:groups){

    post_string[[i]] <-   paste0("post_samps[[", i, "]][,",  1:edges, "]", sep = "")

    prior_string[[i]] <-   paste0("prior_samps[[", i, "]][,",  1:edges, "]", sep = "")

  }

  groups_as_words <- BGGM:::numbers2words(1:groups)

  hyp <- paste( groups_as_words, sep = " ", collapse = "=")

  framed <- BGGM:::framer(hyp)

  mats <- BGGM:::create_matrices(framed = framed, varnames = groups_as_words)


  post_string <- do.call(rbind, post_string)
  prior_string <- do.call(rbind, prior_string)

  BF <- NA
  for(i in 1:edges){

    # temporary string
    temp_string <- paste("cov(cbind(", paste(post_string[,i], sep = " ", collapse = ","), "))")

    # posterior covariance
    cov_post <- eval(parse( text = temp_string))

    # temporary string
    temp_string <- paste("colMeans(cbind(", paste(post_string[,i], sep = " ", collapse = ","), "))")

    # posterior mean
    post_mean <- eval(parse( text = temp_string))

    # temporary string
    temp_string <- paste("cov(cbind(", paste(prior_string[,i], sep = " ", collapse = ","), "))")

    # prior covariance
    cov_prior <- eval(parse( text = temp_string))

    # tranformed posterior
    mu1 <- mats$R_e %*% post_mean
    s1 <- mats$R_e %*% cov_post %*% t(mats$R_e)

    # transformed prior
    mu0 <- mats$R_e %*% rep(0, groups)
    s0 <- mats$R_e %*% cov_prior %*% t(mats$R_e)

    # bayes factor
    log_BF <- mvnfast::dmvn(X = t(mats$r_e), mu = mu1, sigma = s1, log = TRUE) -
      mvnfast::dmvn(X = t(mats$r_e), mu = mu0, sigma = s0, log = TRUE)

    BF[i] <- exp(log_BF)
  }
  BF_01 <- matrix(0, p, p)

  BF_01[upper.tri(BF_01)] <- BF

  BF_01 <- BGGM:::symmteric_mat(BF_01)

  returned_object <- list(BF_01 = BF_01,
                          p = p,
                          info = info,
                          iter = iter,
                          call = match.call(),
                          delta = delta)

  class(returned_object) <- "ggm_compare_bf"
  returned_object
}



