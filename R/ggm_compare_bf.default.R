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
#'
#'# assume null is true
#' Y1 <- MASS::mvrnorm(500, rep(0, 16), Sigma = diag(16))
#' Y2 <- MASS::mvrnorm(500, rep(0, 16), Sigma = diag(16))
#' Y3 <- MASS::mvrnorm(500, rep(0, 16), Sigma = diag(16))
#'
#' bf_ggm <- ggm_compare_bf(Y1, Y2, Y3, prior_sd = .5, iter = 500, cores = 2)

ggm_compare_bf <- function(..., prior_sd = 0.35, iter = 25000, cores = 2){

  priorprob <- 1

  info <- Y_combine(...)

  groups <- length(info$dat)

  if(groups < 2){
    stop("must have (at least) two groups")
  }

  delta <- delta_solve(prior_sd)

  fits <- lapply(1:groups, function(x)   sampling(info$dat[[x]], nu = 10000,
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

    post_samps[[i]] <- do.call(rbind.data.frame,  lapply(1:cores, function(x)   fits[[i]][[x]]$fisher_z_post[,1:edges]))

  }




    post_string <- list()
    prior_string <- list()



    mu_diff <- list()
    sd_diff <- list()

    for(i in 1:groups){

      post_string[[i]] <-   paste0("post_samps[[", i, "]][,",  1:edges, "]", sep = "")

      prior_string[[i]] <-   paste0("prior_samps[[", i, "]][,",  1:edges, "]", sep = "")

    }

    groups_as_words <- numbers2words(1:groups)

    hyp <- paste(groups_as_words, sep = " ", collapse = "=")

    framed <- framer(hyp)

    mats <- create_matrices(framed = framed, varnames = groups_as_words)


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

      if(groups == 2){
        mu_diff[[i]] <-  z2r(post_mean)[1] - z2r(post_mean)[2]
        sd_diff[[i]] <- sqrt(s1)
        }

    }
    BF_01 <- matrix(0, p, p)

    BF_01[upper.tri(BF_01)] <- BF

    BF_01 <- symmteric_mat(BF_01)


    returned_object <- list(BF_01 = BF_01,
                            p = p,
                            info = info,
                            iter = iter,
                            prior_sd = prior_sd,
                            call = match.call(),
                            delta = delta,
                            groups = groups,
                            mu_diff = mu_diff,
                            post_samps = post_samps,
                            prior_samps = prior_samps,
                            re = mats$R_e)

    class(returned_object) <- "ggm_compare_bf"
    returned_object
}

#' @name print.ggm_compare_bf
#' @title  Print method for \code{ggm_compare_bf} objects
#'
#' @param x An object of class \code{ggm_compare_bf}
#' @param ... currently ignored
#' @seealso \code{\link{ggm_compare_bf.default}}
#' @export
print.ggm_compare_bf <- function(x,...) {
  cat("BGGM: Bayesian Gaussian Graphical Models \n")
  cat("--- \n")
  cat("Type: GGM Compare with Bayesian Hypothesis Testing \n")
  p <- x$info$dat_info$p[1]
  cat("Posterior Samples:", x$iter, "\n")
  cat("Observations: \n")
  groups <- length(x$info$dat)
  for(i in 1:groups){
    cat("  Group", paste( i, ":", sep = "") , x$info$dat_info$n[[i]], "\n")
  }
  cat("Variables (p):", p, "\n")
  cat("Edges:", .5 * (p * (p-1)), "\n")
  cat("Groups:", nrow(x$info$dat_info), "\n")
  cat("Prior SD:", x$prior_sd, "\n")
  cat("--- \n")
  if(is.null(x$hypotheses)){
    cat("Call: \n")
    print(x$call)
    cat("--- \n")
    cat("Date:", date(), "\n")
  }
  if(isFALSE(is.null(x$hypotheses))){
    cat("Date:", date(), "\n")
  }
}

#' @name summary.ggm_compare_bf
#' @title Summary method for \code{ggm_compare_bf.default} objects
#'
#' @param object An object of class \code{ggm_compare_bf}
#' @param ... currently ignored
#' @seealso \code{\link{ggm_compare_bf.default}}
#' @return A list containing the summarized posterior hypothesis probabilities
#' @export
summary.ggm_compare_bf <- function(object, ...) {
  name_temp <- matrix(0, object$p, object$p)

  name_temp[] <-
    unlist(lapply(1:object$p , function(x)
      paste(1:object$p, x, sep = "--")))

  edge_name <- name_temp[upper.tri(name_temp)]

  groups <- length(object$info$dat)

  BF_10 <- 1 /  object$BF_01[upper.tri(object$BF_01)]
  prob_H1 <- BF_10 / (1 + BF_10)
  prob_H0 <- 1 - prob_H1

  if (groups == 2) {
    pcor_diff <- apply(object$post_samps[[1]], 2,  BGGM::fisher_z2r) -
      apply(object$post_samps[[2]], 2,  BGGM::fisher_z2r)

    sd_diff <- apply(pcor_diff, 2, sd)
    results <- data.frame(
      Edge = edge_name,
      Estimate = round(unlist(object$mu_diff), 3),
      Est.Error = round(sd_diff, 3),

      Pr.H0 = round(prob_H0, 3),
      Pr.H1 = round(prob_H1, 3)
    )


  } else {
    results <- data.frame(
      Edge = edge_name,

      Pr.H0 = round(prob_H0, 3),
      Pr.H1 = round(prob_H1, 3)
    )

  }
  returned_object <- list(results = results,
                          object = object)
  class(returned_object) <- "summary.ggm_compare_bf"
  return(returned_object)
}


#' @title Summary method for \code{summary.ggm_compare_bf} objects
#' @name print.summary.ggm_estimate_bf

#' @param x An object of class \code{summary.ggm_compare_bf}
#' @param ... currently ignored
#' @seealso \code{\link{ggm_compare_bf.default}}
#' @export
print.summary.ggm_compare_bf <- function(x, ...) {
  cat("BGGM: Bayesian Gaussian Graphical Models \n")
  cat("--- \n")
  cat("Type: GGM Compare with Bayesian Hypothesis Testing \n")
  p <- x$object$info$dat_info$p[1]
  cat("Posterior Samples:", x$object$iter, "\n")
  cat("Observations: \n")
  groups <- length(x$object$info$dat)
  for (i in 1:groups) {
    cat("  Group",
        paste(i, ":", sep = "") ,
        x$object$info$dat_info$n[[i]],
        "\n")
  }
  cat("Variables (p):", p, "\n")
  cat("Edges:", .5 * (p * (p - 1)), "\n")
  cat("Groups:", nrow(x$object$info$dat_info), "\n")
  cat("Prior SD:", x$object$prior_sd, "\n")
  cat("--- \n")
  if (is.null(x$hypotheses)) {
    cat("Call: \n")
    print(x$object$call)
    cat("--- \n")
  }
  cat("Hypotheses:\n")
  cat("H0: rho = 0\n")
  cat("H1: rho != 0\n")
  cat("--- \n")
  cat("Estimates:\n")
  print(x$results, row.names = F, ...)
  cat("--- \n")
}
