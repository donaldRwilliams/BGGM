#' GGM Compare: Confirmatory Hypothesis Testing
#'
#' @name ggm_compare_confirm
#'
#' @param ... at least two matrices (or data frame) of dimensions \emph{n} (observations) by  \emph{p} (variables).
#'
#' @param hypothesis character string. hypothesis (or hypotheses) to be tested. See notes for futher details.
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
#' @param iter number of iterations (posterior samples; defaults to 5000).
#'
#' @param seed The seed for random number generation (default set to \code{1}).
#'
#' @return
#' @export
#'
#' @examples
ggm_compare_confirm <- function(...,
                                hypothesis,
                                formula = NULL,
                                data = NULL,
                                type = "continuous",
                                mixed_type = NULL,
                                prior_sd = 0.20,
                                iter = 5000){

  # set seed
  # set.seed(seed)

  # prior prob
  priorprob <- 1

  # delta parameter
  delta <- BGGM:::delta_solve(prior_sd)

  # group info (e.g., n, p, etc)
  info <- BGGM:::Y_combine(...)

  # number of variables
  p <- info$dat_info$p[1]

  # groups
  groups <- length(info$dat)

  # number of edges
  pcors <- 0.5 * (p * (p - 1))

  # identity matrix
  I_p <- diag(p)

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

  prior_samp <- lapply(1:groups, function(x) {

    Y <- as.matrix(info$dat[[x]])

    .Call(
      '_BGGM_Theta_continuous',
      PACKAGE = 'BGGM',
      Y = Y,
      iter = 10000,
      delta = delta,
      epsilon = 0.1,
      prior_only = 1,
      explore = 0
    )$fisher_z
  })

  # colnames: post samples
  col_names <- BGGM:::numbers2words(1:p)
  mat_names <- lapply(1:groups, function(x) paste0("g", BGGM:::numbers2words(x),
               sapply(col_names, function(x) paste(col_names, x, sep = ""))[upper.tri(I_p)]))

  # posterior start group (one)
  post_group <- matrix(post_samp[[1]]$fisher_z[, , 51:(iter + 50)][upper.tri(I_p)],
                       iter, pcors, byrow = TRUE)

  # prior start group (one)
  prior_group <-  matrix(prior_samp[[1]][ , ,][upper.tri(I_p)],
                         nrow = iter,
                         ncol = pcors,
                         byrow = TRUE)

  # post group
  for(j in 2:(groups)){

    post_group <-  cbind(post_group,
                       matrix(post_samp[[j]]$fisher_z[, , 51:(iter+50)][upper.tri(I_p)],
                              nrow = iter, ncol = pcors,
                              byrow = TRUE))

    prior_group <-  cbind(prior_group,
                        matrix(prior_samp[[j]][ , ,][upper.tri(I_p)], iter, pcors, byrow = TRUE))
  }

  posterior_samples <- post_group
  colnames(posterior_samples) <- unlist(mat_names)

  prior_samples <- prior_group
  colnames(prior_samples) <- unlist(mat_names)

  prior_mu <- colMeans(prior_samples)

  prior_cov <- cov(prior_samples)

  post_mu <- colMeans(posterior_samples)

  post_cov <- cov(posterior_samples)

  BFprior <- BF(prior_mu,
              Sigma = prior_cov,
              hypothesis = group_hyp_helper(hypothesis, x = info$dat[[1]]),
              n = 1)

  BFpost <- BF(post_mu,
             Sigma = post_cov,
             hypothesis = group_hyp_helper(hypothesis, x = info$dat[[1]]),
             n = 1)

  # number of hypotheses
  n_hyps <- nrow(BFpost$BFtable_confirmatory)

  # BF against unconstrained
  BF_tu <- NA

  for (i in seq_len(n_hyps)) {
    # BF tu
    BF_tu[i] <-
      prod(BFpost$BFtable_confirmatory[i, 3:4] / BFprior$BFtable_confirmatory[i, 3:4])

  }

  # posterior hyp probs
  out_hyp_prob <- BF_tu * priorprob / sum(BF_tu * priorprob)

  # BF matrix
  BF_matrix <- matrix(rep(BF_tu, length(BF_tu)),
                      ncol = length(BF_tu),
                      byrow = TRUE)

  BF_matrix[is.nan(BF_matrix)] <- 0
  diag(BF_matrix) <- 1

  BF_matrix <- t(BF_matrix) / (BF_matrix)

  row.names(BF_matrix) <- row.names(BFpost$BFtable_confirmatory)
  colnames(BF_matrix) <- row.names(BFpost$BFtable_confirmatory)

  returned_object <- list(
    BF_matrix = BF_matrix,
    out_hyp_prob = out_hyp_prob,
    info = BFpost,
    groups = groups,
    info_dat = info,
    type = type,
    call = match.call(),
    hypothesis = hypothesis,
    iter = iter,
    p = p,
    posterior_samples = posterior_samples,
    post_group = post_group,
    delta = delta
  )

class(returned_object) <- c("BGGM", "confirm", "ggm_compare_confirm")
returned_object

}
