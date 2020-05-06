#' GGM Compare: Confirmatory Hypothesis Testing
#'
#' @name ggm_compare_confirm
#'
#' @param ... at least two matrices (or data frame) of dimensions \emph{n} (observations) by  \emph{p} (nodes).
#'
#' @param hypothesis character string. hypothesis (or hypotheses) to be tested. See notes for futher details.
#'
#' @param formula an object of class \code{\link[stats]{formula}}. This allows for including
#' control variables in the model (i.e., \code{~ gender}).
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
                                type = "continuous",
                                mixed_type = NULL,
                                prior_sd = 0.20,
                                iter = 5000){
  # prior prob
  priorprob <- 1

  # delta parameter
  delta <- delta_solve(prior_sd)

  # combine data
  dat_list <- list(...)

  # combine data
  info <- Y_combine(...)

  # groups
  groups <- length(info$dat)

  if(type == "continuous"){

    if(is.null(formula)){

      post_samp <- lapply(1:groups, function(x) {

        message("BGGM: Posterior Sampling ", "(Group ",x ,")")

        # data
        Y <- as.matrix(scale(dat_list[[x]], scale = F))

        Y <- na.omit(Y)

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

      # formula
    } else {

      post_samp <- lapply(1:groups, function(x) {

        message("BGGM: Posterior Sampling ", "(Group ",x ,")")

        control_info <- remove_predictors_helper(list(as.data.frame(dat_list[[x]])),
                                                 formula = formula)

        # data
        Y <- as.matrix(scale(control_info$Y_groups[[1]], scale = F))

        # nodes
        p <- ncol(Y)

        # observations
        n <- nrow(Y)

        # model matrix
        X <- as.matrix(control_info$model_matrices[[1]])

        # posterior sample
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

    } else if(type == "binary"){

    # intercept only
    if (is.null(formula)) {


      post_samp <- lapply(1:groups, function(x) {
        message("BGGM: Posterior Sampling ", "(Group ",x ,")")

        # data
        Y <- as.matrix(na.omit(dat_list[[x]]))

        # obervations
        n <- nrow(Y)

        # nodes
        p <- ncol(Y)

        X <- matrix(1, n, 1)

        # posterior sample
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

    } else {

      post_samp <- lapply(1:groups, function(x) {

        message("BGGM: Posterior Sampling ", "(Group ",x ,")")
        control_info <- remove_predictors_helper(list(as.data.frame(dat_list[[x]])),
                                                 formula = formula)

        # data
        Y <-  as.matrix(control_info$Y_groups[[1]])

        # observations
        n <- nrow(Y)

        # nodes
        p <- ncol(Y)

        # model matrix
        X <- as.matrix(control_info$model_matrices[[1]])

        # posterior sample
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
    }

  } else if(type == "ordinal"){

    if(is.null(formula)){

      post_samp <- lapply(1:groups, function(x) {

        message("BGGM: Posterior Sampling ", "(Group ",x ,")")
        # data
        Y <- as.matrix(na.omit(dat_list[[x]]))

        # obervations
        n <- nrow(Y)

        # nodes
        p <- ncol(Y)

        X <- matrix(1, n, 1)

        # categories
        K <- max(apply(Y, 2, function(x) { length(unique(x)) } ))

        # posterior sample
        # call c ++
         .Call(
          "_BGGM_mv_ordinal_albert",
          Y = Y,
          X = X,
          iter = iter + 50,
          delta = delta,
          epsilon = 0.1,
          K = K)
      })

      } else {

        post_samp <- lapply(1:groups, function(x) {

          message("BGGM: Posterior Sampling ", "(Group ", x ,")")

          control_info <- remove_predictors_helper(list(as.data.frame(dat_list[[x]])),
                                                   formula = formula)

          # data
          Y <-  as.matrix(control_info$Y_groups[[1]])

          # observations
          n <- nrow(Y)

          # nodes
          p <- ncol(Y)

          # model matrix
          X <- as.matrix(control_info$model_matrices[[1]])

          # categories
          K <- max(apply(Y, 2, function(x) { length(unique(x)) } ))

          # posterior sample
          # call c ++
          .Call(
            "_BGGM_mv_ordinal_albert",
            Y = Y,
            X = X,
            iter = iter + 50,
            delta = delta,
            epsilon = 0.1,
            K = K)
        })
  }

} else if(type == "mixed") {

    if(!is.null(formula)){

      warning("formula ignored for mixed data at this time")

      post_samp <- lapply(1:groups, function(x) {

        message("BGGM: Posterior Sampling ", "(Group ",x ,")")

        control_info <- remove_predictors_helper(list(as.data.frame(dat_list[[x]])),
                                                 formula = formula)

        # data
        Y <-  as.matrix(control_info$Y_groups[[1]])

        Y <- na.omit(Y)

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
        })

    } else {


      post_samp <- lapply(1:groups, function(x) {

        message("BGGM: Posterior Sampling ", "(Group ",x ,")")

        Y <- na.omit(dat_list[[x]])

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
      })
      }

  } else {
    stop("'type' not supported: must be continuous, binary, ordinal, or mixed.")
}

  message("BGGM: Finished")

  # sample prior
  if(is.null(formula)){

    Yprior <- as.matrix(dat_list[[1]])

    prior_samp <- lapply(1:groups, function(x) {
      .Call(
        '_BGGM_sample_prior',
        PACKAGE = 'BGGM',
        Y = Yprior,
        iter = 10000,
        delta = delta,
        epsilon = 0.1,
        prior_only = 1,
        explore = 0
      )$fisher_z
    })

  } else {

    control_info <- remove_predictors_helper(list(as.data.frame(dat_list[[1]])),
                                             formula = formula)

    Yprior <- as.matrix(scale(control_info$Y_groups[[1]], scale = F))

    prior_samp <- lapply(1:groups, function(x) {

      .Call(
        '_BGGM_sample_prior',
        PACKAGE = 'BGGM',
        Y = Yprior,
        iter = 10000,
        delta = delta,
        epsilon = 0.1,
        prior_only = 1,
        explore = 0
      )$fisher_z
    })
  }

  # nodes
  p <- ncol(Yprior)

  # number of pcors
  pcors <- 0.5 * (p * (p - 1))

  # identity matrix
  I_p <- diag(p)

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
  out_hyp_prob <- (BF_tu * priorprob) / sum(BF_tu * priorprob)

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
    delta = delta,
    formula = formula,
    dat_list = dat_list,
    post_samp = post_samp

  )

class(returned_object) <- c("BGGM",
                            "confirm",
                            "ggm_compare_confirm")
returned_object

}



print_ggm_confirm <- function(x, ...){
  groups <- x$groups
  info <- x$info_dat
  cat("BGGM: Bayesian Gaussian Graphical Models \n")

  cat("Type:",  x$type ,  "\n")

  cat("--- \n")

  cat("Posterior Samples:", x$iter, "\n")

  for(i in 1:groups){
    cat("  Group", paste( i, ":", sep = "") , info$dat_info$n[[i]], "\n")
  }
  # number of variables
  cat("Variables (p):", x$p, "\n")
  # number of edges
  cat("Relations:", .5 * (x$p * (x$p-1)), "\n")
  cat("Delta:", x$delta, "\n")
  cat("--- \n")

  cat("Call:\n")

  print(x$call)

  cat("--- \n")

  cat("Hypotheses: \n\n")

  hyps <- strsplit(x$hypothesis, ";")

  n_hyps <- length(hyps[[1]])

  x$info$hypotheses[1:n_hyps] <- hyps[[1]]

  n_hyps <- length(x$info$hypotheses)

  for (h in seq_len(n_hyps)) {
    cat(paste0("H", h,  ": ", gsub(" ", "", x$info$hypotheses[h])  ))
    cat("\n")
  }

  cat("--- \n")

  cat("Posterior prob: \n\n")

  for(h in seq_len(n_hyps)){

    cat(paste0("p(H",h,"|data) = ", round(x$out_hyp_prob[h], 3 )  ))
    cat("\n")
  }

  cat("--- \n")

  cat('Bayes factor matrix: \n')

  print(round(x$BF_matrix, 3))

  cat("--- \n")

  cat("note: equal hypothesis prior probabilities")
}


