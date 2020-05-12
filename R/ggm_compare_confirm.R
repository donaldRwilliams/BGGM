#' GGM Compare: Confirmatory Hypothesis Testing
#'
#' #' @description Confirmatory hypothesis testing for comparing GGMs. Hypotheses are expressed as equality
#' and/or ineqaulity contraints on the partial correlations of interest. Here the focus is \emph{not}
#' on determining the graph (see \code{\link{explore}}) but testing specific hypotheses related to
#' the conditional (in)dependence structure. These methods were introduced in
#' \insertCite{Williams2019_bf;textual}{BGGM} and in \insertCite{williams2020comparing;textual}{BGGM}
#'
#' @name ggm_compare_confirm
#'
#' @param ... At least two matrices (or data frame) of dimensions \emph{n} (observations) by  \emph{p} (nodes).
#'
#' @param hypothesis Character string. The hypothesis (or hypotheses) to be tested. See notes for futher details.
#'
#' @param formula an object of class \code{\link[stats]{formula}}. This allows for including
#' control variables in the model (i.e., \code{~ gender}).
#'
#' @param prior_sd Numeric. The scale of the prior distribution (centered at zero),
#'                 in reference to a beta distribtuion (defaults to 0.25).
#'
#' @param type Character string. Which type of data for \strong{Y} ? The options include \code{continuous},
#' \code{binary}, or \code{ordinal}. See the note for further details.
#'
#' @param mixed_type numeric vector. An indicator of length p for which varibles should be treated as ranks.
#' (1 for rank and 0 to assume normality). The default is currently (dev version) to treat all integer variables
#' as ranks when \code{type = "mixed"} and \code{NULL} otherwise. See note for further details.
#'
#' @param iter  Number of iterations (posterior samples; defaults to 25,000).
#'
#' @param progress Logical. Should a progress bar be included (defaults to \code{TRUE}) ?
#'
#' @param seed An integer for the random seed.
#'
#' @return The returned object of class \code{confirm} contains a lot of information that
#'         is used for printing and plotting the results. For users of \strong{BGGM}, the following
#'         are the useful objects:
#'
#' \itemize{
#'
#' \item \code{out_hyp_prob} Posterior hypothesis probabilities.
#'
#' \item \code{info} An object of class \code{BF} from the R package \strong{BFpack}.
#'
#' }
#'
#' @details
#' The hypotheses can be written either with the respective column names or numbers.
#' For example, \code{g1_1--2} denotes the relation between the variables in column 1 and 2 for group 1.
#' The \code{g1_1} is required and the only difference from \code{\link{confirm}} (one group).
#' Note that these must correspond to the upper triangular elements of the correlation
#' matrix. This is accomplished by ensuring that the first number is smaller than the second number.
#' This also applies when using column names (i.e,, in reference to the column number).
#'
#' \strong{One Hypothesis}:
#'
#' \itemize{
#'
#' \item  To test whether a relation in larger in one group, while both are expected
#'        to be positive,  this can be written as
#'
#'  \code{hyp <-  c(g1_1--2 > g2_1--2 > 0)},
#'
#'  This is then compared to the complement.
#'}
#'
#' \strong{More Than One Hypothesis}:
#'
#' \itemize{
#'
#' \item The above hypothesis can also be compared to, say, a null model by using ";"
#'       to seperate the hypotheses, for example,
#'
#' \code{hyp <-  c(g1_1--2 > g2_1--2 > 0; g1_1--2 = g2_1--2 = 0)}.
#'
#' Any number of hypotheses can be compared this way.
#'}
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
#'  continuous or \emph{only} discrete variables \insertCite{hoff2007extending}{BGGM}. This is based on the
#'  ranked likelihood which requires sampling the ranks for each variable (i.e., the data is not merely
#'  transformed to ranks). This is computationally expensive when there are many levels. For example,
#'  with continuous data, there are as many ranks as data points!
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
#' @note
#'
#' \strong{"Default" Prior}:
#'
#'  In Bayesian statistics, a default Bayes factor needs to have several properties. I refer
#'  interested users to \insertCite{@section 2.2 in @dablander2020default;textual}{BGGM}. In
#'  \insertCite{Williams2019_bf;textual}{BGGM}, some of these propteries were investigated (e.g.,
#'  model selection consistency). That said, we would not consider this a "default" or "automatic"
#'  Bayes factor and thus we encourage users to perform sensitivity analyses by varying the scale of the prior
#'  distribution.
#'
#'  Furthermore, it is important to note there is no "correct" prior and, also, there is no need
#'  to entertain the possibility of a "true" model. Rather, the Bayes factor can be interpreted as
#'  which hypothesis best (relative to each other) predicts the observed data
#'  \insertCite{@Section 3.2 in @Kass1995}{BGGM}.
#'
#' \strong{Interpretation of Conditional (In)dependence Models for Latent Data}:
#'
#'  See \code{\link{BGGM-package}} for details about interpreting GGMs based on latent data
#' (i.e, all data types besides \code{"continuous"})
#'
#' @examples
#' @export
ggm_compare_confirm <- function(...,
                                hypothesis,
                                formula = NULL,
                                type = "continuous",
                                mixed_type = NULL,
                                prior_sd = 0.20,
                                iter = 25000,
                                progress = TRUE,
                                seed = 1){



  old <- .Random.seed

  set.seed(seed)

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

      set.seed(x)

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

  .Random.seed <<- old

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


