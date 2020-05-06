#' Confirmatory Hypothesis Testing
#'
#' @description Confirmatory hypothesis testing in GGMs, for, say, comparing formalized models encoding
#' scienctific expectations. Hypotheses are expressed as equality and/or ineqaulity contraints on the
#' partial correlations of interest. Here the focus is \emph{not} on determining the graph (see \code{\link{explore}}) but testing specific hypotheses related to
#' the conditional (in)dependence structure. These methods were introduced in \insertCite{Williams2019_bf;textual}{BGGM}.
#'
#' @param Y  matrix (or data frame) of dimensions \emph{n} (observations) by  \emph{p} (variables).
#'
#' @param hypothesis character string. hypothesis (or hypotheses) to be tested. See notes for futher details.
#'
#' @param prior_sd hypothesized standard deviation of the prior distribution.
#'
#' @param formula an object of class \code{\link[stats]{formula}}. This allows for including
#' control variables in the model (e.g.,, \code{~ gender * education}).
#'
#' @param type character string. Which type of data for \strong{Y} ? The options include \code{continuous},
#' \code{binary}, \code{ordinal}, or \code{mixed}. See the note for further details.
#'
#' @param mixed_type numeric vector of length \emph{p}. An indicator for which varibles should be treated as ranks.
#' (1 for rank and 0 to assume normality). The default is currently (dev version) to treat all integer variables
#' as ranks when \code{type = "mixed"} and \code{NULL} otherwise. See note for further details.
#'
#' @param iter posterior and prior samples. 25,000 is the default, as it results in a more stable Bayes factor than
#' using, say, 5,000.
#'
#'
#' @references
#' \insertAllCited{}
#'
#' @return list of class \code{confirm}:
#'
#' \itemize{
#' \item \code{BF_matrix} matrix of Bayes factors for each hypothesis. Also includes the compliment
#' \item \code{post_prob} posterior hypothesis probabilities
#' \item \code{hypotheses} \code{hypothesis}
#' \item \code{call} match.call()
#' \item \code{p} number of variables
#' \item \code{n} number of observations
#' \item \code{iter} number of posterior samples
#' \item \code{delta} hyperparameter of matrix-F prior distribution (corresponds to prior_sd)
#' \item \code{parcors_mat} partial correlation matrix
#' \item \code{returned_mats} contrast matrices
#' }
#'
#' @importFrom MASS ginv
#' @importFrom BFpack BF
#' @importFrom stats cov rbeta

#'
#' @note
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
#' \strong{Hypothesis Examples:}
#'
#' The hypotheses can be written either with the respective column names or numbers. For example, \code{1--2} denotes the
#' relation between variables in column 1 and 2. Note that these must correspond to the upper triangular elements of the correlation
#' matrix. This is accomplished by ensuring that the first number is smaller than the second number. This also applies when using
#' column names (e.g., in reference to the column number).
#'
#'One hypothesis:
#' \itemize{
#'
#' \item  To test whether some relations are larger than others, while others are expected to be equal, this can be writting as
#'
#'  \code{hyp <-  c(1--2 > 1--3  = 1--4 > 0)},
#'
#'  where there is an addition additional contraint that all effects are expected to be positive. This is then compared to the complement.
#'}
#'
#'More than one hypothesis:
#'\itemize{
#'
#'\item The above hypothesis can also be compared to, say, a null model by using ";" to seperate the hypotheses, for example,
#'
#' \code{hyp <-  c(1--2 > 1--3  = 1--4 > 0; 1--2 = 1--3  = 1--4 = 0)}.
#'
#' Any number of hypotheses can be tested.
#'}
#'
#'
#'
#' see \code{methods(class = "confirm")}
#'
#' @examples
#' \donttest{
#'
#' # p = 10
#' Y <- BGGM::bfi[,1:10]
#'
#' # hypothesis
#' hypothesis <- c("1--2 > 1--3 > 1--4 > 1--5")
#'
#' # test inequality contraint
#' test_order <-  confirm(Y = Y, hypothesis  = hypothesis,
#'                       prior_sd = 0.5, iter = 50000,
#'                       cores = 2)
#' # summary
#' summary(test_order)
#'
#'
#'# test hypothesized directions
#'
#'# hypothesis
#'hypothesis <- c("(1--2, 1--3, 1--4)  <  0 < (1--6)")
#'
#'# test directions
#' test_directions <-  confirm(Y = Y, hypothesis  = hypothesis,
#'                       prior_sd = 0.5, iter = 50000,
#'                       cores = 2)
#'# summary
#' test_directions
#'
#'#############################
#'######## Binary #############
#'#############################
#'
#'
#'
#'
#'}
#'@export
confirm <- function(Y, hypothesis,
                    prior_sd = 0.25,
                    formula = NULL,
                    type = "continuous",
                    mixed_type = NULL,
                    iter = 25000){

  priorprob <- 1

  # hyperparameter
  delta <- delta_solve(prior_sd)

  p <- ncol(Y)

  I_p <- diag(p)

  # number of edges
  pcors <- p*(p-1)*0.5

  message("BGGM: Posterior Sampling")

  # continuous
  if(type == "continuous"){




    # no control
    if(is.null(formula)){

      # na omit
      Y <- as.matrix(na.omit(Y))

      # scale Y
      Y <- scale(Y, scale = F)

      # nodes
      p <- ncol(Y)

      n <- nrow(Y)

      # posterior sample
      post_samp <- .Call(
        '_BGGM_Theta_continuous',
        PACKAGE = 'BGGM',
        Y = Y,
        iter = iter + 50,
        delta = delta,
        epsilon = 0.1,
        prior_only = 0,
        explore = 1
      )

      # control for variables
    } else {

      control_info <- remove_predictors_helper(list(as.data.frame(Y)),
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
      post_samp <- .Call(
        "_BGGM_mv_continuous",
        Y = Y,
        X = X,
        delta = delta,
        epsilon = 0.1,
        iter = iter + 50
      )

    } # end control

    # binary
    } else if (type == "binary"){


      # intercept only
      if (is.null(formula)) {

        # data
        Y <- as.matrix(na.omit(Y))

        # obervations
        n <- nrow(Y)

        # nodes
        p <- ncol(Y)

        X <- matrix(1, n, 1)

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

      }

  # posterior sample
    post_samp <-  .Call(
      "_BGGM_mv_binary",
      Y = Y,
      X = X,
      delta = delta,
      epsilon = 0.1,
      iter = iter + 50,
      beta_prior = 0.0001,
      cutpoints = c(-Inf, 0, Inf)
    )

    # ordinal
  } else if(type == "ordinal"){

    # intercept only
    if(is.null(formula)){

      # data
      Y <- as.matrix(na.omit(Y))

      # obervations
      n <- nrow(Y)

      # nodes
      p <- ncol(Y)

      # intercept only
      X <- matrix(1, n, 1)

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

      }

    # categories
    K <- max(apply(Y, 2, function(x) { length(unique(x))   } ))

    # call c ++
    post_samp <- .Call(
      "_BGGM_mv_ordinal_albert",
      Y = Y,
      X = X,
      iter = iter + 50,
      delta = delta,
      epsilon = 0.1,
      K = K
    )

  } else if(type == "mixed"){

    # no control variables allowed
    if (!is.null(formula)) {

      warning("formula ignored for mixed data at this time")

      control_info <- remove_predictors_helper(list(as.data.frame(Y)),
                                               formula = formula)

      # data
      Y <-  as.matrix(control_info$Y_groups[[1]])

      formula <- NULL

    } else {

      Y <- na.omit(Y)
    }

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

  } else {

    stop("'type' not supported: must be continuous, binary, ordinal, or mixed.")
  }

  message("BGGM: Finished")

 # sample prior
  prior_samp <- .Call(
    '_BGGM_sample_prior',
    PACKAGE = 'BGGM',
    Y = Y,
    iter = 10000,
    delta = delta,
    epsilon = 0.01,
    prior_only = 1,
    explore = 0
  )$fisher_z

  col_names <- BGGM:::numbers2words(1:p)

  mat_name <- sapply(col_names, function(x) paste(col_names,x, sep = ""))[upper.tri(I_p)]

  posterior_samples <- matrix(post_samp$fisher_z[, , 51:(iter+50)][upper.tri(I_p)],
                              iter, pcors,
                              byrow = TRUE)

  prior_samples <- matrix(prior_samp[,,][upper.tri(I_p)], 10000, pcors, byrow = TRUE)


  colnames(posterior_samples) <- mat_name

  colnames(prior_samples) <- mat_name

  prior_mu <- colMeans(prior_samples)

  prior_cov <- cov(prior_samples)

  post_mu <- colMeans(posterior_samples)

  post_cov <- cov(posterior_samples)

  BFprior <- BF(prior_mu,
                Sigma = prior_cov,
                hypothesis = BGGM:::convert_hyps(hypothesis, Y = Y),
                n = 1)


  BFpost <- BF(post_mu,
               Sigma = post_cov,
               hypothesis = BGGM:::convert_hyps(hypothesis, Y = Y),
               n = 1)


  # number of hypotheses
  n_hyps <- nrow(BFpost$BFtable_confirmatory)

  # BF against unconstrained
  BF_tu <- NA

  for(i in seq_len(n_hyps)){

    # BF tu
    BF_tu[i] <- prod(BFpost$BFtable_confirmatory[i,3:4] / BFprior$BFtable_confirmatory[i,3:4])

    }

  # posterior hyp probs
  out_hyp_prob <- BF_tu*priorprob / sum(BF_tu*priorprob)

  # BF matrix
  BF_matrix <- matrix(rep(BF_tu, length(BF_tu)),
                      ncol = length(BF_tu),
                      byrow = TRUE)

  BF_matrix[is.nan(BF_matrix)] <- 0
  diag(BF_matrix) <- 1

  BF_matrix <- t(BF_matrix) / (BF_matrix)

  row.names(BF_matrix) <- row.names( BFpost$BFtable_confirmatory)
  colnames(BF_matrix) <- row.names( BFpost$BFtable_confirmatory)

  returned_object <- list(BF_matrix = BF_matrix,
                          out_hyp_prob = out_hyp_prob,
                          info = BFpost,
                          dat = Y,
                          type = type,
                          call = match.call(),
                          hypothesis = hypothesis,
                          iter = iter, p = p,
                          delta = delta)

  class(returned_object) <- c("BGGM", "confirm")
  returned_object
}
