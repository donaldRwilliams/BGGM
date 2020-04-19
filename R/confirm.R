#' Confirmatory Hypothesis Testing in GGMs
#'
#' @description Traditionally, Gaussian graphical models (GGM) are inherently exploratory. That is, automated model selection is performed. A key aspect of \strong{BGGM} is the ability to extend inference beyond exploratory and to
#' confirmatory hypothesis testing. This is accomplished by testing equality and/or inequality constraints for sets of edges (partial correlations).
#'
#' @param Y  matrix (or data frame) of dimensions \emph{n} (observations) by  \emph{p} (variables).
#'
#' @param hypothesis character string. hypothesis (or hypotheses) to be tested. See notes for futher details.
#'
#' @param prior_sd hypothesized standard deviation of the prior distribution
#'
#' @param formula an object of class \code{\link[stats]{formula}}. This allows for including
#' control variables in the model (e.g.,, \code{~ gender * education}).
#'
#' @param data an optional data frame, list or environment (or an object coercible by \code{\link[base]{as.data.frame}})
#' to a data frame containing the variables in \code{formula}. This is required when controlling for variables.
#'
#' @param type character string. Which type of data for \strong{Y} ? The options include \code{continuous},
#' \code{binary}, or \code{ordinal}. See the note for further details.
#'
#' @param iter posterior and prior samples. 25,000 is the default, as it results in a more stable Bayes factor than
#' using, say, 5,000.
#'
#' @param cores integer. How many cores for parallel computing ? (default is 2).
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
#' @export
#'
#' @note Currently inequality and equality restrictions can be tested. The former is an ordering the respective edge sizes,
#' whereas the latter allows for testing whether certain edges are exactly the same.
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
confirm <- function(Y, hypothesis,
                    prior_sd = 0.25,
                    formula = NULL,
                    data = NULL,
                    type = "continuous",
                    iter = 25000,
                    cores = 2){

  # set prior prob to 1--i.e., equal
  priorprob = 1

  # hyperparameter
  delta <- delta_solve(prior_sd)

  # variables
  p <- ncol(Y)

  # number of edges
  edges <- p*(p-1)*0.5

  if(type == "continuous"){
  fit <- explore(Y,  prior_sd = prior_sd,
                 iter = iter,  cores = cores)


  # posterior
  posterior_samples <- do.call(rbind.data.frame,
                               lapply(1:fit$cores, function(z)
                                 fit$samples[[z]]$fisher_z_post))[1:edges]

  # prior samples
  prior_samples <-  do.call(rbind.data.frame,
                               lapply(1:fit$cores, function(z)
                                 fit$samples[[z]]$fisher_z_prior))[1:edges]

  } else if( type == "binary"){

    samples <- sampling_helper_poly(Y, delta, iter)

    posterior_samples <- samples$fisher_z_post[,1:edges]

    prior_samples <- samples$fisher_z_prior[,1:edges]

  }

  prior_mu <- colMeans(prior_samples)

  prior_cov <- cov(prior_samples)

  post_mu <- colMeans(posterior_samples)

  post_cov <- cov(posterior_samples)

  BFprior <- BF(prior_mu,
                Sigma = prior_cov,
                hypothesis = convert_hyps(hypothesis, Y = Y),
                n = 1)


  BFpost <- BF(post_mu,
               Sigma = post_cov,
               hypothesis = convert_hyps(hypothesis, Y = Y),
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
  BF_matrix

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
