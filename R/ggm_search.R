#' @title Search Algorithm
#'
#' @param x Either a \code{n} by \code{p} data matrix or a \code{p} by \code{p} correlation matrix.
#'
#' @param n Numeric. Sample size. This is required if \code{x} is a correlation matrix.
#'
#' @param method Character string. Which search algorithm should be used? The default is
#'               MCMC model composition \insertCite{madigan1994model}{BGGM}.
#'
#' @param prior_prob Numeric. Prior inclusion probability for an individual edge. The default is
#'        \code{0.1} which favors graphs with less edges. If \code{prior_prob = 0.5}, this results
#'        in equal prior model probabilities.
#'
#' @param iter Numeric. The number of iterations (defaults to 5,000). Note that this is not the
#'             number of posterior samples, but the number of proposed graphs.
#'
#' @param stop_early Numeric. Should the sampler be stopped early when the graph does not change? The
#'                   default is \code{1000} (i.e., if 1000 graphs have been
#'                   rejected in a row, the search stops early).
#'
#' @param bma_mean Logical. Should the Bayesian model averaged mean be computed (defaults to \code{TRUE})?
#'
#' @param progress Logical. Should a progress bar be included (defaults to \code{TRUE}) ?
#'
#' @param seed An integer for the random seed (defaults to 1).
#'
#' @param ... Currently ignored.
#'
#' @references
#' \insertAllCited{}
#'
#'
#' @details The algorithm does not scale particularly well with \code{p} (the number of nodes).
#'
#' @return A list of class \code{ggm_search}. For users of \strong{BGGM}, the following
#'         are the useful objects:
#'
#' \itemize{
#' \item \code{pcor_adj} Selected partial correlation matrix (weighted adjacency). This corresponds
#'                       to the most probable model.
#'
#' \item \code{Theta_map} Selected precision matrix. This corresponds
#'                       to the most probable model.
#'
#' \item \code{Theta_bma} Bayesian model averaged precision matrix
#'
#' \item \code{pcor_bma} Bayesian model averaged partial correlation matrix.
#'
#'
#' }
#'
#' @examples
#' \donttest{
#'
#' # data
#' Y <- ptsd
#'
#' # polychoric partials
#' S <- psych::polychoric(Y)$rho
#'
#' # fit model
#' fit <- ggm_search(x = S,
#'                   n  = nrow(Y),
#'                   stop_early = 1000,
#'                   progress = FALSE)
#'}
#' @export
ggm_search <- function(x, n = NULL,
                       method = "mc3",
                       prior_prob = 0.1,
                       iter = 5000,
                       stop_early = 1000,
                       bma_mean = TRUE,
                       seed = 1,
                       progress = TRUE, ...){

  set.seed(seed)

  if (base::isSymmetric(as.matrix(x))) {
    S <- x
  } else {
    S <- cor(x)
    n <- nrow(x)
  }

  p <- ncol(S)

  if(method == "mc3"){

    if(is.null(stop_early)){

      stop_early <- iter

    }

    pcors <- -cov2cor( solve(S) ) + diag(2, p)

    # test full vs missing one edge
    BF_01 <- exp(-0.5 *  (BGGM:::tstat(r = pcors, n = n, p - 2)^2 - log(n)))
    BF_10 <- 1/BF_01

    adj_start <- ifelse(BF_10 > 1, 1, 0)

    old <- hft_algorithm(Sigma = S,
                         adj = adj_start,
                         tol = 1e-10,
                         max_iter = 100)

    bic_old <- bic_fast(Theta = old$Theta,
                        S = S,
                        n = n,
                        prior_prob = prior_prob)

    if(isTRUE(progress)){
      message(paste0("BGGM: Sampling Graphs", ...))
    }

    fit <- search(S = S,
                  iter = iter,
                  old_bic = bic_old,
                  start_adj = adj_start,
                  gamma = prior_prob,
                  n = n,
                  stop_early = stop_early,
                  progress = progress)

    if(isTRUE(progress)){
      message("BGGM: Finished")
    }

    # accepted
    acc <- fit$acc

    # first matrix (starting values)
    fit$adj[,,1] <- adj_start

    # approximate marginal likelihood
    approx_marg_ll <- fit$bics

    # starting bic
    approx_marg_ll[1] <- bic_old

    if(!is.null(stop_early)){
      approx_marg_ll <- approx_marg_ll[which(approx_marg_ll != 0)]
      fit$adj <- fit$adj[,, which(approx_marg_ll != 0)]
    }

    adj_path <- fit$adj

    selected <- which.min(approx_marg_ll)

    if(acc == 0){
      adj <- fit$adj
    } else {
      adj <-  fit$adj[,,selected]
    }

    # BFs vs mpm
    delta <- approx_marg_ll - min(approx_marg_ll)
    probs <- exp(-0.5 * delta) / sum( exp(-0.5 * delta))

    Theta_map <- hft_algorithm(Sigma = S,
                               adj = adj,
                               tol = 1e-10,
                               max_iter = 100)

    pcor_adj <- -cov2cor(Theta_map$Theta) + diag(2, p)

  }


  if(bma_mean){

    graph_ids <- which(duplicated(approx_marg_ll) == 0)[-1]

    delta <- (approx_marg_ll[graph_ids] - min(approx_marg_ll[graph_ids])) * (6 / (2*sqrt(2*n)))

    probs <- exp(- 0.5 * delta) / sum(exp(- 0.5 * delta))

    graphs <- adj_path[,,graph_ids]

    Theta_bma <- lapply(1:length(probs), function(x)

      hft_algorithm(Sigma = S,
                    adj = graphs[,,x],
                    tol =1e-10,
                    max_iter = 1000)$Theta * probs[x]
    )

    Theta_bma <- Reduce("+", Theta_bma)
    pcor_bma <- -cov2cor(Theta_bma) + diag(2, p)

  } else {

    Theta_bma <- NULL
    pcor_bma <- NULL

  }



  returned_object <- list(pcor_adj = pcor_adj,
                          Theta_map = Theta_map,
                          Theta_bma = Theta_bma,
                          pcor_bma = pcor_bma,
                          adj = adj,
                          adj_start = adj_start,
                          probs = probs,
                          approx_marg_ll = approx_marg_ll,
                          selected = selected,
                          BF_start = BF_10,
                          adj_path = adj_path,
                          acc = acc)

  rm(.Random.seed, envir=.GlobalEnv)

  return( returned_object )
}
