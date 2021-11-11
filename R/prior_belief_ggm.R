#' Prior Belief Gaussian Graphical Model
#'
#' @description Incorporate prior information into the estimation of the
#' conditional dependence structure. This prior information is expressed as
#' the prior odds that each relation should be included in the graph.
#'
#' @param Y Matrix (or data frame) of dimensions \emph{n} (observations) by
#'          \emph{p} (variables/nodes).
#'
#' @param prior_ggm Matrix of dimensions \emph{p} by \emph{p}, encoding the prior
#'                  odds for including each relation in the graph (see '\code{Details}')
#'
#' @param post_odds_cut Numeric. Threshold for including an edge (defaults to 3).
#'                      Note \code{post_odds} refers to posterior odds.
#'
#' @param ... Additional arguments passed to \code{\link{explore}}.
#'
#'
#' @details Technically, the prior odds is not for including an edge in the graph,
#' but for (H1)/p(H0), where H1 captures the hypothesized edge size and H0 is the
#' null model  \insertCite{@see Williams2019_bf}{BGGM}. Accordingly, setting an
#' entry in \code{prior_ggm} to, say, 10, encodes a prior belief that H1 is 10 times
#' more likely than H0. Further, setting an entry in \code{prior_ggm} to 1 results
#' in equal prior odds (the default in \code{\link{select.explore}}).
#'
#'
#' @return An object including:
#'
#' \itemize{
#'
#' \item{\strong{adj}}: Adjacency matrix
#'
#' \item{\strong{post_prob}}: Posterior probability for the
#'                            alternative hypothesis.
#'
#' }
#'
#' @examples
#' \donttest{
#' # Assume perfect prior information
#' # synthetic ggm
#' p <- 20
#' main <- gen_net()
#'
#' # prior odds 10:1, assuming graph is known
#' prior_ggm <- ifelse(main$adj == 1, 10, 1)
#'
#' # generate data
#' y <- MASS::mvrnorm(n = 200,
#'                    mu = rep(0, 20),
#'                    Sigma = main$cors)
#'
#' # prior est
#' prior_est <- prior_belief_ggm(Y = y,
#'                               prior_ggm = prior_ggm,
#'                               progress = FALSE)
#'
#' # check scores
#' BGGM:::performance(Estimate = prior_est$adj,
#'                    True = main$adj)
#'
#' # default in BGGM
#' default_est <- select(explore(y, progress = FALSE))
#'
#' # check scores
#' BGGM:::performance(Estimate = default_est$Adj_10,
#'                    True = main$adj)
#'
#' }
#'
#' @export
prior_belief_ggm <- function(Y,
                             prior_ggm,
                             post_odds_cut = 3,
                             ...){

  if (any(dim(prior_ggm) != ncol(Y))) {
    stop("prior_ggm must be a square matrix")
  }

  check_symmetric <-
    all.equal(prior_ggm[lower.tri(prior_ggm)],
          t(prior_ggm)[lower.tri(prior_ggm)])

  if (isFALSE(check_symmetric)){
    stop("prior_ggm must be symmetric")
  }

  if(any(prior_ggm == 0)){
    stop("zeros are not allowed in prior_ggm")
  }

  fit <- explore(Y, ...)

  sel <- select(fit)

  post_odds <- sel$BF_10 * prior_ggm

  adj <- ifelse(post_odds > post_odds_cut, 1, 0)

  post_prob <- post_odds / (1 + post_odds)

  returned_object <- list(adj = adj,
                          post_prob = post_prob)

  return(returned_object)
}
