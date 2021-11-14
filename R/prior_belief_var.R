#' Prior Belief Graphical VAR
#'
#' @param Y Matrix (or data frame) of dimensions \emph{n}
#'          (observations) by \emph{p} (variables/nodes).
#'
#' @param prior_temporal Matrix of dimensions \emph{p} by \emph{p},
#'                       encoding the prior odds for including each relation
#'                       in the temporal graph (see '\code{Details}'). If null
#'                       a matrix of 1's is used, resulting in equal prior odds.
#'
#' @param post_odds_cut Numeric. Threshold for including an edge (defaults to 3).
#'                      Note \code{post_odds} refers to posterior odds.
#'
#' @param est_ggm Logical. Should the contemporaneous network be estimated
#'                (defaults to \code{TRUE})?
#'
#' @param prior_ggm Matrix of dimensions \emph{p} by \emph{p}, encoding the prior
#'                  odds for including each relation in the graph
#'                  (see '\code{Details}'). If null a matrix of 1's is used,
#'                  resulting in equal prior odds.
#'
#' @param progress Logical. Should a progress bar be included
#'                 (defaults to \code{TRUE}) ?
#'
#' @param ... Additional arguments passed to \code{\link{explore}}. Ignored
#'            if \code{prior_ggm = FALSE}.
#'
#' @details Technically, the prior odds is not for including an edge in the graph,
#' but for (H1)/p(H0), where H1 captures the hypothesized edge size and H0 is the
#' null model  \insertCite{@see Williams2019_bf}{BGGM}. Accordingly, setting an
#' entry in \code{prior_ggm} to, say, 10, encodes a prior belief that H1 is 10 times
#' more likely than H0. Further, setting an entry in \code{prior_ggm} or
#' \code{prior_var} to 1 results in equal prior odds
#' (the default in \code{\link{select.explore}}).
#'
#' @return An object including (\code{est_ggm = FALSE}):
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
#' An object including (\code{est_ggm = TRUE}):
#' \itemize{
#'
#' \item{\strong{adj_temporal}}: Adjacency matrix for the temporal network.
#'
#' \item{\strong{post_prob_temporal}}: Posterior probability for the
#'                                     alternative hypothesis (temporal edge)
#'
#' \item{\strong{adj_ggm}}: Adjacency matrix for the contemporaneous
#'                          network (ggm).
#'
#' \item{post_prob_ggm}: Posterior probability for the
#'                       alternative hypothesis (contemporaneous edge)
#' }
#'
#'
#' @export
#' @importFrom stats lm residuals
#' @importFrom BFpack BF
#' @examples
#' \donttest{
#' # affect data from 1 person
#' # (real data)
#' y <- na.omit(subset(ifit, id == 1)[,2:7])
#' p <- ncol(y)
#'
#' # random prior graph
#' # (dont do this in practice!!)
#' prior_var = matrix(sample(c(1,10),
#'                    size = p^2, replace = TRUE),
#'                    nrow = p, ncol = p)
#'
#' # fit model
#' fit <- prior_belief_var(y,
#'                         prior_temporal = prior_var,
#'                         post_odds_cut = 3)
#'}
prior_belief_var <- function(Y,
                             prior_temporal = NULL,
                             post_odds_cut,
                             est_ggm = TRUE,
                             prior_ggm = NULL,
                             progress = TRUE, ...){

  y <- Y

  prior_var <- prior_temporal

  p <- ncol(y)

  n <- nrow(y)

  y <- scale(y)

  colnames(y) <- NULL

  y_t <- as.matrix(y[-1,,drop=FALSE])

  y_t_1 <- as.matrix(y[-nrow(y),,drop=FALSE])

  colnames(y_t) <- 1:p

  BF_mat <- matrix(data = 0,
                   nrow = p,
                   ncol = p)

  if(est_ggm) {
    resid_mat <- matrix(data = 0,
                        nrow = n - 1,
                        ncol = p)
  }

  if(is.null(prior_var)) {
    prior_var <- matrix(1, p, p)
  }

  if(is.null(prior_ggm)){
    prior_ggm <- matrix(1, p, p)
  }

  if(any(prior_ggm == 0)){
    stop("zeros are not allowed in prior_ggm")
  }

  if(any(prior_var == 0)){
    stop("zeros are not allowed in prior_ggm")
  }

  message("testing temporal relations")

  if(progress){
    pb <- utils::txtProgressBar(min = 0,
                                 max = p,
                                 style = 3)
  }

  for(i in 1:p) {

    fit_i <- lm(y_t[, i] ~ 0 + y_t_1)

    if(est_ggm){
      resid_mat[,i] <- residuals(fit_i)
    }

    BF_10 <-
      lapply(names(coef(fit_i)), function(x) {
        1 / BF(fit_i, hypothesis = paste0(x, "=0"))$BFmatrix_confirmatory[1, 2]
      })

    BF_mat[, i] <- unlist(BF_10)

    if(progress){
      utils::setTxtProgressBar(pb, i)
    }

  }

  if (est_ggm) {
    message("\n\ntesting contemporanenous relations")

    fit_ggm <- prior_belief_ggm(Y = resid_mat,
                                prior_ggm = prior_ggm,
                                post_odds_cut = post_odds_cut,
                                ...)
  }

  post_odds <- t(BF_mat) * prior_var

  adj <- ifelse(post_odds > post_odds_cut, 1, 0)

  post_prob <- post_odds/(1 + post_odds)

  if(est_ggm) {
    returned_object <- list(
      adj_temporal = adj,
      post_prob_temporal = post_prob,
      adj_ggm = fit_ggm$adj,
      post_prob_ggm = fit_ggm$post_prob
    )
  } else {
    returned_object <- list(adj = adj,
                            post_prob = post_prob)
  }

  class(returned_object) <- c("BGGM", "prior_var")

  return(returned_object)
}

print_prior_var <- function(x, ...){
  cat("BGGM: Bayesian Gaussian Graphical Models \n")
  cat("Prior Belief Graphical VAR\n")
  cat("--- \n")
  cat("Date:", date(), "\n")
}
