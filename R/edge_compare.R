#' Title
#'
#' @param X1
#' @param X2
#' @param ROPE
#' @param decision_rule
#' @param chains
#' @param iter
#' @param burnin
#'
#' @return
#' @export
#'
#' @examples
edge_compare <- function(X1, X2, ROPE, decision_rule, chains = 2, iter = 1000, burnin = 100){
  rope_helper <- function(x, ROPE){
    mean(x > -ROPE & x < ROPE)

  }
  p <- ncol(X1)
  mat_BF <- mat_selected<- matrix(0, p, p)
  prior_samps <- replicate(chains * iter, cov2cor(rWishart(n = 1, df = p, diag(2))[,,1])[1,2] -  cov2cor(rWishart(n = 1, df = p, diag(2))[,,1])[1,2] )
  fit1 <- BGGM(X1, chains = chains, iter = iter, burnin = burnin)
  fit2 <- BGGM(X2, chains = chains, iter = iter, burnin = burnin)

  difference <- fit1$posterior_samples - fit2$posterior_samples

  pcs <- difference[,c(p*p+1):ncol(difference)]

  prob_pcs <- apply(pcs, 2, rope_helper, ROPE)
  prob_prior <- rope_helper(prior_samps, ROPE)

  BF_01 <- (prob_pcs / prob_prior) * (prob_prior / ROPE)
  mat_BF[] <- BF_01
  mat_selected <- ifelse(mat_BF > decision_rule, 1, 0)
  list(mat_BF = mat_BF, mat_selected = mat_selected)

}
