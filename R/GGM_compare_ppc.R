#' GGM comparison based on the predictive distribution
#'
#' @param g1 data from group 1 (n $\times$ p)
#' @param g2 data from group 1 (n $\times$ p)
#' @param samples numer of posterior and predictive samples
#'
#' @return
#' @export
#'
#' @examples
GGM_compare_ppc <- function(g1, g2, samples){
  # g1: network from group 1
  # g2: network from group 2
  # samples: number of samples

  # g1 sample size
  n <- nrow(g1)

  # g2 sample size
  p <- ncol(g1)

  # n2 sample size
  n2 <- nrow(g2)

  # g1 posterior distribution
  fit_g1 <- BGGM::bayes_estimate(g1, samples = samples)

  # g1 posterior mean
  inv_g1 <- fit_g1$inv_mat

  # g2 scatter matrix
  S_g2 <- t(g2) %*% g2

  # g2 precision martrix
  inv_g2 <- (n2 - 1) * solve(S_g2)

  # get names for matching
  col_names <-  colnames(fit_g1$posterior_samples)

  # posterior sample for precision matrix
  inv_cov <- fit_g1$posterior_samples[,grep(pattern = "cov_inv", x = col_names)]

  # compute predictive risk
  loss_fun <- lapply(1:samples, function(x) ppc_helper(x, inv_g1 = inv_g1, inv_cov, n =  n2, p))

  # retrieve predictive distributions
  # KL divergence
  KLD_ppc <- as.numeric(lapply(1:samples,function(x)  loss_fun[[x]]$KLD))

  # Jensen Shannon Divergence
  JSD_ppc <- as.numeric(lapply(1:samples,function(x)  loss_fun[[x]]$JSD))

  # Quadratic loss
  QL_ppc <- as.numeric(lapply(1:samples,function(x)  loss_fun[[x]]$QL))

  # frobenius loss
  FL_ppc <- as.numeric(lapply(1:samples,function(x)  loss_fun[[x]]$FL))

  # compute g1 vs g2 risk
  # KL divergence
  KLD_groups <- KL(Theta = inv_g1, hatTheta = inv_g2)

  # Jensen Shannon Divergence
  JSD_groups <- 0.5 *  KL(Theta = inv_g1, hatTheta = inv_g2) + 0.5 * KL(Theta = inv_g2, hatTheta = inv_g1)

  # Quadratic loss
  QL_groups <- QL(Theta = inv_g1, hatTheta = inv_g2)

  # frobenius loss
  FL_groups <- sum((inv_g1 - inv_g2)^2)


  # posterior predictive p-values
  # KL divergence
  KLD_pvalue <-  mean(KLD_ppc > KLD_groups)

  # Jensen Shannon Divergence
  JSD_pvalue <-  mean(JSD_ppc > JSD_groups)

  # Quadratic loss
  QL_pvalue  <-  mean(QL_ppc > QL_groups)

  # Quadratic loss
  FL_pvalue   <- mean(FL_ppc > FL_groups)

returned <- list(p_values = data.frame(loss = c("KL", "JSD", "QL","FL"),
                             p_value = c(KLD_pvalue,
                                         JSD_pvalue,
                                         QL_pvalue,
                                         FL_pvalue)),

                 score = data.frame(loss = c("KL", "JSD", "QL","FL"),
                                    score = c(KLD_groups,
                                              JSD_groups,
                                              QL_groups,
                                              FL_groups)),

                 predictive_dists = list(KLD_ppc = KLD_ppc,
                                         JSD_ppc = JSD_ppc,
                                         QL_ppc  = QL_ppc,
                                         FL_ppc = FL_ppc))

}

