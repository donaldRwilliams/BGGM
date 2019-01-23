#' Title
#'
#' @param g1
#' @param g2
#' @param samples
#'
#' @return
#' @export
#'
#' @examples
GGM_compare_ppc <- function(g1, g2, samples){
  # g1: network from group 1
  #
  n <- nrow(g1)
  p <- ncol(g1)

  n2 <- nrow(g2)


  fit_g1 <- BGGM::bayes_estimate(g1, samples = samples)

  S_g2 <- t(g2) %*% g2
  inv_g2 <- (n2 - 1) * solve(S_g2)


  inv_g1 <- fit_g1$inv_mat

  col_names <-  colnames(fit_g1$posterior_samples)


  inv_cov <- fit_g1$posterior_samples[,grep(pattern = "cov_inv", x = col_names)]


  loss_fun <- lapply(1:samples, function(x) ppc_helper(x, inv_g1 = inv_g1, inv_cov ,n =  n2, p))

  KLD_ppc <- as.numeric(lapply(1:samples,function(x)  loss_fun[[x]]$KLD))
  JSD_ppc <- as.numeric(lapply(1:samples,function(x)  loss_fun[[x]]$JSD))
  QL_ppc <- as.numeric(lapply(1:samples,function(x)  loss_fun[[x]]$QL))
  FL_ppc <- as.numeric(lapply(1:samples,function(x)  loss_fun[[x]]$FL))

  KLD_groups <- KL(Theta = inv_g1, hatTheta = inv_g2)
  JSD_groups <- 0.5 *  KL(Theta = inv_g1, hatTheta = inv_g2) + 0.5 * KL(Theta = inv_g2, hatTheta = inv_g1)
  QL_groups <- QL(Theta = inv_g1, hatTheta = inv_g2)
  FL_groups <- sum((inv_g1 - inv_g2)^2)

  KLD_pvalue <-  mean(KLD_ppc > KLD_groups)
  JSD_pvalue <-  mean(JSD_ppc > JSD_groups)
  QL_pvalue  <-  mean(QL_ppc > QL_groups)
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
                               FL_ppc = FL_ppc)
       )

  }

