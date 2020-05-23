## ---- eval = FALSE------------------------------------------------------------
#  # need the developmental version
#  if (!requireNamespace("remotes")) {
#    install.packages("remotes")
#  }
#  
#  # install from github
#  remotes::install_github("donaldRwilliams/BGGM")
#  

## ---- warning =FALSE, message=FALSE-------------------------------------------
# need these packages
library(BGGM)
library(ggplot2)
library(assortnet)
library(networktools)
library(MASS)

# group 1
Yg1 <- MASS::mvrnorm(n = 926, 
                     mu = rep(0, 16), 
                     Sigma = ptsd_cor3, 
                     empirical = TRUE)

# group 2
Yg2 <- MASS::mvrnorm(n = 956, 
                     mu = rep(0, 16), 
                     Sigma = ptsd_cor4, 
                     empirical = TRUE)

## -----------------------------------------------------------------------------
f <- function(Yg1, Yg2){
  # number of nodes
  p <- ncol(Yg1)
  
  # index of off-diagonal
  indices <- upper.tri( diag(p))
  
  # group 1:
  # fit model
  g1_fit <- estimate(Yg1, analytic = TRUE)
  # pcors
  g1_pcors <- pcor_mat(g1_fit)[indices]
  
  # group 2
  # fit model
  g2_fit <- estimate(Yg2, analytic = TRUE)
  # pcors
  g2_pcors <- pcor_mat(g2_fit)[indices]
  
  # test-statistic
  cor(g1_pcors, g2_pcors)
  }

## -----------------------------------------------------------------------------
obs <- f(Yg1, Yg2)

# observed
obs

## ---- message=FALSE, results='hide'-------------------------------------------
ppc <- BGGM::ggm_compare_ppc(Yg1, Yg2, 
                       FUN = f, 
                       custom_obs = obs, 
                       iter = 1000, 
                       loss = FALSE)

## -----------------------------------------------------------------------------
ppc

## ---- eval=FALSE--------------------------------------------------------------
#  plot(ppc)

## ---- echo=FALSE, message=FALSE, warning=FALSE--------------------------------
plot(ppc, col_critical = "lightblue", 
     col_noncritical = "lightblue")[[1]] +
  xlab("Predictive Correlation")


## -----------------------------------------------------------------------------
f <- function(Yg1, Yg2){
  # nodes
  p <- ncol(Yg1)
  
  # index of off-diagonal
  indices <- upper.tri( diag(p))

  # fit models
  fit1 <-  BGGM::estimate(Yg1, analytic = TRUE)
  fit2 <-  BGGM::estimate(Yg2, analytic = TRUE)
  
  # select graphs
  sel1 <- BGGM::select(fit1)
  sel2 <- BGGM::select(fit2)
  
  # hamming distance
  sum((sel1$adj[indices] - sel2$adj[indices]) ^ 2)
}

## -----------------------------------------------------------------------------
obs <- f(Yg1, Yg2)

# observed
obs

## ---- message=FALSE, results='hide'-------------------------------------------
ppc <- BGGM::ggm_compare_ppc(Yg1, Yg2, 
                       FUN = f, 
                       custom_obs = obs, 
                       iter = 1000)

## -----------------------------------------------------------------------------
ppc

## ---- message=FALSE, warning=FALSE--------------------------------------------
plot(ppc)

## -----------------------------------------------------------------------------
f <- function(Yg1, Yg2){
  # nodes
  p <- ncol(Yg1)
  
  # index of off-diagonal
  indices <- upper.tri( diag(p))

  # fit models
  fit1 <-  BGGM::estimate(Yg1, analytic = TRUE)
  fit2 <-  BGGM::estimate(Yg2, analytic = TRUE)
  
  pcor1 <- BGGM::pcor_mat(fit1) 
  pcor2 <- BGGM::pcor_mat(fit2)
  
  # CDM for partial correlations
  # note: numerator is the trace; denominator is the Frobenius norm
  1 - (sum(diag(pcor1 %*% pcor2)) / (norm(pcor1, type = "f") * norm(pcor2, type = "f")))
}

## -----------------------------------------------------------------------------
obs <- f(Yg1, Yg2)

# observed
obs

## ---- message=FALSE, results='hide'-------------------------------------------
ppc <- BGGM::ggm_compare_ppc(Yg1, Yg2, 
                       FUN = f, 
                       custom_obs = obs, 
                       iter = 1000)

## -----------------------------------------------------------------------------
ppc

## -----------------------------------------------------------------------------
hist(ppc$predictive_custom, 
     xlim = c(0, obs), 
     main = "Partial Correlation Matrix Distance")
abline(v = obs)

## -----------------------------------------------------------------------------
# clusters based on DSM-5
comms <- c(
  rep("A", 4),
  rep("B", 7),
  rep("C", 5)
)

f <- function(Yg1, Yg2){

  fit1 <-  BGGM::estimate(Yg1, analytic = TRUE)
  fit2 <-  BGGM::estimate(Yg2, analytic = TRUE)

  pcor1 <- BGGM::pcor_mat(fit1)
  pcor2 <- BGGM::pcor_mat(fit2)

  assort1 <- assortnet::assortment.discrete(pcor1, types = comms,
                                         weighted = TRUE,
                                         SE = FALSE, M = 1)$r

  assort2 <- assortnet::assortment.discrete(pcor2, types = comms,
                                          weighted = TRUE,
                                          SE = FALSE, M = 1)$r
  (assort1 - assort2)
}

## -----------------------------------------------------------------------------
obs <- f(Yg1, Yg2)

# observed
obs

## ---- message=FALSE, results='hide'-------------------------------------------
ppc <- BGGM::ggm_compare_ppc(Yg1, Yg2, 
                       FUN = f, 
                       custom_obs = obs, 
                       iter = 1000)

## -----------------------------------------------------------------------------
ppc

## -----------------------------------------------------------------------------
plot(ppc)

## -----------------------------------------------------------------------------
f <- function(Yg1, Yg2){

  fit1 <-  BGGM::estimate(Yg1, analytic = TRUE)
  fit2 <-  BGGM::estimate(Yg2, analytic = TRUE)

  pcor1 <- BGGM::pcor_mat(fit1)
  pcor2 <- BGGM::pcor_mat(fit2)

  ei1 <- networktools::expectedInf(pcor1)$step1
 
  ei2 <- networktools::expectedInf(pcor2)$step1
   sum((ei1 - ei2)^2)
}

## -----------------------------------------------------------------------------
obs <- f(Yg1, Yg2)

# observed
obs

## ---- message=FALSE, results='hide'-------------------------------------------
ppc <- BGGM:::ggm_compare_ppc(Yg1, Yg2, 
                       FUN = f, 
                       custom_obs = obs, 
                       iter = 1000)

## -----------------------------------------------------------------------------
ppc

## -----------------------------------------------------------------------------
hist(ppc$predictive_custom, 
    xlim = c(0, obs),
     main = "Expected Influence\n Sum of Squared Error")
abline(v = obs)

