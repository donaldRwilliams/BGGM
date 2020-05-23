## ---- eval = FALSE------------------------------------------------------------
#  # need the developmental version
#  if (!requireNamespace("remotes")) {
#    install.packages("remotes")
#  }
#  
#  # install from github
#  remotes::install_github("donaldRwilliams/BGGM")

## ---- warning =FALSE, message=FALSE-------------------------------------------
# need these packages
library(BGGM)
library(ggplot2)
library(assortnet)
library(networktools)

# data
Y <- ptsd[,1:7]

## ---- message=FALSE, warning=FALSE, eval=FALSE--------------------------------
#  library(BGGM)
#  
#  # copula ggm
#  fit <- estimate(Y, type = "mixed", iter = 1000)

## -----------------------------------------------------------------------------
# define function
f <- function(x,...){
  networktools::expectedInf(x,...)$step1
}

## ---- eval = FALSE, message=FALSE, results='hide'-----------------------------
#  # iter = 250 for demonstrative purposes
#  # (but note even 1000 iters takes less than 1 second)
#  # compute
#  net_stat <- roll_your_own(object = fit,
#                            FUN = f,
#                            select = FALSE,
#                            iter = 250)
#  # print
#  net_stat
#  
#  #> BGGM: Bayesian Gaussian Graphical Models
#  #> ---
#  #> Network Stats: Roll Your Own
#  #> Posterior Samples: 250
#  #> ---
#  #> Estimates:
#  #>
#  #>  Node Post.mean Post.sd Cred.lb Cred.ub
#  #>     1     0.701   0.099   0.508   0.871
#  #>     2     0.912   0.113   0.722   1.179
#  #>     3     0.985   0.112   0.742   1.199
#  #>     4     1.056   0.105   0.851   1.247
#  #>     5     1.056   0.116   0.862   1.288
#  #>     6     0.491   0.092   0.329   0.679
#  #>     7     0.698   0.098   0.521   0.878
#  #> ---

## ---- eval = FALSE, results='hide'--------------------------------------------
#  net_stat <- roll_your_own(object = fit,
#                            FUN = f,
#                            select = TRUE,
#                            iter = 250)
#  
#  # print
#  net_stat
#  
#  #> BGGM: Bayesian Gaussian Graphical Models
#  #> ---
#  #> Network Stats: Roll Your Own
#  #> Posterior Samples: 250
#  #> ---
#  #> Estimates:
#  #>
#  #>  Node Post.mean Post.sd Cred.lb Cred.ub
#  #>     1     0.636   0.136   0.386   0.874
#  #>     2     0.792   0.113   0.580   0.996
#  #>     3     0.777   0.122   0.544   1.001
#  #>     4     0.910   0.121   0.667   1.143
#  #>     5     0.525   0.104   0.331   0.727
#  #>     6     0.484   0.110   0.270   0.686
#  #>     7     0.247   0.081   0.088   0.412
#  #> ---

## ---- message=FALSE, eval=FALSE-----------------------------------------------
#  plot(net_stat)

## ---- eval = FALSE, message=FALSE, results='hide'-----------------------------
#  # clusters
#  communities <- substring(colnames(Y), 1, 1)
#  
#  # function is slow
#  f <- function(x, ...){
#  networktools::bridge(x, ...)$`Bridge Strength`
#  }
#  
#  
#  # compute
#  net_stat <- roll_your_own(object = fit,
#                            FUN = f,
#                            communities = communities,
#                            iter = 250)
#  
#  # print
#  net_stat
#  
#  #> BGGM: Bayesian Gaussian Graphical Models
#  #> ---
#  #> Network Stats: Roll Your Own
#  #> Posterior Samples: 250
#  #> ---
#  #> Estimates:
#  #>
#  #>  Node Post.mean Post.sd Cred.lb Cred.ub
#  #>     1     0.162   0.082   0.035   0.347
#  #>     2     0.250   0.113   0.061   0.501
#  #>     3     0.180   0.104   0.049   0.480
#  #>     4     0.280   0.098   0.090   0.480
#  #>     5     0.375   0.093   0.196   0.558
#  #>     6     0.617   0.166   0.339   1.002
#  #>     7     0.628   0.166   0.400   1.025
#  #> ---

## ---- message = FALSE, eval=FALSE---------------------------------------------
#  plot(net_stat,
#       fill = "lightblue") +
#    ggtitle("Bridge Strength") +
#    xlab("Score")

## ---- eval = FALSE, message=FALSE, results='hide'-----------------------------
#  # clusters
#  communities <- substring(colnames(Y), 1, 1)
#  
#  # define function
#  f <- function(x,...){
#    assortnet::assortment.discrete(x, ...)$r
#  }
#  
#  net_stat <- roll_your_own(object = fit,
#                            FUN = f,
#                            types = communities,
#                            weighted = TRUE,
#                            SE = FALSE, M = 1,
#                            iter = 250)
#  
#  # print
#  net_stat
#  
#  #> BGGM: Bayesian Gaussian Graphical Models
#  #> ---
#  #> Network Stats: Roll Your Own
#  #> Posterior Samples: 250
#  #> ---
#  #> Estimates:
#  #>
#  #>  Post.mean Post.sd Cred.lb Cred.ub
#  #>      0.261   0.124   -0.01   0.469
#  #> ---

## ---- eval=FALSE--------------------------------------------------------------
#  hist(net_stat$results, main = "Assortment")

