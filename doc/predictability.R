## ---- eval = FALSE, message=FALSE---------------------------------------------
#  # need the developmental version
#  if (!requireNamespace("remotes")) {
#    install.packages("remotes")
#  }
#  
#  # install from github
#  remotes::install_github("donaldRwilliams/BGGM")
#  library(BGGM)

## ---- eval=FALSE--------------------------------------------------------------
#  # binary data
#  Y <- women_math
#  
#  # fit model
#  fit <- estimate(Y, type = "binary")

## ---- eval=FALSE--------------------------------------------------------------
#  r2 <- predictability(fit)
#  
#  # print
#  r2
#  
#  #> BGGM: Bayesian Gaussian Graphical Models
#  #> ---
#  #> Metric: Bayes R2
#  #> Type: binary
#  #> ---
#  #> Estimates:
#  #>
#  #>  Node Post.mean Post.sd Cred.lb Cred.ub
#  #>     1     0.016   0.012   0.002   0.046
#  #>     2     0.103   0.023   0.064   0.150
#  #>     3     0.155   0.030   0.092   0.210
#  #>     4     0.160   0.021   0.118   0.201
#  #>     5     0.162   0.022   0.118   0.202
#  #>     6     0.157   0.028   0.097   0.208
#  #> ---

## ---- message=FALSE, eval=FALSE-----------------------------------------------
#  plot(r2,
#       type = "error_bar",
#       size = 4,
#       cred = 0.90)

## ---- message=FALSE, eval=FALSE-----------------------------------------------
#  plot(r2,
#       type = "ridgeline",
#       cred = 0.50)

## ---- eval=FALSE--------------------------------------------------------------
#  Y <- ptsd
#  
#  fit <- estimate(Y + 1, type = "ordinal")

## ---- eval=FALSE--------------------------------------------------------------
#  r2 <- predictability(fit)
#  
#  # print
#  r2
#  
#  #> BGGM: Bayesian Gaussian Graphical Models
#  #> ---
#  #> Metric: Bayes R2
#  #> Type: ordinal
#  #> ---
#  #> Estimates:
#  #>
#  #>  Node Post.mean Post.sd Cred.lb Cred.ub
#  #>     1     0.487   0.049   0.394   0.585
#  #>     2     0.497   0.047   0.412   0.592
#  #>     3     0.509   0.047   0.423   0.605
#  #>     4     0.524   0.049   0.441   0.633
#  #>     5     0.495   0.047   0.409   0.583
#  #>     6     0.297   0.043   0.217   0.379
#  #>     7     0.395   0.045   0.314   0.491
#  #>     8     0.250   0.042   0.173   0.336
#  #>     9     0.440   0.048   0.358   0.545
#  #>    10     0.417   0.044   0.337   0.508
#  #>    11     0.549   0.048   0.463   0.648
#  #>    12     0.508   0.048   0.423   0.607
#  #>    13     0.504   0.047   0.421   0.600
#  #>    14     0.485   0.043   0.411   0.568
#  #>    15     0.442   0.045   0.355   0.528
#  #>    16     0.332   0.039   0.257   0.414
#  #>    17     0.331   0.045   0.259   0.436
#  #>    18     0.423   0.044   0.345   0.510
#  #>    19     0.438   0.044   0.354   0.525
#  #>    20     0.362   0.043   0.285   0.454
#  #> ---

## ---- eval=FALSE--------------------------------------------------------------
#  plot(r2)

## ---- eval=FALSE--------------------------------------------------------------
#  # fit model
#  fit <- estimate(Y)
#  
#  # predictability
#  r2 <- predictability(fit)

