## ---- eval = FALSE, message=FALSE---------------------------------------------
#  # need the developmental version
#  if (!requireNamespace("remotes")) {
#    install.packages("remotes")
#  }
#  
#  # install from github
#  remotes::install_github("donaldRwilliams/BGGM")
#  library(BGGM)

## ---- echo=FALSE, message=FALSE-----------------------------------------------
library(BGGM)

## ---- eval=FALSE--------------------------------------------------------------
#  # data
#  Y <- ptsd[,1:10]
#  
#  # fit model
#  # + 1 makes first category a 1
#  fit <- estimate(Y + 1, type = "ordinal")

## ---- eval=FALSE--------------------------------------------------------------
#  convergence(fit, print_names = TRUE)
#  
#  #>  [1] "B1--B2"         "B1--B3"         "B2--B3"         "B1--B4"         "B2--B4"         "B3--B4"         "B1--B5"
#  #>  [8] "B2--B5"         "B3--B5"         "B4--B5"         "B1--C1"         "B2--C1"         "B3--C1"         "B4--C1"
#  #> [15] "B5--C1"         "B1--C2"         "B2--C2"         "B3--C2"         "B4--C2"         "B5--C2"         "C1--C2"
#  #> [22] "B1--D1"         "B2--D1"         "B3--D1"         "B4--D1"         "B5--D1"         "C1--D1"         "C2--D1"
#  #> [29] "B1--D2"         "B2--D2"         "B3--D2"         "B4--D2"         "B5--D2"         "C1--D2"         "C2--D2"
#  #> [36] "D1--D2"         "B1--D3"         "B2--D3"         "B3--D3"         "B4--D3"         "B5--D3"         "C1--D3"
#  #> [43] "C2--D3"         "D1--D3"         "D2--D3"         "B1_(Intercept)" "B2_(Intercept)" "B3_(Intercept)" "B4_(Intercept)"
#  #> [50] "B5_(Intercept)" "C1_(Intercept)" "C2_(Intercept)" "D1_(Intercept)" "D2_(Intercept)" "D3_(Intercept)"

## ---- eval=FALSE--------------------------------------------------------------
#  convergence(fit, param = "B1--B2", type = "acf")

## ---- eval=FALSE--------------------------------------------------------------
#  # sim time series
#  ts.sim <- arima.sim(list(order = c(1,1,0), ar = 0.7), n = 200)
#  
#  acf(ts.sim)

## ---- eval=FALSE--------------------------------------------------------------
#  # extract samples
#  samps <- fit$post_samp$pcors
#  
#  # iterations
#  iter <- fit$iter
#  
#  # thinning interval
#  thin <-  5
#  
#  # save every 5th (add 50 which is the burnin)
#  new_iter <- length(seq(1,to = iter + 50 , by = thin))
#  
#  # replace (add 50 which is the burnin)
#  fit$post_samp$pcors <- samps[,,seq(1,to = iter + 50, by = thin)]
#  
#  # replace iter
#  fit$iter <- new_iter - 50
#  
#  # check thinned
#  convergence(fit, param = "B1--B2", type = "acf")

## ---- eval=FALSE--------------------------------------------------------------
#  convergence(fit, param = "B1--B2", type = "trace")

