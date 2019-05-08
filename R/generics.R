# generics
####################################
########## generics ################
####################################
select <- function(x,...){
  UseMethod("select", x)
}

edge_compare <- function(x, ...){
  UseMethod("edge_compare", x)
}

compare <- function(x, ...){
  UseMethod("compare", x)
}

loocv <- function(x, ...) UseMethod("loocv", x)

estimate <- function(x, ...) UseMethod("estimate")

explore  <- function(x, ...) UseMethod("explore")


confirm  <- function(x, ...) UseMethod("confirm")

loocv <- function(x, ...) UseMethod("loocv")

####################################
############ generics ##############
####################################
summary.compare.predict <- function(x, ...){
  cat("BGGM: Bayesian Gaussian Graphical Models \n")
  cat("--- \n")
  if(is.null(x$test_data)){
    cat("Type: In-sample predictive accuracy \n")
  } else{
    cat("Type: Out-of-sample predictive accuracy \n")
  }
  if(x$measure == "R2"){
    measure <-  "Variance Explained (R2) \n"
  } else{
    measure <-  "Mean Squared Error (MSE) \n"
  }
  cat("Measure:", measure)
  cat("Constrasts: Pairwise \n")
  cat("--- \n")
  cat("Call:\n")
  print(x$call)
  cat("--- \n")
  cat("Posterior Estimates: \n\n")
  print(x$summary_error, ...)

}


print.compare.predict <- function(x, ...){
  cat("BGGM: Bayesian Gaussian Graphical Models \n")
  cat("--- \n")
  if(is.null(x$test_data)){
    cat("Type: In-sample predictive accuracy contrasts \n")
  } else{
    cat("Type: Out-of-sample predictive accuracy contrasts \n")
  }
  if(x$measure == "R2"){
    measure <-  "Variance Explained (R2) \n"
  } else{
    measure <-  "Mean Squared Error (MSE) \n"
  }
  cat("Measure:", measure)
  cat("Contrasts: All pairwise \n")
  cat("--- \n")
  cat("Call:\n")
  print(x$call)
}



summary.loocv <- function(x,...){
  cat("BGGM: Bayesian Gaussian Graphical Models \n")
  cat("--- \n")
  if(is.numeric(x$samples)){
    cat("Type: Leave-One-Out Prediction Error (Bayesian) \n")
  } else {
    cat("Type: Leave-One-Out Prediction Error (Analytic) \n")
  }
  cat("--- \n")
  cat("Call:\n")
  print(x$call)
  cat("--- \n")
  cat("Estimates: \n\n ")
  temp <- x$returned_object
  rownames(temp) <- c()
  print(temp,  row.names = FALSE, ...)
  cat("--- \n")
  cat("Note: rss = residual sum of squares")

}


print.loocv <- function(x,...){
  cat("BGGM: Bayesian Gaussian Graphical Models \n")
  cat("--- \n")

  if(is.numeric(x$samples)){
    cat("Type: Leave-One-Out Prediction Error (Bayesian) \n")
  } else {
    cat("Type: Leave-One-Out Prediction Error (Analytic) \n")
  }
  cat("--- \n")
  cat("Call:\n")
  print(x$call)

}






analytic_solve <- function(X){
  # sample size
  n <- nrow(X)
  p <- ncol(X)
  # centererd mat
  X <- scale(X, scale = T)
  # scale matrix
  S <- t(X) %*% X
  inv_mu <-  solve(S + diag(10^-5,  p)) * (n)
  inv_var <-  (n + p + 1)*(solve(S + diag(10^-5, p) )^2 + tcrossprod(diag(solve(S + diag(10^-5, p)))))


  inv_cor <- diag( 1 / sqrt((diag(inv_mu)))) %*% inv_mu %*% diag( 1 / sqrt((diag(inv_mu))) )
  partials <- inv_cor * -1 + diag(2, p)
  list(inv_mu = inv_mu, inv_var = inv_var, partial = partials)
}








coef.estimate <- function(fit, node, ci_width,  samples = 1000){
  test <- beta_summary(fit, node = node, ci_width = ci_width, samples = samples)
  sums <- list()
  for(i in 1:max(node)){

    sums[[i]] <- test[[i]][[1]]
    names(sums)[[i]] <-  paste("Predicting node", node[i])
  }
  cat("BGGM: Bayesian Gaussian Graphical Models \n")
  cat("--- \n")
  cat("Type: Inverse to Regression \n")
  cat("--- \n")
  cat("Call: \n")
  print(test$call)
  cat("--- \n")
  cat("Posterior Estimates: \n \n")
  sums
}

summary.predict <- function(x,  ...){
  cat("BGGM: Bayesian Gaussian Graphical Models \n")
  cat("--- \n")
  if(is.null(x$test_data)){
    cat("Type: In-sample predictive accuracy \n")
  } else{
    cat("Type: Out-of-sample predictive accuracy \n")

  }
  if(x$measure == "R2"){
    measure <-  "Variance Explained (R2) \n"
  } else{
    measure <-  "Mean Squared Error (MSE) \n"

  }
  cat("Measure:", measure)
  cat("--- \n")
  cat("Call:\n")
  print(x$call)
  cat("--- \n")
  cat("Posterior Estimates: \n\n")
  temp <- cbind.data.frame(Node = 1:nrow(x$summary_error), x$summary_error)
  rownames(temp) <- c()
  print(temp,  row.names = FALSE, ...)
  cat("--- \n")
}


head.predict <- function(x,  ...){
  cat("BGGM: Bayesian Gaussian Graphical Models \n")
  cat("--- \n")
  if(is.null(x$test_data)){
    cat("Type: In-sample predictive accuracy \n")
  } else{
    cat("Type: Out-of-sample predictive accuracy \n")

  }
  if(x$measure == "R2"){
    measure <-  "Variance Explained (R2) \n"
  } else{
    measure <-  "Mean Squared Error (MSE) \n"

  }
  cat("Measure:", measure)
  cat("--- \n")
  cat("Call:\n")
  print(x$call)
  cat("--- \n")
  cat("Posterior Estimates: \n\n")
  temp <- cbind.data.frame(Node = 1:nrow(x$summary_error), x$summary_error)
  rownames(temp) <- c()
  print(temp[1:2,],  row.names = FALSE, ...)
  cat("--- \n")
}







print.predict <- function(x,  ...){
  cat("BGGM: Bayesian Gaussian Graphical Models \n")
  cat("--- \n")
  if(is.null(x$test_data)){
    cat("Type: In-sample predictive accuracy \n")
  } else{
    cat("Type: Out-of-sample predictive accuracy \n")

  }
  if(x$measure == "R2"){
    measure <-  "Variance Explained (R2) \n"
  } else{
    measure <-  "Mean Squared Error (MSE) \n"

  }
  cat("Measure:", measure)
  cat("--- \n")
  cat("Call:\n")
  print(x$call)
  # cat("--- \n")
  # cat("Posterior Estimates: \n\n")
  # print(x$summary_error, ...)
}

summary.estimate <- function(x){
  print(x)

}


print.estimate <- function(x){
  cat("BGGM: Bayesian Gaussian Graphical Models \n")
  cat("--- \n")
  if(!isFALSE( x$analytic)){
    cat("Type: Estimation (Analytic Solution) \n")
  }
  if(isFALSE( x$analytic)){
    cat("Type: Estimation (Sampling) \n")
  }
  cat("Posterior Samples:", x$iter, "\n")
  cat("Observations (n):", nrow(x$dat), "\n")
  cat("Variables (p):", x$p, "\n")
  cat("Edges:", .5 * (x$p * (x$p-1)), "\n")
  cat("--- \n")
  cat("Call: \n")
  print(x$call)
  cat("--- \n")
  cat("Date:", date(), "\n")
}


summary.select.estimate <- function(x, ...){
  cat("BGGM: Bayesian Gaussian Graphical Models \n")
  cat("--- \n")
  if(!isFALSE(x$analytic)){
    cat("Type: Selected Graph (Analytic Solution) \n")
  } else{
    cat("Type: Selected Graph (Sampling) \n")

  }
  if(is.null(x$rope)){
    cat("Credible Interval:", gsub("^.*\\.","", x$ci), "% \n")
    cat("Connectivity:", round(mean(x$adjacency[upper.tri(x$adjacency)]) * 100, 1), "% \n")
    cat("--- \n")
    cat("Call:\n")
    print(x$call)
    cat("--- \n")
    cat("Selected:\n \n")
    colnames( x$partials) <- 1:ncol(x$partials)
    row.names( x$partials) <- 1:ncol(x$partials)
    colnames( x$adjacency) <- 1:ncol(x$partials)
    row.names( x$adjacency) <- 1:ncol(x$partials)
    cat("Partial correlations \n \n")
    print(x$partials, digits = 2)
    cat("--- \n \n")
    cat("Adjacency \n \n")
    print(x$adjacency)
    cat("--- \n")

  } else{
    cat("Probability:", x$prob, "\n")
    cat("Region of Practical Equivalence:", "[", -1 * x$rope, ", ", x$rope, "]", "\n", sep = "")
    cat("Connectivity:", round(mean(x$adjacency_non_zero[upper.tri(x$adjacency_non_zero)]) * 100, 1), "% \n")
    cat("--- \n")
    cat("Call:\n")
    print(x$call)
    cat("--- \n")
    cat("Selected:\n \n")
    colnames(x$partials_non_zero) <- 1:ncol(x$partials_non_zero)
    row.names(x$partials_non_zero) <- 1:ncol(x$partials_non_zero)
    cat("Partial correlations \n \n")
    print(x$partials_non_zero, digits = 2)
    cat("--- \n \n")
    cat("Adjacency non-zero \n \n")
    colnames(x$adjacency_non_zero) <- 1:ncol(x$partials_non_zero)
    rownames(x$adjacency_non_zero) <- 1:ncol(x$partials_non_zero)
    print(x$adjacency_non_zero)
    cat("--- \n \n")
    cat("Adjacency zero \n \n")
    colnames(x$adjacency_zero) <- 1:ncol(x$partials_non_zero)
    rownames(x$adjacency_zero) <- 1:ncol(x$partials_non_zero)
    print(x$adjacency_zero)
  }
}



print.select.estimate <- function(x, ...){
  cat("BGGM: Bayesian Gaussian Graphical Models \n")
  cat("--- \n")
  if(is.numeric(x$rope)){
    cat("Type: Selected Graph (Sampling) \n")
  } else{
    cat("Type: Selected Graph (Analytic Solution) \n")

  }
  if(is.null(x$rope)){
    cat("Credible Interval:", gsub("^.*\\.","", x$ci), "% \n")
    cat("Connectivity:", round(mean(x$adjacency[upper.tri(x$adjacency)]) * 100, 1), "% \n")
    cat("--- \n")
    cat("Call:\n")
    print(x$call)
    cat("--- \n")


  } else{
    cat("Probability:", x$prob, "\n")
    cat("Region of Practical Equivalence:", "[", -1 * x$rope, ", ", x$rope, "]", "\n", sep = "")
    cat("Connectivity:", round(mean(x$adjacency_non_zero[upper.tri(x$adjacency_non_zero)]) * 100, 1), "% \n")
    cat("--- \n")
    cat("Call:\n")
    print(x$call)
    cat("--- \n")

  }
}

head.compare.estimate <- function(x){
  cat("BGGM: Bayesian Gaussian Graphical Models \n")
  cat("--- \n")
  cat("Type: Edge comparison(s) \n")
  cat("Credible Interval:", gsub("^.*\\.","", x$ci), "% \n")

  if(is.numeric(x$rope)){
    cat("Region of Practical Equivalence:", "[", -1 * x$rope, ", ", x$rope, "]", "\n", sep = "")
  }
  cat("--- \n")
  cat("Call:\n")
  print(x$call)
  cat("--- \n")
  print(edge_difference$returned_object[1,])
  cat("---")

}
print.compare.estimate <- function(x, ...){
  cat("BGGM: Bayesian Gaussian Graphical Models \n")
  cat("--- \n")
  cat("Type: Edge comparison(s) \n")
  cat("Credible Interval:", gsub("^.*\\.","", x$ci), "% \n")

  if(is.numeric(x$rope)){
    cat("Region of Practical Equivalence:", "[", -1 * x$rope, ", ", x$rope, "]", "\n", sep = "")
  }
  cat("--- \n")
  cat("Call:\n")
  print(x$call)
  cat("--- \n")
}


summary.compare.estimate <- function(x, ...){
  cat("BGGM: Bayesian Gaussian Graphical Models \n")
  cat("--- \n")
  cat("Type: Edge comparison(s) \n")
  cat("Credible Interval:", gsub("^.*\\.","", x$ci), "% \n")

  if(is.numeric(x$rope)){
    cat("Region of Practical Equivalence:", "[", -1 * x$rope, ", ", x$rope, "]", "\n", sep = "")
  }
  cat("--- \n")
  cat("Call:\n")
  print(x$call)
  cat("--- \n")
  cat("Posterior Estimates: \n\n")
  print(x$returned_object, ...)
}



print.explore <- function(x){
  cat("BGGM: Bayesian Gaussian Graphical Models \n")
  cat("--- \n")
  cat("Type: Hypothesis Testing (Exploratory) \n")
  cat("Posterior Samples:", x$iter, "\n")
  cat("Observations (n):", nrow(x$dat), "\n")
  cat("Variables (p):", x$p, "\n")
  cat("Edges:", .5 * (x$p * (x$p-1)), "\n")
  cat("Delta:", x$delta, "\n")
  cat("--- \n")
  cat("Call: \n")
  print(x$call)
  cat("--- \n")
  cat("Date:", date(), "\n")
}
summary.explore <- function(x){
  cat("BGGM: Bayesian Gaussian Graphical Models \n")
  cat("--- \n")
  cat("Type: Hypothesis Testing (Exploratory) \n")
  cat("Posterior Samples:", x$iter, "\n")
  cat("Observations (n):", nrow(x$dat), "\n")
  cat("Variables (p):", x$p, "\n")
  cat("Edges:", .5 * (x$p * (x$p-1)), "\n")
  cat("Delta:", x$delta, "\n")
  cat("--- \n")
  cat("Call: \n")
  print(x$call)
  cat("--- \n")
  cat("Date:", date(), "\n")
}
