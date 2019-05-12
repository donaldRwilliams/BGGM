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
summary.edge_compare.predict <- function(x, ...){
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
  cat("pr_out: post prob outside of rope \n")
  cat("pr_in: post prob inside of rope \n")
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








# coef.estimate <- function(fit, node, ci_width,  samples = 1000){
#   test <- beta_summary(fit, node = node, ci_width = ci_width, samples = samples)
#   sums <- list()
#   for(i in 1:max(node)){
#
#     sums[[i]] <- test[[i]][[1]]
#     names(sums)[[i]] <-  paste("Predicting node", node[i])
#   }
#   cat("BGGM: Bayesian Gaussian Graphical Models \n")
#   cat("--- \n")
#   cat("Type: Inverse to Regression \n")
#   cat("--- \n")
#   cat("Call: \n")
#   print(test$call)
#   cat("--- \n")
#   cat("Posterior Estimates: \n \n")
#   name(sums) <- ""
#   sums
# }

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
  cat("Estimates: \n\n")
  temp <- cbind.data.frame(node = 1:nrow(x$summary_error), x$summary_error)
  rownames(temp) <- c()
  print(temp,  row.names = FALSE, ...)
  cat("--- \n")
}


head.predict <- function(x, nrow,  ...){
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
  cat("Estimates: \n\n")
  temp <- cbind.data.frame(node = 1:nrow(x$summary_error), x$summary_error)
  rownames(temp) <- c()
  print(temp[1:nrow,],  row.names = FALSE, ...)
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


summary.select.estimate <- function(x, summarize = F, ...){
  cat("BGGM: Bayesian Gaussian Graphical Models \n")
  cat("--- \n")
  if(!is.null(x$analytic)){
    cat("Type: Selected Graph (Analytic Solution) \n")
  } else{
    cat("Type: Selected Graph (Sampling) \n")

  }

  if(isFALSE(summarize)){
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
  if(isTRUE(summarize)){
    if(isTRUE(x$analytic)){
      stop("summary not available for the analytic solution")
    }


    if(is.null(x$rope)){
      p <- ncol(x$partials)
      mat_names <- mu_mat <- ci_low <- ci_up <- mat_temp <- matrix(0, p, p)
      mat_names[] <-  unlist(lapply(1:p, function(z) paste(1:p, z, sep = "--")))


      low <- (1 - x$ci) / 2
      up  <-  1 - low

      mu_mat[] <-  colMeans(x$pcor_samples)
      sd <-  x$pcor_sd[upper.tri(x$pcor_sd)]
      cis <- apply(x$pcor_samples, 2, quantile, c(low, up))
      ci_low[] <- cis[1,]
      ci_up[] <- cis[2,]

      summ <- data.frame(egde = mat_names[upper.tri(mat_names)],
                         post_mean = mu_mat[upper.tri(mu_mat)],
                         post_sd = x$pcor_sd[upper.tri(x$pcor_sd)],
                         temp1 = ci_low[upper.tri(ci_low)],
                         temp2 = ci_up[upper.tri(ci_up)],
                         check.names = F)

      colnames(summ)[4] <-  row.names(cis)[1]
      colnames(summ)[5] <- row.names(cis)[2]
      cat("Credible Interval:", gsub("^.*\\.","", x$ci), "% \n")
      cat("Connectivity:", round(mean(x$adjacency[upper.tri(x$adjacency)]) * 100, 1), "% \n")
      cat("--- \n")
      cat("Call:\n")
      print(x$call)
      cat("--- \n")
      cat("Estimates: \n \n")
      print(summ, row.names = F,...)
      cat("--- \n")

    }else{
      cat("Probability:", x$prob, "\n")
      cat("Region of Practical Equivalence:", "[", -1 * x$rope, ", ", x$rope, "]", "\n", sep = "")
      cat("Connectivity:", round(mean(x$adjacency_non_zero[upper.tri(x$adjacency_non_zero)]) * 100, 1), "% \n")
      cat("--- \n")
      cat("Call:\n")
      print(x$call)
      cat("--- \n")
      cat("pr_out: post prob outside of rope \n")
      cat("pr_in: post prob inside of rope \n")
      cat("--- \n")


      p <- ncol(x$partials)
      mat_names <- mu_mat <- rope_in  <- matrix(0, p, p)
      mat_names[] <-  unlist(lapply(1:p, function(z) paste(1:p, z, sep = "--")))


      low <- (1 - x$ci) / 2
      up  <-  1 - low

      mu_mat[] <-  colMeans(x$pcor_samples)
      sd <-  x$pcor_sd[upper.tri(x$pcor_sd)]

      rope_in[] <- x$in_rope


      cat("Estimates: \n \n")
      summ <- data.frame(egde = mat_names[upper.tri(mat_names)],
                         post_mean = mu_mat[upper.tri(mu_mat)],
                         post_sd = x$pcor_sd[upper.tri(x$pcor_sd)],
                         "pr_out" = 1 - rope_in[upper.tri(rope_in)],
                         "pr_in" = rope_in[upper.tri(rope_in)],
                         check.names = F)




      print(summ, row.names = F,...)
      cat("--- \n")
    }
  }
}











head.select.estimate <- function(x, summarize = F, nrow = 2, ...){
  cat("BGGM: Bayesian Gaussian Graphical Models \n")
  cat("--- \n")
  if(!is.null(x$analytic)){
    cat("Type: Selected Graph (Analytic Solution) \n")
  } else{
    cat("Type: Selected Graph (Sampling) \n")

  }

  if(isFALSE(summarize)){
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
  if(isTRUE(summarize)){
    if(isTRUE(x$analytic)){
      stop("summary not available for the analytic solution")
    }


    if(is.null(x$rope)){
      p <- ncol(x$partials)
      mat_names <- mu_mat <- ci_low <- ci_up <- mat_temp <- matrix(0, p, p)
      mat_names[] <-  unlist(lapply(1:p, function(z) paste(1:p, z, sep = "--")))


      low <- (1 - x$ci) / 2
      up  <-  1 - low

      mu_mat[] <-  colMeans(x$pcor_samples)
      sd <-  x$pcor_sd[upper.tri(x$pcor_sd)]
      cis <- apply(x$pcor_samples, 2, quantile, c(low, up))
      ci_low[] <- cis[1,]
      ci_up[] <- cis[2,]

      summ <- data.frame(egde = mat_names[upper.tri(mat_names)],
                         post_mean = mu_mat[upper.tri(mu_mat)],
                         post_sd = x$pcor_sd[upper.tri(x$pcor_sd)],
                         temp1 = ci_low[upper.tri(ci_low)],
                         temp2 = ci_up[upper.tri(ci_up)],
                         check.names = F)

      colnames(summ)[4] <-  row.names(cis)[1]
      colnames(summ)[5] <- row.names(cis)[2]
      cat("Credible Interval:", gsub("^.*\\.","", x$ci), "% \n")
      cat("Connectivity:", round(mean(x$adjacency[upper.tri(x$adjacency)]) * 100, 1), "% \n")
      cat("--- \n")
      cat("Call:\n")
      print(x$call)
      cat("--- \n")
      cat("Estimates: \n \n")
      print(summ[1:nrow], row.names = F,...)
      cat("--- \n")

    }else{
      cat("Probability:", x$prob, "\n")
      cat("Region of Practical Equivalence:", "[", -1 * x$rope, ", ", x$rope, "]", "\n", sep = "")
      cat("Connectivity:", round(mean(x$adjacency_non_zero[upper.tri(x$adjacency_non_zero)]) * 100, 1), "% \n")
      cat("--- \n")
      cat("Call:\n")
      print(x$call)
      cat("--- \n")
      cat("pr_out: post prob outside of rope \n")
      cat("pr_in: post prob inside of rope \n")
      cat("--- \n")


      p <- ncol(x$partials)
      mat_names <- mu_mat <- rope_in  <- matrix(0, p, p)
      mat_names[] <-  unlist(lapply(1:p, function(z) paste(1:p, z, sep = "--")))


      low <- (1 - x$ci) / 2
      up  <-  1 - low

      mu_mat[] <-  colMeans(x$pcor_samples)
      sd <-  x$pcor_sd[upper.tri(x$pcor_sd)]

      rope_in[] <- x$in_rope


      cat("Estimates: \n \n")
      summ <- data.frame(egde = mat_names[upper.tri(mat_names)],
                         post_mean = mu_mat[upper.tri(mu_mat)],
                         post_sd = x$pcor_sd[upper.tri(x$pcor_sd)],
                         "pr_out" = 1 - rope_in[upper.tri(rope_in)],
                         "pr_in" = rope_in[upper.tri(rope_in)],
                         check.names = F)




      print(summ[1:nrow,], row.names = F,...)
      cat("--- \n")
    }
  }
}



# summary.select.estimate <- function(x, ...){
#   cat("BGGM: Bayesian Gaussian Graphical Models \n")
#   cat("--- \n")
#   if(!isFALSE(x$analytic)){
#     cat("Type: Selected Graph (Analytic Solution) \n")
#   } else{
#     cat("Type: Selected Graph (Sampling) \n")
#
#   }
#   if(is.null(x$rope)){
#     cat("Credible Interval:", gsub("^.*\\.","", x$ci), "% \n")
#     cat("Connectivity:", round(mean(x$adjacency[upper.tri(x$adjacency)]) * 100, 1), "% \n")
#     cat("--- \n")
#     cat("Call:\n")
#     print(x$call)
#     cat("--- \n")
#     cat("Selected:\n \n")
#     colnames( x$partials) <- 1:ncol(x$partials)
#     row.names( x$partials) <- 1:ncol(x$partials)
#     colnames( x$adjacency) <- 1:ncol(x$partials)
#     row.names( x$adjacency) <- 1:ncol(x$partials)
#     cat("Partial correlations \n \n")
#     print(x$partials, digits = 2)
#     cat("--- \n \n")
#     cat("Adjacency \n \n")
#     print(x$adjacency)
#     cat("--- \n")
#
#   } else{
#     cat("Probability:", x$prob, "\n")
#     cat("Region of Practical Equivalence:", "[", -1 * x$rope, ", ", x$rope, "]", "\n", sep = "")
#     cat("Connectivity:", round(mean(x$adjacency_non_zero[upper.tri(x$adjacency_non_zero)]) * 100, 1), "% \n")
#     cat("--- \n")
#     cat("Call:\n")
#     print(x$call)
#     cat("--- \n")
#     cat("Selected:\n \n")
#     colnames(x$partials_non_zero) <- 1:ncol(x$partials_non_zero)
#     row.names(x$partials_non_zero) <- 1:ncol(x$partials_non_zero)
#     cat("Partial correlations \n \n")
#     print(x$partials_non_zero, digits = 2)
#     cat("--- \n \n")
#     cat("Adjacency non-zero \n \n")
#     colnames(x$adjacency_non_zero) <- 1:ncol(x$partials_non_zero)
#     rownames(x$adjacency_non_zero) <- 1:ncol(x$partials_non_zero)
#     print(x$adjacency_non_zero)
#     cat("--- \n \n")
#     cat("Adjacency zero \n \n")
#     colnames(x$adjacency_zero) <- 1:ncol(x$partials_non_zero)
#     rownames(x$adjacency_zero) <- 1:ncol(x$partials_non_zero)
#     print(x$adjacency_zero)
#   }
# }
#
#
#

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

head.edge_compare.estimate <- function(x, nrow = 2, ...){
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
  cat("Estimates: \n \n")
  print(edge_difference$returned_object[1:nrow,-c(4:5)], row.names = FALSE, ...)
  cat("---")

}
print.edge_compare.estimate <- function(x, ...){
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


summary.edge_compare.estimate <- function(x, ...){
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
  print(x$returned_object[,-c(4:5)], row.names = FALSE,  ...)
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



print.select.explore <- function(x, hyp = "H1",  log = TRUE, summarize = FALSE, ...){

  # name of package
  cat("BGGM: Bayesian Gaussian Graphical Models \n")
  cat("--- \n")

  # hypothesis testing
  cat("Type: Hypothesis Testing \n")

  # alternative argument
  cat("Alternative:", sel$alternative, "\n")

  # if exhaustive hypothesis testing
  if(x$alternative == "exhaustive"){

    # posterior probability
    cat("Posterior probability:", x$prob, "\n")
    cat("--- \n")
    cat("Call:\n")
    print(x$call)
    cat("--- \n")
  }
  if(x$alternative == "greater"){
    cat("Bayes Factor:", x$BF_cut, "\n")
    if(isFALSE(summarize)){
      cat("Connectivity:", round(mean(x$Adj_20[upper.tri(x$Adj_20)]) * 100, 1), "% \n")

    }
    cat("Call:\n")
    print(x$call)
    cat("--- \n")
  }

  # if less than zero
  if(x$alternative == "less"){
    cat("Bayes Factor:", x$BF_cut, "\n")
    if(isFALSE(summarize)){
      cat("Connectivity:", round(mean(x$Adj_20[upper.tri(x$Adj_20)]) * 100, 1), "% \n")

    }
    cat("Call:\n")
    print(x$call)
    cat("--- \n")
  }


}

summary.select.explore <- function(x, hyp = "H1",  log = TRUE, summarize = FALSE, ...){
  cat("BGGM: Bayesian Gaussian Graphical Models \n")
  cat("--- \n")
  cat("Type: Hypothesis Testing \n")
  cat("Alternative:", x$alternative, "\n")

  # exhaustive
  if(x$alternative == "exhaustive"){
    cat("Posterior probability:", x$prob, "\n")
    cat("--- \n")
    cat("Call:\n")
    print(x$call)
    cat("--- \n")
    if(!isFALSE(summarize)){

      summ <-  cbind.data.frame(edge = x$post_prob$edge,
                                post_mean = x$pcor_mat[upper.tri(x$pcor_mat)],
                                post_sd = x$pcor_sd[upper.tri(x$pcor_sd)],
                                "p(H0|Y)" = x$post_prob[,2],
                                "p(H1|Y)" = x$post_prob[,3],
                                "p(H2|Y)" = x$post_prob[,4])

      cat("Hypotheses: \n")
      cat("H0: rho = 0\nH1: rho > 0\nH2: rho < 0", "\n")
      cat("--- \n")
      cat("Estimates: \n \n ")
      print(summ, row.names = FALSE, ...)
      cat("--- \n")
    } else{
      if(hyp == "H1"){
        cat("Hypothesis: \n")
        cat("H1: rho > 0 \n")
        cat("--- \n")
        p <- ncol(x$pos_mat)
        colnames(x$pos_mat) <- 1:p
        row.names(x$pos_mat) <- 1:p
        cat("Partial Correlations \n \n")
        print(x$pcor_mat * ifelse(x$pos_mat > x$prob, 1, 0), ...)
        cat("--- \n \n")
        cat("Adjancency \n \n")
        print(ifelse(x$pos_mat > x$prob, 1, 0), ...)
        cat("--- \n")

      }

      if(hyp == "H0"){
        cat("Hypothesis: \n")
        cat("H0: rho = 0 \n")
        cat("--- \n")
        p <- ncol(x$null_mat)
        colnames(x$null_mat) <- 1:p
        row.names(x$null_mat) <- 1:p
        cat("Partial Correlations \n \n")
        print(x$pcor_mat * ifelse(x$null_mat > x$prob, 1, 0), ...)
        cat("--- \n \n")
        cat("Adjancency \n \n")
        print(ifelse(x$null_mat > x$prob, 1, 0), ...)
        cat("--- \n")

      }

      if(hyp == "H2"){
        cat("Hypothesis: \n")
        cat("H0: rho < 0 \n")
        cat("--- \n")
        p <- ncol(x$neg_mat)
        colnames(x$neg_mat) <- 1:p
        row.names(x$neg_mat) <- 1:p
        cat("Partial Correlations \n \n")
        print(x$pcor_mat * ifelse(x$neg_mat > x$prob, 1, 0), ...)
        cat("--- \n \n")
        cat("Adjancency \n \n")
        print(ifelse(x$neg_mat > x$prob, 1, 0), ...)
        cat("--- \n")
      }
    }
  }

  # greater than
  if(x$alternative == "greater"){
    cat("Bayes Factor:", x$BF_cut, "\n")
    if(isFALSE(summarize)){
      cat("Connectivity:", round(mean(x$Adj_20[upper.tri(x$Adj_20)]) * 100, 1), "% \n")

    }
    cat("Call:\n")
    print(x$call)
    cat("--- \n")

    if(!isFALSE(summarize)){
      p <- ncol(x$pcor_mat)
      mat_names <- matrix(0, p, p)

      edge_name <-

        mat_names[] <-  unlist(lapply(1:p, function(z) paste(1:p, z, sep = "--")))

      if(log == TRUE){
        summ <-  cbind.data.frame(edge = mat_names[upper.tri(mat_names)],
                                  post_mean = x$pcor_mat[upper.tri(x$pcor_mat)],
                                  post_sd = x$pcor_sd[upper.tri(x$pcor_sd)],
                                  "BF 20" = log(x$BF_20[upper.tri(x$BF_20)]),
                                  "BF 01" = log(x$BF_01[upper.tri(x$BF_01)]))
      } else{

        summ <-  cbind.data.frame(edge = mat_names[upper.tri(mat_names)],
                                  post_mean = x$pcor_mat[upper.tri(x$pcor_mat)],
                                  post_sd = x$pcor_sd[upper.tri(x$pcor_sd)],
                                  "BF 20" = x$BF_20[upper.tri(x$BF_20)],
                                  "BF 01" = x$BF_01[upper.tri(x$BF_01)])
      }

      cat("Hypotheses: \n")
      cat("H0: rho = 0\nH1: rho > 0", "\n")
      cat("--- \n")
      cat("Estimates: \n \n ")
      print(summ, row.names = FALSE, ...)
      cat("--- \n")
      cat("note: BF 20 is a one-sided Bayes factor for H1 \n")
      cat("--- \n")
    } else{
      if(hyp == "H1"){

        cat("Hypothesis: \n")
        cat("H1: rho > 0 \n")
        cat("--- \n")
        p <- ncol(x$partials_positive)
        colnames(x$partials_positive) <- 1:p
        row.names(x$partials_positive) <- 1:p
        cat("Partial Correlations \n \n")
        print(x$partials_positive, ...)
        cat("--- \n \n")
        cat("Adjancency (positive) \n \n")
        colnames(x$Adj_20) <- 1:p
        row.names(x$Adj_20) <- 1:p
        print(x$Adj_20, ...)
        cat("--- \n")
      }


      if(hyp == "H0"){

        cat("Hypothesis: \n")
        cat("H0: rho = 0 \n")
        cat("--- \n")
        p <- ncol(x$pcor_mat)
        colnames(x$pcor_mat) <- 1:p
        row.names(x$pcor_mat) <- 1:p
        cat("Partial Correlations \n \n")
        print(x$pcor_mat * x$Adj_01, ...)
        cat("--- \n \n")
        cat("Adjancency (null) \n \n")
        colnames(x$Adj_01) <- 1:p
        row.names(x$Adj_01) <- 1:p
        print(x$Adj_01, ...)
        cat("--- \n")
      }
    }


  }


  # less than
  if(x$alternative == "less"){
    cat("Bayes Factor:", x$BF_cut, "\n")
    if(isFALSE(summarize)){
      cat("Connectivity:", round(mean(x$Adj_20[upper.tri(x$Adj_20)]) * 100, 1), "% \n")

    }
    cat("--- \n")
    cat("Call:\n")
    print(x$call)
    cat("--- \n")


    if(!isFALSE(summarize)){
      p <- ncol(x$pcor_mat)
      mat_names <- matrix(0, p, p)

      edge_name <-

        mat_names[] <-  unlist(lapply(1:p, function(z) paste(1:p, z, sep = "--")))

      if(log == TRUE){
        summ <-  cbind.data.frame(edge = mat_names[upper.tri(mat_names)],
                                  post_mean = x$pcor_mat[upper.tri(x$pcor_mat)],
                                  post_sd = x$pcor_sd[upper.tri(x$pcor_sd)],
                                  "BF 20" = log(x$BF_20[upper.tri(x$BF_20)]),
                                  "BF 01" = log(x$BF_01[upper.tri(x$BF_01)]))
      } else{

        summ <-  cbind.data.frame(edge = mat_names[upper.tri(mat_names)],
                                  post_mean = x$pcor_mat[upper.tri(x$pcor_mat)],
                                  post_sd = x$pcor_sd[upper.tri(x$pcor_sd)],
                                  "BF 20" = x$BF_20[upper.tri(x$BF_20)],
                                  "BF 01" = x$BF_01[upper.tri(x$BF_01)])
      }

      cat("Hypotheses: \n")
      cat("H0: rho = 0\nH1: rho < 0", "\n")
      cat("--- \n")
      cat("Estimates: \n \n ")
      print(summ, row.names = FALSE, ...)
      cat("--- \n")
      cat("note: BF 20 is a one-sided Bayes factor for H1 \n")
      cat("--- \n")
    }

    else{
      if(hyp == "H1"){

        cat("Hypothesis: \n")
        cat("H1: rho < 0 \n")
        cat("--- \n")
        p <- ncol(x$partials_negative)
        colnames(x$partials_negative) <- 1:p
        row.names(x$partials_negative) <- 1:p
        cat("Partial Correlations \n \n")
        print(x$partials_negative, ...)
        cat("--- \n \n")
        cat("Adjancency (positive) \n \n")
        colnames(x$Adj_20) <- 1:p
        row.names(x$Adj_20) <- 1:p
        print(x$Adj_20, ...)
        cat("--- \n")
      }


      if(hyp == "H0"){

        cat("Hypothesis: \n")
        cat("H0: rho = 0 \n")
        cat("--- \n")
        p <- ncol(x$pcor_mat)
        colnames(x$pcor_mat) <- 1:p
        row.names(x$pcor_mat) <- 1:p
        cat("Partial Correlations \n \n")
        print(x$pcor_mat * x$Adj_01, ...)
        cat("--- \n \n")
        cat("Adjancency (null) \n \n")
        colnames(x$Adj_01) <- 1:p
        row.names(x$Adj_01) <- 1:p
        print(x$Adj_01, ...)
        cat("--- \n")
      }
    }

  }

  if(x$alternative == "two.sided"){
    cat("Bayes Factor:", x$BF_cut, "\n")
    if(isFALSE(summarize)){
      cat("Connectivity:", round(mean(x$Adj_10[upper.tri(x$Adj_10)]) * 100, 1), "% \n")

    }
    cat("--- \n")
    cat("Call:\n")
    print(x$call)
    cat("--- \n")


    if(!isFALSE(summarize)){
      p <- ncol(x$pcor_mat)
      mat_names <- matrix(0, p, p)



      mat_names[] <-  unlist(lapply(1:p, function(z) paste(1:p, z, sep = "--")))

      if(log == TRUE){
        summ <-  cbind.data.frame(edge = mat_names[upper.tri(mat_names)],
                                  post_mean = x$pcor_mat[upper.tri(x$pcor_mat)],
                                  post_sd = x$pcor_sd[upper.tri(x$pcor_sd)],
                                  "BF 10" = log(x$BF_10[upper.tri(x$BF_10)]),
                                  "BF 01" = log(x$BF_01[upper.tri(x$BF_01)]))
      } else{

        summ <-  cbind.data.frame(edge = mat_names[upper.tri(mat_names)],
                                  post_mean = x$pcor_mat[upper.tri(x$pcor_mat)],
                                  post_sd = x$pcor_sd[upper.tri(x$pcor_sd)],
                                  "BF 20" = x$BF_10[upper.tri(x$BF_10)],
                                  "BF 01" = x$BF_01[upper.tri(x$BF_01)])
      }


      cat("Hypotheses: \n")
      cat("H0: rho = 0\nH1: rho != 0", "\n")
      cat("--- \n")
      cat("Estimates: \n \n ")
      print(summ, row.names = FALSE, ...)
      cat("--- \n")
    }
    else{

      if(hyp == "H1"){

        cat("Hypothesis: \n")
        cat("H1: rho != 0 \n")
        cat("--- \n")
        p <- ncol(x$partials_non_zero)
        colnames(x$partials_non_zero) <- 1:p
        row.names(x$partials_non_zero) <- 1:p
        cat("Partial Correlations \n \n")
        print(x$partials_non_zero)
        cat("--- \n \n")
        cat("Adjancency (non-zero) \n \n")
        colnames(x$Adj_10) <- 1:p
        row.names(x$Adj_10) <- 1:p
        print(x$Adj_10)
        cat("--- \n")
      }


      if(hyp == "H0"){

        cat("Hypothesis: \n")
        cat("H0: rho = 0 \n")
        cat("--- \n")
        p <- ncol(x$pcor_mat)
        colnames(x$pcor_mat) <- 1:p
        row.names(x$pcor_mat) <- 1:p
        cat("Partial Correlations \n \n")
        print(x$pcor_mat * x$Adj_01, ...)
        cat("--- \n \n")
        cat("Adjancency (null) \n \n")
        colnames(x$Adj_01) <- 1:p
        row.names(x$Adj_01) <- 1:p
        print(x$Adj_01, ...)
        cat("--- \n")
      }
    }




  }

}


