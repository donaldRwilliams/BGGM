
ci_helper <- function(x, ci_width){
  low <- (1 - ci_width) / 2
  up  <-  1 - low
  interval <-  quantile(x, c(low, up))
  as.numeric(ifelse(interval[1] < 0 & interval[2] > 0, 0, 1))
}

inverse_2_beta <- function(fit, samples = 500){

  if(samples >= fit$iter){
    stop("Samples used to compute R2 cannot be greater than the number used for fitting the model")

  }
  # get posterior estimate for precision matrix
  inv <-  fit$posterior_samples[,  grep("cov_inv", colnames(fit$posterior_samples))]

  # seperate estimates by row
  node_wise_elements <- lapply(split(colnames(inv), 1:fit$p), function(x) inv[,x])

  # convert elemetns to regression counterparts
  betas <- lapply(1:fit$p, function(x) node_wise_elements[[x]][-x] *  as.matrix((1/ node_wise_elements[[x]][x]) ) * -1)
  sigmas <- lapply(1:fit$p, function(x) as.matrix((1/ node_wise_elements[[x]][x])))

  betas <- lapply(betas, function(x)  x[1:samples,])
  sigmas <- lapply(sigmas, function(x) sqrt(x[1:samples,]))
  returned_object <- list(betas = betas,  sigmas = sigmas, p = fit$p, data = fit$dat)
  class(returned_object) <- "inverse_2_beta"

  returned_object
}

summary_beta_helper <- function(x, node, ci_width){
  # index for row_names
  row_names <- 1:x$p

  # lower and upper ci
  low <- (1 - ci_width) / 2
  up <- 1 - low

  # beta posterior mean
  beta <- apply(x$betas[[node]], 2, mean)

  # beta poster sd
  post_sd <- apply(x$betas[[node]], 2, sd)

  # lower and upper of posterior
  beta_ci <- t(apply(x$betas[[node]], 2, quantile, probs = c(low, up)))

  # sd of the outcome
  sd_y <- sd(x$data[,node])

  # sd of the predictors
  sd_x <- apply(x$data[,-node], 2, sd)

  # standardized (std) beta
  beta_std_temp <- x$betas[[node]] * (sd_x / sd_y)

  # beta std posterior mean
  beta_std <- apply(beta_std_temp, 2, mean)

  # beta std posterior sd
  beta_std_post_sd <- apply(beta_std_temp, 2, sd)

  # lower and upper of posterior
  beta_std_ci <- t(apply(beta_std_temp, 2,  quantile, probs = c(low, up)))

  returned_object <- round(cbind(node =  row_names[-node],
                                 beta,
                                 post_sd,
                                 beta_ci,
                                 beta_std,
                                 post_sd = beta_std_post_sd,
                                 beta_std_ci),3)

  # remove row names
  row.names(returned_object ) <- NULL

  # make list for naming purposes
  returned_object <- list(as.data.frame(returned_object))

  # list elemenet name as the response
  names(returned_object) <- paste("predicting node", node)

  returned_object
}


R2_helper <- function(ypred, y, ci_width) {
  low <- (1 - ci_width) / 2
  up  <-  1 - low
  e <- -1 * sweep(ypred, 2, y)
  var_ypred <- apply(ypred, 1, var)
  var_e <- apply(e, 1, var)
  r2 <- unlist(var_ypred / (var_ypred + var_e))
  ci <- quantile(r2, prob = c(low, up) )
  mu_r2 <- mean(r2)
  sd_r2 <- sd(r2)
  summary_r2 <- c(post_mean = mu_r2, post_sd = sd_r2, ci)
  list(summary_r2 = summary_r2, R2 = r2)
}


name_helper <-  function(x){

  x <-  gsub("[A-z].*,", replacement = "", x)
  col_names <- gsub("[]]", "", x)
  col_names
}
