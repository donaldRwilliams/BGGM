
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
  row_names <- 1:x$p
  low <- (1 - ci_width) / 2
  up <- 1 - low
  beta <- apply(x$betas[[node]], 2, mean)
  post_sd <- apply(x$betas[[node]], 2, sd)
  beta_ci <- t(apply(x$betas[[node]], 2, quantile, probs = c(low, up)))
  sd_y <- sd(x$data[,node])
  sd_x <- apply(x$data[,-node], 2, sd)
  beta_std_temp <- x$betas[[node]] * (sd_x / sd_y)
  beta_std <- apply(beta_std_temp, 2, mean)
  beta_std_post_sd <- apply(beta_std_temp, 2, sd)
  beta_std_ci <- t(apply(beta_std_temp, 2,  quantile, probs = c(low, up)))






  returned_object <- round(cbind(node =  row_names[-node],
                                 beta,
                                 post_sd,
                                 beta_ci,
                                 beta_std,
                                 post_sd = beta_std_post_sd,
                                 beta_std_ci),3)

  row.names(returned_object ) <- NULL

  returned_object <- list(as.data.frame(returned_object))

  names(returned_object) <- paste("predicting node", node)

  returned_object
}

