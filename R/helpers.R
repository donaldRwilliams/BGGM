#' @importFrom stats coef cov2cor dnorm lm na.omit pnorm quantile rWishart sd qnorm
#' @importFrom utils combn
#' @importFrom foreach %dopar% foreach
#' @import ggplot2
compare_predict_helper <- function(x, ci_width){
  post_mean <- mean(x)
  post_sd <- stats::sd(x)
  low <- (1 - ci_width) / 2
  up  <-  1 - low
  interval <-  t(stats::quantile(x, c(low, up)))
  summ <-  round(cbind.data.frame(post_mean = post_mean,
                                  post_sd = post_sd,
                                  interval), 3)
}

# delta give prior_sd
delta_solve = function(x){
  (x^2)^-1 - 1
}

# fisher z to r
z2r <- function (z) {
  (exp(2 * z) - 1)/(1 + exp(2 * z))
}

# lower triangle of matrix
get_lower_tri<-function(cormat){
  cormat[upper.tri(cormat)] <- NA
  return(cormat)
}

# analytic solution
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
  list(inv_mu = inv_mu,
       inv_var = inv_var,
       partial = partials)
}

# summarize coefficients
beta_summary <- function(x, node, ci_width, samples){

  # convert inverse to beta
  x <- inverse_2_beta(x, samples = samples)

  # stop if not the correct class
  if(class(x) != "inverse_2_beta"){
    stop("class must be inverse_2_beta")
  }

  # check ci_width is allowed
  if(ci_width >= 1 | ci_width <= 0){
    stop("ci_width must be between 0 and 1")
  }
  returned_object <- lapply(node, function(y) summary_beta_helper(node =  y,
                                                                  x = x,
                                                                  ci_width))
  class(returned_object) <- "beta_summary"
  returned_object$betas <- x$betas[[node]]
  returned_object$sigma <- x$sigma[[node]]
  returned_object$call <- match.call()
  returned_object
}

rope_helper <- function(x, rope){
  mean(- rope < x & x < rope )

}

ci_helper <- function(x, ci_width){
  low <- (1 - ci_width) / 2
  up  <-  1 - low
  interval <-  stats::quantile(x, c(low, up))
  as.numeric(ifelse(interval[1] < 0 & interval[2] > 0, 0, 1))
}

Mo_risk_help_node <- function(x, post,  n1, n2, p){

  inv_mat <- post[,,x]

  Y_rep1 <- mvnfast::rmvn(n = n1,  mu = rep(0, p), sigma = stats::cov2cor(solve(inv_mat)))

  Y_rep2 <-  mvnfast::rmvn(n = n2, mu = rep(0, p), sigma = stats::cov2cor(solve(inv_mat)))

  jsd_node <- unlist(lapply(1:ncol(Y_rep1), function(z) node_jsd_help(z, Y_rep1, Y_rep2)))

  jsd_node



}


node_jsd_help <- function(x, Y_rep1, Y_rep2){

  Y_rep1 <- scale(Y_rep1)

  Y_rep2 <- scale(Y_rep2)

  pred1 <- Y_rep1[,-x]  %*% beta_helper(Y_rep1, x)

  pred2 <- Y_rep2[,-x]  %*%  beta_helper(Y_rep2, x)

  jsd_node <- (kl_func(stats::var(pred1), stats::var(pred2)) +
               kl_func(stats::var(pred2), stats::var(pred1))) * .5

  jsd_node

}


beta_helper <- function(x, which_one){
  y <- x[,which_one]
  X <- x[,-which_one]

  fit <- lm(y ~ 0 + X)
  coef(fit)

}

inverse_2_beta <- function(fit, samples = 500){
  fit <- fit[1:6]

  # check number of samples
  if(samples > fit$iter){
    stop("Samples used to compute R2 cannot be greater than the number used for fitting the model")

  }

  # get posterior estimate for precision matrix
  inv <-  fit$posterior_samples[,  grep("cov_inv", colnames(fit$posterior_samples))]

  # seperate estimates by row
  node_wise_elements <- lapply(split(colnames(inv), 1:fit$p), function(x) inv[,x])

  # convert off-diagonals to betas
  betas <- lapply(1:fit$p, function(x) node_wise_elements[[x]][-x] *  as.matrix((1/ node_wise_elements[[x]][x]) ) * -1)

  # convert diagonals to residual variance
  sigmas <- lapply(1:fit$p, function(x) as.matrix((1/ node_wise_elements[[x]][x])))


  # betas: select number of posterior samples
  betas <- lapply(betas, function(x)  x[1:samples,])

  # sigma: select number of posterior samples (sd scale)
  sigmas <- lapply(sigmas, function(x) sqrt(x[1:samples,]))


  # returned object
  returned_object <- list(betas = betas,
                          sigmas = sigmas,
                          p = fit$p,
                          data = fit$dat)

  # assign class
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
  post_sd <- apply(x$betas[[node]], 2, stats::sd)

  # lower and upper of posterior
  beta_ci <- t(apply(x$betas[[node]], 2, quantile, probs = c(low, up)))

  # sd of the outcome
  sd_y <- stats::sd(x$data[,node])

  # sd of the predictors
  sd_x <- apply(x$data[,-node], 2, stats::sd)

  # standardized (std) beta
  beta_std_temp <- x$betas[[node]] * (sd_x / sd_y)

  # beta std posterior mean
  beta_std <- apply(beta_std_temp, 2, mean)

  # beta std posterior sd
  beta_std_post_sd <- apply(beta_std_temp, 2, stats::sd)

  # lower and upper of posterior
  beta_std_ci <- t(apply(beta_std_temp, 2,  quantile, probs = c(low, up)))

  returned_object <- round(cbind(Node =  row_names[-node],
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
  returned_object$call <- match.call()
  returned_object
}


R2_helper <- function(ypred, y, ci_width) {
  low <- (1 - ci_width) / 2
  up  <-  1 - low
  e <- -1 * sweep(ypred, 2, y)
  var_ypred <- apply(ypred, 1, stats::var)
  var_e <- apply(e, 1, stats::var)
  r2 <- unlist(var_ypred / (var_ypred + var_e))
  ci <- quantile(r2, prob = c(low, up) )
  mu_r2 <- mean(r2)
  sd_r2 <- stats::sd(r2)
  summary_r2 <- c(post_mean = mu_r2, post_sd = sd_r2, ci)
  list(summary_r2 = summary_r2, R2 = r2)
}


MSE_helper <- function(ypred, y, ci_width){
  low <- (1 - ci_width) / 2
  up  <-  1 - low
  mse <- apply(ypred, MARGIN = 1, function(x){mean((x - y)^2)})
  ci <- quantile(mse, prob = c(low, up) )
  mu_mse <- mean(mse)
  sd_mse <- stats::sd(mse)
  summary_mse <- c(post_mean = mu_mse, post_sd = sd_mse, ci)
  list(summary_mse = summary_mse, MSE = mse)
}


name_helper <-  function(x){

  x <-  gsub("[A-z].*,", replacement = "", x)
  col_names <- gsub("[]]", "", x)
  col_names
}

error_helper <- function(ypred, y, ci_width, measure, sigmas =  NULL) {

  low <- (1 - ci_width) / 2

  up  <-  1 - low


  all_residual <- sweep(ypred, 2, y)
  if(measure == "mse"){
    out <- rowMeans(all_residual^2)
  }
  if(measure == "mae"){
    out <- rowMeans(abs(all_residual))
  }
  if(measure == "kl"){
    out <- kl_func(stats::var(y), sigmas^2)
    }
  ci <- quantile(out, prob = c(low, up) )
  mu_out <- mean(out)
  sd_out <- stats::sd(out)
  summary <- c(post_mean = mu_out, post_sd = sd_out, ci)
  list(summary = summary, error = out)
}

kl_func <- function(sigma_1, sigma_2){

  log(sqrt(sigma_2) / sqrt(sigma_1)) + (sigma_1 / (2 * sigma_2)) - .5

}

ppc_helper <- function(x, inv_g1, inv_cov, n, p){

  inv_mat <- matrix(0, p , p)
  inv_mat[,] <- as.numeric(inv_cov[x,])


  y_rep <- mvnfast::rmvn(n, mu = rep(0, p),sigma =  solve(inv_mat))

  S_rep <- t(y_rep) %*% y_rep

  theta_rep <- (n - 1) * solve(S_rep)

  KLD <- KL(Theta = inv_g1, hatTheta = theta_rep)

  JSD <- 0.5 * KL(Theta = inv_g1, hatTheta = theta_rep) + 0.5 * KL(hatTheta = theta_rep, Theta = inv_g1)

  QL <-  QL(Theta = inv_g1, hatTheta = theta_rep)

  FL <- sum((stats::cov2cor(inv_g1) *-1 - stats::cov2cor((theta_rep) * -1)^2))

  return <- list(KLD = KLD, JSD = JSD, QL = QL, FL = FL)
}

contrast_helper <- function(x){
  temp <- unlist(regmatches(x, gregexpr("[[:digit:]]+", x)))
  paste("Y", temp, sep = "_g", collapse = "_vs_")
}


axis_ticks_helper <- function(x){

  paste(stringr::str_sub(x, start = 1, end = 3),
        stringr::str_sub(x, start = 4, end = 5),
        stringr::str_sub(x, start = 6, end = 8))
}



KL = function(Theta,hatTheta){

  # Kuismin, M., & Sillanpaa, M. J. (2016). Use of Wishart prior and simple extensions for
  # sparse precision matrix estimation. PloS one, 11(2), e0148171.
  p = ncol(Theta)

  invTheta = solve(Theta,diag(1,p))

  kl  = 0.5 * (sum(diag(invTheta%*%hatTheta)) - log(det(invTheta%*%hatTheta)) - p)

  return(kl)

}

QL = function(Theta,hatTheta){


  # Kuismin, M., & Sillanpaa, M. J. (2016). Use of Wishart prior and simple extensions for
  # sparse precision matrix estimation. PloS one, 11(2), e0148171.

  p = ncol(Theta)
  I = diag(1,p)

  invTheta = solve(Theta,I)

  osa = sum(diag(invTheta%*%hatTheta - I))
  tulos = osa^2

  return(tulos)

}


unbiased_cov <- function(x){
  x <- scale(x)
  n <- nrow(x) - 1
  mle_cov <- n^-1 * t(x) %*% x
  stats::cov2cor(solve(mle_cov))
}




Mo_risk_help <- function(x, post, n1, n2, p){
  inv_mat <- post[,,x]
  Y_rep1 <-  mvnfast::rmvn(n = n1,  mu = rep(0, p), sigma = stats::cov2cor(solve(inv_mat)))
  Y_rep2 <-  mvnfast::rmvn(n = n2, mu = rep(0, p), sigma =  stats::cov2cor(solve(inv_mat)))

  jsd <- 0.5 *  KL(unbiased_cov(Y_rep1), unbiased_cov(Y_rep2)) +
         0.5 *  KL(unbiased_cov(Y_rep2), unbiased_cov(Y_rep1))

  jsd
}



Y_combine <- function(...){

  dat <- list(...)

  dat <- lapply(1:length(dat), function(x) na.omit(dat[[x]]))

  dat_info <- lapply(1:length(dat), function(x) {
    p <- ncol(dat[[x]])

    n <- nrow(dat[[x]])

    data.frame(p = p, n = n)
  })

  list(dat = dat, dat_info =  do.call(rbind, dat_info),
       pairwise = t(combn(1:length(dat), 2)))
}



approx_sd <- function(r, n, k){
  sqrt((1-r^2)/(n  - k - 2))
}

positive_helper <- function(pcor, post_sd, BF_null){
  dens_greater  <- (1  - pnorm(0, pcor, post_sd)) * 2
  BF_null * dens_greater
}

negative_helper <- function(pcor, post_sd, BF_null){
  dens_less  <- pnorm(0, pcor, post_sd) * 2
  BF_null * dens_less
}


exhaustive_helper <- function(BF_null, BF_positive, BF_negative){
  c(BF_null, BF_positive, BF_negative) /  sum(BF_null, BF_positive, BF_negative)
}

symmteric_mat <- function(x){
  x[lower.tri(x)] <- t(x)[lower.tri(x)]
  x
}

colnames_helper <- function(x, col_names){
  colnames(x) <- col_names
}

sampling_helper = function(X,  nu, delta,  n_samples){
  X <- as.matrix(X)
  # number of variables
  p <- ncol(X)
  # number of observations
  n <- nrow(X)
  # number of partial correlations
  pcors <- (p * (p - 1)) / 2
  # names for the partial correlations
  col_names <- numbers2words(1:p)

  mat_name <- matrix(unlist(lapply(col_names, function(x) paste(col_names,x, sep = ""))), p , p)
  mat_name_up <- mat_name[upper.tri(mat_name)]
  mat_name_low <- mat_name[lower.tri(mat_name)]

  # center the data
  Xhat <- X - rep(1,n)%*%t(apply(X,2,mean))

  # scatter matrix
  S <- t(Xhat)%*%Xhat

  # storage
  pcor_store_up <- pcor_store_low <- prior_store_up <- prior_store_low <- matrix(NA, nrow = n_samples, ncol = pcors)
  inv_cov_store <-  array(NA, c(p, p, n_samples))

  # initial values
  Psi <- b_inv <- diag(p)

  for(i in 1:n_samples){
    # draw from posterior
    post <- post_helper(S = S, n = n, nu = nu, p = p, delta = delta, Psi = Psi, b_inv = b_inv * 10000)
    # store partials
    pcor_store_up[i,] <- post$pcors_post_up
    pcor_store_low[i,] <- post$pcors_post_low
    # store the inverse
    inv_cov_store[,,i] <- post$sigma_inv
    # draw from prior and store
    prior_samps <- prior_helper(nu = nu, delta = delta, p = p)
    prior_store_up[i,] <- prior_samps$pcors_prior_up
    prior_store_low[i, ] <- prior_samps$pcors_prior_up
    # Psi
    Psi <- post$Psi
  }

  # transform posterior samples
  fisher_z_post_up <- apply(pcor_store_up, 2, fisher_z)
  fisher_z_post_low <- apply(pcor_store_low, 2, fisher_z)

  fisher_z_prior_up <- apply(prior_store_up, 2, fisher_z)
  fisher_z_prior_low <- apply(prior_store_low, 2, fisher_z)


  colnames(fisher_z_prior_up) <- mat_name_up
  colnames(fisher_z_post_up) <- mat_name_up
  colnames(pcor_store_up) <- mat_name_up


  colnames(fisher_z_prior_low) <- mat_name_low
  colnames(fisher_z_post_low) <- mat_name_low
  colnames(pcor_store_low) <- mat_name_low



  # returned list
  list(fisher_z_post = cbind(fisher_z_post_up, fisher_z_post_low),
       pcor_post = cbind(pcor_store_up, pcor_store_low),
       inv_cov_post = inv_cov_store,
       pcor_prior = cbind(prior_store_up, prior_store_low),
       fisher_z_prior =cbind(fisher_z_prior_up, fisher_z_prior_low))
}

prior_helper <- function(nu, p, delta){
  # sample from inverse Wishart
  inv_wish_prior <- solve(rWishart(1, df =  delta + p - 1, diag(p) * 1000)[,,1], tol = 1e-20)

  # sample from Wishart
  sigma_inv_prior <- rWishart(1, df = nu, inv_wish_prior)[,,1]

  # partical correlation matrix
  pcor_mat_prior <- - diag(1/sqrt(diag(sigma_inv_prior)))%*%sigma_inv_prior%*%diag(1/sqrt(diag(sigma_inv_prior)))
  pcors_prior_up <- pcor_mat_prior[upper.tri(pcor_mat_prior)]
  pcors_prior_low <- pcor_mat_prior[lower.tri(pcor_mat_prior)]

  list(pcors_prior_up = pcors_prior_up, pcors_prior_low = pcors_prior_low)
}



post_helper <- function(S, n, nu, p, delta, Psi, b_inv){
  # precision matrix
  sigma_inv <- rWishart(1, delta + n - 1, solve(Psi+S,  tol =  1e-20))[,,1]

  # Psi
  Psi <- rWishart(1, nu + delta + p - 2, solve(sigma_inv + b_inv,  tol  = 1e-20))[,,1]

  # partial correlation matrix
  pcor_mat <- - diag(1/sqrt(diag(sigma_inv)))%*%sigma_inv%*%diag(1/sqrt(diag(sigma_inv)))
  pcors_post_up = pcor_mat[upper.tri(pcor_mat)]
  pcors_post_low = pcor_mat[lower.tri(pcor_mat)]

  # returned list
  list(pcors_post_up = pcors_post_up, pcors_post_low = pcors_post_low, sigma_inv = sigma_inv, Psi = Psi)
}

sampling <- function(X, nu, delta, n_samples = 20000, cores = 4){
  # register parallel
  cl <- parallel::makeCluster(cores)
  doSNOW::registerDoSNOW(cl)

  # samples for each "chain"
  samps <- rep(round(n_samples / cores), cores)
  chains <- cores

  # global variable
  i <- 1

  # sample from priors and posteriors
  samples <- foreach::foreach(i = 1:chains,
                              .export = c("fisher_z", "sampling_helper",
                                          "numbers2words", "prior_helper", "post_helper")) %dopar% {
                                                sampling_helper(X = X, nu = nu,
                                                                delta = delta,
                                                                n_samples = samps[i])
                                            }
  # stop cluster
  parallel::stopCluster(cl)

  return(samples)
}

fisher_z <- function(rho){
  .5 * log(( 1 + rho )/ ( 1 - rho ))
}

sd_helper <- function(post_samples, prior_at_zero){
  prior_at_zero /  dnorm(0, mean(post_samples), stats::sd(post_samples))
}

pcor_name_helper <- function(x){
  keep_vars <-  unlist(strsplit(gsub("[^[:alnum:] ]", "", x), " +"))
  keep_vars
}

framer <- function(x){
  pos_comparisons <- unlist(gregexpr("[<>=]", x))
  leftside <- rep(NA, length(pos_comparisons) + 1)
  rightside <- rep(NA, length(pos_comparisons) + 1)
  pos1 <- c(-1, pos_comparisons)
  pos2 <- c(pos_comparisons, nchar(x) + 1)
  for(i in seq_along(pos1)){
    leftside[i] <- substring(x, pos1[i] + 1, pos1[i+1] - 1)
    rightside[i] <- substring(x, pos2[i] + 1, pos2[i+1] - 1)
  }
  leftside <- leftside[-length(leftside)]
  rightside <- rightside[-length(rightside)]
  comparisons <- substring(x, pos_comparisons, pos_comparisons)
  data.frame(left = leftside,
             comp = comparisons,
             right = rightside,
             stringsAsFactors = FALSE)

}



create_matrices <- function(framed, varnames){
  k <- length(varnames)
  if(any(grepl(",", framed$left)) || any(grepl(",", framed$right))){
    if(nrow(framed) > 1){
      for(r in 1:(nrow(framed)-1)){
        if(all.equal(framed$right[r], framed$left[r+1])){
          if(substring(framed$right[r], 1, 1) == "(") {
            framed$right[r] <- sub("),.+", ")", framed$right[r])
            framed$left[r+1] <- sub(".+),", "", framed$left[r +1])
          } else{
            framed$right[r] <- sub(",.+", "", framed$right[r])
            framed$left[r+1] <- sub("[^,]+,", "", framed$left[r+1])
          }
        }
      }
    }

    commas_left <- framed$left[grep(",", framed$left)]
    commas_right <- framed$right[grep(",", framed$right)]
    if(isTRUE(any(!grepl("\\(.+)", commas_left))) || isTRUE(any(!grepl("\\(.+)", commas_right))) ||
       isTRUE(any(grepl(").+", commas_left))) || isTRUE(any(grepl(").+", commas_right))) ||
       isTRUE(any(grepl(".+\\(", commas_left))) || isTRUE(any(grepl(".+\\(", commas_right)))) {
      stop("Incorrect hypothesis syntax or extra character, check specification")
    }

    framed$left <- gsub("[()]", "", framed$left)
    framed$right <- gsub("[()]", "", framed$right)
    commas <- unique(c(grep(",", framed$left), grep(",", framed$right)))

    if(length(commas) > 0){
      multiples <- vector("list", length = length(commas))

      for(r in seq_along(commas)){
        several <- framed[commas,][r, ]

        if(several$comp == "="){

          several <- c(several$left, several$right)
          separate <- unlist(strsplit(several, split = ","))
          if(any(grepl("^$", several))) stop("Misplaced comma in hypothesis")
          converted_equality <- paste(separate, collapse = "=")
          multiples[[r]] <- framer(converted_equality)

        } else{
          leftvars <- unlist(strsplit(several$left, split = ","))
          rightvars <- unlist(strsplit(several$right, split = ","))
          if(any(grepl("^$", leftvars)) || any(grepl("^$", rightvars))) stop("Misplaced comma in hypothesis")

          left <- rep(leftvars, length.out = length(rightvars)*length(leftvars))
          right <- rep(rightvars, each = length(leftvars))
          comp <- rep(several$comp, length(left))

          multiples[[r]] <- data.frame(left = left, comp = comp, right = right, stringsAsFactors = FALSE)
        }
      }

      framed <- framed[-commas,]
      multiples <- do.call(rbind, multiples)
      framed <- rbind(multiples, framed)
    }
  }

  equality <- framed[framed$comp == "=",]
  inequality <- framed[!framed$comp == "=",]

  #****Equality part string-to-matrix
  if(nrow(equality) == 0) {
    R_e <- r_e <- NULL
  } else{
    outcomes <- suppressWarnings(apply(equality[, -2], 2, as.numeric))
    outcomes <- matrix(outcomes, ncol = 2, byrow = FALSE)
    if(any(rowSums(is.na(outcomes)) == 0)) stop("Value compared with value rather than variable, e.g., '2 = 2', check hypotheses")
    rows <- which(rowSums(is.na(outcomes)) < 2)
    specified <- t(outcomes[rows,])
    specified <- specified[!is.na(specified)]
    r_e <- ifelse(rowSums(is.na(outcomes)) == 2, 0, specified)
    r_e <- matrix(r_e)

    var_locations <- apply(equality[, -2], 2, function(x) ifelse(x %in% varnames, match(x, varnames), 0))
    var_locations <- matrix(var_locations, ncol = 2)

    R_e <- matrix(rep(0, nrow(equality)*length(varnames)), ncol = length(varnames))

    for(i in seq_along(r_e)){
      if(!all(var_locations[i, ] > 0)){
        R_e[i, var_locations[i,]] <- 1
      } else{
        R_e[i, var_locations[i,]] <- c(1, -1)
      }
    }
  }


  #****Inequality part string-to-matrix
  if(nrow(inequality) == 0) {
    R_i <- r_i <- NULL
  } else{
    outcomes <- suppressWarnings(apply(inequality[, -2], 2, as.numeric))
    outcomes <- matrix(outcomes, ncol = 2, byrow = FALSE)
    if(any(rowSums(is.na(outcomes)) == 0)) stop("Value compared with value rather than variable, e.g., '2 > 2', check hypotheses")
    cols <- which(rowSums(is.na(outcomes)) < 2)
    specified <- t(outcomes[cols,])
    specified <- specified[!is.na(specified)]
    r_i <- ifelse(rowSums(is.na(outcomes)) == 2, 0, specified)
    r_i <- matrix(r_i)

    leq <- which(inequality$comp == "<")
    var_locations <- apply(inequality[, -2], 2, function(x) ifelse(x %in% varnames, match(x, varnames), 0))
    var_locations <- matrix(var_locations, ncol = 2)

    R_i <- matrix(rep(0, nrow(inequality)*length(varnames)), ncol = length(varnames))

    for(i in seq_along(r_i)){
      if(!all(var_locations[i, ] > 0)){

        if(var_locations[i, 1] == 0){
          if(i %in% leq){
            value <-  1
          } else{
            r_i[i] <- r_i[i]*-1
            value <- -1
          }
        } else{
          if(i %in% leq){
            r_i[i] <- r_i[i]*-1
            value <-  -1
          } else{
            value <- 1
          }
        }

        R_i[i, var_locations[i,]] <- value

      } else{
        value <- if(i %in% leq) c(-1, 1) else c(1, -1)
        R_i[i, var_locations[i,]] <- value
      }
    }
  }

  #3)check comparisons----------------
  if(is.null(R_i)){
    comparisons <- "only equality"
  } else if(is.null(R_e)){
    comparisons <- "only inequality"
  } else{
    comparisons <- "both comparisons"
  }

  #set prior mean
  R_ei <- rbind(R_e,R_i)
  r_ei <- rbind(r_e,r_i)
  Rr_ei <- cbind(R_ei,r_ei)
  beta_zero <- MASS::ginv(R_ei)%*%r_ei

  if(nrow(Rr_ei) > 1){
    rref_ei <- pracma::rref(Rr_ei)
    nonzero <- rref_ei[,k+1]!=0
    if(max(nonzero)>0){
      row1 <- max(which(nonzero==T))
      if(sum(abs(rref_ei[row1,1:k]))==0){
        stop("Default prior mean cannot be constructed from constraints.")
      }
    }
  }


  list(R_i = R_i,
       r_i = r_i,
       R_e = R_e,
       r_e = r_e,
       R_ei = R_ei,
       Rr_ei = Rr_ei,
       r_ei = r_ei,
       beta_zero = beta_zero,
       comparisons = comparisons)

}

word2num <- function(word){
  wsplit <- strsplit(tolower(word)," ")[[1]]
  one_digits <- list(zero=0, one=1, two=2, three=3, four=4, five=5,
                     six=6, seven=7, eight=8, nine=9)
  teens <- list(eleven=11, twelve=12, thirteen=13, fourteen=14, fifteen=15,
                sixteen=16, seventeen=17, eighteen=18, nineteen=19)
  ten_digits <- list(ten=10, twenty=20, thirty=30, forty=40, fifty=50,
                     sixty=60, seventy=70, eighty=80, ninety=90)
  doubles <- c(teens,ten_digits)
  out <- 0
  i <- 1
  while(i <= length(wsplit)){
    j <- 1
    if(i==1 && wsplit[i]=="hundred")
      temp <- 100
    else if(i==1 && wsplit[i]=="thousand")
      temp <- 1000
    else if(wsplit[i] %in% names(one_digits))
      temp <- as.numeric(one_digits[wsplit[i]])
    else if(wsplit[i] %in% names(teens))
      temp <- as.numeric(teens[wsplit[i]])
    else if(wsplit[i] %in% names(ten_digits))
      temp <- (as.numeric(ten_digits[wsplit[i]]))
    if(i < length(wsplit) && wsplit[i+1]=="hundred"){
      if(i>1 && wsplit[i-1] %in% c("hundred","thousand"))
        out <- out + 100*temp
      else
        out <- 100*(out + temp)
      j <- 2
    }
    else if(i < length(wsplit) && wsplit[i+1]=="thousand"){
      if(i>1 && wsplit[i-1] %in% c("hundred","thousand"))
        out <- out + 1000*temp
      else
        out <- 1000*(out + temp)
      j <- 2
    }
    else if(i < length(wsplit) && wsplit[i+1] %in% names(doubles)){
      temp <- temp*100
      out <- out + temp
    }
    else{
      out <- out + temp
    }
    i <- i + j
  }
  return(list(word,out))
}

numbers2words <- function(x){
  ## Function by John Fox found here:
  ## http://tolstoy.newcastle.edu.au/R/help/05/04/2715.html
  ## Tweaks by AJH to add commas and "and"
  helper <- function(x){

    digits <- rev(strsplit(as.character(x), "")[[1]])
    nDigits <- length(digits)
    if (nDigits == 1) as.vector(ones[digits])
    else if (nDigits == 2)
      if (x <= 19) as.vector(teens[digits[1]])
    else trim(paste(tens[digits[2]],
                    Recall(as.numeric(digits[1]))))
    else if (nDigits == 3) trim(paste(ones[digits[3]], "hundred and",
                                      Recall(makeNumber(digits[2:1]))))
    else {
      nSuffix <- ((nDigits + 2) %/% 3) - 1
      if (nSuffix > length(suffixes)) stop(paste(x, "is too large!"))
      trim(paste(Recall(makeNumber(digits[
        nDigits:(3*nSuffix + 1)])),
        suffixes[nSuffix],"," ,
        Recall(makeNumber(digits[(3*nSuffix):1]))))
    }
  }
  trim <- function(text){
    #Tidy leading/trailing whitespace, space before comma
    text=gsub("^\ ", "", gsub("\ *$", "", gsub("\ ,",",",text)))
    #Clear any trailing " and"
    text=gsub(" and$","",text)
    #Clear any trailing comma
    gsub("\ *,$","",text)
  }
  makeNumber <- function(...) as.numeric(paste(..., collapse=""))
  #Disable scientific notation
  opts <- options(scipen=100)
  on.exit(options(opts))
  ones <- c("", "one", "two", "three", "four", "five", "six", "seven",
            "eight", "nine")
  names(ones) <- 0:9
  teens <- c("ten", "eleven", "twelve", "thirteen", "fourteen", "fifteen",
             "sixteen", " seventeen", "eighteen", "nineteen")
  names(teens) <- 0:9
  tens <- c("twenty", "thirty", "forty", "fifty", "sixty", "seventy", "eighty",
            "ninety")
  names(tens) <- 2:9
  x <- round(x)
  suffixes <- c("thousand", "million", "billion", "trillion")
  if (length(x) > 1) return(trim(sapply(x, helper)))
  helper(x)
}

samps_inv_helper <- function(x, p){
  inv <- paste("cov_inv", paste(paste("[", paste( 1:p, x, sep = ","), sep = ""), "]", sep = ""), sep = "")
  inv
}

samps_pcor_helper <- function(x, p){

  pcors <- paste("pcors", paste(paste("[", paste( 1:p, x, sep = ","), sep = ""), "]", sep = ""), sep = "")
  pcors
}

hyp_converter <- function(x){

  hyp_converted <- x

  extract_numbers <- unlist(stringr::str_extract_all(hyp_converted, "\\d+"))

  extract_numbers <- extract_numbers[unlist(extract_numbers) != 0 ]
  words <- NA
  for(i in 1:length(extract_numbers)){

    temp <- noquote(extract_numbers[i])
    words[i] <- numbers2words(as.numeric(temp))
    hyp_converted <- sub(temp, numbers2words(as.numeric(temp)), hyp_converted)


  }

  hyp_converted <- stringr::str_remove_all(hyp_converted, "--")

  list(hyp_converted = hyp_converted, words = words)
}

performance <- function(Estimate, True){

  True <- as.matrix(True)
  Estimate <- as.matrix(Estimate)

  # True Negative
  TN <- ifelse(True[upper.tri(True)] == 0 & Estimate[upper.tri(Estimate)] == 0, 1, 0); TN <- sum(TN)
  # False Positive
  FP <- ifelse(True[upper.tri(True)] == 0 & Estimate[upper.tri(Estimate)] != 0, 1, 0); FP <- sum(FP)
  # True Positive
  TP <- ifelse(True[upper.tri(True)] != 0 & Estimate[upper.tri(Estimate)] != 0, 1, 0); TP <- sum(TP)
  # False Negatives
  FN <- ifelse(True[upper.tri(True)] != 0 & Estimate[upper.tri(Estimate)] == 0, 1, 0); FN <- sum(FN)

  Specificity <- TN/(TN + FP)
  Sensitivity <- TP/(TP + FN)
  Precision <- TP/(TP + FP)

  Recall <- TP / (TP + FN)

  F1_score <- 2 * ((Precision * Recall) / (Precision + Recall))

  MCC <- (TP*TN - FP*FN)/sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN))

  results <- c(Specificity, Sensitivity, Precision, Recall,  F1_score, MCC)
  results_name <- c("Specificity", "Sensitivity", "Precision", "Recall",  "F1_score", "MCC")
  results <- cbind.data.frame(measure = results_name, score = results)
  list(results = results)


}
