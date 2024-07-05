#' @importFrom stats coef cov2cor var dnorm lm
#' na.omit pnorm quantile rWishart runif dnorm rnorm
#' sd qnorm residuals fitted density weighted.mean
#' @importFrom utils combn
#' @importFrom graphics plot
#' @importFrom Rdpack reprompt
#' @importFrom MASS ginv
#' @import ggplot2
#' @importFrom stats model.matrix terms cor


tstat <- function(r, n, k){
  r * sqrt(((n - 2 - k) / (1-r^2)))
}

hyp_converter <- function (x) {

  hyp_converted <- x

  extract_numbers <-strsplit(  gsub("[^\\d]+", "",
                                    hyp_converted, perl = TRUE),
                               split = "" )[[1]]

  extract_numbers <- extract_numbers[unlist(extract_numbers) >= 1]

  words <- NA

  for (i in 1:length(extract_numbers)) {
    temp <- noquote(extract_numbers[i])
    words[i] <- numbers2words(as.numeric(temp))
    hyp_converted <- sub(temp, numbers2words(as.numeric(temp)),
                         hyp_converted)
  }

  list(hyp_converted = hyp_converted, words = words)
}


remove_predictors_helper <- function(Y_groups, formula){

  # number of groups
  groups <- length(Y_groups)

  # Y groups
  Y_groups <- lapply(seq_len(groups), function(x)  na.omit(Y_groups[[x]]) )

  model_matrices <- lapply(seq_len(groups) , function(x) {
    stats::model.matrix(formula, Y_groups[[x]])
  })

  # model matrix terms
  mm_terms <- attr(terms(formula), "term.labels")

  if(length(mm_terms) == 0){

    Y_groups <- Y_groups

  } else {

    Y_groups <- lapply(seq_len(groups), function(x){

      # check for factors
      factor_pred <- which(paste0("as.factor(", colnames(Y_groups[[x]]), ")") %in% mm_terms)

      # check for scaled
      scale_pred <- which(paste0("scale(", colnames(Y_groups[[x]]), ")") %in% mm_terms)

      # check for non factors
      cont_pred <- which(colnames(Y_groups[[x]]) %in% mm_terms)

      # remove predictors
      Y_groups[[x]][,-c(factor_pred, cont_pred, scale_pred)]

    })
  }

  list(Y_groups = Y_groups, model_matrices = model_matrices)

  }



binary_latent_helper <- function(x){
  # obervations
  n <- nrow(x)

  # variables
  p <- ncol(x)

  # thresholds
  thresholds <- c(-Inf, 0, Inf)

  # latent data
  latent_data <- sapply(1:p, function(z) qnorm(runif(n, min = pnorm(thresholds[x[,z]],
                                                                    mean = 0,
                                                                    sd = 1),
                                                     max = pnorm(thresholds[x[,z]+1],
                                                                 mean = 0,
                                                                 sd = 1)),
                                               mean = 0, sd = 1))
  # latent data (sd = 1)
  latent_data <- scale(latent_data)

  # latent data
  latent_data
}


ordinal_latent_helper <- function(x, thresholds){

  # observations
  n <- nrow(x)

  # variables
  p <- ncol(x)

  # mean of thresholds
  thresholds <- t(sapply(1:p, function(x) colMeans(thresholds[,,x])))

  latent_data <- sapply(1:p, function(z)  qnorm(runif(n, min = pnorm(thresholds[z, x[,z]],
                                                                     mean = 0,
                                                                     sd = 1),
                                                         max = pnorm(thresholds[z, x[,z]+1],
                                                                     mean = 0,
                                                                     sd = 1)),
                                                mean = 0, sd = 1))
  # latent data
  latent_data <- scale(latent_data)

  # latent data
  latent_data

  }

rank_helper <- function(Y){

  # adapted from hoff (2008). See documentation.
  p <- ncol(Y)

  levels  <- apply(Y, 2, function(x) {match(x, sort(unique(x)))})

  K <-  apply(levels, 2, max, na.rm = TRUE)

  ranks <- apply(Y, 2, rank, ties.method = "max", na.last = "keep")

  n_complete <- apply(!is.na(ranks), 2, sum)

  U <- t(t(ranks)/(n_complete + 1))

  Z <- qnorm(U)

  S <- cov(Z)

  list(K = K,
       levels = levels,
       Sigma_start = S,
       z0_start = Z)

}



group_hyp_helper <- function(hypothesis, x){

  hyp <- gsub(hyp_converter(convert_hyps(hypothesis = hypothesis, cbind(x)))$hyp_converted,
              pattern = "_", replacement = "")
  hyp
}


# convert hypothesis to words
convert_hyps <- function(hypothesis, Y){

  p <- ncol(Y)

  col_names <- numbers2words(1:p)

  mat_name_num <- sapply(1:p, function(x) paste0(1:p, "--", x, sep = ""))

  mat_name_n2w <- sapply(col_names, function(x) paste(col_names,x, sep = ""))

  mat_name_actual <- sapply(colnames(Y), function(x) paste0(colnames(Y), "--", x, sep = ""))

  n_off_diag <- p *( p - 1)*0.5

  where <- lapply(1:n_off_diag, function(x) grep(mat_name_actual[upper.tri(diag(p))][x], hypothesis ))

  where <- which(where == 1)

  if(any(where)){

    for (i in seq_along(where)) {
      hypothesis <-  gsub(mat_name_actual[upper.tri(diag(p))][where[i]],
                     replacement =  mat_name_n2w[upper.tri(diag(p))][where[i]], x = hypothesis )
    }

  } else if (!any(where)){

    where <- lapply(1:n_off_diag, function(x) grep(mat_name_num[upper.tri(diag(p))][x], hypothesis ))
    where <- which(where == 1)


    for (i in seq_along(where)) {
      hypothesis <-  gsub(mat_name_num[upper.tri(diag(p))][where[i]],
                          replacement =  mat_name_n2w[upper.tri(diag(p))][where[i]], x = hypothesis )
    }

  } else {

    stop("error in edge specification\n
             Hints:\n
             1) check that the first is smaller than the second (denoting the upper-triangular elements)-e.g., 1--2 and not 2--1\n
             2) alternatively, the variable might not exist--e.g., p = 10 but 1--11 (11 > 10)")
  }

  hypothesis

}



prior_helper_2 <- function(p, delta , epsilon){

  k <- p

  nu <- 1/epsilon

  parcorMat <- corMat <- ThetaMat <- SigmaMat <- array(0, dim=c(1e4 , k , k))

  for(s in 1:dim(parcorMat)[1]){

    Psi <- rWishart(1,df = nu, Sigma = 1/nu*diag(k))[,,1]

    Sigma <- rWishart(1,df=k-1+delta, Sigma=solve(Psi))[,,1]

    Theta <- solve(Sigma)

    ThetaMat[s,,] <- Theta
    SigmaMat[s,,] <- Sigma

    parcorMat[s,,] <- -diag(1/sqrt(diag(Theta)))%*%Theta%*%diag(1/sqrt(diag(Theta)))

  }
  parcorMat
}

print_summary_explore <- function(x,...){
  summary(x$dat_results, summarize = TRUE)
}


print_post_pred <- function(x,...){
  if(length(x) != 2){
    class(x) <- ""
    x <- round(x, 3)
    print(x)
  } else {
    cat("'summary = FALSE' not printed. See object contents")
  }
}

print_map <- function(x,...){
  cat("BGGM: Bayesian Gaussian Graphical Models \n")
  cat("--- \n")
  cat("Method: Maximum A Posteriori \n")
  cat("--- \n")
  print(x$pcor)
}


print_fitted <- function(x,...){
  if(length(x) == 2){

    cat("'summary = FALSE' not printed. See object contents")

  } else {

    class(x) <- ""
    x <- round(x, 3)
    print(x)
  }
}

print_predict <- function(x,...){
  if(length(x) == 2){

    cat("'summary = FALSE' not printed. See object contents")

  } else {

    class(x) <- ""
    x <- round(x, 3)
    print(x)
  }
}

post_prob <- function(data){
  p1 <- sum(data>0)/length(data)
  p2 <- sum(data<0)/length(data)
  p<- 1 - min(p1,p2)
  return(p)
}



convert_colnames <- function(hyp, Y){
  names_temp <- unlist(strsplit( strsplit(hyp, " ")[[1]], "--"))
  names_temp <- paste(names_temp, collapse = " ")
  names_temp <- unique(strsplit(gsub("[^[:alnum:] ]", "", names_temp), " +")[[1]])

  if(any(names_temp == "0")){
  names_temp <- names_temp[-which(names_temp == "0" )]
  }
  if(!all(names_temp %in% colnames(Y))){
    stop("node names not found in the data")
  }
  for(i in 1:length(names_temp)){
    id <- which(names_temp[[i]]  == colnames(Y))
    hyp <- gsub(x = hyp, pattern = names_temp[[i]],  replacement = id)
  }
  hyp
}

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

## delta_solve takes prior_sd
## prior_sd is in refernce to the stretched sd of the beta distribution with Beta(delta/2, delta/2) * 2-1
## Its SD = sqrt((d/2)*(d/2) / (((d/2)+(d/2))^2*((d/2)+(d/2)+1))) * 2-1;
##     SD = 1/sqrt( d+1 ) - 1
## Solving for d (delta):
##      d = 1/SD^2 -1
## Problem: with SD greater than sqrt(1/2) approx 0.71, d goes below 1 and beta density accumulates at the margins
##          With SD greater than 1, d is negative
##          With SD = sqrt(1/3), d = 2, seems ok as default for model estimation
##          With SD = sqrt(1/4) = .5, d = 3, seems more reasonable for hypothesis testing
## Solution: limit user input to [0,sqrt(1/2)]

delta_solve = function(x){
  if(x <= 0 || x > sqrt(1/2) ) stop("Error: \nPrior_sd must be between 0 and sqrt(1/2) approx. 0.7, to ensure that delta is not less than 1.\nFor delta = 1, set prior_sd to sqrt(1/2)\nFor delta = 2, set prior_sd to sqrt(1/3).")
  1/x^2 - 1
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

  # variables
  p <- ncol(X)

  # centererd mat
  X <- scale(X, scale = T)

  # scatter matrix
  S <- t(X) %*% X

  # degrees of freedom
  df = p

  # map estimate
  inv_map <-  solve(S + diag(10^-5,  p)) * (n + df - p - 1)

  # posterior variane
  inv_var <-  (n + df + 1) * (solve(S + diag(0.1^5, p) )^2 + tcrossprod(diag(solve(S + diag(0.1^5, p)))))

  # inverse correlation matrix
  inv_cor <- diag( 1 / sqrt((diag(inv_map)))) %*% inv_map %*% diag( 1 / sqrt((diag(inv_map))) )

  pcor_mat <- -(inv_cor - diag(p))

  list(inv_map = inv_map,
       inv_var = inv_var,
       pcor_mat = pcor_mat)
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





KL = function(Theta,hatTheta){

  # Kuismin, M., & Sillanpaa, M. J. (2016). Use of Wishart prior and simple extensions for
  # sparse precision matrix estimation. PloS one, 11(2), e0148171.
  p = ncol(Theta)

  invTheta = solve(Theta,diag(1,p))

  kl  = 0.5 * (sum(diag(invTheta%*%hatTheta)) - log(det(invTheta%*%hatTheta)) - p)

  return(kl)

}

QL = function(Theta,hatTheta){
  # Kuismin, M., & Sillanpaa, M. J. (2016). Use of Wishart prior
  # and simple extensions for
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

symmetric_mat <- function(x){
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
  Psi <- b_inv <- sigma_inv <- diag(p)

  for(i in 1:n_samples){

    # draw from posterior
    post <- post_helper(S = S, n = n,
                        nu = nu, p = p,
                        delta = delta,
                        Psi = Psi,
                        sigma_inv = sigma_inv)
    # store partials
    pcor_store_up[i,] <- post$pcors_post_up
    pcor_store_low[i,] <- post$pcors_post_low
    # store the inverse
    inv_cov_store[,,i] <- post$sigma_inv
    # draw from prior and store
    prior_samps <- prior_helper(delta = delta, p = p)
    prior_store_up[i,] <- prior_samps$pcors_prior_up
    prior_store_low[i, ] <- prior_samps$pcors_prior_up
    # Psi
    Psi <- post$Psi
    sigma_inv <- post$sigma_inv
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

prior_helper <- function(p, delta){

  I_p <- diag(p)

  nu <- 1/0.001

  Psi <-  rWishart(1, df = nu, I_p * 0.001)[,,1]

  sigma_inv_prior <- solve(rWishart(1, df = p - 1 + delta, solve(Psi))[,,1])

  pcor_mat_prior <- - diag(1/sqrt(diag(sigma_inv_prior)))%*%sigma_inv_prior%*%diag(1/sqrt(diag(sigma_inv_prior)))

  pcors_prior_up <- pcor_mat_prior[upper.tri(pcor_mat_prior)]
  pcors_prior_low <- pcor_mat_prior[lower.tri(pcor_mat_prior)]

  list(pcors_prior_up = pcors_prior_up,
       pcors_prior_low = pcors_prior_low)
}



post_helper <- function(S, n, nu, p, delta, Psi, b_inv, sigma_inv){

  nu = 1 / 0.001
  delta = delta
  B <- diag(p) * 0.001
  nuMP <- delta + delta - 1
  deltaMP <- nu - p + 1

  BMP <- solve(B)

  BMPinv <- solve(BMP)
  # Psi
  Psi <- rWishart(1, nuMP + deltaMP + p - 1, solve(sigma_inv + BMPinv,  tol  = 1e-20))[,,1]

  # precision matrix
  sigma_inv <- rWishart(1, (deltaMP + p - 1) + (n - 1), solve(Psi+S,  tol =  1e-20))[,,1]



  # partial correlation matrix
  pcor_mat <- - diag(1/sqrt(diag(sigma_inv)))%*%sigma_inv%*%diag(1/sqrt(diag(sigma_inv)))
  pcors_post_up = pcor_mat[upper.tri(pcor_mat)]
  pcors_post_low = pcor_mat[lower.tri(pcor_mat)]

  # returned list
  list(pcors_post_up = pcors_post_up,
       pcors_post_low = pcors_post_low, sigma_inv = sigma_inv, Psi = Psi)
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



create_matrices <- function(framed, varnames) {

  k <- length(varnames)

  if (any(grepl(",", framed$left)) ||
      any(grepl(",", framed$right))) {
    if (nrow(framed) > 1) {
      for (r in 1:(nrow(framed) - 1)) {
        if (all.equal(framed$right[r], framed$left[r + 1])) {
          if (substring(framed$right[r], 1, 1) == "(") {
            framed$right[r] <- sub("),.+", ")", framed$right[r])
            framed$left[r + 1] <- sub(".+),", "", framed$left[r + 1])

          } else{
            framed$right[r] <- sub(",.+", "", framed$right[r])
            framed$left[r + 1] <-
              sub("[^,]+,", "", framed$left[r + 1])
          }
        }
      }
    }

    commas_left <- framed$left[grep(",", framed$left)]
    commas_right <- framed$right[grep(",", framed$right)]
    if (isTRUE(any(!grepl("\\(.+)", commas_left))) ||
        isTRUE(any(!grepl("\\(.+)", commas_right))) ||
        isTRUE(any(grepl(").+", commas_left))) ||
        isTRUE(any(grepl(").+", commas_right))) ||
        isTRUE(any(grepl(".+\\(", commas_left))) ||
        isTRUE(any(grepl(".+\\(", commas_right)))) {
      stop("Incorrect hypothesis syntax or extra character, check specification")
    }

    framed$left <- gsub("[()]", "", framed$left)
    framed$right <- gsub("[()]", "", framed$right)
    commas <-
      unique(c(grep(",", framed$left), grep(",", framed$right)))

    if (length(commas) > 0) {
      multiples <- vector("list", length = length(commas))

      for (r in seq_along(commas)) {
        several <- framed[commas, ][r,]

        if (several$comp == "=") {
          several <- c(several$left, several$right)
          separate <- unlist(strsplit(several, split = ","))
          if (any(grepl("^$", several)))
            stop("Misplaced comma in hypothesis")
          converted_equality <- paste(separate, collapse = "=")
          multiples[[r]] <- framer(converted_equality)

        } else{
          leftvars <- unlist(strsplit(several$left, split = ","))
          rightvars <- unlist(strsplit(several$right, split = ","))
          if (any(grepl("^$", leftvars)) ||
              any(grepl("^$", rightvars)))
            stop("Misplaced comma in hypothesis")

          left <-
            rep(leftvars, length.out = length(rightvars) * length(leftvars))
          right <- rep(rightvars, each = length(leftvars))
          comp <- rep(several$comp, length(left))

          multiples[[r]] <-
            data.frame(
              left = left,
              comp = comp,
              right = right,
              stringsAsFactors = FALSE
            )
        }
      }

      framed <- framed[-commas, ]
      multiples <- do.call(rbind, multiples)
      framed <- rbind(multiples, framed)
    }
  }

  equality <- framed[framed$comp == "=", ]
  inequality <- framed[!framed$comp == "=", ]

  #****Equality part string-to-matrix
  if (nrow(equality) == 0) {
    R_e <- r_e <- NULL
  } else{
    outcomes <- suppressWarnings(apply(equality[,-2], 2, as.numeric))
    outcomes <- matrix(outcomes, ncol = 2, byrow = FALSE)
    if (any(rowSums(is.na(outcomes)) == 0))
      stop("Value compared with value rather than variable, e.g., '2 = 2', check hypotheses")
    rows <- which(rowSums(is.na(outcomes)) < 2)
    specified <- t(outcomes[rows, ])
    specified <- specified[!is.na(specified)]
    r_e <- ifelse(rowSums(is.na(outcomes)) == 2, 0, specified)
    r_e <- matrix(r_e)

    var_locations <-
      apply(equality[,-2], 2, function(x)
        ifelse(x %in% varnames, match(x, varnames), 0))
    var_locations <- matrix(var_locations, ncol = 2)

    R_e <-
      matrix(rep(0, nrow(equality) * length(varnames)), ncol = length(varnames))

    for (i in seq_along(r_e)) {
      if (!all(var_locations[i,] > 0)) {
        R_e[i, var_locations[i, ]] <- 1
      } else{
        R_e[i, var_locations[i, ]] <- c(1,-1)
      }
    }
  }


  #****Inequality part string-to-matrix
  if (nrow(inequality) == 0) {
    R_i <- r_i <- NULL
  } else{
    outcomes <- suppressWarnings(apply(inequality[,-2], 2, as.numeric))
    outcomes <- matrix(outcomes, ncol = 2, byrow = FALSE)
    if (any(rowSums(is.na(outcomes)) == 0))
      stop("Value compared with value rather than variable, e.g., '2 > 2', check hypotheses")
    cols <- which(rowSums(is.na(outcomes)) < 2)
    specified <- t(outcomes[cols, ])
    specified <- specified[!is.na(specified)]
    r_i <- ifelse(rowSums(is.na(outcomes)) == 2, 0, specified)
    r_i <- matrix(r_i)

    leq <- which(inequality$comp == "<")
    var_locations <-
      apply(inequality[,-2], 2, function(x)
        ifelse(x %in% varnames, match(x, varnames), 0))
    var_locations <- matrix(var_locations, ncol = 2)

    R_i <-
      matrix(rep(0, nrow(inequality) * length(varnames)), ncol = length(varnames))

    for (i in seq_along(r_i)) {
      if (!all(var_locations[i,] > 0)) {
        if (var_locations[i, 1] == 0) {
          if (i %in% leq) {
            value <-  1
          } else{
            r_i[i] <- r_i[i] * -1
            value <- -1
          }
        } else{
          if (i %in% leq) {
            r_i[i] <- r_i[i] * -1
            value <-  -1
          } else{
            value <- 1
          }
        }

        R_i[i, var_locations[i, ]] <- value

      } else{
        value <- if (i %in% leq)
          c(-1, 1)
        else
          c(1,-1)
        R_i[i, var_locations[i, ]] <- value
      }
    }
  }

  #3)check comparisons----------------
  if (is.null(R_i)) {
    comparisons <- "only equality"
  } else if (is.null(R_e)) {
    comparisons <- "only inequality"
  } else{
    comparisons <- "both comparisons"
  }

  #set prior mean
  R_ei <- rbind(R_e, R_i)
  r_ei <- rbind(r_e, r_i)
  Rr_ei <- cbind(R_ei, r_ei)
  # beta_zero <- MASS::ginv(R_ei) %*% r_ei

  if (nrow(Rr_ei) > 1) {
    # rref_ei <- pracma::rref(Rr_ei)
    nonzero <- rref_ei[, k + 1] != 0
    if (max(nonzero) > 0) {
      row1 <- max(which(nonzero == T))
      if (sum(abs(rref_ei[row1, 1:k])) == 0) {
        stop("Default prior mean cannot be constructed from constraints.")
      }
    }
  }
  # beta_zero = beta_zero,
list(
    R_i = R_i,
    r_i = r_i,
    R_e = R_e,
    r_e = r_e,
    R_ei = R_ei,
    Rr_ei = Rr_ei,
    r_ei = r_ei,
    comparisons = comparisons
  )

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


# make names for inverse
samps_inv_helper <- function(x, p){
  inv <- paste("cov_inv", paste(paste("[", paste( 1:p, x, sep = ","),
                                      sep = ""), "]", sep = ""), sep = "")
  inv
}

# make names for partials
samps_pcor_helper <- function(x, p){

  pcors <- paste("pcors", paste(paste("[", paste( 1:p, x, sep = ","),
                                      sep = ""), "]", sep = ""), sep = "")
  pcors
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

  results <- c(Specificity,
               Sensitivity,
               Precision,
               Recall,
               F1_score,
               MCC)

  results_name <- c("Specificity",
                    "Sensitivity",
                    "Precision",
                    "Recall",
                    "F1_score",
                    "MCC")

  results <- cbind.data.frame(measure = results_name, score = results)
  list(results = results)


}

################
# taken and modified from the package orddata
# under the GPL-2 licence
rmvord_naiv <-  function(n, probs, Cors, empirical) {

  q = length(probs)
  categ_probs = 0
  cumul_probs = list(0)
  quant_probs = list(0)
  means = 0
  vars = 0

  var.wt = function(x, w) {
    m = weighted.mean(x = x, w = w)
    sum((x[1:length(x)] - m)^2 * w[1:length(x)])
  }

  for (i in 1:q) {
    categ_probs[i] = length(probs[[i]])
    cumul_probs[[i]] = cumsum(1:categ_probs[i]/10^12 + probs[[i]])
    cumul_probs[[i]][categ_probs[i]] = 1
    quant_probs[[i]] = qnorm(p = cumul_probs[[i]], mean = 0,
                             sd = 1)
  }

  retval = MASS::mvrnorm(n = n, mu = rep(0,q),
                         Sigma = Cors,
                         empirical = empirical)

  for (i in 1:q) {
    retval[, i] = cut(x = retval[, i], breaks = c(-1/0, quant_probs[[i]]),
                      right = FALSE)
  }
  retval
}



csws_labels <- ifelse(1:35 %in% c(7,10,16,24,29),
                      "Family Support",
                      ifelse(1:35 %in% c(3,12,20,25,35),
                             "Competition",
                             ifelse(1:35 %in% c(1,4,17,21,30),
                                    "Appearence",
                                    ifelse(1:35%in%c(2,8,18,26,31),
                                           "God's Love",
                                           ifelse(1:35 %in% c(13, 19, 22, 27,  33),
                                                  "Academic Competence",
                                                  ifelse(1:35 %in% c(5, 11, 14, 28, 34),
                                                         "Virtue", "Approval From Others"))))))

tas_labels <- ifelse(1:20 %in% c(1,3,6,7,9,13,14),
                     "Difficulty\nIdentifying Feelings",
                     ifelse(1:20 %in% c(2,4,11,12,17),
                            "Difficulty\nDescribing Feelings",
                            "Externally\nOriented Feelings"))

iri_labels <- ifelse(1:28 %in% c(3, 8, 11, 15, 21, 25, 28),
                     "Perspective Taking",
                     ifelse(1:28 %in% c(2, 4, 9, 14, 18, 20, 22),
                            "Empathic Concern",
                            ifelse(1:28 %in% c(1, 5, 7, 12, 16, 23, 26), "Fantasy",
                                   "Personal Distress")))

rsa_labels <- ifelse(1:33 %in% c(1, 4, 5, 32),
                     "Planned Future",
                     ifelse(1:33 %in% c(2, 11, 17, 25, 31, 33),
                            "Perception of Self",
                            ifelse(1:33 %in% c(3, 7, 13, 16, 24, 29),
                                   "Family Cohesion",
                                   ifelse(1:33  %in% c(6, 9, 10, 12, 15, 19, 27),
                                          "Social Resources",
                                          ifelse(1:33 %in% c(8, 14, 18, 21, 22, 26),
                                                 "Social Competence", "Structured Style")))))

globalVariables(c('Y1','Y2',
                  'X1', 'X2',
                  'contrast',
                  '..quantile..',
                  'value',
                  'node',
                  'BF',
                  'Edge',
                  'Estimate',
                  'selected',
                  'probability',
                  'cred.lb',
                  'sig',
                  'hyp',
                  'label',
                  'color',
                  'fit',
                  'post_mean',
                  'Error',
                  'density', 'Node',
                  'Post.mean',
                  'L1', 'lag', 'acf',
                  'iteration',
                  '.imp',
                  'estimate',
                  'rref_ei', 'explore',
                  'print_coef',
                  'print_confirm',
                  'print_estimate',
                  'print_explore',
                  'print_ggm_compare',
                  'print_ggm_compare_bf',
                  'print_ggm_compare_ppc',
                  'print_ggm_confirm',
                  'print_roll_your_own',
                  'print_select_explore',
                  'print_select_ggm_compare_estimate',
                  'print_summary_coef',
                  'print_summary_estimate',
                  'print_summary_ggm_compare_bf',
                  'print_summary_ggm_estimate_compare',
                  'print_summary_metric',
                  'print_summary_select_explore',
                  '..', 'ppc', 'rope', 'y',
                  'Relation', 'Pr.H1',
                  "Y_missing",
                  "na_indices",
                  'Y'))


gen_pcors <-
  function (p = 20,
            edge_prob = 0.3,
            lb = 0.05,
            ub = 0.3) {
    d <- -1
    trys <- 0
    while (d < 0) {
      trys <- trys + 1
      effects <- p * (p - 1) * 0.5
      mat <- matrix(1, p, p)
      prob_zero <- 1 - edge_prob
      pool <- c(rep(0, effects * prob_zero),
                runif(effects *
                        edge_prob, lb, ub))
      if (length(pool) != effects) {
        pool <- c(0, pool)
      }
      mat[upper.tri(mat)] <- sample(pool, size = effects)
      pcs <- symmetric_mat(mat)
      pcs <- -pcs
      diag(pcs) <- -diag(pcs)
      d <- det(pcs)
    }
    cors <- cov2cor(solve(pcs))
    inv <- solve(cors)
    pcors <- cov2cor(inv) * -1
    diag(pcors) <- 1
    adj <- ifelse(pcs == 0, 0, 1)
    returned_object <- list(
      pcors = pcors,
      cors = cors,
      trys = trys,
      pcs = pcs,
      adj = adj
    )
    returned_object
  }


symm_mat <- function (x) {
  x[lower.tri(x)] <- t(x)[lower.tri(x)]
  x
}
