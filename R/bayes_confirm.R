#' Title
#'
#' @param X 
#' @param nu 
#' @param delta 
#' @param type 
#' @param hyp 
#' @param n_samples 
#' @param cores 
#'
#' @return
#' @export
#'
#' @examples
bayes_confirm <- function(X, nu, delta, type, hyp = NULL, n_samples, cores){
  
  p <- ncol(X)
  off_diag <- (p * (p-1)) / 2
  
  col_names <- numbers2words(1:p)
  pcor_names <- pcor_name_helper(hyp)
  
  pcor_names <- pcor_names[pcor_names != ""]
  
  if(anyNA(suppressWarnings(as.numeric(pcor_names)))) {
    rm_numeric <- which(!is.na(suppressWarnings(as.numeric(pcor_names))))
    pcor_names <- pcor_names[-rm_numeric]
    }
  
  hyp <- stringr::str_remove_all(hyp, "_")
  
    #col_indicator <- unique(unlist(lapply(pcor_names_temp, function(x) word2num(x)[[2]]) ))
  #pcor_names <- colnames(X)[col_indicator]
 
  

  
  pcor_mat <- matrix(0, p, p)
  if(is.null(colnames(X))){
    stop("The variables must be named. These names cannot included numbers")
    
  }
  if(type == "inequality"){
    if(is.null(hyp)){
      stop("inequality constrained hypotheses must be specified")
    }
    
    
    #nu = 16;delta = 2; hyp = hyp;cores = 4; n_samples = 50000
    
    samples <- sampling(X, nu = nu, delta = delta, n_samples = n_samples, cores)
    
    ERr1 <- create_matrices(varnames = pcor_names, hyp = hyp)
    
    
    
    prior <- do.call(rbind, lapply(samples, function(x)  x$fisher_z_prior[,pcor_names]))
    post <-  do.call(rbind, lapply(samples, function(x)  x$fisher_z_post[,pcor_names]))
    
    Sigma0 <- cov(prior) 
    Sigma1 <- cov(post)
    
    mu0 <- ERr1$inequality$R_i %*% colMeans(prior)
    mu1 <- ERr1$inequality$R_i %*% colMeans(post)
    
    s0 <- ERr1$inequality$R_i %*% Sigma0 %*% t(ERr1$inequality$R_i)
    s1 <- ERr1$inequality$R_i %*% Sigma1 %*% t(ERr1$inequality$R_i)
    
    prior_prob_h <- mvtnorm::pmvnorm(lower = as.vector(ERr1$inequality$r_i), upper = Inf,mean = as.vector(mu0),  sigma = s0)[1] 
    post_prob_h <- mvtnorm::pmvnorm(lower = as.vector(ERr1$inequality$r_i), upper = Inf, mean = as.vector(mu1),  sigma = s1)[1] 
    
    prior_prob_not_h <- 1 - prior_prob_h 
    post_prob_not_h <- 1 - post_prob_h
    
    
     BF_1u <- post_prob_h / prior_prob_h
     BF_2u <- post_prob_not_h / prior_prob_not_h
     BF_12 <- BF_1u / BF_2u
    
    
    pcor_mat[upper.tri(pcor_mat)] <- colMeans(samples[[1]]$pcor_post[,1:off_diag])
    pcor_mat[lower.tri(pcor_mat)] <- t(pcor_mat)[lower.tri(pcor_mat)]
    
    colnames(pcor_mat) <- colnames(X)
    
    results <- list(BF_1u = BF_1u, BF_2u = BF_2u, BF_12 = BF_12, pcor_mat = pcor_mat)
   }
  results
}

