
ci_helper <- function(x, ci_width){
  low <- (1 - ci_width) / 2
  up  <-  1 - low
  interval <-  quantile(x, c(low, up))
  as.numeric(ifelse(interval[1] < 0 & interval[2] > 0, 0, 1))
}

inverse_2_beta <- function(fit, samples = 500){

  # check number of samples
  if(samples >= fit$iter){
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
    out <- BGGM:::kl_func(var(y), sigmas^2)
    }
  ci <- quantile(out, prob = c(low, up) )
  mu_out <- mean(out)
  sd_out <- sd(out)
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

  FL <- sum((inv_g1 - theta_rep)^2)

  return <- list(KLD = KLD, JSD = JSD, QL = QL, FL = FL)
}




KL = function(Theta,hatTheta){

  # Kuismin, M., & Sillanpää, M. J. (2016). Use of Wishart prior and simple extensions for
  # sparse precision matrix estimation. PloS one, 11(2), e0148171.


  p = ncol(Theta)

  invTheta = solve(Theta,diag(1,p))

  kl  = 0.5 * (sum(diag(invTheta%*%hatTheta)) - log(det(invTheta%*%hatTheta)) - p)

  return(kl)

}

QL = function(Theta,hatTheta){


  # Kuismin, M., & Sillanpää, M. J. (2016). Use of Wishart prior and simple extensions for
  # sparse precision matrix estimation. PloS one, 11(2), e0148171.

  p = ncol(Theta)
  I = diag(1,p)

  invTheta = solve(Theta,I)

  osa = sum(diag(invTheta%*%hatTheta - I))
  tulos = osa^2

  return(tulos)

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
    post <- post_helper(S = S, n = n, nu = nu, p = p, delta = delta, Psi = Psi, b_inv = b_inv)
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
  inv_wish_prior <- solve(rWishart(1, df =  delta + p - 1, diag(p))[,,1], tol = 1e-20)
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
  chains = cores
  # sample from priors and posteriors
  samples <- foreach(i = 1:chains, .export = c("fisher_z", "sampling_helper",
                                               "numbers2words",
                                               "prior_helper", "post_helper")) %dopar%{
                                                 sampling_helper(X = X, nu = nu, delta = delta,  n_samples = samps[i])

                                               }
  parallel::stopCluster(cl)
  samples
}






fisher_z <- function(rho){
  .5*log((1+rho)/(1-rho))
}

sd_helper <- function(post_samples, prior_at_zero){
  prior_at_zero /  dnorm(0, mean(post_samples), sd(post_samples))
}

equality_test <- function(X, nu, delta, type, pcor_names, hyp = NULL, n_samples, cores){
  p <- ncol(X)
  pcor_mat <- matrix(0, p, p)
  if(is.null(colnames(X))){
    stop("The variables must be named. These names cannot included numbers")

  }
  if(type == "inequality"){
    if(is.null(hyp)){
      stop("inequality constrained hypotheses must be specified")
    }


    # X <- as.matrix(dat_ptsd)
    # nu = 21
    # delta = 1
    # n_samples = 1000
    # cores = 4
    samples <- sampling(X, nu = nu, delta = delta, n_samples = n_samples, cores)

    ERr1 <- create_matrices(varnames = pcor_names, hyp = hyp)

    prior <- do.call(rbind, lapply(samples, function(x)  x$fisher_z_prior[,pcor_names]))
    post <- do.call(rbind, lapply(samples, function(x)  x$fisher_z_post[,pcor_names]))

    Sigma0 <- cov(prior)
    Sigma1 <- cov(post)

    mu0 <- ERr1$inequality$R_i %*% colMeans(prior)
    mu1 <- ERr1$inequality$R_i %*% colMeans(post)

    s0 <- ERr1$inequality$R_i %*% Sigma0 %*% t(ERr1$inequality$R_i)
    s1 <- ERr1$inequality$R_i %*% Sigma1 %*% t(ERr1$inequality$R_i)

    prior_prob <- mvtnorm::pmvnorm(lower = as.vector(ERr1$inequality$r_i), upper = Inf,mean = as.vector(mu0),  sigma = s0)[1]
    post_prob <- mvtnorm::pmvnorm(lower = as.vector(ERr1$inequality$r_i), upper = Inf,mean = as.vector(mu1),  sigma = s1)[1]


    pcor_mat[lower.tri(pcor_mat)] <- colMeans(samples[[1]]$pcor_post)
    pcor_mat[upper.tri(pcor_mat)] <- t(pcor_mat)[upper.tri(pcor_mat)]

    colnames(pcor_mat) <- colnames(X)
    results <- list(BF = post_prob / prior_prob)
    results <- list(results = results,  pcor_mat = pcor_mat)
  }
  results
}



pcor_name_helper <- function(x){
  keep_vars <-  unlist(strsplit(gsub("[^[:alnum:] ]", "", x), " +"))
  keep_vars
}




node_direct_helper_1 <- function(pcor_mat, adj_mat){
  # function to replicate direction of nodes
  # pcor_mat: full pcor mat from first network
  # adj_mat: ae

  # ensure the diagonal
  diag(adj_mat) <- 0

  # list to store hypotheses
  hyp <- list()

  for(i in 1:ncol(pcor_mat) ){

    ith_row <- pcor_mat[i,] * adj_mat[i,]
    ith_row_sign <- sign(ith_row)

    # if only positive
    if(sum(ith_row_sign == 1) != 0 & sum(ith_row_sign == -1) == 0){

      # which is positive
      which_positve <- which(ith_row_sign ==  1)

      # change numbers to words
      change_2_words <- numbers2words(which_positve)

      # join i row to positive edges
      join_ith_row <- paste(numbers2words(i), change_2_words, sep = "_")

      # store
      hyp[[i]] <-  paste("(", paste(join_ith_row, collapse = ", "), ")", "> 0")
    }

    if(sum(ith_row_sign == -1) != 0 & sum(ith_row_sign == 1) == 0){

      # which is negative
      which_negative <- which(ith_row_sign ==  -1)

      # change numbers to words
      change_2_words <- numbers2words(which_negative)

      # join i row to positive edges
      join_ith_row <- paste(numbers2words(i), change_2_words , sep = "_")

      # store
      hyp[[i]] <-  paste("(", paste(join_ith_row, collapse = ", "), ")", "> 0")

    }



    if(sum(ith_row_sign == -1) != 0 & sum(ith_row_sign == 1) != 0){

      # which is positive
      which_positve <- which(ith_row_sign ==  1)

      # change numbers to words
      change_2_words_pos <- numbers2words(which_positve)

      # join i row to positive edges
      join_ith_row_pos <- paste(numbers2words(i), change_2_words_pos,  sep = "_")

      # store
      hyp_pos <-  paste("(", paste(join_ith_row_pos, collapse = ", "), ")", "> 0")



      # which is negative
      which_negative <- which(ith_row_sign ==  -1)

      # change numbers to words
      change_2_words_neg <- numbers2words(which_negative)

      # join i row to positive edges
      join_ith_row_neg <- paste(numbers2words(i), change_2_words_neg , sep = "_")



      hyp_neg <- paste(">", " ( ", paste(join_ith_row_neg, collapse = ", "), " )", sep = "")

      hyp[[i]] <-  paste(hyp_pos, hyp_neg)

    }

  }

  hyp
}





#################################
######## create matrices ########
#################################
create_matrices <- function(varnames, hyp){

  #varnames <- variable.names(object) #provides the variable names of the linear model object, including intercept
  if(is.null(varnames)) stop("Please input proper linear model object")
  varnames <- gsub("(\\(Intercept\\))", "Intercept", varnames) #remove parentheses around intercept so these don't conflict later

  hyp2 <- gsub("[ \n]", "", hyp) #removes all whitespace
  hyp2 <- gsub("(\\(Intercept\\))", "Intercept", hyp2) #If hyp re. intercept input with surrounding parentheses, remove
  if(!grepl("^[0-9a-zA-Z><=,().-]+$", hyp2)) stop("Impermissable characters in hypotheses") #Self-explanatory. NEW parentehese
  if(grepl("[><=]{2,}", hyp2)) stop("Do not use combined comparison signs e.g., '>=' or '=='")

  step1 <- unlist(strsplit(hyp2, split = "[<>=,()]")) #split by special characters and unlist
  input_vars <- step1[grep("[a-zA-Z]+", step1)] #extract subunits that contain at least one letter
  #if(!all(input_vars %in% varnames)) stop("Hypothesis variable(s) not in object, check spelling") #Checks if input variables exist in lm-object

  framer <- function(x){ #As function because same code used once more later
    pos_comparisons <- unlist(gregexpr("[<>=]", x)) #Gives the positions of all comparison signs
    leftside <- rep(NA, length(pos_comparisons) + 1) #empty vector for loop below
    rightside <- rep(NA, length(pos_comparisons) + 1) #empty vector for loop below
    pos1 <- c(-1, pos_comparisons) #positions to extract data to the leftside of comparisons
    pos2 <- c(pos_comparisons, nchar(x) + 1) #positions to extract data to the rightside of comparisons
    for(i in seq_along(pos1)){
      leftside[i] <- substring(x, pos1[i] + 1, pos1[i+1] - 1) #Extract all variables or outcomes to the leftside of a comparison sign
      rightside[i] <- substring(x, pos2[i] + 1, pos2[i+1] - 1) #Extract all variables or outcomes to the rightside of a comparison sign
    }
    leftside <- leftside[-length(leftside)] #remove last element which is a NA due to loop formatting
    rightside <- rightside[-length(rightside)] #remove last element which is a NA due to loop formatting
    comparisons <- substring(x, pos_comparisons, pos_comparisons) #Extract comparison signs
    data.frame(left = leftside, comp = comparisons, right = rightside, stringsAsFactors = FALSE) #hypotheses as a dataframe
  }

  framed <- framer(hyp2) #hypotheses as a dataframe

  if(any(grepl(",", framed$left)) || any(grepl(",", framed$right))){ #Larger loop that deals with commas if the specified hypothesis contains any
    if(nrow(framed) > 1){
      for(r in 1:(nrow(framed)-1)){ #If a hypothesis has been specified with commas e.g., "X1 > 0, X2 > 0" or "(X1, X2) > X3"
        if(all.equal(framed$right[r], framed$left[r+1])){ #The right hand side of the hypothesis df will be equal to the next row left side
          if(substring(framed$right[r], 1, 1) == "(") { #If the first row begins with a ( as when "X1 > (X2, X3)" and opposed to "(X2, X3) > X1"
            framed$right[r] <- sub("),.+", ")", framed$right[r])#If so, remove everything to the right of the parenthesis on the right hand side
            framed$left[r+1] <- sub(".+),", "", framed$left[r +1])#and everything to the left of the parenthesis on the left hand side to correct the df
          } else{
            framed$right[r] <- sub(",.+", "", framed$right[r]) #else, remove everything to the right of the comma on the right hand side
            framed$left[r+1] <- sub("[^,]+,", "", framed$left[r+1]) #and everything to the left of the comma on the left hand side to correct the df
          }
        }
      }
    }

    commas_left <- framed$left[grep(",", framed$left)] #At this point all remaining elements that contain commas should also have parentheses, check this
    commas_right <- framed$right[grep(",", framed$right)] #Necessary to use is isTRUE below in case one of these contains no commas, and 'any' for several rows
    if(isTRUE(any(!grepl("\\(.+)", commas_left))) || isTRUE(any(!grepl("\\(.+)", commas_right))) || #Check so rows contain parenthesis
       isTRUE(any(grepl(").+", commas_left))) || isTRUE(any(grepl(").+", commas_right))) || #Check so parentheses are not followed by anything
       isTRUE(any(grepl(".+\\(", commas_left))) || isTRUE(any(grepl(".+\\(", commas_right)))) { #chekc so parentheses are not preceded by anything
      stop("Incorrect hypothesis syntax or extra character, check specification")
    }


    framed$left <- gsub("[()]", "", framed$left) #drop remaining parentheses
    framed$right <- gsub("[()]", "", framed$right)
    commas <- unique(c(grep(",", framed$left), grep(",", framed$right))) #Gives us the unique rows that still contain commas (multiple comparisons) from left or right columns

    if(length(commas) > 0){ #If there are any multiple comparisons e.g., (X1, X2) below loop separates these in
      multiples <- vector("list", length = length(commas)) #Empty vector to store results for each row in loop below

      for(r in seq_along(commas)){ #for each row containing commas
        several <- framed[commas,][r, ] #select row r

        if(several$comp == "="){ #If e.g., (X1, X2) = X3, convert to X1 = X2 = X3

          several <- c(several$left, several$right)
          separate <- unlist(strsplit(several, split = ",")) #split by special characters and unlist
          if(any(grepl("^$", several))) stop("Misplaced comma in hypothesis") #if empty element
          converted_equality <- paste(separate, collapse = "=") #convert to X1 = X2 = X3 shape
          multiples[[r]] <- framer(converted_equality) #hypotheses as a dataframe

        } else{ #If inequality comparison
          leftvars <- unlist(strsplit(several$left, split = ",")) #separate left hand var
          rightvars <- unlist(strsplit(several$right, split = ",")) #separate right hand vars
          if(any(grepl("^$", leftvars)) || any(grepl("^$", rightvars))) stop("Misplaced comma in hypothesis") #if empty element

          left <- rep(leftvars, length.out = length(rightvars)*length(leftvars)) #repeat each leftvars the number of rightvars
          right <- rep(rightvars, each = length(leftvars)) #complement for rightvars
          comp <- rep(several$comp, length(left)) #repeat the comparison a corresponding number of times

          multiples[[r]] <- data.frame(left = left, comp = comp, right = right, stringsAsFactors = FALSE) #save as df to be able to combine with 'framed'
        }
      }

      framed <- framed[-commas,] #remove old unfixed rows with commas
      multiples <- do.call(rbind, multiples) #make list into dataframe
      framed <- rbind(multiples, framed) #recombine into one dataframe
    }
  } #end comma loop

  equality <- framed[framed$comp == "=",]
  inequality <- framed[!framed$comp == "=",]

  #********Equality
  if(nrow(equality) == 0) { #If there are no '=' comparisons set to NULL
    list_equality <- NULL
  } else{
    outcomes <- suppressWarnings(apply(equality[, -2], 2, as.numeric)) #Convert left/right to numeric, non-numeric values (variables) coerced to NA
    outcomes <- matrix(outcomes, ncol = 2, byrow = FALSE) #Conversion to matrix in case there was only one row in outcomes
    if(any(rowSums(is.na(outcomes)) == 0)) stop("Value compared with value rather than variable, e.g., '2 = 2', check hypotheses")
    rows <- which(rowSums(is.na(outcomes)) < 2) #which rows contain a numeric value (comparing variable to value), that is not two NA-values
    specified <- t(outcomes[rows,]) #transpose so that specified comparison values are extracted in correct order below, e.g, in case when "X1 = 0, 2 = X2"
    specified <- specified[!is.na(specified)] #extract specified comparison values
    r_e <- ifelse(rowSums(is.na(outcomes)) == 2, 0, specified) #If variable = variable -> 0, if variable = value -> value
    r_e <- matrix(r_e) #convert to matrix

    var_locations <- apply(equality[, -2], 2, function(x) ifelse(x %in% varnames, match(x, varnames), 0)) #convert non-variables to 0 and others are given their locations in varnames
    var_locations <- matrix(var_locations, ncol = 2) #Necessary if only one comparison row

    R_e <- matrix(rep(0, nrow(equality)*length(varnames)), ncol = length(varnames)) #Create empty variable matrix

    for(i in seq_along(r_e)){ # for each row i in R_e, replace the columns specified in var_locations row i
      if(!all(var_locations[i, ] > 0)){ #If only one variable is specified (i.e., other one is set to zero)
        R_e[i, var_locations[i,]] <- 1 #Set this variable to 1 in R_e row i
      } else{ #If two variables specified
        R_e[i, var_locations[i,]] <- c(1, -1) #Set one column to 1 and the other to -1 in R_e row i
      }
    }
    list_equality <- list(R_e = R_e, r_e = r_e) #Note column 1 in R_e is for intercept
  }

  #********Inequality
  if(nrow(inequality) == 0) { #If there are no '>' or '<' comparisons set to NULL
    list_inequality <- NULL
  } else{
    outcomes <- suppressWarnings(apply(inequality[, -2], 2, as.numeric)) #Convert left/right to numeric, non-numeric values (variables) coerced to NA
    outcomes <- matrix(outcomes, ncol = 2, byrow = FALSE) #Conversion to matrix in case there was only one row in outcomes
    if(any(rowSums(is.na(outcomes)) == 0)) stop("Value compared with value rather than variable, e.g., '2 > 2', check hypotheses")
    cols <- which(rowSums(is.na(outcomes)) < 2) #which columns contain a numeric value (comparing variable to value), that is not two NA-values
    specified <- t(outcomes[cols,]) #transpose so that specified comparison values are extracted in correct order below
    specified <- specified[!is.na(specified)] #extract specified comparison values
    r_i <- ifelse(rowSums(is.na(outcomes)) == 2, 0, specified) #If variable = variable -> 0, if variable = value -> value
    r_i <- matrix(r_i) #convert to matrix

    leq <- which(inequality$comp == "<") #gives the rows that contain '<' (lesser or equal) comparisons
    var_locations <- apply(inequality[, -2], 2, function(x) ifelse(x %in% varnames, match(x, varnames), 0)) #convert non-variables to 0 and others are given their locations
    var_locations <- matrix(var_locations, ncol = 2) #Necessary if only one comparison row

    R_i <- matrix(rep(0, nrow(inequality)*length(varnames)), ncol = length(varnames)) #Create empty variable matrix

    for(i in seq_along(r_i)){ # for each row i in R_i, replace the columns specified in var_locations row i
      if(!all(var_locations[i, ] > 0)){ #If only one variable is specified (i.e., other one is set to zero)

        if(var_locations[i, 1] == 0){ #If first value is not the variable (i.e, a comparison value)
          if(i %in% leq){#Then if comparison is 'lesser or equal'
            value <-  1  #set variable value to 1
          } else{ #else if comparison 'larger or equal'
            r_i[i] <- r_i[i]*-1 #invert comparison value
            value <- -1 # set variable value to -1
          }
        } else{ #If first value is the variable (i.e., the second is a comparison value)
          if(i %in% leq){ #then if comparison is 'lesser or equal'
            r_i[i] <- r_i[i]*-1 #invert comparison value
            value <-  -1  #set variable value to -1
          } else{
            value <- 1 #else if comparison 'larger or equal' set to 1
          }
        }

        R_i[i, var_locations[i,]] <- value #Set this variable to 1/-1 in R_i row i

      } else{ #If two variables specified
        value <- if(i %in% leq) c(-1, 1) else c(1, -1) #If comparison is 'leq' take var2 - var1, if 'larger or equal' take var1 - var2
        R_i[i, var_locations[i,]] <- value #Set one column to 1 and the other to -1 in R_i row i
      }
    }

    list_inequality<- list(R_i = R_i, r_i = r_i) #Note column 1 in R_i is for intercept
  }

  matrices <- list(equality = list_equality, inequality = list_inequality); matrices #final output

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


