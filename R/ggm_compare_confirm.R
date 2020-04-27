

# sample posterior
post_samp <- lapply(1:groups, function(x) {

  Y <- as.matrix(scale(info$dat[[x]], scale = FALSE))

  .Call(
    '_BGGM_Theta_continuous',
    PACKAGE = 'BGGM',
    Y = Y,
    iter = iter + 50,
    delta = delta,
    epsilon = 0.1,
    prior_only = 0,
    explore = 1
  )
})

# sample prior
prior_samp <- lapply(1:groups, function(x) {

   Y <- info$dat[[x]]

  .Call(
    '_BGGM_Theta_continuous',
    PACKAGE = 'BGGM',
    Y = Y,
    iter = 10000,
    delta = delta,
    epsilon = 0.1,
    prior_only = 1,
    explore = 0
  )$fisher_z
})


# variables
p <- ncol(Y1)

# I_p
I_p <- diag(p)

# number of edges
pcors <- p*(p-1)*0.5

I_p <- diag(p)

col_names <- BGGM:::numbers2words(1:p)

mat_name <- sapply(col_names, function(x) paste(col_names,x, sep = ""))[upper.tri(I_p)]


temp_names <- lapply(1:groups, function(x) paste0("g", BGGM:::numbers2words(x),  mat_name))

post_group <- matrix(post_samp[[1]]$fisher_z[, , 51:(iter+50)][upper.tri(I_p)],
                            iter, pcors,
                            byrow = TRUE)

prior_group <-  matrix(prior_samp[[1]][ , ,][upper.tri(I_p)], iter, pcors, byrow = TRUE)


# combined groups
for(j in 2:(groups)){
  post_group <-  cbind(post_group,
                       matrix(post_samp[[1]]$fisher_z[, , 51:(iter+50)][upper.tri(I_p)],
                              iter, pcors,
                              byrow = TRUE))
  prior_group <-  cbind(prior_group,
                        matrix(prior_samp[[1]][ , ,][upper.tri(I_p)], iter, pcors, byrow = TRUE))
}


posterior_samples <- post_group
colnames(posterior_samples) <- unlist(temp_names)

prior_samples <- prior_group
colnames(prior_samples) <- unlist(temp_names)



prior_mu <- colMeans(prior_samples)

prior_cov <- cov(prior_samples)

post_mu <- colMeans(posterior_samples)

post_cov <- cov(posterior_samples)

Y <- cbind(Y1, Y2)

hypothesis <- c("g1_A--B > g1_A--C")
colnames(Y1) <- LETTERS[1:10]
colnames(Y2) <- LETTERS[1:10]

hypothesis <- BGGM:::convert_hyps(hypothesis = hypothesis, cbind(Y1))
hypothesis

hyp <- gsub(BGGM:::hyp_converter(hypothesis)$hyp_converted, pattern = "_", replacement = "")
hyp


BFprior <- BF(prior_mu,
              Sigma = prior_cov,
              hypothesis = hyp,
              n = 1)


BFpost <- BF(post_mu,
             Sigma = post_cov,
             hypothesis = hyp,
             n = 1)




