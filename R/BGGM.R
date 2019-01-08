#' Title
#'
#' @param x data matrix
#' @param chains number of chains
#' @param iter number of posterior samples for each chain
#' @param burnin discarded samples during warmup
#'
#' @return
#' @export
#'
#' @examples
#' fit <- BGGM(X)

BGGM <- function(x, chains = 2, iter = 1000, burnin = 100){
  x <- as.matrix(x)
  ## model
  mod = "model{
  for ( i in 1:n) {
  x[i,1:p] ~ dmnorm( mu[1:p] , cov_inv[1:p,1:p] )
  }
  for ( i in 1:p){
  mu[i] <- 0
  }
  cov_inv ~ dwish(r_mat[1:p,1:p], wish_df)


  for(i in 1:p) {
  inv_sd[i] <- sqrt(cov_inv[i,i])
  }

  for ( i in 1:p ) {
  for ( j in 1:p ) {
  pcors[i,j] <- ( (-1 * cov_inv[i,j]) / (inv_sd[i]*inv_sd[j]) )
  }
  }
}
"
dat <- scale(x)
## data list
data_list = list(x = dat,
                 n = nrow(x),
                 wish_df = ncol(x),
                 p = ncol(x),
                 r_mat = diag(x=1,nrow=ncol(x)))




jags_model <- rjags::jags.model(textConnection(mod),
                                data = data_list,
                                n.adapt = burnin,
                                n.chains = chains,
                                quiet = T)

mcmc_samples <- rjags::coda.samples(model =  jags_model,variable.names = c("cov_inv", "pcors"),
                                    n.iter = iter, data = data_list)

parcor_mat <- matrix(0, ncol = ncol(x), ncol(x))
inv_mat <- prob_mat <- parcor_mat

# posterior samples
df_samps <- do.call(rbind.data.frame, mcmc_samples)


parcor_mat[] <- colMeans(df_samps[,  grep("pcors", colnames(df_samps))])
diag(parcor_mat) <- 0


prob_mat[] <- apply(df_samps[,  grep("pcors", colnames(df_samps))], 2,
                    FUN = function(x) min(mean(x > 0), mean(x < 0)))


inv_mat[]   <- colMeans(df_samps[,  grep("cov_inv", colnames(df_samps))])

list(parcors_mat = parcor_mat, inv_mat = inv_mat, posterior_samples = df_samps,
     post_prob = prob_mat, p = ncol(x), dat = dat, iter = chains * iter)

}
