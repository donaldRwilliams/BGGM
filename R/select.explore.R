#' Title
#' @param x data matrix
#' @param BF_cut evidentiary threshold
#' @param alternative type of hypothesis (see notes)
#' @param hyp_prob posterior probability threshold (see notes)
#'
#' @return
#' @export
#'
#' @examples
select.explore <- function(x,
                           BF_cut = 3,
                           alternative = "two.sided",
                           hyp_prob = NULL){




  posterior_samples <- do.call(rbind.data.frame,
                               lapply(1:x$cores, function(z)  x$samples[[z]]$fisher_z_post))

  prior_samples <- unlist(do.call(rbind.data.frame,
                                  lapply(1:x$cores, function(z)  x$samples[[z]]$fisher_z_prior)))

  if(alternative == "two.sided"){

    BF_10_mat <- BF_01_mat <- matrix(0, x$p, x$p)



    prior_dens <- dnorm(0, mean(prior_samples), sd(prior_samples))

    BF_10 <- apply(posterior_samples, MARGIN = 2, FUN = function(z){prior_dens / dnorm(0, mean(z), sd(z))} )
    BF_01 <- 1/ BF_10

    BF_10_mat[upper.tri(BF_10_mat)] <- BF_10[1:x$edge]
    BF_10_mat <- BGGM:::symmteric_mat(BF_10_mat)

    BF_01_mat <- 1 / BF_10_mat


    Adj_10 <- ifelse(BF_10_mat > BF_cut, 1, 0)
    Adj_01 <- ifelse(BF_10_mat < 1 / BF_cut, 1, 0)
    diag(Adj_01) <- 0
    diag(BF_01_mat) <- 0
    returned_object = list(partials_non_zero = x$parcors_mat * Adj_10,
                           pcor_mat = x$parcors_mat,
                           pcor_sd = x$parcors_sd,
                           Adj_10 = Adj_10,
                           Adj_01 = Adj_01,
                           BF_10 = BF_10_mat,
                           BF_01 = BF_01_mat,
                           BF_cut = BF_cut,
                           alternative = alternative,
                           call = match.call())
  }
  if(alternative == "greater"){

    prior_dens <- dnorm(0, mean(prior_samples), sd(prior_samples))

    BF_20_mat <- BF_01_mat <- matrix(0, x$p, x$p)

    BF_10 <- apply(posterior_samples, MARGIN = 2, FUN = function(z){prior_dens / dnorm(0, mean(z), sd(z)) } )

    dens_greater <-  apply(posterior_samples, MARGIN = 2, FUN = function(z){(1  - pnorm(0, mean(z), sd(z))) * 2} )

    BF_20 <- dens_greater * (BF_10)

    BF_01 <- 1 / BF_10

    BF_20_mat[upper.tri(BF_20_mat)] <- BF_20[1:x$edge]

    BF_20_mat <- BGGM:::symmteric_mat(BF_20_mat)

    BF_01_mat[upper.tri(BF_01_mat)] <- BF_01[1:x$edge]

    BF_01_mat <- BGGM:::symmteric_mat(BF_01_mat)

    Adj_20 <- ifelse(BF_20_mat > BF_cut, 1, 0)
    Adj_01 <- ifelse(BF_01_mat > BF_cut, 1, 0)

    diag(Adj_01) <- 0

    diag(BF_01_mat) <- 0


    returned_object = list(partials_positive = x$parcors_mat * Adj_20,
                           pcor_mat = x$parcors_mat,
                           pcor_sd = x$parcors_sd,
                           Adj_01 = Adj_01,
                           Adj_20 =  Adj_20,
                           BF_20 = BF_20_mat,
                           BF_01 = BF_01_mat,
                           BF_cut = BF_cut,
                           alternative = alternative,
                           call = match.call())


  }
  if(alternative == "less")  {

    prior_dens <- dnorm(0, mean(prior_samples), sd(prior_samples))

    BF_20_mat <- BF_01_mat <- matrix(0, x$p, x$p)

    BF_10 <- apply(posterior_samples, MARGIN = 2, FUN = function(z){prior_dens / dnorm(0, mean(z), sd(z))} )

    dens_less <-  apply(posterior_samples, MARGIN = 2, FUN = function(z){(pnorm(0, mean(z), sd(z))) * 2} )

    BF_20 <- dens_less * (BF_10)

    BF_01 <- 1 / BF_10

    BF_20_mat[upper.tri(BF_20_mat)] <- BF_20[1:x$edge]

    BF_20_mat <- BGGM:::symmteric_mat(BF_20_mat)

    BF_01_mat[upper.tri(BF_01_mat)] <- BF_01[1:x$edge]

    BF_01_mat <- BGGM:::symmteric_mat(BF_01_mat)

    Adj_20 <- ifelse(BF_20_mat > BF_cut, 1, 0)
    Adj_01 <- ifelse(BF_10_mat < 1 / BF_cut, 1, 0)

    diag(Adj_01) <- 0

    diag(BF_01_mat) <- 0


    returned_object = list(partials_negative = x$parcors_mat * Adj_20,
                           pcor_mat = x$parcors_mat,
                           pcor_sd = x$parcors_sd,
                           Adj_01 = Adj_01,
                           Adj_20 = Adj_20,
                           BF_20 = BF_20_mat,
                           BF_01 = BF_01_mat,
                           BF_cut = BF_cut,
                           alternative = alternative,
                           call = match.call())

  }
  if(alternative == "exhaustive"){

    if(is.null(hyp_prob)){
      stop("posterior probability must be specificed \n for exhaustive hypothesis testing")

    }

    mat_names <- null_mat <- pos_mat <- neg_mat <- matrix(0, x$p, x$p)

    mat_names[] <-  unlist(lapply(1:x$p, function(z) paste(1:x$p, z, sep = "--")))

    prior_dens <- dnorm(0, mean(prior_samples), sd(prior_samples))



    BF_10 <- apply(posterior_samples, MARGIN = 2, FUN = function(z){prior_dens / dnorm(0, mean(z), sd(z))} )

    dens_less <-  apply(posterior_samples, MARGIN = 2, FUN = function(z){(pnorm(0, mean(z), sd(z))) * 2} )

    BF_less <- dens_less * (BF_10)

    BF_01 <- 1 / BF_10


    dens_greater <-  apply(posterior_samples, MARGIN = 2, FUN = function(z){(1  - pnorm(0, mean(z), sd(z))) * 2} )

    BF_greater <- dens_greater * (BF_10)


    BF_mat <- data.frame(BF_01 = BF_01,  BF_greater = BF_greater, BF_less = BF_less)
    post_prob <-  data.frame( t(apply(BF_mat, 1, FUN = function(x) {x / sum(x) })) )
    colnames(post_prob) <- c("prob_zero", "prob_greater", "prob_less")
    row.names(post_prob) <- c()
    post_prob <-  round(post_prob, 3)
    post_prob <- cbind.data.frame(edge = mat_names[upper.tri(mat_names)], post_prob)[1:x$edge,]


    null_mat[upper.tri(null_mat)] <- ifelse(post_prob$prob_zero > hyp_prob, post_prob$prob_zero, 0)
    null_mat <- BGGM:::symmteric_mat(null_mat)

    pos_mat[upper.tri(pos_mat)] <-  ifelse(post_prob$prob_greater >  hyp_prob, post_prob$prob_greater, 0)
    pos_mat <- BGGM:::symmteric_mat(pos_mat)

    neg_mat[upper.tri(neg_mat)]  <-  ifelse(post_prob$prob_less > hyp_prob, post_prob$prob_less , 0)
    neg_mat <- BGGM:::symmteric_mat(neg_mat)


    returned_object <- list(post_prob = post_prob,
                            neg_mat = neg_mat,
                            pos_mat = pos_mat,
                            null_mat = null_mat,
                            alternative = alternative,
                            pcor_mat = x$parcors_mat,
                            pcor_sd = x$parcors_sd,
                            call = match.call(),
                            prob = hyp_prob)
}

  class(returned_object) <- "select.explore"
  returned_object
}

