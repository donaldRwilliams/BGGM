#' Edge (Partial Correaltion) Contrasts with the Bayes Factor
#'
#' @description Allows for comparing partial correlations \emph{within} the same GGM--e.g., to determine the largest edge for each node. The differences are tested
#' with the Bayes factor. The test can be either directional, for example whether one edge is larger or smaller than another, or two sided testing is also possible. Further,
#' there is also an exhaustive option that provides the posterior probability of the null, greater than zero, and less than zero.
#'
#'
#' @param x object from \code{explore}
#' @param contrast partial correlations to compare
#' @param alternative \code{greater}, \code{less}, \code{two.sided}, or \code{exhaustive}
#'
#' @return list of class \code{edge_compare.explore}:
#' \itemize{
#' \item \code{results} data.frame summary of the edge contrasts
#' \item \code{alternative} \code{greater}, \code{less}, \code{two.sided}, or \code{exhaustive}
#' \item \code{call} match.call()
#' }
#'
#'
#' @export
#'
#' @examples
#'
#' # p = 10
#' Y <- BGGM::bfi[,1:10]
#'
#' # sample from posterior
#' fit_bf <- explore(Y, prior_sd = 0.5,
#'               iter = 5000,
#'               cores = 2)
#'
#' # edge compare
#'  edge_comp <- edge_compare(fit_bf,
#'                          contrast = c("1--5 - 1--3",
#'                                       "1--2 - 1--6",
#'                                       "1--4 - 1--7",
#'                                       "1--5 - 1--10",
#'                                       "1--2 - 1--9"),
#'                          alternative = "exhaustive")
#' # summary
#' summary(edge_comp)
edge_compare.explore <- function(x, contrast, alternative){

  p <- x$p


  # names
  mat_names <- null_mat <- pos_mat <- neg_mat <- matrix(0, x$p, x$p)

  mat_names[] <-  unlist(lapply(1:x$p, function(z) paste(1:x$p, z, sep = "--")))


  # tranformed
  posterior_samples_z <- do.call(rbind.data.frame,
                                 lapply(1:x$cores, function(z)  x$samples[[z]]$fisher_z_post))
  prior_samples_z <-    do.call(rbind.data.frame,
                                lapply(1:x$cores, function(z)  x$samples[[z]]$fisher_z_prior))

  colnames(posterior_samples_z) <- c(mat_names[upper.tri(mat_names)],mat_names[lower.tri(mat_names)] )
  colnames(prior_samples_z) <- c(mat_names[upper.tri(mat_names)],mat_names[lower.tri(mat_names)] )


  # not transformed
  posterior_samples <- do.call(rbind.data.frame,
                               lapply(1:x$cores, function(z)  x$samples[[z]]$pcor_post))
  prior_samples <-    do.call(rbind.data.frame,
                              lapply(1:x$cores, function(z)  x$samples[[z]]$pcor_prior))

  colnames(posterior_samples) <- c(mat_names[upper.tri(mat_names)],mat_names[lower.tri(mat_names)] )
  colnames(prior_samples) <- c(mat_names[upper.tri(mat_names)],mat_names[lower.tri(mat_names)] )

  dat <- list()

  for(i in 1:length(contrast)){

    one <- sub(".* - ","", contrast[i])

    two <- gsub(" .*","", contrast[i])

    if(anyNA(match(c(one, two), colnames(posterior_samples)))){

      stop("rewrite the contrasts. must be 1--2 - 1--3")

    }

    post_diff_z <- posterior_samples_z[,two] - posterior_samples_z[,one]
    prior_diff_z <- prior_samples_z[,two] - prior_samples_z[,one]


    post_diff <- posterior_samples[,two] - posterior_samples[,one]
    prior_diff <- prior_samples[,two] - prior_samples[,one]




    if(alternative == "two.sided"){
      BF <- dnorm(0, mean(prior_diff_z), sd(prior_diff_z)) / dnorm(0, mean(post_diff_z), sd(post_diff_z))
      prob <- BF / (BF + 1)

      dat[[i]] <- data.frame(contrast = contrast[i],
                             post_mean = mean(post_diff),
                             post_sd = sd(post_diff),
                             "BF 10" = BF,
                             "p(H0|Y)" = 1 - prob,
                             "p(H1|Y)" = prob ,
                             check.names = F)

    }



    if(alternative == "greater"){

      dens_greater <- (1  - pnorm(0, mean(post_diff_z), sd(post_diff_z))) * 2

      BF <- dnorm(0, mean(prior_diff_z), sd(prior_diff_z)) / dnorm(0, mean(post_diff_z), sd(post_diff_z))

      BF_20 <- BF * dens_greater
      prob <- BF_20 / (BF_20 + 1)

      dat[[i]] <- data.frame(contrast = contrast[i],
                             post_mean = mean(post_diff),
                             post_sd = sd(post_diff),
                             "BF 10" = BF_20,
                             "p(H0|Y)" = 1 - prob,
                             "p(H1|Y)" = prob ,
                             check.names = F)
    }


    if(alternative == "less"){

      dens_less <- (pnorm(0, mean(post_diff_z), sd(post_diff_z))) * 2

      BF <- dnorm(0, mean(prior_diff_z), sd(prior_diff_z)) / dnorm(0, mean(post_diff_z), sd(post_diff_z))

      BF_20 <- BF * dens_less
      prob <- BF_20 / (BF_20 + 1)

      dat[[i]] <- data.frame(contrast = contrast[i],
                             post_mean = mean(post_diff),
                             post_sd = sd(post_diff),
                             "BF 10" = BF_20,
                             "p(H0|Y)" = 1 - prob,
                             "p(H1|Y)" = prob ,
                             check.names = F)
    }



    if(alternative == "exhaustive"){

      BF <- dnorm(0, mean(prior_diff_z), sd(prior_diff_z)) / dnorm(0, mean(post_diff_z), sd(post_diff_z))

      dens_greater <- (1  - pnorm(0, mean(post_diff_z), sd(post_diff_z))) * 2

      BF_greater <- BF * dens_greater

      dens_less <- (pnorm(0, mean(post_diff_z), sd(post_diff_z))) * 2

      BF_less <- BF * dens_less


      BF_mat <- data.frame(BF_01 = 1/ BF,  BF_greater = BF_greater, BF_less = BF_less)

      post_prob <-  data.frame( t(apply(BF_mat, 1, FUN = function(x) {x / sum(x) })) )


      colnames(post_prob) <- c("p(H0|Y)", "p(H1|Y)", "p(H2|Y)")
      row.names(post_prob) <- c()
      post_prob <-  round(post_prob, 3)
      dat[[i]] <- cbind.data.frame(contrast = contrast[i],
                                   post_mean = mean(post_diff),
                                   post_sd = sd(post_diff),
                                   post_prob)

    }


  }
  returned_object <- list(results = do.call(rbind, dat),
                          alternative = alternative,
                          call = match.call())
  class(returned_object) <- "edge_compare.explore"
  returned_object
}
