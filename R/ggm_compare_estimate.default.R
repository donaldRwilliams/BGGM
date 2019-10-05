#' Compare Edges Between GGMs with the Posterior Distribution
#'
#' @param ... data matrices (\emph{n} by  \emph{p}). Requires at least two.
#' @param iter number of posterior samples
#' @param cred credible interval width used for the decision rule
#' @param rope region of practical equivalence
#'
ggm_compare_estimate.default <- function(..., iter = 5000,
                                         cred = NULL,
                                         rope = NULL){

  ci_width <- cred
  if(is.null( ci_width) & is.null(rope)) stop("ci_width or rope must be specified")
  low <- (1 - ci_width) / 2

  up  <-  1 - low

  info <- Y_combine(...)

  p <- info$dat_info$p[1]

  groups <- length(info$dat)

  inv_mat <- lapply(1:groups, function(x) {Y <- info$dat[[x]];
                                Y <- scale(Y, scale = T);
                                S <- t(Y) %*% Y;
                                n <- nrow(Y);
                                samps <- rWishart(iter, n - 1, solve(S))
                                })


  pcors <- list()
  for(i in 1:groups){

  temp <- lapply(1:iter, function(x) {pcors <-  cov2cor(inv_mat[[i]][,,x]);
                                 pcors[upper.tri(pcors)] })

  pcors[[i]] <- do.call(rbind, temp)

  }

  pcors_diffs <- list()

  for(i in 1:nrow(info$pairwise)){

    temp <- info$pairwise[i,]

    diff <- list(pcors[[temp[1]]] - pcors[[temp[2]]] )

    names(diff) <- paste("Y_g", temp, sep = "", collapse = " vs ")

    pcors_diffs[[i]] <- diff

  }

  name_temp <- matrix(0, p, p)

  name_temp[] <-  unlist(lapply(1:p , function(x) paste( 1:p, x, sep = "--")  ) )


  if(is.null(rope)){

    dat_results <- list()

   for(i in 1:nrow(info$pairwise)){

   diff_ci <-  t(apply(pcors_diffs[[i]][[1]], MARGIN = 2, quantile, c(low, up)))

   diff_mu <-    apply(pcors_diffs[[i]][[1]], MARGIN = 2, mean)

   diff_sd <-    apply(pcors_diffs[[i]][[1]], MARGIN = 2, sd)


   results_temp <- list(cbind.data.frame(data.frame(edge = name_temp[upper.tri(name_temp)],
                                               post_mean =  diff_mu,
                                               post_sd = diff_sd),
                                    diff_ci))

   names(results_temp) <-  names(pcors_diffs[[i]])

   dat_results[[i]] <- results_temp


   }

  }

  if(is.numeric(rope)){

    dat_results <- list()

    for(i in 1:nrow(info$pairwise)){

      pr_in <-  apply(pcors_diffs[[i]][[1]], MARGIN = 2, rope_helper, rope)

      pr_out <- 1 - pr_in

      diff_mu <-    apply(pcors_diffs[[i]][[1]], MARGIN = 2, mean)

      diff_sd <-    apply(pcors_diffs[[i]][[1]], MARGIN = 2, sd)


      results_temp <- list(data.frame(edge = name_temp[upper.tri(name_temp)],
                                                       post_mean =  diff_mu,
                                                       post_sd = diff_sd,
                                                       pr_out = pr_out,
                                                       pr_in  = pr_in))

      names(results_temp) <-  names(pcors_diffs[[i]])

      dat_results[[i]] <- results_temp
      }

  }

  returned_object <- list(dat_results = dat_results,
                          p = p,
                          info = info,
                          iter = iter,
                          call = match.call(),
                          rope = rope,
                          ci_width = ci_width)

  class(returned_object) <- "ggm_compare_estimate"

  returned_object

}


