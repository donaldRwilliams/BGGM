#' Compare Partial Correlations with the Posterior Distribution
#' @name ggm_compare_estimate
#'
#' @description Compare edges (partial correlations) that are estimated from groups to, say, detect a differences or equivalence.
#'
#' @param ... matrices (or data frame) of dimensions \emph{n} (observations) by  \emph{p} (variables).
#' Requires at least two.
#'
#' @param type character string. Which type of data for \strong{Y} ? The options include \code{continuous},
#' \code{binary}, or \code{ordinal}. See the note for further details.
#'
#' @param iter number of iterations (posterior samples; defaults to 5000).
#'
#' @param analytic logical. Should the analytic solution be computed (default is \code{FALSE}) ?
#'
#' @return
#' A list of class \code{ggm_compare_estimate} containing:
#'  \itemize{
#'  \item \code{pcor_diffs} partial correlation differences (posterior distribution)
#'  \item \code{p} number of variable
#'  \item \code{info} list containing information about each group (e.g., sample size, etc.)
#'  \item \code{iter} number of posterior samples
#'  \item \code{call} \code{match.call}
#'  }
#'
#' @note The work flow for most functions in \strong{BGGM} is to first fit the model
#' and then select the graph  (in this case the differences) with  \code{\link{select}}.
#'
#' @seealso \code{\link{select.ggm_compare_estimate}}
#'
#' @examples
#' # data
#' Y1 <- BGGM::bfi[1:500,1:5]
#' Y2 <- BGGM::bfi[501:1000, 1:5]
#'
#' # fit model
#' fit <- ggm_compare_estimate(Y1, Y2)
#'
#' # posterior summary of differences
#' summary(fit)
#'
#' # select (threshold) with credible intervals
#' sel <- select(fit)
#'
#' # summary
#' summary(sel)
#'
#'# selected differences
#' sel$mat_pcor
#'
#' # adjacency matrix
#' sel$mat_adj
#' @export
ggm_compare_estimate <- function(...,
                                 type = "continuous",
                                 analytic = FALSE,
                                 prior_sd = 0.25,
                                 iter = 5000){

  if(type == "continuous"){

    info <- Y_combine(...)

    p <- info$dat_info$p[1]

    groups <- length(info$dat)

    comparisons <- nrow(info$pairwise)

  if(groups < 2){

    stop("must have (at least) two groups")

    }


  if(!analytic){

    delta <- delta_solve(prior_sd)

    post_samp <- lapply(1:groups, function(x){

    Y <- info$dat[[x]]

     .Call('_BGGM_Theta_continuous',
                        PACKAGE = 'BGGM',
                        Y = Y,
                        iter = iter + 50,
                        delta = delta,
                        epsilon = 0.1,
                        prior_only = 0,
                        explore = 1)$pcors
    })

   diff <- lapply(1:comparisons, function(x) {

     contrast <- info$pairwise[x,]

     post_samp[[contrast[[1]]]][,,(51:iter + 50)] - post_samp[[contrast[[2]]]][,,(51:iter + 50)]

     })


    names(diff)  <- sapply(1:comparisons,function(x)  paste("Y_g",
                                          info$pairwise[x,],
                                          sep = "",
                                          collapse = " - "))


    analytic <- FALSE
    type <- "continuous"

    # returned object
    returned_object <- list(
      diff = diff,
      p = p,
      info = info,
      iter = iter,
      analytic = analytic,
      type = type,
      call = match.call()
    )

    } else {

      z_stat <- lapply(1:comparisons, function(x){

        contrast <- info$pairwise[x,]

        g1 <- analytic_solve(info$dat[[contrast[[1]]]])
        g2 <- analytic_solve(info$dat[[contrast[[2]]]])

        z_stat <- abs((g1$inv_map - g2$inv_map) /   sqrt(g1$inv_var + g2$inv_var))

        })


      diff <- lapply(1:comparisons, function(x){

        contrast <- info$pairwise[x,]

        g1 <- analytic_solve(info$dat[[contrast[[1]]]])
        g2 <- analytic_solve(info$dat[[contrast[[2]]]])

          (g1$pcor_mat - g2$pcor_mat)

      })



      names(z_stat)  <- sapply(1:comparisons,function(x)  paste("Y_g",
                                                              info$pairwise[x,],
                                                              sep = "",
                                                              collapse = " - "))


      returned_object <- list(
        z_stat = z_stat,
        diff = diff,
        p = p,
        info = info,
        iter = iter,
        type = type,
        analytic = analytic,
        call = match.call())
    }
  } else {

    stop("type not supported. ordinal, binary, and mixed data will be implemented soon.")

  }

    class(returned_object) <- c("BGGM",
                              "ggm_compare_estimate",
                              "estimate")
  returned_object

}


#' @name summary.ggm_compare_estimate
#' @title Summary method for \code{ggm_compare_estimate} objects
#'
#' @param object An object of class \code{ggm_compare_estimate}
#' @param cred credible interval width
#' @param ... currently ignored
#' @seealso \code{\link{ggm_compare_estimate}}
#' @return A list containing the summarized posterior distributions
#' @examples
#' # data
#' Y1 <- BGGM::bfi[1:500,1:5]
#' Y2 <- BGGM::bfi[501:1000, 1:5]
#'
#' # fit model
#' fit <- ggm_compare_estimate(Y1, Y2)
#'
#' # posterior summary of differences
#' summary(fit)
#'
#' @export
summary.ggm_compare_estimate <- function(object, cred = 0.95,...) {

  lb <- (1 - cred) / 2
  ub <- 1 - lb
  name_temp <- matrix(0, object$p, object$p)

  name_temp[] <-
    unlist(lapply(1:object$p , function(x)
      paste(1:object$p, x, sep = "--")))

  dat_results <- list()

  for (i in 1:nrow(object$info$pairwise)) {

    ci <- apply(object$pcors_diffs[[i]][[1]], MARGIN = 2,
                FUN = function(x){ quantile(x, probs = c(lb, ub)) })
    diff_mu <-
      apply(object$pcors_diffs[[i]][[1]], MARGIN = 2, mean)

    diff_sd <-
      apply(object$pcors_diffs[[i]][[1]], MARGIN = 2, sd)

    results_temp <-
      data.frame(
        edge = name_temp[upper.tri(name_temp)],
        post_mean =  round(diff_mu, 3),
        post_sd = round(diff_sd, 3),
        ci = round(t(ci), 3)
      )

    colnames(results_temp) <- c(
      "Edge",
      "Post.mean",
      "Post.sd",
      "Cred.lb",
      "Cred.ub"
      )
    dat_results[[i]] <- results_temp
  }
  returned_object <- list(dat_results = dat_results,
                          object = object)
  class(returned_object) <- c("BGGM",
                              "summary",
                              "ggm_compare_estimate",
                               "estimate")
  returned_object
}


