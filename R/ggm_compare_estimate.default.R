#' GGM Compare: Estimate
#'
#' @name ggm_compare_estimate
#'
#' @description Compare edges (partial correlations) that are estimated from groups to, say, detect a differences or equivalence.
#'
#' @param ... matrices (or data frame) of dimensions \emph{n} (observations) by  \emph{p} (variables).
#' Requires at least two.
#'
#' @param formula an object of class \code{\link[stats]{formula}}. This allows for including
#' control variables in the model (i.e., \code{~ gender}).
#'
#' @param data an optional data frame, list or environment (or an object coercible by \code{\link[base]{as.data.frame}})
#' to a data frame containing the variables in \code{formula}. This is required when controlling for variables.
#'
#' @param prior_sd The scale of the prior distribution (centered at zero), in reference to a beta distribtuion.
#' The `default` is 0.25. See note for further details.
#'
#' @param type character string. Which type of data for \strong{Y} ? The options include \code{continuous},
#' \code{binary}, or \code{ordinal}. See the note for further details.
#'
#' @param mixed_type numeric vector. An indicator of length p for which varibles should be treated as ranks.
#' (1 for rank and 0 to assume normality). The default is currently (dev version) to treat all integer variables
#' as ranks when \code{type = "mixed"} and \code{NULL} otherwise. See note for further details.
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
                                 formula = NULL,
                                 data = NULL,
                                 type = "continuous",
                                 mixed_type = NULL,
                                 analytic = FALSE,
                                 prior_sd = 0.25,
                                 iter = 5000,
                                 seed = 1){


  # set seed
  set.seed(seed)

  # combine data
  info <- Y_combine(...)

  # number of variables
  p <- info$dat_info$p[1]

  # number of observation
  n = info$dat_info$n[1]

  # number of groups
  groups <- length(info$dat)

  # number of comparisons
  comparisons <- nrow(info$pairwise)

  # delta rho ~ beta(delta/2, delta/2)
  delta <- delta_solve(prior_sd)

  if(groups < 2){
    stop("must have (at least) two groups")
    }

  # sample
  if(!analytic){
    # continuous
    if(type == "continuous"){
      # no control variables
      if(is.null(formula)){

        # posterior sample
        post_samp <- lapply(1:groups, function(x) {

        Y <- as.matrix(scale(info$dat[[x]], scale = F))

        .Call(
          '_BGGM_Theta_continuous',
          PACKAGE = 'BGGM',
          Y = Y,
          iter = iter + 50,
          delta = delta,
          epsilon = 0.1,
          prior_only = 0,
          explore = 1
        )$pcors
        })

        # control for variables
        } else {

          # model matrix
          X <- model.matrix(formula, data)

          post_samp <- lapply(1:groups, function(x) {

            Y <- as.matrix(scale(info$dat[[x]], scale = F))

            .Call(
              "_BGGM_mv_continuous",
              Y = Y,
              X = X,
              delta = delta,
              epsilon = 0.1,
              iter = iter + 50
            )$pcors

          })

          }

      } else if(type == "binary") {


        if(is.null(formula)){

          X <- matrix(1, n, 1)

          } else {

            # model matrix
            X <- model.matrix(formula, data)

            }

        post_samp <- lapply(1:groups, function(x) {

          Y <- as.matrix(info$dat[[x]])

          .Call(
            "_BGGM_mv_binary",
            Y = Y,
            X = X,
            delta = delta,
            epsilon = 0.1,
            iter = iter + 50,
            beta_prior = 0.0001,
            cutpoints = c(-Inf, 0, Inf)
          )$pcors

        })

      } else if(type == "ordinal"){

        # intercept only
        if(is.null(formula)){

          X <- matrix(1, n, 1)

        } else {

          # model matrix
          X <- model.matrix(formula, data)

        }

        post_samp <- lapply(1:groups, function(x) {

          Y <- as.matrix(info$dat[[x]])

          K <- max(apply(Y, 2, function(x) { length(unique(x))   } ))

          .Call("_BGGM_mv_ordinal_albert",
              Y = Y + 1,
              X = X,
              iter = iter + 50,
              delta = delta,
              epsilon = 0.1,
              K = K)$pcors
        })

      } else if(type == "mixed"){


        if(!is.null(formula)){

          warning("formula ignored for mixed data at this time")

          formula <- NULL

          }



        post_samp <- lapply(1:groups, function(x) {

          Y <- as.matrix(info$dat[[x]])

          # default for ranks
          if(is.null(mixed_type)) {

            idx = colMeans(round(Y) == Y)
            idx = ifelse(idx == 1, 1, 0)

            # user defined
          } else {

            idx = mixed_type

          }


          # rank following hoff (2008)
          rank_vars <- rank_helper(Y)

          .Call("_BGGM_copula",
                z0_start = rank_vars$z0_start,
                levels = rank_vars$levels,
                K = rank_vars$K,
                Sigma_start = rank_vars$Sigma_start,
                iter = iter + 50,
                delta = delta,
                epsilon = 0.1,
                idx = idx)$pcors
          })

        } else {
          stop("'type' not supported: must be continuous, binary, ordinal, or mixed.")
          }




      diff <- lapply(1:comparisons, function(x) {
        contrast <- info$pairwise[x, ]

        post_samp[[contrast[[1]]]][, , 51:(iter + 50)] - post_samp[[contrast[[2]]]][, , 51:(iter + 50)]

      })


      names(diff)  <- sapply(1:comparisons, function(x)
        paste("Y_g",
              info$pairwise[x, ],
              sep = "",
              collapse = " - "))




      # returned object
      returned_object <- list(
        diff = diff,
        p = p,
        info = info,
        iter = iter,
        analytic = analytic,
        type = type,
        formula = formula,
        call = match.call()
      )



  # analytic
    } else {

      if(type != "continuous"){
        warning("analytic solution only available for 'type = continuous'")
        type <- "continuous"
      }

      formula <- NULL

      z_stat <- lapply(1:comparisons, function(x) {
        contrast <- info$pairwise[x, ]

        g1 <- analytic_solve(info$dat[[contrast[[1]]]])
        g2 <- analytic_solve(info$dat[[contrast[[2]]]])

        z_stat <-
          abs((g1$inv_map - g2$inv_map) /   sqrt(g1$inv_var + g2$inv_var))

      })


      diff <- lapply(1:comparisons, function(x) {
        contrast <- info$pairwise[x, ]

        g1 <- analytic_solve(info$dat[[contrast[[1]]]])
        g2 <- analytic_solve(info$dat[[contrast[[2]]]])

        (g1$pcor_mat - g2$pcor_mat)

      })

      names(diff)  <- sapply(1:comparisons, function(x)
        paste("Y_g",
              info$pairwise[x, ],
              sep = "",
              collapse = " - "))

      names(z_stat)  <-
        sapply(1:comparisons, function(x)
          paste("Y_g",
                info$pairwise[x, ],
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
        call = match.call()
      )

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





