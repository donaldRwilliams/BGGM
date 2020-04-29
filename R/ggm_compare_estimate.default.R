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
#' @param seed The seed for random number generation (default set to \code{1}).
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
                                 iter = 5000){

  # set seed
  # set.seed(seed)

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

          # data Y_gj
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

          # posterior sample
          post_samp <- lapply(1:groups, function(x) {

            # data Y_gj
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

      # binary data
      } else if(type == "binary") {

        # intercept only
        if(is.null(formula)){

          X <- matrix(1, n, 1)

          # model matrix
          } else {

            X <- model.matrix(formula, data)

            }

        # posterior sample
        post_samp <- lapply(1:groups, function(x) {

          # Y_gj
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

        # posterior sample
        post_samp <- lapply(1:groups, function(x) {

          # Y_gj
          Y <- as.matrix(info$dat[[x]])

          # categories
          K <- max(apply(Y, 2, function(x) { length(unique(x))   } ))

          # call c ++
          .Call("_BGGM_mv_ordinal_albert",
              Y = Y + 1,
              X = X,
              iter = iter + 50,
              delta = delta,
              epsilon = 0.1,
              K = K)$pcors
          })

        # mixed data
        } else if(type == "mixed"){

          # no control variables allowed
          if(!is.null(formula)){

            warning("formula ignored for mixed data at this time")
            formula <- NULL
            }

          # posterior samples
          post_samp <- lapply(1:groups, function(x) {

            # Y_gj
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
#'
#' @title Summary method for \code{ggm_compare_estimate} objects
#'
#' @param object an object of class \code{ggm_compare_estimate}
#'
#' @param col_names logical. Should the summary include the column names (default is \code{TRUE})?
#'                  Setting to \code{FALSE} includes the column numbers (e.g., \code{1--2}).
#'
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
summary.ggm_compare_estimate <- function(object,
                                         col_names = TRUE,
                                         cred = 0.95,...) {

  # nodes
  p <- object$p

  # identity matrix
  I_p <- diag(p)

  # lower bound
  lb <- (1 - cred) / 2

  # upper bound
  ub <- 1 - lb

  # relation names
  name_mat <- matrix(0, p, p)

  # number of comparisons
  comparisons <- length(names(object$diff))

  # column names
  cn <-  colnames(object$info$dat[[1]])


  if(col_names | is.null(cn)){

    mat_names <-  sapply(cn , function(x) paste(cn, x, sep = "--"))[upper.tri(I_p)]

    } else {

      mat_names <- sapply(1:p , function(x) paste(1:p, x, sep = "--"))[upper.tri(I_p)]

      }

  dat_results <- list()

  # summary for comparison i
  for(i in seq_len(comparisons)){

    post_mean <- round(apply( object$diff[[i]], 1:2, mean), 3)[upper.tri(I_p)]

    post_sd  <- round(apply( object$diff[[i]], 1:2, sd), 3)[upper.tri(I_p)]

    post_lb <- round(apply( object$diff[[i]], 1:2, quantile, lb), 3)[upper.tri(I_p)]

    post_ub <- round(apply( object$diff[[i]], 1:2, quantile, ub), 3)[upper.tri(I_p)]


    results_i <-
      data.frame(
        relation = mat_names,
        post_mean =  post_mean,
        post_sd = post_sd,
        post_lb = post_lb,
        post_ub = post_ub
      )

    colnames(results_i) <- c(
      "Relation",
      "Post.mean",
      "Post.sd",
      "Cred.lb",
      "Cred.ub"
    )


    dat_results[[i]] <- results_i
  }

  returned_object <- list(dat_results = dat_results,
                          object = object)
  class(returned_object) <- c("BGGM",
                              "summary", "summary.ggm_compare_estimate",
                              "ggm_compare_estimate",
                              "estimate")
  returned_object
}

# print ggm compare
print_ggm_compare <- function(x, ...){
  cat("BGGM: Bayesian Gaussian Graphical Models \n")
  cat("--- \n")
  cat("Type:",  x$type, "\n")
  cat("Analytic:", x$analytic, "\n")
  cat("Formula:", paste(as.character(fit$formula), collapse = " "), "\n")
  # number of iterations
  cat("Posterior Samples:", x$iter, "\n")
  # number of observations
  cat("Observations (n):\n")
  groups <- length(x$info$dat)
  for(i in 1:groups){
    cat("  Group", paste( i, ":", sep = "") , x$info$dat_info$n[[i]], "\n")
  }
  # number of variables
  cat("Nodes (p):", x$p, "\n")
  # number of edges
  cat("Relations:", .5 * (x$p * (x$p-1)), "\n")
  cat("--- \n")
  cat("Call: \n")
  print(x$call)
  cat("--- \n")
  cat("Date:", date(), "\n")
}

# print summary
print_summary_ggm_estimate_compare <- function(x,...){

  cat("BGGM: Bayesian Gaussian Graphical Models \n")
  cat("--- \n")
  cat("Type:",  x$object$type, "\n")
  cat("Analytic:", x$object$analytic, "\n")
  cat("Formula:", paste(as.character(fit$formula), collapse = " "), "\n")
 # number of iterations
  cat("Posterior Samples:", x$object$iter, "\n")
  # number of observations
  cat("Observations (n):\n")
  groups <- length(x$object$info$dat)
  for (i in 1:groups) {
    cat("  Group",
        paste(i, ":", sep = "") ,
        x$object$info$dat_info$n[[i]],
        "\n")
  }
  # number of variables
  cat("Nodes (p):", x$object$p, "\n")
  # number of edges
  cat("Relations:", .5 * (x$object$p * (x$object$p - 1)), "\n")
  cat("--- \n")
  cat("Call: \n")
  print(x$object$call)
  cat("--- \n")
  cat("Estimates:\n")
  for (i in 1:nrow(x$object$info$pairwise)) {
    cat("\n", names(x$object$pcors_diffs[[i]]), "\n")

    print(x$dat_results[[i]], right = FALSE, row.names = FALSE,...)

  }
  cat("--- \n")
}





#' Plot \code{summary.ggm_compare_estimate} Objects
#'
#' @param x an object of class \code{estimate} or \code{ggm_compare_estimate}
#' @param color color of error bar
#' @param width width of error bar cap
#' @param ... currently ignored
#'
#' @return an object of class \code{ggplot}
#' @export
plot.summary.ggm_compare_estimate <- function(x, color = "black",
                                  size = 2,
                                  width = 0, ...){

  n_plt  <- length(x$dat_results)

  # plots
  lapply(1:seq_len(n_plt), function(i){

    dat_temp <- x$dat_results[[i]][order(x$dat_results[[i]]$Post.mean,
                                         decreasing = F), ]

    dat_temp$Relation <-
      factor(dat_temp$Relation,
             levels = dat_temp$Relation,
             labels = dat_temp$Relation)


    ggplot(dat_temp,
           aes(x = Relation,
               y = Post.mean)) +

      geom_errorbar(aes(ymax = dat_temp[, 4],
                        ymin = dat_temp[, 5]),
                    width = width,
                    color = color) +
      geom_point(size = size) +
      xlab("Index") +
      theme(axis.text.x = element_text(
        angle = 90,
        vjust = 0.5,
        hjust = 1
      )) +
      ggtitle(paste(names(x$object$diff)))
  })
}

