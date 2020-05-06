#' GGM Compare: Estimate
#'
#' @name ggm_compare_estimate
#'
#' @description Compare partial correlations that are estimated from any number of groups. This method works for
#' continuous, binary, ordinal, and mixed data (a combination of categorical and continuous variables).
#' The approach (i.e., a difference between posterior distributions) was
#' described in  \insertCite{Williams2019;textual}{BGGM}.
#'
#' @param ... matrices (or data frame) of dimensions \emph{n} (observations) by  \emph{p} (variables).
#' Requires at least two.
#'
#' @param formula an object of class \code{\link[stats]{formula}}. This allows for including
#' control variables in the model (i.e., \code{~ gender}). See the note for further details.
#'
#'
#' @param prior_sd The scale of the prior distribution (centered at zero), in reference to a beta distribtuion.
#' The `default` is 0.50. See note for further details.
#'
#' @param type character string. Which type of data for \strong{Y} ? The options include \code{continuous},
#' \code{binary}, \code{ordinal}, or \code{continuous}. See the note for further details.
#'
#' @param mixed_type numeric vector. An indicator of length \emph{p} for which varibles should be treated as ranks.
#' (1 for rank and 0 to use the 'empirical' or observed distribution). The default is currently to treat all integer variables
#' as ranks when \code{type = "mixed"} and \code{NULL} otherwise. See note for further details.
#'
#' @param iter number of iterations (posterior samples; defaults to 5000).
#'
#' @param analytic logical. Should the analytic solution be computed (default is \code{FALSE})? This is only available
#'                 for continous data. Note that if \code{type = "mixed"} and \code{analytic = TRUE}, the data will
#'                 automatically be treated as continuous.
#'
#' @references
#' \insertAllCited{}
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
#' @note This function can be used to compare the partial correlations for any number of groups.
#' This is accomplished with pairwise comparisons for each relation. In the case of three groups,
#' for example, group 1 and group 2 are compared, then group 1 and group 3 are compared, and then
#' group 2 and group 3 are compared. There is a full distibution for each difference that can be
#' summarized (i.e., \code{\link{summary.ggm_compare_estimate}}) and then visualized
#' (i.e., \code{\link{plot.summary.ggm_compare_estimate}}).
#'
#' \strong{Selecting the Graph of Differences:}
#'
#' Selecting the differences is implemented in \code{\link{select.ggm_compare_estimate}}. Currently the only available
#' option is to detect differences with  credible interval exclusion of zero. This
#' corresponds to a directional posterior probability \insertCite{Marsman2017a}{BGGM}. For example,
#' the probability (conditional on the model) of a difference  is at least 97.5 $\%$ when the 95 $\%$
#' credible interval excludes zero.
#'
#'
#' \strong{Interpretation of conditional (in)dependence models for latent data:}
#'
#' A  tetrachoric correlation (binary data) is a special case of a polychoric correlation (ordinal data).
#' Both relations are between "theorized normally distributed continuous latent variables"
#' (\href{https://en.wikipedia.org/wiki/Polychoric_correlation}{Wikipedia}). In both instances,
#' the correpsonding partial correlation between observed variables is conditioned
#' on the remaining variables in the \emph{latent} space. This implies that interpration
#' is much the same as for continuous data, but with respect to latent variables.
#' We refer interested reader to \insertCite{@page 2364, section 2.2, in  @webb2008bayesian;textual}{BGGM}.
#'
#'
#' \strong{Mixed Data:}
#'
#' The mixed data approach was introduced  \insertCite{@in @hoff2007extending;textual}{BGGM}
#' (our paper describing an extension to Bayesian hypothesis testing if forthcoming).
#' This is a semi-paramateric copula model based on the ranked likelihood. This is computationally
#' expensive when treating continuous data as ranks. The current default is to treat only integer data as ranks.
#' This should of course be adjusted for continous data that is skewed. This can be accomplished with the
#' argument \code{mixed_type}. A \code{1} in the numeric vector of length \emph{p}indicates to treat that
#' respective node as a rank (corresponding to the column number) and a zero indicates to use the observed
#' (or "emprical") data.
#'
#'
#' It is also important to note that \code{type = "mixed"} is not restricted to mixed data (containing a combination of
#' categorical and continuous): all the nodes can be ordinal or continuous (but again this will take some time).
#'
#' \strong{Prior Distribution:}
#'
#'  We use the novel matrix-F prior distribution that was recently introduces in X. In the context of GGMs, this prior
#'  was first used in X. The key advantage is that the implied prior distrbution for each partial correlation is (1)
#'  invariant to the network size (see section X in X) and (2) approximately beta distributed (Equation X in X). Accordingly,
#'  the arguement \code{prior_sd} corresponds the standard deviation of a beta distrubtion (with a mean of zero). The function
#'  \code{\link{plot_prior}} can be used to visualize the prior distibution.
#'
#' \strong{Additional GGM Compare Methods}
#'
#' Bayesian hypothesis testing is implemented in \code{\link{ggm_compare_explore}} and
#' \code{\link{ggm_compare_confirm}} \insertCite{Williams2019_bf}{BGGM}. The latter allows for confirmatory
#' hypothesis testing.  An approach based on a posterior predictive check is implemented in \code{\link{ggm_compare_ppc}}
#' \insertCite{williams2020comparing}{BGGM}. This provides  a 'global' test for comparing the entire GGM and a 'nodewise'
#' test for comparing each variable in the network. Testing for nodewise differences in predictabilty is implemented in
#' \code{\link{test.R2}} \insertCite{Williams2019;textual}{BGGM}.
#'
#' @examples
#' # comparing two groups (males vs. females) in a personality network
#'
#'
#'
#'
#' @export
ggm_compare_estimate <- function(...,
                                 formula = NULL,
                                 type = "continuous",
                                 mixed_type = NULL,
                                 analytic = FALSE,
                                 prior_sd = 0.50,
                                 iter = 5000){
  # combine data
  dat_list <- list(...)

  # combine data
  info <- Y_combine(...)

  # # number of variables
  # p <- info$dat_info$p[1]
  #
  # # number of observation
  # n = info$dat_info$n[1]

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

    post_samp <- lapply(1:groups, function(x) {

      Y <- dat_list[[x]]

      # call estimate
      estimate(Y, formula = formula,
               type = type,
               prior_sd = prior_sd,
               iter = iter,
               mixed_type = mixed_type,
               ... = paste0("(Group ", x, ")"))

    })

    # number of variables
    p <- post_samp[[1]]$p

    # compute difference
    diff <- lapply(1:comparisons, function(x) {

      contrast <- info$pairwise[x, ]
      post_samp[[contrast[[1]]]]$post_samp$pcors[, , 51:(iter + 50)] - post_samp[[contrast[[2]]]]$post_samp$pcors[, , 51:(iter + 50)]

      })

    # name posterior (differences) array
    names(diff)  <- sapply(1:comparisons, function(x)
        paste("Y_g",
              info$pairwise[x, ],
              sep = "",
              collapse = " - "))

    # pcor_mats
    pcor_mats <- lapply(1:length(diff), function(x) {
       round(apply(diff[[x]], 1:2, mean), 3)

       })

      # name pcor_mats
      names(pcor_mats) <- names(diff)

      # returned object
      returned_object <- list(
        pcor_mats = pcor_mats,
        diff = diff,
        p = p,
        info = info,
        iter = iter,
        analytic = analytic,
        type = type,
        formula = formula,
        call = match.call(),
        post_samp= post_samp
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
  cn <-  colnames(object$post_samp[[1]]$Y)


  if(is.null(cn) | isFALSE(col_names)){

    mat_names <- sapply(1:p , function(x) paste(1:p, x, sep = "--"))[upper.tri(I_p)]

    } else {

      mat_names <-  sapply(cn , function(x) paste(cn, x, sep = "--"))[upper.tri(I_p)]

      }

  dat_results <- list()

  # summary for comparison i
  for(i in seq_len(comparisons)){

    post_mean <- round(apply(object$diff[[i]], 1:2, mean), 3)[upper.tri(I_p)]

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

    if(isTRUE( object$analytic)){

      results_i <- results_i[,1:2]
    }

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
  lapply(1:n_plt, function(i){

    dat_temp <- x$dat_results[[i]][order(x$dat_results[[i]]$Post.mean,
                                         decreasing = F), ]

    dat_temp$Relation <-
      factor(dat_temp$Relation,
             levels = dat_temp$Relation,
             labels = dat_temp$Relation)


   plt <- ggplot(dat_temp,
           aes(x = Relation,
               y = Post.mean)) +
      geom_point(size = size) +
      xlab("Index") +
      theme(axis.text.x = element_text(
        angle = 90,
        vjust = 0.5,
        hjust = 1
      )) +
      ggtitle(paste(names(x$object$diff)[i]))

      if(isFALSE( x$object$analytic)){

      plt <- plt +  geom_errorbar(aes(ymax = dat_temp[, 4],
                          ymin = dat_temp[, 5]),
                      width = width,
                      color = color)
      }

      plt
  })
}

