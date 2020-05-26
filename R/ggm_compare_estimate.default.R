#' GGM Compare: Estimate
#'
#' @name ggm_compare_estimate
#'
#' @description Compare partial correlations that are estimated from any number of groups. This method works for
#' continuous, binary, ordinal, and mixed data (a combination of categorical and continuous variables).
#' The approach (i.e., a difference between posterior distributions) was
#' described in  \insertCite{Williams2019;textual}{BGGM}.
#'
#' @param ... Matrices (or data frames) of dimensions \emph{n} (observations) by  \emph{p} (variables).
#' Requires at least two.
#'
#' @param formula An object of class \code{\link[stats]{formula}}. This allows for including
#' control variables in the model (i.e., \code{~ gender}). See the note for further details.
#'
#'
#' @param prior_sd The scale of the prior distribution (centered at zero), in reference to a beta distribtuion
#'                 (defaults to 0.50).
#'                 See note for further details.
#'
#' @param type Character string. Which type of data for \strong{Y} ? The options include \code{continuous},
#' \code{binary}, \code{ordinal}, or \code{continuous}. See the note for further details.
#'
#' @param mixed_type Numeric vector. An indicator of length \emph{p} for which varibles should be treated as ranks.
#' (1 for rank and 0 to use the 'empirical' or observed distribution). The default is currently to treat all integer variables
#' as ranks when \code{type = "mixed"} and \code{NULL} otherwise. See note for further details.
#'
#' @param iter Number of iterations (posterior samples; defaults to 5000).
#'
#' @param analytic Logical. Should the analytic solution be computed (default is \code{FALSE})? This is only available
#'                 for continous data. Note that if \code{type = "mixed"} and \code{analytic = TRUE}, the data will
#'                 automatically be treated as continuous.
#'
#' @param progress Logical. Should a progress bar be included (defaults to \code{TRUE}) ?
#'
#' @param seed An integer for the random seed.
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
#' @details
#' This function can be used to compare the partial correlations for any number of groups.
#' This is accomplished with pairwise comparisons for each relation. In the case of three groups,
#' for example, group 1 and group 2 are compared, then group 1 and group 3 are compared, and then
#' group 2 and group 3 are compared. There is a full distibution for each difference that can be
#' summarized (i.e., \code{\link{summary.ggm_compare_estimate}}) and then visualized
#' (i.e., \code{\link{plot.summary.ggm_compare_estimate}}). The graph of difference is selected with
#' \code{\link{select.ggm_compare_estimate}}).
#'
#'
#' \strong{Controlling for Variables}:
#'
#' When controlling for variables, it is assumed that \code{Y} includes \emph{only}
#' the nodes in the GGM and the control variables. Internally, \code{only} the predictors
#' that are included in \code{formula} are removed from \code{Y}. This is not behavior of, say,
#' \code{\link{lm}}, but was adopted to ensure  users do not have to write out each variable that
#' should be included in the GGM. An example is provided below.
#'
#' \strong{Mixed Type}:
#'
#'  The term "mixed" is somewhat of a misnomer, because the method can be used for data including \emph{only}
#'  continuous or \emph{only} discrete variables. This is based on the ranked likelihood which requires sampling
#'  the ranks for each variable (i.e., the data is not merely transformed to ranks). This is computationally
#'  expensive when there are many levels. For example, with continuous data, there are as many ranks
#'  as data points!
#'
#'  The option \code{mixed_type} allows the user to determine  which variable should be treated as ranks
#'  and the "emprical" distribution is used otherwise. This is accomplished by specifying an indicator
#'  vector of length \emph{p}. A one indicates to use the ranks, whereas a zero indicates to "ignore"
#'  that variable. By default all integer variables are handled as ranks.
#'
#' \strong{Dealing with Errors}:
#'
#' An error is most likely to arise when \code{type = "ordinal"}. The are two common errors (although still rare):
#'
#' \itemize{
#'
#' \item The first is due to sampling the thresholds, especially when the data is heavily skewed.
#'       This can result in an ill-defined matrix. If this occurs, we recommend to first try
#'       decreasing \code{prior_sd} (i.e., a more informative prior). If that does not work, then
#'       change the data type to \code{type = mixed} which then estimates a copula GGM
#'       (this method can be used for data containing \strong{only} ordinal variable). This should
#'       work without a problem.
#'
#' \item  The second is due to how the ordinal data are categorized. For example, if the error states
#'        that the index is out of bounds, this indicates that the first category is a zero. This is not allowed, as
#'        the first category must be one. This is addressed by adding one (e.g., \code{Y + 1}) to the data matrix.
#'
#' }
#'
#'
#' @note
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
#'
#' \strong{Interpretation of Conditional (In)dependence Models for Latent Data}:
#'
#' See \code{\link{BGGM-package}} for details about interpreting GGMs based on latent data
#' (i.e, all data types besides \code{"continuous"})
#'
#'
#'
#' \strong{Additional GGM Compare Methods}
#'
#' Bayesian hypothesis testing is implemented in \code{\link{ggm_compare_explore}} and
#' \code{\link{ggm_compare_confirm}} \insertCite{Williams2019_bf}{BGGM}. The latter allows for confirmatory
#' hypothesis testing.  An approach based on a posterior predictive check is implemented in \code{\link{ggm_compare_ppc}}
#' \insertCite{williams2020comparing}{BGGM}. This provides  a 'global' test for comparing the entire GGM and a 'nodewise'
#' test for comparing each variable in the network \insertCite{Williams2019;textual}{BGGM}.
#'
#' @examples
#' \donttest{
#' # note: iter = 250 for demonstrative purposes
#'
#' # data
#' Y <- bfi
#'
#' # males and females
#' Ymale <- subset(Y, gender == 1,
#'                    select = -c(gender,
#'                                education))[,1:10]
#'
#' Yfemale <- subset(Y, gender == 2,
#'                      select = -c(gender,
#'                                  education))[,1:10]
#'
#' # fit model
#' fit <- ggm_compare_estimate(Ymale,  Yfemale,
#'                            type = "ordinal",
#'                            iter = 250,
#'                            prior_sd = 0.25,
#'                            progress = FALSE)
#'
#' ###########################
#' ### example 2: analytic ###
#' ###########################
#' # only continuous
#'
#' # fit model
#' fit <- ggm_compare_estimate(Ymale, Yfemale,
#'                             analytic = TRUE)
#'
#' # summary
#' summ <- summary(fit)
#'
#' # plot summary
#' plt_summ <- plot(summary(fit))
#'
#' # select
#' E <- select(fit)
#'
#' # plot select
#' plt_E <- plot(select(fit))
#'
#' }
#'
#' @export
ggm_compare_estimate <- function(...,
                                 formula = NULL,
                                 type = "continuous",
                                 mixed_type = NULL,
                                 analytic = FALSE,
                                 prior_sd = 0.50,
                                 iter = 5000,
                                 progress = TRUE,
                                 seed = 1){
  # combine data
  dat_list <- list(...)

  # combine data
  info <- Y_combine(...)

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
               seed = x,
               progress = progress,
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

      p <- ncol(diff[[1]])

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


#' @title Summary method for \code{ggm_compare_estimate} objects
#'
#' @description Summarize the posterior distribution of each partial correlation
#' difference with the posterior mean and standard deviation.
#'
#' @name summary.ggm_compare_estimate
#'
#' @param object An object of class \code{ggm_compare_estimate}.
#'
#' @param col_names Logical. Should the summary include the column names (default is \code{TRUE})?
#'                  Setting to \code{FALSE} includes the column numbers (e.g., \code{1--2}).
#'
#' @param cred Numeric. The credible interval width for summarizing the posterior
#'             distributions (defaults to 0.95; must be between 0 and 1).
#'
#' @param ... Currently ignored.
#'
#' @seealso \code{\link{ggm_compare_estimate}}
#'
#' @return A list containing the summarized posterior distributions.
#'
#' @examples
#' \donttest{
#' # note: iter = 250 for demonstrative purposes
#' # data
#' Y <- bfi
#'
#' # males and females
#' Ymale <- subset(Y, gender == 1,
#'                 select = -c(gender,
#'                             education))[,1:5]
#'
#' Yfemale <- subset(Y, gender == 2,
#'                   select = -c(gender,
#'                               education))[,1:5]
#'
#' # fit model
#' fit <- ggm_compare_estimate(Ymale,  Yfemale,
#'                             type = "ordinal",
#'                             iter = 250,
#'                             prior_sd = 0.25,
#'                             progress = FALSE)
#'
#' summary(fit)
#' }
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

   if(isFALSE(object$analytic )){

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

   } else {

     post_mean <- round(object$diff[[i]][upper.tri(I_p)], 3)

     results_i <-
       data.frame(
         relation = mat_names,
         post_mean =  post_mean
       )


     colnames(results_i) <- c(
       "Relation",
       "Post.mean"
     )

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
  cat("Formula:", paste(as.character(x$formula), collapse = " "), "\n")
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
  cat("Formula:", paste(as.character(x$object$formula), collapse = " "), "\n")
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
#' @description Visualize the posterior distribution differences.
#'
#' @param x An object of class \code{ggm_compare_estimate}.
#'
#' @param size Numeric. The size of the points (defaults to 2).
#'
#' @param color Character string. The color of the points
#' (defaults to \code{"black"}).
#'
#' @param width Numeric. The width of error bar ends (defaults to \code{0}).
#'
#' @param ... Currently ignored.
#'
#' @return An object of class \code{ggplot}
#'
#' @seealso \code{\link{ggm_compare_estimate}}
#'
#' @examples
#' \donttest{
#' # note: iter = 250 for demonstrative purposes
#' # data
#' Y <- bfi
#'
#' # males and females
#' Ymale <- subset(Y, gender == 1,
#'                 select = -c(gender,
#'                             education))[,1:5]
#'
#' Yfemale <- subset(Y, gender == 2,
#'                   select = -c(gender,
#'                               education))[,1:5]
#'
#' # fit model
#' fit <- ggm_compare_estimate(Ymale,  Yfemale,
#'                             type = "ordinal",
#'                             iter = 250,
#'                             prior_sd = 0.25,
#'                             progress = FALSE)
#'
#' plot(summary(fit))
#' }
#'
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

