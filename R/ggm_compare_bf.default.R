#' GGM Compare: Exploratory Hypothesis Testing
#'
#' @name ggm_compare_explore
#'
#' @description Compare Gaussian graphical models with exploratory hypothesis testing using the matrix-F prior
#' distribution \insertCite{Mulder2018}{BGGM}. A test for each partial correlation in the model for any number
#' of groups. This provides evidence for the null hypothesis of no difference and the alternative hypothesis
#' of difference. With more than two groups, the test is for \emph{all} groups simultaneously (i.e., the relation
#' is the same or different in all groups). This method was introduced in \insertCite{williams2020comparing;textual}{BGGM}.
#' For confirmatory hypothesis testing see \code{confirm_groups}.
#'
#' @param ... At least two matrices (or data frame) of dimensions \emph{n} (observations) by  \emph{p} (variables).
#'
#' @param formula An object of class \code{\link[stats]{formula}}. This allows for including
#' control variables in the model (i.e., \code{~ gender}).
#'
#' @param prior_sd Numeric. The scale of the prior distribution (centered at zero), in reference to a beta distribtuion.
#' The `default` is 0.25. See note for further details.
#'
#' @param type Character string. Which type of data for \code{Y} ? The options include \code{continuous},
#' \code{binary}, or \code{ordinal}. See the note for further details.
#'
#' @param mixed_type Numeric vector. An indicator of length p for which varibles should be treated as ranks.
#' (1 for rank and 0 to assume normality). The default is currently (dev version) to treat all integer variables
#' as ranks when \code{type = "mixed"} and \code{NULL} otherwise. See note for further details.
#'
#' @param analytic logical. Should the analytic solution be computed (default is \code{FALSE}) ? See note for details.
#'
#' @param iter number of iterations (posterior samples; defaults to 5000).
#'
#' @param progress Logical. Should a progress bar be included (defaults to \code{TRUE}) ?
#'
#' @param seed An integer for the random seed.
#'
#' @references
#' \insertAllCited{}
#'
#'
#' @return The returned object of class \code{ggm_compare_explore} contains a lot of information that
#'         is used for printing and plotting the results. For users of \strong{BGGM}, the following
#'         are the useful objects:
#'
#' \itemize{
#'
#' \item \code{BF_01} A \emph{p} by \emph{p} matrix including
#'                     the Bayes factor for the null hypothesis.
#'
#' \item \code{pcor_diff} A \emph{p} by \emph{p} matrix including
#'                        the difference in partial correlations (only for two groups).
#'
#' \item \code{samp} A list containing the fitted models (of class \code{explore}) for each group.
#'
#' }
#' @details
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
#' @note
#'
#' \strong{"Default" Prior}:
#'
#'  In Bayesian statistics, a default Bayes factor needs to have several properties. I refer
#'  interested users to \insertCite{@section 2.2 in @dablander2020default;textual}{BGGM}. In
#'  \insertCite{Williams2019_bf;textual}{BGGM}, some of these propteries were investigated, such
#'  model selection consistency. That said, we would not consider this a "default" Bayes factor and
#'  thus we encourage users to perform sensitivity analyses by varying the scale of the prior
#'  distribution.
#'
#'  Furthermore, it is important to note there is no "correct" prior and, also, there is no need
#'  to entertain the possibility of a "true" model. Rather, the Bayes factor can be interpreted as
#'  which hypothesis best (relative to each other) predicts the observed data
#'  \insertCite{@Section 3.2 in @Kass1995}{BGGM}.
#'
#' \strong{Interpretation of Conditional (In)dependence Models for Latent Data}:
#'
#' See \code{\link{BGGM-package}} for details about interpreting GGMs based on latent data
#' (i.e, all data types besides \code{"continuous"})
#'
#'
#' @examples
#'
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
#' #############################
#' ### example 1: continuous ###
#' #############################
#'
#' # fit model
#' fit <- ggm_compare_explore(Ymale, Yfemale,
#'                            iter = 250,
#'                            type = "continuous")
#'
#' # summary
#' summary(fit)
#'
#' # plot summary
#' plot(summary(fit))
#'
#' # select graph
#' select(fit)
#'
#' # plot graph
#' plot(select(fit))
#'
#' ##########################
#' ### example 2: ordinal ###
#' ##########################
#'
#' # fit model
#' fit <- ggm_compare_explore(Ymale,  Yfemale,
#'                            type = "ordinal",
#'                            iter = 250)
#'
#' # summary
#' summary(fit)
#'
#' # plot summary
#' plot(summary(fit))
#'
#' # select graph
#' select(fit)
#'
#'
#' #########################
#' ### example 3: mixed  ###
#' #########################
#'
#' # fit model
#' fit <- ggm_compare_explore(Ymale, Yfemale,
#'                            type = "mixed",
#'                            iter = 250)
#'
#' # summary
#' summary(fit)
#'
#' # plot summary
#' plot(summary(fit))
#' }
#'
#' @export
ggm_compare_explore <- function(...,
                           formula = NULL,
                           type = "continuous",
                           mixed_type = NULL,
                           analytic = FALSE,
                           prior_sd = 0.20,
                           iter = 5000,
                           progress = TRUE,
                           seed = 1){

  # combine data
  dat_list <- list(...)

  # combine data
  info <- Y_combine(...)

  # groups
  groups <- length(info$dat)

  delta <- delta_solve(prior_sd)

  # check at least two groups
  if(groups < 2){

    stop("must have (at least) two groups")

    }

  # sample
  if(!analytic){

  samp <- lapply(1:groups, function(x) {

    # mixed
    # message("BGGM: Posterior Sampling ", "(Group ",x ,")")
    Y <- dat_list[[x]]


    # call estimate
    explore(Y, formula = formula,
             type = type,
             prior_sd =  prior_sd,
             iter = iter,
             mixed_type = mixed_type,
             progress = progress,
             seed = x,
             ... = paste0("(Group ", x, ")"))


    })

  post_samp <- lapply(1:groups, function(x) samp[[x]]$post_samp )

  prior_samp <-  lapply(1:groups, function(x) samp[[x]]$prior_samp)

  # p with predictors removed
  p <- samp[[1]]$p


  # store pcor diff
  pcor_diff <- BF_01_mat <- matrix(0, p, p)

  # upper triangular elements
  indices <- which(upper.tri(diag(p)), arr.ind = TRUE )

  # make contrast matrices
  ## words for compatability
  groups_as_words <- numbers2words(1:groups)

  ## hypotheses
  hyp <- paste(groups_as_words, sep = " ", collapse = "=")

  ## `framed` hypotheses
  framed <- framer(hyp)

  ## contrast matrices
  mats <- create_matrices(framed = framed,
                          varnames = groups_as_words)


  # loop through upper triangular
  for(i in seq_len(nrow(indices))){

    rho_ij <- indices[i,]

    # start
    post_group <-  post_samp[[1]]$fisher_z[ rho_ij[1], rho_ij[2], (51:(iter + 50))]
    prior_group <-  prior_samp[[1]]$fisher_z[ 1, 2,]

    # combined groups
    for(j in 2:(groups)){
      post_group <-  cbind(post_group,  post_samp[[j]]$fisher_z[ rho_ij[1], rho_ij[2], (51:(iter + 50))])
      prior_group <-  cbind(prior_group,  prior_samp[[j]]$fisher_z[1, 2,])
    }

    # posterior covariance
    cov_post <- cov(post_group)

    # prior covariance
    cov_prior <- cov(prior_group)

    # posterior mean
    post_mean <- colMeans(post_group)

    # tranformed posterior
    mu_post <- mats$R_e %*% post_mean
    s_post <- mats$R_e %*% cov_post %*% t(mats$R_e)

    # transformed prior
    mu_prior <- mats$R_e %*% rep(0, groups)

    s_prior <- mats$R_e %*% cov_prior %*% t(mats$R_e)

    # bayes factor
    log_BF <- mvnfast::dmvn(X = t(mats$r_e),
                            mu = mu_post,
                            sigma = s_post,
                            log = TRUE) -
      mvnfast::dmvn(X = t(mats$r_e),
                    mu = mu_prior,
                    sigma = s_prior,
                    log = TRUE)

    BF_01_mat[ rho_ij[1], rho_ij[2] ] <- exp(log_BF)

    if(groups == 2){
      pcor_diff[ rho_ij[1], rho_ij[2] ] <-  (z2r(post_mean)[1] - z2r(post_mean)[2])
        }
  }

  BF_01 <-  symmetric_mat(BF_01_mat)

  pcor_diff <- symmetric_mat(pcor_diff)

  returned_object <- list(BF_01 = BF_01,
                          info = info,
                          iter = iter,
                          prior_sd = prior_sd,
                          call = match.call(),
                          delta = delta,
                          groups = groups,
                          pcor_diff = pcor_diff,
                          samp = samp,
                          type = type,
                          p = p)

  # analytic solution
  } else {

    stop("analytic not currently implemented")

  }

  class(returned_object) <- c("BGGM",
                              "ggm_compare_explore",
                              "explore")
  returned_object

  }

print_summary_ggm_compare_bf <- function(x, ...){
  groups <- x$object$groups
  cat("BGGM: Bayesian Gaussian Graphical Models \n")
  cat("--- \n")
  cat("Type:",  x$object$type, "\n")
  cat("Formula:", paste(as.character(x$formula), collapse = " "), "\n")
  # number of iterations
  cat("Posterior Samples:", x$object$iter, "\n")
  # number of observations
  cat("Observations (n):\n")
  groups <- length(x$object$info$dat)
  for(i in 1:groups){
    cat("  Group", paste( i, ":", sep = "") , x$object$info$dat_info$n[[i]], "\n")
  }
  # number of variables
  cat("Variables (p):", x$object$p, "\n")
  # number of edges
  cat("Relations:", .5 * (x$object$p * (x$object$p-1)), "\n")
  cat("Delta:", x$object$delta, "\n")
  cat("--- \n")
  cat("Call: \n")
  print(x$object$call)
  cat("--- \n")
  cat("Hypotheses:\n")
  cat("H0:", paste0("rho_g", 1:groups, collapse = " = "), "\n")
  cat("H1:", paste0("rho_g", 1:groups, collapse = " - "), " = 0\n")
  cat("--- \n\n")

  print(x$results, right = FALSE, row.names = FALSE)
  cat("--- \n")
}

print_ggm_compare_bf <- function(x, ...){
  cat("BGGM: Bayesian Gaussian Graphical Models \n")
  cat("--- \n")
  cat("Type:",  x$type, "\n")
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
  cat("Variables (p):", x$p, "\n")
  # number of edges
  cat("Relations:", .5 * (x$p * (x$p-1)), "\n")
  cat("Delta:", x$delta, "\n")
  cat("--- \n")
  cat("Call: \n")
  print(x$call)
  cat("--- \n")
  cat("Hypotheses:\n")
  cat("H0:", paste0("rho_g", 1:groups, collapse = " = "), "\n")
  cat("H1:", paste0("rho_g", 1:groups, collapse = " - "), " = 0\n")
  cat("--- \n")
  cat("Date:", date(), "\n")
}

#' @title Summary Method for \code{ggm_compare_explore} Objects
#'
#' @description Summarize the posterior hypothesis probabilities
#'
#' @name summary.ggm_compare_explore
#'
#' @param object An object of class \code{ggm_compare_explore}.
#'
#' @param col_names Logical. Should the summary include the column names (default is \code{TRUE})?
#'                  Setting to \code{FALSE} includes the column numbers (e.g., \code{1--2}).
#'
#' @param ... Currently ignored.
#'
#' @return An object of class \code{summary.ggm_compare_explore}
#'
#' @seealso \code{\link{ggm_compare_explore}}
#'
#' @export
summary.ggm_compare_explore <- function(object,
                                        col_names = TRUE,
                                        ...){

  # nodes
  p <- object$p

  # identity matrix
  I_p <- diag(p)

  # prob null
  prob_H0 <-   round(object$BF_01 / (object$BF_01 + 1), 3)

  # prob h1
  prob_H1 <-   round(1 - prob_H0, 3)

  # column names
  cn <-  colnames(object$samp[[1]]$Y)


  if(is.null(cn)){

    mat_names <- sapply(1:p , function(x) paste(1:p, x, sep = "--"))[upper.tri(I_p)]

  } else {


    mat_names <-  sapply(cn , function(x) paste(cn, x, sep = "--"))[upper.tri(I_p)]

  }

  if(object$groups == 2){

    post_mean <- round(object$pcor_diff[upper.tri(I_p)], 3)

    post_sd <- round(apply(object$samp[[1]]$post_samp$pcors -
                          object$samp[[2]]$post_samp$pcors, 1:2, sd)[upper.tri(I_p)], 3)

    results <- data.frame(Relation = mat_names,
                                Post.mean = post_mean,
                                Post.sd = post_sd,
                                Pr.H0 = prob_H0[upper.tri(I_p)],
                                Pr.H1 = prob_H1[upper.tri(I_p)])

    } else {

    results <- data.frame(Relation = mat_names,
                          Pr.H0 = prob_H0[upper.tri(I_p)],
                          Pr.H1 = prob_H1[upper.tri(I_p)])

    }

  returned_object <- list(results = results,
                          object = object)

  class(returned_object) <- c("BGGM",
                              "ggm_compare_explore",
                              "summary.ggm_compare_explore",
                              "explore")
  returned_object
}


#' @title Plot \code{summary.ggm_compare_explore} Objects
#'
#' @description Visualize the posterior hypothesis probabilities.
#'
#' @name plot.summary.ggm_compare_explore
#'
#' @param x An object of class \code{summary.ggm_compare_explore}
#'
#' @param size Numeric. The size of the points (defaults to 2).
#'
#' @param color Character string. The color of the points
#' (defaults to \code{"black"}).
#'
#' @param ... Currently ignored.
#'
#' @return A \code{ggplot} object
#'
#' @seealso \code{\link{ggm_compare_explore}}
#'
#' @export
plot.summary.ggm_compare_explore <- function(x,
                                             size = 2,
                                             color = "black", ...){


  dat_temp <- x$results[order(x$results$Pr.H1,
                              decreasing = F), ]

  dat_temp$Relation <-
    factor(dat_temp$Relation,
           levels = dat_temp$Relation,
           labels = dat_temp$Relation)


  ggplot(dat_temp,
         aes(x = Relation,
             y = Pr.H1)) +
    geom_point(size = size, color = color) +

    theme(axis.text.x = element_text(
      angle = 90,
      vjust = 0.5,
      hjust = 1
    )) +
    coord_flip()

}


