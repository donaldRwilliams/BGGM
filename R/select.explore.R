#' Select Graphical Structure with the Bayes Factor
#' @name select.explore
#' @description This allows for not only estimating the conditional dependence structure, that is non-zero edges, but also the conditional \strong{in}dependence
#' structure (evidence for no relation).
#'
#' @param object object of class \code{explore.default}
#' @param BF_cut evidentiary threshold
#' @param alternative type of hypothesis (see notes)
#' @param ... currently not used
#'
#' @note The \code{alternative} can be either \code{greater}, \code{less}, \code{two.sided}, or \code{exhaustive}. The argument \code{hyp_prob} is used only
#' when \code{alternative = hypothesis}. \code{greater} and \code{less} test directional hypotheses, and thus, the graphical structure will
#' only included edges in that direction (i.e., positive or negative). \code{two.sided} is the customary approach, and test for the presence or
#' absence of an edge. \code{exhaustive} tests negative vs. positive  vs. zero. Here \code{hyp_prob} is the posterior probability threshold for the respective
#' hypotheses.
#'
#'
#' @return list of class \code{select.explore}:
#'
#' \code{alternative = "two.sided"}:
#' \itemize{
#'  \item \code{partials_non_zero} selected partial correlation matrix
#'  \item \code{pcor_mat} partial correlation matrix (non set to zero)
#'  \item \code{pcor_sd} partial correlation standard deviations
#'  \item \code{Adj_10} adjacency matrix for the selected edges (in favor of the \code{alternative})
#'  \item \code{Adj_01} adjacency matrix for the null hypothesis  (conditional independence)
#'  \item \code{BF_10} Bayes factors for \code{alternative}
#'  \item \code{BF_01} Bayes factors for the null hypothesis
#'  \item \code{BF_cut} evidentiary threshold
#'  \item \code{alternative} \code{"two.sided"}
#'  \item \code{code} \code{match.call()}
#' }
#'
#' \code{alternative = "greater"}:
#' \itemize{
#'  \item \code{partials_positive} selected partial correlation matrix
#'  \item \code{pcor_mat} partial correlation matrix (none set to zero)
#'  \item \code{pcor_sd} partial correlation standard deviations
#'  \item \code{Adj_20} adjacency matrix for the selected edges (in favor of the \code{alternative})
#'  \item \code{Adj_01} adjacency matrix for the null hypothesis  (conditional independence)
#'  \item \code{BF_20} Bayes factors for \code{alternative}
#'  \item \code{BF_01} Bayes factors for the null hypothesis
#'  \item \code{BF_cut} evidentiary threshold
#'  \item \code{alternative} \code{"greater"}
#'  \item \code{code} \code{match.call()}
#' }
#'
#' \code{alternative = "less"}:
#' \itemize{
#'  \item \code{partials_negative} selected partial correlation matrix
#'  \item \code{pcor_mat} partial correlation matrix (none set to zero)
#'  \item \code{pcor_sd} partial correlation standard deviations
#'  \item \code{Adj_20} adjacency matrix for the selected edges (in favor of the \code{alternative})
#'  \item \code{Adj_01} adjacency matrix for the null hypothesis  (conditional independence)
#'  \item \code{BF_20} Bayes factors for \code{alternative}
#'  \item \code{BF_01} Bayes factors for the null hypothesis
#'  \item \code{BF_cut} evidentiary threshold
#'  \item \code{alternative} \code{"less"}
#'  \item \code{code} \code{match.call()}
#' }
#'
#' \code{alternative = "exhaustive"}
#' \itemize{
#' \item \code{post_prob} data.frame with posterior probabilities for each edge
#' \item \code{neg_mat} adjacency matrix for negative edges
#' \item \code{post_mat} adjacency matrix for positive edges
#' \item \code{null_mat} adjacency matrix for zero  (conditional independence)
#' \item \code{"alternative"} "exhaustive"
#' \item \code{pcor_mat} partial correlation matrix (non set to zero)
#' \item \code{pcor_sd} partial correlation standard deviations
#' \item \code{code} \code{match.call()}
#' \item \code{prob} \code{hyp_prob}
#' }
#'
#' @examples
#' \donttest{
#' # p = 10
#' Y <- BGGM::bfi[,1:10]
#'
#' # sample posterior
#' fit <- explore(Y, iter = 5000)
#'
#' # select E
#' E <- select(fit, BF_cut = 3)
#'
#' # summarize
#' summary(E)
#'
#' # non-zero edges
#' E$partials_non_zero
#'
#' # adjacency matrix
#' E$Adj_10
#'
#' # null adjacency matrix
#' E$Adj_01
#' }
#' @export
select.explore <- function(object,
                           BF_cut = 3,
                           alternative = "two.sided",
                           ...){
  # rename
  x <- object

  # hyp probability
  hyp_prob <- BF_cut / (BF_cut + 1)




    # posterior samples
    post_samp <- x$post_samp

    # prior samples
    prior_samp <- x$prior_samp



  # two sided testing
  if(alternative == "two.sided"){

    # posterior
    post_sd <- apply(post_samp$fisher_z[,,(51:x$iter)], 1:2, sd)
    post_mean  <- x$pcor_mat
    post_dens <- dnorm(0, post_mean, post_sd )

    # prior
    prior_sd <- apply(prior_samp$fisher_z[,,(51:x$iter)], 1:2, sd)
    prior_dens <- dnorm(0, 0, mean(prior_sd))

    # BF
    BF_10_mat <- prior_dens / post_dens
    BF_01_mat <- 1 / BF_10_mat
    diag(BF_01_mat) <- 0
    diag(BF_10_mat) <- 0

    # selected: alternative
    Adj_10 <- ifelse(BF_10_mat > BF_cut, 1, 0)

    # selected: null
    Adj_01 <- ifelse(BF_10_mat < 1 / BF_cut, 1, 0)
    diag(Adj_01) <- 0
    diag(Adj_10) <- 0

    # returned object
    returned_object = list(pcor_mat_zero = post_mean * Adj_10,
                           pcor_mat = round(post_mean, 3),
                           pcor_sd = round(post_sd, 3),
                           Adj_10 = Adj_10,
                           Adj_01 = Adj_01,
                           BF_10 = BF_10_mat,
                           BF_01 = BF_01_mat,
                           BF_cut = BF_cut,
                           alternative = alternative,
                           call = match.call(),
                           type = x$type,
                           formula = x$formula,
                           analytic = x$analytic,
                           object = object
                           )

    # one sided greater
    } else if(alternative == "greater"){

      # posterior
      post_sd <- apply(post_samp$fisher_z[,,(51:x$iter)], 1:2, sd)
      post_mean  <- x$pcor_mat
      post_dens <- dnorm(0, post_mean, post_sd )

      # prior
      prior_sd <- apply(prior_samp$fisher_z[,,(51:x$iter)], 1:2, sd)
      prior_dens <- dnorm(0, 0, mean(prior_sd))

      # BF (two sided)
      BF_10_mat <- prior_dens / post_dens

      # BF one sided
      BF_20_mat <-  BF_10_mat * ((1 - pnorm(0, post_mean, post_sd)) * 2)

      # BF null
      BF_02_mat <- 1 / BF_20_mat
      diag(BF_02_mat) <- 0
      diag(BF_20_mat) <- 0

      # selected edges (alternative)
      Adj_20 <- ifelse(BF_20_mat > BF_cut, 1, 0)

      # selected edges (null)
      Adj_02 <- ifelse(BF_02_mat > BF_cut, 1, 0)
      diag(Adj_02) <- 0
      diag(Adj_20) <- 0

      # returned object
      returned_object = list(
        pcor_mat_zero = post_mean * Adj_20,
        pcor_mat = round(post_mean, 3),
        pcor_sd = round(post_sd, 3),
        Adj_20 = Adj_20,
        Adj_02 = Adj_02,
        BF_20 = BF_20_mat,
        BF_02 = BF_02_mat,
        BF_cut = BF_cut,
        alternative = alternative,
        call = match.call(),
        type = x$type,
        formula = x$formula,
        analytic = x$analytic,
        object = object
      )

    # one side less
    } else if(alternative == "less")  {

      # posterior
      post_sd <- apply(post_samp$fisher_z[,,(51:x$iter)], 1:2, sd)
      post_mean  <- x$pcor_mat
      post_dens <- dnorm(0, post_mean, post_sd )

      # prior
      prior_sd <- apply(prior_samp$fisher_z[,,(51:x$iter)], 1:2, sd)
      prior_dens <- dnorm(0, 0, mean(prior_sd))

      # BF (two sided)
      BF_10_mat <- prior_dens / post_dens

      # BF one sided
      BF_20_mat <- BF_10_mat * (pnorm(0, post_mean, post_sd) * 2)

      # BF null
      BF_02_mat <- 1 / BF_20_mat
      diag(BF_02_mat) <- 0
      diag(BF_20_mat) <- 0

      # selected edges (alternative)
      Adj_20 <- ifelse(BF_20_mat > BF_cut, 1, 0)

      # selected edges (null)
      Adj_02 <- ifelse(BF_02_mat > BF_cut, 1, 0)
      diag(Adj_02) <- 0
      diag(Adj_20) <- 0

      # returned object
      returned_object = list(
        pcor_mat_zero = post_mean * Adj_20,
        pcor_mat = round(post_mean, 3),
        pcor_sd = round(post_sd, 3),
        Adj_20 = Adj_20,
        Adj_02 = Adj_02,
        BF_20 = BF_20_mat,
        BF_02 = BF_02_mat,
        BF_cut = BF_cut,
        alternative = alternative,
        call = match.call(),
        type = x$type,
        formula = x$formula,
        analytic = x$analytic,
        object = object
      )

      # exhaustive testing
      } else if (alternative == "exhaustive")

        if(alternative == "exhaustive"){

          if(is.null(hyp_prob)){

            stop("posterior probability must be specificed \n for exhaustive hypothesis testing")

            }



          # column names
          cn <-  colnames(x$Y)

          p  <- ncol(x$pcor_mat)

          I_p <- diag(p)

          if(is.null(cn)){

            mat_names <- sapply(1:p , function(z) paste(1:p, z, sep = "--"))[upper.tri(I_p)]

          } else {

            mat_names <-  sapply(cn , function(z) paste(cn, z, sep = "--"))[upper.tri(I_p)]

          }

          # posterior
          post_sd <- apply(post_samp$fisher_z[,,(51:x$iter)], 1:2, sd)
          post_mean  <- x$pcor_mat
          post_dens <- dnorm(0, post_mean, post_sd)

          # prior
          prior_sd <- apply(prior_samp$fisher_z[,,(51:x$iter)], 1:2, sd)
          prior_dens <- dnorm(0, 0, mean(prior_sd))

          # BF (two sided)
          BF_10_mat <- prior_dens / post_dens

          # BF less
          BF_less <- BF_10_mat  * (pnorm(0, post_mean, post_sd) * 2)
          BF_greater <-  BF_10_mat * ((1 - pnorm(0, post_mean, post_sd)) * 2)

          # BF null
          BF_null <- 1 / BF_10_mat

          prob_null <-  BF_null / (BF_null + BF_greater + BF_less)
          prob_greater <-  BF_greater / (BF_null + BF_greater + BF_less)
          prob_less <-  BF_less / (BF_null + BF_greater + BF_less)

          prob_mat <-  prob_null + prob_greater + prob_less

          prob_dat = data.frame(edge = mat_names,
                                prob_zero = prob_null[upper.tri(prob_null)],
                                prob_greater = prob_greater[upper.tri(prob_greater)],
                                prob_less = prob_less[upper.tri(prob_less)])

          # no rownames
          row.names(prob_dat) <- c()

          # selected  (null)
          null_mat <- ifelse(prob_null > hyp_prob, 1, 0)

          # selected (positive)
          pos_mat <-  ifelse(prob_greater > hyp_prob, 1, 0)


          # selected  (negative)
          neg_mat <-  ifelse(prob_less > hyp_prob, 1, 0)

          # negative edges
          returned_object <- list(
            post_prob = prob_dat,
            neg_mat = neg_mat,
            pos_mat = pos_mat,
            null_mat = null_mat,
            alternative = alternative,
            pcor_mat = round(post_mean, 3),
            pcor_sd = round(post_sd, 3),
            call = match.call(),
            prob = hyp_prob,
            type = x$type,
            formula = x$formula,
            analytic = x$analytic,
            object = object
          )
        } else {

          stop("alternative not supported. see documentation")
    }

  class(returned_object) <- c("BGGM",
                              "select.explore",
                              "explore",
                              "select")
  returned_object
}




print_select_explore <- function(x,
                                 ...){

  p <- ncol(x$pcor_mat_zero)
  cat("BGGM: Bayesian Gaussian Graphical Models \n")
  cat("--- \n")
  cat("Type:", x$type, "\n")
  cat("Analytic:", x$analytic, "\n")
  cat("Formula:", paste(as.character(fit$formula), collapse = " "), "\n")
  cat("Alternative:", x$alternative, "\n")
  if(x$alternative == "two.sided"){
    cat("Bayes Factor:", x$BF_cut, "\n")

  }
  cat("--- \n")
  cat("Call:\n")
  print(x$call)
  cat("--- \n")
  cat("Hypotheses: \n")

  if(x$alternative == "two.sided"){

    cat("H0: rho = 0\nH1: rho != 0", "\n")
    cat("--- \n")
    colnames(x$Adj_10) <- 1:p
    row.names(x$Adj_10) <- 1:p
    colnames( x$pcor_mat_zero) <- 1:p
    row.names(x$pcor_mat_zero) <- 1:p
    cat("Partial Correlations:\n\n")
    print(round(x$pcor_mat_zero, 2))
    cat("--- \n")
    cat("Adjacency:\n\n")
    print(x$Adj_10)
    cat("--- \n")
  } else if (x$alternative == "greater"){

    cat("H0: rho = 0\nH1: rho > 0", "\n")
    cat("--- \n")
    colnames(x$Adj_20) <- 1:p
    row.names(x$Adj_20) <- 1:p
    colnames( x$pcor_mat_zero) <- 1:p
    row.names(x$pcor_mat_zero) <- 1:p
    cat("Partial Correlations:\n\n")
    print(round(x$pcor_mat_zero, 2))
    cat("--- \n")
    cat("Adjacency:\n\n")
    print(x$Adj_20)
    cat("--- \n")

  } else if (x$alternative == "less"){

    cat("H0: rho = 0\nH1: rho < 0", "\n")
    cat("--- \n")
    colnames(x$Adj_20) <- 1:p
    row.names(x$Adj_20) <- 1:p
    colnames( x$pcor_mat_zero) <- 1:p
    row.names(x$pcor_mat_zero) <- 1:p
    cat("Partial Correlations:\n\n")
    print(round(x$pcor_mat_zero, 2))
    cat("--- \n")
    cat("Adjacency:\n\n")
    print(x$Adj_20)
    cat("--- \n")
  } else {

    cat("H0: rho = 0\nH1: rho > 0\nH2: rho < 0", "\n")
    cat("--- \n")
    cat("Summary:\n\n")
    dat <- x$post_prob
    dat$prob_zero <- round(dat$prob_zero, 3)
    dat$prob_greater <- round(dat$prob_greater, 3)
    dat$prob_less <- round(dat$prob_less, 3)
    colnames(dat) <- c("Relation", "Pr.H0", "Pr.H1", "Pr.H2")
    print(dat, row.names = FALSE, right = FALSE)
    cat("--- \n")
  }
}



#' Summary Method for \code{select.explore} Objects
#'
#' @param object object of class \code{select.explore}.
#'
#' @return a data frame including the posterior mean, standard deviation,
#' and posterior hypothesis probabilities for each relation.
#' @export
summary.select.explore <- function(object, col_names = TRUE){

  x <- object

  p <- ncol(x$pcor_mat)

  I_p <- diag(p)

  # column names
  cn <-  colnames(object$object$Y)


  if(!isTRUE(col_names) | is.null(cn)){

    mat_names <- sapply(1:p , function(x) paste(1:p, x, sep = "--"))[upper.tri(I_p)]

  } else {


    mat_names <-  sapply(cn , function(x) paste(cn, x, sep = "--"))[upper.tri(I_p)]

  }



  if(x$alternative == "two.sided"){

    post_mean <- x$pcor_mat[upper.tri(x$pcor_mat)]
    post_sd <-  x$pcor_sd[upper.tri(x$pcor_sd)]
    prob_H1 <- x$BF_10[upper.tri(x$BF_10)] / (x$BF_10[upper.tri(x$BF_10)] + 1)
    prob_H0 <- 1 - prob_H1
    summ <-  data.frame(
      Relation = mat_names,
      Post.mean = post_mean,
      Post.sd = post_sd,
      Pr.H0 = round(prob_H0, 3),
      Pr.H1 = round(prob_H1, 3)
    )

  } else if (x$alternative == "greater"){

    post_mean <- x$pcor_mat[upper.tri(x$pcor_mat)]
    post_sd <-  x$pcor_sd[upper.tri(x$pcor_sd)]
    prob_H1 <- x$BF_20[upper.tri(x$BF_20)] / (x$BF_20[upper.tri(x$BF_20)] + 1)
    prob_H0 <- 1 - prob_H1
    summ <-  data.frame(
      Relation = mat_names,
      Post.mean = post_mean,
      Post.sd = post_sd,
      Pr.H0 = round(prob_H0, 3),
      Pr.H1 = round(prob_H1, 3)
    )



  } else if (x$alternative == "less" | x$alternative == "greater"){

    post_mean <- x$pcor_mat[upper.tri(x$pcor_mat)]
    post_sd <-  x$pcor_sd[upper.tri(x$pcor_sd)]
    prob_H1 <- x$BF_20[upper.tri(x$BF_20)] / (x$BF_20[upper.tri(x$BF_20)] + 1)
    prob_H0 <- 1 - prob_H1
    summ <-  data.frame(
      Relation = mat_names[upper.tri(mat_names)],
      Post.mean = post_mean,
      Post.sd = post_sd,
      Pr.H0 = round(prob_H0, 3),
      Pr.H1 = round(prob_H1, 3)
    )



  } else {

    summ <- cbind.data.frame( x$post_prob[,1],
                              x$pcor_mat[upper.tri(x$pcor_mat)],
                              x$pcor_sd[upper.tri(x$pcor_sd)],
                              round(x$post_prob[,2:4], 3))

    colnames(summ) <- c("Relation",
                        "Post.mean",
                        "Post.sd",
                        "Pr.H0",
                        "Pr.H1",
                        "Pr.H2")


  }

  returned_object <- list(summary = summ, object = object)

  class(returned_object) <- c("BGGM", "summary.select.explore",
                              "explore", "select.explore",
                              "summary")
  returned_object


}



print_summary_select_explore <- function(x,...){

  cat("BGGM: Bayesian Gaussian Graphical Models \n")
  cat("--- \n")
  cat("Type:", x$object$type, "\n")
  cat("Alternative:", x$object$alternative, "\n")
  cat("--- \n")
  cat("Call:\n")
  print(x$object$call)
  cat("--- \n")
  cat("Hypotheses: \n")

  if(x$object$alternative == "two.sided"){

    cat("H0: rho = 0\nH1: rho != 0", "\n")

  } else if (x$object$alternative == "greater"){

    cat("H0: rho = 0\nH1: rho > 0", "\n")

  } else if (x$object$alternative == "less"){

    cat("H0: rho = 0\nH1: rho < 0", "\n")

  } else {

    cat("H0: rho = 0\nH1: rho > 0\nH2: rho < 0", "\n")

  }

  cat("--- \n\n")

  print(x$summary, right = FALSE, row.names = FALSE)


}


#' GGM: Plot \code{summary.select.explore] Objects}
#'
#' @description Visualize posterior hypothesis probabilities.
#'
#' @param x object of class \code{summary.select.explore}
#' @param size numeric. Size of the points.
#' @param color string. Color of the points.
#'
#' @return
#' @export
plot.summary.select.explore <- function(x, size = 2, color = "black"){


  dat_temp <- x$summary[order(x$summary$Pr.H1,
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
    ))

}
