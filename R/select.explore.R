#' Select Graphical Structure with the Bayes Factor
#' @name select.explore
#' @description This allows for not only estimating the conditional dependence structure, that is non-zero edges, but also the conditional \strong{in}dependence
#' structure (evidence for no relation).
#'
#' @param object object of class \code{explore.default}
#' @param BF_cut evidentiary threshold
#' @param alternative type of hypothesis (see notes)
#' @param hyp_prob posterior probability threshold (see notes)
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
#' @export
select.explore <- function(object,
                           BF_cut = 3,
                           alternative = "two.sided",
                           hyp_prob = NULL, ...){
  # rename
  x <- object

  # posterior samples
  posterior_samples <- do.call(rbind.data.frame,
                               lapply(1:x$cores, function(z)
                                 x$samples[[z]]$fisher_z_post))

  # prior samples
  prior_samples <- unlist(do.call(rbind.data.frame,
                                  lapply(1:x$cores, function(z)
                                  x$samples[[z]]$fisher_z_prior)))

  if(alternative == "two.sided"){

    # matrices for storage
    BF_10_mat <- BF_01_mat <- matrix(0, x$p, x$p)

    # prior density
    prior_dens <- dnorm(0, mean(prior_samples), sd(prior_samples))

    # BF for alternative
    BF_10 <- apply(posterior_samples, MARGIN = 2,
                   FUN = function(z){ prior_dens / dnorm(0, mean(z), sd(z)) })

    # BF for null
    BF_01 <- 1 / BF_10

    # alternative BF mat
    BF_10_mat[upper.tri(BF_10_mat)] <- BF_10[1:x$edge]
    BF_10_mat <- symmteric_mat(BF_10_mat)

    # null BF mat
    BF_01_mat <- 1 / BF_10_mat
    diag(BF_01_mat) <- 0

    # selected edges (alternative)
    Adj_10 <- ifelse(BF_10_mat > BF_cut, 1, 0)

    # selected edges (null)
    Adj_01 <- ifelse(BF_10_mat < 1 / BF_cut, 1, 0)
    diag(Adj_01) <- 0

    # returned object
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

    # prior density
    prior_dens <- dnorm(0, mean(prior_samples), sd(prior_samples))

    # matrices for storage
    BF_20_mat <- BF_01_mat <- matrix(0, x$p, x$p)

    # BF for alternative
    BF_10 <- apply(posterior_samples, MARGIN = 2,
                   FUN = function(z){ prior_dens / dnorm(0, mean(z), sd(z)) })

    # density greater than zero
    dens_greater <-  apply(posterior_samples, MARGIN = 2,
                           FUN = function(z){ (1  - pnorm(0, mean(z), sd(z))) * 2 })

    # one sided BF
    BF_20 <- dens_greater * (BF_10)

    # BF null
    BF_01 <- 1 / BF_10

    # alternative BF mat
    BF_20_mat[upper.tri(BF_20_mat)] <- BF_20[1:x$edge]
    BF_20_mat <- symmteric_mat(BF_20_mat)

    # null BF mat
    BF_01_mat[upper.tri(BF_01_mat)] <- BF_01[1:x$edge]
    BF_01_mat <- symmteric_mat(BF_01_mat)
    diag(BF_01_mat) <- 0

    # selected edges (alternative)
    Adj_20 <- ifelse(BF_20_mat > BF_cut, 1, 0)

    # selected edges (null)
    Adj_01 <- ifelse(BF_01_mat > BF_cut, 1, 0)
    diag(Adj_01) <- 0

    # returned object
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

    # prior density
    prior_dens <- dnorm(0, mean(prior_samples), sd(prior_samples))

    # matrices for storage
    BF_20_mat <- BF_01_mat <- matrix(0, x$p, x$p)

    # BF for alternative
    BF_10 <- apply(posterior_samples, MARGIN = 2,
                   FUN = function(z){ prior_dens / dnorm(0, mean(z), sd(z)) })

    # density less than zero
    dens_less <-  apply(posterior_samples, MARGIN = 2,
                        FUN = function(z){ (pnorm(0, mean(z), sd(z))) * 2 })

    # one sided BF
    BF_20 <- dens_less * (BF_10)

    # BF null
    BF_01 <- 1 / BF_10

    # alternative BF mat
    BF_20_mat[upper.tri(BF_20_mat)] <- BF_20[1:x$edge]
    BF_20_mat <- symmteric_mat(BF_20_mat)

    # null BF mat
    BF_01_mat[upper.tri(BF_01_mat)] <- BF_01[1:x$edge]
    BF_01_mat <- symmteric_mat(BF_01_mat)
    diag(BF_01_mat) <- 0

    # selected edges (alternative)
    Adj_20 <- ifelse(BF_20_mat > BF_cut, 1, 0)

    # selected edges (null)
    Adj_01 <- ifelse(BF_01_mat > BF_cut, 1, 0)
    diag(Adj_01) <- 0

    # returned object
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

    # matrices for storage
    mat_names <- null_mat <- pos_mat <- neg_mat <- matrix(0, x$p, x$p)

    # matrix names
    mat_names[] <-  unlist(lapply(1:x$p,
                                  FUN = function(z) paste(1:x$p, z, sep = "--")))

    # prior density
    prior_dens <- dnorm(0, mean(prior_samples), sd(prior_samples))

    # two sided BF
    BF_10 <- apply(posterior_samples, MARGIN = 2,
                   FUN = function(z){ prior_dens / dnorm(0, mean(z), sd(z)) })

    # density less than zero
    dens_less <-  apply(posterior_samples, MARGIN = 2,
                        FUN = function(z){ (pnorm(0, mean(z), sd(z))) * 2 })

    # BF less than zero
    BF_less <- dens_less * (BF_10)

    # BF null
    BF_01 <- 1 / BF_10

    # density greater than zero
    dens_greater <-  apply(posterior_samples, MARGIN = 2,
                           FUN = function(z){ (1  - pnorm(0, mean(z), sd(z))) * 2 })

    # BF greater than zero
    BF_greater <- dens_greater * (BF_10)

    # combine BFs
    BF_mat <- data.frame(BF_01 = BF_01,
                         BF_greater = BF_greater,
                         BF_less = BF_less)

    # posterior hypothesis probabilities
    post_prob <-  data.frame(t(apply(BF_mat, 1, FUN = function(x) { x / sum(x) })))

    # column names
    colnames(post_prob) <- c("prob_zero", "prob_greater", "prob_less")

    # no rownames
    row.names(post_prob) <- c()

    # round
    post_prob <-  round(post_prob, 3)

    # combine probabilities and edge names
    post_prob <- cbind.data.frame(edge = mat_names[upper.tri(mat_names)], post_prob)[1:x$edge,]

    # selected edges (null)
    null_mat[upper.tri(null_mat)] <- ifelse(post_prob$prob_zero > hyp_prob, post_prob$prob_zero, 0)
    null_mat <- symmteric_mat(null_mat)

    # selected edges (positive)
    pos_mat[upper.tri(pos_mat)] <-  ifelse(post_prob$prob_greater >  hyp_prob, post_prob$prob_greater, 0)
    pos_mat <- symmteric_mat(pos_mat)

    # selected edges (negative)
    neg_mat[upper.tri(neg_mat)]  <-  ifelse(post_prob$prob_less > hyp_prob, post_prob$prob_less , 0)
    neg_mat <- symmteric_mat(neg_mat)

    # negative edges
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

#' @name print.select.explore
#' @title  Print method for \code{select.explore} objects
#'
#' @param x An object of class \code{select.explore}
#' @param summarize summarize connectivity
#' @param ... currently ignored
#' @seealso \code{\link{select.explore}}
#' @export
print.select.explore <- function(x, summarize = FALSE, ...){
  # name of package
  cat("BGGM: Bayesian Gaussian Graphical Models \n")
  cat("--- \n")
  # hypothesis testing
  cat("Type: Hypothesis Testing \n")
  # alternative argument
  cat("Alternative:", x$alternative, "\n")
  # if exhaustive hypothesis testing
  if(x$alternative == "exhaustive"){
    # posterior probability
    cat("Posterior probability:", x$prob, "\n")
    cat("--- \n")
    cat("Call:\n")
    print(x$call)
    cat("--- \n")
  }
  if(x$alternative == "greater"){
    cat("Bayes Factor:", x$BF_cut, "\n")
    if(isFALSE(summarize)){
      cat("Connectivity:", round(mean(x$Adj_20[upper.tri(x$Adj_20)]) * 100, 1), "% \n")
    }
    cat("--- \n")
    cat("Call:\n")
    print(x$call)
    cat("--- \n")
  }
  # if less than zero
  if(x$alternative == "less"){
    cat("Bayes Factor:", x$BF_cut, "\n")
    if(isFALSE(summarize)){
      cat("Connectivity:", round(mean(x$Adj_20[upper.tri(x$Adj_20)]) * 100, 1), "% \n")
    }
    cat("--- \n")
    cat("Call:\n")
    print(x$call)
    cat("--- \n")
  }
  if (x$alternative == "two.sided") {
    if(isFALSE(summarize)){
      cat("Connectivity:", round(mean(x$Adj_10[upper.tri(x$Adj_10)]) * 100, 1), "% \n")
    }
    cat("--- \n")
    cat("Call:\n")
    print(x$call)
    cat("--- \n")
  }
}

#' @name summary.select.explore
#' @title Summary method for \code{select.explore} objects
#'
#' @param object An object of class \code{select.explore}
#' @param hyp hypothesis to summarize (default is the alternative H1)
#' @param log log scale for the Bayes factors (default is \code{TRUE})
#' @param summarize if \code{TRUE} partial correlations and credible intervals are provided
#' @param ... currently ignored
#' @seealso \code{\link{select.explore}}
#' @export
summary.select.explore <- function(object, hyp = "H1",
                                   log = TRUE, summarize = FALSE, ...){
  x <- object
  cat("BGGM: Bayesian Gaussian Graphical Models \n")
  cat("--- \n")
  cat("Type: Hypothesis Testing \n")
  cat("Alternative:", x$alternative, "\n")

  # exhaustive
  if(x$alternative == "exhaustive"){
    cat("Posterior probability:", x$prob, "\n")
    cat("--- \n")
    cat("Call:\n")
    print(x$call)
    cat("--- \n")
    if(!isFALSE(summarize)){

      summ <-  cbind.data.frame(edge = x$post_prob$edge,
                                post_mean = x$pcor_mat[upper.tri(x$pcor_mat)],
                                post_sd = x$pcor_sd[upper.tri(x$pcor_sd)],
                                "p(H0|Y)" = x$post_prob[,2],
                                "p(H1|Y)" = x$post_prob[,3],
                                "p(H2|Y)" = x$post_prob[,4])

      cat("Hypotheses: \n")
      cat("H0: rho = 0\nH1: rho > 0\nH2: rho < 0", "\n")
      cat("--- \n")
      cat("Estimates: \n \n ")
      print(summ, row.names = FALSE, ...)
      cat("--- \n")
    } else{
      if(hyp == "H1"){
        cat("Hypothesis: \n")
        cat("H1: rho > 0 \n")
        cat("--- \n")
        p <- ncol(x$pos_mat)
        colnames(x$pos_mat) <- 1:p
        row.names(x$pos_mat) <- 1:p
        cat("Partial Correlations \n \n")
        print(x$pcor_mat * ifelse(x$pos_mat > x$prob, 1, 0), ...)
        cat("--- \n \n")
        cat("Adjancency \n \n")
        print(ifelse(x$pos_mat > x$prob, 1, 0), ...)
        cat("--- \n")

      }

      if(hyp == "H0"){
        cat("Hypothesis: \n")
        cat("H0: rho = 0 \n")
        cat("--- \n")
        p <- ncol(x$null_mat)
        colnames(x$null_mat) <- 1:p
        row.names(x$null_mat) <- 1:p
        cat("Partial Correlations \n \n")
        print(x$pcor_mat * ifelse(x$null_mat > x$prob, 1, 0), ...)
        cat("--- \n \n")
        cat("Adjancency \n \n")
        print(ifelse(x$null_mat > x$prob, 1, 0), ...)
        cat("--- \n")

      }

      if(hyp == "H2"){
        cat("Hypothesis: \n")
        cat("H0: rho < 0 \n")
        cat("--- \n")
        p <- ncol(x$neg_mat)
        colnames(x$neg_mat) <- 1:p
        row.names(x$neg_mat) <- 1:p
        cat("Partial Correlations \n \n")
        print(x$pcor_mat * ifelse(x$neg_mat > x$prob, 1, 0), ...)
        cat("--- \n \n")
        cat("Adjancency \n \n")
        print(ifelse(x$neg_mat > x$prob, 1, 0), ...)
        cat("--- \n")
      }
    }
  }

  # greater than
  if(x$alternative == "greater"){
    cat("Bayes Factor:", x$BF_cut, "\n")
    if(isFALSE(summarize)){
      cat("Connectivity:", round(mean(x$Adj_20[upper.tri(x$Adj_20)]) * 100, 1), "% \n")

    }
    cat("--- \n")
    cat("Call:\n")
    print(x$call)
    cat("--- \n")

    if(!isFALSE(summarize)){
      p <- ncol(x$pcor_mat)
      mat_names <- matrix(0, p, p)

      edge_name <-

        mat_names[] <-  unlist(lapply(1:p, function(z) paste(1:p, z, sep = "--")))

      if(log == TRUE){
        summ <-  cbind.data.frame(edge = mat_names[upper.tri(mat_names)],
                                  post_mean = x$pcor_mat[upper.tri(x$pcor_mat)],
                                  post_sd = x$pcor_sd[upper.tri(x$pcor_sd)],
                                  "BF 20" = log(x$BF_20[upper.tri(x$BF_20)]),
                                  "BF 01" = log(x$BF_01[upper.tri(x$BF_01)]))
      } else{

        summ <-  cbind.data.frame(edge = mat_names[upper.tri(mat_names)],
                                  post_mean = x$pcor_mat[upper.tri(x$pcor_mat)],
                                  post_sd = x$pcor_sd[upper.tri(x$pcor_sd)],
                                  "BF 20" = x$BF_20[upper.tri(x$BF_20)],
                                  "BF 01" = x$BF_01[upper.tri(x$BF_01)])
      }

      cat("Hypotheses: \n")
      cat("H0: rho = 0\nH1: rho > 0", "\n")
      cat("--- \n")
      cat("Estimates: \n \n ")
      print(summ, row.names = FALSE, ...)
      cat("--- \n")
      cat("note: BF 20 is a one-sided Bayes factor for H1 \n")
      cat("--- \n")
    } else{
      if(hyp == "H1"){

        cat("Hypothesis: \n")
        cat("H1: rho > 0 \n")
        cat("--- \n")
        p <- ncol(x$partials_positive)
        colnames(x$partials_positive) <- 1:p
        row.names(x$partials_positive) <- 1:p
        cat("Partial Correlations \n \n")
        print(x$partials_positive, ...)
        cat("--- \n \n")
        cat("Adjancency (positive) \n \n")
        colnames(x$Adj_20) <- 1:p
        row.names(x$Adj_20) <- 1:p
        print(x$Adj_20, ...)
        cat("--- \n")
      }


      if(hyp == "H0"){

        cat("Hypothesis: \n")
        cat("H0: rho = 0 \n")
        cat("--- \n")
        p <- ncol(x$pcor_mat)
        colnames(x$pcor_mat) <- 1:p
        row.names(x$pcor_mat) <- 1:p
        cat("Partial Correlations \n \n")
        print(x$pcor_mat * x$Adj_01, ...)
        cat("--- \n \n")
        cat("Adjancency (null) \n \n")
        colnames(x$Adj_01) <- 1:p
        row.names(x$Adj_01) <- 1:p
        print(x$Adj_01, ...)
        cat("--- \n")
      }
    }
}
  # less than
  if(x$alternative == "less"){
    cat("Bayes Factor:", x$BF_cut, "\n")
    if(isFALSE(summarize)){
      cat("Connectivity:", round(mean(x$Adj_20[upper.tri(x$Adj_20)]) * 100, 1), "% \n")
    }
    cat("--- \n")
    cat("Call:\n")
    print(x$call)
    cat("--- \n")

    if(!isFALSE(summarize)){
      p <- ncol(x$pcor_mat)
      mat_names <- matrix(0, p, p)

      # edge_name <- unlist(lapply(1:p, function(z) paste(1:p, z, sep = "--")))

        mat_names[] <-  unlist(lapply(1:p, function(z) paste(1:p, z, sep = "--")))

      if(log == TRUE){
        summ <-  cbind.data.frame(edge = mat_names[upper.tri(mat_names)],
                                  post_mean = x$pcor_mat[upper.tri(x$pcor_mat)],
                                  post_sd = x$pcor_sd[upper.tri(x$pcor_sd)],
                                  "BF 20" = log(x$BF_20[upper.tri(x$BF_20)]),
                                  "BF 01" = log(x$BF_01[upper.tri(x$BF_01)]))
      } else{

        summ <-  cbind.data.frame(edge = mat_names[upper.tri(mat_names)],
                                  post_mean = x$pcor_mat[upper.tri(x$pcor_mat)],
                                  post_sd = x$pcor_sd[upper.tri(x$pcor_sd)],
                                  "BF 20" = x$BF_20[upper.tri(x$BF_20)],
                                  "BF 01" = x$BF_01[upper.tri(x$BF_01)])
      }

      cat("Hypotheses: \n")
      cat("H0: rho = 0\nH1: rho < 0", "\n")
      cat("--- \n")
      cat("Estimates: \n \n ")
      print(summ, row.names = FALSE, ...)
      cat("--- \n")
      cat("note: BF 20 is a one-sided Bayes factor for H1 \n")
      cat("--- \n")
    }

    else{
      if(hyp == "H1"){
        cat("Hypothesis: \n")
        cat("H1: rho < 0 \n")
        cat("--- \n")
        p <- ncol(x$partials_negative)
        colnames(x$partials_negative) <- 1:p
        row.names(x$partials_negative) <- 1:p
        cat("Partial Correlations \n \n")
        print(x$partials_negative, ...)
        cat("--- \n \n")
        cat("Adjancency (negative) \n \n")
        colnames(x$Adj_20) <- 1:p
        row.names(x$Adj_20) <- 1:p
        print(x$Adj_20, ...)
        cat("--- \n")
      }

      if(hyp == "H0"){
        cat("Hypothesis: \n")
        cat("H0: rho = 0 \n")
        cat("--- \n")
        p <- ncol(x$pcor_mat)
        colnames(x$pcor_mat) <- 1:p
        row.names(x$pcor_mat) <- 1:p
        cat("Partial Correlations \n \n")
        print(x$pcor_mat * x$Adj_01, ...)
        cat("--- \n \n")
        cat("Adjancency (null) \n \n")
        colnames(x$Adj_01) <- 1:p
        row.names(x$Adj_01) <- 1:p
        print(x$Adj_01, ...)
        cat("--- \n")
      }
    }
  }

  if(x$alternative == "two.sided"){
    cat("Bayes Factor:", x$BF_cut, "\n")
    if(isFALSE(summarize)){
      if(hyp == "H1"){
        cat("Connectivity:", round(mean(x$Adj_10[upper.tri(x$Adj_10)]) * 100, 1), "% \n")
      } else{
        cat("Connectivity:", round(mean(x$Adj_01[upper.tri(x$Adj_01)]) * 100, 1), "% \n")
      }

    }
    cat("--- \n")
    cat("Call:\n")
    print(x$call)
    cat("--- \n")
    if (!isFALSE(summarize)) {
      p <- ncol(x$pcor_mat)
      mat_names <- matrix(0, p, p)
      mat_names[] <-  unlist(lapply(1:p, function(z) paste(1:p, z, sep = "--")))
      if(log == TRUE){
        summ <-  cbind.data.frame(edge = mat_names[upper.tri(mat_names)],
                                  post_mean = x$pcor_mat[upper.tri(x$pcor_mat)],
                                  post_sd = x$pcor_sd[upper.tri(x$pcor_sd)],
                                  "BF 10" = log(x$BF_10[upper.tri(x$BF_10)]))
      } else {
        summ <-  cbind.data.frame(edge = mat_names[upper.tri(mat_names)],
                                  post_mean = x$pcor_mat[upper.tri(x$pcor_mat)],
                                  post_sd = x$pcor_sd[upper.tri(x$pcor_sd)],
                                  "BF 10" = x$BF_10[upper.tri(x$BF_10)])
      }
      cat("Hypotheses: \n")
      cat("H0: rho = 0\nH1: rho != 0", "\n")
      cat("--- \n")
      cat("Estimates: \n \n ")
      print(summ, row.names = FALSE, ...)
      cat("--- \n")
      cat("note: BF_10 is evidence in favor of H1")
  } else {
    if(hyp == "H1"){
      cat("Hypothesis: \n")
      cat("H1: rho != 0 \n")
      cat("--- \n")
      p <- ncol(x$partials_non_zero)
      colnames(x$partials_non_zero) <- 1:p
      row.names(x$partials_non_zero) <- 1:p
      cat("Partial Correlations \n \n")
      print(x$partials_non_zero)
      cat("--- \n \n")
      cat("Adjancency (non-zero) \n \n")
      colnames(x$Adj_10) <- 1:p
      row.names(x$Adj_10) <- 1:p
      print(x$Adj_10)
      cat("--- \n")
    }
    if(hyp == "H0"){
      cat("Hypothesis: \n")
      cat("H0: rho = 0 \n")
      cat("--- \n")
      p <- ncol(x$pcor_mat)
      colnames(x$pcor_mat) <- 1:p
      row.names(x$pcor_mat) <- 1:p
      cat("Partial Correlations \n \n")
      print(x$pcor_mat * x$Adj_01, ...)
      cat("--- \n \n")
      cat("Adjancency (null) \n \n")
      colnames(x$Adj_01) <- 1:p
      row.names(x$Adj_01) <- 1:p
      print(x$Adj_01, ...)
      cat("--- \n")
      cat("note: connectivity reflects conditionally independent relations")
    }
  }
}
}

