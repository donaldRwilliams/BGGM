#' Summarize \code{select.explore}
#'
#' @param x object of class \code{select.explore}
#' @param hyp hypothesis to summarize-e.g., \code{"H0"}, \code{"H1"}, etc.
#' @param log log tranform the Bayes factors
#' @param summarize \code{summarize = TRUE} provides a summary output of the edges, Bayes factor,
#'  and posterior probabilities. default is \code{FALSE}, which outputs the matrix for the respective hypothesis (i.e., \code{hyp})
#' @param ... \code{digits}
#'
#' @return
#' @export
#'
#' @examples
#'
#'Y <- BGGM::bfi[,1:10]
#'
#'# sample from posterior
#' fit_bf <- explore(Y, prior_sd = 0.5,
#'                   iter = 5000,
#'                   cores = 2)
#'
#'# select the graph
#' E <- select(fit_bf,
#'            hyp_prob = 0.90,
#'            alternative = "exhaustive")
#'
#'# hyp = "H0"
#'summary(E, hyp = "H0", summarize = F)
#'
summary.select.explore <- function(x, hyp = "H1",  log = TRUE, summarize = FALSE, ...){
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


    if(!isFALSE(summarize)){
      p <- ncol(x$pcor_mat)
      mat_names <- matrix(0, p, p)



      mat_names[] <-  unlist(lapply(1:p, function(z) paste(1:p, z, sep = "--")))

      if(log == TRUE){
        summ <-  cbind.data.frame(edge = mat_names[upper.tri(mat_names)],
                                  post_mean = x$pcor_mat[upper.tri(x$pcor_mat)],
                                  post_sd = x$pcor_sd[upper.tri(x$pcor_sd)],
                                  "BF 10" = log(x$BF_10[upper.tri(x$BF_10)]))
      } else{

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
    }
    else{

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
