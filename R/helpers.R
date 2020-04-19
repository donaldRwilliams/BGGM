#' @importFrom stats coef cov2cor var dnorm lm
#' na.omit pnorm quantile rWishart
#' sd qnorm residuals fitted density
#' @importFrom utils combn
#' @importFrom foreach %dopar% foreach
#' @importFrom graphics plot
#' @importFrom Rdpack reprompt
#' @import ggplot2



print_confirm <- function(x, ...){
  cat("BGGM: Bayesian Gaussian Graphical Models \n")
  cat("Type: Confirmatory Hypothesis Testing \n")
  cat("--- \n")
  cat("Call:\n")
  print(x$call)
  cat("--- \n")
  cat("Hypotheses: \n")

  if(length(x$hypotheses) == length(x$post_prob)){

    hyps <- data.frame( t(t(names(x$post_prob))), c(t(x$hypotheses)))
    colnames(hyps) <- NULL
    print(hyps, row.names = F)
  }



  if(length(x$hypotheses) != length(x$post_prob)){

    if(length(x$hypotheses) > 1){
      hyps <- data.frame( t(t(names(x$post_prob))), c(t(x$hypotheses),
                                                      paste("'not ", "H1-",
                                                            length(x$hypotheses), "'", sep = "")))
      colnames(hyps) <- NULL
      print(hyps, row.names = FALSE, right = F)

    } else{
      hyps <- data.frame( t(t(names(x$post_prob))), c(t(x$hypotheses),
                                                      paste("'not ", "H1", "'", sep = "")))
      colnames(hyps) <- NULL
      print(hyps, row.names = FALSE )

    }

  }

  cat("--- \n")
  cat("Posterior prob: \n")

  temp <- data.frame( (paste("p(", names(x$post_prob), "|Y) = ", round(x$post_prob, 4),  sep = "")))
  colnames(temp) <- ""

  print(temp, row.names = F, right = F)
  cat("--- \n")
  cat('Bayes factor matrix: \n')
  print(t(x$BF_matrix))
  cat("--- \n")
  cat("note: equal hypothesis prior probabilities")
}


print_ggm_compare_bf <- function(object, ...) {
  x <- object
  name_temp <- matrix(0, object$p, object$p)

  name_temp[] <-
    unlist(lapply(1:object$p , function(x)
      paste(1:object$p, x, sep = "--")))

  edge_name <- name_temp[upper.tri(name_temp)]

  groups <- length(object$info$dat)

  BF_10 <- 1 /  object$BF_01[upper.tri(object$BF_01)]
  prob_H1 <- BF_10 / (1 + BF_10)
  prob_H0 <- 1 - prob_H1

  if (groups == 2) {
    pcor_diff <- apply(object$post_samps[[1]], 2,  BGGM::fisher_z2r) -
      apply(object$post_samps[[2]], 2,  BGGM::fisher_z2r)

    sd_diff <- apply(pcor_diff, 2, sd)
    results <- data.frame(
      Edge = edge_name,
      Estimate = round(unlist(object$mu_diff), 3),
      Est.Error = round(sd_diff, 3),

      Pr.H0 = round(prob_H0, 3),
      Pr.H1 = round(prob_H1, 3)
    )


  } else {
    results <- data.frame(
      Edge = edge_name,

      Pr.H0 = round(prob_H0, 3),
      Pr.H1 = round(prob_H1, 3)
    )

  }
  x <- list(results = results,
                          object = object)
  # class(returned_object) <- "summary.ggm_compare_bf"
  # return(returned_object)

  cat("BGGM: Bayesian Gaussian Graphical Models \n")
  cat("--- \n")
  cat("Type: GGM Compare with Bayesian Hypothesis Testing \n")
  p <- x$object$info$dat_info$p[1]
  cat("Posterior Samples:", x$object$iter, "\n")
  cat("Observations: \n")
  groups <- length(x$object$info$dat)
  for (i in 1:groups) {
    cat("  Group",
        paste(i, ":", sep = "") ,
        x$object$info$dat_info$n[[i]],
        "\n")
  }
  cat("Variables (p):", p, "\n")
  cat("Edges:", .5 * (p * (p - 1)), "\n")
  cat("Groups:", nrow(x$object$info$dat_info), "\n")
  cat("Prior SD:", x$object$prior_sd, "\n")
  cat("--- \n")
  if (is.null(x$hypotheses)) {
    cat("Call: \n")
    print(x$object$call)
    cat("--- \n")
  }
  cat("Hypotheses:\n")
  cat("H0: rho = 0\n")
  cat("H1: rho != 0\n")
  cat("--- \n")
  cat("Estimates:\n")
  print(x$results, row.names = F, ...)
  cat("--- \n")

}



print_select_ggm_compare_bf <- function(x,...){
  object <- x
  cat("BGGM: Bayesian Gaussian Graphical Models \n")
  cat("--- \n")
  cat("Type: GGM Compare with Bayesian Hypothesis Testing \n")
  # number of iterations
  p <- object$object$info$dat_info$p[1]
  cat("Posterior Samples:", object$object$iter, "\n")
  cat("Observations: \n")
  groups <- length(object$object$info$dat)
  for (i in 1:groups) {
    cat("  Group",
        paste(i, ":", sep = "") ,
        object$object$info$dat_info$n[[i]],
        "\n")
  }
  # number of variables
  cat("Variables (p):", object$object$p, "\n")
  # number of edges
  cat("Edges:", .5 * (object$object$p * (object$object$p-1)), "\n")
  cat("--- \n")
  cat("Call: \n")
  print(object$call)
  cat("--- \n")
  cat("Selected:\n\n")
  cat("Adjacency non-zero \n \n")
  colnames(object$adj_10) <- 1:p
  row.names(object$adj_10) <- 1:p
  print(object$adj_10)
  cat("\n")
  cat("Adjacency zero \n \n")
  colnames(object$adj_01) <- 1:p
  row.names(object$adj_01) <- 1:p
  print(object$adj_01)
  cat("--- \n")
  cat("note: matrices (e.g., selected partial correlations) are in the select object")
}





print_select_explore <- function(x, hyp = "H1",
                                   log = TRUE, summarize = FALSE, ...){

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

      summ <-  cbind.data.frame(Edge = x$post_prob$edge,
                                Estimate = x$pcor_mat[upper.tri(x$pcor_mat)],
                                Est.Error = x$pcor_sd[upper.tri(x$pcor_sd)],
                                "Pr.H0" = x$post_prob[,2],
                                "Pr.H1" = x$post_prob[,3],
                                "Pr.H2" = x$post_prob[,4])

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
        summ <-  cbind.data.frame(Edge = mat_names[upper.tri(mat_names)],
                                  Estimate = x$pcor_mat[upper.tri(x$pcor_mat)],
                                  Est.Error = x$pcor_sd[upper.tri(x$pcor_sd)],
                                  "BF 20" = log(x$BF_20[upper.tri(x$BF_20)]),
                                  "BF 01" = log(x$BF_01[upper.tri(x$BF_01)]))
      } else{

        summ <-  cbind.data.frame(Edge = mat_names[upper.tri(mat_names)],
                                  Estimate = x$pcor_mat[upper.tri(x$pcor_mat)],
                                  Est.Error = x$pcor_sd[upper.tri(x$pcor_sd)],
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
        summ <-  cbind.data.frame(Edge = mat_names[upper.tri(mat_names)],
                                  Estimate = x$pcor_mat[upper.tri(x$pcor_mat)],
                                  Est.Error = x$pcor_sd[upper.tri(x$pcor_sd)],
                                  "BF 20" = log(x$BF_20[upper.tri(x$BF_20)]),
                                  "BF 01" = log(x$BF_01[upper.tri(x$BF_01)]))
      } else{

        summ <-  cbind.data.frame(Edge = mat_names[upper.tri(mat_names)],
                                  Estimate = x$pcor_mat[upper.tri(x$pcor_mat)],
                                  Est.Error = x$pcor_sd[upper.tri(x$pcor_sd)],
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
        summ <-  cbind.data.frame(Edge = mat_names[upper.tri(mat_names)],
                                  Estimate = x$pcor_mat[upper.tri(x$pcor_mat)],
                                  Est.Error = x$pcor_sd[upper.tri(x$pcor_sd)],
                                  "BF 10" = log(x$BF_10[upper.tri(x$BF_10)]))
      } else {
        summ <-  cbind.data.frame(Edge = mat_names[upper.tri(mat_names)],
                                  Estimate = x$pcor_mat[upper.tri(x$pcor_mat)],
                                  Est.Error = x$pcor_sd[upper.tri(x$pcor_sd)],
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


print_summary_explore <- function(x,...){
  summary(x$dat_results, summarize = TRUE)
}


print_explore <- function(x,...){
  cat("BGGM: Bayesian Gaussian Graphical Models \n")
  cat("--- \n")
  cat("Type: Hypothesis Testing (Exploratory) \n")
  cat("Posterior Samples:", x$iter, "\n")
  cat("Observations (n):", nrow(x$dat), "\n")
  cat("Variables (p):", x$p, "\n")
  cat("Edges:", .5 * (x$p * (x$p-1)), "\n")
  cat("Delta:", x$delta, "\n")
  cat("--- \n")
  cat("Call: \n")
  print(x$call)
  cat("--- \n")
  cat("Date:", date(), "\n")
}




print_select_ggm_compare_estimate <- function(x,...){
  object <- x
  cat("BGGM: Bayesian Gaussian Graphical Models \n")
  cat("--- \n")
  cat("Type: GGM Compare with the Posterior Distribution\n")
  # number of iterations
  cat("Posterior Samples:", object$object$iter, "\n")
  if (is.null(object$rope)) {

    cat("Credible Interval:",  gsub("*0.","", formatC( round(object$cred, 4), format='f', digits=2)), "% \n")

  } else {

    cat("Probability:", object$prob, "\n")
    cat("Region of Practical Equivalence:", "[", -1 * object$rope, ", ", object$rope, "]", "\n", sep = "")
  }
  # number of observations
  cat("Observations (n):\n")
  groups <- length(object$object$info$dat)
  for(i in 1:groups){
    cat("  Group", paste( i, ":", sep = "") , object$object$info$dat_info$n[[i]], "\n")
  }
  # number of variables
  cat("Variables (p):", object$object$p, "\n")
  # number of edges
  cat("Edges:", .5 * (object$object$p * (object$object$p-1)), "\n")
  cat("--- \n")
  cat("Call: \n")
  print(object$call)
  cat("--- \n")
  cat("Selected:\n\n")
  cat("Partial Correlations \n\n")
  if(is.null(object$rope)){
    for(i in 1:length(object$mat_adj)){
      cat(names(object$mat_adj)[[i]], "\n")

      print(object$mat_pcor[[i]])
      cat("\n")
    }
  } else {
    for(i in 1:length(object$mat_in_adj)){
      cat(names(object$mat_in_adj)[[i]], "\n")

      print(object$mat_out_pcor[[i]])
      cat("\n")
    }
    cat("--- \n")
    cat("note: null matrices in the select object")
  }
}


print_select_estimate <- function(x, summarize = FALSE, ...){
  # x <- object
  p <- ncol(x$partials_non_zero)
  cat("BGGM: Bayesian Gaussian Graphical Models \n")
  cat("--- \n")
  if(!is.null(x$analytic)){
    cat("Type: Selected Graph (Analytic Solution) \n")
  } else{
    cat("Type: ", x$type, "\n")

  }

  if(isFALSE(summarize)){
    if(is.null(x$rope)){
      cat("Credible Interval:",  gsub("*0.","", formatC( round(x$ci, 4), format='f', digits=2)), "% \n")
      cat("Connectivity:", round(mean(x$adjacency_non_zero[upper.tri(x$adjacency_non_zero)]) * 100, 1), "% \n")
      cat("--- \n")
      cat("Call:\n")
      print(x$call)
      cat("--- \n")
      cat("Selected:\n \n")
      colnames( x$partials_non_zero)  <- 1:p
      row.names( x$partials_non_zero) <- 1:p
      colnames( x$partials_non_zero)  <- 1:p
      row.names( x$partials_non_zero) <- 1:p
      cat("Partial correlations \n \n")
      print(x$partials_non_zero, digits = 2)
      cat("--- \n \n")
      cat("Adjacency \n \n")
      colnames(x$adjacency_non_zero) <- 1:p
      row.names(x$adjacency_non_zero) <- 1:p
      print(x$adjacency_non_zero)
      cat("--- \n")

    } else{
      cat("Probability:", x$prob, "\n")
      cat("Region of Practical Equivalence:", "[", -1 * x$rope, ", ", x$rope, "]", "\n", sep = "")
      cat("Connectivity:", round(mean(x$adjacency_non_zero[upper.tri(x$adjacency_non_zero)]) * 100, 1), "% \n")
      cat("--- \n")
      cat("Call:\n")
      print(x$call)
      cat("--- \n")
      cat("Selected:\n \n")
      colnames(x$partials_non_zero) <- 1:p
      row.names(x$partials_non_zero) <- 1:p
      cat("Partial correlations \n \n")
      print(x$partials_non_zero, digits = 2)
      cat("--- \n \n")
      cat("Adjacency non-zero \n \n")
      colnames(x$adjacency_non_zero) <- 1:p
      rownames(x$adjacency_non_zero) <- 1:p
      print(x$adjacency_non_zero)
      cat("--- \n \n")
      cat("Adjacency zero \n \n")
      colnames(x$adjacency_zero) <- 1:p
      rownames(x$adjacency_zero) <- 1:p
      print(x$adjacency_zero)
    }
  }
  if(isTRUE(summarize)){
    if(isTRUE(x$analytic)){
      stop("summary not available for the analytic solution")
    }
    if(is.null(x$rope)){
      p <- ncol(x$partials_non_zero)
      mat_names <- mu_mat <- ci_low <- ci_up <- mat_temp <- matrix(0, p, p)
      mat_names[] <-  unlist(lapply(1:p, function(z) paste(1:p, z, sep = "--")))

      low <- (1 - x$ci) / 2
      up  <-  1 - low

      mu_mat[] <-  colMeans(x$pcor_samples)
      sd <-  x$pcor_sd[upper.tri(x$pcor_sd)]
      cis <- apply(x$pcor_samples, 2, quantile, c(low, up))
      ci_low[] <- cis[1,]
      ci_up[] <- cis[2,]

      summ <- data.frame(edge = mat_names[upper.tri(mat_names)],
                         post_mean = mu_mat[upper.tri(mu_mat)],
                         post_sd = x$pcor_sd[upper.tri(x$pcor_sd)],
                         temp1 = ci_low[upper.tri(ci_low)],
                         temp2 = ci_up[upper.tri(ci_up)],
                         check.names = F)

      colnames(summ) <- c("Edge", "Estimate", "Est.Error",  paste(c("lb.", "ub."),
                                                                  gsub("*0.","", formatC( round(x$ci, 4), format='f', digits=2)), "%", sep = ""))
      cat("Credible Interval:", gsub("^.*\\.","", x$ci), "% \n")
      cat("Connectivity:", round(mean(x$adjacency_non_zero[upper.tri(x$adjacency_non_zero)]) * 100, 1), "% \n")
      cat("--- \n")
      cat("Call:\n")
      print(x$call)
      cat("--- \n")
      cat("Estimates: \n \n")
      print(summ, row.names = F,...)
      cat("--- \n")

    }else{
      cat("Probability:", x$prob, "\n")
      cat("Region of Practical Equivalence:", "[", -1 * x$rope, ", ", x$rope, "]", "\n", sep = "")
      cat("Connectivity:", round(mean(x$adjacency_non_zero[upper.tri(x$adjacency_non_zero)]) * 100, 1), "% \n")
      cat("--- \n")
      cat("Call:\n")
      print(x$call)
      cat("--- \n")
      cat("Pr.out: post prob outside of rope \n")
      cat("Pr.in: post prob inside of rope \n")
      cat("--- \n")

      p <- ncol(x$partials_non_zero)
      mat_names <- mu_mat <- rope_in  <- matrix(0, p, p)
      mat_names[] <-  unlist(lapply(1:p, function(z) paste(1:p, z, sep = "--")))

      low <- (1 - x$ci) / 2
      up  <-  1 - low

      mu_mat[] <-  colMeans(x$pcor_samples)
      sd <-  x$pcor_sd[upper.tri(x$pcor_sd)]

      rope_in[] <- x$in_rope

      cat("Estimates: \n \n")
      summ <- data.frame(edge = mat_names[upper.tri(mat_names)],
                         post_mean = mu_mat[upper.tri(mu_mat)],
                         post_sd = x$pcor_sd[upper.tri(x$pcor_sd)],
                         "pr_out" = 1 - rope_in[upper.tri(rope_in)],
                         "pr_in" = rope_in[upper.tri(rope_in)],
                         check.names = F)

      colnames(summ) <- c("Edge", "Estimate",
                          "Est.Error",  "Pr.out", "Pr.in")
      print(summ, row.names = F,...)
      cat("--- \n")
    }
  }
}
# print_select_estimate <- function(x, ...){
#   cat("BGGM: Bayesian Gaussian Graphical Models \n")
#   cat("--- \n")
#   if(is.numeric(x$rope)){
#     cat("Type: Selected Graph (Sampling) \n")
#   } else{
#     cat("Type: Selected Graph (Analytic Solution) \n")
#   }
#   if(is.null(x$rope)){
#     cat("Credible Interval:", gsub("^.*\\.","", x$ci), "% \n")
#     cat("Connectivity:", round(mean(x$adjacency[upper.tri(x$adjacency)]) * 100, 1), "% \n")
#     cat("--- \n")
#     cat("Call:\n")
#     print(x$call)
#     cat("--- \n")
#   } else{
#     cat("Probability:", x$prob, "\n")
#     cat("Region of Practical Equivalence:", "[", -1 * x$rope, ", ", x$rope, "]", "\n", sep = "")
#     cat("Connectivity:", round(mean(x$adjacency_non_zero[upper.tri(x$adjacency_non_zero)]) * 100, 1), "% \n")
#     cat("--- \n")
#     cat("Call:\n")
#     print(x$call)
#     cat("--- \n")
#   }
# }


print_summary_estimate <- function(x, ...) {
  # analytic == TRUE
  if (isTRUE(x$object$analytic)) {
    cat(print(x$object), "\n")
    cat("note: posterior summary not available for analytic solution")
  } else {
    cat("BGGM: Bayesian Gaussian Graphical Models \n")
    cat("--- \n")
    # number of iterations
    cat("Posterior Samples:", x$object$iter, "\n")
    # number of observations
    cat("Observations (n):", nrow(x$object$dat), "\n")
    # number of variables
    cat("Variables (p):", x$object$p, "\n")
    # number of edges
    cat("Edges:", .5 * (x$object$p * (x$object$p - 1)), "\n")
    cat("--- \n")
    cat("Call: \n")
    print(x$object$call)
    cat("--- \n")
    cat("Estimates:\n\n")
    print(x$dat_results, row.names = F)
    cat("--- \n")

  }
}



print_post_pred <- function(x,...){
  if(length(x) != 2){
    class(x) <- ""
    x <- round(x, 3)
    print(x)
  } else {
    cat("'summary = FALSE' not printed. See object contents")
  }
}




print_summary_metric <- function(x, digits = 2,...){
  cat("BGGM: Bayesian Gaussian Graphical Models \n")
  cat("--- \n")
  if(x$metric == "bayes_R2"){
    cat("Metric:", "Bayes R2\n")
  } else if(x$metric == "bayes_R2_diff"){
    cat("Metric:", "Bayes R2 Difference \n")
  } else {
    cat("Metric:", x$metric, "\n")

  }
  cat("Type:", x$type, "\n")
  cat("Credible Interval:", x$cred, "\n")
  cat("--- \n")
  cat("Estimates:\n\n")
  dat <- x$summary
  colnames(dat) <- c(colnames(dat)[1:3], "Cred.lb", "Cred.ub")
  print(as.data.frame( sapply(dat , round, digits)),
        row.names = FALSE)
}






print_summary_ggm_compare_ppc <- function(x, ...){
  cat("BGGM: Bayesian Gaussian Graphical Models \n")
  cat("--- \n")
  if(x$object$type == "nodewise"){
    cat("Type: GGM Comparison (Nodewise Predictive Check) \n")
  } else{
    cat("Type: GGM Comparison (Global Predictive Check) \n")
  }
  p <- x$object$info$dat_info$p[1]
  cat("Posterior Samples:", x$object$iter, "\n")

  groups <- length(x$object$info$dat)
  for (i in 1:groups) {
    cat("  Group",
        paste(i, ":", sep = "") ,
        x$object$info$dat_info$n[[i]],
        "\n")
  }
  cat("Variables (p):", p, "\n")
  cat("Edges:", .5 * (p * (p-1)), "\n")
  cat("--- \n")
  cat("Call: \n")
  print(x$object$call)
  cat("--- \n")
  cat("Estimates: \n \n")
  if(x$object$type == "global"){
    print(x$results, right = T, row.names = F,...)
    cat("--- \n")
    cat("note: \np_value = p(T(Y_rep) > T(y)|Y)\nKLD = (symmetric) Kullback-Leibler divergence")
  }
  if(x$object$type == "nodewise"){
    for(i in 1:length(x$object$obs_jsd)){
      cat(do.call(rbind, x$object$names)[[i]], "\n")
      print(   x$results[[i]],  row.names = F, ...)
      cat("\n")
    }

    cat("--- \n")
    cat("note: \np_value = p(T(Y_rep) > T(y)|Y)\nKLD = (symmetric) Kullback-Leibler divergence")
  }
}


print_summary_ggm_estimate_compare <- function(x,...){

  cat("BGGM: Bayesian Gaussian Graphical Models \n")
  cat("--- \n")
  cat("Type: GGM Compare with the Posterior Distribution\n")
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
  cat("Variables (p):", x$object$p, "\n")
  # number of edges
  cat("Edges:", .5 * (x$object$p * (x$object$p - 1)), "\n")
  cat("--- \n")
  cat("Call: \n")
  print(x$object$call)
  cat("--- \n")
  cat("Estimates:\n")
  for (i in 1:nrow(x$object$info$pairwise)) {
    cat("\n", names(x$object$pcors_diffs[[i]]), "\n")

    print(x$dat_results[[i]], row.names = FALSE,...)

  }
  cat("--- \n")
}







summary_ggm_compare_ppc <- function(object, ...){

  p <- object$info$dat_info$p[1]

  if (object$type == "global") {
    results <- data.frame(
      contrast = do.call(rbind, object$names),
      KLD =  do.call(rbind, object$obs_jsd),
      p_value = object$pvalue
    )
  }
  if (object$type == "nodewise") {
    results <- list()
    for (i in 1:length(object$obs_jsd)) {
      results[[i]] <-
        data.frame(
          node = 1:p ,
          KLD =  round(do.call(rbind, object$obs_jsd[[i]]), 3),
          p_value = unlist(object$pvalue[[i]])
        )
      names(results)[[i]] <- object$names[[i]]
    }
  }

  returned_object <- list(results = results,
                          object = object)
}


print_ggm_compare_ppc <- function(x,...){

  print_summary_ggm_compare_ppc((summary_ggm_compare_ppc(x)))
}




print_ggm_compare <- function(x, ...){
  cat("BGGM: Bayesian Gaussian Graphical Models \n")
  cat("--- \n")
  cat("Type: GGM Compare with the Posterior Distribution\n")
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
  cat("Edges:", .5 * (x$p * (x$p-1)), "\n")
  cat("--- \n")
  cat("Call: \n")
  print(x$call)
  cat("--- \n")
  cat("Date:", date(), "\n")
}



print_coef <- function(x,...){

  res_sigma <- x$inv_2_beta$sigma

  lb <- (1 - x$cred) / 2

  ub <- 1 - lb

  cred_int <- stats::quantile(res_sigma, prob = c(lb, ub))

  res_sigma_summ <- data.frame(Estimate = mean(res_sigma),
                               Est.Error = sd(res_sigma),
                               t(cred_int))

  # R2
  ypred <- t(apply(as.matrix(x$inv_2_beta$betas)[1:x$iter,], 1,
                   function(z)  z %*% t(as.matrix(x$data[,- x$node]))))

  r2 <- R2_helper(ypred, x$data[,x$node], ci_width = x$cred)
  cred_in <- stats::quantile(r2$R2, prob = c(lb, ub))

  res_r2_summ <- data.frame(Estimate = mean(r2$R2), Est.Error = sd(r2$R2), t(cred_in))

  colnames(res_sigma_summ) <- c("Estimate", "Est.Error", "CrI.lb", "CrI.ub")

  colnames(res_r2_summ) <- c("Estimate", "Est.Error", "CrI.lb", "CrI.ub")

  colnames(x$summary_inv_2_beta) <- c("Node", "Estimate",
                                      "Est.Error", "Cred.lb",
                                      "Cred.ub")
  cat("BGGM: Bayesian Gaussian Graphical Models \n")
  cat("--- \n")
  cat("Type: Inverse to Regression \n")
  cat("Credible Interval:",  gsub("*0.","", formatC( round(x$cred, 4), format='f', digits=2)), "% \n")
  cat("Node:", x$node, "\n")
  cat("--- \n")
  cat("Call: \n")
  print(x$call)
  cat("--- \n")
  cat("Coefficients: \n \n")
  summary_inv_2_beta <- data.frame(x$summary_inv_2_beta,
                                   check.names = F)
  print(summary_inv_2_beta, row.names = FALSE, ...)
  cat("--- \n")
  cat("Sigma:\n\n")
  print(round(res_sigma_summ, 3), row.names = FALSE, ...)
  cat("--- \n")
  cat("Bayesian R2:\n\n")
  print(round(res_r2_summ, 3), row.names = FALSE, ...)
  cat("--- \n")
}


print_map <- function(x,...){
  cat("BGGM: Bayesian Gaussian Graphical Models \n")
  cat("--- \n")
  cat("Method: Maximum A Posteriori \n")
  cat("--- \n")
  print(x$pcor)
}



print_fitted <- function(x,...){
  if(length(x) == 2){

    cat("'summary = FALSE' not printed. See object contents")

  } else {

    class(x) <- ""
    x <- round(x, 3)
    print(x)
  }
}

print_predict <- function(x,...){
  if(length(x) == 2){

    cat("'summary = FALSE' not printed. See object contents")

  } else {

    class(x) <- ""
    x <- round(x, 3)
    print(x)
  }
}


print_estimate <- function(x, ...){
  cat("BGGM: Bayesian Gaussian Graphical Models \n")
  cat("--- \n")
  # analytic == TRUE
  if(!isFALSE( x$analytic)){
    cat("Type: Estimation (Analytic Solution) \n")
  }
  # analytic  == FALSE
  if(isFALSE( x$analytic)){
    cat("Type:", x$type, "\n")
  }
  # number of iterations
  cat("Posterior Samples:", x$iter, "\n")
  # number of observations
  cat("Observations (n):", nrow(x$dat), "\n")
  # number of variables
  cat("Variables (p):", x$p, "\n")
  # number of edges
  cat("Edges:", .5 * (x$p * (x$p-1)), "\n")
  cat("--- \n")
  cat("Call: \n")
  print(x$call)
  cat("--- \n")
  cat("Date:", date(), "\n")
}



post_prob <- function(data){
  p1 <- sum(data>0)/length(data)
  p2 <- sum(data<0)/length(data)
  p<- 1 - min(p1,p2)
  return(p)
}


R2_ppc <- function(fit, betas, adj, which_one, sims){

  # data
  dat <- fit$dat

  # sample size
  n <- nrow(dat)

  # selected betas
  beta <- sweep(as.matrix(betas$betas[[which_one]]) , MARGIN=2,
                adj[which_one,-which_one], `*`)

  # sigmas
  sigma <- betas$sigmas[[which_one]]

  # posterior predictive
  ppc <- t(sapply(1:sims, function(s) rnorm(n = n,
                                            mean = dat[,-which_one] %*% beta[s,],
                                            sd = sigma[s])))
  # fitted values
  pred <- t(sapply(1:sims, function(s) dat[,-which_one] %*% beta[s,]))
  # r <- rowSums(pred^2)/ rowSums(ppc^2)

  # Bayes R2
  apply(pred, 1, var) / apply(ppc, 1, var)

}

convert_colnames <- function(hyp, Y){
  names_temp <- unlist(strsplit( strsplit(hyp, " ")[[1]], "--"))
  names_temp <- paste(names_temp, collapse = " ")
  names_temp <- unique(strsplit(gsub("[^[:alnum:] ]", "", names_temp), " +")[[1]])

  if(any(names_temp == "0")){
  names_temp <- names_temp[-which(names_temp == "0" )]
  }
  if(!all(names_temp %in% colnames(Y))){
    stop("node names not found in the data")
  }
  for(i in 1:length(names_temp)){
    id <- which(names_temp[[i]]  == colnames(Y))
    hyp <- gsub(x = hyp, pattern = names_temp[[i]],  replacement = id)
  }
  hyp
}

compare_predict_helper <- function(x, ci_width){
  post_mean <- mean(x)
  post_sd <- stats::sd(x)
  low <- (1 - ci_width) / 2
  up  <-  1 - low
  interval <-  t(stats::quantile(x, c(low, up)))
  summ <-  round(cbind.data.frame(post_mean = post_mean,
                                  post_sd = post_sd,
                                  interval), 3)
}

# delta give prior_sd
delta_solve = function(x){
  (x^2)^-1 - 1
}

# fisher z to r
z2r <- function (z) {
  (exp(2 * z) - 1)/(1 + exp(2 * z))
}

# lower triangle of matrix
get_lower_tri<-function(cormat){
  cormat[upper.tri(cormat)] <- NA
  return(cormat)
}

# analytic solution
analytic_solve <- function(X){
  # sample size
  n <- nrow(X)
  p <- ncol(X)
  # centererd mat
  X <- scale(X, scale = T)
  # scale matrix
  S <- t(X) %*% X
  inv_mu <-  solve(S + diag(10^-5,  p)) * (n)
  inv_var <-  (n + p + 1)*(solve(S + diag(10^-5, p) )^2 + tcrossprod(diag(solve(S + diag(10^-5, p)))))
  inv_cor <- diag( 1 / sqrt((diag(inv_mu)))) %*% inv_mu %*% diag( 1 / sqrt((diag(inv_mu))) )
  partials <- inv_cor * -1 + diag(2, p)
  list(inv_mu = inv_mu,
       inv_var = inv_var,
       partial = partials)
}

# summarize coefficients
beta_summary <- function(x, node, ci_width, samples){

  # convert inverse to beta
  x <- inverse_2_beta(x, samples = samples)

  # stop if not the correct class
  if(class(x) != "inverse_2_beta"){
    stop("class must be inverse_2_beta")
  }

  # check ci_width is allowed
  if(ci_width >= 1 | ci_width <= 0){
    stop("ci_width must be between 0 and 1")
  }
  returned_object <- lapply(node, function(y) summary_beta_helper(node =  y,
                                                                  x = x,
                                                                  ci_width))
  class(returned_object) <- "beta_summary"
  returned_object$betas <- x$betas[[node]]
  returned_object$sigma <- x$sigma[[node]]
  returned_object$call <- match.call()
  returned_object
}

rope_helper <- function(x, rope){
  mean(- rope < x & x < rope )

}

ci_helper <- function(x, ci_width){
  low <- (1 - ci_width) / 2
  up  <-  1 - low
  interval <-  stats::quantile(x, c(low, up))
  as.numeric(ifelse(interval[1] < 0 & interval[2] > 0, 0, 1))
}

Mo_risk_help_node <- function(x, post,  n1, n2, p){

  inv_mat <- post[,,x]

  Y_rep1 <- mvnfast::rmvn(n = n1,  mu = rep(0, p), sigma = stats::cov2cor(solve(inv_mat)))

  Y_rep2 <-  mvnfast::rmvn(n = n2, mu = rep(0, p), sigma = stats::cov2cor(solve(inv_mat)))

  jsd_node <- unlist(lapply(1:ncol(Y_rep1), function(z) node_jsd_help(z, Y_rep1, Y_rep2)))

  jsd_node



}


node_jsd_help <- function(x, Y_rep1, Y_rep2){

  Y_rep1 <- scale(Y_rep1)

  Y_rep2 <- scale(Y_rep2)

  pred1 <- Y_rep1[,-x]  %*% beta_helper(Y_rep1, x)

  pred2 <- Y_rep2[,-x]  %*%  beta_helper(Y_rep2, x)

  jsd_node <- (kl_func(stats::var(pred1), stats::var(pred2)) +
               kl_func(stats::var(pred2), stats::var(pred1))) * .5

  jsd_node

}


beta_helper <- function(x, which_one){
  y <- x[,which_one]
  X <- x[,-which_one]

  fit <- lm(y ~ 0 + X)
  coef(fit)

}

inverse_2_beta <- function(fit, samples = 500){
  fit <- fit[1:6]

  # check number of samples
  if(samples > fit$iter){
    stop("Samples used to compute R2 cannot be greater than the number used for fitting the model")

  }

  # get posterior estimate for precision matrix
  inv <-  fit$posterior_samples[,  grep("cov_inv", colnames(fit$posterior_samples))]

  # seperate estimates by row
  node_wise_elements <- lapply(split(colnames(inv), 1:fit$p), function(x) inv[,x])

  # convert off-diagonals to betas
  betas <- lapply(1:fit$p, function(x) node_wise_elements[[x]][-x] *  as.matrix((1/ node_wise_elements[[x]][x]) ) * -1)

  # convert diagonals to residual variance
  sigmas <- lapply(1:fit$p, function(x) as.matrix((1/ node_wise_elements[[x]][x])))


  # betas: select number of posterior samples
  betas <- lapply(betas, function(x)  x[1:samples,])

  # sigma: select number of posterior samples (sd scale)
  sigmas <- lapply(sigmas, function(x) sqrt(x[1:samples,]))


  # returned object
  returned_object <- list(betas = betas,
                          sigmas = sigmas,
                          p = fit$p,
                          data = fit$dat)

  # assign class
  class(returned_object) <- "inverse_2_beta"

  returned_object
}

summary_beta_helper <- function(x, node, ci_width){
  # index for row_names
  row_names <- 1:x$p

  # lower and upper ci
  low <- (1 - ci_width) / 2
  up <- 1 - low

  # beta posterior mean
  beta <- apply(x$betas[[node]], 2, mean)

  # beta poster sd
  post_sd <- apply(x$betas[[node]], 2, stats::sd)

  # lower and upper of posterior
  beta_ci <- t(apply(x$betas[[node]], 2, quantile, probs = c(low, up)))

  # sd of the outcome
  sd_y <- stats::sd(x$data[,node])

  # sd of the predictors
  sd_x <- apply(x$data[,-node], 2, stats::sd)

  # standardized (std) beta
  beta_std_temp <- x$betas[[node]] * (sd_x / sd_y)

  # beta std posterior mean
  beta_std <- apply(beta_std_temp, 2, mean)

  # beta std posterior sd
  beta_std_post_sd <- apply(beta_std_temp, 2, stats::sd)

  # lower and upper of posterior
  beta_std_ci <- t(apply(beta_std_temp, 2,  quantile, probs = c(low, up)))

  returned_object <- round(cbind(Node =  row_names[-node],
                                 beta,
                                 post_sd,
                                 beta_ci,
                                 beta_std,
                                 post_sd = beta_std_post_sd,
                                 beta_std_ci),3)

  # remove row names
  row.names(returned_object ) <- NULL

  # make list for naming purposes
  returned_object <- list(as.data.frame(returned_object))

  # list elemenet name as the response
  names(returned_object) <- paste("predicting node", node)
  returned_object$call <- match.call()
  returned_object
}


R2_helper <- function(ypred, y, ci_width) {
  low <- (1 - ci_width) / 2
  up  <-  1 - low
  e <- -1 * sweep(ypred, 2, y)
  var_ypred <- apply(ypred, 1, stats::var)
  var_e <- apply(e, 1, stats::var)
  r2 <- unlist(var_ypred / (var_ypred + var_e))
  ci <- quantile(r2, prob = c(low, up) )
  mu_r2 <- mean(r2)
  sd_r2 <- stats::sd(r2)
  summary_r2 <- c(post_mean = mu_r2, post_sd = sd_r2, ci)
  list(summary_r2 = summary_r2, R2 = r2)
}


MSE_helper <- function(ypred, y, ci_width){
  low <- (1 - ci_width) / 2
  up  <-  1 - low
  mse <- apply(ypred, MARGIN = 1, function(x){mean((x - y)^2)})
  ci <- quantile(mse, prob = c(low, up) )
  mu_mse <- mean(mse)
  sd_mse <- stats::sd(mse)
  summary_mse <- c(post_mean = mu_mse, post_sd = sd_mse, ci)
  list(summary_mse = summary_mse, MSE = mse)
}


name_helper <-  function(x){

  x <-  gsub("[A-z].*,", replacement = "", x)
  col_names <- gsub("[]]", "", x)
  col_names
}

error_helper <- function(ypred, y, ci_width, measure, sigmas =  NULL) {

  low <- (1 - ci_width) / 2

  up  <-  1 - low


  all_residual <- sweep(ypred, 2, y)

  if(measure == "mse"){
    out <- rowMeans(all_residual^2)
  }
  if(measure == "mae"){
    out <- rowMeans(abs(all_residual))
  }
  if(measure == "kl"){
    out <- kl_func(stats::var(y), sigmas^2)
    }
  ci <- quantile(out, prob = c(low, up) )
  mu_out <- mean(out)
  sd_out <- stats::sd(out)
  summary <- c(post_mean = mu_out, post_sd = sd_out, ci)
  list(summary = summary, error = out)
}

kl_func <- function(sigma_1, sigma_2){

  log(sqrt(sigma_2) / sqrt(sigma_1)) + (sigma_1 / (2 * sigma_2)) - .5

}

ppc_helper <- function(x, inv_g1, inv_cov, n, p){

  inv_mat <- matrix(0, p , p)
  inv_mat[,] <- as.numeric(inv_cov[x,])


  y_rep <- mvnfast::rmvn(n, mu = rep(0, p),sigma =  solve(inv_mat))

  S_rep <- t(y_rep) %*% y_rep

  theta_rep <- (n - 1) * solve(S_rep)

  KLD <- KL(Theta = inv_g1, hatTheta = theta_rep)

  JSD <- 0.5 * KL(Theta = inv_g1, hatTheta = theta_rep) + 0.5 * KL(hatTheta = theta_rep, Theta = inv_g1)

  QL <-  QL(Theta = inv_g1, hatTheta = theta_rep)

  FL <- sum((stats::cov2cor(inv_g1) *-1 - stats::cov2cor((theta_rep) * -1)^2))

  return <- list(KLD = KLD, JSD = JSD, QL = QL, FL = FL)
}

contrast_helper <- function(x){
  temp <- unlist(regmatches(x, gregexpr("[[:digit:]]+", x)))
  paste("Y", temp, sep = "_g", collapse = "_vs_")
}


axis_ticks_helper <- function(x){

  paste(stringr::str_sub(x, start = 1, end = 3),
        stringr::str_sub(x, start = 4, end = 5),
        stringr::str_sub(x, start = 6, end = 8))
}



KL = function(Theta,hatTheta){

  # Kuismin, M., & Sillanpaa, M. J. (2016). Use of Wishart prior and simple extensions for
  # sparse precision matrix estimation. PloS one, 11(2), e0148171.
  p = ncol(Theta)

  invTheta = solve(Theta,diag(1,p))

  kl  = 0.5 * (sum(diag(invTheta%*%hatTheta)) - log(det(invTheta%*%hatTheta)) - p)

  return(kl)

}

QL = function(Theta,hatTheta){
  # Kuismin, M., & Sillanpaa, M. J. (2016). Use of Wishart prior
  # and simple extensions for
  # sparse precision matrix estimation. PloS one, 11(2), e0148171.

  p = ncol(Theta)
  I = diag(1,p)

  invTheta = solve(Theta,I)

  osa = sum(diag(invTheta%*%hatTheta - I))
  tulos = osa^2

  return(tulos)

}

unbiased_cov <- function(x){
  x <- scale(x)
  n <- nrow(x) - 1
  mle_cov <- n^-1 * t(x) %*% x
  stats::cov2cor(solve(mle_cov))
}

Mo_risk_help <- function(x, post, n1, n2, p){
  inv_mat <- post[,,x]
  Y_rep1 <-  mvnfast::rmvn(n = n1,  mu = rep(0, p), sigma = stats::cov2cor(solve(inv_mat)))
  Y_rep2 <-  mvnfast::rmvn(n = n2, mu = rep(0, p), sigma =  stats::cov2cor(solve(inv_mat)))

  jsd <- 0.5 *  KL(unbiased_cov(Y_rep1), unbiased_cov(Y_rep2)) +
         0.5 *  KL(unbiased_cov(Y_rep2), unbiased_cov(Y_rep1))

  jsd
}



Y_combine <- function(...){

  dat <- list(...)

  dat <- lapply(1:length(dat), function(x) na.omit(dat[[x]]))

  dat_info <- lapply(1:length(dat), function(x) {
    p <- ncol(dat[[x]])

    n <- nrow(dat[[x]])

    data.frame(p = p, n = n)
  })

  list(dat = dat, dat_info =  do.call(rbind, dat_info),
       pairwise = t(combn(1:length(dat), 2)))
}



approx_sd <- function(r, n, k){
  sqrt((1-r^2)/(n  - k - 2))
}

positive_helper <- function(pcor, post_sd, BF_null){
  dens_greater  <- (1  - pnorm(0, pcor, post_sd)) * 2
  BF_null * dens_greater
}

negative_helper <- function(pcor, post_sd, BF_null){
  dens_less  <- pnorm(0, pcor, post_sd) * 2
  BF_null * dens_less
}


exhaustive_helper <- function(BF_null, BF_positive, BF_negative){
  c(BF_null, BF_positive, BF_negative) /  sum(BF_null, BF_positive, BF_negative)
}

symmteric_mat <- function(x){
  x[lower.tri(x)] <- t(x)[lower.tri(x)]
  x
}

colnames_helper <- function(x, col_names){
  colnames(x) <- col_names
}

sampling_helper = function(X,  nu, delta,  n_samples){
  X <- as.matrix(X)
  # number of variables
  p <- ncol(X)
  # number of observations
  n <- nrow(X)
  # number of partial correlations
  pcors <- (p * (p - 1)) / 2
  # names for the partial correlations
  col_names <- numbers2words(1:p)

  mat_name <- matrix(unlist(lapply(col_names, function(x) paste(col_names,x, sep = ""))), p , p)
  mat_name_up <- mat_name[upper.tri(mat_name)]
  mat_name_low <- mat_name[lower.tri(mat_name)]

  # center the data
  Xhat <- X - rep(1,n)%*%t(apply(X,2,mean))

  # scatter matrix
  S <- t(Xhat)%*%Xhat

  # storage
  pcor_store_up <- pcor_store_low <- prior_store_up <- prior_store_low <- matrix(NA, nrow = n_samples, ncol = pcors)
  inv_cov_store <-  array(NA, c(p, p, n_samples))

  # initial values
  Psi <- b_inv <- diag(p)

  for(i in 1:n_samples){
    # draw from posterior
    post <- post_helper(S = S, n = n, nu = nu, p = p, delta = delta, Psi = Psi, b_inv = b_inv * 10000)
    # store partials
    pcor_store_up[i,] <- post$pcors_post_up
    pcor_store_low[i,] <- post$pcors_post_low
    # store the inverse
    inv_cov_store[,,i] <- post$sigma_inv
    # draw from prior and store
    prior_samps <- prior_helper(nu = nu, delta = delta, p = p)
    prior_store_up[i,] <- prior_samps$pcors_prior_up
    prior_store_low[i, ] <- prior_samps$pcors_prior_up
    # Psi
    Psi <- post$Psi
  }

  # transform posterior samples
  fisher_z_post_up <- apply(pcor_store_up, 2, fisher_z)
  fisher_z_post_low <- apply(pcor_store_low, 2, fisher_z)

  fisher_z_prior_up <- apply(prior_store_up, 2, fisher_z)
  fisher_z_prior_low <- apply(prior_store_low, 2, fisher_z)


  colnames(fisher_z_prior_up) <- mat_name_up
  colnames(fisher_z_post_up) <- mat_name_up
  colnames(pcor_store_up) <- mat_name_up


  colnames(fisher_z_prior_low) <- mat_name_low
  colnames(fisher_z_post_low) <- mat_name_low
  colnames(pcor_store_low) <- mat_name_low



  # returned list
  list(fisher_z_post = cbind(fisher_z_post_up, fisher_z_post_low),
       pcor_post = cbind(pcor_store_up, pcor_store_low),
       inv_cov_post = inv_cov_store,
       pcor_prior = cbind(prior_store_up, prior_store_low),
       fisher_z_prior =cbind(fisher_z_prior_up, fisher_z_prior_low))
}

prior_helper <- function(nu, p, delta){
  # sample from inverse Wishart
  inv_wish_prior <- solve(rWishart(1, df =  delta + p - 1, diag(p) * 1000)[,,1], tol = 1e-20)

  # sample from Wishart
  sigma_inv_prior <- rWishart(1, df = nu, inv_wish_prior)[,,1]

  # partical correlation matrix
  pcor_mat_prior <- - diag(1/sqrt(diag(sigma_inv_prior)))%*%sigma_inv_prior%*%diag(1/sqrt(diag(sigma_inv_prior)))
  pcors_prior_up <- pcor_mat_prior[upper.tri(pcor_mat_prior)]
  pcors_prior_low <- pcor_mat_prior[lower.tri(pcor_mat_prior)]

  list(pcors_prior_up = pcors_prior_up, pcors_prior_low = pcors_prior_low)
}



post_helper <- function(S, n, nu, p, delta, Psi, b_inv){
  # precision matrix
  sigma_inv <- rWishart(1, delta + n - 1, solve(Psi+S,  tol =  1e-20))[,,1]

  # Psi
  Psi <- rWishart(1, nu + delta + p - 2, solve(sigma_inv + b_inv,  tol  = 1e-20))[,,1]

  # partial correlation matrix
  pcor_mat <- - diag(1/sqrt(diag(sigma_inv)))%*%sigma_inv%*%diag(1/sqrt(diag(sigma_inv)))
  pcors_post_up = pcor_mat[upper.tri(pcor_mat)]
  pcors_post_low = pcor_mat[lower.tri(pcor_mat)]

  # returned list
  list(pcors_post_up = pcors_post_up, pcors_post_low = pcors_post_low, sigma_inv = sigma_inv, Psi = Psi)
}

sampling <- function(X, nu, delta, n_samples = 20000, cores = 4){
  # register parallel
  cl <- parallel::makeCluster(cores)
  doParallel::registerDoParallel(cl)

  # samples for each "chain"
  samps <- rep(round(n_samples / cores), cores)
  chains <- cores

  # global variable
  i <- 1

  # sample from priors and posteriors
  samples <- foreach::foreach(i = 1:chains,
                              .export = c("fisher_z", "sampling_helper",
                                          "numbers2words", "prior_helper", "post_helper")) %dopar% {
                                                sampling_helper(X = X, nu = nu,
                                                                delta = delta,
                                                                n_samples = samps[i])
                                            }
  # stop cluster
  parallel::stopCluster(cl)

  return(samples)
}

fisher_z <- function(rho){
  .5 * log(( 1 + rho )/ ( 1 - rho ))
}

sd_helper <- function(post_samples, prior_at_zero){
  prior_at_zero /  dnorm(0, mean(post_samples), stats::sd(post_samples))
}

pcor_name_helper <- function(x){
  keep_vars <-  unlist(strsplit(gsub("[^[:alnum:] ]", "", x), " +"))
  keep_vars
}

framer <- function(x){
  pos_comparisons <- unlist(gregexpr("[<>=]", x))
  leftside <- rep(NA, length(pos_comparisons) + 1)
  rightside <- rep(NA, length(pos_comparisons) + 1)
  pos1 <- c(-1, pos_comparisons)
  pos2 <- c(pos_comparisons, nchar(x) + 1)
  for(i in seq_along(pos1)){
    leftside[i] <- substring(x, pos1[i] + 1, pos1[i+1] - 1)
    rightside[i] <- substring(x, pos2[i] + 1, pos2[i+1] - 1)
  }
  leftside <- leftside[-length(leftside)]
  rightside <- rightside[-length(rightside)]
  comparisons <- substring(x, pos_comparisons, pos_comparisons)
  data.frame(left = leftside,
             comp = comparisons,
             right = rightside,
             stringsAsFactors = FALSE)

}



create_matrices <- function(framed, varnames){
  k <- length(varnames)
  if(any(grepl(",", framed$left)) || any(grepl(",", framed$right))){
    if(nrow(framed) > 1){
      for(r in 1:(nrow(framed)-1)){
        if(all.equal(framed$right[r], framed$left[r+1])){
          if(substring(framed$right[r], 1, 1) == "(") {
            framed$right[r] <- sub("),.+", ")", framed$right[r])
            framed$left[r+1] <- sub(".+),", "", framed$left[r +1])
          } else{
            framed$right[r] <- sub(",.+", "", framed$right[r])
            framed$left[r+1] <- sub("[^,]+,", "", framed$left[r+1])
          }
        }
      }
    }

    commas_left <- framed$left[grep(",", framed$left)]
    commas_right <- framed$right[grep(",", framed$right)]
    if(isTRUE(any(!grepl("\\(.+)", commas_left))) || isTRUE(any(!grepl("\\(.+)", commas_right))) ||
       isTRUE(any(grepl(").+", commas_left))) || isTRUE(any(grepl(").+", commas_right))) ||
       isTRUE(any(grepl(".+\\(", commas_left))) || isTRUE(any(grepl(".+\\(", commas_right)))) {
      stop("Incorrect hypothesis syntax or extra character, check specification")
    }

    framed$left <- gsub("[()]", "", framed$left)
    framed$right <- gsub("[()]", "", framed$right)
    commas <- unique(c(grep(",", framed$left), grep(",", framed$right)))

    if(length(commas) > 0){
      multiples <- vector("list", length = length(commas))

      for(r in seq_along(commas)){
        several <- framed[commas,][r, ]

        if(several$comp == "="){

          several <- c(several$left, several$right)
          separate <- unlist(strsplit(several, split = ","))
          if(any(grepl("^$", several))) stop("Misplaced comma in hypothesis")
          converted_equality <- paste(separate, collapse = "=")
          multiples[[r]] <- framer(converted_equality)

        } else{
          leftvars <- unlist(strsplit(several$left, split = ","))
          rightvars <- unlist(strsplit(several$right, split = ","))
          if(any(grepl("^$", leftvars)) || any(grepl("^$", rightvars))) stop("Misplaced comma in hypothesis")

          left <- rep(leftvars, length.out = length(rightvars)*length(leftvars))
          right <- rep(rightvars, each = length(leftvars))
          comp <- rep(several$comp, length(left))

          multiples[[r]] <- data.frame(left = left, comp = comp, right = right, stringsAsFactors = FALSE)
        }
      }

      framed <- framed[-commas,]
      multiples <- do.call(rbind, multiples)
      framed <- rbind(multiples, framed)
    }
  }

  equality <- framed[framed$comp == "=",]
  inequality <- framed[!framed$comp == "=",]

  #****Equality part string-to-matrix
  if(nrow(equality) == 0) {
    R_e <- r_e <- NULL
  } else{
    outcomes <- suppressWarnings(apply(equality[, -2], 2, as.numeric))
    outcomes <- matrix(outcomes, ncol = 2, byrow = FALSE)
    if(any(rowSums(is.na(outcomes)) == 0)) stop("Value compared with value rather than variable, e.g., '2 = 2', check hypotheses")
    rows <- which(rowSums(is.na(outcomes)) < 2)
    specified <- t(outcomes[rows,])
    specified <- specified[!is.na(specified)]
    r_e <- ifelse(rowSums(is.na(outcomes)) == 2, 0, specified)
    r_e <- matrix(r_e)

    var_locations <- apply(equality[, -2], 2, function(x) ifelse(x %in% varnames, match(x, varnames), 0))
    var_locations <- matrix(var_locations, ncol = 2)

    R_e <- matrix(rep(0, nrow(equality)*length(varnames)), ncol = length(varnames))

    for(i in seq_along(r_e)){
      if(!all(var_locations[i, ] > 0)){
        R_e[i, var_locations[i,]] <- 1
      } else{
        R_e[i, var_locations[i,]] <- c(1, -1)
      }
    }
  }


  #****Inequality part string-to-matrix
  if(nrow(inequality) == 0) {
    R_i <- r_i <- NULL
  } else{
    outcomes <- suppressWarnings(apply(inequality[, -2], 2, as.numeric))
    outcomes <- matrix(outcomes, ncol = 2, byrow = FALSE)
    if(any(rowSums(is.na(outcomes)) == 0)) stop("Value compared with value rather than variable, e.g., '2 > 2', check hypotheses")
    cols <- which(rowSums(is.na(outcomes)) < 2)
    specified <- t(outcomes[cols,])
    specified <- specified[!is.na(specified)]
    r_i <- ifelse(rowSums(is.na(outcomes)) == 2, 0, specified)
    r_i <- matrix(r_i)

    leq <- which(inequality$comp == "<")
    var_locations <- apply(inequality[, -2], 2, function(x) ifelse(x %in% varnames, match(x, varnames), 0))
    var_locations <- matrix(var_locations, ncol = 2)

    R_i <- matrix(rep(0, nrow(inequality)*length(varnames)), ncol = length(varnames))

    for(i in seq_along(r_i)){
      if(!all(var_locations[i, ] > 0)){

        if(var_locations[i, 1] == 0){
          if(i %in% leq){
            value <-  1
          } else{
            r_i[i] <- r_i[i]*-1
            value <- -1
          }
        } else{
          if(i %in% leq){
            r_i[i] <- r_i[i]*-1
            value <-  -1
          } else{
            value <- 1
          }
        }

        R_i[i, var_locations[i,]] <- value

      } else{
        value <- if(i %in% leq) c(-1, 1) else c(1, -1)
        R_i[i, var_locations[i,]] <- value
      }
    }
  }

  #3)check comparisons----------------
  if(is.null(R_i)){
    comparisons <- "only equality"
  } else if(is.null(R_e)){
    comparisons <- "only inequality"
  } else{
    comparisons <- "both comparisons"
  }

  #set prior mean
  R_ei <- rbind(R_e,R_i)
  r_ei <- rbind(r_e,r_i)
  Rr_ei <- cbind(R_ei,r_ei)
  beta_zero <- MASS::ginv(R_ei)%*%r_ei

  if(nrow(Rr_ei) > 1){
    rref_ei <- pracma::rref(Rr_ei)
    nonzero <- rref_ei[,k+1]!=0
    if(max(nonzero)>0){
      row1 <- max(which(nonzero==T))
      if(sum(abs(rref_ei[row1,1:k]))==0){
        stop("Default prior mean cannot be constructed from constraints.")
      }
    }
  }


  list(R_i = R_i,
       r_i = r_i,
       R_e = R_e,
       r_e = r_e,
       R_ei = R_ei,
       Rr_ei = Rr_ei,
       r_ei = r_ei,
       beta_zero = beta_zero,
       comparisons = comparisons)

}

word2num <- function(word){
  wsplit <- strsplit(tolower(word)," ")[[1]]
  one_digits <- list(zero=0, one=1, two=2, three=3, four=4, five=5,
                     six=6, seven=7, eight=8, nine=9)
  teens <- list(eleven=11, twelve=12, thirteen=13, fourteen=14, fifteen=15,
                sixteen=16, seventeen=17, eighteen=18, nineteen=19)
  ten_digits <- list(ten=10, twenty=20, thirty=30, forty=40, fifty=50,
                     sixty=60, seventy=70, eighty=80, ninety=90)
  doubles <- c(teens,ten_digits)
  out <- 0
  i <- 1
  while(i <= length(wsplit)){
    j <- 1
    if(i==1 && wsplit[i]=="hundred")
      temp <- 100
    else if(i==1 && wsplit[i]=="thousand")
      temp <- 1000
    else if(wsplit[i] %in% names(one_digits))
      temp <- as.numeric(one_digits[wsplit[i]])
    else if(wsplit[i] %in% names(teens))
      temp <- as.numeric(teens[wsplit[i]])
    else if(wsplit[i] %in% names(ten_digits))
      temp <- (as.numeric(ten_digits[wsplit[i]]))
    if(i < length(wsplit) && wsplit[i+1]=="hundred"){
      if(i>1 && wsplit[i-1] %in% c("hundred","thousand"))
        out <- out + 100*temp
      else
        out <- 100*(out + temp)
      j <- 2
    }
    else if(i < length(wsplit) && wsplit[i+1]=="thousand"){
      if(i>1 && wsplit[i-1] %in% c("hundred","thousand"))
        out <- out + 1000*temp
      else
        out <- 1000*(out + temp)
      j <- 2
    }
    else if(i < length(wsplit) && wsplit[i+1] %in% names(doubles)){
      temp <- temp*100
      out <- out + temp
    }
    else{
      out <- out + temp
    }
    i <- i + j
  }
  return(list(word,out))
}

numbers2words <- function(x){
  ## Function by John Fox found here:
  ## http://tolstoy.newcastle.edu.au/R/help/05/04/2715.html
  ## Tweaks by AJH to add commas and "and"
  helper <- function(x){

    digits <- rev(strsplit(as.character(x), "")[[1]])
    nDigits <- length(digits)
    if (nDigits == 1) as.vector(ones[digits])
    else if (nDigits == 2)
      if (x <= 19) as.vector(teens[digits[1]])
    else trim(paste(tens[digits[2]],
                    Recall(as.numeric(digits[1]))))
    else if (nDigits == 3) trim(paste(ones[digits[3]], "hundred and",
                                      Recall(makeNumber(digits[2:1]))))
    else {
      nSuffix <- ((nDigits + 2) %/% 3) - 1
      if (nSuffix > length(suffixes)) stop(paste(x, "is too large!"))
      trim(paste(Recall(makeNumber(digits[
        nDigits:(3*nSuffix + 1)])),
        suffixes[nSuffix],"," ,
        Recall(makeNumber(digits[(3*nSuffix):1]))))
    }
  }
  trim <- function(text){
    #Tidy leading/trailing whitespace, space before comma
    text=gsub("^\ ", "", gsub("\ *$", "", gsub("\ ,",",",text)))
    #Clear any trailing " and"
    text=gsub(" and$","",text)
    #Clear any trailing comma
    gsub("\ *,$","",text)
  }
  makeNumber <- function(...) as.numeric(paste(..., collapse=""))
  #Disable scientific notation
  opts <- options(scipen=100)
  on.exit(options(opts))
  ones <- c("", "one", "two", "three", "four", "five", "six", "seven",
            "eight", "nine")
  names(ones) <- 0:9
  teens <- c("ten", "eleven", "twelve", "thirteen", "fourteen", "fifteen",
             "sixteen", " seventeen", "eighteen", "nineteen")
  names(teens) <- 0:9
  tens <- c("twenty", "thirty", "forty", "fifty", "sixty", "seventy", "eighty",
            "ninety")
  names(tens) <- 2:9
  x <- round(x)
  suffixes <- c("thousand", "million", "billion", "trillion")
  if (length(x) > 1) return(trim(sapply(x, helper)))
  helper(x)
}

samps_inv_helper <- function(x, p){
  inv <- paste("cov_inv", paste(paste("[", paste( 1:p, x, sep = ","), sep = ""), "]", sep = ""), sep = "")
  inv
}

samps_pcor_helper <- function(x, p){

  pcors <- paste("pcors", paste(paste("[", paste( 1:p, x, sep = ","), sep = ""), "]", sep = ""), sep = "")
  pcors
}

hyp_converter <- function(x){

  hyp_converted <- x

  extract_numbers <- unlist(stringr::str_extract_all(hyp_converted, "\\d+"))

  extract_numbers <- extract_numbers[unlist(extract_numbers) != 0 ]
  words <- NA
  for(i in 1:length(extract_numbers)){

    temp <- noquote(extract_numbers[i])
    words[i] <- numbers2words(as.numeric(temp))
    hyp_converted <- sub(temp, numbers2words(as.numeric(temp)), hyp_converted)


  }

  hyp_converted <- stringr::str_remove_all(hyp_converted, "--")

  list(hyp_converted = hyp_converted, words = words)
}

performance <- function(Estimate, True){

  True <- as.matrix(True)
  Estimate <- as.matrix(Estimate)

  # True Negative
  TN <- ifelse(True[upper.tri(True)] == 0 & Estimate[upper.tri(Estimate)] == 0, 1, 0); TN <- sum(TN)
  # False Positive
  FP <- ifelse(True[upper.tri(True)] == 0 & Estimate[upper.tri(Estimate)] != 0, 1, 0); FP <- sum(FP)
  # True Positive
  TP <- ifelse(True[upper.tri(True)] != 0 & Estimate[upper.tri(Estimate)] != 0, 1, 0); TP <- sum(TP)
  # False Negatives
  FN <- ifelse(True[upper.tri(True)] != 0 & Estimate[upper.tri(Estimate)] == 0, 1, 0); FN <- sum(FN)

  Specificity <- TN/(TN + FP)
  Sensitivity <- TP/(TP + FN)
  Precision <- TP/(TP + FP)

  Recall <- TP / (TP + FN)

  F1_score <- 2 * ((Precision * Recall) / (Precision + Recall))

  MCC <- (TP*TN - FP*FN)/sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN))

  results <- c(Specificity, Sensitivity, Precision, Recall,  F1_score, MCC)
  results_name <- c("Specificity", "Sensitivity", "Precision", "Recall",  "F1_score", "MCC")
  results <- cbind.data.frame(measure = results_name, score = results)
  list(results = results)


}

csws_labels <- ifelse(1:35 %in% c(7,10,16,24,29),
                      "Family Support",
                      ifelse(1:35 %in% c(3,12,20,25,35),
                             "Competition",
                             ifelse(1:35 %in% c(1,4,17,21,30),
                                    "Appearence",
                                    ifelse(1:35%in%c(2,8,18,26,31),
                                           "God's Love",
                                           ifelse(1:35 %in% c(13, 19, 22, 27,  33),
                                                  "Academic Competence",
                                                  ifelse(1:35 %in% c(5, 11, 14, 28, 34),
                                                         "Virtue", "Approval From Others"))))))

tas_labels <- ifelse(1:20 %in% c(1,3,6,7,9,13,14),
                     "Difficulty\nIdentifying Feelings",
                     ifelse(1:20 %in% c(2,4,11,12,17),
                            "Difficulty\nDescribing Feelings",
                            "Externally\nOriented Feelings"))

iri_labels <- ifelse(1:28 %in% c(3, 8, 11, 15, 21, 25, 28),
                     "Perspective Taking",
                     ifelse(1:28 %in% c(2, 4, 9, 14, 18, 20, 22),
                            "Empathic Concern",
                            ifelse(1:28 %in% c(1, 5, 7, 12, 16, 23, 26), "Fantasy",
                                   "Personal Distress")))

rsa_labels <- ifelse(1:33 %in% c(1, 4, 5, 32),
                     "Planned Future",
                     ifelse(1:33 %in% c(2, 11, 17, 25, 31, 33),
                            "Perception of Self",
                            ifelse(1:33 %in% c(3, 7, 13, 16, 24, 29),
                                   "Family Cohesion",
                                   ifelse(1:33  %in% c(6, 9, 10, 12, 15, 19, 27),
                                          "Social Resources",
                                          ifelse(1:33 %in% c(8, 14, 18, 21, 22, 26),
                                                 "Social Competence", "Structured Style")))))

globalVariables(c('Y1','Y2',
                  'X1', 'X2',
                  'contrast',
                  '..quantile..',
                  'value',
                  'node',
                  'BF',
                  'Edge',
                  'Estimate',
                  'selected',
                  'probability',
                  'cred.lb',
                  'sig',
                  'hyp',
                  'label',
                  'color',
                  'fit',
                  'post_mean',
                  'Error',
                  'density', 'Node',
                  'Post.mean',
                  'L1', 'lag', 'acf',
                  'iteration'))
