#' Summarize \code{select.ggm_compare_bf}
#'
#' @param x
#' @param type "adj" for adjacency matrix (0's and 1's). "BF" for matrix of selected Bayes factors
#'
#' @return
#' @export
#'
#' @examples
summary.select.ggm_compare_bf <- function(x, type = "adj"){

  edges <- .5 * (x$p * (x$p-1))

  prop_alt <- sum( x$adj_10[upper.tri( x$adj_10 )]) / edges

  prop_null <- sum( x$adj_01[upper.tri( x$adj_01 )]) / edges

  incon <- 1 - (prop_alt + prop_null)

  group <- nrow( x$info$dat_info )

  cat("BGGM: Bayesian Gaussian Graphical Models \n")
  cat("--- \n")
  # hypothesis testing
  cat("Type: GGM Comparison (Bayesian Hypothesis Testing) \n")
  cat("Bayes Factor: 3 \n")
  cat("--- \n")
  cat("Call: \n")
  print(x$call)
  cat("--- \n")
  cat("Hypotheses: \n")
  cat("H0",paste("rho_g", 1:group, "_ij", sep = "", collapse = " = "), "\n")
  cat("H1", "'not H0'", "\n")
  cat("--- \n")
  cat("Evidence for H0:", round(prop_null * 100, 2), "%\n")
  cat("Evidence for H1:", round(prop_alt * 100, 2), "%\n")
  cat("Inconclusive:", round(incon * 100, 2), "%\n")

  if(type == "adj"){
    cat("--- \n\n")
    cat("Adjacency Matrix (H0) \n \n")
    temp_01 <- as.data.frame( x$adj_01)
    colnames(temp_01) <- 1:ncol(temp_01)
    print(temp_01)
    cat("--- \n \n")
    cat("Adjacency Matrix (H1) \n \n")
    temp_10 <- as.data.frame( x$adj_10)
    colnames(temp_10) <- 1:ncol(temp_10)
    print(temp_10)
    cat("--- \n")
  }
  if(type == "BF"){
    cat("--- \n\n")
    cat("Bayes Factor Matrix (H0) \n \n")
    temp_01 <- as.data.frame( round(x$BF_01_adj,2))
    colnames(temp_01) <- 1:ncol(temp_01)
    print(temp_01)
    cat("--- \n \n")
    cat("Bayes Factor Matrix (H1) \n \n")
    temp_10 <- as.data.frame(round( x$BF_10_adj,2))
    colnames(temp_10) <- 1:ncol(temp_10)
    print(temp_10)
    cat("--- \n")
  }
  if(type == "none"){
    cat("--- \n")
  }


}
