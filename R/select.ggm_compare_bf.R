#' Select Graphical Structure with the Bayes Factor
#'
#' @description Compare GGMs with the Bayes factor. This method allows for
#' assessing (relative) evidence for edge equality or edges differences across any number of groups. Further, confirmatory hypothesis testing
#' can be used to test predictions or expectations regarding difference or similarities in different groups (e.g., male vs. female).
#'
#'
#' @param object object of class \code{ggm_compare_bf}
#' @param BF_cut evidentiary threshold
#' @param ... currently ignored
#' @note The test provides relative evidence for whether all groups have the same edge strength (for each edge in the model). It is possible to test whether
#' many groups are all the same. In this case, if there is evidence for the alternative (not equal), information is not immediately available for which group differs
#' from the others. Thus pairwise group contrasts could be used after testing the equality of more than two groups.
#'
#' @return list of class \code{select.ggm_compare_bf}
#'
#' \itemize{
#' \item \code{BF_10} Bayes factors for the alternative ("not equal")
#' \item \code{BF_01} Bayes factors for the null hypothesis
#' \item \code{BF_10_adj} Bayes factor adjacency matrix for the alternative ("not equal")
#' \item \code{BF_01_adj} Bayes factor adjacency matrix for the null hypothesis
#' \item \code{adj_10} adjacency matrix for the selected edges (in favor of the "not equal")
#' \item \code{adj_01} adjacency matrix for the selected edges (in favor of the null hypothesis)
#' \item \code{pcor_mat_10} partial correlation matrix for the alternative ("not equal")
#' }
#'
#' @export
#'
#' @examples
#'
#' Y1 <- MASS::mvrnorm(500, rep(0, 16), Sigma = diag(16))
#' Y2 <- MASS::mvrnorm(500, rep(0, 16), Sigma = diag(16))
#' Y3 <- MASS::mvrnorm(500, rep(0, 16), Sigma = diag(16))
#'
#' ggm_bf <- ggm_compare_bf(Y1, Y2, Y3,
#'                         prior_sd = .20,
#'                         iter = 500,
#'                         cores = 2)
#'# select with BF_cut = 3
#' ggm_bf_sel <- select(ggm_bf, BF_cut = 3)
#'
#' # summary
#' summary(ggm_bf_sel)

select.ggm_compare_bf <- function(object, BF_cut = 3,...){

  x <- object
  BF_10 <- 1/ x$BF_01

  diag(BF_10) <- 0

  adj_10 <- ifelse(BF_10 > BF_cut, 1, 0)

  adj_01 <- ifelse(x$BF_01 > BF_cut, 1, 0)

  BF_01_adj <- adj_01 * x$BF_01

  BF_10_adj <- adj_10 * BF_10

  pcor_mat <- matrix(0, x$p, x$p)
  if(x$groups == 2){
    pcor_mat[upper.tri(pcor_mat)] <- unlist(x$mu_diff)
    pcor_mat <- symmteric_mat(pcor_mat)
    pcor_mat <- adj_10 * pcor_mat

    }

  returned_object <- list(BF_10 = BF_10,
                          BF_01 = x$BF_01,
                          BF_01_adj = BF_01_adj,
                          BF_10_adj = BF_10_adj,
                          adj_10 = adj_10,
                          adj_01 = adj_01,
                          call = match.call(),
                          p = ncol(BF_10),
                          iter = x$iter,
                          info = x$info,
                          BF = BF_cut,
                          pcor_mat_10 = pcor_mat,
                          object = object)

  class(returned_object) <- "select.ggm_compare_bf"
  returned_object

}


#' @name print.select.ggm_compare_bf
#' @title  Print method for \code{select.ggm_compare_bf} objects
#'
#' @param x An object of class \code{select.ggm_compare_bf}
#' @param ... currently ignored
#' @seealso \code{\link{select.ggm_compare_bf}}
#' @export
print.select.ggm_compare_bf <- function(x,...){
print(x$object)
}

#' @name summary.select.ggm_compare_bf
#' @title  Print method for \code{select.ggm_compare_bf} objects
#'
#' @param object An object of class \code{select.ggm_compare_bf}
#' @param ... currently ignored
#' @seealso \code{\link{select.ggm_compare_bf}}
#' @export
summary.select.ggm_compare_bf <- function(object,...){
  cat("BGGM: Bayesian Gaussian Graphical Models \n")
  cat("--- \n")
  cat("Type: GGM Compare with the Posterior Distribution\n")
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



#' @title Plot \code{select.ggm_compare_bf} Objects
#'
#' @description This function plots the selected graph when comparing GGMs with Bayesian hypothesis
#' testing. There are two heatmap plots in total. The first includes egdes for which there was
#' evidence for a diffirence. The second includes edges for which there is evidence for the null
#' hypothesis of equality. The tiles in the heatmap correspond to the Bayes factor
#'
#'
#' @param x object of class \code{select.ggm_compare_bf}
#' @param H0_low tile color for the smallest Bayes factors in the null hypothesis heatmap
#' @param H0_high tile color for the largest Bayes factors in the null hypothesis heatmap
#' @param H1_low tile color for the smallest Bayes factors in the alternative hypothesis heatmap
#' @param H1_high  tile color for the largest Bayes factors in the alternative hypothesis heatmap
#' @param ... currently ignored
#'
#' @return list containing two \code{ggplot} objects
#' @export
#'
#' @examples
#' # group 1
#' Y1 <- MASS::mvrnorm(500, mu = rep(0,16),
#'                     Sigma =  BGGM::ptsd_cor2,
#'                     empirical = FALSE)
#'
#' # group 2
#' Y2 <- MASS::mvrnorm(500, mu = rep(0,16),
#'                     Sigma =  BGGM::ptsd_cor2,
#'                     empirical = FALSE)
#' # fit model
#' fit <- ggm_compare_bf(Y1, Y2,
#'                      prior_sd = 0.20,
#'                      iter = 1000,
#'                      cores = 2)
#' # select E
#' E <- select(fit, BF_cut = 3)
#'
#' # plot E
#' plot(E)
plot.select.ggm_compare_bf <- function(x, H0_low = "lightblue",
                                       H0_high = "purple",
                                       H1_low = "yellow",
                                       H1_high = "red", ...){

  # make alternative matrix
  alt_df <- reshape::melt(x$BF_10_adj)
  alt_df$BF <- ifelse(alt_df$value == 0, NA, alt_df$value)

  min_scale <- round(log(x$BF),2)

  max_scale <- round(max(log(na.omit(alt_df$BF))),2)

  alt_mat <- ggplot(data = alt_df, aes(x=as.factor(X1),
                                       y = as.factor(X2),
                                       fill=log(BF))) +
    geom_tile() +
    ylim(rev(levels(as.factor(alt_df$X1)))) +
    scale_fill_continuous(low= H1_low,
                          high= H1_high,
                          guide="colorbar",
                          na.value="white",
                          limits = c(min_scale, max_scale),
                          breaks = c(min_scale, max_scale)) +
    ylab("") +
    xlab("Alternative Hypothesis Matrix")

  # make alternative matrix
  null_df <- reshape::melt(x$BF_01_adj)
  null_df$BF <- ifelse(null_df$value == 0, NA, null_df$value)

  min_scale <- round(log(x$BF),2)

  max_scale <- round(max(log(na.omit(null_df$BF))),2)

  null_mat <- ggplot(data = null_df, aes(x=as.factor(X1),
                                         y = as.factor(X2),
                                         fill=log(BF))) +
    geom_tile() +
    ylim(rev(levels(as.factor(alt_df$X1)))) +
    scale_fill_continuous(low = H0_low,
                          high = H0_high,
                          guide="colorbar",
                          na.value="white",
                          limits = c(min_scale, max_scale),
                          breaks = c(min_scale, max_scale)) +
    ylab("") +
    xlab("Null Hypothesis Matrix")

  returned_object <- list(alt_mat = alt_mat,
                          null_mat = null_mat)

  return(returned_object)
}

