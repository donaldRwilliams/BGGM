#' Select Partial Correlattion Differences for \code{ggm_compare_bf} Objects
#'
#' @param object object of class \code{ggm_compare_bf}.
#'
#' @param post_prob numeric. Posterior `inclusion` probability.  The default is set to 0.50.
#'
#' @param ... not currently used.
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
#' \donttest{
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
#' }

select.ggm_compare_bf <- function(object, post_prob = 0.50,...){

  x <- object

  BF_cut = (post_prob) / (1 - post_prob)

  BF_10 <- 1/ x$BF_01

  diag(BF_10) <- 0

  adj_10 <- ifelse(BF_10 > BF_cut, 1, 0)

  adj_01 <- ifelse(x$BF_01 > BF_cut, 1, 0)

  BF_01_adj <- adj_01 * x$BF_01

  BF_10_adj <- adj_10 * BF_10

  pcor_mat <- matrix(0, x$p, x$p)

  if(x$groups == 2){

    pcor_mat <- adj_10 * x$pcor_diff

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
                          post_prob = post_prob,
                          pcor_mat_10 = pcor_mat,
                          object = object)

  class(returned_object) <- c("BGGM",
                              "explore",
                              "select.ggm_compare_bf")
  returned_object

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
#' \donttest{
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
#'                      iter = 50,
#'                      cores = 2)
#' # select E
#' E <- select(fit, BF_cut = 3)
#'
#' # plot E
#' plot(E)
#' }
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

