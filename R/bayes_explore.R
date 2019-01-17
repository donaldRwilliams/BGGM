#' bayes_explore
#'
#' @param X
#' @param exhaustive
#' @param adjustment
#' @param sample
#' @param prior
#'
#' @return
#' @export
#'
#' @examples
bayes_explore <- function(X,
                          exhaustive = FALSE,
                          adjustment = NULL,
                          sample = NULL,
                          prior = NULL){

  # ensure input is a matrix
  X <- na.omit(as.matrix(X))

  # dimensions of data
  n = dim(X)[1]
  p <- dim(X)[2]

  # temporary storage for pcor name
  # (fill then select only the upper triangle)
  mat_name_temp <- matrix(0, p, p)
  mat_name_temp[] <- unlist(lapply(colnames(X), function(x) paste(x, colnames(X), sep = "_")))

  # stop: high dimensional data not allowed
  if(p >= n) stop("The number of observations (n) must be greater than the number of variables (p)")

  # stop: check for column names
  if(is.null(colnames(X))) stop("Variables must have names")

  #if(p < 10 & sample == FALSE) stop("The normal apporx requires more than 10 variables. Set sample = TRUE")

  # scatter matrix
  S = t(X) %*% X

  # precision matrix
  precision_mat <- (n + p) * solve(S + diag(p))

  # partial correlation matrix ("-1" reverse the direction)
  partial_mat <- diag(1/sqrt(diag( precision_mat ))) %*%  precision_mat  %*% diag(1/sqrt(diag( precision_mat ))) * - 1

  # partial correlations in the upper triangle
  pcors <- partial_mat[upper.tri(partial_mat)]

  # approximate posterior standard deviation
  pcors_sd <- approx_sd(pcors, n, p - 1)

  # Bayes factor for null hypothesis
  BF_null <- dnorm(0, pcors, pcors_sd) / dnorm(0, 0, p^(-1/2))

  if(is.null(adjustment)){

    # matrix for storage
    BF_01 <- matrix(0, p, p)

    # Bayes factors for null hypothesis into matrix
    BF_01[upper.tri(BF_01)] <- BF_null
    BF_01 <- symmteric_mat(BF_01)

    # returned object
    results <- list(BF_01 = BF_01,
                    exhaustive = FALSE,
                    p = p,
                    dat = X,
                    pcors = pcors,
                    pcors_sd = pcors_sd,
                    partial_mat = partial_mat,
                    adjustment = FALSE)
    }
  # adjusted bayes factor
  if(!is.null(adjustment)){
    if(!is.numeric(adjustment)) stop("adjustment factor must be numeric")

    # matrix for storage
    BF_01 <- matrix(0, p, p)

    # Bayes factors for null hypothesis into matrix
    BF_01[upper.tri(BF_01)] <- BF_null * (dnorm(0, 0, p^(-1/2)) / adjustment)
    BF_01 <- symmteric_mat(BF_01)

    # returned object
    results <- list(BF_01 = BF_01,
                    exhaustive = FALSE,
                    p = p,
                    dat = X,
                    pcors = pcors,
                    pcors_sd = pcors_sd,
                    partial_mat = partial_mat,
                    adjustment = TRUE)
    }

  # exhaustive approach
  if(exhaustive == TRUE){

    # bayes factor: negative
    BF_negative <- negative_helper(pcors, pcors_sd, 1 / BF_null)

    # bayes factor: positive
    BF_positive <- positive_helper(pcors, pcors_sd, 1 / BF_null)

    # exhaustive results
    exhaustive_results <- as.data.frame(cbind(pcors, t(mapply(exhaustive_helper, BF_null, BF_positive, BF_negative))))

    # name columns
    colnames(exhaustive_results) <- c("pcor", "null_prob", "positive_prob", "negative_prob")

    # name rows
    row.names(exhaustive_results) <- mat_name_temp[upper.tri(mat_name_temp)]

    # returned object
    results <- list(exhaustive_results = exhaustive_results,
                    exhaustive = TRUE,
                    p = p,
                    dat = X,
                    partial_mat = partial_mat)
  }

  class(results) <- "Bayes_explore_object"
  results

}

