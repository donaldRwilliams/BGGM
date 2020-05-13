#' @title Zero-Order Correlations
#'
#'
#' @name zero_order_cors
#'
#' @description Estimate zero-order correlations for any type of data. Note zero-order refers to the fact that
#' no variables are controlled for (i.e., bivariate correlations). To our knowledge, this is the only Bayesian
#' implementation in \code{R} that can estiamte Pearson's,  tetrachoric (binary), polychoric
#' (ordinal with more than two cateogries), and rank based correlation coefficients.
#'
#' @param Y  matrix (or data frame) of dimensions \emph{n} (observations) by  \emph{p} (variables).

#' @param type character string. Which type of data for \strong{Y} ? The options include \code{continuous},
#' \code{binary}, \code{ordinal}, or \code{mixed}. See the note for further details.
#'
#' @param mixed_type numeric vector. An indicator of length p for which varibles should be treated as ranks.
#' (1 for rank and 0 to assume normality). The default is currently to treat all integer variables as ranks
#' when \code{type = "mixed"} and \code{NULL} otherwise. See note for further details.
#'
#' @param iter number of iterations (posterior samples; defaults to 5000).
#'
#'
#' @return
#'
#' \itemize{
#'
#' \item \code{R} An array including the correlation matrices
#'               (of dimensions \emph{p} by \emph{p} by \emph{iter})
#'
#' \item \code{R_mean} Posterior mean of the correlations (of dimensions \emph{p} by \emph{p})
#' }
#'

#' @export
zero_order_cors <- function(Y,  type = "continuous",
                            iter = 5000,
                            mixed_type = NULL){

  fit <- estimate(Y,
                  type = type,
                  iter = iter,
                  mixed_type = mixed_type)

  cors <- pcor_to_cor(fit)

  return(cors)

}


