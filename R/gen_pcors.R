#' Simulate a Partial Correlation Matrix
#'
#' @param p number of variables (nodes)
#'
#' @param edge_prob connectivity
#'
#' @param lb lower bound for the partial correlations
#'
#' @param ub upper bound for the partial correlations
#'
#' @note The function checks for a valid matrix (positive definite),
#' but sometimes this will still fail. For example, for
#' larger \code{p}, to have large partial correlations this
#' requires a sparse GGM
#' (accomplished by setting \code{edge_prob}
#' to a small value).
#'
#' @return A list containing the following:
#'
#' \itemize{
#'
#' \item{\strong{pcor}}: Partial correlation matrix, encoding
#' the conditional (in)dependence structure.
#'
#' \item{\strong{cors}}: Correlation matrix.
#'
#' \item{\strong{adj}}: Adjacency matrix.
#'
#' \item{\strong{trys}}: Number of attempts to obtain a
#' positive definite matrix.
#'
#' }
#'
#' @export
#'
#' @importFrom stats runif
#'
#' @examples
#'
#' true_net <- gen_net(p = 10)
gen_net <- function(p = 20,
                    edge_prob = 0.3,
                    lb = 0.05,
                    ub = 0.3) {

  # negative determinant
  d <- -1

  # number of trys
  trys <- 0

  # until d is positive
  while (d < 0) {

    trys <- trys + 1

    effects <- p * (p - 1) * 0.5

    mat <- matrix(1, p, p)

    prob_zero <- 1 - edge_prob

    pool <- c(rep(0, effects * prob_zero),
              runif(effects * edge_prob, lb, ub))

    if (length(pool) != effects) {
      pool <- c(0, pool)
    }

    mat[upper.tri(mat)] <- sample(pool, size = effects)

    pcs <- symm_mat(mat)

    pcs <- -pcs

    diag(pcs) <- -diag(pcs)

    d <- det(pcs)

  }

  cors <- cov2cor(solve(pcs))

  inv <- solve(cors)

  pcors <- cov2cor(inv) * -1

  diag(pcors) <- 1

  adj <- ifelse(pcs == 0, 0, 1)

  returned_object <- list(
    pcors = pcors * adj,
    cors = cors,
    adj = adj,
    trys = trys
  )

  return(returned_object)

}
