#' Maximum A Posteriori Precision Matrix
#'
#'
#' @param Y Matrix (or data frame) of dimensions \emph{n} (observations) by  \emph{p} (variables).
#'
#' @return An object of class \code{map}, including the precision matrix,
#'         partial correlation matrix, and regression parameters.
#'
#' @export
#'
#' @examples
#' Y <- BGGM::bfi[, 1:5]
#'
#' # map
#' map <- map(Y)
#' map
map <- function(Y){

  Y <- na.omit(Y)

  p <- ncol(Y)

  fit <- analytic_solve(Y)

  map <- fit$inv_map

  pcor <- fit$pcor_mat

  betas <- lapply(1:p, function(z) -1 * (map[z,-z] / map[z,z]) )
  betas <- do.call(rbind, betas)

  returned_object <- list(precision = round(map, 3),
                          pcor = round(pcor, 3),
                          betas = betas,
                          dat = Y)

  class(returned_object) <- c("BGGM",
                              "estimate",
                              "map")

  return(returned_object)
}

