#' Maximum A Posteriori Precision Matrix
#'
#' @param Y data matrix (\emph{n} by  \emph{p})
#' @return object of class \code{map}
#' @export
#'
#' @examples
#' # p = 20
#' Y <- BGGM::bfi[, 1:5]
#'
#' # map
#' map <- map(Y)
#' map
map <- function(Y){

  x <- na.omit(Y)
  x <- scale(x, scale = FALSE)
  p <- ncol(x)
  n <- nrow(x)

  S <- t(x) %*% x

  map <- (n - p - 1) * solve(S)

  pcor <- -1 * (cov2cor(map) - diag(p))

  betas <- lapply(1:p, function(z) -1 * (map[z,-z] / map[z,z]) )
  betas <- do.call(rbind, betas)

  returned_object <- list(precision = map,
                          pcor = pcor,
                          betas = betas,
                          dat = x)

  class(returned_object) <- "map"

  return(returned_object)
}


#' Print Method for \code{map} Objects
#'
#' @param x object of class map
#' @param ... currently ignored
#' @export
print.map <- function(x,...){
  cat("BGGM: Bayesian Gaussian Graphical Models \n")
  cat("--- \n")
  cat("Method: Maximum A Posteriori")
  print(x$pcor)
}
