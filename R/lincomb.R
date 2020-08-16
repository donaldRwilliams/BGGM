#' Test a linear combination of posterior samples
#'
#' @param lin_comb A string specifying a linear combination of variables
#'
#' @param object An object of class \code{BGGM}
#'
#' @param cred Numeric. The credible interval width
#'                  (defaults to 0.95; must be between 0 and 1).
#'
#' @param rope Specify a region of practical equivalence (i.e., ROPE). This is Optional.
#'
#' @return An object of class \code{lin_comb}
#'
#' @examples
#' \donttest{
#' # data
#' Y <- BGGM::ptsd
#'
#' # names
#' colnames(Y) <- letters[1:20]
#'
#' # estimate model
#' est <- BGGM::estimate(Y)
#'
#' # test
#' bggm_comb <- lincomb("a--c + a--d > b--c + b--d",
#'                        object = est,
#'                        cred = 0.90)
#'
#' # print
#' bggm_comb
#' }
#' @export
lincomb <- function(lin_comb, object, cred = 0.95, rope = NULL){

  out <- bayeslincom::lin_comb(lin_comb = lin_comb,
                               obj = object,
                               cri_level = cred,
                               rope = rope)

  class(out) <- c("BGGM", "bayeslincom")
  return(out)

}

print_lincomb <- function(x, ...){
  cat("BGGM: Bayesian Gaussian Graphical Models \n")
  cat("--- \n")
  bayeslincom::print.bayeslincom(x)
}
