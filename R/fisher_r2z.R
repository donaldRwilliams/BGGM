#' Fisher Z Transformation
#'
#' @description Tranform correlations to Fisher's Z
#' @param r correlation (can be a vector)
#'
#' @return  Fisher Z transformed correlation(s)
#' @export
#'
#' @examples
#' fisher_r2z(0.5)
fisher_r_to_z <- function(r){
  z <- fisher_z(r)
  return(z)
}


