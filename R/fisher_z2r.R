#' Fisher Z Back Transformation
#' @description Back tranform Fisher's Z to correlations
#' @param z Fisher Z
#'
#' @return Correlation (s) (backtransformed)
#' @export
#'
#' @examples
#' fisher_r2z(0.5)
fisher_z_to_r <- function(z){
  r <- z2r(z)
  return(r)
}
