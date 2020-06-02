#' Graph Selection for \code{VAR.estimate} Object
#'
#' @param object An object of class \code{VAR.estimate}.
#'
#' @param cred Numeric. The credible interval width for selecting the graph
#'  (defaults to 0.95; must be between 0 and 1).
#'
#' @param alternative A character string specifying the alternative hypothesis. It
#'                    must be one of "two.sided" (default), "greater"  or "less".
#'                    See note for futher details.
#'
#' @return An object of class \code{select.VAR_estimate}
#'
#' @export
select.VAR_estimate <- function(object,
                                cred = 0.95,
                                alternative = "two.sided"){
}
