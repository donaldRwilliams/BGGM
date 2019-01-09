#' Title
#'
#' @param x
#' @param node
#' @param ci_width
#' @param samples
#'
#' @return
#' @export
#'
#' @examples
beta_summary <- function(x, node, ci_width, samples){
  # x: object "bayes_estimate"
  # node: which nodes to summarize
  # ci_width: credible interval probability
  # samples: number of posterior samples

  # convert inverse to beta
  x <- BGGM:::inverse_2_beta(x, samples = samples)

  # stop if not the correct class
  if(class(x) != "inverse_2_beta"){
    stop("class must be inverse_2_beta")
  }

  # check ci_width is allowed
  if(ci_width >= 1 | ci_width <= 0){
    stop("ci_width must be between 0 and 1")
  }
  returned_object <- lapply(node, function(y) summary_beta_helper(node =  y,
                                                                  x = x,
                                                                  ci_width))
  returned_object
}
