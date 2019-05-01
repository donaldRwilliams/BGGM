#' Title
#'
#' @param x
#' @param ci_width
#'
#' @return
#' @export
#'
#' @examples
compare.predict <- function(x, ci_width){

  pairwise_comb <- t(combn(1:x$p, 2))
  node_names <- names(x$post_samples)
  diff <- list()
  names_temp <- list()
  for(i in 1:nrow(pairwise_comb)){
    # paste()
    temp <-  x$post_samples[[node_names[pairwise_comb[i,1]]]] - x$post_samples[[node_names[pairwise_comb[i,2]]]]
    names_temp[[i]]  <- paste(which(names( x$post_samples) == node_names[pairwise_comb[i,1]]),
                              which(names( x$post_samples) == node_names[pairwise_comb[i,2]]), sep = " - " )
    diff[[i]] <- temp

  }


  names(diff) <- names_temp


  res <- lapply(1:nrow(pairwise_comb), function(x)   BGGM:::compare_predict_helper(diff[[x]] , ci_width)       )


  summ <- cbind.data.frame(contrast = names(diff), do.call(rbind.data.frame, res))

  retuned_object <- list(summary_error = summ, test_data = x$test_data,  measure = x$measure, call = match.call())

  class(retuned_object) <- "compare.predict"

  return(retuned_object)

}
