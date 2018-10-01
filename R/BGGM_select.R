#' Title
#'
#' @param x BGGM fitted object
#' @param BF_directional Decision rule based on Bayes factors
#' @param probability Decision rule based on Posterior probabiities
#'
#' @return
#' @export
#'
#' @examples
BGGM_select  <- function(x, BF_directional = NULL, probability = NULL){
  pcors <-as.matrix( x$parcors_mat)
  x <- x$post_prob

  if(is.null(BF_directional) & is.null(probability)){
    stop("BF_directional or probability must be specified ")
  }

  if(is.numeric(probability)){

    selected <- ifelse(x < (1-probability), 1, 0)
    diag(selected) <- 0
    list_return <- list(adjacency_mat =   selected, pcors_selected = selected   * pcors)
  }

  if(is.numeric(BF_directional)){

    selected <- ifelse(x > BF_directional | x < (1 / BF_directional), 1, 0)
    list_return <- list(adjacency_mat =   selected, pcors_selected = selected   * pcors)


    }
return(list_return)

}

