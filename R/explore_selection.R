#' Title
#'
#' @param X
#' @param threshold
#' @param type
#'
#' @return
#' @export
#'
#' @examples
explore_selection <- function(X, threshold = 3, type = NULL){

  # check the class
  if(class(X) != "Bayes_explore_object") stop("X must be from Bayes_explore")

  if(X$exhaustive == FALSE & sum(type == c("null_vs_positive", "null_vs_negative")) != 0) {
    stop("These hypothese are not allowed with exhaustive = FALSE")
  }



  if(X$exhaustive == FALSE){
    # ensure type != NULL
    if(is.null(type)) stop("Must specificy hypothesis")

    # ensure type is permitted
    if(sum(type == c("two_sided", "greater_than", "less_than")) == 0) stop("Must be two_sided, greater_than, or less_than")

    # two-sided hyp testing
    if(type == "two_sided"){
      # select based on threshold
      BF_null_adj <- ifelse(X$BF_01 > threshold, 1, 0)
      BF_alt_adj  <-  ifelse(1 / X$BF_01 > threshold, 1, 0)

      # returned object
      results <-  list(BF_null_adj = BF_null_adj,
                       BF_alt_adj = BF_alt_adj,
                       BF_null = X$BF_01,
                       BF_alt = 1 / X$BF_01,
                       partial_mat = X$partial_mat  * BF_alt_adj)
    }

    # one-sided (greater_than) hyp testing
    if(type == "greater_than"){
      BF_pos_adj <- BF_pos <- matrix(0, ncol = X$p, nrow = X$p)

      BF <- positive_helper(pcor = X$pcors, post_sd = X$pcors_sd, BF_null = 1/X$BF_01[upper.tri(X$BF_01)])
      BF_pos[upper.tri(BF_pos)] <- BF
      BF_pos <- symmteric_mat(BF_pos)

      BF_pos_adj <- ifelse(BF_pos > threshold, 1, 0)

      results <- list(BF_pos_adj = BF_pos_adj,
                      BF_pos = BF_pos,
                      partial_mat = X$partial_mat * BF_pos_adj)
    }





  } else {



    if(is.null(type)) stop("Must specificy hypothesis")

    if(type == "null_vs_positive"){
      # matrices for storage
      BF_null_adj <- BF_positive_adj <- BF_null <- BF_positive <- matrix(0, ncol = X$p, nrow = X$p)

      # evidence for null vs positive
      BF_null_adj[upper.tri(BF_null_adj)] <- ifelse(X$exhaustive_results$null_prob / X$exhaustive_results$positive_prob > threshold &
                                                      X$exhaustive_results$null_prob / X$exhaustive_results$negative_prob > threshold, 1, 0)


      # evidnece for postive vs null
      BF_positive_adj[upper.tri(BF_positive_adj)] <-
        # simple ifelse (broken up to reduce width)
        ifelse(1/(X$exhaustive_results$null_prob /
                    X$exhaustive_results$positive_prob) > threshold
               & # check complement
                 X$exhaustive_results$positive_prob /
                 X$exhaustive_results$negative_prob > threshold, 1, 0)



      BF_null_adj <- symmteric_mat(BF_null_adj)
      BF_positive_adj <- symmteric_mat(BF_positive_adj)

      BF_null[upper.tri(BF_null)] <- X$exhaustive_results$null_prob / X$exhaustive_results$positive_prob
      BF_null <- symmteric_mat(BF_null)
      BF_null <- BF_null * BF_null_adj
      BF_positive[upper.tri(BF_positive)] <-  1/(X$exhaustive_results$null_prob / X$exhaustive_results$positive_prob)
      BF_positive <- symmteric_mat(BF_positive)
      BF_positive <- BF_positive * BF_positive_adj

      # name the columns based on the data
      colnames(BF_null_adj) <-
        colnames(BF_positive_adj) <-
        colnames(BF_null) <-
        colnames(BF_positive) <-
        colnames(X$dat)


      # name the columns based on the data
      row.names(BF_null_adj) <-
        row.names(BF_positive_adj) <-
        row.names(BF_null) <-
        row.names(BF_positive) <-
        colnames(X$dat)

      results <- list(BF_null_adj = BF_null_adj,
                      BF_positive_adj = BF_positive_adj,
                      BF_positive = BF_positive,
                      BF_null = BF_null)
    }

    if(type == "null_vs_negative"){
      BF_null_adj <- BF_negative_adj <- BF_null <- BF_negative <- matrix(0, ncol = X$p, nrow = X$p)

      BF_null_adj[upper.tri(BF_null_adj)] <- ifelse(X$exhaustive_results$null_prob / X$exhaustive_results$negative_prob > threshold &
                                                      X$exhaustive_results$null_prob / X$exhaustive_results$positive_prob > threshold, 1, 0)



      BF_negative_adj[upper.tri(BF_negative_adj)] <- ifelse(1 / (X$exhaustive_results$null_prob /
                                                                 X$exhaustive_results$negative_prob) > threshold &
                                                                 X$exhaustive_results$negative_prob /
                                                                 X$exhaustive_results$positive_prob > threshold, 1, 0)



      BF_null_adj <- symmteric_mat(BF_null_adj)
      BF_negative_adj <- symmteric_mat(BF_negative_adj)

      prob_null[upper.tri(prob_null)] <- X$exhaustive_results$null_prob
      prob_null <- symmteric_mat(prob_null)
      prob_null <- prob_null * BF_null_adj
      prob_negative[upper.tri(prob_negative)] <-  X$exhaustive_results$negative_prob
      prob_negative <- symmteric_mat(prob_negative)
      prob_negative <- prob_negative * BF_negative_adj

      colnames(BF_null_adj) <- colnames(BF_negative_adj) <- colnames(prob_null) <- colnames(prob_negative) <-  colnames(X$dat)
      row.names(BF_null_adj) <- row.names(BF_negative_adj) <-  row.names(prob_null) <- row.names(prob_negative) <- colnames(X$dat)

      results <- list(BF_null_adj = BF_null_adj,
                      BF_negative_adj = BF_negative_adj,
                      prob_negative = prob_negative,
                      prob_null = prob_null)
    }




  }
  results

}
