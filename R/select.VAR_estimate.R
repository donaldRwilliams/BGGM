#' Graph Selection for \code{var.estimate} Object
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
select.var_estimate <- function(object,
                                cred = 0.95,
                                alternative = "two.sided"){


  pcors <- object$fit$pcors[,,51:(object$iter +50)]

  pcor_mat <- apply(pcors, 1:2, mean)

  beta <- object$fit$beta[,,51:(object$iter +50)]

  beta_mat <- apply(beta, 1:2, mean)


  if(alternative == "two.sided"){

    lb <- (1 - cred) / 2

    ub <- 1 - lb

    pcor_adj <- ifelse(apply(pcors, 1:2, quantile, lb) < 0 &
                       apply(pcors, 1:2, quantile, ub) > 0, 0, 1)


    beta_adj <- ifelse(apply(beta, 1:2, quantile, lb) < 0 &
                       apply(beta, 1:2, quantile, ub) > 0, 0, 1)

  } else if(alternative == "greater") {

    lb <- (1 - cred)

    pcor_adj <- ifelse(apply(pcors, 1:2, quantile, lb) > 0, 1, 0)

    beta_adj <- ifelse(apply(beta, 1:2, quantile, lb) > 0, 1, 0)


  } else {

    ub <- cred

    pcor_adj <- ifelse(apply(pcors, 1:2, quantile, ub) < 0, 1, 0)

    beta_adj <- ifelse(apply(beta, 1:2, quantile, ub) < 0, 1, 0)

  }

  beta_weighted_adj <-  beta_adj * beta_mat
  pcor_weighted_adj <- pcor_adj * pcor_mat


  returned_object <- list(
    pcor_adj = pcor_adj,
    beta_adj = beta_adj,
    beta_weighted_adj =  beta_weighted_adj,
    pcor_weighted_adj =  pcor_weighted_adj,
    beta_mu = beta_mat,
    pcor_mu = pcor_mat,
    alternative = alternative,
    cred = cred,
    object = object
  )

  class(returned_object) <- c("BGGM",
                              "select.var_estimate",
                              "VAR_estimate",
                              "select")
  return(returned_object)

}
