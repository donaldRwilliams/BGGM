#'  Select Partial Correlattion Differences for \code{ggm_compare_estimate} Objects
#' @name select.ggm_compare_estimate
#'
#' @param object object of class \code{ggm_compare_estimate}
#'
#' @param cred credible interval width used for the decision rule
#'
#' @param rope region of practical equivalence
#'
#' @param ... not currently used
#'
#' @return A list of class \code{select.ggm_compare_estimate}:
#' \code{rope} = NULL:
#' \itemize{
#' \item \code{mat_adj} adjacency matrix (one for each contrast)
#' \item \code{mat_pcor} selected partial correlations (one for each contrast)
#' \item \code{call} \code{match.call()}
#' \item \code{object} object of class \code{ggm_compare_estimate}
#' \item \code{rope} region of practical equivalence
#' \item \code{cred} credible interval
#' \item \code{prob} posterior probability
#' }
#' @examples
#' # data
#' Y1 <- BGGM::bfi[1:500,1:5]
#' Y2 <- BGGM::bfi[501:1000, 1:5]
#'
#' # fit model
#' fit <- ggm_compare_estimate(Y1, Y2)
#'
#' # posterior summary of differences
#' summary(fit)
#'
#' # select (threshold) with credible intervals
#' sel <- select(fit)
#'
#' # summary
#' summary(sel)
#'
#'# selected differences
#' sel$mat_pcor
#'
#' # adjacency matrix
#' sel$mat_adj
#' @export
select.ggm_compare_estimate <- function(object,
                                        cred = 0.95,
                                        rope = NULL,
                                        ...) {

  # number of contrasts
  contrasts <- nrow(object$info$pairwise)

  if(is.null(rope)){

    if(isFALSE(object$analytic)){

      lb <- (1 - cred) / 2
      ub <- 1 - lb

      mean_diff <- lapply(1:contrasts, function(x){

      apply(object$diff[[x]]  , 1:2, mean)

    })


    adj <- lapply(1:contrasts, function(x){

      ifelse(apply(object$diff[[x]], 1:2, quantile, lb) < 0 &
             apply(object$diff[[x]], 1:2, quantile, ub) > 0, 0, 1)
    })


    pcor_adj <- lapply(1:contrasts, function(x){
      mean_diff[[x]] * adj[[x]]
      })

    returned_object <- list(mean_diff = mean_diff,
                            pcor_adj = pcor_adj,
                            adj = adj,
                            call = match.call(),
                            object = object,
                            rope = rope,
                            cred = cred)
    } else {

      # analytic
      critical <- abs(qnorm((1 - cred) / 2))

      adj <- lapply(1:contrasts, function(x)  {

        ifelse(object$z_stat[[x]] > critical, 1, 0)

        })

      pcor_adj <- lapply(1:contrasts, function(x){

        object$diff[[x]] * adj[[x]]

        })

      returned_object <- list(adj = adj,
                              pcor_adj = pcor_adj,
                              adj = adj,
                              call = match.call(),
                              object = object,
                              rope = rope,
                              cred = cred)

      } # end analytic

    } else { # for rope. future direction

      stop("rope is not currently implemented.")

      }

  class(returned_object) <- c("BGGM",
                              "select.ggm_compare_estimate",
                              "estimate", "select")
  return(returned_object)

}



print_select_ggm_compare_estimate <- function(x,...){
  object <- x
  comparisons <- length(object$pcor_adj)
  p <- ncol(object$pcor_adj[[1]])
  cat("BGGM: Bayesian Gaussian Graphical Models \n")
  cat("--- \n")
  cat("Type:", object$object$type, "\n")
  cat("Analytic:", object$object$analytic, "\n")
  cat("Posterior Samples:", object$object$iter, "\n")
  cat("Credible Interval:",  gsub("*0.","", formatC( round(object$cred, 4), format='f', digits=2)), "% \n")
  cat("--- \n")
  cat("Call: \n")
  print(object$object$call)
  cat("--- \n")
  cat("Selected:\n\n")
  for(i in 1:comparisons){

    cat(names(object$object$diff)[i], "\n")
    mat <- object$pcor_adj[[i]]
    colnames(mat) <- 1:p
    row.names(mat) <- 1:p
    print(round(mat, 3))
    cat("--- \n\n")
  }
}


