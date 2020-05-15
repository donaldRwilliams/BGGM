#' @title Graph Selection for \code{ggm_compare_estimate} Objects

#' @description Provides the selected graph (of differences) based on credible intervals for
#' the partial correlations that did not contain zero
#' \insertCite{Williams2019}{BGGM}.
#'
#' @name select.estimate
#'
#' @param object An object of class \code{estimate.default}.
#'
#' @param cred Numeric. The credible interval width for selecting the graph
#'  (defaults to 0.95; must be between 0 and 1).
#'
#' @param ... not currently used
#'
#' @return The returned object of class \code{select.ggm_compare_estimate} contains a lot of information that
#'         is used for printing and plotting the results. For users of \strong{BGGM}, the following
#'         are the useful objects:
#'
#'
#' \itemize{
#'
#' \item \code{mean_diff} A list of matrices for each group comparsion (partial correlation differences).
#'
#' \item \code{pcor_adj} A list of weighted adjacency matrices for each group comparsion.
#'
#' \item \code{adj} A list of adjacency matrices for each group comparsion.
#'
#' }
#'
#' @examples
#' \donttest{
#'
#' ##################
#' ### example 1: ###
#' ##################
#' data
#' Y <- bfi
#'
#' # males and females
#' Ymale <- subset(Y, gender == 1,
#'                select = -c(gender,
#'                            education))
#' # fit model
#' fit <- ggm_compare_estimate(Ymale, Yfemale,
#'                            type = "continuous")
#'
#'
#' E <- select(fit)
#'
#' }
#' @export
select.ggm_compare_estimate <- function(object,
                                        cred = 0.95,
                                        ...) {


  # rope removed, but will add to minor realeses
  rope = NULL

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


