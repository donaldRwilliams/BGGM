#' @title GGM Select: Estimation
#'
#' @description Provides the selected graph based on credible intervals for the partial correlations
#' that did not contain of zero \inserCite{Williams2019}{BGGM}.
#'
#' @name select.estimate
#'
#' @param object An object of class \code{estimate.default}.
#'
#' @param cred Numeric. The credible interval width for selecting the graph
#'  (defaults to 0.95; must be between 0 and 1).
#'
#' @param alternative A character string specifying the alternative hypothesis. It
#'                    must be one of "two.sided" (default), "greater"  or "less".
#'                    See note for futher details.

#' @param ... Not currently used.
#'
#'
#' @return The returned object of class \code{select.estimate} contains a lot of information that
#'         is used for printing and plotting the results. For users of \strong{BGGM}, the following
#'         are the useful objects:
#'
#' \itemize{
#' \item \code{pcor_adj} Selected partial correlation matrix (weighted adjacency).
#' \item \code{adj} Adjacency matrix for the selected edges
#' \item \code{object} An object of class \code{estimate} (the fitted model).
#'
#' }
#'
#' @seealso estimate ggm_compare_estimate
#'
#'
#' @details
#'
#' This package was built for the social-behavioral sciences in particular. In these applications, there is
#' strong theory that expects \emph{all} effects to be positive. This is known as a "positive manifold" and
#' this notion has a rich tradition in psychometrics. Hence, this can be incorportated into the graph with
#' \code{alternative = "greater"}. This results in the estimted structure including only positive edges. Further
#' details can be found at the blog "Dealing with Negative (Red) Edges in Psychological Networks: Frequentist Edition"
#' (\href{https://donaldrwilliams.github.io/2020/03/29/dealing-with-negative-red-edges-in-psychological-networks-frequentist-edition/}{link})
#'
#'
#' @examples
#' \donttest{
#' # data
#' Y <- bfi[,1:25]
#'
#' # estimate
#' fit <- estimate(Y, iter = 250)
#'
#'
#' # select edge set
#' E <- select(fit)
#'
#' }
#'
#' @export
select.estimate <- function(object,
                            cred = 0.95,
                            alternative = "two.sided",
                            ...){

  if(isFALSE(object$analytic)){

    pcors <- object$post_samp$pcors[,,51:(object$iter +50)]

    if(alternative == "two.sided"){

      lb <- (1 - cred) / 2
      ub <- 1 - lb

      adj <- ifelse(apply(pcors, 1:2, quantile, lb) < 0 &
                    apply(pcors, 1:2, quantile, ub) > 0, 0, 1)

      } else if(alternative == "greater") {

        lb <- (1 - cred)
        adj <- ifelse(apply(pcors, 1:2, quantile, lb) > 0, 1, 0)

        } else {

          ub <- cred
          adj <- ifelse(apply(pcors, 1:2, quantile, ub) < 0, 1, 0)

          }

    # analytic
    } else {

      if(alternative == "two.sided"){

        lb <- (1 - cred) / 2
        ub <- 1 - lb

        z_stat <- abs(object$analytic_fit$inv_map /  sqrt(object$analytic_fit$inv_var))

        adj <- ifelse(z_stat >  qnorm(ub), 1, 0)

      } else if (alternative == "greater") {

        ub <- 1 - cred

        z_stat <- (-object$analytic_fit$inv_map) /  sqrt(object$analytic_fit$inv_var)

        adj <- ifelse(z_stat > qnorm(ub, lower.tail = FALSE), 1, 0)



      } else if(alternative == "less"){


        ub <- 1 - cred

        z_stat <- (object$analytic_fit$inv_map) /  sqrt(object$analytic_fit$inv_var)

        adj <- ifelse(z_stat > qnorm(ub, lower.tail = FALSE), 1, 0)

      }

      }

  pcor_adj <- adj * object$pcor_mat

  returned_object <- list(
    pcor_adj = pcor_adj,
    adj = adj,
    alternative = alternative,
    cred = cred,
    object = object
  )

  class(returned_object) <- c("BGGM",
                              "select.estimate",
                              "estimate",
                              "select")
  returned_object
}

#' @title S3 \code{select} method
#' @name select
#' @description S3 select method
#' @param object object of class \code{estimate} or\code{explore}
#' @param ... not currently used
#' @return \code{select} works with the following methods:
#' \itemize{
#' \item \code{\link{select.estimate}}
#' \item \code{\link{select.explore}}
#' \item \code{\link{select.ggm_compare_estimate}}
#' }
#' @export
select <- function(object,...){
  UseMethod("select", object)
}


print_select_estimate <- function(x, ...){
  object <- x
  p <- ncol(object$pcor_adj)
  cat("BGGM: Bayesian Gaussian Graphical Models \n")
  cat("--- \n")
  cat("Type:", object$object$type, "\n")
  cat("Analytic:", object$object$analytic, "\n")
  cat("Formula:", paste(as.character(x$formula), collapse = " "), "\n")
  cat("Posterior Samples:", object$object$iter, "\n")
  cat("Credible Interval:",  gsub("*0.","", formatC( round(object$cred, 4), format='f', digits=2)), "% \n")
  cat("--- \n")
  cat("Call: \n")
  print(object$object$call)
  cat("--- \n")
  cat("Selected:\n\n")
  mat <- object$pcor_adj
  colnames(mat) <- 1:p
  row.names(mat) <- 1:p
  print(round(mat, 3))
  cat("--- \n")
}
