#' GGM Selection: Estimation
#' @name select.estimate
#'
#' @description Select the edge set for \code{estimate} objects. The graph is determined with
#' the posterior distribution, in particular credible interval exclusion of zero. This
#' corresponds to a directional posterior probability \insertCite{Marsman2017a}{BGGM}. For example,
#' the probability (conditional on the model) of a positive edge is at least 97.5% when the 95%
#' credible interval excludes zero.
#'
#' @param object object of class \code{estimate.default}
#' @param cred credible interval width used for the decision rule
#'
#' @param alternative a character string specifying the alternative hypothesis,
#'                    must be one of "two.sided", "greater" (default) or "less".
#'                    See note for futher details.

#' @param ... not currently used
#'
#' @references
#' \insertAllCited{}
#'
#' @return An object of class \code{select.estimate}:
#'
#' \itemize{
#' \item \code{pcor_adj} selected partial correlation matrix
#' \item \code{adj} adjacency matrix for the selected edges
#' \item \code{ci} credible interval width
#' \item \code{object} object of class \code{estimate}
#' }
#'

#'
#'
#' @note
#'
#' \strong{Alternative Default:}
#'
#' This package was built for the social-behavioral sciences in particular. In these applications, there is
#' strong theory that expects \emph{all} effects to be positive. This is known as a "positive manifold" and
#' this notion has a rich tradition in psychometrics. Hence, by default this expecation is included into graph selection
#' (\code{alternative = "greater"}). This results in the estimted structure including only positive edges. Further
#' details can be found at the blog "Dealing with Negative (Red) Edges in Psychological Networks: Frequentist Edition"
#' (\href{https://donaldrwilliams.github.io/2020/03/29/dealing-with-negative-red-edges-in-psychological-networks-frequentist-edition/}{link})
#'
#' @examples
#'
#' # Analytic = TRUE
#'# p = 5
#' Y <- BGGM::bfi[,1:5]
#'
#' # analytic solution
#' fit_analytic <- estimate(Y, analytic = TRUE)
#'
#' # select E
#' E <- select(fit_analytic, ci_width = 0.95)
#'
#' # non-zero partial correlations
#' E$partials_non_zero
#'
#' # adjacency matrix
#' E$adjacency_non_zero
#' @export
select.estimate <- function(object,
                            cred = 0.95,
                            alternative = "greater"){

  if(isFALSE(object$analytic)){

    pcors <- object$post_samp$pcors[,,51:(fit$iter +50)]

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

        adj <- ifelse( z_stat > qnorm(ub, lower.tail = FALSE), 1, 0)



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

#' @title S3 select method
#' @name select
#' @description S3 select method
#' @param object object of class \code{estimate}, \code{explore}, or ..
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
  cat("Formula:", paste(as.character(fit$formula), collapse = " "), "\n")
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
