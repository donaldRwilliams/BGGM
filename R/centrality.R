#' Bayesian Centrality
#'
#' @description Compute Bayesian centrality. Currently the only option is for strength (see note for details).
#' For each variable, this provides a posterior distribution for the sum of the respective partial
#' correlations. Testing strength (e.g., differences or confirmatory testing of order contraints) is
#' implemented in \code{test.centrality}.
#'
#' @name centrality.BGGM
#'
#' @param object an object of class \code{estimate} or \code{explore}. All data types (e.g., ordinal, binary, etc.)
#'               are supported.
#'
#' @param type character string. Which type of centrality? Currently the only options are
#'             \code{strength} and \code{bridge_strength}. See note for further details.
#'
#' @param cluster charcter string. Which clusters or communities do the nodes belong to? This is
#'                required when \code{type = bridge_stength}.
#'
#' @param select logical. Should strength be condtional on the selected model (default set to \code{NULL})? Note
#'                 that is might seem desirable to compute strength from only the selected edges. However, this
#'                 presents challenges for making statistical inference. See note for further details.
#'
#' @param cred numeric. Credible interval width used for selecting the graph with an \code{estimate}
#'             object (default is \code{0.95}).
#'
#' @param BF_cut numeric. Evidentiary threshold used for selected the graph with an \code{explore}
#'               object (default is \code{3})
#'
#'
#'@references
#'\insertAllCited{}
#'
#' @note Centrality is a new feature to \strong{BGGM} (> 1.0). At this time, only centrality indices based
#' on the partial correlations (e.g., strength) are implemented and this will not likely change
#' \insertCite{@see critiques on centrality in @bringmann2019centrality; @hallquist2019problems}{BGGM}. This is because
#' strength is \emph{merely} the sum of partial correlations, which have a known (approximate) distribution
#' (\href{https://en.wikipedia.org/wiki/Partial_correlation}{Wikipedia}). This
#' translates into the sum also having those same properties. This is ideal for making statistical inference. In fact,
#' it is not uncommon to test the sum of coefficients with classical methods
#' \insertCite{@page 96, section 4.5.2, in  @baum2006introduction;textual}{BGGM}.
#'
#' We are open to all suggestions for implementing strength related measures: perhaps there is a particular sum
#' of theorectial and/or substantative interest.
#'
#' \strong{Beware of Model Selection Bias:}
#'
#' The default setting is \code{select = FALSE}. This means that \emph{all} of the partial correlation for a given
#' node are summed and not only those that were determined to be non-zero. Although it may seem that including only
#' the selected edges leads to desirable "inference," this is not the case. This would create a problem known
#' as model selection bias that arises when \emph{conditioning} on the selected model \insertCite{@leeb2006can}{BGGM} .
#'
#' To make sense of this, note that:
#'
#' \enumerate{
#' \item To be selected, the partial correlations must far enough way from zero.
#' \item As a result of 1, this means that edges will often be upward biased when they are only considered
#'          when they are selected \insertCite{@more specifically a truncated sampling distribution, p. 5 in  @leeb2006can}{BGGM}
#' \item This then can translate into strength also having issues when only computing it from the selected edges.
#' }
#'
#' Overcoming this issue is an active area of research that is often called 'selective' or 'post-selection' inference
#' (see the references in \href{http://bactra.org/notebooks/post-model-selection-inference.html}{weblink}).
#' However, a viable and simple approach is to \strong{not condition} strength on the selected model (i.e., \code{select = FALSE}).
#' This then readily allows for comparing the sum of partial correlations, say, differences in bridge strength. In future
#' versions (> 2.0), the option \code{select = FALSE} may be removed altogether due to the innate proclivity to compare estimates
#' visually (\code{plot.centrality}).
#'
#'
#'
#' @return object of class \code{BGGM} and \code{centrality}
#' @export
centrality.BGGM <- function(object,
                            type = "strength",
                            cluster = NULL,
                            select = FALSE,
                            cred = NULL,
                            BF_cut = NULL){

  # object from estimate methods
  if(is(object, "estimate")){

    # no model selection (preferred)
    if(isFALSE(select)){

      p <- object$p

      pcor_samples <- object$posterior_samples[,  grep("pcors", colnames(object$posterior_samples))]

      # names of partrials in matrix
      mat_names <- matrix(0, nrow = p, ncol = p)
      mat_names[] <- colnames(pcor_samples)

      # compute bayesian strength
      strength <-  lapply(1:p, FUN = function(x) {
        edges <- pcor_samples[,mat_names[x,][-x]]
        rowSums(edges *  apply(edges, 2, sign))
        })


      # warnging: conditional on selected model !
  } else {

   warning(paste("\n`select = TRUE`. please do not attempt to make",
                 "statistical infernece from the `selected` edges",
                 "(e.g., noting differences). See note for futher details.", sep = "\n"))


    object <- select(object, cred = cred)

    # selected edges
    adj <- object$adjacency_non_zero
    diag(adj) <- 0
    p <- ncol(adj)

    mat_names <- matrix(0, nrow = p, ncol = p)
    mat_names[] <- colnames(object$pcor_sample)

    strength <-  lapply(1:p, FUN = function(x) {
      edges <- object$pcor_samples[,mat_names[x, adj[x,] == 1 ]]
      rowSums(edges *  apply(edges, 2, sign))
      })
    } # end of selection

    # end selecte.estimate
  } else {

    stop("object not supported. must be of class 'select.estimate' or 'select.explore'")
  }

  returned_object <- list(strength = strength, object = object)
  class(returned_object) <- c("BGGM", "centrality")
  return(returned_object)
}






