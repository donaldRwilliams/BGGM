#' @title Extract the Partial Correlation Matrix
#'
#' @description Extract the partial correlation matrix (posterior mean)
#' from \code{\link{estimate}}, \code{\link{explore}}, \code{\link{ggm_compare_estimate}},
#' and \code{\link{ggm_compare_explore}} objects. It is also possible to extract the
#' partial correlation differences for \code{\link{ggm_compare_estimate}} and
#' \code{\link{ggm_compare_explore}} objects.
#'
#' @param object A model estimated with \strong{BGGM}. All classes are supported, assuming
#' there is matrix to be extracted.
#'
#' @param difference Logical. Should the difference be returned (defaults to \code{FALSE}) ? Note
#' that this assumes there is a difference (e.g., an object of class \code{ggm_compare_estimate})
#' and ignored otherwise.
#'
#' @param ... Currently ignored.
#'
#' @return The estimated partial correlation matrix.
#'
#' @examples
#' \donttest{
#' # note: iter = 250 for demonstrative purposes
#'
#' # data
#' Y <- ptsd[,1:5] + 1
#'
#' # ordinal
#' fit <- estimate(Y, type = "ordinal",
#'                 iter = 250)
#'
#' pcor_mat(fit)
#' }
#' @export
pcor_mat <- function(object, difference = FALSE, ...){

  # estimate default
  if(all(c("estimate", "default") %in% class(object))){

    cn <- colnames(object$Y)

    if(is.null(cn)) {
      cn <- 1:object$p
    }

    colnames(object$pcor_mat) <- cn
    rownames(object$pcor_mat) <- cn

    pcor_mat <- round(object$pcor_mat, 3)

    # return
    pcor_mat

    # explore default
  } else if(all(c("explore", "default") %in% class(object))){

    cn <- colnames(object$Y)

    if(is.null(cn)) {
      cn <- 1:object$p
    }

    colnames(object$pcor_mat) <- cn
    rownames(object$pcor_mat) <- cn

    pcor_mat <- round(object$pcor_mat, 3)

    pcor_mat

    # ggm compare estimate
  } else if(is(object, "ggm_compare_estimate")){

    # analytic is false
    if(isFALSE(object$analytic)){

      cn <- colnames(object$post_samp$Y)
      comparisons <- length(object$pcor_mats)

      if(is.null(cn)) {
        cn <- 1:object$p
      }

      # difference ?
      if(isTRUE(difference)){

        # name matrices rows and columns
        for(i in seq_len(comparisons)){

          colnames(object$pcor_mats[[i]]) <- cn
          rownames(object$pcor_mats[[i]]) <- cn

        }

        pcor_mat <- object$pcor_mats


      } else {

        pcor_mat <- list()

        for(i in seq_len(comparisons)){

          pcor_mat[[i]] <-  round(object$post_samp[[i]]$pcor_mat, 3)

        }

        # name for clarity
        names(pcor_mat) <- paste0("Y_g", seq_len(comparisons))

      }

      # return
      pcor_mat

      # analytic
    } else {

      # difference
      if(isTRUE(difference)){

        # names
        cn <- colnames(object$info$dat[[1]])

        # comparisons
        comparisons <- length(object$diff)

        # name matrices rows and columns
        for(i in seq_len(comparisons)){
          colnames(object$diff[[i]]) <- cn
          rownames(object$diff[[i]]) <- cn
        }

        pcor_mat <- lapply(seq_len(comparisons), function(x) {

          round(object$diff[[x]], 3)

          })

        names(pcor_mat) <- names(object$diff)

        pcor_mat

        # no difference
      } else {

        # groups
        groups <- length(object$info$dat)

        pcor_mat <- lapply(seq_len(groups), function(x) {

          round(analytic_solve( object$info$dat[[x]])$pcor_mat, 3)

          })

        # names
        names(pcor_mat) <- paste0("Y_g", seq_len(groups))

        # return
        pcor_mat
      }
    } # end analytic

  } else if(is(object, "ggm_compare_explore")){

    cn <- colnames(object$info$dat[[1]])

    if(is.null(cn)){
      cn <- 1:object$p
    }

    if(isTRUE(difference)){

      if(object$groups > 2){

        stop("difference only available with two groups. see 'estimate'.")

      }

      # pcor mat
      pcor_mat <- round(object$pcor_diff, 3)

      # names
      colnames(pcor_mat) <- cn
      rownames(pcor_mat) <- cn

      # return
      pcor_mat

    } else {

      # groups
      groups <- length(object$info$dat)

      pcor_mat <- lapply(seq_len(groups), function(x) {

        round(analytic_solve( object$info$dat[[x]])$pcor_mat, 3)

        })

      # names
      names(pcor_mat) <- paste0("Y_g", seq_len(groups))

      # pcor mat
      pcor_mat

    }

  } else {
    stop("partial correlation matrix not found.")
  }

}
