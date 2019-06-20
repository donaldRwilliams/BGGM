#'  Edge Differences and (Practical) Equivalence Between GGMs
#'
#' @param x object of class \code{ggm_compare_estimate}
#' @param type decision rule to use. either \code{"ci"} or \code{"rope"}
#' @param prob posterior probability. only used for \code{type = "rope"}
#'
#' @return list of class \code{select.ggm_compare_estimate}:
#'
#' \code{type = "ci"}:
#'  \itemize{
#'  \item \code{mats_diff} adjacency matrices for edges differences (one for each group contrast)
#'  \item  \code{type} \code{ci} (used internally for summarizing and plotting)
#'  \item  \code{prob} \code{NULL} (used internally for summarizing and plotting)
#'  }
#' \code{type = "rope"}:
#'  \itemize{
#'  \item \code{mats_diff} adjacency matrices for edge differences (one for each group contrast)
#'  \item  \code{mats_null} adjacency matrices for edge (practical) equivalence (one for each group contrast)
#'  \item  \code{type} \code{ci} (used internally for summarizing and plotting)
#'  \item  \code{prob} \code{NULL} (used internally for summarizing and plotting)
#'  }
#' @export
#'
#' @examples
#' # Assume null is true for all edges
#' Y1 <- MASS::mvrnorm(500, rep(0, 16), Sigma = diag(16))
#' Y2 <- MASS::mvrnorm(500, rep(0, 16), Sigma = diag(16))
#' Y3 <- MASS::mvrnorm(500, rep(0, 16), Sigma = diag(16))
#'
#' ggm_est <- ggm_compare_estimate(Y1, Y2, Y3, iter = 10000, rope = .1)
#' sel <- select(ggm_est, type = "rope", prob = .95)
#'
#' # practically equivalent
#' sel$mats_null

select.ggm_compare_estimate <- function(x, type, prob = NULL){

      contrasts <- length(x$dat_results)

      mats <- replicate(contrasts, expr =   matrix(0,  x$p, x$p))

      mats_adj <- list()

    if(type ==  "ci"){

        if(is.numeric(x$rope)) stop("rope cannot be used. please use type = 'rope' or use ci_width instead of the rope when estimating the models")

        mats_adj <- list()

     for(i in 1:contrasts){


          mat_temp <- matrix(0, x$p, x$p)


          mat_temp[upper.tri(mat_temp)] <- ifelse( x$dat_results[[i]][[1]][,4] < 0 & x$dat_results[[i]][[1]][,5] > 0, 0, 1)

          mat_temp <- list(BGGM:::symmteric_mat(mat_temp))

          names(mat_temp) <- names(x$dat_results[[i]])

          mats_adj[[i]] <- mat_temp

        }

        returned_object <- list(mats_diff = mats_adj,
                                type = type,
                                prob = prob,
                                call = match.call(),
                                p = x$p,
                                rope = x$rope,
                                ci = x$ci_width)


        }

      if(type == "rope"){

        if(is.numeric(x$ci_width)) stop("ci_width cannot be used. please use type = 'ci' or use the rope instead of ci_width= when estimating the models")

        mats_adj_null <- mats_adj_diff <- list()


        if(is.null(prob)) stop("probability must be specified for type = 'rope'")

        for(i in 1:contrasts){

          mat_temp <- matrix(0, x$p, x$p)


          mat_temp[upper.tri(mat_temp)] <- ifelse(x$dat_results[[i]][[1]]$pr_out > prob, 1, 0)

          mat_temp <- list(BGGM:::symmteric_mat(mat_temp))

          names(mat_temp) <- names(x$dat_results[[i]])

          mats_adj_diff[[i]] <- mat_temp

          mat_temp <- matrix(0, x$p, x$p)
          mat_temp[upper.tri(mat_temp)] <- ifelse(x$dat_results[[1]][[1]]$pr_in > prob, 1, 0)

          mat_temp <- list(BGGM:::symmteric_mat(mat_temp))

          names(mat_temp) <- names(x$dat_results[[i]])

          mats_adj[[i]] <- mat_temp


          }


         returned_object <- list(mats_null = mats_adj,
                                 mats_diff = mats_adj_diff,
                                 type = type,
                                 prob = prob,
                                 call = match.call(),
                                 p = x$p,
                                 rope = x$rope,
                                 ci = x$ci_width)
 }

  class(returned_object) <- "select.ggm_compare_estimate"

  returned_object

}


