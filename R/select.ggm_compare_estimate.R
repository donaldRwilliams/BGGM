#'  Edge Differences and (Practical) Equivalence Between GGMs
#' @name select.ggm_compare_estimate
#' @param object object of class \code{ggm_compare_estimate}
#' @param cred credible interval width used for the decision rule
#' @param rope region of practical equivalence
#' @param prob posterior probability (see notes)
#' @param ... not currently used
#'
#' @note \code{prob} (posterior probability) is the decision rule for the rope. For example, with \code{rope = 0.1}
#' and \code{prob = 0.95}, differences require that 95 \% of the posterior exlcudes +- 0.1,
#' whereas equivalence requires that 95 \% of the posterior is within  +- 0.1.
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
                                        prob = 0.95, ...) {

  # number of contrasts
  contrasts <- nrow(object$info$pairwise)

  # matrices for storage
  mats <- replicate(contrasts, expr =   matrix(0,  object$p, object$p))

  # no rope
  if(is.null(rope)){
    message("no rope specified. prob ignored.")

    # summary
    summ <- summary(object, cred = cred)

    # storage matrices
    mat_adj <- mat_pcor <- list()

    for(i in 1:contrasts){
      # temp matrix
      mat_temp_adj <-  mat_temp_pcor <- matrix(0, object$p, object$p)

      # select
      mat_temp_adj[upper.tri(mat_temp_adj)]   <-
        ifelse(summ$dat_results[[i]][,4] < 0 &
                 summ$dat_results[[i]][,5] > 0,
               0,
               1)

      # difference mean
      mat_temp_pcor[upper.tri(mat_temp_pcor)] <-
        summ$dat_results[[i]]$Post.mean

      # selected pcors
      mat_temp_pcor <- mat_temp_pcor * mat_temp_adj

      # symmetric matrix
      mat_temp_adj <- symmteric_mat(mat_temp_adj)

      #  column and row names
      colnames(mat_temp_adj) <- 1:object$p
      row.names(mat_temp_adj) <- 1:object$p
      colnames(mat_temp_pcor) <- 1:object$p
      row.names(mat_temp_pcor) <- 1:object$p

      # symmetric matrix
      mat_temp_pcor <- symmteric_mat(mat_temp_pcor)

      # store in list
      mat_adj[[i]] <- mat_temp_adj
      mat_pcor[[i]] <- mat_temp_pcor
      names(mat_adj)[[i]] <- names(object$pcors_diffs[[i]])
      names(mat_pcor)[[i]] <- names(object$pcors_diffs[[i]])
    }

    returned_object <- list(mat_adj = mat_adj,
                            mat_pcor = mat_pcor,
                            call = match.call(),
                            object = object,
                            rope  = rope,
                            cred = cred,
                            prob = prob)

  } else {


    mat_in_adj <- mat_in_pcor <- mat_out_adj <- mat_out_pcor <-  list()


    for(i in 1:contrasts){

      mat_temp_in_adj <-
        mat_temp_pcor <- mat_temp_out_adj <- matrix(0, object$p, object$p)

      in_rope <-
        apply(object$pcors_diffs[[i]][[1]], 2, rope_helper, rope)
      out_rope <- 1 - in_rope


      mat_temp_in_adj[upper.tri(mat_temp_in_adj)] <-
        ifelse(in_rope > prob, 1, 0)

      mat_temp_out_adj[upper.tri(mat_temp_in_adj)] <-
        ifelse(out_rope > prob, 1, 0)

      mat_temp_pcor[upper.tri(mat_temp_pcor)] <-
        colMeans(object$pcors_diffs[[i]][[1]])

      mat_temp_out_pcor <- mat_temp_pcor * mat_temp_out_adj
      mat_temp_out_pcor <- symmteric_mat(mat_temp_out_pcor)
      mat_temp_out_adj  <- symmteric_mat(mat_temp_out_adj)

      mat_temp_in_pcor <- mat_temp_pcor  * mat_temp_in_adj
      mat_temp_in_pcor <- symmteric_mat(mat_temp_in_pcor)
      mat_temp_in_adj  <- symmteric_mat(mat_temp_in_adj)

      # names row and columns (pcors)
      colnames(mat_temp_out_pcor) <- 1:object$p
      row.names(mat_temp_out_pcor) <- 1:object$p
      colnames(mat_temp_in_pcor) <- 1:object$p
      row.names(mat_temp_in_pcor) <- 1:object$p

      # names rows and columns (adj)
      colnames(mat_temp_in_adj) <- 1:object$p
      row.names(mat_temp_in_adj) <- 1:object$p
      colnames(mat_temp_out_adj) <- 1:object$p
      row.names(mat_temp_out_adj) <- 1:object$p


      mat_in_adj[[i]] <- mat_temp_in_adj
      mat_out_adj[[i]] <- mat_temp_out_adj

      mat_in_pcor[[i]] <- round(mat_temp_in_pcor, 3)
      mat_out_pcor[[i]] <- round(mat_temp_out_pcor, 3)

      names(mat_in_adj)[[i]] <- names(object$pcors_diffs[[i]])
      names(mat_out_adj)[[i]] <- names(object$pcors_diffs[[i]])

      names(mat_in_pcor)[[i]] <- names(object$pcors_diffs[[i]])
      names(mat_out_pcor)[[i]] <- names(object$pcors_diffs[[i]])

    }

    returned_object <- list(mat_in_adj = mat_in_adj,
                            mat_out_adj = mat_out_adj,
                            mat_in_pcor = mat_in_pcor,
                            mat_out_pcor = mat_out_pcor,
                            call = match.call(),
                            object = object,
                            cred = cred,
                            rope = rope,
                            prob = prob)


  }

  class(returned_object) <- "select.ggm_compare_estimate"
  return(returned_object)

}

#' @name print.select.ggm_compare_estimate
#' @title  Print method for \code{select.ggm_compare_estimate} objects
#'
#' @param x An object of class \code{select.ggm_compare_estimate}
#' @param ... currently ignored
#' @seealso \code{\link{select.ggm_compare_estimate}}
#' @export
print.select.ggm_compare_estimate <- function(x,...){
  print(x$object)
}


#' @name summary.select.ggm_compare_estimate
#' @title Summary method for \code{select.ggm_compare_estimate} objects
#'
#' @param object An object of class \code{select.ggm_compare_estimate}
#' @param ... currently ignored
#' @seealso \code{\link{select.ggm_compare_estimate}}
#' @export
summary.select.ggm_compare_estimate <- function(object,...){
  cat("BGGM: Bayesian Gaussian Graphical Models \n")
  cat("--- \n")
  cat("Type: GGM Compare with the Posterior Distribution\n")
  # number of iterations
  cat("Posterior Samples:", object$object$iter, "\n")
  if (is.null(object$rope)) {

    cat("Credible Interval:",  gsub("*0.","", formatC( round(object$cred, 4), format='f', digits=2)), "% \n")

  } else {

    cat("Probability:", object$prob, "\n")
    cat("Region of Practical Equivalence:", "[", -1 * object$rope, ", ", object$rope, "]", "\n", sep = "")
  }
  # number of observations
  cat("Observations (n):\n")
  groups <- length(object$object$info$dat)
  for(i in 1:groups){
    cat("  Group", paste( i, ":", sep = "") , object$object$info$dat_info$n[[i]], "\n")
  }
  # number of variables
  cat("Variables (p):", object$object$p, "\n")
  # number of edges
  cat("Edges:", .5 * (object$object$p * (object$object$p-1)), "\n")
  cat("--- \n")
  cat("Call: \n")
  print(object$call)
  cat("--- \n")
  cat("Selected:\n\n")
  cat("Partial Correlations \n\n")
  if(is.null(object$rope)){
    for(i in 1:length(object$mat_adj)){
      cat(names(object$mat_adj)[[i]], "\n")

      print(object$mat_pcor[[i]])
      cat("\n")
    }
  } else {
    for(i in 1:length(object$mat_in_adj)){
      cat(names(object$mat_in_adj)[[i]], "\n")

      print(object$mat_out_pcor[[i]])
      cat("\n")
    }
    cat("--- \n")
    cat("note: null matrices in the select object")
  }
}

