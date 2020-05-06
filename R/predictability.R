#' Predictability: Bayesian Variance Explained (R2)
#'
#' @name predictability
#'
#' @description  Compute nodewise predictability or  Bayesian variance explained \insertCite{@R2 @gelman_r2_2019}{BGGM}.
#'               In the context of GGMs, this method was described in \insertCite{Williams2019;textual}{BGGM}.
#'               Currently, continuous and mixed data are supported.
#'
#'
#' @param object object of class \code{estimate}
#'
#' @param select logical. Should the graph be selected ? The default is currently \code{FALSE}.
#'
#' @param cred numeric. credible interval between 0 and 1  (default is 0.95) that is used for selecting the graph.
#'
#' @param iter interger. iterations (posterior samples) used for computing R2.
#'
#' @param ... currently ignored.
#'
#' @return object of classes \code{bayes_R2} and \code{metric}
#'
#' @note
#'
#' \strong{Mixed Data}:
#'
#' The mixed data approach is somewhat ad-hoc \insertCite{@see for example p. 277 in  @hoff2007extending;textual}{BGGM}. This
#' is becaue uncertainty in the ranks is not incorporated, which means that variance explained is computed from
#' the 'empirical' \emph{CDF}.
#'
#' \strong{Model Selection}:
#'
#' Currently the default to include all nodes in the model when computing R2. This can be changed (i.e., \code{select = TRUE}), which
#' then sets those edges not detected to zero. This is accomplished by subsetting the correlation matrix according to each neighborhood
#' of relations.
#'
#' @references
#' \insertAllCited{}
#'
#' @examples
#' \donttest{
#' #####################
#' #### continuous #####
#' #####################
#' # data
#' Y <- BGGM::ptsd
#'
#' # estimate the model
#' fit <- estimate(Y, iter = 1000)
#'
#' # predictability
#' r2 <- predictability(fit,
#'                      select = TRUE,
#'                      iter = 1000)
#'
#' # print summary
#' r2
#'
#' # plot
#' plot(r2)
#'
#'
#' #####################
#' ####### mixed #######
#' #####################
#'
#' # estimate the model
#' fit <- estimate(Y,
#'                 type = "mixed",
#'                 iter = 1000)
#'
#' # predictability
#' r2 <- predictability(fit,
#'                      select = TRUE,
#'                      iter = 1000)
#'
#' # print summary
#' r2
#'
#' # plot
#' plot(r2)
#' }
#' @export
predictability <- function(object,
                           select = FALSE,
                           cred = 0.95,
                           iter = 1000){


  if(object$type == "binary" |
     object$type == "ordinal"){

    stop("data type currently not available. must be type = 'continuous' or type = 'mixed'")
  }

  if(object$type == "mixed"){

  Y <- rank_helper(object$Y)$z0_start

  } else {

    # scale data
    Y <- as.matrix(scale(object$Y))

  }
  # nodes
  p <- ncol(Y)

  # observations
  n <- nrow(Y)

  # correlations
  cors <- pcor_to_cor(object)$R[,,1:iter]

  # not conditional on selected model
  if(isFALSE(select)){

    r2 <-  lapply(1:p,   function(x)  {

      # computed from selected model
      .Call(
        "_BGGM_predictability_helper",
        Y[, -x],
        y =  Y[, x],
        XX =  cors[-x,-x, ],
        Xy =  cors[x, -x, ],
        n = n,
        iter = iter
      )$r2
    })

  } else {

    # select model
    sel <- select(object, cred = cred)

    # adjacency
    adj <- sel$adj

    # R2
    r2 <- lapply(1:p, function(x)  {

      # non selected return zero
      if(sum(adj[x,])  == 0 ){
        0

        # a neighborhood exists
      } else {

        # neighborhood
        selected <- which(adj[x,] == 1)

        # check length 1 (return correlation)
        if(length(selected) == 1){

          cors[x, selected,]

          # more than one relation: call c++
        } else {

          # computed from selected model
          .Call(
            "_BGGM_predictability_helper",
            Y[, selected],
            y =  Y[, x],
            XX =  cors[selected, selected, ],
            Xy =  cors[x, selected, ],
            n = n,
            iter = iter
          )$r2
        }
      }
    })

  }

  # R2
  scores <- lapply(1:p, function(x) {

    r2_new <- r2[[x]]
    if(length(r2_new) > 0){

      r2_new[r2_new > 0]
    }

    })

  # returned object
  returned_object <- list(scores = scores,
                          type = "post.pred",
                          metric = "bayes_R2",
                          cred = cred,
                          Y = Y)
  class(returned_object) <- c("BGGM",
                              "metric",
                              "R2",
                              "estimate")
  return(returned_object)
}

