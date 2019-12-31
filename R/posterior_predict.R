#' Posterior Predictive Distribution for \code{estimate} Objects
#'
#' @inheritParams predict.estimate
#'
#' @return
#' @export
#'
#' @examples
posterior_predict <- function(object, iter = 500,
                              cred = 0.95, newdata = NULL,
                              summary = TRUE,...){

  # check for object class
  if(class(object) != "estimate"){
    stop("object must be of class estimate")
  }

  # check iteration exist
  if(object$iter < iter){
    stop("iter cannot be larger than the number of samples used for fitting")
  }

  # variables in original fitting
  p <- object$p

  lb <- (1 - cred) / 2
  ub <- 1 - lb

  # no newdata
  if(is.null(newdata)){

    # data
    dat <- object$dat

    # sample size
    n <- nrow(dat)

    } else {

      # newdata
      dat <- newdata

      # sample size
      n <- nrow(dat)

      # check newdata has same number of columns
      if(ncol(dat) != p){
        stop("new data not allowed (number of variables differs)")
      }
    }

  # betas
  betas <- BGGM:::inverse_2_beta(object, samples = iter)

  # summarize ?
  if(isTRUE(summary)){

    ppc_array <- array(0, dim = c(n, 4, p ),
                       dimnames = list(1:n,
                                       c("Post.mean",
                                         "Post.sd",
                                         "Cred.lb",
                                         "Cred.ub"),
                                       paste0("node_", 1:p)))

    for(i in 1:p){

      # selected betas
      beta <- betas$betas[[i]]

      # sigmas
      sigma <- betas$sigmas[[i]]

      # posterior predictive
      ppc <- t(sapply(1:iter, function(s) rnorm(n = n,
                                          mean = dat[,-i] %*% t(beta[s,]),
                                          sd = sigma[s])))

      ppc_array[,1,i] <- apply(ppc, 2, mean)
      ppc_array[,2,i] <- apply(ppc, 2, sd)
      ppc_array[,3:4,i] <- t(apply(ppc, 2, quantile, prob = c(lb, ub)))

      }

    returned_object <- ppc_array

    } else {

      # storage for non-summarized
      pred <- list()

      for(i in 1:p){

        # selected betas
        beta <- betas$betas[[i]]

        # sigmas
        sigma <- betas$sigmas[[i]]

        # posterior predictive
        pred[[i]] <- t(sapply(1:iter, function(s) rnorm(n = n,
                                              mean = dat[,-i] %*% t(beta[s,]),
                                              sd = sigma[s])))
        }

  returned_object <- list(pred = pred, dat = dat)

  }
  class(returned_object) <- "post.pred"
  return(returned_object)

  }



