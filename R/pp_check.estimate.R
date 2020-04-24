#' Posterior Predictive Checks for \code{BGGM} Objects
#' @name pp_check.estimate
#' @aliases pp_check
#' @description Perform posterior predictive checks with the help
#' of the \pkg{bayesplot} package (with code taken from \pkg{brms}).
#'
#' @param object object of of class \code{estimate}
#' @param iter number of posterior samples used
#' @param type type of ppc plot. (\code{type = "xyz"}
#' provides a list of suppored types)
#' @param ... currently ignored
#'
#' @details further details are provided here \code{\link[bayesplot:PPC-overview]{PPC}}
#' @examples
#' # data
#' Y <- bfi[,1:5]
#'
#' # fit model
#' fit <- estimate(Y)
#'
#' pp_check(fit, iter = 50,
#' type = "stat")
#' @importFrom bayesplot pp_check
#' @export pp_check
#' @export
pp_check.estimate <- function(object, iter,
                              type = "stat",
                              ...){


  dots <- list(stat = "sd")
  p <- object$p
  dat <- object$dat

  if(missing(iter)){

    iter <- object$iter

    }

  # taken from brms
  valid_types <- as.character(bayesplot::available_ppc(""))
  valid_types <- sub("^ppc_", "", valid_types)

if (!type %in% valid_types) {
  stop("Type '", type, "' is not a valid ppc type. ",
        "Valid types are:\n",
       paste0(valid_types[-grep("grouped", x = valid_types)],
              collapse = ", "))
}

ppc_fun <- get(paste0("ppc_", type),
               asNamespace("bayesplot"))


plts <- list()

pred <- posterior_predict(object,
                          iter = iter,
                          summary = F)$pred

for(i in 1:p){

  ppc_args <- c(list(y = dat[,i],
                   yrep = pred[[i]][1:iter,]), dots)

  plt <- do.call(ppc_fun, ppc_args)

  plts[[i]] <- plt + ggtitle(paste("Node", i))

  }

plts

}

