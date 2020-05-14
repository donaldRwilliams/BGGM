#' MCMC Convergence
#'
#' Monitor convergence of the MCMC algorithms.
#'
#' @param object An object of class \code{estimate} or \code{explore}
#'
#' @param param Character string. Names of parameters for which to monitor MCMC convergence.
#'
#' @param type Character string. Which type of convergence plot ? The current
#'             options are \code{trace} (default) and \code{acf}.
#'
#' @param print_names Logical. Should the parameter names be printed (defaults to \code{FALSE})? This
#'                   can be used to first determine the parameter names to specify in \code{type}.
#'
#' @return A list of \code{ggplot} objects.
#'
#' @note An overview of MCMC diagnostics can be found \href{http://sbfnk.github.io/mfiidd/mcmc_diagnostics.html}{here}.
#'
#' @importFrom stats acf
#'
#' @examples
#'
#' \donttest{
#' # note: iter = 250 for demonstrative purposes
#'
#' # data
#' Y <- ptsd
#'
#' #########################
#' ###### continuous #######
#' #########################
#' fit <- estimate(Y, iter = 250)
#'
#' # print names first
#' convergence(fit, print_names = TRUE)
#'
#' # trace plots
#' convergence(fit, type = "trace",
#'             param = c("B1--B2", "B1--B3"))
#'
#' # acf plots
#' convergence(fit, type = "acf",
#'             param = c("B1--B2", "B1--B3"))
#'
#' #########################
#' ######## mixed ##########
#' #########################
#' # copula
#'
#' fit <- estimate(Y, type = "mixed",
#'                 iter = 250)
#'
#' # print names first
#' convergence(fit, print_names = TRUE)
#'
#' # trace plots
#' convergence(fit, type = "trace",
#'             param = c("B1--B2", "B1--B3"))
#'
#' # acf plots
#' convergence(fit, type = "acf",
#'             param = c("B1--B2", "B1--B3"))
#'
#' #########################
#' ######## ordinal ########
#' #########################
#' fit <- estimate(Y + 1, type = "ordinal",
#'                 iter = 250)
#'
#' # print names first
#' convergence(fit, print_names = TRUE)
#'
#' # trace plots
#' convergence(fit, type = "trace",
#'             param = c("B1--B2", "B1--B3"))
#'
#' # acf plots
#' convergence(fit, type = "acf",
#'             param = c("B1--B2", "B1--B3"))
#' }
#' @export
convergence <- function(object,
                        param = NULL,
                        type = "trace",
                        print_names = FALSE){

  # posterior samples
  samps <- posterior_samples(object)

  # print names ?
  if(!isFALSE(print_names)){

    print(colnames(samps))

  } else {

    # trace plot
    if(type == "trace"){

      # number of params
      params <- length(param)

      plts <- lapply(1:params, function(x){

        dat <- as.data.frame( samps[,param[x]])

        dat$iteration <- 1:nrow(dat)

        ggplot(data = dat,
               mapping = aes(x = iteration,
                             y = dat[,1])) +
          geom_line(alpha = 0.75) +
          geom_hline(yintercept = mean(dat[,1]),
                     color = "red")+
          ggtitle(param[x]) +
          ylab("Estimate")
      })

    } else if(type == "acf"){

      params <- length(param)

      plts <- lapply(1:params, function(x) {

        dat <- with(acf(samps[,param[x]],
                        plot = FALSE),
                    data.frame(lag, acf));

        ggplot(data = dat,
               mapping = aes(x = lag,
                             y = acf)) +
          geom_hline(aes(yintercept = 0)) +
          geom_segment(mapping = aes(xend = lag,
                                     yend = 0)) +
          ggtitle(param[x])

      })

    } else {

      stop("type not supported. must be 'trace' or 'acf'")

    }

    plts
  }
}
