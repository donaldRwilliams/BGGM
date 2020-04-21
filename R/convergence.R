#' MCMC Convergence Plots
#'
#' Monitor converge of the MCMC sampler
#'
#' @param x object of class \code{estimate} or \code{explore}
#' @param type \code{acf} or \code{trace} plot
#' @param param edge name(s) (e.g., "1--2" or c("1--2", "1--3"))
#'
#' @return \code{ggplots}
#' @importFrom stats acf
#' @export
#'
#' @examples
#' # plot
#' Y <- bfi[,1:5]
#'
#' # fit model
#' fit <- estimate(Y)
#'
#' # convergence plots
#' convergence(fit, type = "trace")
convergence <- function(x, type = "acf", param = "1--2"){


  if(class(x) == 'estimate'){

    samps <- x$posterior_samples
    p <-  x$p

    mat <- matrix(0, p, p)
    mat[] <- unlist(lapply(1:p, function(x) paste(1:p, x, sep = "--")))
    samps <- samps[, grep("pcor", colnames(samps))]
    colnames(samps) <-  unlist(as.list(mat))


    if(type == "acf"){
      params <- length(param)
      plts <- lapply(1:params, function(x) {

        dat <- with(acf(samps[,param[x]], plot = FALSE),
                    data.frame(lag, acf));

        ggplot(data = dat, mapping = aes(x = lag, y = acf)) +
          geom_hline(aes(yintercept = 0)) +
          geom_segment(mapping = aes(xend = lag, yend = 0)) +
          ggtitle(param[x])
      } )

    } else if (type == "trace"){

      params <- length(param)
      plts <- lapply(1:params, function(x) {

        dat <- as.data.frame( samps[,param[x]])
        dat$iteration <- 1:nrow(dat)
        ggplot(data = dat, mapping = aes(x = iteration, y = dat[,1])) +
          geom_line(alpha = 0.75) +
          geom_hline(yintercept = mean(dat[,1]), color = "red")+
          ggtitle(param[x]) +
          ylab("Estimate")
      } )


    } else {

      stop("type not supported")
    }

  } else {
    stop("class not currently supported")
  }

  plts
}

