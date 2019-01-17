#' Sample from prior distributions
#'
#' @param distribution "wishart" or "matrix-f"
#' @param p
#' @param delta
#' @param samples
#'
#' @return
#' @export
#'
#' @note This function can be used to determined the prior distribution for hypothesis testing, where it will be very important
#' for determing the graphical structure. For estimation the prior is not as important
#'
#' @examples
sample_prior <- function(distribution, p = NULL, delta = NULL, samples = 500){


  if(distribution == "wishart"){
    if(!is.numeric(p)) stop("p should be the dimensions of the network and numeric")
    prior_samples <- replicate(samples,  cov2cor(rWishart(1, p, Sigma = diag(2))[,,1])[1,2] )

  }


  if(distribution == "matrix_f"){
    if(!is.numeric(delta)) stop("delta must be numeric and positive")

    prior_samples <- extraDistr::rnsbeta(samples, shape1 = delta / 2, shape2 = delta / 2, min = -1, max = 1)

  }

  plt <- ggplot(data.frame(prior_samples), aes(x = prior_samples)) +
         geom_histogram(fill = "lightblue",
                        color = "white") + theme_classic()
    returned_object <- list(plt = plt, prior_samples = prior_samples)

returned_object
    }



