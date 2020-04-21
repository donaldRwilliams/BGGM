#' Robustness Check for \code{ggm_compare_bf} Objects
#'
#' @param ... objects of class \code{ggm_compare_bf} (at least two)
#' @param BF_cut evidentiary threshold (default is 3)
#'
#' @return object of class \code{robust_check}
#' @export
#'
#' @examples
#' # all equal
#' Y1 <- MASS::mvrnorm(500, rep(0, 16), Sigma = diag(16))
#' Y2 <- MASS::mvrnorm(500, rep(0, 16), Sigma = diag(16))
#' Y3 <- MASS::mvrnorm(500, rep(0, 16), Sigma = diag(16))
#'
#' # two priors
#' prior_10 <- ggm_compare_bf(Y1, Y2, Y3, prior_sd = .10, iter = 250, cores = 2)
#' prior_20 <- ggm_compare_bf(Y1, Y2, Y3, prior_sd = .20, iter = 250, cores = 2)
#'
#' # robustness check
#' rc <- robust_check(prior_10, prior_20)
#'
#' # print results
#' rc
#' # plot results
#' plot(rc)
robust_check <- function(..., BF_cut = 3){


  models <- list(...)

  cls <- class(models[[1]])
  if(cls != "ggm_compare_bf"){
    stop("object must be of class ggm_compare_bf")

  }

  n_models <- length(models)
  props <- lapply(1:n_models, function(x) {

    sel <- select(models[[x]], BF_cut = BF_cut);
    prior_sd <- models[[x]]$prior_sd;
    null <- sel$BF_01[upper.tri( sel$BF_01)];
    null <- mean(null > BF_cut);
    alt <- sel$BF_10[upper.tri( sel$BF_10)];
    alt <- mean(alt > BF_cut);
    incon <- 1 - (null + alt)
    data.frame(measure = c("null", "alt", "incon"),
               props = c(null, alt, incon),
               prior_sd = prior_sd)
  })


  results <- do.call(rbind.data.frame, props)

  class(results) <- c("robust_check", "data.frame")

   return(results)
}

#' @title  Print method for \code{robust_check} objects
#'
#' @param x object of class \code{robust_check}
#' @param ... currently ignored
#' @export
print.robust_check <- function(x, ... ){
  cat("BGGM: Bayesian Gaussian Graphical Models \n")
  cat("Type: GGM Compare with Bayesian Hypothesis Testing \n")
  cat("--- \n")
  cat("Robustness Check \n\n")
  sds <- unique( as.data.frame(x)$prior_sd)
  for(i in seq(sds)){

    cat("Prior SD:", sds[i], "\n")
    res <- subset(as.data.frame(x), prior_sd == sds[i])[,-3]
    cat("  Equivalent:",
        paste(round(res$props[1] * 100, 2), "%\n"))
    cat("  Different:",
        paste(round(res$props[2] * 100, 2), "%\n"))
    if(i < length(sds)){
    cat("  Inconclusive:",
        paste(round(res$props[3] * 100, 2), "%\n\n"))
    } else {
      cat("  Inconclusive:",
          paste(round(res$props[3] * 100, 2), "%\n"))

  }
  }
    cat("--- \n")
    cat("note: percentage of the total number of edges")

}

#' @title Plot \code{robust_check} Objects
#'
#' @param x object of class \code{robust_check}
#' @param size geom point size
#' @param ... currently ignored
#'
#' @return \code{ggplot} object
#' @export
#'
#' @examples
#' # all equal
#' Y1 <- MASS::mvrnorm(500, rep(0, 16), Sigma = diag(16))
#' Y2 <- MASS::mvrnorm(500, rep(0, 16), Sigma = diag(16))
#' Y3 <- MASS::mvrnorm(500, rep(0, 16), Sigma = diag(16))
#'
#' # two priors
#' prior_10 <- ggm_compare_bf(Y1, Y2, Y3, prior_sd = .10, iter = 250, cores = 2)
#' prior_20 <- ggm_compare_bf(Y1, Y2, Y3, prior_sd = .20, iter = 250, cores = 2)
#'
#' # robustness check
#' rc <- robust_check(prior_10, prior_20)
#'
#' # print results
#' rc
#' # plot results
#' plot(rc)
plot.robust_check <- function(x, size = 2,...){

  dat <- as.data.frame(x)
  dat$measure <- factor(dat$measure,
                        levels = c("alt", "null", "incon"),
                        labels = c("Different",
                                  "Equivalent",
                                   "Inconclusive"))

  ggplot(dat, aes(x = as.factor(prior_sd),
                  y = props,
                  group = measure)) +
  geom_line(aes(linetype = measure)) +
  geom_point(aes(color = measure), size = size) +
  theme_bw() +
    theme(legend.position = "top",
        legend.title=element_blank()) +
  scale_y_continuous(limits = c(0,1),
                     breaks = seq(0,1, 0.25),
                     labels = paste(seq(0,1, 0.25) * 100, "%")) +
  ylab("Percentage of Edges") +
  xlab("Prior SD")


 }
