#' Bayesian Variance Explained (R2)
#'
#' @name bayes_R2
#' @description  Compute Bayesian R2. In contrast to the functions
#' \code{\link{mse}}, \code{\link{mae}}, etc., this can be used to
#' compare predictabiltiy between nodes within a network or between
#' networks. Also, only posterior predictive R2 is implemented.
#'
#' @param object object of class \code{estimate}
#' @param cred credible interval width used for selecting the network
#' @param iter iterations used for computing R2
#' @param cores number of cores for parallel computing
#' @param ... currently ignored
#'
#' @return object of classes \code{bayes_R2} and \code{metric}
#' @export
#'
#' @examples
#' \donttest{
#' # data
#' Y <- subset(tas, gender == "M")[,-ncol(tas)]
#'
#' # fit model
#' fit <- estimate(Y)
#'
#' # bayes R2
#' r2 <- bayes_R2(fit, iter = 50)
#'
#' # print summary
#' r2
#'
#' # plot
#' plot(r2)
#' }
bayes_R2 <- function(object,
                           cred = 0.95,
                           iter = 1000,
                           cores = 2,...){

  p <- object$p

  cl <- parallel::makeCluster(cores)

  doParallel::registerDoParallel(cl)

  adj <-  select(object, cred = cred)$adjacency_non_zero

  betas <- inverse_2_beta(object, samples = iter)

 scores <- parallel::parLapply(cl = cl,
                              X = 1:p, function(x) R2_ppc(fit = object,
                                                          betas = betas,
                                                          adj = adj,
                                                          which_one = x,
                                                          sims = iter))

 parallel::stopCluster(cl)

 # returned object
 returned_object <- list(scores = scores,
                         type = "post.pred",
                         metric = "bayes_R2",
                         cred = cred)

 class(returned_object) <- c("metric", "R2")

 return(returned_object)

}

#' Assess Predictability
#'
#' @name test.R2
#'
#' @description  Compare nodes within networks or between networks. Currently the only
#' option available is Bayesian R2.
#'
#' @param ... object(s) of class \code{R2}
#'
#' @return object of class \code{metric}
#' @export
#'
#' @examples
#' \donttest{
#' # two groups (males vs females)
#' Y_m <- subset(tas, gender == "M")[,-ncol(tas)]
#' Y_f <- subset(tas, gender == "F")[,-ncol(tas)]
#'
#' # fit models
#' fit_m <- estimate(Y_m)
#' fit_f <- estimate(Y_f)
#'
#' # r2
#' r2_m <- bayes_R2(fit_m, iter = 50)
#' r2_f <- bayes_R2(fit_f, iter = 50)
#'
#' # assess predictability
#' assess_pred  <- assess_predictability(r2_m, r2_f)
#'
#' # summary
#' assess_pred
#'
#' # plot (small rope value removes it from the plot )
#' plot(assess_pred, rope = 0.001)
#' }
assess_predictability <-  function(...){

  temp <- list(...)

  if(!all(sapply(temp, class)[2,] == "R2")){
    stop("all object must be of class R2 (bayes_R2)")
}

  groups <- length(temp)

  if(groups != 2){
    stop("currently only comparing groups nodewise is possible")
  }

  g1 <- temp[[1]]
  g2 <- temp[[2]]

  p <- length(g1$scores)

  diff <- lapply(1:p, function(x) g1$scores[[x]] - g2$scores[[x]])
  returned_object <- list(scores = diff,
                          type = "post.pred",
                          metric = "bayes_R2_diff")
  class(returned_object) <- "metric"
  returned_object
}
