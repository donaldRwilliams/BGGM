#' Compute Custom Network Statistics
#'
#' This function allows for computing custom network statistics for
#' weighted adjacency matrices (partial correlations). The statistics are computed for
#' each of the sampled matrices, resulting in a distribution.
#'
#' @param object An object of class \code{estimate}.
#'
#' @param FUN A custom function for computing the statistic. The first argument must be
#'            a partial correlation matrix.
#'
#' @param iter  Number of iterations (posterior samples; defaults to the number in the object).
#'
#' @param select Logical. Should the graph be selected ? The default is currently \code{FALSE}.
#'
#' @param cred Numeric. Credible interval between 0 and 1  (default is 0.95) that is used for selecting the graph.
#'
#' @param ... Arguments passed to the function.
#'
#' @return An object defined by \code{FUN}.
#'
#' @details
#'
#' The user has complete control of this function. Hence, care must be taken as to what \code{FUN} returns
#' and in what format.
#'
#'
#' @examples
#' \donttest{
#'
#' ####################################
#' ###### example 1: assortment #######
#' ####################################
#' library(assortnet)
#'
#' Y <- BGGM::bfi[,1:10]
#' membership <- c(rep("a", 5), rep("c", 5))
#'
#' # fit model
#' fit <- estimate(Y = Y,
#'                 analytic = FALSE,
#'                 iter = 1000)
#'
#' # list of columns belowinging in each group
#' e.g., first 5 are "a", last 5 are "c"
#'
#' membership <- c(rep("a", 5), rep("c", 5))
#'
#'f <- function(x,...){
#' assortment.discrete(x, ...)$r
#'}
#'
#'
#' net_stat <- roll_your_own(object = fit,
#'                           FUN = f,
#'                           types = membership,
#'                           weighted = TRUE,
#'                           SE = FALSE, M = 1)
#'
#' hist(net_stat)
#'
#'
#' ############################################
#' ###### example 2: expected influence #######
#' ############################################
#' library(networkTools)
#'
#' # data
#' Y <- depression
#'
#' # fit model
#' fit <- estimate(Y = Y, iter = 5000)
#'
#' # define function
#' f <- function(x,...){
#'   expectedInf(x,...)$step1
#' }
#'
#' # compute
#' net_stat <- roll_your_own(object = fit,
#'                           FUN = f,
#'                           iter = 1000)
#'
#' colmeans
#' colMeans(t(net_stat))
#'
#' # full distribution (node 1)
#' hist(net_stat[1,])
#' }
#' @export
roll_your_own <- function(object,
                          FUN,
                          iter = NULL,
                          select = FALSE,
                          cred = 0.95, ...){

  if(! all( c("estimate", "default") %in% class(fit)) ){
    stop("class must be 'estimate'")
  }

  if(!isFALSE(select)){

    sel <- select(fit, cred = cred)
    adj <- sel$adj

  } else {

    p <- ncol(object$pcor_mat)
    adj <- matrix(1, p, p)

  }

  if(is.null(iter)){

    iter <- object$iter

  }

  pcors <- object$post_samp$pcors[, , 51:(iter + 50)]

  pb <- utils::txtProgressBar(min = 0, max = iter, style = 3)

  results <- sapply(1:iter, function(x) {

    pcors_s <- pcors[, , x] * adj

    est <- FUN(pcors_s, ...)

    utils::setTxtProgressBar(pb, x)

    est
  })

  returned_object <- list(results = results, iter = iter)

  class(returned_object) <- c("BGGM",
                              "roll_your_own")

  return(returned_object)
  }


print_roll_your_own <- function(x, cred = 0.95, ...) {
  cat("BGGM: Bayesian Gaussian Graphical Models \n")
  cat("--- \n")
  cat("Network Stats: Roll Your Own\n")
  cat("Posterior Samples:", x$iter, "\n")
  cat("--- \n")
  cat("Estimates: \n\n")
  lb <- (1-cred) / 2
  ub <- 1 - lb
  dims <- dim(x$results)
  if(is.null(dims)){
    mu <- mean(x$results)
    scale <- sd(x$results)
    res <- data.frame(Post.mean = round(mean(x$results), 3),
                      Post.sd =    round(sd(x$results), 3),
                      Cred.lb = round(quantile(x$results, probs = lb), 3),
                      Cred.ub = round(quantile(x$results, probs = lb), 3) )
  } else {

    mu <-  apply( x$results, 1, mean)
    p <- length(mu)
    scale <- apply( x$results, 1, sd)
    ci_lb <- apply( x$results, 1, quantile, lb)
    ci_ub <- apply( x$results, 1, quantile, ub)

    res<- data.frame(Node = 1:p, Post.mean = round(mu, 3),
                     Post.sd = round(scale, 3),
                     Cred.lb = round(ci_lb, 3),
                     Cred.ub = round(ci_ub, 3))
    }

  print(res, row.names = FALSE)
  cat("--- \n")
}


#' Plot \code{roll_your_own} Objects
#'
#' @param x An object of class \code{roll_your_own}
#'
#' @param fill Character string specifying the color for the ridges.
#'
#' @param alpha Numeric. Transparancey of the ridges
#'
#' @param ... Currently ignored
#'
#' @return An object of class \code{ggplot}
#'
#' @importFrom ggridges stat_density_ridges
#'
#' @examples
#' \donttest{
#' ####################################
#' ###### example 1: assortment #######
#' ####################################
#' library(assortnet)
#'
#' Y <- BGGM::bfi[,1:10]
#' membership <- c(rep("a", 5), rep("c", 5))
#'
#' # fit model
#' fit <- estimate(Y = Y,
#'                 analytic = FALSE,
#'                 iter = 1000)
#'
#' # list of columns belowinging in each group
#' e.g., first 5 are "a", last 5 are "c"
#'
#' membership <- c(rep("a", 5), rep("c", 5))
#'
#'f <- function(x,...){
#' assortment.discrete(x, ...)$r
#'}
#'
#'
#' net_stat <- roll_your_own(object = fit,
#'                           FUN = f,
#'                           types = membership,
#'                           weighted = TRUE,
#'                           SE = FALSE, M = 1)
#'
#' # print
#' net_stat
#'
#' # plot
#' plot(net_stat)
#'
#'
#' ############################################
#' ###### example 2: expected influence #######
#' ############################################
#' library(networkTools)
#'
#' # data
#' Y <- depression
#'
#' # fit model
#' fit <- estimate(Y = Y, iter = 5000)
#'
#' # define function
#' f <- function(x,...){
#'   expectedInf(x,...)$step1
#' }
#'
#' # compute
#' net_stat <- roll_your_own(object = fit,
#'                           FUN = f,
#'                           iter = 1000)
#'
#' # print
#' net_stat
#'
#' #plot
#' plot(net_stat)
#' }
#'
#' @export
plot.roll_your_own <- function(x, fill = "#CC79A7", alpha = 0.5, ...){

  dims <- dim(x$results)

  if(is.null(dims)){

    dat <- data.frame(x= x$results, y = 1)

    plt <- ggplot(dat, aes(x = x,
                           y = as.factor(y))) +
      geom_density_ridges(fill = fill,
                          alpha = alpha)

  } else {

    dat <- reshape::melt(t(x$results))

    mus <- tapply(dat$value, dat$X2, mean)

    dat$order <- factor(dat$X2, levels =  unique(dat$X2)[order(mus)],
                        labels = unique(dat$X2)[order(mus)] )

    plt <- ggplot(dat, aes(x = value,
                           group = order,
                           y = as.factor(order))) +
      geom_density_ridges(fill = fill,
                          alpha = alpha) +
      ylab("")


  }
  plt
}
