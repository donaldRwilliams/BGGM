#' VAR: Estimation
#'
#' @param Y Matrix (or data frame) of dimensions \emph{n} (observations) by  \emph{p} (variables).
#'
#' @param rho_sd Numeric. Scale of the prior distribution for the partial correlations,
#' approximately the standard deviation of a beta distribution
#' (defaults to 0.50).
#'
#' @param beta_sd Numeric. Standard deviation of the prior distribution for the regression coefficients
#'        (defaults to 1).
#'
#' @param iter Number of iterations (posterior samples; defaults to 5000).
#'
#' @return An object of class \code{var_estimate}
#' @export
var_estimate <- function(Y, rho_sd = 0.50,
                         beta_sd = 1,
                         iter = 5000) {


  Y <- scale(na.omit(Y))
  # number of nodes
  p <- ncol(Y)

  # number of obs
  n <- nrow(Y)

  # Y lagged: add NA
  Y_lag <-   rbind(NA, Y)

  # Y lagged: column names
  colnames(Y_lag) <- paste0(colnames(Y), ".l1")

  # combine all data
  Y_all <- na.omit(cbind.data.frame(rbind(Y, NA), Y_lag))

  # nodes in GGM
  Y <- as.matrix(Y_all[,1:p])

  # predictors (lagged effects)
  X <- as.matrix(Y_all[,(p+1):(p*2)])

  # delta: rho ~ beta(delta/2, delta/2)
  delta <- BGGM:::delta_solve(rho_sd)

  # prior variance
  beta_var <- beta_sd^2

  message(paste0("BGGM: Posterior Sampling "))

  fit <-.Call(
      "_BGGM_var",
      Y =  as.matrix(Y),
      X = as.matrix(X),
      delta = delta,
      epsilon = 0.001,
      beta_prior = diag(p) * (1 / beta_var),
      iter = iter + 50,
      start = solve(cor(Y)),
      progress = TRUE
    )
  message("BGGM: Finished")

  pcor_mu <- round(
    apply(fit$pcors[,,51:(iter + 50)], 1:2, mean),
    digits = 3)

  beta_mu <- round(
    apply(fit$beta[,,51:(iter + 50)], 1:2, mean),
    digits = 3)

  colnames(pcor_mu) <- colnames(Y)
  rownames(pcor_mu) <- colnames(Y)
  colnames(beta_mu) <- colnames(Y)
  row.names(beta_mu) <- colnames(X)

  returned_object <- list(fit = fit,
                          iter = iter,
                          beta_mu = beta_mu,
                          pcor_mu = pcor_mu,
                          p = p,
                          n = n,
                          Y = Y,
                          X = X,
                          call = match.call())

  class(returned_object) <- c("BGGM",
                              "var_estimate",
                              "default")
  return(returned_object)
}


print_var_estimate <- function(x, ...){
  cat("BGGM: Bayesian Gaussian Graphical Models \n")
  cat("--- \n")
  cat("Vector Autoregressive Model (VAR) \n")
  cat("--- \n")
  # number of iterations
  cat("Posterior Samples:", x$iter, "\n")
  # number of observations
  cat("Observations (n):", x$n,"\n")
  # number of variables
  cat("Nodes (p):", x$p, "\n")
  # number of edges
  cat("--- \n")
  cat("Call: \n")
  print(x$call)
  cat("--- \n")
  cat("Partial Correlations: \n\n")
  print(x$pcor_mu)
  cat("--- \n")
  cat("Coefficients: \n\n")
  print(x$beta_mu)
  cat("--- \n")
  cat("Date:", date(), "\n")
}


#' @title Summary Method for \code{var_estimate} Objects
#'
#' @name summary.var_estimate
#'
#' @description Summarize the posterior distribution of each partial correlation
#' and regression coefficient with the posterior mean, standard deviation, and
#' credible intervals.
#'
#' @param object An object of class \code{var_estimate}
#'
#' @param cred Numeric. The credible interval width for summarizing the posterior
#' distributions (defaults to 0.95; must be between 0 and 1).
#'
#' @param ... Currently ignored.
#'
#' @seealso \code{\link{var_estimate}}
#'
#' @return A dataframe containing the summarized posterior distributions,
#' including both the partial correlations and the regression coefficients.
#' @export
summary.var_estimate <- function(object,
                                 cred = 0.95,
                                 ...){

  # nodes
  p <- object$p

  # identity matrix
  I_p <- diag(p)

  # lower bound
  lb <- (1 - cred) / 2

  # upper bound
  ub <- 1 - lb

  # column names
  cn <-  colnames(object$Y)

  if(is.null(cn)){

    mat_names <- sapply(1:p , function(x) paste(1:p, x, sep = "--"))[upper.tri(I_p)]

  } else {

    mat_names <-  sapply(cn , function(x) paste(cn, x, sep = "--"))[upper.tri(I_p)]

  }

  pcor_mean <- round(
    apply(object$fit$pcors[, , 51:(object$iter + 50)], 1:2, mean),
    3)[upper.tri(I_p)]

  pcor_sd  <- round(
    apply(object$fit$pcors[,, 51:(object$iter + 50) ], 1:2, sd),
    digits = 3)[upper.tri(I_p)]

  pcor_lb <- round(
    apply( object$fit$pcors[,, 51:(object$iter + 50) ], 1:2, quantile, lb),
    digits = 3)[upper.tri(I_p)]

  pcor_ub <- round(
    apply(object$fit$pcors[,, 51:(object$iter + 50) ], 1:2, quantile, ub),
    digits =  3)[upper.tri(I_p)]

  beta_mean <- round(
    apply(object$fit$beta[,, 51:(object$iter + 50) ], 1:2, mean),
    digits = 3)

  beta_sd  <- round(
    apply(object$fit$beta[,, 51:(object$iter + 50) ], 1:2, sd),
    digits = 3)

  beta_lb <- round(
    apply( object$fit$beta[,, 51:(object$iter + 50) ], 1:2, quantile, lb),
    digits =  3)

  beta_ub <- round(
    apply(object$fit$beta[,, 51:(object$iter + 50) ], 1:2, quantile, ub),
    digits = 3)

  pcor_results <-
    data.frame(
      relation = mat_names,
      post_mean =  pcor_mean,
      post_sd = pcor_sd,
      post_lb = pcor_lb,
      post_ub = pcor_ub
    )

  colnames(pcor_results) <- c(
    "Relation",
    "Post.mean",
    "Post.sd",
    "Cred.lb",
    "Cred.ub")

  beta_results <-
    lapply(1:p, function (x) {
      res_p <- data.frame(
        relation = colnames(object$X),
        post_mean =  beta_mean[, x],
        post_sd = beta_sd[, x],
        post_lb = beta_lb[, x],
        post_ub = beta_ub[, x]
      )

      colnames(res_p) <- c("Relation",
                           "Post.mean",
                           "Post.sd",
                           "Cred.lb",
                           "Cred.ub")
      res_p
    })

  names(beta_results) <- colnames(object$Y)

  returned_object <- list(pcor_results = pcor_results,
                          beta_results = beta_results)

  class(returned_object) <- c("BGGM",
                              "var_estimate",
                              "summary.var_estimate")
  return(returned_object)


}

print_summary_var_estimate <- function(x, param = "all", ...){
  p <- nrow(x$beta_results[[1]])
  cn <- gsub("\\..*","" , x$beta_results[[1]]$Relation)
  cat("BGGM: Bayesian Gaussian Graphical Models \n")
  cat("--- \n")
  cat("Vector Autoregressive Model (VAR) \n")
  cat("--- \n")
  if(param == "all" | param == "pcor"){
    cat("Partial Correlations: \n\n")
    print(x$pcor_results, row.names = FALSE)
    cat("--- \n\n")
  }
  if(param == "all" | param == "beta") {
    cat("Coefficients: \n\n")
    for (i in seq_len(p)) {
      # print outcome
      cat(paste0(cn[i], " \n\n"))
      # coefs for node i
      coef_i <- x$beta_results[[i]]
      # predictor names
      # colnames(coef_i) <- cn[-i]
      # print coefs
      print(coef_i, row.names = FALSE)
      cat("---\n")
    }
  }
}




#' Plot \code{summary.var_estimate} Objects
#'
#' @description Visualize the posterior distributions of each partial correlation and
#' regression coefficient.
#'
#' @param x An object of class \code{summary.var_estimate}
#'
#' @param color Character string. The color for the error bars.
#' (defaults to \code{"black"}).
#'
#' @param size  Numeric. The size for the points (defaults to \code{2}).
#'
#' @param width Numeric. The width of error bar ends (defaults to \code{0}).
#'
#' @param param Character string. Which parameters should be plotted ? The options
#' are \code{pcor}, \code{beta}, or \code{all} (default).
#'
#' @param order Logical. Should the relations be ordered by size (defaults to \code{TRUE}) ?
#'
#' @param ... Currently ignored
#'
#' @return A list of \code{ggplot} objects.
#'
#' @export
plot.summary.var_estimate <- function(x,
                                      color = "black",
                                      size = 2,
                                      width = 0,
                                      param = "all",
                                      order = TRUE,
                                      ...){
  if(param == "all" |
     param == "pcor"){

    dat_temp <- x$pcor_results

    if(isTRUE(order)){
      dat_temp <-  dat_temp[order(dat_temp$Post.mean,
                                  decreasing = FALSE), ]
    }

    dat_temp$Relation <-
      factor(dat_temp$Relation,
             levels = dat_temp$Relation,
             labels = dat_temp$Relation)

    pcor_plt <- ggplot(dat_temp,
                       aes(x = Relation,
                           y = Post.mean)) +
      geom_errorbar(aes(ymax = dat_temp[, 4],
                        ymin = dat_temp[, 5]),
                    width = width,
                    color = color) +
      geom_point(size = size) +
      xlab("Index") +
      theme(axis.text.x = element_text(
        angle = 90,
        vjust = 0.5,
        hjust = 1
      )) +
      ggtitle("Partial Correlations")

    if(param == "pcor"){
      beta_plt <- NULL
    }

  }
  if (param == "all" | param == "beta") {

    cn <-  names(x$beta_results)
    p <- nrow(x$beta_results[[1]])

    beta_plt <- lapply(1:p, function(i) {

      dat_temp <- x$beta_results[[i]]

      if(isTRUE(order)){
        dat_temp <-  dat_temp[order(dat_temp$Post.mean,
                                    decreasing = FALSE),]
      }

      dat_temp$Relation <-
        factor(dat_temp$Relation,
               levels = dat_temp$Relation,
               labels = dat_temp$Relation)

      ggplot(dat_temp,
             aes(x = Relation,
                 y = Post.mean)) +
        geom_errorbar(aes(ymax = dat_temp[, 4],
                          ymin = dat_temp[, 5]),
                      width = width,
                      color = color) +
        geom_point(size = size) +
        xlab("Index") +
        theme(axis.text.x = element_text(
          angle = 90,
          vjust = 0.5,
          hjust = 1
        )) +
        ggtitle(paste0("Coefficients: ", cn[i]))

    })

    names(beta_plt) <- cn

    if (param == "beta") {
      pcor_plt <- NULL
    }
  }

  list(pcor_plt = pcor_plt, beta_plt = beta_plt)
}
