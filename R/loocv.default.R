#' Nodewise Leave-One-Out Cross-Validation
#'
#' @description  Assesses the predictability of each node in the "network" of GGM. Currently there are two options avaliable. The first is Bayesian
#' leave-one-out cross-validation and has a measure of uncertainty, whereas the latter is based only on the point estimates and is known as PRESS
#' (predicted residual sums of squares).
#'
#' @param x fitted model
#' @param ci_width credible interval width for selected edges
#' @param samples number of samples (\code{analytic = FALSE})
#' @param n_cl number of clusters (see notes)
#'
#' @return list of class \code{loocv}:
#'
#' \itemize{
#'
#' \item \code{returned_object} data.frame summary of predictability
#' \item \code{call} \code{match.call()}
#' \item \code{samples} number of samples used for computing Bayesian loo
#' }
#'
#'
#'
#' @note This function can be used to compute leave-one-out prediction error analytically for each node in the selected graph. This uses the point estimates
#' for the coefficients, and thus does not provide a measure of uncertainty. Becuase this solution is fast it allows for computing, say, prediction error
#' while adjusting the \emph{p} to \emph{n} ratio. This can provide insight into overfitting. See here for the derivations:
#'  \href{https://robjhyndman.com/hyndsight/loocv-linear-models/}{analytic CV}.
#'
#'  Ideally, Bayesian (approximate) leave-one-out (loo) cross-validation should be computed. This requires posterior samples, which provides a measure
#'  of uncertainty--i.e., a standard error. This is computationally more involved
#'  than the analytic expression, but still much (much) faster than prediciting individual data points. \code{n_cl} allows for parallel computation.
#'  \code{samples} number of posterior samples used for computing loo (cannot exceed the number of samples in the fitted model).
#'
#'  This approach should \strong{not} be used to explicitly compare nodes, as each is fit to a different outcome. However, because the variables are on
#'  the same scale (standardized in advance), it is seems that one could infer which node in the graph is the easiest to predict.
#'
#'  methods(class = "loocv")
#'
#' @export
#'
#' @references
#' Allen, D. M. (1971). Mean square error of prediction as a criterion for selecting variables. Technometrics, 13(3), 469-475.
#'
#' Gelman, A., Hwang, J., & Vehtari, A. (2014). Understanding predictive information criteria for Bayesian models.
#' Statistics and computing, 24(6), 997-1016.
#'
#' Vehtari, A., Gelman, A., & Gabry, J. (2017). Practical Bayesian model evaluation using leave-one-out cross-validation
#' and WAIC. Statistics and Computing, 27(5), 1413-1432.
#'
#'
#'
#'
#' @examples
#'
#' X <- BGGM::bfi[,1:5]
#'
#'###########################
#'##### Bayesian loocv ######
#'###########################
#' fit <- estimate(X,
#'                 analytic = FALSE,
#'                samples = 500)
#'
#' nodewise_loo <- loocv(fit,
#'                       ci_width = 0.99,
#'                       samples = 250,
#'                       n_cl = detectCores() - 1)
#'# summary
#' summary(nodewise_loo)
#'
#'# plot
#' plot(nodewise_loo,
#'     size = 4,
#'      color = "brown4") +
#'   theme_bw() +
#'   theme(panel.grid.minor = element_blank(),
#'         panel.grid.major.x = element_blank()) +
#'   ggtitle("Bayesian Leave-One-Out Cross-Validation") +
#'   ylab("Leave-One-Out Prediction Error \n (lower is better)") +
#'   xlab("Node")
#'
#'
#'###########################
#'##### analytic loocv ######
#'###########################
#' fit <- estimate(X,
#'                 analytic = TRUE)
#'
#' nodewise_loo <- loocv(fit,
#'                       ci_width = 0.99)
#'# summary
#' summary(nodewise_loo)
#'
#'# plot
#' plot(nodewise_loo,
#'      size = 6,
#'      color = "pink") +
#'   theme_bw() +
#'   theme(panel.grid.minor = element_blank(),
#'         panel.grid.major = element_blank()) +
#'   ggtitle("Leave-One-Out Cross-Validation") +
#'   ylab("Leave-One-Out Prediction Error \n (lower is better)") +
#'   xlab("Node") +
#'   scale_y_continuous(limits = c(0, 2500),
#'                      expand = c(0, 0))

loocv.default <- function(x, ci_width = 0.95, samples = 100, n_cl = 2){

  # select graph (error conditional on the selected edges)
  selected <- select(x, ci_width = ci_width)$adjacency

  # remove na'a
  X <- na.omit(x$dat)

  # if analytic solution
  if(!isFALSE(x$analytic)){

  # convert to coefficients
  coefs <- lapply(1:x$p, function(z)     selected[z,-z] * (x$fit$inv_mu[z,-z] / x$fit$inv_mu[z,z]) * -1)

  # store loo score
  loo_score  <- list()

  # store residual sum of squares
  rss <- list()

  for(i in 1:x$p) {

    # subset selected columns
    Xnew <- X[, which(selected[i,] == 1)]

    # compute the hat matrix
    hat <-  diag(Xnew %*%  solve((t(Xnew) %*% Xnew)) %*% t(Xnew))

    # predict y
    yhat <- as.matrix(Xnew) %*% coefs[[i]][which(selected[i,-i] == 1)]

    # residuals
    resid <- X[,i] - yhat

    # so-called PRESS statistic
    press <- (resid / (1 - hat))

    # store loo scores
    loo_score[[i]] <- sum(press^2)

    # store rss
    rss[[i]] <- sum(resid ^ 2)
  }
  # returned object
  returned_object <- data.frame(
                                # nodes 1:p
                                node = 1:x$p,
                                # loo scores
                                loo = unlist(loo_score),
                                # rss
                                rss =  unlist(rss))

  returned_object <- list(returned_object = returned_object,
                          call = match.call())
  }
  # if samples are present
  if(isFALSE(x$analytic)){

    # make cluster
    cl <-  parallel::makeCluster(n_cl)

    # register cluster
    doSNOW::registerDoSNOW(cl)

    # ensure data is scaled
    X <- scale(X, scale = T)

    # list for storage
    loo_list <- list()

    # convert inverse to regression
    post_samps <- BGGM:::inverse_2_beta(x, samples = samples)

    # subset betas
    betas <- post_samps$betas

    # subset sigmas
    sigmas <- post_samps$sigmas

    loos <- foreach(i  = 1:x$p) %dopar% {

      # selected row
      row_select <- selected[i, -i]

      # if none were selected
      if(sum(row_select) == 0){
        loo_list[[i]] <- 0
      }

      # beta selected
      beta_select <- as.matrix(t(apply(betas[[i]], 1, function(z) z * row_select)))

      # colnames of the selected
      col_names <- BGGM:::name_helper(colnames(betas[[i]])[row_select == 1])

      # selected betas
      beta_select <- beta_select[,which(colMeans(as.matrix(beta_select)) != 0)]

      # subset selected columns
      dat_select <- X[, as.numeric(col_names)]

      # predict y
      ypred <- t(apply(as.matrix(beta_select), 1, function(z)  z %*% t(as.matrix(dat_select))))

      # compute log-likelihood matrix
      log_lik <- sapply(1:samples, function(j) {dnorm(X[,i] , ypred[j,] , sigmas[[i]][j] , log=TRUE)})

      # store log-likelihood
      loo_list[[i]] <- suppressWarnings(loo::loo(t(log_lik)))

    }

    parallel::stopCluster(cl)

    # summarize loo
    summary_loo <- do.call(rbind.data.frame, lapply(loos, function(z)
                           data.frame(looic = z$estimates[3], looic_se  = z$estimates[6])))

    # row names
    row.names(summary_loo) <- colnames(x$dat)

    # list names
    names(loos) <- colnames(fit$dat)

    # returned object
    returned_object <- data.frame(
                                  # nodes 1:p
                                  node = 1:x$p,
                                  # loo score
                                  loo = summary_loo$looic,
                                  # loo se
                                  loo_se =  summary_loo$looic_se)

    returned_object <- list(returned_object = returned_object,
                            call = match.call(),
                            samples = samples)
    }
  # assign class
  class(returned_object) <- "loocv"
  return(returned_object)
}

