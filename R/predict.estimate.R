#' Nodewise In- and Out-of-Sample Predictive Accuracy
#'
#' @description Bayesian predictive error. The measures are computed with respect to the posterior distributions.
#' This provides uncertainty for variance explained (Bayesian R2) and mean squared error.
#' @param fit fitted object of class estimate
#' @param test_data option test data set
#' @param ci_width width of the measure
#' @param samples number of samples used to
#' @param measure Bayesian R2 (R2) or mean squared error (MSE)
#'
#' @return summary_error posterior mean, standard deviation, and credible interval for each variable
#' @return  posterior_samples list containing the computed measure for each posterior draw
#' @export
#'
#' @references Gelman, A., Goodrich, B., Gabry, J., & Vehtari, A. (2017). R-squared for Bayesian regression models.
#'             The American Statistician, (just-accepted), 1-6. \href{http://www.stat.columbia.edu/~gelman/research/unpublished/bayes_R2.pdf}{pre-print}
#'
#' @examples
#'
#' dat <- BGGM::ptsd
#' fit <- estimate(dat, iter = 5000)
#'
#' # select graph
#' selected <- select(fit, ci_width = 0.95)$adjacency_mat
#'
#' error <- predict(fit,
#'               measure = "MSE",
#'               test_data  = NULL,
#'               ci_width = 0.95,
#'               samples = 1000)
#'
#' summary(error)
#'
#' plot(error)
#'
#' ################################
#' ## training/test data example ##
#' ################################
#'
#' train_dat <- BGGM::ptsd[1:200,]
#' test_dat <- BGGM::ptsd[201:221,]
#' fit_train <- estimate(train_dat, iter = 5000)
#'
#' selected <- select(fit_train, ci_width = 0.95)$adjacency_mat
#'
#' r2 <- predict(fit = fit_train,
#'               measure = "R2",
#'               test_data  = test_dat,
#'               ci_width = 0.95,
#'               samples = 1000)
#'
#' summary(r2)
#'
#' plot(r2)
#'
#' # plot test and training error
#' plot(x2 = r2, x1 = error)
#'

predict.estimate <- function(fit, test_data = NULL, ci_width, samples = 1000, measure = c("R2", "MSE")){

    selected <- select(fit, ci_width = ci_width)$adjacency

    # check class
    if(class(fit) != "estimate"){
      stop("must be of class estimate")
    }

    # check samples
    if(samples > length(fit$posterior_samples[,1])){
      stop("samples must be equal to or smaller than the posterior samples")
    }


   # ensure data is scaled
    dat <- scale(fit$dat, scale = T)

    # new data
    if(!is.null(test_data)) {
      if(ncol(dat) != ncol(test_data)) {
        stop("the dimensions of the training and test data must be the same")
      }
      # ensure data is scaled
      dat <- na.omit(scale(as.matrix(test_data), scale = T))
    }

    # lists for storate
    summary <- post_samples <- list()

    # compute regression coefficients
    betas <- BGGM:::inverse_2_beta(fit, samples = samples)$betas

    # predicted values for each regression model
    for(i in 1:fit$p){

      # selected betas for row (outcome)
      row_select <- selected[i, -i]

      # here no edges were selected
      if(sum(row_select) == 0){
        summary[[i]] <- 0
        post_samples[[i]] <- 0
      } else{

        # selected betas (i.e., mulitpled by 0 or 1)
        beta_select <- as.matrix(t(apply(betas[[i]], 1, function(x) x * row_select)))

        # column number as selected names
        col_names <- BGGM:::name_helper(colnames(betas[[i]])[row_select == 1])

        # select the betas
        beta_select <- beta_select[,which(colMeans(as.matrix(beta_select)) != 0)]

        # selected predictors
        dat_select <- dat[, as.numeric(col_names)]

        # predictions
        ypred <- t(apply(as.matrix(beta_select), 1, function(x)  x %*% t(as.matrix(dat_select))))

        if(measure == "R2"){
        # compute error measure
        r2 <- BGGM:::R2_helper(ypred = ypred, y = dat[,i], ci_width = 0.95)
        # store the sumamaries
        summary[[i]] <- t(data.frame(r2$summary_r2))
        # store the posterior samples
        post_samples[[i]] <- r2$R2
        }
        if(measure == "MSE"){
          mse <- BGGM:::MSE_helper(ypred = ypred, y = dat[,i], ci_width = 0.95)
          summary[[i]] <- t(data.frame(mse$summary_mse))
          post_samples[[i]] <- mse$MSE
        }



        }
    }

    # rbind all summaries
    summary <- do.call(rbind.data.frame, summary)

    # row names = column names
    row.names(summary) <- colnames(dat)

    # names elements of poserior samples list
    names(post_samples) <- colnames(dat)

    # returned object
    returned_object <- list(summary_error = summary,
                            post_samples = post_samples,
                            p = fit$p,
                            measure = measure,
                            test_data = test_data)
    returned_object$call <- match.call()

    # assign class (used with other functions)
    class(returned_object) <- "predict"








    # returned object
    returned_object

  }




#' Title
#'
#' @param x
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
summary.predict <- function(x,  ...){
  cat("BGGM: Bayesian Gaussian Graphical Models \n")
  cat("--- \n")
  if(is.null(x$test_data)){
    cat("Type: In-sample predictive accuracy \n")
  } else{
    cat("Type: Out-of-sample predictive accuracy \n")

  }
  if(x$measure == "R2"){
    measure <-  "Variance Explained (R2) \n"
  } else{
    measure <-  "Mean Squared Error (MSE) \n"

  }
  cat("Measure:", measure)
  cat("--- \n")
  cat("Call:\n")
  print(x$call)
  cat("--- \n")
  cat("Estimates: \n\n")
  temp <- cbind.data.frame(node = 1:nrow(x$summary_error), x$summary_error)
  rownames(temp) <- c()
  print(temp,  row.names = FALSE, ...)
  cat("--- \n")
}


#' Title
#'
#' @param x
#' @param nrow
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
head.predict <- function(x, nrow,  ...){
  cat("BGGM: Bayesian Gaussian Graphical Models \n")
  cat("--- \n")
  if(is.null(x$test_data)){
    cat("Type: In-sample predictive accuracy \n")
  } else{
    cat("Type: Out-of-sample predictive accuracy \n")

  }
  if(x$measure == "R2"){
    measure <-  "Variance Explained (R2) \n"
  } else{
    measure <-  "Mean Squared Error (MSE) \n"

  }
  cat("Measure:", measure)
  cat("--- \n")
  cat("Call:\n")
  print(x$call)
  cat("--- \n")
  cat("Estimates: \n\n")
  temp <- cbind.data.frame(node = 1:nrow(x$summary_error), x$summary_error)
  rownames(temp) <- c()
  print(temp[1:nrow,],  row.names = FALSE, ...)
  cat("--- \n")
}

#' Title
#'
#' @param x
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
print.predict <- function(x,  ...){
  cat("BGGM: Bayesian Gaussian Graphical Models \n")
  cat("--- \n")
  if(is.null(x$test_data)){
    cat("Type: In-sample predictive accuracy \n")
  } else{
    cat("Type: Out-of-sample predictive accuracy \n")

  }
  if(x$measure == "R2"){
    measure <-  "Variance Explained (R2) \n"
  } else{
    measure <-  "Mean Squared Error (MSE) \n"

  }
  cat("Measure:", measure)
  cat("--- \n")
  cat("Call:\n")
  print(x$call)
}
