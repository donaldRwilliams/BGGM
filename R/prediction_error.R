#' Title
#'
#' @param fit fitted model from \code{bayes_estimate}
#' @param measure mae = mean absolute error; mse = mean squared error, kl = KL-divergence
#' @param newdata out-of-sample prediction
#' @param selected the selected graphical structure
#' @param ci_width credible interval to summarise measure
#' @param samples number of posterior samples used
#'
#' @return summary as data frame that can be used for plotting with \code{predictive plot}
#' @export
#'
#' @examples
#'
#' # fit the model
#' fit <- bayes_estimate(X)
#'
#' # model selection
#' select <- estimate_select(fit, ci_width = .99)
#'
#' # mean squared error
#' pred_error <- prediction_error(fit, measure = "mse", selected = select$adjacency_mat)
#'
#' # summary information
#' pred_error$summary_error
#'
#' # plot
#' plt <- predictive_plot(pred_error, size = 2, color = "blue")


prediction_error <- function(fit, measure, newdata = NULL, selected, ci_width = .95, samples = 500){



  # ensure data is centered
  dat <- scale(fit$dat, scale = F)

  # new data
  if(!is.null(newdata)) {
    dat <- scale(as.matrix(test), scale = F)
  }

  # sigma for KL Divergence
  if(measure == "kl"){
  sigmas <- BGGM:::inverse_2_beta(fit, samples = samples)$sigmas

  }

  # list for storage
  post_samples <- summary_list <- list()

  for(i in 1:fit$p){

    # selected betas for row (outcome)
    row_select <- selected[i, -i]


    # here no edges were selected
    if(sum(row_select) == 0){

      # mae from intercept only model
      if(measure == "mae") {
        summary_list[[i]] <- mean(abs(dat[,i] - mean(dat[,i])))
        post_samples[[i]] <- 0
        }

      # mse from intercept only model
      if(measure == "mse"){
        summary_list[[i]] <- mean((dat[,i] - mean(dat[,i])^2))
        post_samples[[i]] <- 0
         }

      # zero for KL, but note zero will have divergence of apprx 0
      if(measure == "kl"){
        summary_list[[i]] <- 0
        post_samples[[i]] <- 0
      }
    }else{

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

      # compute error measure
      out <- BGGM:::error_helper(ypred = ypred,
                                 y = dat[,i],
                                 ci_width = ci_width,
                                 measure = measure,
                                 sigmas = sigmas[[i]]
                                 )
      # store the sumamaries
      summary_list[[i]] <- t(data.frame(out$summary))

      # store the posterior samples
      post_samples[[i]] <- out$error
    }
  }

  # rbind all summaries
  summary_error <- do.call(rbind.data.frame, summary_list)

  # row names = column names
  row.names(summary_error) <- colnames(dat)

  # names elements of poserior samples list
  names(post_samples) <- colnames(dat)

  # returned object
  returned_object <- list(summary_error = summary_error,
                          post_samples = post_samples,
                          p = fit$p,
                          measure = measure)

  # assign class (used with other functions)
  class(returned_object) <- "prediction_error"

  # return
  returned_object
  }


