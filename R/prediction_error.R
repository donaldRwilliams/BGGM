#' Title
#'
#' @param fit
#' @param measure
#' @param newdata
#' @param selected
#' @param ci_width
#' @param samples
#'
#' @return
#' @export
#'
#' @examples
prediction_error <- function(fit, measure = c("KL"), newdata = NULL, selected, ci_width = .95, samples = 500){



  # ensure data is centered
  dat <- scale(fit$dat, scale = F)

  # new data
  if(!is.null(newdata)) {
    dat <- scale(as.matrix(test), scale = F)
  }



  betas <- BGGM:::inverse_2_beta(fit, samples = samples)$betas

  if(measure == "kl"){


  sigmas <- BGGM:::inverse_2_beta(fit, samples = samples)$sigmas
  dat <- scale(fit$dat, scale = F)
  }
  dat <- fit$dat

  post_samples <- summary_list <- list()

  for(i in 1:fit$p){

    # selected betas for row (outcome)
    row_select <- selected[i, -i]


    #
    if(sum(row_select) == 0){
      #
      if(measure == "mae") {
        summary_list[[i]] <- mean(abs(dat[,i] - mean(dat[,i])))
        post_samples[[i]] <- 0
        }
        if(measure == "mse"){
          summary_list[[i]] <- mean((dat[,i] - mean(dat[,i])^2))
          post_samples[[i]] <- 0
         }
       if(measure == "kl"){
      summary_list[[i]] <- 0
      post_samples[[i]] <- 0
      }
     }
      else{

      beta_select <- as.matrix(t(apply(betas[[i]], 1, function(x) x * row_select)))

      col_names <- BGGM:::name_helper(colnames(betas[[i]])[row_select == 1])

      beta_select <- beta_select[,which(colMeans(as.matrix(beta_select)) != 0)]
      dat_select <- dat[, as.numeric(col_names)]






      ypred <- t(apply(as.matrix(beta_select), 1, function(x)  x %*% t(as.matrix(dat_select))))




      out <- BGGM:::error_helper(ypred = ypred,
                                 y = dat[,i],
                                 ci_width = ci_width,
                                 measure = measure,
                                  sigmas = sigmas[[i]]
                                 )
      summary_list[[i]] <- t(data.frame(out$summary))
      post_samples[[i]] <- out$error
    }}


  summary_error <- do.call(rbind.data.frame, summary_list)
  row.names(summary_error) <- colnames(dat)
  names(post_samples) <- colnames(dat)

  returned_object <- list(summary_error = summary_error, post_samples = post_samples, p = fit$p, measure = measure)

  class(returned_object) <- "prediction_error"

  returned_object



  }


