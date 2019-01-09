#' Title
#'
#' @param fit
#' @param newdata
#' @param selected
#' @param ci_width
#' @param samples
#'
#' @return
#' @export
#'
#' @examples
bayes_R2 <- function(fit, newdata = NULL, selected, ci_width = .95, samples = 500){

  # ensure data is centered
  dat <- scale(fit$dat, scale = F)

  # new data
  if(!is.null(newdata)) {
    dat <- scale(as.matrix(test), scale = F)
  }


  # lists for storate
  summary_r2 <- post_samples <- list()

  # compute regression coefficients
  betas <- BGGM:::inverse_2_beta(fit, samples = samples)$betas

  # predicted values for each regression model
  for(i in 1:fit$p){
    dat_temp <- dat[,-i]
    row_select <- selected[i, -i]
    if(sum(row_select) == 0){
      summary_r2[[i]] <- 0
      post_samples[[i]] <- 0
    } else{

      beta_select <- as.matrix(t(apply(betas[[i]], 1, function(x) x * row_select)))

      col_names <- name_helper(colnames(betas[[i]])[row_select == 1])

      beta_select <- beta_select[,which(colMeans(as.matrix(beta_select)) != 0)]
      dat_select <- dat[, as.numeric(col_names)]


      ypred <- t(apply(as.matrix(beta_select), 1, function(x)  x %*% t(as.matrix(dat_select))))
      r2 <- R2_helper(ypred, dat[,i],ci_width)
      summary_r2[[i]] <- t(data.frame(r2$summary_r2))
      post_samples[[i]] <- r2$R2
    }
  }
  summary_r2 <- do.call(rbind.data.frame, summary_r2)
  row.names(summary_r2) <- colnames(dat)
  names(post_samples) <- colnames(dat)
  list(summary_r2 = summary_r2, post_samples = post_samples)
}
