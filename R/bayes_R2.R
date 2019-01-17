#' Bayesian R2
#'
#' @param fit object of class bayes_estimate
#' @param newdata out-of-sample option
#' @param selected the selected adjacency matrix
#' @param ci_width credible interval width for Bayes R2
#' @param samples Number of samples for compute Bayes R2
#'
#' @return
#' @export
#'
#' @examples
#'X <- BGGM::ptsd
#'
#' # X is the data matrix of dimensions n by p
#' fit <- bayes_estimate(X)
#'
#' # select the graphical structure
#' select <- estimate_select(fit, ci_width = .95)
#'
#' # compute Bayesian R2
#' r2 <-  bayes_R2(fit, selected = select$adjacency_mat, ci_width = .95, samples = 1000)
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

    # selected betas for row (outcome)
    row_select <- selected[i, -i]

    # here no edges were selected
    if(sum(row_select) == 0){
      summary_r2[[i]] <- 0
      post_samples[[i]] <- 0
    } else{

      # selected betas (i.e., mulitpled by 0 or 1)
      beta_select <- as.matrix(t(apply(betas[[i]], 1, function(x) x * row_select)))

      # column number as selected names
      col_names <- name_helper(colnames(betas[[i]])[row_select == 1])

      # select the betas
      beta_select <- beta_select[,which(colMeans(as.matrix(beta_select)) != 0)]

      # selected predictors
      dat_select <- dat[, as.numeric(col_names)]

      # predictions
      ypred <- t(apply(as.matrix(beta_select), 1, function(x)  x %*% t(as.matrix(dat_select))))

      # compute error measure
      r2 <- BGGM:::R2_helper(ypred, dat[,i],ci_width)

      # store the sumamaries
      summary_r2[[i]] <- t(data.frame(r2$summary_r2))

      # store the posterior samples
      post_samples[[i]] <- r2$R2
    }
  }

  # rbind all summaries
  summary_r2 <- do.call(rbind.data.frame, summary_r2)

  # row names = column names
  row.names(summary_r2) <- colnames(dat)

  # names elements of poserior samples list
  names(post_samples) <- colnames(dat)

  # returned object
  returned_object <- list(summary_error = summary_r2,
                          post_samples = post_samples,
                          p = fit$p,
                          measure = "bayes_R2")

  # assign class (used with other functions)
  class(returned_object) <- "bayes_R2"

  # returned object
  returned_object

}
