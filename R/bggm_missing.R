#' GGM: Missing Data
#'
#' Estimation and exploratory hypothesis testing with missing data.
#'
#' @param x An object of class \code{mid} \code{\link[mice]{mice}}.
#'
#' @param method Character string. Which method should be used (default set to \code{estimate})? The current
#'               options are \code{"estimate"} and \code{"explore"}.
#'
#' @param ... Additional arguments passed to either
#'            \code{\link{estimate}} or \code{\link{explore}}.
#'
#' @return an object of class \code{estimate} or \code{explore}
#' @export
#'
#' @note Currently, \strong{BGGM} is compatible with the package \code{\link[mice]{mice}} for handling
#'       the missing data. This is accomplished by fitting a model for each imputed dataset
#'       (i.e., more than one to account for uncertainty in the imputation step) and then pooling
#'       the estimates.
#'
#'       In a future version, an additional option will be added that allows for
#'       imputing the missing values during model fitting. This option will be incorporated directly into
#'       the \code{\link{estimate}} or \code{\link{explore}} functions, such that \code{bggm_missing} will
#'       always support missing data with \code{\link[mice]{mice}}.
#'
#' \strong{Support}:
#'
#'  There is limited support for missing data. As of version \code{2.0.0}, it is possible to
#'  determine the graphical structure with either  \code{\link{estimate}} or \code{\link{explore}}, in addition
#'  to plotting the graph with \code{\link{plot.select}}. All data types \emph{are} currently supported.

#' @examples
#' \donttest{
#' # note: iter = 250 for demonstrative purposes
#'
#' library(mice)
#'
#' # data
#' Y <- ptsd
#'
#' # matrix for indices
#' mat <- matrix(0, nrow = 221, ncol = 20)
#'
#' # indices
#' indices <- which(mat == 0, arr.ind = TRUE)
#'
#' # 50 NAs
#' Y[indices[sample(1:nrow(indices), 50),]] <- NA
#'
#' # impute
#' x <- mice(Y, m = 5, print = FALSE)
#'
#' #########################
#' #######   copula    #####
#' #########################
#' # rank based parital correlations
#'
#' # estimate the model
#' fit_est <-  bggm_missing(x, type =  "mixed",
#'                          iter = 250)
#'
#' # select edge set
#' E <- select(fit_est)
#'
#' # plot E
#' plt_E <- plot(E)
#'
#'}
bggm_missing <- function(x, method = "estimate", ...){

  # check for mice
  if(!requireNamespace("mice", quietly = TRUE)) {
    stop("Please install the '", "mice", "' package.")
  }
  # check for abind
  if(!requireNamespace("abind", quietly = TRUE)) {
    stop("Please install the '", "abind", "' package.")
  }

  # combine data in long format
  data_sets <- mice::complete(x, action = "long")

  # number of data sets
  n_data_sets <- length(unique(data_sets$.imp))

  # remove row id
  Y <- data_sets[,-c(2)]

  if(method == "explore"){

    # fit the models
    fits <- lapply(1:n_data_sets, function(x) explore(subset(Y, .imp == x)[,-1], ...))

    # iterations
    iter <- fits[[1]]$iter

    # partial correlations
    post_start_pcors <-  fits[[1]]$post_samp$pcors[,,]

    # fisher z
    post_start_fisher <- fits[[1]]$post_samp$fisher_z

    # regression (for multivariate)
    if(!is.null( fits[[1]]$formula)){

      post_start_beta <- fits[[1]]$post_samp$beta
    }

    # combinate the imputations
    samps <- for(i in 2:n_data_sets) {

      post_start_pcors <-  abind::abind(post_start_pcors ,
                                        fits[[i]]$post_samp$pcors[,,])

      post_start_fisher <-  abind::abind(post_start_fisher,
                                         fits[[i]]$post_samp$fisher_z[,,])

      # multivarate
     if(!is.null(fits[[1]]$formula)){

       post_start_beta <-  abind::abind(post_start_beta,
                                        fits[[i]]$post_samp$beta[,,])
       }
    }

    # dimensions
   dims <- dim(post_start_pcors)

   # replace samples
   fits[[1]]$post_samp$pcors <- post_start_pcors[,,sample(1:dims[3], iter + 50)]
   fits[[1]]$post_samp$fisher_z <- post_start_fisher[,,sample(1:dims[3], iter + 50)]

   if(!is.null( fits[[1]]$formula)){

      fits[[1]]$post_samp$beta <- post_start_beta
     }

  }

  # estimate models
  if(method == "estimate"){

    # fit the models
    fits <- lapply(1:n_data_sets, function(x) estimate(subset(Y, .imp == x)[,-1], ...))
    iter <- fits[[1]]$iter
    post_start_pcors <-  fits[[1]]$post_samp$pcors[,,]
    post_start_fisher <- fits[[1]]$post_samp$fisher_z

    if(!is.null( fits[[1]]$formula)){

      post_start_beta <- fits[[1]]$post_samp$beta

    }

    samps <- for(i in 2:n_data_sets) {

      post_start_pcors <-  abind::abind(post_start_pcors,
                                        fits[[i]]$post_samp$pcors[,,])

      post_start_fisher <-  abind::abind(post_start_fisher,
                                         fits[[i]]$post_samp$fisher_z[,,])

      if(!is.null( fits[[1]]$formula)){
        post_start_beta <-  abind::abind(post_start_beta,
                                         fits[[i]]$post_samp$beta[,,])
      }
    }

    dims <- dim(post_start_pcors)

    fits[[1]]$post_samp$pcors <- post_start_pcors[,,sample(1:dims[3], iter + 50)]
    fits[[1]]$post_samp$fisher_z <- post_start_fisher

    if(!is.null( fits[[1]]$formula)){
      fits[[1]]$post_samp$beta <- post_start_beta
    }
  }

  fit <- fits[[1]]
  fit
}
