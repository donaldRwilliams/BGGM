#' Mean Squared Error
#' @description Compute mean squared error for either the observed data or future data. The former
#' is computed by plugging in the observed y (the predicted node). The latter is computed from
#' replicated data sets (posterior predictive y), which results in the posterior predictive
#' mean squared error.  Both provide a measure of uncertainty,  as the error is computed from
#' the posterior samples. However, the posterior predictive approach fully captures uncertainty.
#'
#' @name mse
#' @param object object of class \code{post.pred} or  \code{predict.estimate}
#' @param ... currently ignored
#'
#' @return object of class \code{metric}
#' @export
#'
#' @examples
#' # data
#' Y <- subset(tas, gender == "M")[,-ncol(tas)]
#'
#' # fit model
#' fit <- estimate(Y)
#'
#' # predict (note summary = FALSE)
#' pred <- predict(fit, iter = 50, summary = FALSE)
#'
#' mse(pred)
mse <- function(object, ...){

  # check summary is false
  if(length(object) != 2){
    stop("summary must be set to false")
  }

  # data
  dat <- object$dat

  # number of variables
  p <- ncol(object$dat)

  # predictions
  pred <- object$pred

  if(class(object) == "predict.estimate" |
     class(object) == "fitted.estimate"){

    type <- class(object)

    # scores
    scores <- lapply(1:p, function(x) colMeans((t(pred[[x]]) - dat[,x])^2))

  } else if (class(object) == "post.pred"){

    type <- class(object)

    scores <- lapply(1:p, function(x) colMeans((t(pred[[x]]) - dat[,x])^2))

  } else {

    stop("object class not supported (must be predict.estimate or post.pred)")

  }


  # returned object
  returned_object <- list(scores = scores,
                          type = type,
                          metric = "mse")

  class(returned_object) <- "metric"

  return(returned_object)

}

#' Mean Absolute Error
#' @name mae
#' @inheritParams mse
#' @return object of class \code{metric}
#' @export
#'
#' @examples
#' # data
#' Y <- subset(tas, gender == "M")[,-ncol(tas)]
#'
#' # fit model
#' fit <- estimate(Y)
#'
#' # predict (note summary = FALSE)
#' pred <- predict(fit, iter = 50, summary = FALSE)
#'
#' mae(pred)
mae <- function(object, ...){

  # data
  dat <- object$dat

  # number of variables
  p <- ncol(object$dat)

  # predictions
  pred <- object$pred


  if(is(object, "predict")|
     is(object, "fitted")){

    type <- class(object)

    # scores
    scores <- lapply(1:p, function(x) colMeans(abs(t(pred[[x]]) - dat[,x])))

  } else if (is(object, "post.pred")){

    type <- class(object)

    # scores
    scores <- lapply(1:p, function(x) colMeans(abs(t(pred[[x]]) - dat[,x])))

  } else {

    stop("object class not supported (must be predict.estimate or post.pred)")

  }

  returned_object <- list(scores = scores,
                          metric = "mae",
                          type = type)

  class(returned_object) <- c("BGGM", "metric", "estimate")
  return(returned_object)
}

#' Root Mean Squared Error
#' @name rmse
#' @inheritParams mse
#' @return object of class \code{metric}
#' @export
#'
#' @examples
#' # data
#' Y <- subset(tas, gender == "M")[,-ncol(tas)]
#'
#' # fit model
#' fit <- estimate(Y)
#'
#' # predict (note summary = FALSE)
#' pred <- predict(fit, iter  = 50, summary = FALSE)
#'
#' rmse(pred)
rmse <- function(object, ...){

  # data
  dat <- object$dat

  # number of variables
  p <- ncol(object$dat)

  # predictions
  pred <- object$pred

  if(class(object) == "predict.estimate"|
     class(object) == "fitted.estimate"){

    type <- class(object)

    # scores
    scores <- lapply(1:p, function(x) sqrt(colMeans((t(pred[[x]]) - dat[,x])^2)))

  } else if (class(object) == "post.pred"){

    type <- class(object)

    # scores
    scores <- lapply(1:p, function(x) sqrt(colMeans((t(pred[[x]]) - dat[,x])^2)))

  } else {

    stop("object class not supported (must be predict.estimate or post.pred)")

  }

  # returned object
  returned_object <- list(scores = scores,
                          metric = "rmse",
                          type = type)
  class(returned_object) <- "metric"
  return(returned_object)

}

#' Mean Absolute Percentage Error
#' @name mape
#' @inheritParams mse
#'
#' @return object of class \code{metric}
#' @export
#'
#' @examples
#' # data
#' Y <- subset(tas, gender == "M")[,-ncol(tas)]
#'
#' # fit model
#' fit <- estimate(Y)
#'
#' # predict (note summary = FALSE)
#' pred <- predict(fit, iter = 50, summary = FALSE)
#'
#' mape(pred)
mape <- function(object, ...){

  # data
  dat <- object$dat

  # number of variables
  p <- ncol(object$dat)

  # predictions
  pred <- object$pred


  if(class(object) == "predict.estimate" |
     class(object) == "fitted.estimate"){

    type <- class(object)

    # scores
    scores <- lapply(1:p, function(x) colMeans((abs((dat[,x] - t(pred[[x]]))/dat[,x]))))

  } else if (class(object) == "post.pred"){

    type <- class(object)

    # scores
    scores <- lapply(1:p, function(x) colMeans((abs((dat[,x] - t(pred[[x]]))/dat[,x]))))

  } else {

    stop("object class not supported (must be predict.estimate or post.pred)")

  }

  # returned object
  returned_object <- list(scores = scores,
                          metric = "mape",
                          type = type)

  class(returned_object) <- c("BGGM",  "metric", "estimate")
  return(returned_object)

}


#' Summary Method for \code{metric} Objects
#'
#' @param object object of class \code{metric}
#' @param cred credible interval
#' @param ... currently ignored
#' @export
summary.metric <- function(object, cred = 0.95, ...){

  lb <- (1 - cred) / 2

  ub <- 1 - lb

  p <- length(object$scores)

  # identity matrix
  I_p <- diag(p)

 # column names
  cn <-  colnames(object$Y)
  if(is.null(cn)){
    mat_names <- 1:p
  } else {
    mat_names <- cn

  }


  iter <- length(object$scores[[1]])

  dat_summ <- data.frame(Node = mat_names,
                         Post.mean  = round(sapply(object$scores, mean), 3) ,
                         Post.sd = round(sapply(object$scores, sd),3),
                         Cred = round(t(sapply(object$scores,
                                         quantile,
                                         c(lb, ub))), 3))

  dat_summ[is.na(dat_summ)] <- 0

  returned_object <- list(summary = dat_summ,
                          metric = object$metric,
                          type = object$type,
                          iter = iter,
                          cred = cred)

  class(returned_object) <- c("BGGM", "metric",
                              "estimate",
                              "summary",
                              "data.frame")
  returned_object
}



#' Plot \code{metric} Objects
#' @param x object of class \code{metric}
#' @param type \code{"error_bar"} or \code{"ridgeline"}
#' @param cred credible interval
#' @param width width of error bar end (\code{type = "error_bar"})
#' @param size point size (\code{type = "error_bar"})
#' @param color point (\code{type = "error_bar"}) or
#' tail region (\code{type = "ridgeline"} ) color
#' @param alpha transparency of tail region (\code{type = "ridgeline"})
#' @param scale overlap of densities (\code{type = "ridgeline"})
#' @param rope region of practical equivalence (only for Bayes R2 difference)
#' @param ... currently ignored
#'
#' @return \code{ggplot}
#' @importFrom reshape2 melt
#' @importFrom ggridges stat_density_ridges
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
#' # posterior predictions
#' pred <- posterior_predict(fit, iter = 500,
#'                           summary = FALSE)
#'
#' # prediction error
#' error <- mse(pred)
#'
#' # plot
#' plot(error)
#' }
plot.metric <- function(x, type = "error_bar",
                        cred = 0.95, alpha = 0.5,
                        scale = 1, width = 0,
                        size = 1, color = "blue",
                        rope = 0.1,
                        ...){

  if(type == "error_bar"){
    # summary
    summ <- summary(x, cred = cred)

    # temporary dat
    temp <- summ$summary

    if(x$metric == "bayes_R2" | x$metric == "bayes_R2_diff"){

      temp$Post.mean <- ifelse(is.nan(temp$Post.mean), 0, temp$Post.mean)

      #add ordered levelse
      temp$Node <- factor(temp$Node,
                          levels = temp$Node[(order(temp$Post.mean))],
                          labels = temp$Node[(order(temp$Post.mean))])

    } else {

      # add ordered levels
      temp$Node <- factor(temp$Node,
                          levels = (order(temp$Post.mean)),
                          labels =(order(temp$Post.mean)))

    }
    # plot
    plt <- ggplot(temp, aes(x = Node,
                            y = Post.mean))

    if(x$metric == "bayes_R2_diff"){

      plt <- plt + annotate("rect", xmin = -Inf,
                            xmax = Inf, ymin = -rope,
                            ymax =rope,
                            alpha = .1)


    }

    plt <- plt + geom_errorbar(aes(ymin =  temp[,4],
                                   ymax =  temp[,5]),
                               width = width) +
      geom_point(size = size,
                 color = color) +
      coord_flip() +
      xlab("Node") +
      ylab(x$metric) +
      ggtitle(x$type)



  } else if (type == "ridgeline"){

    lb <- (1 - cred) / 2
    ub <- 1 - lb

    dat <- reshape2::melt(x$scores)
    dat$L1 <- as.factor(dat$L1)

    dat$L1 <- factor(dat$L1,
                     labels = rev(order(tapply(dat$value,
                                               dat$L1,
                                               mean))),
                     levels = rev(order(tapply(dat$value,
                                               dat$L1,
                                               mean))))

    color <- grDevices::adjustcolor(color,
                                    alpha.f = alpha)


    plt <- ggplot(dat, aes(x = value,
                           y = as.factor(L1),
                           fill=factor(..quantile..)))

    if(x$metric == "bayes_R2_diff"){

      plt <- plt + annotate("rect", ymin = -Inf,
                            ymax = Inf, xmin = -rope,
                            xmax =rope,
                            alpha = .1)


    }

    plt <- plt + stat_density_ridges(rel_min_height = 0.01,
                                     scale = scale,
                                     geom = "density_ridges_gradient",
                                     calc_ecdf = TRUE,
                                     quantiles = c(lb, ub)) +
      scale_fill_manual(name = "Probability",
                        values = c(color,
                                   "#A6A6A680",
                                   color)) +
      theme(legend.position = "none") +
      ylab("Node") +
      xlab(x$metric) +
      ggtitle(x$type)

  }
  plt
}
