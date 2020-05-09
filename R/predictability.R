#' Predictability: Bayesian Variance Explained (R2)
#'
#' @name predictability
#'
#' @description  Compute nodewise predictability or  Bayesian variance explained \insertCite{@R2 @gelman_r2_2019}{BGGM}.
#'               In the context of GGMs, this method was described in \insertCite{Williams2019;textual}{BGGM}.
#'
#'
#' @param object object of class \code{estimate} or \code{explore}
#'
#' @param select logical. Should the graph be selected ? The default is currently \code{FALSE}.
#'
#' @param cred numeric. credible interval between 0 and 1  (default is 0.95) that is used for selecting the graph.
#'
#' @param BF_cut numeric. evidentiary threshold (default is 3).
#'
#' @param iter interger. iterations (posterior samples) used for computing R2.
#'
#' @param ... currently ignored.
#'
#' @return object of classes \code{bayes_R2} and \code{metric}
#'
#' @note
#'
#'
#' \strong{Binary and Ordinal Data}:
#'
#' R2 is computed from the latent data.
#'
#'
#' \strong{Mixed Data}:
#'
#' The mixed data approach is somewhat ad-hoc \insertCite{@see for example p. 277 in  @hoff2007extending;textual}{BGGM}. This
#' is becaue uncertainty in the ranks is not incorporated, which means that variance explained is computed from
#' the 'empirical' \emph{CDF}.
#'
#' \strong{Model Selection}:
#'
#' Currently the default to include all nodes in the model when computing R2. This can be changed (i.e., \code{select = TRUE}), which
#' then sets those edges not detected to zero. This is accomplished by subsetting the correlation matrix according to each neighborhood
#' of relations.
#'
#' @references
#' \insertAllCited{}
#'
#' @examples
#' \donttest{
#' #####################
#' #### continuous #####
#' #####################
#' # data
#' Y <- BGGM::ptsd
#'
#' # estimate the model
#' fit <- estimate(Y, iter = 1000)
#'
#' # predictability
#' r2 <- predictability(fit,
#'                      select = TRUE,
#'                      iter = 1000)
#'
#' # print summary
#' r2
#'
#' # plot
#' plot(r2)
#'
#'
#' #####################
#' ####### mixed #######
#' #####################
#'
#' # estimate the model
#' fit <- estimate(Y,
#'                 type = "mixed",
#'                 iter = 1000)
#'
#' # predictability
#' r2 <- predictability(fit,
#'                      select = TRUE,
#'                      iter = 1000)
#'
#' # print summary
#' r2
#'
#' # plot
#' plot(r2)
#' }
#' @export
predictability <- function(object,
                           select = FALSE,
                           cred = 0.95,
                           BF_cut = 3,
                           iter = NULL){


  if(object$type == "continuous"){

    Y <- as.matrix(scale(object$Y))

    } else if(object$type == "binary"){

      Y <- binary_latent_helper(object$Y+1)

    } else if(object$type == "ordinal"){

      # latent data
      Y <- ordinal_latent_helper(object$Y, object$post_samp$thresh)

      } else {

        # latent data
        Y <- rank_helper(object$Y)$z0_start
    }

   # nodes
  p <- ncol(Y)

  # observations
  n <- nrow(Y)


  if(is.null(iter)){
    iter <- 1000
    }

  # correlations
  cors <- pcor_to_cor(object)$R[,,1:iter]

  # not conditional on selected model
  if(isFALSE(select)){

    # progress bar
    pb <- utils::txtProgressBar(min = 0, max = p, style = 3)


    r2 <-  lapply(1:p,   function(x)  {

      # computed from selected model
      r2_p <- .Call(
        "_BGGM_predictability_helper",
        Y[, -x],
        y =  Y[, x],
        XX =  cors[-x,-x, ],
        Xy =  cors[x, -x, ],
        n = n,
        iter = iter
      )$r2

      utils::setTxtProgressBar(pb, x)

      r2_p

    })



  } else {


    if(is(object, "estimate") & is(object, "default")){

      # select model
      sel <- select(object, cred = cred)
      # adjacency
      adj <- sel$adj

    } else if(is(object, "explore") &  is(object, "default")){

      sel <- select(object, BF_cut = BF_cut)

      adj <- sel$Adj_10

    }

    # progress bar
    pb <- utils::txtProgressBar(min = 0, max = p, style = 3)

    # R2
    r2 <- lapply(1:p, function(x)  {

      # non selected return zero
      if(sum(adj[x,])  == 0 ){

       r2_p <- 0

       r2_p

        # a neighborhood exists
      } else {

        # neighborhood
        selected <- which(adj[x,] == 1)

        # check length 1 (return correlation)
        if(length(selected) == 1){

          r2_p <- cors[x, selected,]

          r2_p


          # more than one relation: call c++
        } else {

          # computed from selected model
          r2_p <- .Call(
            "_BGGM_predictability_helper",
            Y[, selected],
            y =  Y[, x],
            XX =  cors[selected, selected, ],
            Xy =  cors[x, selected, ],
            n = n,
            iter = iter
          )$r2
        }
      }

      utils::setTxtProgressBar(pb, x)

      r2_p

    })
  }

  # R2
  scores <- lapply(1:p, function(x) {

    r2_new <- r2[[x]]
    if(length(r2_new) > 0){

      r2_new[r2_new > 0]
    }

    })

  # returned object
  returned_object <- list(scores = scores,
                          type = "post.pred",
                          metric = "bayes_R2",
                          cred = cred,
                          BF_cut = BF_cut,
                          data_type = object$type,
                          Y = Y)

  class(returned_object) <- c("BGGM",
                              "predictability",
                              "metric",
                              "R2",
                              "estimate")
  return(returned_object)
}




#' Summary Method for \code{predictability} Objects
#'
#' @param object an object of class \code{predictability}
#'
#' @param cred  credible interval
#'
#' @param ... currently ignored
#'
#' @export

summary.predictability <- function(object, cred = 0.95, ...){

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
                          data_type = object$data_type,
                          cred = cred)

  class(returned_object) <- c("BGGM",
                              "predictability",
                              "metric",
                              "estimate",
                              "summary",
                              "data.frame"
                              )
  returned_object
}




print_summary_metric <- function(x, digits = 2,...){
  cat("BGGM: Bayesian Gaussian Graphical Models \n")
  cat("--- \n")
  if(x$metric == "bayes_R2"){
    cat("Metric:", "Bayes R2\n")
  } else if(x$metric == "bayes_R2_diff"){
    cat("Metric:", "Bayes R2 Difference \n")
  } else {
    cat("Metric:", x$metric, "\n")

  }
  cat("Type:", x$data_type, "\n")
  # cat("Credible Interval:", x$cred, "\n")
  cat("--- \n")
  cat("Estimates:\n\n")
  dat <- x$summary
  colnames(dat) <- c(colnames(dat)[1:3], "Cred.lb", "Cred.ub")
  print(as.data.frame(dat),
        row.names = FALSE)
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
plot.predictability <- function(x, type = "error_bar",
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
      ylab(x$metric)



  } else if (type == "ridgeline"){

    # intervals
    lb <- (1 - cred) / 2

    ub <- 1 - lb

    summ <- summary(x)

    dat <- reshape2::melt(x$scores)

    dat$L1 <- as.factor(dat$L1)

    dat$L1 <- factor(dat$L1,
                     labels = (order(tapply(dat$value,
                                            dat$L1,
                                            mean))),
                     levels = (order(tapply(dat$value,
                                            dat$L1,
                                            mean))))

    color <- grDevices::adjustcolor(color,
                                    alpha.f = alpha)

    plt <- ggplot(dat, aes(x = as.numeric(value),
                           y = as.factor(L1),
                           fill=factor(stat(quantile)))) +
      stat_density_ridges(rel_min_height = 0.01,
                          geom = "density_ridges_gradient",
                          calc_ecdf = TRUE,
                          quantiles = c(lb, ub),
                          scale = scale) +
      scale_fill_manual(name = "Probability",
                        values = c(color,
                                  "#A6A6A680",
                                   color)) +
      theme(legend.position = "none") +
      ylab("Node") +
      xlab(x$metric) +
      scale_y_discrete(labels = summ$summary$Node[order(summ$summary$Post.mean)])


  } else {

    stop("type not supported. must be 'error_bar' or 'ridgeline'.")
  }
  plt
}

