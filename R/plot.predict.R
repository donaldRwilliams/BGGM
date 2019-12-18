#' Predictive Plots for In- and Out-of-Sample Data
#'
#' @param x1 object of class \code{predict}
#' @param x2 optional object of class \code{predict} (typically using \code{test_data})
#' @param size geom_point size
#' @param color geom_point color
#' @param width geom_errorbar width
#'
#' @return A baseline \code{ggplot} object. Further customization is possible. An example is provided below.
#' @export
#' @note x1 can be a \code{predict} object that assesses predictive accuracy on in- or out-of-sample (\code{test_data}) data. x2 can be used to
#' visually compare in- and out-of sample predictive accuracy. This is accomplished by providing the respective \code{predict} objects to
#' x1 and x2. An example is provided below.
#' @examples
#'
#' dat <- BGGM::ptsd
#' fit <- estimate(dat, samples)
#'
#' select graph
#' selected <- select(fit, ci_width = 0.95)$adjacency_mat
#'
#' error <- predict(fit,
#'               selected = selected,
#'               measure = "MSE",
#'               test_data  = NULL,
#'               ci_width = 0.95,
#'               samples = 1000)
#'
#' plt <- plot(error, size = 5, color = "red", width = 1)
#'
#' # plot can be changed further with ggplot--e.g., changing theme
#'
#' plt + theme_classic()
plot.predict <- function(x1,  x2 = NULL, size = 2, color = "red", width = .1, order = NULL){

  if(is.null(x2)){

  x <- x1

  # column for nodes
  x$summary_error$node <- 1:x$p
  # add ordered levels
  x$summary_error$node <- factor(x$summary_error$node,
                               levels = order(x$summary_error$post_mean),
                               labels = order(x$summary_error$post_mean))

  if(x$measure == "R2"){
  # column for nodes
  x$summary_error$node <- 1:x$p

  # add ordered levels
  x$summary_error$node <- factor(x$summary_error$node,
                                 levels = order(x$summary_error$post_mean),
                                 labels = order(x$summary_error$post_mean))

  plt <- ggplot(x$summary_error, aes(x = node,
                                     y = post_mean)) +
    geom_errorbar(aes(ymin =  x$summary_error[,3],
                      ymax = x$summary_error[,4]), width = width) +
    geom_point(size = size,
               color = color) +
    coord_flip() +
    xlab("Node") +
    ylab(expression("Bayesian  "* R^2)) +
    theme_classic()
}
  if(x$measure == "MSE"){
    # column for nodes
    x$summary_error$node <- 1:x$p

    # add ordered levels
    x$summary_error$node <- factor(x$summary_error$node,
                                 levels = rev(order(x$summary_error$post_mean)),
                                 labels = rev(order(x$summary_error$post_mean)))


    plt <-  ggplot(x$summary_error, aes(x = node,
                                      y = post_mean)) +
    geom_errorbar(aes(ymin =  x$summary_error[,3],
                      ymax =  x$summary_error[,4]), width = width) +
    geom_point(size = size,
               color = color) +
    coord_flip() +
    xlab("Node") +
    ylab("Mean Squared Error")
    }
  }
   # for test and training plot
   if(!is.null(x2)){

     # if R2 do the following
     if(x1$measure == "R2"){

       # check that x2 is out-of-sample
       if(is.null(x2$test_data)) stop("x2 must be out-of-sample predictive accuracy")

       x1$summary_error$node <- 1:x1$p
       x1$summary_error$Error <- "train"

       x2$summary_error$node <- 1:x2$p
       x2$summary_error$Error <- "test"

       if(is.null(order) || order == "train") {
       x1$summary_error$node <- factor(x1$summary_error$node,
                                    levels = order(x1$summary_error$post_mean),
                                    labels = order(x1$summary_error$post_mean))

       temp_dat <- rbind.data.frame(x1$summary_error, x2$summary_error)
       } else {
         x2$summary_error$node <- factor(x2$summary_error$node,
                                         levels = order(x2$summary_error$post_mean),
                                         labels = order(x2$summary_error$post_mean))

         temp_dat <- rbind.data.frame(x2$summary_error, x1$summary_error)
         }





   plt <-  ggplot(temp_dat, aes(x = node,
                                         y = post_mean,
                                  group = Error)) +

       geom_errorbar(aes(ymin =  temp_dat[,3],
                         ymax =  temp_dat[,4]),
                     width = width,
                     position = position_dodge(1),
                     color = "grey75") +
       geom_point(position = position_dodge(1),
                  aes(color = Error),
                  size = 4) +
       coord_flip() +
       xlab("Node") +
       ylab(expression("Bayesian  "* R^2)) +
       theme_classic() +
       theme(panel.grid = element_blank())

     }

     if(x1$measure == "MSE"){

       # check that x2 is out-of-sample
       if(is.null(x2$test_data)) stop("x2 must be out-of-sample predictive accuracy")

       x1$summary_error$node <- 1:x1$p
       x1$summary_error$Error <- "Train"

       x2$summary_error$node <- 1:x2$p
       x2$summary_error$Error <- "test"

       if(is.null(order) || order == "train") {
         x1$summary_error$node <- factor(x1$summary_error$node,
                                         levels = order(x1$summary_error$post_mean),
                                         labels = order(x1$summary_error$post_mean))

         temp_dat <- rbind.data.frame(x1$summary_error, x2$summary_error)
       } else {
         x2$summary_error$node <- factor(x2$summary_error$node,
                                         levels = order(x2$summary_error$post_mean),
                                         labels = order(x2$summary_error$post_mean))

         temp_dat <- rbind.data.frame(x2$summary_error, x1$summary_error)
       }





       plt <-  ggplot(temp_dat, aes(x = node,
                                    y = post_mean,
                                    group = Error)) +

         geom_errorbar(aes(ymin =  temp_dat[,3],
                           ymax =  temp_dat[,4]),
                       width = width,
                       position = position_dodge(1),
                       color = "grey75") +
         geom_point(position = position_dodge(1),
                    aes(color = Error),
                    size = 4) +
         coord_flip() +
         xlab("Node") +
         ylab(expression("Mean Squared Error")) +
         theme_classic() +
         theme(panel.grid = element_blank())

     }
  }
plt
}



