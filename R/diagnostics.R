#' Regression Diagnostic Plots for \code{estimate} Objects
#'
#' @description GGMs have a direct correspondence to multiple regression. Hence this function
#' provides diagnostic plots for inspecting the fitted regression models. This allows for
#' visually inspecting assumptions of the model (e.g., normality of the residuals, etc.)
#'
#' @param object object of class \code{estimate}
#' @param iter  iterations used for computing residuals and fitted values
#' @param ... currently ignored
#' @importFrom cowplot plot_grid
#'
#' @return list of \code{ggplot} objects (a plot for each node)
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
#' # diagnostic plot (iter = 10 as an example)
#' diagnostics(fit, iter = 10)
#'}
diagnostics <- function(object, iter = 500,...){


  oldw <- getOption("warn")
  options(warn = -1)

  if(class(object) != "estimate"){
    stop("object must be of class estimate")
  }

   p <- object$p
  .resid <- residuals(object, iter = iter, summary = TRUE)
  .fitted <- fitted(object, iter = iter, summary = TRUE)

  plot <- list()
  for(i in 1:p){

    dat <- data.frame(.resid = .resid[,1,i], .fitted = .fitted[,1,i])

    plot1 <- ggplot(dat, aes(.fitted, .resid)) +
      geom_point() +
      stat_smooth(method="loess", se = FALSE) +
      geom_hline(yintercept = 0, col = "red", linetype = "dashed") +
      xlab("Fitted values") +
      ylab("Residuals") +
      theme_bw() +
      ggtitle("Residuals vs. Fitted")



    plot2 <- ggplot(dat, aes(sample=.resid/sd(.resid)))+stat_qq() +
      stat_qq_line() +
      theme_bw() +
      xlab("Theoretical Quantiles") +
      ylab("Standardized Residuals") +
      ggtitle("Normal Q-Q")

    suppressMessages(
      plot3 <- ggplot(dat, aes(x = .resid)) +
        geom_histogram(color = "white", aes(y = stat(density))) +
        xlab("Residuals") +
        stat_function(
          fun = dnorm,
          args = list(mean = mean(dat$.resid), sd = sd(dat$.resid)),
          lwd = 1.5,
          col = 'red',
          alpha = 0.85,
        ) +
        theme_bw() +
        ylab("Density") +
        ggtitle("Residual Histogram")
    )


    suppressMessages(
      bottom <- cowplot::plot_grid(plot1, plot2, plot3, nrow = 1)
    )

    suppressMessages(
    plot[[i]] <- cowplot::plot_grid("", bottom, nrow = 2,
                                    labels = paste("Node", i), rel_heights = c(1, 10))
  )

  }
  options(warn = oldw)
  plot
}
